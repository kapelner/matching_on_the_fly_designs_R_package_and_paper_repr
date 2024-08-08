#/c/Program\ Files/R/R-devel/bin/R.exe CMD INSTALL -l ~/AppData/Local/R/win-library/4.5/ SeqExpMatch/
pacman::p_load(PTE, dplyr, ggplot2, gridExtra, profvis)
library(SeqExpMatch)
options(error=recover)
set.seed(1984)

diamonds_subset = ggplot2::diamonds[sample(1 : 20000, 200), ]
data(continuous_example)

finagle_different_responses_from_continuous = function(y_cont){
  list(
    continuous = y_cont,
    incidence =  ifelse(y_cont > median(y_cont), 1, 0),
    proportion = (y_cont - min(y_cont)) / max(y_cont - min(y_cont)),
    count =      round(y_cont - min(y_cont)),
    survival =   y_cont - min(y_cont) + .Machine$double.eps
  )  
}
datasets = list(
  pte_example = list(
    X = continuous_example$X %>% 
      mutate(treatment = as.factor(treatment)) %>% 
      mutate(x4 = as.character(x4)),
    responses = finagle_different_responses_from_continuous(continuous_example$y),
    beta_T = list(
      continuous = 2,
      incidence =  0,
      proportion = 0,
      count =      2,
      survival =   2      
    )
  ),
  diamonds = list(
    X = diamonds_subset %>% 
      select(-price) %>%
      mutate(clarity = as.character(clarity)) %>%
      mutate(color = as.character(color)) %>%
      mutate(cut = as.character(cut)),
    responses = finagle_different_responses_from_continuous(diamonds$price),
    beta_T = list(
      continuous = 2,
      incidence =  0,
      proportion = 0,
      count =      2,
      survival =   2      
    )
  ),
  boston = list(
    X = MASS::Boston[1 : 200, ] %>% select(-medv),
    responses = finagle_different_responses_from_continuous(MASS::Boston[1 : 200, ]$medv),
    beta_T = list(
      continuous = 2,
      incidence =  0,
      proportion = 0,
      count =      2,
      survival =   2      
    )
  )
)

nsim_exact_test = 301
num_cores = 1
missing_data_prop = 0
prob_of_adding_response = 0.5

res = list()

#test all response_types and all designs
for (dataset_name in names(datasets)){
  res[[dataset_name]] = list()
  Xorig = datasets[[dataset_name]]$X
  n = nrow(Xorig)
  
  #add some missing data
  X = Xorig
  for (i in 2 : n){
    for (j in 1 : ncol(Xorig)){
      if (runif(1) <= missing_data_prop){
        X[i, j] = NA
      }
    }
  }
  
  for (response_type in names(datasets[[dataset_name]][["responses"]])){ #"continuous", "incidence", "count", "proportion", "survival"
    y = datasets[[dataset_name]][["responses"]][[response_type]]
    beta_T = datasets[[dataset_name]][["beta_T"]][[response_type]]
    res[[dataset_name]][[response_type]] = list()
    
    importance_plots = list()
    for (d in c("CRD", "iBCRD", "Efron", "Atkinson", "KK14", "KK21", "KK21stepwise")){ #"CRD", "iBCRD", "Efron", "Atkinson", "KK14", "KK21", "KK21stepwise"
      cat("\ndataset_name =", dataset_name, "response_type =", response_type, "d =", d, "\n")
      res[[dataset_name]][[response_type]][[d]] = list()
      
      seq_des_obj = SeqDesign$new(n, d, response_type = response_type, verbose = FALSE)
      
      # profvis({
      response_added = array(FALSE, n)
      for (t in 1 : n){
        w_t = seq_des_obj$add_subject_to_experiment_and_assign(X[t, ])
        if (t %% 50 == 0){
          seq_des_obj$print_current_subject_assignment()
        }
        if (runif(1) < prob_of_adding_response){
          seq_des_obj$add_subject_response(t, y[t] + beta_T * w_t, dead = 1)
          response_added[t] = TRUE
        }
      }    
      # })
      
      # profvis({
      for (t in which(!response_added)){
        suppressWarnings(seq_des_obj$add_subject_response(t, y[t], dead = 0))
      }
      # })
      # cat("experiment completed?", seq_des_obj$check_experiment_completed(), "\n")
      
      if (grepl("KK", d)){
        print(seq_des_obj$matching_statistics())
        if (grepl("KK21", d)){
          importance_plots[[d]] = ggplot(data.frame(x = names(seq_des_obj$covariate_weights), imp = seq_des_obj$covariate_weights)) +
            aes(x = x, y = imp) +
            geom_bar(stat = "identity") +
            ggtitle(paste("Design:", d, "response:", response_type)) + 
            theme(axis.text.x = element_text(angle = 90)) + 
            ylim(0, 0.6)
        }
      }
      
      if (d == "KK21stepwise"){
        grid.arrange(importance_plots[["KK21"]], importance_plots[["KK21stepwise"]], nrow = 2)
      }
      
      # let's do some inference
      for (test_type in c("MLE-or-KM-based", "randomization-exact")){ #"MLE-or-KM-based", "randomization-exact"
        res[[dataset_name]][[response_type]][[d]][[test_type]] = list()
        for (estimate_type in c("mean_difference", "median_difference", "default_regression")){
          if (estimate_type == "median_difference" & response_type != "survival"){
            next
          }
          seq_des_inf_obj = SeqDesignInference$new(seq_des_obj, estimate_type = estimate_type, test_type = test_type, num_cores = num_cores)
          beta_hat_T = seq_des_inf_obj$compute_treatment_estimate()
          pval = seq_des_inf_obj$compute_two_sided_pval_for_treatment_effect(nsim_exact_test = nsim_exact_test)
          cat("  beta_hat_T =", beta_hat_T, "\n")
          # if (test_type == "MLE-or-KM-based" | (test_type == "randomization-exact" & response_type == "continuous")){
          #   cat("  95% CI for betaT = [", paste(seq_des_inf_obj$compute_confidence_interval(nsim_exact_test = nsim_exact_test), collapse = ", "), "]\n")
          # }
          cat("  pval =", pval, "\n")
          
          res[[dataset_name]][[response_type]][[d]][[test_type]][["beta_hat_T"]] = beta_hat_T
          res[[dataset_name]][[response_type]][[d]][[test_type]][["pval"]] = pval
        }
      }
    }
  }
}
