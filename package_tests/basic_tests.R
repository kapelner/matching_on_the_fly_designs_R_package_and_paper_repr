#/c/Program\ Files/R/R-devel/bin/R.exe CMD INSTALL -l ~/AppData/Local/R/win-library/4.5/ SeqExpMatch/
pacman::p_load(PTE, dplyr, ggplot2, profvis)
library(SeqExpMatch)
options(error=recover)
set.seed(1984)

diamonds_subset = ggplot2::diamonds[sample(1 : 20000, 200), ]
data(continuous_example)

datasets = list(
  pte_example = list(
    X = continuous_example$X %>% 
      mutate(treatment = as.factor(treatment)) %>% 
      mutate(x4 = as.character(x4)),
    y = continuous_example$y
  ),
  diamonds = list(
    X = diamonds_subset %>% select(-price),
    y = diamonds_subset$price
  ),
  boston = list(
    X = MASS::Boston[1 : 200, ] %>% select(-medv),
    y = MASS::Boston[1 : 200, ]$medv
  )
)
  


nsim_exact_test = 301
num_cores = 1
missing_data_prop = 0.002
prob_of_adding_response = 0.5


#test all response_types and all designs
for (dataset_name in names(datasets)){
  
  Xorig = datasets[[dataset_name]]$X
  n = nrow(Xorig)
  
  for (response_type in c("continuous", "incidence", "count", "proportion", "survival")){ #"continuous", "incidence", "count", "proportion", "survival"
    for (d in c("CRD", "iBCRD", "Efron", "Atkinson", "KK14", "KK21", "KK21stepwise")){ #"CRD", "iBCRD", "Efron", "Atkinson", "KK14", "KK21", "KK21stepwise"
      #add some missing data
      X = Xorig
      for (i in 2 : n){
        for (j in 1 : ncol(Xorig)){
          if (runif(1) <= missing_data_prop){
            X[i, j] = NA
          }
        }
      }
      
      cat("\ndataset_name =", dataset_name, "response_type =", response_type, "d =", d, "\n")
      
      seq_des_obj = SeqDesign$new(n, d, response_type = response_type, verbose = FALSE)
      
      y_cont = datasets[[dataset_name]]$y
      y_cont_pos = y_cont - min(y_cont) + .Machine$double.eps
      y = switch(response_type,
            continuous = y_cont,
            incidence = ifelse(y_cont > median(y_cont), 1, 0),
            proportion = y_cont_pos / (max(y_cont_pos) + 2 * .Machine$double.eps),
            count = round(y_cont_pos),
            survival = y_cont_pos
      )
      # profvis({
      response_added = array(FALSE, n)
      for (t in 1 : n){
        seq_des_obj$add_subject_to_experiment_and_assign(X[t, ])
        if (t %% 50 == 0){
          seq_des_obj$print_current_subject_assignment()
        }
        if (runif(1) < prob_of_adding_response){
          seq_des_obj$add_subject_response(t, y[t], dead = 1)
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
          plot(ggplot(data.frame(x = 1 : length(seq_des_obj$covariate_weights), imp = seq_des_obj$covariate_weights)) +
                 aes(x = x, y = imp) +
                 geom_bar(stat = "identity"))
        }
      }
      
      # let's do some inference
      for (test_type in c("MLE-or-KM-based")){ #"MLE-or-KM-based", "randomization-exact"
        for (estimate_type in c("mean_difference", "median_difference", "default_regression")){
          if (estimate_type == "median_difference" & response_type != "survival"){
            next
          }
          seq_des_inf_obj = SeqDesignInference$new(seq_des_obj, estimate_type = estimate_type, test_type = test_type, num_cores = num_cores)
          cat("  beta_hat_T =", seq_des_inf_obj$compute_treatment_estimate(), "\n")
          # if (test_type == "MLE-or-KM-based" | (test_type == "randomization-exact" & response_type == "continuous")){
          #   cat("  95% CI for betaT = [", paste(seq_des_inf_obj$compute_confidence_interval(nsim_exact_test = nsim_exact_test), collapse = ", "), "]\n")
          # }
          cat("  pval =", seq_des_inf_obj$compute_two_sided_pval_for_treatment_effect(nsim_exact_test = nsim_exact_test), "\n")
        }
      }
    }
  }
}
