#/c/Program\ Files/R/R-devel/bin/R.exe CMD INSTALL -l ~/AppData/Local/R/win-library/4.5/ SeqExpMatch/
pacman::p_load(PTE, datasets, mlbench, dplyr, ggplot2, gridExtra, profvis)
library(SeqExpMatch)
options(error=recover)
options(warn=2)
# set.seed(1984)


nsim_exact_test = 301
num_cores = 1
missing_data_prop = 0.000
prob_of_adding_response = 0.3
max_n_dataset = 200

data(continuous_example)
data(Glass)
data(Sonar)
data(Soybean)
data(BreastCancer)
airquality_subset = rbind(airquality[1, ], airquality %>% slice_sample(n = max_n_dataset - 1, replace = TRUE))
diamonds_subset = ggplot2::diamonds %>% slice_sample(n = max_n_dataset, replace = TRUE)
boston_subset = MASS::Boston %>%        slice_sample(n = max_n_dataset, replace = TRUE)
cars_subset = MASS::Cars93 %>%          slice_sample(n = max_n_dataset, replace = TRUE) %>% select(-Make, -Model)
glass_subset = Glass %>%                slice_sample(n = max_n_dataset, replace = TRUE) %>% mutate(Type = as.numeric(Type) - 1)
pima_subset = MASS::Pima.tr2 %>%        slice_sample(n = max_n_dataset, replace = TRUE) %>% mutate(type = ifelse(type == "Yes", 1, 0))
sonar_subset = Sonar %>%                slice_sample(n = max_n_dataset, replace = TRUE) %>% mutate(Class = ifelse(Class == "M", 1, 0))
soybean_subset = Soybean %>%            slice_sample(n = max_n_dataset, replace = TRUE) %>% mutate(Class = ifelse(Class == "brown-spot", 1, 0))
cancer_subset = BreastCancer %>%        slice_sample(n = max_n_dataset, replace = TRUE) %>% select(-Id) %>% mutate(Class = ifelse(Class == "malignant", 1, 0))

finagle_different_responses_from_continuous = function(y_cont){
  list(
    continuous = y_cont,
    incidence =  ifelse(y_cont > median(y_cont), 1, 0),
    proportion = (y_cont - min(y_cont) + .Machine$double.eps) / max(y_cont - min(y_cont) + 2 * .Machine$double.eps),
    count =      round(y_cont - min(y_cont)),
    survival =   y_cont - min(y_cont) + .Machine$double.eps
  )
}

datasets = list(
  # airquality = list(
  #   X = data.table(airquality_subset %>% select(-Wind)),
  #   responses = finagle_different_responses_from_continuous(airquality_subset$Wind)
  # ),
  # cars = list(
  #   X = data.table(cars_subset %>% select(-Price)),
  #   responses = finagle_different_responses_from_continuous(cars_subset$Price)
  # ),
  glass = list(
    X = data.table(glass_subset %>% select(-Type)),
    responses = list(count = glass_subset$Type)
  ),
  pima = list(
    X = data.table(pima_subset %>% select(-type)),
    responses = list(incidence = pima_subset$type)
  ),
  sonar = list(
    X = data.table(sonar_subset %>% select(-Class)),
    responses = list(incidence = sonar_subset$Class)
  ),
  soybean = list(
    X = data.table(soybean_subset %>% select(-Class)),
    responses = list(incidence = soybean_subset$Class)
  ),
  cancer = list(
    X = data.table(cancer_subset %>% select(-Class)),
    responses = list(incidence = cancer_subset$Class)
  ),
  pte_example = list(
    X = data.table(continuous_example$X %>% 
      mutate(treatment = as.factor(treatment)) %>% 
      mutate(x4 = as.character(x4))),
    responses = finagle_different_responses_from_continuous(continuous_example$y),
    beta_T = list(
      continuous = 0,
      incidence =  0,
      proportion = 0,
      count =      0,
      survival =   0      
    )
  ),
  diamonds = list(
    X = data.table(diamonds_subset %>%
      select(-price) %>%
      mutate(clarity = as.character(clarity)) %>%
      mutate(color = as.character(color)) %>%
      mutate(cut = as.character(cut))),
    responses = finagle_different_responses_from_continuous(diamonds$price),
    beta_T = list(
      continuous = 0,
      incidence =  0,
      proportion = 0,
      count =      0,
      survival =   0
    )
  ),
  boston = list(
    X = data.table(boston_subset %>% select(-medv)),
    responses = finagle_different_responses_from_continuous(boston_subset$medv),
    beta_T = list(
      continuous = 0,
      incidence =  0,
      proportion = 0,
      count =      0,
      survival =   0
    )
  )
)


res = data.frame()

#test all response_types and all designs
# profvis({
for (dataset_name in names(datasets)){
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
    # if (response_type != "survival"){
    #   next
    # }
    
    y = datasets[[dataset_name]][["responses"]][[response_type]]
    beta_T =  if (is.null(datasets[[dataset_name]][["beta_T"]])){
                0
              } else if (is.null(datasets[[dataset_name]][["beta_T"]][[response_type]])){
                0
              } else {
                datasets[[dataset_name]][["beta_T"]][[response_type]]
              }
    
    
    importance_plots = list()
    for (d in c("CRD", "iBCRD", "Efron", "Atkinson", "KK21", "KK21stepwise")){ #"CRD", "iBCRD", "Efron", "Atkinson", "KK14", "KK21", "KK21stepwise"
      cat("\ndataset_name =", dataset_name, "response_type =", response_type, "d =", d, "\n")
      
      seq_des_obj = SeqDesign$new(n, d, response_type = response_type, verbose = FALSE)
      
      response_added = array(FALSE, n)
      
      # profvis({
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
      
      for (t in which(!response_added)){
        suppressWarnings(seq_des_obj$add_subject_response(t, y[t] + beta_T * w_t, dead = 0))
      }
      # cat("experiment completed?", seq_des_obj$check_experiment_completed(), "\n")
      
      if (grepl("KK", d)){
        stats = seq_des_obj$matching_statistics()
        print(paste("  KK stats: prop matches:", stats$prop_subjects_matched, "# remaining in reservoir:", stats$num_subjects_remaining_in_reservoir))
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
        tryCatch({
          plot(grid.arrange(importance_plots[["KK21"]], importance_plots[["KK21stepwise"]], nrow = 2))  
        }, error = function(e){})
      }
      
      # let's do some inference
      for (test_type in c("MLE-or-KM-based")){ #"MLE-or-KM-based", "randomization-exact"
        for (estimate_type in c(
          ########################################### CONTINUOUS
          "continuous_simple_mean_difference",
          "continuous_regression_with_covariates",
          "continuous_KK_compound_mean_difference",  	
          "continuous_KK_compound_multivariate_regression",
          "continuous_KK_regression_with_covariates_with_matching_dummies",
          "continuous_KK_regression_with_covariates_with_random_intercepts",
          ########################################### INCIDENCE
          "incidence_simple_mean_difference",
          "incidence_simple_log_odds",	
          "incidence_logistic_regression",
          #"incidence_KK_compound_univariate_logistic_regression",
          #"incidence_KK_compound_multivariate_logistic_regression",	
          "incidence_KK_multivariate_logistic_regression_with_matching_dummies",	
          "incidence_KK_multivariate_logistic_regression_with_random_intercepts_for_matches",
          ########################################### PROPORTION
          "proportion_simple_mean_difference",
          "proportion_simple_logodds_regression",
          #"proportion_beta_regression",
          #"proportion_KK_compound_univariate_beta_regression",
          #"proportion_KK_compound_multivariate_beta_regression",
          #"proportion_KK_multivariate_beta_regression_with_matching_dummies",
          ########################################### COUNT
          "count_simple_mean_difference",
          "count_univariate_negative_binomial_regression",
          "count_multivariate_negative_binomial_regression",
          #"count_KK_compound_univariate_negative_binomial_regression",	
          #"count_KK_compound_multivariate_negative_binomial_regression",
          "count_KK_multivariate_negative_binomial_regression_with_matching_dummies",
          "count_KK_multivariate_negative_binomial_regression_with_random_intercepts_for_matches",
          ########################################### SURVIVAL
          "survival_simple_median_difference",	
          "survival_simple_restricted_mean_difference",
          "survival_univariate_weibull_regression",	
          "survival_multivariate_weibull_regression",
          # "survival_KK_compound_univariate_weibull_regression",	
          # "survival_KK_compound_multivariate_weibull_regression",
          #"survival_KK_multivariate_weibull_regression_with_matching_dummies",
          "survival_univariate_coxph_regression",	
          "survival_multivariate_coxph_regression",		
          "survival_KK_multivariate_coxph_regression_with_matching_dummies",		
          "survival_KK_multivariate_coxph_regression_with_random_intercepts_for_matches"				
        )){
          #only do inference by appropriate response type
          if (strsplit(estimate_type, "_")[[1]][1] != response_type){
            next
          }
          #only do inference for KK if it's a KK design
          if (grepl("KK", estimate_type) & !grepl("KK", d)){
            next
          }
          seq_des_inf_obj = SeqDesignInference$new(seq_des_obj, estimate_type = estimate_type, test_type = test_type, num_cores = num_cores)
          beta_hat_T = seq_des_inf_obj$compute_treatment_estimate()
          cat("  beta_hat_T =", beta_hat_T, "\n")
          pval = seq_des_inf_obj$compute_two_sided_pval_for_treatment_effect(nsim_exact_test = nsim_exact_test)
          cat("  pval =", pval, "\n")
          
          if (test_type == "MLE-or-KM-based"){
            ci = seq_des_inf_obj$compute_confidence_interval(nsim_exact_test = nsim_exact_test)
            if (test_type == "MLE-or-KM-based" | (test_type == "randomization-exact" & response_type == "continuous")){
              cat("  95% CI for betaT = [", paste(round(ci, 4), collapse = ", "), "] width =", round(ci[2] - ci[1], 4), "\n")
            }
          } else {
            ci = c(NA, NA)
          }
          
          res = rbind(res, data.frame(
            dataset_name = dataset_name,
            response_type = response_type,
            design = d,
            test_type = test_type,
            estimate_type = estimate_type,
            beta_hat_T = beta_hat_T,
            ci_a = ci[1],
            ci_b = ci[2],
            p_val = pval
          ))
        }
      }
    }
  }
}
# })



#####comparisons
###optimal PM on covariates seeing all x's up front
###optimal PM on mu_i's