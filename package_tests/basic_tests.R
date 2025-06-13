#/c/Program\ Files/R/R-devel/bin/R.exe CMD INSTALL -l ~/AppData/Local/R/win-library/4.5/ SeqExpMatch/
pacman::p_load(PTE, datasets, mlbench, AppliedPredictiveModeling, dplyr, ggplot2, gridExtra, profvis, data.table, profvis)
library(SeqExpMatch)
options(error=recover)
options(warn=2)
# set.seed(1984)


nsim_exact_test = 301
num_cores = 1
missing_data_prop = 0
prob_of_adding_response = 0.4
prob_of_uncensored_survival_observation = 0.3
max_n_dataset = 200

data(continuous_example)
data(Glass)
data(Sonar)
data(Soybean)
data(BreastCancer)
data(iris)
data(Ionosphere)
data(abalone)
data(FuelEconomy)

airquality_subset = airquality      %>% slice_sample(n = max_n_dataset, replace = TRUE)
diamonds_subset = ggplot2::diamonds %>% slice_sample(n = max_n_dataset, replace = TRUE)
boston_subset = MASS::Boston %>%        slice_sample(n = max_n_dataset, replace = TRUE)
cars_subset = MASS::Cars93 %>%          slice_sample(n = max_n_dataset, replace = TRUE) %>% select(-Make, -Model)
glass_subset = Glass %>%                slice_sample(n = max_n_dataset, replace = TRUE) %>% mutate(Type = as.numeric(Type) - 1)
pima_subset = MASS::Pima.tr2 %>%        slice_sample(n = max_n_dataset, replace = TRUE) %>% mutate(type = ifelse(type == "Yes", 1, 0))
sonar_subset = Sonar %>%                slice_sample(n = max_n_dataset, replace = TRUE) %>% mutate(Class = ifelse(Class == "M", 1, 0))
soybean_subset = Soybean %>%            slice_sample(n = max_n_dataset, replace = TRUE) %>% mutate(Class = ifelse(Class == "brown-spot", 1, 0)) %>%
  transform(plant.stand = as.numeric(plant.stand), precip = as.numeric(precip), temp = as.numeric(temp), germ = as.numeric(germ), leaf.size = as.numeric(leaf.size))
cancer_subset = BreastCancer %>%        slice_sample(n = max_n_dataset, replace = TRUE) %>% select(-Id) %>% mutate(Class = ifelse(Class == "malignant", 1, 0)) %>%
  transform(Cl.thickness = as.numeric(Cl.thickness), Cell.size = factor(Cell.size, ordered = FALSE), Cell.shape = factor(Cell.shape, ordered = FALSE), Marg.adhesion = factor(Marg.adhesion, ordered = FALSE), Epith.c.size = as.numeric(Epith.c.size))
ionosphere_subset = Ionosphere %>%      slice_sample(n = max_n_dataset, replace = TRUE) %>% select(-V1, -V2)
abalone_subset = abalone %>%            slice_sample(n = max_n_dataset, replace = TRUE)
fuel_subset = cars2012 %>%              slice_sample(n = max_n_dataset, replace = TRUE)
rm(max_n_dataset)


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
  airquality = list(
    X = data.table(airquality_subset %>% select(-Wind)),
    responses = finagle_different_responses_from_continuous(airquality_subset$Wind),
    beta_T = list(
      continuous = 1,
      incidence =  1,
      proportion = 1,
      count =      1,
      survival =   1
    )
  ),
  cars = list(
    X = data.table(cars_subset %>% select(-Price)),
    responses = finagle_different_responses_from_continuous(cars_subset$Price),
    beta_T = list(
      continuous = 1,
      incidence =  1,
      proportion = 1,
      count =      1,
      survival =   1
    )
  ),
  glass = list(
    X = data.table(glass_subset %>% select(-Type)),
    responses = list(count = glass_subset$Type),
    beta_T = list(
      count =  1
    )
  ),
  pima = list(
    X = data.table(pima_subset %>% select(-type)),
    responses = list(incidence = pima_subset$type),
    beta_T = list(
      incidence =  1
    )
  ),
  sonar = list(
    X = data.table(sonar_subset %>% select(-Class)),
    responses = list(incidence = sonar_subset$Class),
    beta_T = list(
      incidence =  1
    )
  ),
  soybean = list(
    X = data.table(soybean_subset %>% select(-Class)),
    responses = list(incidence = soybean_subset$Class),
    beta_T = list(
      incidence =  1
    )
  ),
  cancer = list(
    X = data.table(cancer_subset %>% select(-Class)),
    responses = list(incidence = cancer_subset$Class),
    beta_T = list(
      incidence =  1
    )
  ),
  pte_example = list(
    X = data.table(continuous_example$X %>%
      mutate(treatment = as.factor(treatment)) %>%
      mutate(x4 = as.character(x4))),
    responses = finagle_different_responses_from_continuous(continuous_example$y),
    beta_T = list(
      continuous = 1,
      incidence =  1,
      proportion = 1,
      count =      1,
      survival =   1
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
      continuous = 1,
      incidence =  1,
      proportion = 1,
      count =      1,
      survival =   1
    )
  ),
  boston = list(
    X = data.table(boston_subset %>% select(-medv)),
    responses = finagle_different_responses_from_continuous(boston_subset$medv),
    beta_T = list(
      continuous = 1,
      incidence =  1,
      proportion = 1,
      count =      1,
      survival =   1
    )
  ),
  iris = list(
    X = data.table(iris %>% select(-Sepal.Width)),
    responses = finagle_different_responses_from_continuous(iris$Sepal.Width),
    beta_T = list(
      continuous = 1,
      incidence =  1,
      proportion = 1,
      count =      1,
      survival =   1
    )
  ),
  ionosphere = list(
    X = data.table(ionosphere_subset %>% select(-V3)),
    responses = finagle_different_responses_from_continuous(ionosphere_subset$V3),
    beta_T = list(
      continuous = 1,
      incidence =  1,
      proportion = 1,
      count =      1,
      survival =   1
    )
  ),
  abalone = list(
    X = data.table(abalone_subset %>% select(-LongestShell)),
    responses = finagle_different_responses_from_continuous(abalone_subset$LongestShell),
    beta_T = list(
      continuous = 1,
      incidence =  1,
      proportion = 1,
      count =      1,
      survival =   1
    )
  ),
  fuel = list(
    X = data.table(fuel_subset %>% select(-EngDispl)),
    responses = finagle_different_responses_from_continuous(fuel_subset$EngDispl),
    beta_T = list(
      continuous = 1,
      incidence =  1,
      proportion = 1,
      count =      1,
      survival =   1
    )
  )
)
rm(continuous_example,Glass,Sonar,Soybean,BreastCancer,iris,Ionosphere,abalone,cars2010,cars2011,cars2012)
rm(airquality_subset,diamonds_subset,boston_subset,cars_subset,glass_subset,pima_subset,sonar_subset,soybean_subset,cancer_subset,ionosphere_subset,abalone_subset,fuel_subset)

res = data.frame()

log_odds = function(p){log(p / (1 - p))}
expit = function(log_odds){exp(log_odds) / (1 + exp(log_odds))}
compute_response_with_treatment = function(y_t, dead_t, w_t, beta_T, response_type){
  switch(response_type,
         continuous =    y_t + beta_T * w_t,
         incidence =     if (y_t == 0 & beta_T > 0 & w_t == 1){
                           rbinom(1, 1, expit(beta_T))
                         } else if (y_t == 1 & beta_T < 0 & w_t == 1){
                           rbinom(1, 1, 1 - expit(beta_T))
                         } else {
                           y_t
                         },
         proportion =    if (y_t == 1 & beta_T > 0 & w_t == 1){
                           y_t
                         } else if (y_t == 0 & beta_T < 0 & w_t == 1){
                           y_t
                         } else {
                           expit(log_odds((y_t + .Machine$double.eps) / (1 + 2 * .Machine$double.eps)) + beta_T * w_t)
                         },
         count =         round(exp(log(y_t + .Machine$double.eps) + beta_T * w_t)),
         survival =      if (dead_t == 0){
                           exp(log(y_t) + beta_T * w_t) * runif(1) #multiply the real y by a random proportion
                         } else {
                           exp(log(y_t) + beta_T * w_t)
                         }
  )
}


all_exp_settings = expand.grid(
  dataset_name = names(datasets),
  beta_T = c("zero", "nonzero"),
  response_type = c("continuous", "incidence", "count", "proportion", "survival"),
  d = c("KK21stepwise", "KK21", "Atkinson", "Efron", "iBCRD", "CRD")
)
exp_settings = data.frame()

for (i in 1 : nrow(all_exp_settings)){
  exp_setting = all_exp_settings[i, ]
  dataset_name = exp_setting$dataset_name
  beta_T = exp_setting$beta_T
  response_type = exp_setting$response_type
  d = exp_setting$d
  test_type = exp_setting$test_type
  
  if (!(response_type %in% names(datasets[[dataset_name]][["responses"]]))){
    next
  }
  exp_settings = rbind(exp_settings, exp_setting)
}
rm(i, exp_setting, dataset_name, beta_T, d, test_type, all_exp_settings)

#test all response_types and all designs

visualizations_during_testing = TRUE
#we can parallelize this loop after we're done debugging and Rcpp'ing
for (i in 1 : nrow(exp_settings)){
  exp_setting = exp_settings[i, ]
  dataset_name = as.character(exp_setting$dataset_name)
  beta_T = exp_setting$beta_T
  response_type = as.character(exp_setting$response_type)
  d = as.character(exp_setting$d)
  # if(response_type != "incidence"){
  #   next
  # }
  
  Xorig = datasets[[dataset_name]]$X
  X = Xorig
  if (missing_data_prop > 0){
    #add some missing data
    for (i in 2 : n){
      for (j in 1 : ncol(Xorig)){
        if (runif(1) <= missing_data_prop){
          X[i, j] = NA
        }
      }
    }    
  } else {
    X = na.omit(X)
    if (nrow(X) %% 2 == 1){
      X = X[1 : (nrow(X) - 1), ]
    }
  }
  n = nrow(X)
  
  y = datasets[[dataset_name]][["responses"]][[response_type]]
  beta_T =  if (beta_T == "zero"){
              0
            } else {
              datasets[[dataset_name]][["beta_T"]][[response_type]]
            }
  
  
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
      dead_t = as.numeric(runif(1) > prob_of_uncensored_survival_observation)
      seq_des_obj$add_subject_response(t, compute_response_with_treatment(y[t], w_t, dead_t, beta_T, response_type), dead = dead_t)
      response_added[t] = TRUE
    }
  }
  # })
  
  for (t in which(!response_added)){
    dead_t = as.numeric(runif(1) > prob_of_uncensored_survival_observation)
    suppressWarnings(seq_des_obj$add_subject_response(t, compute_response_with_treatment(y[t], w_t, dead_t, beta_T, response_type), dead = dead_t))
  }
  # cat("experiment completed?", seq_des_obj$check_experiment_completed(), "\n")
  
  if (visualizations_during_testing & grepl("KK", d)){
    stats = seq_des_obj$matching_statistics()
    cat(paste("  KK stats: prop matches:", stats$prop_subjects_matched, "# remaining in reservoir:", stats$num_subjects_remaining_in_reservoir))
    if (grepl("KK21", d))
      if (!is.null(seq_des_obj$covariate_weights)){
        ggplot_obj = ggplot(data.frame(x = names(seq_des_obj$covariate_weights), imp = seq_des_obj$covariate_weights)) +
          aes(x = x, y = imp) +
          geom_bar(stat = "identity") +
          ggtitle(paste("Dataset:", dataset_name, " Design:", d, " Response:", response_type)) +
          theme(axis.text.x = element_text(angle = 90)) +
          ylim(0, 0.6)
        print(ggplot_obj)
      }
    }
  }
  
  # let's do some inference
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
    # "incidence_KK_multivariate_logistic_regression_with_matching_dummies",	
    # "incidence_KK_multivariate_logistic_regression_with_random_intercepts_for_matches",
    ########################################### PROPORTION
    "proportion_simple_mean_difference",
    # "proportion_simple_logodds_regression",
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
    # "count_KK_multivariate_negative_binomial_regression_with_matching_dummies",
    # "count_KK_multivariate_negative_binomial_regression_with_random_intercepts_for_matches",
    ########################################### SURVIVAL
    "survival_simple_median_difference",	
    "survival_simple_restricted_mean_difference",
    "survival_univariate_weibull_regression",	
    "survival_multivariate_weibull_regression",
    # "survival_KK_compound_univariate_weibull_regression",	
    # "survival_KK_compound_multivariate_weibull_regression",
    #"survival_KK_multivariate_weibull_regression_with_matching_dummies",
    "survival_univariate_coxph_regression",	
    "survival_multivariate_coxph_regression"	
    # "survival_KK_multivariate_coxph_regression_with_matching_dummies",		
    # "survival_KK_multivariate_coxph_regression_with_random_intercepts_for_matches"				
  )){
    
    for (test_type in c("MLE-or-KM-based")){
      #only do inference by appropriate response type
      if (strsplit(estimate_type, "_")[[1]][1] != response_type){
        next
      }
      #only do inference for KK if it's a KK design
      if (grepl("KK", estimate_type) & !grepl("KK", d)){
        next
      }
      #there are two estimate types defined for all responses
      if (grepl("simple_mean_difference", estimate_type)){
        estimate_type = "simple_mean_difference"
      }
      if (grepl("KK_compound_mean_difference", estimate_type)){
        estimate_type = "KK_compound_mean_difference"
      }
      seq_des_inf_obj = SeqDesignInference$new(seq_des_obj, estimate_type = estimate_type, test_type = test_type, num_cores = num_cores)
      beta_hat_T = seq_des_inf_obj$compute_treatment_estimate()
      cat("  beta_T = ", beta_T, "  beta_hat_T =", beta_hat_T, "\n")
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


res = data.table(res)
res_summary = res[, .(
    beta_T_hat = mean(beta_hat_T, na.rm = TRUE),
    ci_a = mean(ci_a, na.rm = TRUE),
    ci_b = mean(ci_b, na.rm = TRUE),
    p_val = mean(p_val, na.rm = TRUE), 
    size = mean(p_val < 0.05, na.rm = TRUE),
    prop_NA = sum(is.na(p_val)) / .N
  ), by = c("estimate_type", "design")]
res_summary[prop_NA > 0, .(estimate_type, design, prop_NA)]
res_summary[order(size), .(estimate_type, design, size)]
