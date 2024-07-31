#/c/Program\ Files/R/R-devel/bin/R.exe CMD INSTALL -l ~/AppData/Local/R/win-library/4.4/ SeqExpMatch/
pacman::p_load(PTE, dplyr, profvis)
library(SeqExpMatch)

options(error=recover)
set.seed(1984)

data(continuous_example)
X = as_tibble(continuous_example$X) %>% 
  select(-treatment) %>%
  mutate(x1 = if_else(x1 == "YES", 1, 0)) %>%
  mutate(x2 = if_else(x2 == "ON", 1, 0)) %>%
  mutate(x3 = if_else(x3 == "YES", 1, 0)) %>%
  mutate(x4_high = if_else(x4 == "HIGH", 1, 0)) %>%
  mutate(x4_medium = if_else(x4 == "MEDIUM", 1, 0)) %>%
  select(-x4) %>%
  as.matrix

n = nrow(X)
p = ncol(X)
nsim_exact_test = 101

#test all response_types and all designs
for (response_type in c("continuous", "incidence", "proportion", "count", "survival")){
  for (d in c("CRD", "KK21")){ #"CRD", "iBCRD", "Efron", "Atkinson", "KK14", "KK21", "KK21stepwise"
    cat("\n\n\n\n\nresponse_type =", response_type, "d =", d, "\n\n\n\n\n")
    
    seq_des_obj = SeqDesign$new(n, p, d, response_type = response_type, verbose = FALSE)
    
    y_cont = continuous_example$y 
    y_cont_pos = y_cont - min(y_cont) + .Machine$double.eps
    y = switch(response_type,
          continuous = y_cont,
          incidence = ifelse(y_cont > median(y_cont), 1, 0),
          proportion = y_cont_pos / (max(y_cont_pos) + 2 * .Machine$double.eps),
          count = round(y_cont_pos),
          survival = y_cont_pos
    )
    dead = switch(response_type,
            continuous = rep(1, n),
            incidence = rep(1, n),
            proportion = rep(1, n),
            count = rep(1, n),
            survival = rbinom(n, 1, prob = 0.5)
    )
      
    for (t in 1 : n){
      seq_des_obj$add_subject_to_experiment(X[t, ])
      #seq_des_obj$print_current_subject_assignment()
      if (d %in% c("KK21", "KK21stepwise")){
        seq_des_obj$add_subject_response(t, y[t], dead[t])
      }
    }    
  
    if (!(d %in% c("KK21", "KK21stepwise"))){
      seq_des_obj$add_all_subject_responses(y, dead)
    }
    # cat("experiment completed?", seq_des_obj$check_experiment_completed(), "\n")
  
    # if (grepl("KK", d)){
    #   print(seq_des_obj$matching_statistics())
    # }
    
    # let's do some inference
    for (test_type in c("MLE-or-KM-based", "randomization-exact")){ #
      for (estimate_type in c("mean_difference", "median_difference", "default_regression")){
        if (estimate_type == "median_difference" & response_type != "survival"){
          next
        }
        seq_des_inf_obj = SeqDesignInference$new(seq_des_obj, estimate_type = estimate_type, test_type = test_type, num_cores = 6)
        cat("  beta_hat_T =", seq_des_inf_obj$compute_treatment_estimate(), "\n")
        cat("  95% CI for betaT = [", paste(seq_des_inf_obj$compute_confidence_interval(nsim_exact_test = nsim_exact_test), collapse = ", "), "]\n")
        # cat("  pval =", seq_des_inf_obj$compute_two_sided_pval_for_treatment_effect(nsim_exact_test = nsim_exact_test), "\n")
      }
    }
  }
}

