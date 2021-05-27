#R CMD INSTALL -l ~/Documents/R/win-library/3.6/ SeqExpMatch
options(error=recover)
pacman::p_load(SeqExpMatch, PTE, dplyr, profvis)
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
y = continuous_example$y

#test all designs
for (d in c("CRD", "BCRD", "Efron", "Atkinson", "KK14", "KK21", "KK21stepwise")){
# for (d in c("KK14", "KK21", "KK21stepwise")){
# for (d in c("KK21stepwise")){
  seq_des_obj = SeqDesign$new(n, p, d)
  
  for (i in 1 : n){
    seq_des_obj$add_subject_to_experiment(X[i, ])
    seq_des_obj$print_current_subject_assignment()
    if (d %in% c("KK21", "KK21stepwise")){
      seq_des_obj$add_current_subject_response(y[i])
    }
  }    

  if (!(d %in% c("KK21", "KK21stepwise"))){
    seq_des_obj$add_all_subject_responses(y)
  }
  # cat("experiment completed?", seq_des_obj$check_experiment_completed(), "\n")

  # if (d %in% c("KK14", "KK21", "KK21stepwise")){
  #   print(seq_des_obj$matching_statistics())
  # }
  
  #let's do some inference
  
  # seq_des_inf_obj = SeqDesignInference$new(seq_des_obj, test_type = "normal-based", estimate_type = "difference-in-means", num_cores = 10)
  # cat("  beta_hat_T =", seq_des_inf_obj$compute_treatment_estimate(), "\n")
  # cat("  pval =", seq_des_inf_obj$compute_pval_for_no_treatment_effect(), "\n")
  # cat("  95% CI for betaT = [", paste(seq_des_inf_obj$compute_confidence_interval(), collapse = ", "), "]\n")
  # seq_des_inf_obj = SeqDesignInference$new(seq_des_obj, test_type = "normal-based", estimate_type = "OLS", num_cores = 1)
  # cat("  beta_hat_T =", seq_des_inf_obj$compute_treatment_estimate(), "\n")
  # cat("  pval =", seq_des_inf_obj$compute_pval_for_no_treatment_effect(), "\n")
  # cat("  95% CI for betaT = [", paste(seq_des_inf_obj$compute_confidence_interval(), collapse = ", "), "]\n")
  # seq_des_inf_obj = SeqDesignInference$new(seq_des_obj, test_type = "randomization-exact", estimate_type = "difference-in-means", num_cores = 6)
  # cat("  beta_hat_T =", seq_des_inf_obj$compute_treatment_estimate(), "\n")
  # # profvis({seq_des_inf_obj$compute_pval_for_no_treatment_effect(nsim_exact_test = 5)})
  # cat("  pval =", seq_des_inf_obj$compute_pval_for_no_treatment_effect(nsim_exact_test = 501), "\n")
  # seq_des_inf_obj = SeqDesignInference$new(seq_des_obj, test_type = "randomization-exact", estimate_type = "OLS", num_cores = 6)
  # cat("  beta_hat_T =", seq_des_inf_obj$compute_treatment_estimate(), "\n")
  # # profvis({seq_des_inf_obj$compute_pval_for_no_treatment_effect(nsim_exact_test = 50)})
  # cat("  pval =", seq_des_inf_obj$compute_pval_for_no_treatment_effect(nsim_exact_test = 501), "\n")
  
  
  
  for (test_type in c("normal-based", "randomization-exact")){
    for (estimate_type in c("difference-in-means", "OLS")){
      seq_des_inf_obj = SeqDesignInference$new(seq_des_obj, estimate_type = estimate_type, test_type = test_type, num_cores = 6)
      beta_hat_T = seq_des_inf_obj$compute_treatment_estimate()
      cat("  beta_hat_T =", beta_hat_T, "\n")
      pval = seq_des_inf_obj$compute_pval_for_no_treatment_effect(nsim_exact_test = 50)
      cat("  pval =", pval, "\n")

      if (test_type == "normal-based"){
        cat("  95% CI for betaT = [", paste(seq_des_inf_obj$compute_confidence_interval(), collapse = ", "), "]\n")
      }

    }
  }

}

# b_T_sims = seq_des_inf_obj$randomization_inference_samples_for_no_treatment_effect(nsim_exact_test = 50)
