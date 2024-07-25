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

#test all response_types and all designs
for (response_type in c("continuous", "incidence", "proportion", "count_uncensored", "count_censored", "survival_uncensored", "survival_censored")){
  for (d in c("CRD")){ #, "iBCRD", "Efron", "Atkinson", "KK14", "KK21", "KK21stepwise"
    cat("\n\n\n\n\nresponse_type =", response_type, "d =", d, "\n\n\n\n\n")
    
    seq_des_obj = SeqDesign$new(n, p, d, response_type = response_type, verbose = FALSE)
    
    y_cont = continuous_example$y 
    y_cont_pos = y_cont - min(y_cont) + .Machine$double.eps
    y = switch(response_type,
          continuous = y_cont,
          incidence = ifelse(y_cont > median(y_cont), 1, 0),
          proportion = y_cont_pos / (max(y_cont_pos) + 2 * .Machine$double.eps),
          count_uncensored = round(y_cont_pos),
          count_censored = round(y_cont_pos),
          survival_uncensored = y_cont_pos,
          survival_censored = y_cont_pos
    )
    dead = switch(response_type,
            continuous = rep(1, n),
            incidence = rep(1, n),
            proportion = rep(1, n),
            count_uncensored = rep(1, n),
            count_censored = rbinom(n, 1, prob = 0.7),
            survival_uncensored = rep(1, n),
            survival_censored = rbinom(n, 1, prob = 0.7)
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
  
    # if (d %in% c("KK14", "KK21", "KK21stepwise")){
    #   print(seq_des_obj$matching_statistics())
    # }
    
    #let's do some inference
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
}

