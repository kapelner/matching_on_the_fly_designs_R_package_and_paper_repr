pacman::p_load(SeqExpMatch, dplyr)

n = 50
p = 2
mu_x = 1
sigma_x = 1
sigma_e = 1
beta_T = 1
nsim_exact_test = 501
num_cores = 6
Nsim = 1000

#build mvnp covariates
set.seed(1984)
errors = rnorm(n, 0, sigma_e)

all_betas_and_correlations = list()
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0, 	betas = c(1, 1, 1, 0, 0)))) #QUAD EVEN
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0.75, 	betas = c(1, 1, 1, 0, 0)))) #QUAD MORE UNEVEN with CORR
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0, 	betas = c(6, 1, 2, 0, 0)))) #QUAD MORE UNEVEN
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0.75, 	betas = c(6, 1, 2, 0, 0)))) #QUAD MORE UNEVEN with CORR

res = data.frame(
  betas = character(),
  rho = numeric(),
  design = character(),
  estimate_type = character(),
  test_type = character(),
  beta_hat_T = numeric(),
  pval = numeric()
)
for (nsim in 1 : Nsim){
  cat ("nsim:", nsim, "/", Nsim, "\n")
  for (all_betas_and_correlation in all_betas_and_correlations){
    betas = all_betas_and_correlation[["betas"]]
    rho = all_betas_and_correlation[["rho"]] 
    
    Sigma = sigma_x * (matrix(rho, nrow = p, ncol = p) + diag(1 - rho, p))
    X = MASS::mvrnorm(n, rep(mu_x, p), Sigma)
    z = betas[1] * X[, 1] +
      betas[2] * X[, 2] + 
      betas[3] * X[, 1]^2 +
      betas[4] * X[, 2]^2 +
      betas[5] * X[, 1] * X[, 2]
    y = array(NA, n)
    
    #test all designs
    for (d in c("CRD", "BCRD", "Efron", "Atkinson", "KK14", "KK21", "KK21stepwise")){
      seq_des_obj = SeqDesign$new(n, p, d, verbose = FALSE)
      
      for (t in 1 : n){
        seq_des_obj$add_subject_to_experiment(X[t, ])
        w_t = seq_des_obj$w[seq_des_obj$t]
        y[t] = beta_T * w_t + z[t] + errors[t]
        seq_des_obj$add_current_subject_response(y[t])
      }
      
      for (test_type in c("normal-based", "randomization-exact")){
        for (estimate_type in c("difference-in-means", "OLS")){
          seq_des_inf_obj = SeqDesignInference$new(seq_des_obj, estimate_type = estimate_type, test_type = test_type, num_cores = num_cores, verbose = FALSE)
          
          beta_hat_T = seq_des_inf_obj$compute_treatment_estimate()
          pval = seq_des_inf_obj$compute_pval_for_no_treatment_effect(nsim_exact_test = nsim_exact_test)
          
          res = rbind(res, data.frame(
            betas = paste0(betas, collapse=""),
            rho = rho,
            design = d,
            estimate_type = estimate_type,
            test_type = test_type,
            beta_hat_T = beta_hat_T,
            pval = pval
          ))
        }
      }
      
    }
  }
  
}

res_mod = res %>% 
  mutate(sq_err = (beta_hat_T - beta_T)^2, rej = pval < 0.05) %>%
  group_by(betas, rho, design, estimate_type, test_type) %>%
  summarize(mse = mean(sq_err), power = sum(rej) / n())

res_mod %>% filter(betas == "11100" & rho == 0 & estimate_type == "difference-in-means" & test_type == "normal-based") %>% as.data.frame
res_mod %>% filter(betas == "11100" & rho == 0 & estimate_type == "OLS" & test_type == "normal-based") %>% as.data.frame
res_mod %>% filter(betas == "11100" & rho == 0 & estimate_type == "difference-in-means" & test_type == "randomization-exact") %>% as.data.frame
res_mod %>% filter(betas == "11100" & rho == 0 & estimate_type == "OLS" & test_type == "randomization-exact") %>% as.data.frame

res_mod %>% filter(betas == "11100" & rho == 0.75 & estimate_type == "difference-in-means" & test_type == "normal-based") %>% as.data.frame
res_mod %>% filter(betas == "11100" & rho == 0.75 & estimate_type == "OLS" & test_type == "normal-based") %>% as.data.frame
res_mod %>% filter(betas == "11100" & rho == 0.75 & estimate_type == "difference-in-means" & test_type == "randomization-exact") %>% as.data.frame
res_mod %>% filter(betas == "11100" & rho == 0.75 & estimate_type == "OLS" & test_type == "randomization-exact") %>% as.data.frame

res_mod %>% filter(betas == "61200" & rho == 0 & estimate_type == "difference-in-means" & test_type == "normal-based") %>% as.data.frame
res_mod %>% filter(betas == "61200" & rho == 0 & estimate_type == "OLS" & test_type == "normal-based") %>% as.data.frame
res_mod %>% filter(betas == "61200" & rho == 0 & estimate_type == "difference-in-means" & test_type == "randomization-exact") %>% as.data.frame
res_mod %>% filter(betas == "61200" & rho == 0 & estimate_type == "OLS" & test_type == "randomization-exact") %>% as.data.frame

res_mod %>% filter(betas == "61200" & rho == 0.75 & estimate_type == "difference-in-means" & test_type == "normal-based") %>% as.data.frame
res_mod %>% filter(betas == "61200" & rho == 0.75 & estimate_type == "OLS" & test_type == "normal-based") %>% as.data.frame
res_mod %>% filter(betas == "61200" & rho == 0.75 & estimate_type == "difference-in-means" & test_type == "randomization-exact") %>% as.data.frame
res_mod %>% filter(betas == "61200" & rho == 0.75 & estimate_type == "OLS" & test_type == "randomization-exact") %>% as.data.frame







