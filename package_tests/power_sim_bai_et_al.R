pacman::p_load(SeqExpMatch, dplyr, data.table)
rm(list = ls())
options(error = recover)

#n = 50
p = 2
mu_x = 1
sigma_x = 1
sigma_e = 1
nsim_exact_test = 501
num_cores = 20
Nsim = 5000

#build mvnp covariates
set.seed(1984)


all_betas_and_correlations = list()
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0, 	betas = c(1, 1, 1, 0, 0)))) #QUAD EVEN
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0.75, 	betas = c(1, 1, 1, 0, 0)))) #QUAD MORE UNEVEN with CORR
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0, 	betas = c(6, 1, 2, 0, 0)))) #QUAD MORE UNEVEN
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0.75, 	betas = c(6, 1, 2, 0, 0)))) #QUAD MORE UNEVEN with CORR
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0, 	betas = c(4, 1, 2, 3, 2)))) #QUAD MORE UNEVEN
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0.75, 	betas = c(4, 1, 2, 3, 2)))) #QUAD MORE UNEVEN with CORR
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0, 	betas = c(4, 1, 2, 3, 2)))) #QUAD MORE UNEVEN
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0.75, 	betas = c(4, 1, 2, 3, 2)))) #QUAD MORE UNEVEN with CORR

res = data.frame(
  beta_T = numeric(),
  n = numeric(),
  betas = character(),
  rho = numeric(),
  design = character(),
  mor = character(),
  infrence = character(),
  reservoir = numeric(),
  beta_hat_T = numeric(),
  pval = numeric()
)
times = rep(NA, Nsim)
for (nsim in 1 : Nsim){
  start_time = Sys.time()
  cat ("nsim:", nsim, "/", Nsim, "\n")
  
  tryCatch({
    for (beta_T in c(0,1)){
      for (n in c(50, 100, 500, 1000)){
        errors = rnorm(n, 0, sigma_e)
        for (all_betas_and_correlation in all_betas_and_correlations){
          betas = all_betas_and_correlation[["betas"]]
          rho = all_betas_and_correlation[["rho"]] 
          
          Sigma = sigma_x * (matrix(rho, nrow = p, ncol = p) + diag(1 - rho, p))
          X = data.table(MASS::mvrnorm(n, rep(mu_x, p), Sigma))
          z = betas[1] * X[, 1] +
            betas[2] * X[, 2] + 
            betas[3] * X[, 1]^2 +
            betas[4] * X[, 2]^2 +
            betas[5] * X[, 1] * X[, 2]
          y = array(NA, n)
          #test all designs
          for (d in c("KK14", "KK21", "KK21stepwise")){
            des = paste0("SeqDesign", d)
            des_class = get(des)
            seq_des_obj_morrison_n = des_class$new(n = n, prob_T = 0.5, p = p, response_type = "continuous", verbose = FALSE, morrison = TRUE)
            seq_des_obj_morrison = des_class$new(prob_T = 0.5, p = p, response_type = "continuous", verbose = FALSE, morrison = TRUE)
            seq_des_obj = des_class$new(n = n, prob_T = 0.5, p = p, response_type = "continuous", verbose = FALSE, morrison = FALSE)
            for(mor in c("_morrison_n", "_morrison", "")){
              #cat(d, mor, '\n')
              
              cat(d,mor,'\n')
              
              cur_seq_des_obj = get(paste0("seq_des_obj", mor))
              
              for (t in 1 : n){
                cur_seq_des_obj$add_subject_to_experiment_and_assign(X[t, ])
                w_t = cur_seq_des_obj$get_w()[cur_seq_des_obj$get_t()]
                y[t] = beta_T * w_t + z[t] + errors[t]
                cur_seq_des_obj$add_subject_response(t = t, y = as.numeric(y[t]))
              }

              for(infrence in c("AllKKCompoundMeanDiff", "BaiAdjustedT")){
                des = paste0("SeqDesignInference", infrence)
                if(infrence == "BaiAdjustedT"){
                  des = paste0(des, d)
                  des_class = get(des)
                  for(conex in c(TRUE, FALSE)){
                    seq_des_inf_obj = des_class$new(cur_seq_des_obj, num_cores = num_cores, verbose = FALSE, convex = conex)
                    if(conex == TRUE){
                      infrence_specific = paste0(infrence, "_convex")
                    } else {
                      infrence_specific = infrence
                    }
                    beta_hat_T = seq_des_inf_obj$compute_treatment_estimate()
                    pval = seq_des_inf_obj$compute_mle_two_sided_pval_for_treatment_effect()
                    
                    res = rbind(res, data.frame(
                      betas = paste0(betas, collapse=""),
                      beta_T = beta_T,
                      n = n,
                      rho = rho,
                      design = d,
                      mor = mor,
                      infrence = infrence_specific,
                      reservoir = sum(seq_des_obj$get_match_indic() == 0),
                      beta_hat_T = beta_hat_T,
                      pval = pval
                    ))
                  }
                } else {
                  des_class = get(des)
                  seq_des_inf_obj = des_class$new(cur_seq_des_obj, num_cores = num_cores, verbose = FALSE)
                  
                  beta_hat_T = seq_des_inf_obj$compute_treatment_estimate()
                  pval = seq_des_inf_obj$compute_mle_two_sided_pval_for_treatment_effect()
                  
                  res = rbind(res, data.frame(
                    betas = paste0(betas, collapse=""),
                    beta_T = beta_T,
                    n = n,
                    rho = rho,
                    design = d,
                    mor = mor,
                    infrence = infrence,
                    reservoir = sum(seq_des_obj$get_match_indic() == 0),
                    beta_hat_T = beta_hat_T,
                    pval = pval
                  ))
                }
              }
            }
          }
        }
      }
    }
    times[nsim] = Sys.time() - start_time
    print(Sys.time() - start_time)
  }, error = function(e) {
    message(paste("⚠️ Error in nsim =", nsim, ":", conditionMessage(e)))
    times[nsim] = NA  
  })
}
res_mod = res %>%
  mutate(sq_err = (beta_hat_T - beta_T)^2, rej = pval < 0.05) %>%
  group_by(betas, rho, beta_T, n, design, mor, inference) %>%
  summarize(reservoir_size = mean(reservoir), mse = mean(sq_err), percent_reject = sum(rej) / n())