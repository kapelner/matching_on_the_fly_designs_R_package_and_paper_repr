pacman::p_load(SeqExpMatch, dplyr, data.table, doFuture, future, doRNG, foreach, progressr) #doParallel
rm(list = ls())
set.seed(1986)
options(error = recover)

p = 2
mu_x = 1
sigma_x = 1
sigma_e = 1
nsim_exact_test = 501
num_cores = availableCores()-1
Nsim = 10

beta_Ts = c(0, 1)
configs = list()
configs = c(configs, list(setting = list(rho = 0, 	betas = c(1, 1, 1, 0, 0)))) #QUAD EVEN
configs = c(configs, list(setting = list(rho = 0.75, 	betas = c(1, 1, 1, 0, 0)))) #QUAD MORE UNEVEN with CORR
configs = c(configs, list(setting = list(rho = 0, 	betas = c(6, 1, 2, 0, 0)))) #QUAD MORE UNEVEN
configs = c(configs, list(setting = list(rho = 0.75, 	betas = c(6, 1, 2, 0, 0)))) #QUAD MORE UNEVEN with CORR
configs = c(configs, list(setting = list(rho = 0, 	betas = c(4, 1, 2, 3, 2)))) #QUAD MORE UNEVEN
configs = c(configs, list(setting = list(rho = 0.75, 	betas = c(4, 1, 2, 3, 2)))) #QUAD MORE UNEVEN with CORR
designs = c("KK14", "KK21", "KK21stepwise")
morisons = c("_morrison_n", "_morrison", "default")
ns = c(50, 100, 500, 1000)

params = expand.grid(
  i = 1:Nsim,
  beta_T = beta_Ts,
  config = configs,
  n = ns,
  design = designs,
  morison = morisons,
  KEEP.OUT.ATTRS = FALSE
)
params = params %>%
  arrange(i, beta_T, config, n, design, morison)

run_simulation = function(i, beta_T, config, design, morison, n){
  #cat("nsim:", i, ", 1\n")
  res = data.frame(
    i = numeric(),
    beta_T = numeric(),
    n = numeric(),
    betas = character(),
    rho = numeric(),
    design = character(),
    mor = character(),
    inference = character(),
    reservoir = numeric(),
    beta_hat_T = numeric(),
    pval = numeric()
  )
  errors = rnorm(n, 0, sigma_e)
  betas = config$setting$betas
  rho = config$setting$rho
  Sigma = sigma_x * (matrix(rho, nrow = p, ncol = p) + diag(1 - rho, p))
  X = data.table(MASS::mvrnorm(n, rep(mu_x, p), Sigma))
  z = betas[1] * X[, 1] +
    betas[2] * X[, 2] + 
    betas[3] * X[, 1]^2 +
    betas[4] * X[, 2]^2 +
    betas[5] * X[, 1] * X[, 2]
  y = array(NA, n)
  
  
  des = paste0("SeqDesign", design)
  des_class = get(des)
  
  
  if (morison == "_morrison_n") {
    seq_des_obj = des_class$new(n = n, prob_T = 0.5, p = p,
                                 response_type = "continuous", verbose = FALSE, morrison = TRUE)
  } else if (morison == "_morrison") {
    seq_des_obj = des_class$new(prob_T = 0.5, p = p,
                                 response_type = "continuous", verbose = FALSE, morrison = TRUE)
  } else {
    seq_des_obj = des_class$new(n = n, prob_T = 0.5, p = p,
                                 response_type = "continuous", verbose = FALSE, morrison = FALSE)
  }
  
  for (t in 1 : n){
    seq_des_obj$add_subject_to_experiment_and_assign(X[t, ])
    w_t = seq_des_obj$get_w()[seq_des_obj$get_t()]
    y[t] = beta_T * w_t + z[t] + errors[t]
    seq_des_obj$add_subject_response(t = t, y = as.numeric(y[t]))
  }
  rm(errors, Sigma, X, z, y, t, des, des_class)
  gc()
  #cat("nsim:", i, ", 2\n")
  for(inference in c("AllKKCompoundMeanDiff", "BaiAdjustedT")){ 
    #cat(inference, "\n")
    des_inf = paste0("SeqDesignInference", inference)
    if(inference == "BaiAdjustedT"){
      des_inf = paste0(des_inf, design)
      des_inf_class = get(des_inf)
      for(conex in c(TRUE, FALSE)){
        #for (bai_v in c(TRUE, FALSE)){
          seq_des_inf_obj = des_inf_class$new(seq_des_obj, verbose = FALSE, convex = conex) #, var_squ_flag = var_sqr
          if(conex == TRUE){
            inference_specific = paste0(inference, "_convex")
          } else {
            inference_specific = inference
          }
          beta_hat_T = seq_des_inf_obj$compute_treatment_estimate()
          pval = seq_des_inf_obj$compute_mle_two_sided_pval_for_treatment_effect()
          
          res = rbind(res, data.frame(
            betas = paste0(betas, collapse=""),
            i = i,
            beta_T = beta_T,
            n = n,
            rho = rho,
            design = design,
            mor = morison,
            inference = inference_specific,
            reservoir = sum(seq_des_obj$get_match_indic() == 0),
            beta_hat_T = beta_hat_T,
            pval = pval
          ))
        #}
      }
    } else {
      des_inf_class = get(des_inf)
      seq_des_inf_obj = des_inf_class$new(seq_des_obj, verbose = FALSE)
      beta_hat_T = seq_des_inf_obj$compute_treatment_estimate()
      pval = seq_des_inf_obj$compute_mle_two_sided_pval_for_treatment_effect()
      res = rbind(res, data.frame(
        betas = paste0(betas, collapse=""),
        i = i,
        beta_T = beta_T,
        n = n,
        rho = rho,
        design = design,
        mor = morison,
        inference = inference,
        reservoir = sum(seq_des_obj$get_match_indic() == 0),
        beta_hat_T = beta_hat_T,
        pval = pval
      ))
    }
  }
  #cat("nsim:", i, ", 3\n")
  return(res)
}


handlers(global = TRUE)
handlers("txtprogressbar")


registerDoFuture()
plan(multisession, workers = num_cores)

start_time = Sys.time()

with_progress({
  prog = progressor(along = 1:nrow(params))
  
  results = foreach(row = iter(params, by = "row"), .combine = rbind, .packages = c("SeqExpMatch", "data.table", "dplyr")) %dorng% {
    
    i = row$i
    beta_T = row$beta_T
    config = row$config
    n = row$n
    design = row$design
    morison = row$morison
    cat(glue::glue("Running i={i}, n={n}, design={design}, morison={morison}"), '\n')
    res = tryCatch({
      out = run_simulation(i, beta_T, config, design, morison, n)
      cat("Successfully ran simulation")
      prog()
      out
    }, error = function(e) {
      cat(glue::glue("Error in sim i={i}, n={n}, design={design}, mor={morison}: {e$message}"), '\n')
      prog()  # still update progress bar even if it fails
      NULL    # return NULL if failed, will be dropped in rbind
    })
  }
})


plan(sequential)
end_time = Sys.time()


res_mod = results %>%
  mutate(sq_err = (beta_hat_T - beta_T)^2, rej = pval < 0.05) %>%
  group_by(betas, rho, beta_T, n, design, mor, inference) %>%
  summarize(reservoir_size = mean(reservoir), mse = mean(sq_err), percent_reject = sum(rej) / n())