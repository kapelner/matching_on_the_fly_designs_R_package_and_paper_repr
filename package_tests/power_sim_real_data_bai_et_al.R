pacman::p_load(SeqExpMatch, dplyr, data.table, datasets)
pacman::p_load(doParallel, PTE, datasets, qgam, mlbench, AppliedPredictiveModeling, dplyr, ggplot2, gridExtra, profvis, data.table, profvis)
rm(list = ls())
options(error = recover)

#n = 50
nsim_exact_test = 501
num_cores = 20
Nsim = 5000

#build mvnp covariates
set.seed(1984)

max_n_dataset = 25000

source("_dataset_load.R")

for (dataset_name in names(datasets_and_response_models)){
  if (!("continuous" %in% names(datasets_and_response_models[[dataset_name]]$beta_T)) || (dataset_name %in% c("pte_example", "iris"))){
    datasets_and_response_models[[dataset_name]] = NULL
  } else {
    datasets_and_response_models[[dataset_name]]$response_models = list()
    response_type = "continuous"
    y_orig = datasets_and_response_models[[dataset_name]]$y_original[[response_type]]
    datasets_and_response_models[[dataset_name]]$response_models[[response_type]] = lm(y_orig ~ ., datasets_and_response_models[[dataset_name]]$X)
  }
}
rm(y_orig, dataset_name, response_type)


res = data.frame(
  data_set = character(),
  beta_T = numeric(),
  n = numeric(),
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
    for (dataset_name in names(datasets_and_response_models)){
      mod = datasets_and_response_models[[dataset_name]]
      for (beta_T in c(0, mod$beta_T$continuous)){
        for (n in c(50, 100, 500, 1000)){
          
          Xy = data.table(mod$X, y = mod$y_original$continuous)
          Xy_trim_inx = sample(1:nrow(Xy), n)
          Xy_trim = Xy[Xy_trim_inx, ]
          y = Xy_trim[,y.V1]
          X = Xy_trim[,y.V1 := NULL]
          y_filled = rep(NA, n)
          p = ncol(X)
    
          #test all designs
          for (d in c("KK14", "KK21", "KK21stepwise")){
            
            des = paste0("SeqDesign", d)
            des_class = get(des)
            seq_des_obj_morrison_n = des_class$new(n = n, prob_T = 0.5, p = p, response_type = "continuous", verbose = FALSE, morrison = TRUE)
            seq_des_obj_morrison = des_class$new(prob_T = 0.5, p = p, response_type = "continuous", verbose = FALSE, morrison = TRUE)
            seq_des_obj = des_class$new(n = n, prob_T = 0.5, p = p, response_type = "continuous", verbose = FALSE, morrison = FALSE)
            
            for(mor in c("_morrison_n", "_morrison", "")){
              
              cur_seq_des_obj = get(paste0("seq_des_obj", mor))
              
              for (t in 1 : n){
                cur_seq_des_obj$add_subject_to_experiment_and_assign(X[t, ])
                w_t = cur_seq_des_obj$get_w()[cur_seq_des_obj$get_t()]
                y_filled[t] = beta_T * w_t + y[t]
                cur_seq_des_obj$add_subject_response(t = t, y = as.numeric(y_filled[t]))
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
                      data_set = dataset_name,
                      beta_T = beta_T,
                      n = n,
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
                    data_set = dataset_name,
                    beta_T = beta_T,
                    n = n,
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