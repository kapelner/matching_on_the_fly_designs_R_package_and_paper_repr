pacman::p_load(SeqExpMatch, dplyr, data.table, doFuture, future, doRNG, foreach, progressr) #doParallel
pacman::p_load(doParallel, PTE, datasets, qgam, mlbench, AppliedPredictiveModeling, dplyr, ggplot2, gridExtra, profvis, data.table, profvis)
rm(list = ls())
set.seed(1986)
options(error = recover)


nsim_exact_test = 501
num_cores = availableCores()-10
Nsim = 1

max_n_dataset = 5000

source("_dataset_load.R")

for (dataset_name in names(datasets_and_response_models)){
  if (!("continuous" %in% names(datasets_and_response_models[[dataset_name]]$beta_T)) || (dataset_name %in% c("pte_example", "iris"))){
    datasets_and_response_models[[dataset_name]] = NULL
  } else {
    datasets_and_response_models[[dataset_name]]$response_models = list()
    response_type = "continuous"
    y_orig = datasets_and_response_models[[dataset_name]]$y_original[[response_type]]
    datasets_and_response_models[[dataset_name]]$response_models[[response_type]] = lm(y_orig ~ ., datasets_and_response_models[[dataset_name]]$X)
    for (type in setdiff(names(datasets_and_response_models[[dataset_name]]$y_original), 'continuous')){
      datasets_and_response_models[[dataset_name]]$y_original[[type]] = NULL
    }
  }
}
rm(y_orig, dataset_name, response_type)

data_sets = names(datasets_and_response_models)
designs = c("KK14", "KK21", "KK21stepwise")
morisons = c("_morrison_n", "_morrison", "default")
ns = c(50, 100, 500)
beta_Ts = c(0, 0.2)


params = expand.grid(
  i = 1:Nsim,
  data_set = data_sets,
  beta_T = beta_Ts,
  n = ns,
  design = designs,
  morison = morisons,
  KEEP.OUT.ATTRS = FALSE
)
params = params %>%
  arrange(i, data_set, beta_T, n, design, morison)

run_simulation = function(i, data_set, beta_T, design, morison, n){
  res = data.frame(
    i = numeric(),
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
  
  mod = datasets_and_response_models[[data_set]]
  
  Xy = data.table(mod$X, y = mod$y_original$continuous)
  Xy_trim_inx = sample(1:nrow(Xy), n)
  Xy_trim = Xy[Xy_trim_inx, ]
  y = Xy_trim[,y.V1]
  X = Xy_trim[,y.V1 := NULL]
  y_filled = rep(NA, n)
  p = ncol(X)
  
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
    y_filled[t] = beta_T * w_t + y[t]
    seq_des_obj$add_subject_response(t = t, y = as.numeric(y_filled[t]))
  }
  
  #cat("nsim:", i, ", 2\n")
  for(infrence in c("AllKKCompoundMeanDiff", "BaiAdjustedT")){ 
    #cat(infrence, "\n")
    des_inf = paste0("SeqDesignInference", infrence)
    if(infrence == "BaiAdjustedT"){
      des_inf = paste0(des_inf, design)
      des_inf_class = get(des_inf)
      for(conex in c(TRUE, FALSE)){
        seq_des_inf_obj = des_inf_class$new(seq_des_obj, verbose = FALSE, convex = conex) 
        if(conex == TRUE){
          infrence_specific = paste0(infrence, "_convex")
        } else {
          infrence_specific = infrence
        }
        
        beta_hat_T = seq_des_inf_obj$compute_treatment_estimate()
        pval = seq_des_inf_obj$compute_mle_two_sided_pval_for_treatment_effect()
        
        res = rbind(res, data.frame(
          i = i,
          data_set = data_set,
          beta_T = beta_T,
          n = n,
          design = design,
          mor = morison,
          infrence = infrence_specific,
          reservoir = sum(seq_des_obj$get_match_indic() == 0),
          beta_hat_T = beta_hat_T,
          pval = pval
        ))
      }
    } else {
      des_inf_class = get(des_inf)
      seq_des_inf_obj = des_inf_class$new(seq_des_obj, verbose = FALSE)
      beta_hat_T = seq_des_inf_obj$compute_treatment_estimate()
      pval = seq_des_inf_obj$compute_mle_two_sided_pval_for_treatment_effect()
      res = rbind(res, data.frame(
        i = i,
        data_set = data_set,
        beta_T = beta_T,
        n = n,
        design = design,
        mor = morison,
        infrence = infrence,
        reservoir = sum(seq_des_obj$get_match_indic() == 0),
        beta_hat_T = beta_hat_T,
        pval = pval
      ))
    }
  }
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
    data_set = row$data_set
    beta_T = row$beta_T
    n = row$n
    design = row$design
    morison = row$morison
    cat(glue::glue("Running i={i}, n={n}, design={design}, morison={morison}"), '\n')
    res = tryCatch({
      out = run_simulation(i, data_set, beta_T, design, morison, n)
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