rm(list = ls())
setwd(utils::getSrcDirectory(function(){})[1])
source("performance_study_setup.R")
options(error=recover)


nC = 5
plan(callr, workers = nC, gc = TRUE)
# cl = makeCluster(nC)
# registerDoSEQ()
# progress=NULL

# cl <- makeSOCKcluster(nC)
# registerDoSNOW(cl)
# pb <- tkProgressBar(max=nrow(exp_settings))
# progress <- function(i) setTkProgressBar(pb, i, label =
#                           paste0(i, " of ", nrow(exp_settings), " i.e., ", round(i / nrow(exp_settings) * 100, 1), "%"))

# clusterExport(cl, list(
#   "exp_settings", 
#   "X100", 
#   "linear_mean_function", 
#   "response_functions", 
#   "estimands_betaT_one", 
#   "betas", 
#   "betaToverall", 
#   "censoring_mechanism",
#   "survival_mu_multiple", 
#   "survival_k", 
#   "mu_survival_max_to_be_observed",
#   "nsim_exact_test"
# ), envir = environment())






# exp_settings = exp_settings[1:1000,]

# profvis({
# for (n_setting in 1 : nrow(exp_settings)){
# res = foreach(
#     n_setting = 1 : nrow(exp_settings),
#     .inorder = FALSE,
#     .combine = rbind,
#     .packages = 'SeqExpMatch',
#     .options.snow = list(progress = progress)) %dopar% {
run_all_simulations = function(n_settings) {
  progress_bar = progressr::progressor(along = n_settings)
  total_number_of_iterations = max(n_settings)
  future_lapply(n_settings, 
                future.seed = TRUE, 
                future.packages = 'SeqExpMatch',
                future.globals = c(
                  "exp_settings",
                  "X100",
                  "linear_mean_function",
                  "response_functions",
                  "estimands_betaT_one",
                  "betas",
                  "betaToverall",
                  "censoring_mechanism",
                  "survival_mu_multiple",
                  "survival_k",
                  "mu_survival_max_to_be_observed",
                  "nsim_exact_test"
                ),
                FUN = function(n_setting) {
    #get this run's settings
    n =                exp_settings$n[n_setting]
    response_type =    as.character(exp_settings$response_type[n_setting])
    design =           as.character(exp_settings$design[n_setting])
    test_type =        as.character(exp_settings$test_type[n_setting])
    inference_method = as.character(exp_settings$inference_method[n_setting])
    betaT =            exp_settings$betaT[n_setting]
    prob_of_adding_response = exp_settings$prob_of_adding_response[n_setting]
    
    progress_bar(sprintf(paste0("n_setting = %g/", total_number_of_iterations), n_setting),
      message = paste0(
        "#", stringr::str_pad(n_setting, 6, side = "left", "0"), 
        "\n  response_type: ", response_type, 
        "\n  design:", design,
        "\n  test_type:", test_type,
        "\n  inference_method:", inference_method,
        "\n"
      ))
    
    response_function = response_functions[[response_type]]
    dead_function =     response_functions[["dead"]]
    

    seq_des_obj = SeqDesign$new(n, design, response_type, verbose = FALSE)
    
    #now we add the subjects and sometimes assign y_t values
    response_added = array(FALSE, n)
    for (t in 1 : n){
      x_t = X100[t, ]
      w_t = seq_des_obj$add_subject_to_experiment_and_assign(x_t)
      if (runif(1) < prob_of_adding_response){
        y_t = response_function(as.numeric(x_t), w_t, betaT)
        dead_t = ifelse(response_type == "survival", response_functions[["dead"]](y_t), 1)
        seq_des_obj$add_subject_response(t, y_t, dead_t)
        response_added[t] = TRUE
      }
    }
    
    for (t in which(!response_added)){
      y_t = response_function(as.numeric(X100[t, ]), seq_des_obj$w[t], betaT)
      dead_t = ifelse(response_type == "survival", response_functions[["dead"]](y_t), 1)
      seq_des_obj$add_subject_response(t, y_t, dead_t)
    }
    
    seq_des_inf_obj = SeqDesignInference$new(seq_des_obj, estimate_type = inference_method, test_type = test_type, verbose = FALSE)
  
    estimand = ifelse(betaT == 0, 0, estimands_betaT_one[[inference_method]])
    estimate = seq_des_inf_obj$compute_treatment_estimate()
    if (test_type == "MLE-or-KM-based"){
      ci = seq_des_inf_obj$compute_confidence_interval(nsim_exact_test = nsim_exact_test)
      ci_a = ci[1]
      ci_b = ci[2]
    } else {
      ci_a = NA
      ci_b = NA
    }
    p_val = seq_des_inf_obj$compute_two_sided_pval_for_treatment_effect(nsim_exact_test = nsim_exact_test)
    
    
    #if it's sequential
    if (n_setting %% 37 == 0){
      cat("run", n_setting, "of", nrow(exp_settings), "\n")
    }
    
    res_row = list(
      idx =              n_setting,
      nsim =             exp_settings$nsims[n_setting],
      n =                n,
      response_type =    response_type, 
      design =           design, 
      test_type =        test_type, 
      inference_method = inference_method,
      prob_of_adding_response = prob_of_adding_response,
      betaT =            betaT,
      estimand =         estimand, 
      estimate =         estimate,
      ci_a =             ci_a, 
      ci_b =             ci_b,
      p_val =            p_val,
      weights =          paste(seq_des_obj$covariate_weights, collapse = ",")
    )
    
    #much faster
    # data.table::set(res, n_setting, names(res), res_row)
    # save(res_row, file = paste0("res/result_", stringr::str_pad(n_setting, side = "left", 7, pad = "0")))
    res_row
  })
}
res = run_all_simulations(1 : nrow(exp_settings))
# stopCluster(cl)
# rm(cl); gc()

filename = paste0("res/perf", ".RData")
save(res, file = filename)