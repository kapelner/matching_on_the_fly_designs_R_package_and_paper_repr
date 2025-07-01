#/c/Program\ Files/R/R-devel/bin/R.exe CMD INSTALL -l ~/AppData/Local/R/win-library/4.5/ SeqExpMatch/
rm(list = ls())
pacman::p_load(doParallel, PTE, datasets, qgam, mlbench, AppliedPredictiveModeling, dplyr, ggplot2, gridExtra, profvis, data.table, profvis)
library(SeqExpMatch)
# options(error=recover)
# options(warn=2)
set.seed(1)

#options(future.debug = TRUE)
# registerDoFuture()
nC = 18 #num cores
# plan(sequential)

nsim_exact_test = 301
missing_data_prop = 0
prob_of_adding_response = 1 #CHANGE THIS
prob_of_uncensored_survival_observation = 0.3
max_n_dataset = 500
weibull_fixed_scale = 1
visualizations_during_testing = FALSE
ns = c(100)
Nrep = 501

source("_dataset_load.R")

for (dataset_name in names(datasets_and_response_models)){
  datasets_and_response_models[[dataset_name]]$response_models = list()
  for (response_type in names(datasets_and_response_models[[dataset_name]]$beta_T)){
    y_orig = datasets_and_response_models[[dataset_name]]$y_original[[response_type]]
    datasets_and_response_models[[dataset_name]]$response_models[[response_type]] = switch(response_type,
       continuous =  lm(y_orig ~ ., datasets_and_response_models[[dataset_name]]$X),
       incidence =   glm(y_orig ~ ., datasets_and_response_models[[dataset_name]]$X, family = "binomial"),
       proportion =  betareg(y_orig ~ ., datasets_and_response_models[[dataset_name]]$X),
       count =       MASS::glm.nb(y_orig ~ ., datasets_and_response_models[[dataset_name]]$X),
       survival =    survreg(Surv(y_orig, rep(1, length(y_orig))) ~ ., datasets_and_response_models[[dataset_name]]$X, scale = weibull_fixed_scale)
    )
  }
}
rm(y_orig, dataset_name, response_type)
#now we fit all the models for the parametric bootstrap


# log_odds = function(p){log(p / (1 - p))}
# expit = function(log_odds){exp(log_odds) / (1 + exp(log_odds))}
draw_response_with_treatment = function(x_t, dead_t, w_t, beta_T, response_type){
  x_t = c(1, as.numeric(x_t)) #add a 1 for the intercept term
  mod = datasets_and_response_models[[dataset_name]]$response_models[[response_type]]
  b = switch(class(mod)[1],
        betareg = coef(summary(mod))$mu[,1],
        coef(mod)
      )
  
  #ensure we only predict on variables whose coefficients were estimated
  x_t = x_t[!is.na(b)]
  b = b[!is.na(b)]
  eta = b %*% x_t + beta_T * w_t
  exp_eta = exp(eta)
  expit_eta = exp(eta - log1pexp(eta)) #for numerical stability as exp(eta) can be very large

  # stop("boom")
  switch(response_type,
   continuous =  rnorm(1, eta, summary(mod)$sigma),
   incidence =   rbinom(1, 1, expit_eta),
   proportion =  {
     phi = coef(summary(mod))$phi[1]
     rbeta(1, expit_eta * phi, (1 - expit_eta) * phi)
   },
   count =       rnbinom(1, size = mod$theta, mu = exp_eta),
   survival =    rweibull(1, mod$scale, exp_eta)
  )
}


exp_settings = data.table(expand.grid(
  dataset_name = names(datasets_and_response_models),
  beta_T = c("zero", "nonzero"),
  response_type = c("continuous", "incidence", "count", "proportion", "survival"),
  d = c("KK21stepwise", "KK21", "KK14", "Atkinson", "Efron", "iBCRD", "CRD"), #c("KK14", "iBCRD"), #
  n = ns,
  n_rep = 1 : Nrep
))
check_dataset_has_response = function(dataset_name, response_type){
  response_type %in% names(datasets_and_response_models[[dataset_name]]$beta_T)
}
  
exp_settings[, has_response_type := .(Map(check_dataset_has_response, dataset_name, response_type))]
exp_settings = exp_settings[has_response_type == TRUE]
exp_settings[, has_response_type := NULL]

rm(check_dataset_has_response)

#test all response_types and all designs

# res = data.frame()


# handlers(global = TRUE)
# handlers("progress")
# 
# foreach_wrap <- function(n_settings_arr) {
#   progressor_obj <- progressor(along = n_settings_arr)

  
cl = makeCluster(nC)
registerDoParallel(cl)
clusterExport(cl, c('exp_settings', 'datasets_and_response_models', 'missing_data_prop', 'visualizations_during_testing', 'nsim_exact_test', 'prob_of_uncensored_survival_observation'))
  
res = foreach(
    nsim = 1:nrow(exp_settings),
    .inorder = FALSE,
    .combine = rbind,
    .packages = c("data.table", "SeqExpMatch", "qgam")
    # .errorhandling = "pass"
    # .verbose = TRUE,
  ) %dopar% {
#we can parallelize this loop after we're done debugging and Rcpp'ing
# for (nsim in 1 : nrow(exp_settings)){
  set.seed(nsim)
    
    cat("\n", nsim, "/", nrow(exp_settings), "\n")
    
  # exp_setting = exp_settings[nsim, ]
  # print(res)
  dataset_name = as.character(exp_settings$dataset_name[nsim])
  beta_T = exp_settings$beta_T[nsim]
  response_type = as.character(exp_settings$response_type[nsim])
  design = as.character(exp_settings$d[nsim])
  n = exp_settings$n[nsim]
  n_rep = exp_settings$n_rep[nsim]
  
  # if (response_type != "continuous"){
  #   next
  # }
  # if (!(dataset_name %in% c("boston", "airquality"))){
  #   next
  # }
  
  X = datasets_and_response_models[[dataset_name]]$X
  if (n > nrow(X)){
    next
  }
  
  X = copy(X[sample(1 : .N, n)])
  if (missing_data_prop > 0){
    #add some missing data
    for (i in 2 : n){
      for (j in 1 : ncol(Xorig)){
        if (runif(1) <= missing_data_prop){
          X[i, j] = NA
        }
      }
    }    
  }
  
  beta_T = ifelse(beta_T == "zero", 0, datasets_and_response_models[[dataset_name]][["beta_T"]][[response_type]])
  
  # if (!(response_type == "survival")){
  #   next
  # }
  
  cat("\n", nsim, "/", nrow(exp_settings), " dataset_name =", dataset_name, "response_type =", response_type, "design =", design, "\n")
  
  t_0 = proc.time()
  seq_des_obj = SeqDesign$new(n, design, response_type = response_type, verbose = FALSE)
  response_added = array(FALSE, n)
  # stop("boom")
  # profvis({
  for (t in 1 : n){
    w_t = seq_des_obj$add_subject_to_experiment_and_assign(X[t, ])
    # if (t %% 50 == 0){
    #   seq_des_obj$print_current_subject_assignment()
    # }
    if (runif(1) < prob_of_adding_response){
      dead_t = ifelse(response_type == "survival", as.numeric(runif(1) > prob_of_uncensored_survival_observation), 1)
      seq_des_obj$add_subject_response(t, draw_response_with_treatment(X[t], dead_t, w_t, beta_T, response_type), dead = dead_t)
      response_added[t] = TRUE
    }
  }
  # })
  assignment_time = proc.time() - t_0
  
  for (t in which(!response_added)){
    dead_t = ifelse(response_type == "survival", as.numeric(runif(1) > prob_of_uncensored_survival_observation), 1)
    suppressWarnings(seq_des_obj$add_subject_response(t, draw_response_with_treatment(X[t, ], dead_t, seq_des_obj$w[t], beta_T, response_type), dead = dead_t))
  }
  cat("experiment completed?", seq_des_obj$check_experiment_completed(), "\n")
  # mean(log(seq_des_obj$y[seq_des_obj$w == 1])) - mean(log(seq_des_obj$y[seq_des_obj$w == 0]))
  
  prop_subjects_matched = NA
  if (visualizations_during_testing & grepl("KK", design)){
    stats = seq_des_obj$matching_statistics()
    prop_subjects_matched = stats$prop_subjects_matched
    cat(paste("  KK stats: prop matches:", prop_subjects_matched, "# remaining in reservoir:", stats$num_subjects_remaining_in_reservoir, "\n"))
    if (grepl("KK21", design)){
      if (!is.null(seq_des_obj$covariate_weights)){
        wts = seq_des_obj$covariate_weights / sum(seq_des_obj$covariate_weights)
        ggplot_obj = ggplot(data.frame(x = names(wts), imp = wts)) +
          aes(x = x, y = imp) +
          geom_bar(stat = "identity") +
          ggtitle(paste("Dataset:", dataset_name, " Design:", design, " Response:", response_type)) +
          theme(axis.text.x = element_text(angle = 90)) +
          ylim(0, 0.6)
        print(ggplot_obj)
      }
    }
  }
  
  res = data.frame()
  # let's do some inference
  for (test_type in c("MLE-or-KM-based", "randomization-exact")){
    for (estimate_type in c(
      ########################################### CONTINUOUS
      "continuous_simple_mean_difference",
      "continuous_multivariate_regression",
      "continuous_KK_compound_mean_difference",
      "continuous_KK_compound_multivariate_regression",
      # "continuous_KK_regression_with_covariates_with_matching_dummies",
      # "continuous_KK_regression_with_covariates_with_random_intercepts",
      ########################################### INCIDENCE
      "incidence_simple_mean_difference",
      "incidence_simple_log_odds",	
      "incidence_multivariate_logistic_regression",
      "incidence_KK_compound_mean_difference",
      # "incidence_KK_compound_univariate_logistic_regression",
      # "incidence_KK_compound_multivariate_logistic_regression",
      # "incidence_KK_multivariate_logistic_regression_with_matching_dummies",
      # "incidence_KK_multivariate_logistic_regression_with_random_intercepts_for_matches",
      ########################################### PROPORTION
      "proportion_simple_mean_difference",
      "proportion_simple_logodds_regression",
      "proportion_multivariate_beta_regression",
      "proportion_KK_compound_mean_difference",
      # "proportion_KK_compound_univariate_beta_regression",
      # "proportion_KK_compound_multivariate_beta_regression",
      # "proportion_KK_multivariate_beta_regression_with_matching_dummies",
      ########################################### COUNT
      "count_simple_mean_difference",
      "count_KK_compound_mean_difference",
      # "count_univariate_negative_binomial_regression",
      "count_multivariate_negative_binomial_regression",
      #"count_KK_compound_univariate_negative_binomial_regression",	
      #"count_KK_compound_multivariate_negative_binomial_regression",
      # "count_KK_multivariate_negative_binomial_regression_with_matching_dummies",
      # "count_KK_multivariate_negative_binomial_regression_with_random_intercepts_for_matches",
      ########################################### SURVIVAL
      "survival_simple_mean_difference",
      "survival_KK_compound_mean_difference",
      "survival_simple_median_difference",	
      "survival_simple_restricted_mean_difference",
      # "survival_univariate_weibull_regression",
      "survival_multivariate_weibull_regression",
      # "survival_KK_compound_univariate_weibull_regression",	
      # "survival_KK_compound_multivariate_weibull_regression",
      #"survival_KK_multivariate_weibull_regression_with_matching_dummies",
      # "survival_univariate_coxph_regression",
      "survival_multivariate_coxph_regression"
      # "survival_KK_multivariate_coxph_regression_with_matching_dummies",		
      # "survival_KK_multivariate_coxph_regression_with_random_intercepts_for_matches"				
    )){
      
      #only do inference by appropriate response type
      if (strsplit(estimate_type, "_")[[1]][1] != response_type){
        next
      }
      #only do inference for KK if it's a KK design
      if (grepl("KK", estimate_type) & !grepl("KK", design)){
        next
      }
      
      t_0 = proc.time()
      
      seq_des_inf_obj = SeqDesignInference$new(seq_des_obj, estimate_type = estimate_type, test_type = test_type, num_cores = 1)
      beta_hat_T = seq_des_inf_obj$compute_treatment_estimate()
      cat("  beta_T = ", beta_T, "  beta_hat_T =", beta_hat_T, "\n")
      pval = seq_des_inf_obj$compute_two_sided_pval_for_treatment_effect(nsim_exact_test = nsim_exact_test)
      cat("  pval =", pval, "\n")
    
      if (test_type == "MLE-or-KM-based" & response_type == "continuous"){
        ci = seq_des_inf_obj$compute_confidence_interval(nsim_exact_test = nsim_exact_test)
        cat("  95% CI for betaT = [", paste(round(ci, 4), collapse = ", "), "] width =", round(ci[2] - ci[1], 4), "\n")
      } else {
        ci = c(NA, NA)
      }
      
      inference_time = proc.time() - t_0
      
      res = rbind(res, data.frame(
        nsim = nsim,
        n_rep = n_rep,
        dataset_name = dataset_name,
        n = n,
        response_type = response_type,
        design = design,
        prop_subjects_matched = prop_subjects_matched,
        test_type = test_type,
        estimate_type = estimate_type,
        beta_T = beta_T,
        beta_hat_T = beta_hat_T,
        ci_a = ci[1],
        ci_b = ci[2],
        p_val = pval,
        assignment_time = assignment_time[["elapsed"]],
        inference_time = inference_time[["elapsed"]]
      ))
    }
  }
  save(res, file = paste0("seq_exp_match_", nsim, ".RData"))
  # progressor_obj(sprintf("nsim=%g", nsim))
  res
  # }
}

# res = foreach_wrap(1:2)
fwrite(res, file = paste0("seq_exp_match_sims.csv"))



# Nsim = 1e3
# yts = array(NA, Nsim)
# ycs = array(NA, Nsim)
# for (nsim in 1 : Nsim){
#   yts[nsim] = draw_response_with_treatment(X[t, ], 1, dead_t, beta_T, response_type)
#   ycs[nsim] = draw_response_with_treatment(X[t, ], 0, dead_t, beta_T, response_type)
# }
# mean(yts)
# mean(ycs)

# res[response_type == "continuous" & estimate_type == "simple_mean_difference" & beta_T == 0 & design == "CRD"]



# cppFunction("
# 	LogicalVector which_cols_vary_cpp(NumericMatrix X){
# 	  int n = X.nrow();
# 	  int p = X.ncol();
# 	  LogicalVector which_cols_vary(p);
# 	  for (int j = 0; j < p; j++){
#   	  for (int i = 1; i < n; i++){
#   	    if (X(i, j) != X(0, j)){
#   	      which_cols_vary[j] = TRUE;
#   	      break;
#   	    }
#   	  }	    
# 	  }
# 	  return which_cols_vary;
# 	}			
# ")


# covariate_data_matrix = as.matrix(MASS::Boston[, 1:13])
# response_obj = MASS::Boston$medv
# 
# microbenchmark::microbenchmark(
#   regular = {
#     ols_mod = lm(response_obj ~ covariate_data_matrix)
#     abs(coef(suppressWarnings(summary(ols_mod)))[2, 3])
#   },
#   speedup = {
#     Xmat = cbind(1, covariate_data_matrix)
#     qr_decomp = qr(Xmat)
#     Qmat = qr.Q(qr_decomp)
#     Rmat = qr.R(qr_decomp)
#     Qmatt = t(Qmat)
#     Rmatt = t(Rmat)
#     RtRinv = solve(Rmatt %*% Rmat)
#     b = (RtRinv %*% Rmatt %*% Qmatt %*% response_obj)
#     s_sq_e = sum((response_obj - Xmat %*% b)^2) / (length(response_obj) - ncol(Qmat))
#     abs(b[2]) / sqrt(s_sq_e * RtRinv[2, 2])
#   },
#   times = 1e3
# )

rm(list = ls())
pacman::p_load(ggplot2, data.table, magrittr, stringi)
filenames = dir(pattern = "*.RData")
load(filenames[1])
res$prop_subjects_matched = NA_real_
res = data.table(res)
blank = data.table(matrix(NA, 1e6, ncol(res)))
setnames(blank, colnames(res))
all_res = rbind(res, blank)
i = nrow(res) + 1
for (f in 2 : length(filenames)){
  load(filenames[f])
  set(all_res, i : (i + nrow(res) - 1), j = colnames(res), value = res)
  i = i + nrow(res)
}
rm(blank, i, f, res)
all_res = all_res[!is.na(nsim)]
gc()
table(all_res$n_rep)


all_res[, design := factor(design, levels = c("KK21stepwise", "KK21", "KK14", "Atkinson", "Efron", "iBCRD", "CRD"))]
all_res[abs(beta_hat_T) > 100, beta_hat_T := NA]
all_res[is.na(beta_hat_T), ci_a := NA]
all_res[is.na(beta_hat_T) > 100, ci_b := NA]
all_res[, covers := as.numeric(ci_a <= beta_T & beta_T <= ci_b)]
all_res[, design := factor(design, levels = rev(c("KK21stepwise", "KK21", "KK14", "Atkinson", "Efron", "iBCRD", "CRD")))]

all_res[, estimate_type := stri_replace_all_regex(estimate_type, "continuous_regression_with_covariates", "continuous_multivariate_regression")]
all_res[, estimate_type := stri_replace_all_regex(estimate_type, "incidence_logistic_regression", "incidence_multivariate_logistic_regression")]
all_res[, estimate_type := stri_replace_all_regex(estimate_type, "proportion_beta_regression", "proportion_multivariate_beta_regression")]


all_res[, estimate_type := stri_replace_all_regex(estimate_type, "continuous_", "")]
all_res[, estimate_type := stri_replace_all_regex(estimate_type, "incidence_", "")]
all_res[, estimate_type := stri_replace_all_regex(estimate_type, "count_", "")]
all_res[, estimate_type := stri_replace_all_regex(estimate_type, "proportion_", "")]
all_res[, estimate_type := stri_replace_all_regex(estimate_type, "survival_", "")]
all_res[, resp_and_est := paste0(response_type, "_", estimate_type)]

table(all_res$resp_and_est)
all_res[, .(
  avg_assignment_time = mean(assignment_time, na.rm = TRUE)
), by = c("n", "design", "response_type")][order(design, response_type)]

all_res[, .(
  avg_inference_time = mean(inference_time, na.rm = TRUE)
), by = c("n", "test_type", "estimate_type", "response_type")]

res_summary_mle = all_res[test_type == "MLE-or-KM-based", .(
  dataset = first(dataset_name),
  estimate_type = first(estimate_type),
  response_type = first(response_type),
  beta_T_hat = mean(beta_hat_T, na.rm = TRUE),
  ase = mean((beta_hat_T - beta_T)^2, na.rm = TRUE),
  coverage = mean(ci_a <= beta_T & beta_T <= ci_b, na.rm = TRUE),
  p_val = mean(p_val, na.rm = TRUE), 
  avg_pow = mean(p_val < 0.05, na.rm = TRUE),
  p_val_prop_test = prop.test(x = sum(p_val < 0.05, na.rm = TRUE), n = sum(!is.na(p_val)), p = 0.05)$p.value,
  prop_NA = sum(is.na(p_val)) / .N,
  nrep = .N
), by = c("dataset_name", "design", "n", "resp_and_est", "beta_T")]
res_summary_mle[, is_KK := grepl("KK", estimate_type)]
res_summary_mle[, is_multi := grepl("multivariate", estimate_type)]

res_summary_no_effect = res_summary_mle[beta_T == 0]
res_summary_tx_effect = res_summary_mle[beta_T > 0]
res_summary_tx_effect[, p_val_prop_test := NULL]
res_summary_no_effect[, p_val_prop_test := p_val_prop_test < 0.05 / .N]
res_summary_no_effect[p_val_prop_test == TRUE]
sum(res_summary_no_effect$p_val_prop_test) / nrow(res_summary_no_effect)


summary(lm(avg_pow ~ design + response_type + dataset_name, res_summary_tx_effect))
summary(lm(avg_pow ~ design + estimate_type + dataset_name, res_summary_tx_effect))
summary(lm(avg_pow ~ design * resp_and_est + dataset_name, res_summary_tx_effect))
summary(lm(avg_pow ~ design + (is_KK + is_multi) + dataset_name, res_summary_tx_effect))
summary(lm(avg_pow ~ design * (is_KK + is_multi) + dataset_name, res_summary_tx_effect))

res_summary_tx_effect[order(dataset_name, estimate_type, design)][response_type == "incidence"][dataset_name == "boston"]


res_summary_rand = all_res[test_type == "randomization-exact", .(
  dataset = first(dataset_name),
  beta_T_hat = mean(beta_hat_T, na.rm = TRUE),
  p_val = mean(p_val, na.rm = TRUE), 
  avg_pow = mean(p_val < 0.05, na.rm = TRUE),
  p_val_prop_test = prop.test(x = sum(p_val < 0.05, na.rm = TRUE), n = length(!is.na(p_val)), p = 0.05)$p.value,
  prop_NA = sum(is.na(p_val)) / .N,
  nrep = .N
), by = c("dataset_name", "design", "n", "resp_and_est", "estimate_type", "response_type", "beta_T")]




res_summary_no_effect = res_summary_rand[beta_T == 0][, beta_T := NULL]
res_summary_tx_effect = res_summary_rand[beta_T > 0][, beta_T := NULL]
res_summary_tx_effect[, p_val_prop_test := NULL]
res_summary_no_effect[, p_val_prop_test := p_val_prop_test < 0.05 / .N]
res_summary_no_effect[p_val_prop_test == TRUE]
sum(res_summary_no_effect$p_val_prop_test) / nrow(res_summary_no_effect)


res_summary_tx_effect[order(dataset_name, estimate_type, design)][response_type == "survival"][dataset_name == "boston"]
res_summary_no_effect[order(dataset_name, estimate_type, design)][response_type == "survival"][dataset_name == "boston"]
res_summary_tx_effect[order(dataset_name, estimate_type, design)][response_type == "continuous"][dataset_name == "ionosphere"]
res_summary_no_effect[order(dataset_name, estimate_type, design)][response_type == "continuous"][dataset_name == "ionosphere"]
res_summary_tx_effect[, .(avg_power = mean(avg_pow), nrep = sum(nrep)), by = c("response_type", "estimate_type", "design")][order(response_type, estimate_type, design)] %>% data.frame

summary(lm(avg_pow ~ design + response_type + dataset_name, res_summary_tx_effect))
summary(lm(avg_pow ~ design + resp_and_est, res_summary_tx_effect))

summary(lm(avg_pow ~ design + estimate_type + dataset_name, res_summary_tx_effect))


