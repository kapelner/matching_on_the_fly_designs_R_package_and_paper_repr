block_results = data.frame(
  max_std_diff_balance = numeric(),
  max_ks_stat = numeric(),
  beta_hat_T = numeric(),
  Rsq = numeric(),
  Ha_acceptance = numeric(),
  pct_T = numeric(),
  conv_guessing_strategy_pct_correct = numeric(),
  #metrics only for the sequential matching algorithms
  T_stat = numeric(),
  final_reservoir_size = numeric(),
  ssqr_over_sum = numeric(),
  ssqr_eq_ssqdbar_pval = numeric(),
  matched_correlation = numeric(),
  pct_only_matches = numeric(),
  pct_only_reservoir = numeric(),
  true_var_prop_diff = numeric()
)
indic_Ts = matrix(NA, nrow = Nsim_per_block, ncol = n)


# block_results = foreach(nsim = 1 : Nsim_per_block, .inorder = FALSE, .combine = rbind) %dopar% {
#   if (nsim %% 100 == 0){
#     cat(".")
#   }
# 
#   Sigma = sigma_x * (matrix(rho, nrow = p, ncol = p) + diag(1 - rho, p))
#   #build mvnp covariates
#   x_s = MASS::mvrnorm(n, rep(mu_x, p), Sigma)
#   #build errors independent of x (and thus a priori of treatment allocation)
#   errors = rnorm(n, 0, sigma_e)
#   #sys.source("common_crd.R", envir = environment(), toplevel.env = environment())
#   indic_T = rbinom(n, 1, prob_trt)
#   sys.source(paste("create_response_", response_model, ".R", sep = ""), envir = environment(), toplevel.env = environment())
#   #create design matrix
#   Xy = data.frame(cbind(x_s, indic_T, y))
#   colnames(Xy) = c(paste0("x", 1 : p), "indic_T", "y")
# 
#   #pull out yT, yC
#   yTs = Xy[Xy$indic_T == 1, "y"]
#   yCs = Xy[Xy$indic_T == 0, "y"]
#   betas
#   NA
# }


#the doparallel requires "uploading" relevant variables to the nodes' local environments
clusterExport(cl, list("beta_T", "prob_trt", "response_model", "betas", "t_0_matching_pct", "prob_match_cutoff_lambda", "z_test", "Nsim_exact_test", "NUM_BOOT_DSQD_DIST", "stepwise_weights", "find_rsqs"), envir = environment())
#now do the parallelization
block_results = foreach(nsim = 1 : Nsim_per_block, .inorder = FALSE, .combine = rbind) %dopar% {


#for (nsim in 1 : Nsim_per_block){
  if (nsim %% 100 == 0){
    cat(".")
  }
      
      
  # print(c(stepwise_weights, find_rsqs, all_betas, beta_T, betas, colnames_master_results, data_and_err_gen, find_rsqs, m, master_results, matching_algorithm, metrics_for_each_run, n, ns_to_test, Nsim_exact_test, Nsim_per_block, NUM_BOOT_DSQD_DIST, NUM_CORES_R, p, prob_match_cutoff_lambda, prob_match_cutoff_lambdas, prob_trt, randomization_type, randomization_types, response_model, results, results_filepath, sim_types, t_0_matching_pct, t_0_matching_pcts, treatment_effects, z_test, Z_TESTS))
  # stop("BOOOOOM  ", paste(ls(), collapse = ", "))
      
	#generate data
	#make exchangeable covariance matrix
	Sigma = sigma_x * (matrix(rho, nrow = p, ncol = p) + diag(1 - rho, p))
	#build mvnp covariates
	x_s = MASS::mvrnorm(n, rep(mu_x, p), Sigma)
	#build errors independent of x (and thus a priori of treatment allocation)
	errors = rnorm(n, 0, sigma_e)

	#run one run of whatever simulation type
  sys.source(paste("rand_type_", randomization_type, ".R", sep = ""), envir = environment())
  
	#run the convergent guessing strategy
	sys.source("convergent_guessing_strategy.R", envir = environment())
  
  # stop("BOOOOOM  ", paste(ls(), collapse = ", "))
    
	#save indic_Ts for debugging purposes
	indic_Ts[nsim, ] = Xy$indic_T
	
	#calculate balance metrics
	xTs = Xy[Xy$indic_T == 1, 1 : p]
	xCs = Xy[Xy$indic_T == 0, 1 : p]
	
	#compute balance	
	max_std_diff_balance = -.Machine$double.xmax
	for (j in 1 : p){
		std_diff_balance = abs(mean(xTs[, j]) - mean(xCs[, j])) / sqrt(var(xTs[, j]) / length(xTs[, j]) + var(xCs[, j]) / length(xCs[, j]))
		if (std_diff_balance > max_std_diff_balance){
		  max_std_diff_balance = std_diff_balance
		}
	}
	max_ks_stat = -.Machine$double.xmax
	for (j in 1 : p){
		ks_stat = ks.test(xTs[, j], xCs[, j])$statistic
		if (ks_stat > max_ks_stat){
			max_ks_stat = ks_stat
		}
	}
	# cat("   nsim:", nsim, "\n")
#	block_results = rbind(block_results, data.frame(
	 data.frame(
	  max_std_diff_balance = max_std_diff_balance,
	  max_ks_stat = max_ks_stat,
	  beta_hat_T = beta_hat_T,
	  Rsq = is_defined_or_NA(Rsq),
	  Ha_acceptance = ifelse(pval < 0.05, 1, 0),
	  pct_T = nrow(xTs) / n,
	  conv_guessing_strategy_pct_correct = sum(indic_T_guesses == indic_T) / n,
	  #metrics only for the sequential matching algorithms
	  T_stat = is_defined_or_NA(T_stat),
	  final_reservoir_size = is_defined_or_NA(final_reservoir_size),
	  ssqr_over_sum = is_defined_or_NA(ssqr_over_sum),
	  ssqr_eq_ssqdbar_pval = is_defined_or_NA(ssqr_eq_ssqdbar_pval),
	  matched_correlation = is_defined_or_NA(matched_correlation),
	  pct_only_matches = is_defined_or_NA(pct_only_matches),
	  pct_only_reservoir = is_defined_or_NA(pct_only_reservoir),
	  true_var_prop_diff = is_defined_or_NA(true_var_prop_diff),
	  conv_guessing_strategy_pct_correct_after_n_0 = 
	  	ifelse(is_defined(min_t0_based_on_matching_pct),
			sum(indic_T_guesses[min_t0_based_on_matching_pct : n] == indic_T[min_t0_based_on_matching_pct : n]) / (n - min_t0_based_on_matching_pct + 1),
			NA
		)
	 )
#	)
}

cat("\n")

#add to results
all_beta_hats[[randomization_type]] = block_results$beta_hat_T #save all of these so we can try to understand the distr
results["avg_max_std_diff_bal", 1] = mean(block_results$max_std_diff_balance, na.rm = TRUE)
results["avg_beta_T", 1] = mean(block_results$beta_hat_T, na.rm = TRUE)	
results["avg_abs_bias", 1] = mean(abs(block_results$beta_hat_T - beta_T), na.rm = TRUE)
results["avg_max_ks_stat", 1] = mean(block_results$max_ks_stat, na.rm = TRUE)
results["std_err_beta_T", 1] = sd(block_results$beta_hat_T, na.rm = TRUE)
results["power", 1] = mean(block_results$Ha_acceptance, na.rm = TRUE)
results["pct_trt_diff", 1] =  mean(abs(block_results$pct_T - prob_trt), na.rm = TRUE)
results["conv_guessing_strategy_pct_correct_avg", 1] = mean(block_results$conv_guessing_strategy_pct_correct, na.rm = TRUE) 
results["conv_guessing_strategy_pct_correct_sd", 1] = sd(block_results$conv_guessing_strategy_pct_correct, na.rm = TRUE) 

#only save reservoir data strategies that involve a matching algorithm
if (matching_algorithm){
	results["res_end_prop_avg", 1] = mean(block_results$final_reservoir_size)
#	results["ssqr_over_sum", 1] = mean(block_results$ssqr_over_sum, na.rm = TRUE)
#	results["ssqr_eq_ssqd_pval", 1] = mean(block_results$ssqr_eq_ssqdbar_pval, na.rm = TRUE)
	results["pct_only_matches", 1] = mean(block_results$pct_only_matches, na.rm = TRUE)
	results["pct_only_reservoir", 1] = mean(block_results$pct_only_reservoir, na.rm = TRUE)
	results["match_corr", 1] = mean(block_results$matched_correlation, na.rm = TRUE)
#	results["true_var_prop_diff", 1] = mean(block_results$true_var_prop_diff, na.rm = TRUE)
	results["conv_guessing_strategy_pct_correct_after_n_0_avg", 1] = mean(block_results$conv_guessing_strategy_pct_correct_after_n_0, na.rm = TRUE)
	results["conv_guessing_strategy_pct_correct_after_n_0_sd", 1] = sd(block_results$conv_guessing_strategy_pct_correct_after_n_0, na.rm = TRUE)
	
}
