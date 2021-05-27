
for (n in ns_to_test){
	for (beta_T in treatment_effects){
	  for (all_betas_and_correlation in all_betas_and_correlations){

		
		  all_beta_hats = list()
    	  betas = all_betas_and_correlation[["betas"]]
    	  beta_string = paste(betas, collapse = "_")
    	  rho = all_betas_and_correlation[["rho"]] 
		  for (randomization_type in randomization_types){
			for (prob_match_cutoff_lambda in prob_match_cutoff_lambdas){
				for (t_0_matching_pct in t_0_matching_pcts){
					for (z_test in Z_TESTS){
					  num_simulation_setting = num_simulation_setting + 1
						cat("SIM BLOCK:", num_simulation_setting, "/", num_simulation_settings,
							# "response_model =", response_model,
						    "n =", n, 
						    "betaT =", beta_T,
						    "betas =", beta_string, 
						    "rho =", rho, 
						    "alg =", randomization_type, 
							"lambda =", prob_match_cutoff_lambda, 
							"t_0% =", t_0_matching_pct,
							"z_test =", z_test
						) #, "\n"
						
						
						#each prob_match and response gets their own results matrix and it gets saved to a different CSV file
						results = matrix(NA, nrow = length(metrics_for_each_run), ncol = 1)
						rownames(results) = array(NA, length(metrics_for_each_run))
						#add row names to the results matrix
						for (m in 1 : length(metrics_for_each_run)){
							rownames(results) = metrics_for_each_run	
						}
						
						matching_algorithm = FALSE ##default, and overwritten by rand_type code
						#run the simulation series!
						source("inner_simulation.R")
						#save the results iteratively
						
						master_results = rbind(master_results, data.frame(
						    algorithm = randomization_type,
							n = n,
							beta_T = beta_T,
							betas = beta_string,
							rho = rho,
							cutoff = ifelse(matching_algorithm, prob_match_cutoff_lambda, NA),
							t_0_pct = ifelse(matching_algorithm, t_0_matching_pct, NA),
							t(results)
						))
						#now write out one row to the file	
#								stop("boom")
						write.csv(master_results, paste0(results_filepath, "_n_", n, "_betaT_", beta_T, "_betas_", beta_string, "_rho_", rho, ".csv"), row.names = FALSE)		
						
						num_simulation_settings_left = num_simulation_settings - num_simulation_setting
						sec_elapsed = as.numeric(difftime(as.POSIXct(Sys.time()), as.POSIXct(time_began), unit = "sec"))
						sec_per_simulation_setting = sec_elapsed / num_simulation_setting
						sec_estimated_left = sec_per_simulation_setting * num_simulation_settings_left
						
						min_elapsed = round(sec_elapsed / 60)
						# cat("  t diff:", Sys.time() - time_began, "\n")
						cat("  # min elapsed:", min_elapsed, "\n")
						cat("  # min est remaining:", round(sec_estimated_left / 60), "\n")
						#we don't need to run these inner loops many times since they will be the same
						if (!matching_algorithm){
						  break
						}
					}
				  if (!matching_algorithm){
				    break
				  }						  
				}
			  if (!matching_algorithm){
			    break
			  }					  
			}
		  }
	  }
	}
}