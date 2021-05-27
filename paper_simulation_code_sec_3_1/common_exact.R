#the observed value is the sample avg difference
beta_hat_T = mean(yTs) - mean(yCs)

## now we have to monte-carlo the exact test
b_T_sims = array(NA, Nsim_exact_test)
for (nsim_exact_test in 1 : Nsim_exact_test){
	permuted_indic_T = do.call(indic_T_permute_function, indic_T_permute_function_args)
	yTs = Xy[permuted_indic_T == 1, "y"]
	yCs = Xy[permuted_indic_T == 0, "y"]
	
	b_T_sims[nsim_exact_test] = mean(yTs) - mean(yCs)
}

#hist(b_T_sims, br = 100)
#mean(b_T_sims)
#sum(b_T_sims_unc > 1) / Nsim_exact_test
#sum(b_T_sims_cond > 1) / Nsim_exact_test
#hist(b_T_sims_unc, br = 100)
#hist(b_T_sims_cond, br = 100)
#ks.test(b_T_sims_unc, b_T_sims_cond)

#this is the empirical two-sided p-value based on simulation
pval = sum(abs(beta_hat_T) < abs(b_T_sims)) / Nsim_exact_test