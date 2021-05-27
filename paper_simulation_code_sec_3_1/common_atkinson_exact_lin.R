#the observed value is the OLS estimate
b_T_obs = coef(lm(y ~ ., Xy))[p + 2]

Xy_copy = Xy
## now we have to monte-carlo the exact test
b_T_sims = array(NA, Nsim_exact_test)
for (nsim_exact_test in 1 : Nsim_exact_test){
	Xy_copy$indic_T = atkinson_assignment(x_s)
	#rerun regression
	b_T_sims[nsim_exact_test] = coef(lm(y ~ ., Xy_copy))[p + 2]
}

#hist(b_T_sims, br = 100)
#mean(b_T_sims)
#sum(b_T_sims_unc > 1) / Nsim_exact_test
#sum(b_T_sims_cond > 1) / Nsim_exact_test
#hist(b_T_sims_unc, br = 100)
#hist(b_T_sims_cond, br = 100)
#ks.test(b_T_sims_unc, b_T_sims_cond)

#this is the empirical two-sided p-value based on simulation
beta_hat_T = NA
pval = sum(abs(b_T_obs) < abs(b_T_sims)) / Nsim_exact_test