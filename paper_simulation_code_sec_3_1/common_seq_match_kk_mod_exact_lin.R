## now we have to monte-carlo the exact test
Xyleft_copy = Xyleft
Xymatched_copy = Xymatched
b_T_sims = array(NA, Nsim_exact_test)

#sometimes the reservoir just isn't large enough so just use the matches
if (nRT <= 2 || nRC <= 2 || nRT + nRC <= p + 2){
	
	for (nsim_exact_test in 1 : Nsim_exact_test){
		#### unconditional permuting - flip coins for each (equivalent to conditional permuting)
		for (im in 1 : m){
			Xymatched_copy[Xymatched_copy$match_indic == im, "indic_T"] = sample(c(0, 1)) #permute
		}	
		
		XyD = matrix(NA, nrow = m, ncol = ncol(Xymatched))
		colnames(XyD) = colnames(Xymatched)
		
		for (im in 1 : m){
			Xyim = Xymatched_copy[Xymatched_copy$match_indic == im, ]
			XyimT = Xyim[Xyim$indic_T == 1, ]
			XyimC = Xyim[Xyim$indic_T == 0, ]
			XyD[im, ] = as.numeric(XyimT - XyimC)
		}
		linear_mod_matched = lm(y ~ . - match_indic - indic_T - z, data = as.data.frame(XyD))
		coefs = coef(summary(linear_mod_matched))
		
		b_T_sims[nsim_exact_test] = coefs[1, 1]
	}
	
} else if (m == 0){ #sometimes there's no matches
	for (nsim_exact_test in 1 : Nsim_exact_test){
		Xyleft_copy$indic_T = sample(c(rep(1, nRT), rep(0, nRC))) #permute
		linear_mod_reservoir = lm(y ~ . - match_indic - z, data = Xyleft_copy)		
		coefs_reservoir = coef(summary(linear_mod_reservoir))		
		b_T_sims[nsim_exact_test] = coefs_reservoir[p + 2, 1]
	}
} else {
	#compute estimator from matched pairs by regression
	for (nsim_exact_test in 1 : Nsim_exact_test){
		#### unconditional permuting - flip coins for each (equivalent to conditional permuting)
		for (im in 1 : m){
			Xymatched_copy[Xymatched_copy$match_indic == im, "indic_T"] = sample(c(0, 1)) #permute
		}	
		
		XyD = matrix(NA, nrow = m, ncol = ncol(Xymatched))
		colnames(XyD) = colnames(Xymatched)
		
		for (im in 1 : m){
			Xyim = Xymatched_copy[Xymatched_copy$match_indic == im, ] #permute
			XyimT = Xyim[Xyim$indic_T == 1, ]
			XyimC = Xyim[Xyim$indic_T == 0, ]
			XyD[im, ] = as.numeric(XyimT - XyimC)
		}
		linear_mod_matched = lm(y ~ . - match_indic - indic_T - z, data = as.data.frame(XyD))
		coefs_matched = coef(summary(linear_mod_matched))
		beta_match_regression = coefs_matched[1, 1]	
		ssqd_match_regression = coefs_matched[1, 2]^2 #lin mod returns SE not VAR, so square it
				
		#compute estimator reservoir sample std error
		Xyleft_copy$indic_T = sample(c(rep(1, nRT), rep(0, nRC))) #permute
		linear_mod_reservoir = lm(y ~ . - match_indic - z, data = Xyleft_copy)
		coefs_reservoir = coef(summary(linear_mod_reservoir))
		beta_reservoir_regression = coefs_reservoir[p + 2, 1]
		ssqd_reservoir_regression = coefs_reservoir[p + 2, 2]^2 #lin mod returns SE not VAR, so square it
		w_star = ssqd_reservoir_regression / (ssqd_reservoir_regression + ssqd_match_regression) #just a convenience for faster runtime	

		b_T_sims[nsim_exact_test] = w_star * beta_match_regression + (1 - w_star) * beta_reservoir_regression #proper weighting	
	}
}

#this is the empirical two-sided p-value based on simulation
pval = sum(abs(b_T_sims) > abs(beta_hat_T)) / Nsim_exact_test #beta_hat_T comes from common_seq_match_kk_mod_ols.R
