
#sometimes the reservoir just isn't large enough
if (nRT <= 2 || nRC <= 2 || nRT + nRC <= p + 2){
	XyD = matrix(NA, nrow = m, ncol = ncol(Xymatched))
	colnames(XyD) = colnames(Xymatched)
	
	for (im in 1 : m){
		Xyim = Xymatched[Xymatched$match_indic == im, ]
		XyimT = Xyim[Xyim$indic_T == 1, ]
		XyimC = Xyim[Xyim$indic_T == 0, ]
		XyD[im, ] = as.numeric(XyimT - XyimC)
	}
	linear_mod_matched = lm(y ~ . - match_indic - indic_T - z, data = as.data.frame(XyD))
	coefs = coef(summary(linear_mod_matched))
	
	beta_hat_T = coefs[1, 1]
	T_stat = coefs[1, 3]
	pval = coefs[1, 4]
	pct_only_matches = 1
	pct_only_reservoir = 0
	
#and sometimes there's no matches	
} else if (m == 0){ 
	linear_mod_reservoir = lm(y ~ . - match_indic - z, data = Xyleft)
	
	coefs_reservoir = coef(summary(linear_mod_reservoir))

	beta_hat_T = coefs_reservoir[p + 2, 1]
	T_stat = coefs_reservoir[p + 2, 3]
	pval = coefs_reservoir[p + 2, 4]
	pct_only_matches = 0
	pct_only_reservoir = 1

#but most of the time... we have matches and a goodly-sized reservoir
} else {
	#compute estimator from matched pairs by regression
	XyD = matrix(NA, nrow = m, ncol = ncol(Xymatched))
	colnames(XyD) = colnames(Xymatched)
	
	for (im in 1 : m){
		Xyim = Xymatched[Xymatched$match_indic == im, ]
		XyimT = Xyim[Xyim$indic_T == 1, ]
		XyimC = Xyim[Xyim$indic_T == 0, ]
		XyD[im, ] = as.numeric(XyimT - XyimC)
	}
	
	linear_mod_matched = lm(y ~ . - match_indic - indic_T - z, data = as.data.frame(XyD))
	coefs_matched = coef(summary(linear_mod_matched))
	beta_match_regression = coefs_matched[1, 1]	
	ssqd_match_regression = coefs_matched[1, 2]^2 #lin mod returns SE not VAR, so square it
	
	
	#compute estimator reservoir sample std error
	linear_mod_reservoir = lm(y ~ . - match_indic - z, data = Xyleft)
	
	coefs_reservoir = coef(summary(linear_mod_reservoir))
	beta_reservoir_regression = coefs_reservoir[p + 2, 1]
	ssqd_reservoir_regression = coefs_reservoir[p + 2, 2]^2 #lin mod returns SE not VAR, so square it
	
	w_star = ssqd_reservoir_regression / (ssqd_reservoir_regression + ssqd_match_regression) #just a convenience for faster runtime	
	ssqr_over_sum = w_star	
	
	beta_hat_T = w_star * beta_match_regression + (1 - w_star) * beta_reservoir_regression #proper weighting	
	
	
	b_T_est_sd = sqrt(ssqd_match_regression * ssqd_reservoir_regression / (ssqd_match_regression + ssqd_reservoir_regression)) #analagous eq's
	
	T_stat = beta_hat_T / b_T_est_sd #we calculate the T-stat against the null of zero effect	
	
	if (z_test){
		pval = 2 * (1 - pnorm(abs(T_stat))) #approximate by using real Z
	} else {
		pval = 2 * (1 - pt(abs(T_stat), max(m - 1, min(nRT, nRC)))) #approximate by using a T
	}
	
	ssqr_over_sum = ssqd_reservoir_regression / (ssqd_reservoir_regression + ssqd_match_regression)	
	ssqr_eq_ssqdbar_pval = pf(ssqd_reservoir_regression * 2 / ssqd_match_regression, nRT + nRC - 2, m - 1, lower.tail = F)
	pct_only_matches = 0
	pct_only_reservoir = 0	
	yMT = Xymatched[Xymatched$indic_T == 1, ]$y
	yCT = Xymatched[Xymatched$indic_T == 0, ]$y
	matched_correlation = cor(yMT, yCT)	
}
