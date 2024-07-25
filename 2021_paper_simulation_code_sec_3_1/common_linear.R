#get beta_hat_T from a regression
linear_mod = lm(y ~ ., data = Xy)

coefs = coef(summary(linear_mod))
beta_hat_T = coefs[p + 2, 1]
#cat("   s_beta_T_hat:", coefs[p + 2, 2], "\n")
pval = coefs[p + 2, 4]
Rsq = summary(linear_mod)$r.squared
#if (nsim == 500){
#	stop("boooooom")
#}
