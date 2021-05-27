sys.source("common_stratification.R", envir = environment())

#get beta_hat_T from a regression
linear_mod = lm(y ~ . - blocks, data = Xy)

coefs = coef(summary(linear_mod))
beta_hat_T = coefs[p + 2, 1]
pval = coefs[p + 2, 4]
Rsq = summary(linear_mod)$r.squared