#now get beta_hat_T via a regular classic pooled 2-sample t-test
ttest = t.test(yTs, yCs)
beta_hat_T = mean(yTs) - mean(yCs)
pval = ttest$p.value
Rsq = 1 - sum((Xy$y - beta_hat_T * Xy$indic_T)^2) / sum((Xy$y - mean(Xy$y))^2)