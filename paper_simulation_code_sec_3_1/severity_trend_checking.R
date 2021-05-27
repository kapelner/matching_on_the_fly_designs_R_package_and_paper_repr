#save(indic_Ts, file = "indic_Ts.RData")

#z-tests for all entrants
alpha_corr = 0.05
z = qnorm(1 - alpha_corr / 2)
l = 0.5 - z * 0.5 / sqrt(Nsim_per_block)
u = 0.5 + z * 0.5 / sqrt(Nsim_per_block)
l
u
prob_T_t_est = colMeans(indic_Ts)
rejections = prob_T_t_est < l | prob_T_t_est > u
which(rejections)
sum(rejections)
pvals = sort(sapply(prob_T_t_est, function(x){2 * (1 - pnorm(abs(mean(x) - 0.5) / (0.5 / sqrt(Nsim_per_block))))}))
pvals

#logistic regression to test time trend
indic_Ts[1,]
assignment_time_df = data.frame(
		t = rep(1:n, Nsim_per_block),
		w = c(t(indic_Ts))
)
summary(glm(w ~ t, assignment_time_df, family = "binomial"))

assignment_time_df = data.frame(
		t = rep((t_0_matching_pcts[1] * n) : n, Nsim_per_block),
		w = c(t(indic_Ts[, (t_0_matching_pcts[1] * n):n]))
)
summary(glm(w ~ t, assignment_time_df, family = "binomial"))

#C-A test
contingency_table = data.frame(
		Ts = colSums(indic_Ts),
		Cs = Nsim_per_block - colSums(indic_Ts)
)
pacman::p_load(DescTools)
CochranArmitageTest(contingency_table, "increasing")
CochranArmitageTest(contingency_table, "decreasing")



# dose <- matrix(c(10,9,10,7, 0,1,0,3), byrow=TRUE, nrow=2, dimnames=list(resp=0:1, dose=0:3))
# Desc(dose)
# 
# CochranArmitageTest(dose, "increasing")
# CochranArmitageTest(dose)

#KS-test
ks.test(prob_T_t_est[1 : (t_0_matching_pcts[1] * n)], prob_T_t_est[(t_0_matching_pcts[1] * n + 1) : n])

binned_contingency_table = matrix(NA, nrow = 10, ncol = 2)
idx = seq(0, 90, by = 10)
for (i in 1 : length(idx)){
	binned_contingency_table[i, ] = apply(contingency_table[(idx[i] + 1) : (idx[i] + 10), ], 2, sum)
}
binned_contingency_table

pacman::p_load(ggplot2)
alpha_corr = 0.05 / 100
z_crit = qnorm(1 - alpha_corr / 2)
ggplot(data.frame(
						t = 1 : n,
						Design = c(rep("Bernoulli", n * t_0_matching_pcts[1]), rep("Our CARA", n - n * t_0_matching_pcts[1])),
						prob_T_t_est = prob_T_t_est, 
						l = prob_T_t_est - z_crit * sqrt(prob_T_t_est * (1 - prob_T_t_est) / Nsim_per_block), 
						u = prob_T_t_est + z_crit * sqrt(prob_T_t_est * (1 - prob_T_t_est) / Nsim_per_block)
				)) + geom_point(aes(x = t, y = prob_T_t_est, col = Design)) +
		geom_errorbar(aes(x = t, y = prob_T_t_est, ymin = l, ymax = u, col = Design), width=1) +
		geom_hline(yintercept = 0.5, col = "black")

