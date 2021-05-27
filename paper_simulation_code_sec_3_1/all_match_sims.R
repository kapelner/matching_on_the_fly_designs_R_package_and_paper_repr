pacman::p_load(MASS, iterators, tidyverse, lobstr, data.table, magrittr, microbenchmark, doParallel)

source("helper_functions.R")

NUM_CORES_R = 64
cl = makeCluster(NUM_CORES_R)
registerDoParallel(cl)
#registerDoSEQ(cl)


########set simulation parameters here
results_filepath = "~" # "results/results"

#how precise do we want our simulated results?
Nsim_per_block = 10000 #1000
Nsim_exact_test = 501 #501 #501
NUM_BOOT_DSQD_DIST = 500 #500 #500

#the treatment effects: 0 tests the size of test and something \neq 0 tests power
treatment_effects = c(
  1,
  0
)

#do we use the Z test and/or the T test (only when applicable)
Z_TESTS = c(TRUE)

#balanced assignment of treatment and control
prob_trt = 0.5

#how many subjects enter the sequential experiment?
ns_to_test = c( 
	# 50,
	100
	# 200
)

mu_x = 1
sigma_x = 1
sigma_e = 1

#the number of dimensions is fixed by simulation
p = 2

#how should y be generated as a function of the x's and noise?
response_model = "linear_quad_and_interact_medley"
############################BE VERY CAREFUL!!!!!!
#response_model = "linear_quad_and_interact_medley_severity_sorted"

#what kind of models are we running?
sim_types = c("series_of_quadratics_and_interactions")


#cutoff parameter for matching
prob_match_cutoff_lambdas = c(
	#  0.025,
	# 0.05
	# 0.075,
	0.10
	# 0.125,
	# 0.15,
	# 0.175,
	# 0.2
)

t_0_matching_pcts = c(
#	0.05,
	# 0.1,
	# 0.2,
	0.35
	# 0.4,
	# 0.5,
	# 0.6
)

#correlated x's is handled inside of these settings
all_betas_and_correlations = list()
# all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0, betas = c(0, 0, 0, 0, 0)))) #ZERO EFFECTS
# all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0, betas = c(1, 1, 0, 0, 0)))) #LINEAR EVEN
# all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0.75, betas = c(1, 1, 0, 0, 0)))) #LINEAR EVEN with CORR
# all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0, betas = c(2, 1, 0, 0, 0)))) #LINEAR UNEVEN
# all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0, betas = c(2, 1, 1, 0, 0)))) #QUAD MORE UNEVEN
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0, 	betas = c(1, 1, 1, 0, 0)))) #QUAD EVEN
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0.75, 	betas = c(1, 1, 1, 0, 0)))) #QUAD MORE UNEVEN with CORR
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0, 	betas = c(6, 1, 2, 0, 0)))) #QUAD MORE UNEVEN
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0.75, 	betas = c(6, 1, 2, 0, 0)))) #QUAD MORE UNEVEN with CORR

#How do we randomize subjects into treatment or control? 
#Then, how do we analyze the resulting data to obtain an effect size and significance level?
#Everything uncommented becomes part of the simulation
randomization_types = c(
 "crd_ttest",						#completely randomized design (CRD)							#analyzed via t-test
 "crd_lin",						#completely randomized design (CRD)							#analyzed via OLS
 "crd_exact",						#completely randomized design (CRD)							#analyzed via classic permutation test
 "crd_exact_lin",					#completely randomized design (CRD)							#analyzed via classic permutation test with OLS estimate

 "bcrd_ttest",					#completely randomized design (CRD)							#analyzed via t-test
 "bcrd_lin",						#completely randomized design (CRD)							#analyzed via OLS
 "bcrd_exact",					#completely randomized design (CRD)							#analyzed via classic permutation test
 "bcrd_exact_lin",				#completely randomized design (CRD)							#analyzed via classic permutation test with OLS estimate

 "atkinson_ttest",				#Atkinson's optimal design									#analyzed via t-test
 "atkinson_lin",					#Atkinson's optimal design									#analyzed via OLS
 "atkinson_exact",					#Atkinson's optimal design									#analyzed via OLS
 "atkinson_exact_lin",					#Atkinson's optimal design									#analyzed via OLS

 "efron_ttest",					#Efron's biased coin alpha=2/3 design						#analyzed via t-test
 "efron_lin",						#Efron's biased coin alpha=2/3 design						#analyzed via OLS
 "efron_exact",					#Efron's biased coin alpha=2/3 design						#analyzed via classic permutation test
 "efron_exact_lin",				#Efron's biased coin alpha=2/3 design
#   # "strat_ttest",				#stratification by tertiles => 9 blocks design  			#analyzed via t-test
#   # "strat_lin",					#stratification by tertiles => 9 blocks design  			#analyzed via OLS
# # "strat_exact",					#stratification by tertiles => 9 blocks design  			#analyzed via classic permutation test
# # "strat_exact_lin",				#stratification by tertiles => 9 blocks design  			#analyzed via classic permutation test with OLS estimate
 
 "ps_min_ttest",					#Pocock & Simon's minimization design						#analyzed via t-test
 "ps_min_lin",					#Pocock & Simon's minimization design						#analyzed via OLS
#  # "ps_min_exact",				#Pocock & Simon's minimization design						#analyzed via classic permutation test with OLS estimate
#  # "ps_min_exact_lin",			#Pocock & Simon's minimization design						#analyzed via classic permutation test
# 

"seq_match_kk_ttest",					#KK14 design, t-test  							#analyzed via classic test of Section 2.3.1, Equation 4
"seq_match_kk_lin",	    #KK14 design, OLS exact test
"seq_match_kk_exact",		#KK14 design, exact test  					#analyzed via permutation test of Section 2.3.3
"seq_match_kk_exact_lin",		#KK14 design, exact test  					#analyzed via permutation test of Section 2.3.3


"seq_match_weighted_bootstrap_ttest",		#KK20 naive design  					#analyzed via permutation test of Section 2.3.3
"seq_match_weighted_bootstrap_lin",		#KK20 naive design   					#analyzed via permutation test of Section 2.3.3
"seq_match_weighted_bootstrap_exact",		#KK20 naive design  					#analyzed via permutation test of Section 2.3.3
"seq_match_weighted_bootstrap_exact_lin",		#KK20 naive design   					#analyzed via permutation test of Section 2.3.3 with OLS estimate

 "seq_match_weighted_bootstrap_stepwise_ttest", #KK20 design
 "seq_match_weighted_bootstrap_stepwise_lin", #KK20 design
 "seq_match_weighted_bootstrap_stepwise_exact", #KK20 design
 "seq_match_weighted_bootstrap_stepwise_exact_lin" #KK20 design

		
)

##########deprecated

# "seq_match_weighted_bootstrap_exact_orthog",
# "seq_match_weighted_bootstrap_exact_orthog_lin",
# "seq_match_weighted_bootstrap_orthog",
#	"post_match_crd_kk",			#Post Matching Design with CRD assignment					#analyzed via classic test of Section 2.3.1, Equation 4
#	"post_match_strat_kk",			#Post Matching Design with stratified assignment			#analyzed via classic test of Section 2.3.1, Equation 4
# "post_match_crd_kk_lin",		#Post Matching Design with CRD assignment					#analyzed via modified OLS of Section 2.3.2, Equation 6
# "post_match_strat_kk_lin",		#Post Matching Design with stratified assignment			#analyzed via modified OLS of Section 2.3.2, Equation 6
#	"post_match_crd_kk_exact",		#Post Matching Design with CRD assignment					#analyzed via permutation test of Section 2.3.3
#	"post_match_strat_kk_exact"		#Post Matching Design with stratified assignment			#analyzed via permutation test of Section 2.3.3

#what do we measure for each set of simulations?
metrics_for_each_run = c(
	"power",
	"res_end_prop_avg",
	"pct_only_matches",
	"pct_only_reservoir",
	"avg_beta_T", 
	"avg_abs_bias",
	"std_err_beta_T",
	"avg_max_std_diff_bal", 
	"avg_max_ks_stat",
	"pct_trt_diff",
	"conv_guessing_strategy_pct_correct_avg",
	"conv_guessing_strategy_pct_correct_sd",
	"match_corr",
	"ssqr_over_sum",
	"ssqr_eq_ssqd_pval",
	"true_var_prop_diff",
	"conv_guessing_strategy_pct_correct_after_n_0_avg",
	"conv_guessing_strategy_pct_correct_after_n_0_sd"
)

#make master results matrix
colnames_master_results = c(
  "algorithm",
  "n",
  "beta_T",
  "betas",
  "rho",
  "cutoff",
  "t_0_pct",
  metrics_for_each_run
)
master_results = data.frame(matrix(NA, nrow = 0, ncol = length(colnames_master_results)))
colnames(master_results) = colnames_master_results

time_began = Sys.time()
num_simulation_settings = length(treatment_effects) * length(ns_to_test) * length(prob_match_cutoff_lambdas) * length(t_0_matching_pcts) * length(all_betas_and_correlations) * length(randomization_types)
num_simulation_setting = 0
run_sims = function(){
  for (sim_type in sim_types){
    source(paste("sim_type_", sim_type, ".R", sep = ""))
  }
}
run_sims()
stopCluster(cl)
master_results
num_atkinson_assignment_errors

############## convergence guessing checking
#source("convergence_guessing_checking.R")

############## severity trend checking
#source("severity_trend_checking.R")



#############checks on distributions
# pacman::p_load(ggplot2)
# fake_normal = rnorm(Nsim_per_block, mean = 1, sd = sd(all_beta_hats[["bcrd_lin"]]))
# ggplot(data.frame(betaThat = all_beta_hats[["bcrd_lin"]], fake_normal = fake_normal)) + 
#   geom_density(aes(x = betaThat), fill = "red", alpha = 0.5) + 
#   geom_density(aes(x = fake_normal), fill = "green", alpha = 0.5)
# ks.test(all_beta_hats[["bcrd_lin"]], fake_normal)

