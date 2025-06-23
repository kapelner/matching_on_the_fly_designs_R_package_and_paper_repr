pacman::p_load(SeqExpMatch, data.table, stringr, dplyr, ggplot2, gridExtra, profvis, doSNOW, tcltk, future.callr, future.apply, progressr)

betaToverall = 1
prob_of_adding_responses = c(0.5, 1)
betaTs = c(0, betaToverall)
ns = c(100)
Nsim = 500

censoring_mechanism = list(q = 0.8, prob = 0.2)
betas = c(3,-3,1,-1,0,0)
nsim_exact_test = 501




response_types = c("continuous", "incidence", "proportion", "count", "survival")
designs = c("CRD", "iBCRD", "Efron", "Atkinson", "KK14", "KK21", "KK21stepwise")
test_types = c("MLE-or-KM-based", "randomization-exact")
# covariate_distributions = c("unif")



linear_mean_function = function(x_vec, w, betaT){
  x_vec %*% betas + (w - 0.5) * betaT
}

Nplot = 6e4
X = matrix(runif(Nplot * 6, -0.3, 0.3), ncol = 6)
mus = data.table()
mus[, w := rbinom(Nplot, 1, 0.5)]
mus[, continuous := as.numeric(X %*% betas + (w - 0.5) * betaToverall)]
mus[, incidence  := exp(continuous) / (1 + exp(continuous))]
mus[, proportion := incidence]
mus[, count      := exp(continuous)]
mus[, survival   := count]

ggplot(mus) +
  geom_histogram(aes(x = count, fill = as.factor(w)), bins = 1e2, alpha = 0.5)
mus[, lapply(.SD, mean), by=w, .SDcols=setdiff(colnames(mus), "w")]

mu_survival_max_to_be_observed = quantile(mus$survival, censoring_mechanism$q)
rm(mus)

survival_k = 0.8
survival_mu_multiple = gamma(1 + 1 / survival_k)
response_functions = list(
  continuous = function(x_vec, w, betaT){
    rnorm(1, linear_mean_function(x_vec, w, betaT), sd = 1)
  },
  incidence = function(x_vec, w, betaT){
    mu = linear_mean_function(x_vec, w, betaT)
    rbinom(1, 1, exp(mu) / (1 + exp(mu)))
  },
  proportion = function(x_vec, w, betaT){
    mu = linear_mean_function(x_vec, w, betaT)
    mu = exp(mu) / (1 + exp(mu))
    sigsq = mu * (1 - mu) / 2
    const = sigsq + mu^2 - mu
    alpha = -mu * const / sigsq
    beta = (mu - 1) * const / sigsq
    rbeta(1, alpha, beta)
  },
  count =      function(x_vec, w, betaT){
    mu = exp(linear_mean_function(x_vec, w, betaT))
    size = 1
    rnbinom(1, size = size, mu = mu)
  },
  survival =   function(x_vec, w, betaT){
    mu = exp(linear_mean_function(x_vec, w, betaT))
    y = rweibull(1, shape = survival_mu_multiple, scale = survival_mu_multiple * mu)
    ifelse(y > mu_survival_max_to_be_observed, mu_survival_max_to_be_observed, y)
  },
  dead = function(y_surv){
    if (y_surv > mu_survival_max_to_be_observed){
      0
    # } else if (runif(1) < censoring_mechanism$prob){
    #   0
    } else {
      1
    }
  }
)


#now let's calculate the betaT's
X100 = X[1:100, ]
rm(X)

Nsum_res_beta_T = 5000
res_beta_T = data.table()
if (!setequal(betaTs, 0)){
  for (nsim_res in 1 : Nsum_res_beta_T){
    res_beta_T = rbind(res_beta_T, data.table(
      # y_cont_T_beta_T_zero = apply(X100, 1, function(xi){
      #   response_functions$continuous(xi, 1, 0)
      # }),
      # y_cont_C_beta_T_zero = apply(X100, 1, function(xi){
      #   response_functions$continuous(xi, 0, 0)
      # }),
      y_cont_T_beta_T_one = apply(X100, 1, function(xi){
        response_functions$continuous(xi, 1, betaToverall)
      }),
      y_cont_C_beta_T_one = apply(X100, 1, function(xi){
        response_functions$continuous(xi, 0, betaToverall)
      }),    
      
      # y_incid_T_beta_T_zero = apply(X100, 1, function(xi){
      #   response_functions$incidence(xi, 1, 0)
      # }),
      # y_incid_C_beta_T_zero = apply(X100, 1, function(xi){
      #   response_functions$incidence(xi, 0, 0)
      # }),
      y_incid_T_beta_T_one = apply(X100, 1, function(xi){
        response_functions$incidence(xi, 1, betaToverall)
      }),
      y_incid_C_beta_T_one = apply(X100, 1, function(xi){
        response_functions$incidence(xi, 0, betaToverall)
      }),  
      
      # y_prop_T_beta_T_zero = apply(X100, 1, function(xi){
      #   response_functions$proportion(xi, 1, 0)
      # }),
      # y_prop_C_beta_T_zero = apply(X100, 1, function(xi){
      #   response_functions$proportion(xi, 0, 0)
      # }),
      y_prop_T_beta_T_one = apply(X100, 1, function(xi){
        response_functions$proportion(xi, 1, betaToverall)
      }),
      y_prop_C_beta_T_one = apply(X100, 1, function(xi){
        response_functions$proportion(xi, 0, betaToverall)
      }), 
      
      # y_count_T_beta_T_zero = apply(X100, 1, function(xi){
      #   response_functions$count(xi, 1, 0)
      # }),
      # y_count_C_beta_T_zero = apply(X100, 1, function(xi){
      #   response_functions$count(xi, 0, 0)
      # }),
      y_count_T_beta_T_one = apply(X100, 1, function(xi){
        response_functions$count(xi, 1, betaToverall)
      }),
      y_count_C_beta_T_one = apply(X100, 1, function(xi){
        response_functions$count(xi, 0, betaToverall)
      }), 
      
      # y_surv_T_beta_T_zero = apply(X100, 1, function(xi){
      #   response_functions$survival(xi, 1, 0)
      # }),
      # y_surv_C_beta_T_zero = apply(X100, 1, function(xi){
      #   response_functions$survival(xi, 0, 0)
      # }),
      y_surv_T_beta_T_one = apply(X100, 1, function(xi){
        response_functions$survival(xi, 1, betaToverall)
      }),
      y_surv_C_beta_T_one = apply(X100, 1, function(xi){
        response_functions$survival(xi, 0, betaToverall)
      })    
    ))
  }
  
  all_means = as.numeric(res_beta_T[, lapply(.SD, mean), .SDcols = colnames(res_beta_T)])
  all_mean_diffs = all_means[seq(1, length(all_means), by = 2)] - all_means[seq(2, length(all_means), by = 2)]
  rm(all_means)
} else {
  all_mean_diffs = rep(NA, 5)
}

#now calculate estimands
estimands_betaT_one = list(
  continuous = list(
    simple_mean_difference = betaToverall,
    KK_compound_mean_difference = betaToverall,
    "continuous_multivariate_regression" = betaToverall,	
    "continuous_KK_compound_multivariate_regression" = betaToverall
    # "continuous_KK_regression_with_covariates_with_matching_dummies" = betaToverall,
    # "continuous_KK_regression_with_covariates_with_random_intercepts" = betaToverall
  ),
  incidence = list(
    simple_mean_difference = all_mean_diffs[2],
    KK_compound_mean_difference = all_mean_diffs[2],
    "incidence_simple_log_odds" = betaToverall,	
    "incidence_multivariate_logistic_regression" = betaToverall,
    "incidence_KK_compound_univariate_logistic_regression" = betaToverall,
    "incidence_KK_compound_multivariate_logistic_regression" = betaToverall
    # "incidence_KK_multivariate_logistic_regression_with_matching_dummies" = betaToverall,	
    # "incidence_KK_multivariate_logistic_regression_with_random_intercepts_for_matches" = betaToverall
  ),
  proportion = list(
    simple_mean_difference = all_mean_diffs[3],
    KK_compound_mean_difference = all_mean_diffs[3],
    "proportion_simple_logodds_regression" = betaToverall,
    "proportion_multivariate_beta_regression" = betaToverall
    #"proportion_KK_compound_univariate_beta_regression",
    #"proportion_KK_compound_multivariate_beta_regression",
    #"proportion_KK_multivariate_beta_regression_with_matching_dummies",
  ),
  count = list(
    simple_mean_difference = all_mean_diffs[4],
    KK_compound_mean_difference = all_mean_diffs[4],
    "count_univariate_negative_binomial_regression" = betaToverall,
    "count_multivariate_negative_binomial_regression" = betaToverall
    #"count_KK_compound_univariate_negative_binomial_regression",	
    #"count_KK_compound_multivariate_negative_binomial_regression",
    # "count_KK_multivariate_negative_binomial_regression_with_matching_dummies" = betaToverall,
    # "count_KK_multivariate_negative_binomial_regression_with_random_intercepts_for_matches" = betaToverall
    
  ),
  survival = list(
    "survival_simple_median_difference" = median(res_beta_T$y_surv_T_beta_T_one) - median(res_beta_T$y_surv_C_beta_T_one),
    "survival_simple_restricted_mean_difference" = all_mean_diffs[5],
    "survival_univariate_weibull_regression" = betaToverall,
    "survival_multivariate_weibull_regression" = betaToverall,
    #"survival_KK_compound_univariate_weibull_regression",	
    #"survival_KK_compound_multivariate_weibull_regression",
    #"survival_KK_multivariate_weibull_regression_with_matching_dummies",
    "survival_univariate_coxph_regression" = survival_k * betaToverall,	
    "survival_multivariate_coxph_regression" = survival_k * betaToverall	
    # "survival_KK_multivariate_coxph_regression_with_matching_dummies" = survival_k * betaToverall,			
    # "survival_KK_multivariate_coxph_regression_with_random_intercepts_for_matches" = survival_k * betaToverall
  )
)

rm(all_mean_diffs, res_beta_T, Nsum_res_beta_T)

X100 = data.table(X100)

exp_settings = data.table(expand.grid(
  n = ns,
  response_type = response_types,
  design = designs,
  test_type = test_types,
  prob_of_adding_responses = prob_of_adding_responses,
  betaT = betaTs
))
#now we kill illegal inference methods
exp_settings = exp_settings[!(!grepl("KK21", design) & prob_of_adding_responses != 0.5), ]
exp_settings = merge(data.frame(nsims = 1 : Nsim), exp_settings, by = NULL)
exp_settings = data.table(exp_settings)
table(exp_settings$prob_of_adding_responses, exp_settings$design)


rm(ns, response_types, designs, test_types, inference_methods, Nsim, Nplot)

