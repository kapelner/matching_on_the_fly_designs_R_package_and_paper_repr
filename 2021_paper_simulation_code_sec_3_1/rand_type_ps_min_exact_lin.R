sys.source("common_ps_min.R", envir = environment())
indic_T_permute_function = function(){rbinom(n, 1, prob_trt)} #wrong since doesn't respect structure
indic_T_permute_function_args = list()
sys.source("common_exact_lin.R", envir = environment())