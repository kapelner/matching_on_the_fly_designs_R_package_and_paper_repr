sys.source("common_stratification.R", envir = environment())
indic_T_permute_function = function(){rbinom(n, 1, prob_trt)} #wrong since we're not permuting according to the structure of the allocations
indic_T_permute_function_args = list()
sys.source("common_exact.R", envir = environment())