sys.source("common_crd.R", envir = environment())
indic_T_permute_function = function(){rbinom(n, 1, prob_trt)}
indic_T_permute_function_args = list()
sys.source("common_exact.R", envir = environment())