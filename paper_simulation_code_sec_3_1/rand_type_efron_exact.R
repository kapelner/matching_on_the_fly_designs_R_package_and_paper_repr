sys.source("common_efron.R", envir = environment())
indic_T_permute_function = function(){efron_biased_coin_design(n, prob_trt)}
indic_T_permute_function_args = list()
sys.source("common_exact.R", envir = environment())