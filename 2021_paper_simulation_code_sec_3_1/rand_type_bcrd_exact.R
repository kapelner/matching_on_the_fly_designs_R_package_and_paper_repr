sys.source("common_bcrd.R", envir = environment())
indic_T_permute_function = function(){sample(Xy$indic_T)}
indic_T_permute_function_args = list()
sys.source("common_exact.R", envir = environment())