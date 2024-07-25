sys.source("common_atkinson.R", envir = environment())
indic_T_permute_function = function(Xatkinson){atkinson_assignment(Xatkinson)}
indic_T_permute_function_args = list(Xatkinson = x_s)
sys.source("common_exact_lin.R", envir = environment())