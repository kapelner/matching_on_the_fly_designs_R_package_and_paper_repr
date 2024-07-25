
if (!exists("atkinson_assignment")){
	num_atkinson_assignment_errors = 0
	atkinson_assignment = function(Xraw){
		X_with_T_and_int = cbind(rbinom(n, 1, prob_trt), 1, Xraw)
#		prob_Ts = array(NA, n)
#		check_probTs = array(NA, n)
		for (t in (p + 2 + 1) : n){
			X_t_min_1 = X_with_T_and_int[1 : (t - 1), ]	
			
			#if everyone is T or everyone is C, Bernoulli-randomize
			if (length(unique(X_t_min_1[, 1])) == 1){
				X_with_T_and_int[t, 1] = rbinom(1, 1, prob_trt)
			#otherwise do Atkinson
			} else {
				
				tryCatch({
					XtX = t(X_t_min_1) %*% X_t_min_1
					
	#				if (Matrix::rankMatrix(XtX) < t - 1){
	#					print(length(unique(X_t_min_1[1, ])))
	#					print(X_t_min_1)
	#				}
					
					M = (t - 1) * solve(XtX)
					
					
					A = M[1, 2 : (p + 2)] %*% c(1, Xraw[t, ]) 
					s_over_A_plus_one_sq = (M[1, 1] / A + 1)^2
					prob_T = s_over_A_plus_one_sq / (s_over_A_plus_one_sq + 1)
	#			prob_Ts[t] = prob_T
					
	#			n1 = sum(X_with_T_and_int[1 : (t-1), 1])
	#			n2 = (t - 1) - n1
	#			check_probTs[t] = n2^2 / (n1^2 + n2^2)
					X_with_T_and_int[t, 1] = rbinom(1, 1, prob_T)	
				}, error = function(e){
					num_atkinson_assignment_errors = num_atkinson_assignment_errors + 1
					X_with_T_and_int[t, 1] = rbinom(1, 1, prob_trt)
				})
			}
			

		}	
		# prob_Ts
		# check_probTs
		# X_with_T_and_int[, 1] #indicT
		
#		list(indicT = X_with_T_and_int[, 1], prob_Ts = prob_Ts, check_probTs = check_probTs)
		X_with_T_and_int[, 1]
	}
}

indic_T = atkinson_assignment(x_s)
sys.source("common_crd_and_bcrd.R", envir = environment(), toplevel.env = environment())