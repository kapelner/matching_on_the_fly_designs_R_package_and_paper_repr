indic_T_guesses = array(NA, n)
indic_T_guesses[1] = rbinom(1, 1, 0.5)

sum_T = indic_T[1]
for (t in 2 : n){
	t_minus_one_over_two = (t - 1) / 2
	if (sum_T == t_minus_one_over_two){
		#50-50 guess
		indic_T_guesses[t] = rbinom(1, 1, 0.5)
	} else {
		#deterministically guess the least represented tx group
		indic_T_guesses[t] = ifelse(sum_T < t_minus_one_over_two, 1, 0)	
	}
	sum_T = sum_T + indic_T[t]
}
#stop("boom")

