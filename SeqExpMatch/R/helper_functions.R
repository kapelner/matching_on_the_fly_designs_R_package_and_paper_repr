robust_survreg = function(y, dead, cov_matrix_or_vector, dist = "weibull", num_max_iter = 50){
	robust_survreg_with_surv_object(survival::Surv(y, dead), cov_matrix_or_vector, dist = dist, y = y, num_max_iter = num_max_iter)
}

survreg_control_default = survival::survreg.control(
	maxiter = 100,			#default
	rel.tolerance = 1e-9, 	#default
	outer.max = 10			#default
)

robust_survreg_with_surv_object = function(surv_object, cov_matrix_or_vector, dist = "weibull", y = NULL, num_max_iter = 50){
	if (is.null(y)){
		y = as.numeric(surv_object)[1 : length(surv_object)]
	}
	init = lm.fit(cbind(1, cov_matrix_or_vector), log(y))$coefficients
	surv_reg_formula = surv_object ~ .
	cov_matrix_or_vector_data_frame = data.frame(cov_matrix_or_vector)
	num_iter = 1
	repeat {
		mod = suppressWarnings(survival::survreg(
			surv_reg_formula, 
			data = cov_matrix_or_vector_data_frame, 
			dist = dist, 
			init = init,
			control = survreg_control_default
		))
		if (any(is.na(mod$coefficients))){
			if (num_iter == num_max_iter){
				return(mod) #nothing we can do... the survival package is just badly implemented
			}
			#cat(paste("  robust survreg num_iter", num_iter, "init", paste(round(init, 3), collapse = ", "), "\n"))
			init = init + rnorm(length(init))
			num_iter = num_iter + 1
		} else {
			return(mod)
		}
	}
}