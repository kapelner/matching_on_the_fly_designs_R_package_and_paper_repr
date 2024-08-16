robust_survreg = function(y, dead, cov_matrix_or_vector, dist = "weibull", num_max_iter = 50){
	robust_survreg_with_surv_object(survival::Surv(y, dead), cov_matrix_or_vector, dist = dist, y = y, num_max_iter = num_max_iter)
}

survreg_control_default = survival::survreg.control(
	maxiter = 100,			#default
	rel.tolerance = 1e-9, 	#default
	outer.max = 10			#default
)

robust_survreg_with_surv_object = function(surv_object, cov_matrix_or_vector, dist = "weibull", y = NULL, num_max_iter = 50){
#	if (is.null(y)){
#		y = as.numeric(surv_object)[1 : length(surv_object)]
#	}
#	init = lm.fit(cbind(1, cov_matrix_or_vector), log(y))$coefficients
	surv_reg_formula = surv_object ~ .
	cov_matrix_or_vector_data_frame = data.frame(cov_matrix_or_vector)
	init = rep(0, ncol(cov_matrix_or_vector_data_frame) + 1)
	num_iter = 1
	repeat {
		tryCatch({
			mod = suppressWarnings(survival::survreg(
				surv_reg_formula, 
				data = cov_matrix_or_vector_data_frame, 
				dist = dist, 
				init = init,
				control = survreg_control_default
			))
			if (!any(is.na(mod$coefficients))){
				return(mod)
			}	
		}, error = function(e){})
		if (num_iter == num_max_iter){
			break
		}
		#cat(paste("  robust survreg num_iter", num_iter, "init", paste(round(init, 3), collapse = ", "), "\n"))
		init = init + rnorm(length(init))
		num_iter = num_iter + 1
	}
	
	#no more mister nice guy...
	#now we start to eliminate columns that may be causing multicollinearity
	j_killed = c()
	repeat {
		R = cor(cov_matrix_or_vector)
		pairs_to_eliminate_one_of = which(R > .95 & R < 1, arr.ind = TRUE)
		if (nrow(pairs_to_eliminate_one_of) == 0){
			break
		}
		col_indices = c(pairs_to_eliminate_one_of)
		col_indices = col_indices[col_indices > 1] #but we cannot eliminate the first column
		j_kill = sample_mode(col_indices)
		j_killed = c(j_killed, j_kill)
		cov_matrix_or_vector = cov_matrix_or_vector[, -j_kill]
	}
	
	#now we play again
	cov_matrix_or_vector_data_frame = data.frame(cov_matrix_or_vector)
	num_iter = 1
	repeat {
		tryCatch({
			mod = suppressWarnings(survival::survreg(
							surv_reg_formula, 
							data = cov_matrix_or_vector_data_frame, 
							dist = dist, 
							init = init,
							control = survreg_control_default
					))
			if (!any(is.na(mod$coefficients))){
				return(mod)
			}
		}, error = function(e){})
		if (num_iter == num_max_iter){
			stop("could not get survival regression to work!!!") #we give up!!!!
		}
		#cat(paste("  robust survreg num_iter", num_iter, "init", paste(round(init, 3), collapse = ", "), "\n"))
		init = init + rnorm(length(init))
		num_iter = num_iter + 1
	}	
}

robust_betareg = function(form_obj, data_obj){
	repeat {
		tryCatch({
			mod = suppressWarnings(betareg::betareg(form_obj, data = data_obj))
			return(mod)
		}, error = function(e){})
		data_obj = data_obj[, 1 : (ncol(data_obj) - 1)] #chop off one column at a time until it works
		if (ncol(data_obj) == 0){
			break
		}
	}
}

robust_negbinreg = function(form_obj, data_obj){
	repeat {
		tryCatch({
			mod = suppressWarnings(MASS::glm.nb(form_obj, data = data_obj))
			return(mod)
		}, error = function(e){})
		data_obj = data_obj[, 1 : (ncol(data_obj) - 1)] #chop off one column at a time until it works
		if (ncol(data_obj) == 0){
			break
		}
	}
	NA
}

#There's no standard R function for this.
sample_mode = function(data){
	as.numeric(names(sort(-table(data)))[1])
}