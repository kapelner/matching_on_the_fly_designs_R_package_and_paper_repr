assertResponseType = function(response_type, needed_response_type){
	if (response_type != needed_response_type){
		stop("This type of inference is only available for ", needed_response_type, " responses.")
	}
}

assertNoCensoring = function(any_censoring){
	if (any_censoring){
		stop("This type of inference is only available for uncensored responses.")
	}
}

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
	
	if (ncol(cov_matrix_or_vector_data_frame) == 1){ #no tricks will work...
		return(NULL)
	}
	
	#no more mister nice guy...
	#now we start to eliminate columns that may be causing multicollinearity
	j_killed = c()
	repeat {
		R = cor(cov_matrix_or_vector)
		pairs_to_eliminate_one_of = which(abs(R) > .95 & R < 1, arr.ind = TRUE)
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
			return(NULL)
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
		data_obj = data_obj[, 1 : (ncol(data_obj) - 1), drop = FALSE] #chop off one column at a time until it works
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
		data_obj = data_obj[, 1 : (ncol(data_obj) - 1), drop = FALSE] #chop off one column at a time until it works
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

#same as summary.glm except it doesn't calculate the residuals as this takes a long time
summary_glm_lean = function (object, dispersion = NULL, correlation = FALSE, symbolic.cor = FALSE, ...){
	est.disp <- FALSE
	df.r <- object$df.residual
	if (is.null(dispersion)) {
		fam <- object$family
		dispersion <- if (!is.null(fam$dispersion) && !is.na(fam$dispersion)) 
					fam$dispersion
				else if (fam$family %in% c("poisson", "binomial")) 
					1
				else if (df.r > 0) {
					est.disp <- TRUE
					if (any(object$weights == 0)) 
						warning("observations with zero weight not used for calculating dispersion")
					sum((object$weights * object$residuals^2)[object$weights > 
											0])/df.r
				}
				else {
					est.disp <- TRUE
					NaN
				}
	}
	aliased <- is.na(coef(object))
	p <- object$rank
	if (p > 0) {
		p1 <- 1L:p
		Qr <- stats:::qr.lm(object)
		coef.p <- object$coefficients[Qr$pivot[p1]]
		covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
		dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
		covmat <- dispersion * covmat.unscaled
		var.cf <- diag(covmat)
		s.err <- sqrt(var.cf)
		tvalue <- coef.p/s.err
		dn <- c("Estimate", "Std. Error")
		if (!est.disp) {
			pvalue <- 2 * pnorm(-abs(tvalue))
			coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
			dimnames(coef.table) <- list(names(coef.p), c(dn, 
							"z value", "Pr(>|z|)"))
		}
		else if (df.r > 0) {
			pvalue <- 2 * pt(-abs(tvalue), df.r)
			coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
			dimnames(coef.table) <- list(names(coef.p), c(dn, 
							"t value", "Pr(>|t|)"))
		}
		else {
			coef.table <- cbind(coef.p, NaN, NaN, NaN)
			dimnames(coef.table) <- list(names(coef.p), c(dn, 
							"t value", "Pr(>|t|)"))
		}
		df.f <- NCOL(Qr$qr)
	}
	else {
		coef.table <- matrix(, 0L, 4L)
		dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error", 
						"t value", "Pr(>|t|)"))
		covmat.unscaled <- covmat <- matrix(, 0L, 0L)
		df.f <- length(aliased)
	}
	keep <- match(c("call", "terms", "family", "deviance", "aic", 
					"contrasts", "df.residual", "null.deviance", "df.null", 
					"iter", "na.action"), names(object), 0L)
	ans <- c(object[keep], list(
					coefficients = coef.table, 
					aliased = aliased, 
					dispersion = dispersion, 
					df = c(object$rank, df.r, df.f), 
					cov.unscaled = covmat.unscaled, 
					cov.scaled = covmat
				)
			)
	if (correlation && p > 0) {
		dd <- sqrt(diag(covmat.unscaled))
		ans$correlation <- covmat.unscaled/outer(dd, dd)
		ans$symbolic.cor <- symbolic.cor
	}
	class(ans) <- "summary.glm"
	return(ans)
}