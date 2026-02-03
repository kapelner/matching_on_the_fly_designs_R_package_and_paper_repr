#' Logit
#'
#' Calculates the logit i.e., log(p / (1 - p))
#'
#' @param p		The value between 0 and 1 non inclusive
#' @return		Its corresponding logit value as a real number
#' @export
logit = function(p){
	if (any(p <= 0 | p >= 1, na.rm = TRUE)){
		stop("p must be between 0 and 1 non-inclusive")
	}
	log(p / (1 - p))
}

#' Inverse Logit
#'
#' Calculates the inverse logit i.e., 1 / (1 + exp(-x))
#'
#' @param x		Any real number
#' @return		Its corresponding inverse logit value between 0 and 1 non inclusive
#' @export
inv_logit = function(x){
	1 / (1 + exp(-x))
}

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

#' Robust Survival Regression
#'
#' Performs Survival regression from the original response and censoring vectors 
#' that if it fails, keeps trying with a different initialization point until a maximum number of iterations
#'
#' @param y						The response vector
#' @param dead					The censoring vector (1 if dead/uncensored and 0 if censored)
#' @param cov_matrix_or_vector  The model matrix
#' @param dist  				The parametric distribution form (default is Weibull)
#' @param num_max_iter			Maximum # of iterations to repeat (default is 50)
#' @return						The Survival regression model object
#' @export
robust_survreg = function(y, dead, cov_matrix_or_vector, dist = "weibull", num_max_iter = 50){
	robust_survreg_with_surv_object(survival::Surv(y, dead), cov_matrix_or_vector, dist = dist, num_max_iter = num_max_iter)
}

#' Robust Survival Regression
#'
#' Performs Survival regression from a survival::Surv object 
#' that if it fails, keeps trying with a different initialization point until a maximum number of iterations
#'
#' @param surv_object			The survival object (built from the response vector and censoring vector)
#' @param cov_matrix_or_vector  The model matrix
#' @param dist  				The parametric distribution form (default is Weibull)
#' @param num_max_iter			Maximum # of iterations to repeat (default is 50)
#' @return						The Survival regression model object
#' @export
robust_survreg_with_surv_object = function(surv_object, cov_matrix_or_vector, dist = "weibull", num_max_iter = 50){
	surv_reg_formula = surv_object ~ .
	X_mat = as.matrix(cov_matrix_or_vector)
	
	# Eliminate columns that may be causing multicollinearity before attempting model fit
	if (ncol(X_mat) > 1){
		repeat {
			# Use a slightly lower threshold for better robustness
			R = suppressWarnings(stats::cor(X_mat))
			R[is.na(R)] = 1 # Treat constant columns as collinear

			pairs_to_eliminate_one_of = which(abs(R) > .99, arr.ind = TRUE)
			# Filter out diagonal
			pairs_to_eliminate_one_of = pairs_to_eliminate_one_of[pairs_to_eliminate_one_of[,1] < pairs_to_eliminate_one_of[,2], , drop=FALSE]

			if (nrow(pairs_to_eliminate_one_of) == 0){
				break
			}
			# Drop one from each pair, prefer higher index
			j_kill = unique(pairs_to_eliminate_one_of[, 2])
			X_mat = X_mat[, -j_kill, drop = FALSE]
			if (ncol(X_mat) <= 1) break
		}
	}
	
	cov_matrix_or_vector_data_frame = data.frame(X_mat)
	
	# Optimization: Use fast_weibull_regression for initialization if applicable
	if (dist == "weibull") {
		y = surv_object[, 1]
		dead = surv_object[, 2]
		
		# fast_weibull_regression expects X without intercept (it adds it)
		# BUT robust_survreg might already have intercept-like cols? 
		# No, the formula ~ . adds intercept.
		
		res = tryCatch(fast_weibull_regression(y, dead, X_mat), error = function(e) NULL)
		
		if (!is.null(res) && is.finite(res$neg_log_lik)) {
			init_vals = c(res$coefficients, res$log_sigma)
			
			mod = tryCatch({
				suppressWarnings(survival::survreg(
					surv_reg_formula, 
					data = cov_matrix_or_vector_data_frame, 
					dist = dist, 
					init = init_vals,
					control = survival::survreg.control(maxiter = 100, rel.tolerance = 1e-9, outer.max = 10)
				))
			}, error = function(e) NULL)
			
			if (!is.null(mod) && !any(is.na(mod$coefficients))){
				return(mod)
			}
		}
	}

	init = rep(0, ncol(cov_matrix_or_vector_data_frame) + 1)
	num_iter = 1
	repeat {
		tryCatch({
			mod = suppressWarnings(survival::survreg(
				surv_reg_formula, 
				data = cov_matrix_or_vector_data_frame, 
				dist = dist, 
				init = init,
				control = 	survival::survreg.control(
								maxiter = 100,			#default
								rel.tolerance = 1e-9, 	#default
								outer.max = 10			#default
							)
			))
			if (!any(is.na(mod$coefficients))){
				return(mod)
			}	
		}, error = function(e){})
		if (num_iter >= num_max_iter){
			break
		}
		init = init + stats::rnorm(length(init))
		num_iter = num_iter + 1
	}
	
	return(NULL)
}

#' Robust Negative Binomial Regression
#'
#' Performs Negative Binomial regression that if it fails, keeps dropping one column from the model matrix until it works
#'
#' @param form_obj	The formula
#' @param data_obj  The data frame to run Negative Binomial regression on 
#' @return			The Negative Binomial regression model object
#' @export
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

#' Sample Mode from an array
#'
#' Compute sample mode from an array of numbers
#'
#' @param data	The array of numbers
#' @return		The sample mode
#' @export
sample_mode = function(data){
	sample_mode_cpp(data)
}

#' Lean GLM Summary
#'
#' Same as summary.glm except it doesn't calculate the residuals as this takes a long time
#'
#' @param object		The GLM object
#' @param dispersion	See GLM documentation 
#' @param correlation	See GLM documentation 
#' @param symbolic.cor	See GLM documentation 
#' @param ...			Other parameters (currently unused)
#' @return				The summary of the GLM
#' @export
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
	aliased <- is.na(stats::coef(object))
	p <- object$rank
	if (p > 0) {
		p1 <- 1L:p
		Qr <- object$qr
		coef.p <- object$coefficients[Qr$pivot[p1]]
		covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
		dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
		covmat <- dispersion * covmat.unscaled
		var.cf <- diag(covmat)
		s.err <- sqrt(var.cf)
		tvalue <- coef.p/s.err
		dn <- c("Estimate", "Std. Error")
		if (!est.disp) {
			pvalue <- 2 * stats::pnorm(-abs(tvalue))
			coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
			dimnames(coef.table) <- list(names(coef.p), c(dn, 
							"z value", "Pr(>|z|)"))
		}
		else if (df.r > 0) {
			pvalue <- 2 * stats::pt(-abs(tvalue), df.r)
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
		coef.table <- matrix( 0L, 4L)
		dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error", 
						"t value", "Pr(>|t|)"))
		covmat.unscaled <- covmat <- matrix( 0L, 0L)
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

#' Fast Mean Calculation
#'
#' Calculates the mean of a numeric vector using Rcpp for speed.
#'
#' @param x A numeric vector.
#' @return The mean of the vector.
#' @name mean_cpp
#' @rdname mean_cpp
#' @export
#' @examples
#' \dontrun{
#' x <- rnorm(100)
#' mean_cpp(x)
#' }
NULL

#' Fast Variance Calculation
#'
#' Calculates the variance of a numeric vector using Rcpp for speed.
#'
#' @param x A numeric vector.
#' @return The variance of the vector.
#' @name var_cpp
#' @rdname var_cpp
#' @export
#' @examples
#' \dontrun{
#' x <- rnorm(100)
#' var_cpp(x)
#' }
NULL
