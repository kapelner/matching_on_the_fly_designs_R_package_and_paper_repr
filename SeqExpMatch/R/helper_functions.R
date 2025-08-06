
fast_glm_with_var = function(Xmm, y, glm_function){
	mod = glm_function(Xmm, y)
	XtWX = eigen_Xt_times_diag_w_times_X_cpp(Xmm, mod$w)	
	mod$ssq_b_2 = eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(XtWX, 2)
	mod	
}

negative_binomial_link_cache = MASS::negative.binomial(theta = 1)

fast_negbin_regression = function(Xmm, y){
	mod = glm.fit(Xmm, y, family = negative_binomial_link_cache)
	list(b = mod$coefficients, w = mod$weights)
}

fast_negbin_regression_with_var = function(Xmm, y){
	fast_glm_with_var(Xmm, y, fast_negbin_regression)
}

	
#	loglik <- function(n, th, mu, y, w) sum(w * (lgamma(th + 
#        y) - lgamma(th) - lgamma(y + 1) + th * log(th) + y * 
#        log(mu + (y == 0)) - (th + y) * log(th + mu)))
#    link <- log
#    fam0 <- if (missing(init.theta)) 
#        do.call("poisson", list(link = link))
#    else do.call("negative.binomial", list(theta = init.theta, 
#        link = link))
#    mf <- Call <- match.call()
#    m <- match(c("formula", "data", "subset", "weights", "na.action", 
#        "etastart", "mustart", "offset"), names(mf), 0)
#    mf <- mf[c(1, m)]
#    mf$drop.unused.levels <- TRUE
#    mf[[1L]] <- quote(stats::model.frame)
#    mf <- eval.parent(mf)
#    Terms <- attr(mf, "terms")
#    if (method == "model.frame") 
#        return(mf)
#    Y <- model.response(mf, "numeric")
#    X <- if (!is.empty.model(Terms)) 
#        model.matrix(Terms, mf, contrasts)
#    else matrix(, NROW(Y), 0)
#    w <- model.weights(mf)
#    if (!length(w)) 
#        w <- rep(1, nrow(mf))
#    else if (any(w < 0)) 
#        stop("negative weights not allowed")
#    offset <- model.offset(mf)
#    mustart <- model.extract(mf, "mustart")
#    etastart <- model.extract(mf, "etastart")
#    n <- length(Y)
#    if (!missing(method)) {
#        if (!exists(method, mode = "function")) 
#            stop(gettextf("unimplemented method: %s", sQuote(method)), 
#                domain = NA)
#        glm.fitter <- get(method)
#    }
#    else {
#        method <- "glm.fit"
#        glm.fitter <- stats::glm.fit
#    }
#    if (control$trace > 1) 
#        message("Initial fit:")
#    fit <- glm.fitter(x = X, y = Y, weights = w, start = start, 
#        etastart = etastart, mustart = mustart, offset = offset, 
#        family = fam0, control = list(maxit = control$maxit, 
#            epsilon = control$epsilon, trace = control$trace > 
#                1), intercept = attr(Terms, "intercept") > 0)
#    class(fit) <- c("glm", "lm")
#    mu <- fit$fitted.values
#    th <- as.vector(theta.ml(Y, mu, sum(w), w, limit = control$maxit, 
#        trace = control$trace > 2))
#    if (control$trace > 1) 
#        message(gettextf("Initial value for 'theta': %f", signif(th)), 
#            domain = NA)
#    fam <- do.call("negative.binomial", list(theta = th, link = link))
#    iter <- 0
#    d1 <- sqrt(2 * max(1, fit$df.residual))
#    d2 <- del <- 1
#    g <- fam$linkfun
#    Lm <- loglik(n, th, mu, Y, w)
#    Lm0 <- Lm + 2 * d1
#    while ((iter <- iter + 1) <= control$maxit && (abs(Lm0 - 
#        Lm)/d1 + abs(del)/d2) > control$epsilon) {
#        eta <- g(mu)
#        fit <- glm.fitter(x = X, y = Y, weights = w, etastart = eta, 
#            offset = offset, family = fam, control = list(maxit = control$maxit, 
#                epsilon = control$epsilon, trace = control$trace > 
#                  1), intercept = attr(Terms, "intercept") > 
#                0)
#        t0 <- th
#        th <- theta.ml(Y, mu, sum(w), w, limit = control$maxit, 
#            trace = control$trace > 2)
#        fam <- do.call("negative.binomial", list(theta = th, 
#            link = link))
#        mu <- fit$fitted.values
#        del <- t0 - th
#        Lm0 <- Lm
#        Lm <- loglik(n, th, mu, Y, w)
#        if (control$trace) {
#            Ls <- loglik(n, th, Y, Y, w)
#            Dev <- 2 * (Ls - Lm)
#            message(sprintf("Theta(%d) = %f, 2(Ls - Lm) = %f", 
#                iter, signif(th), signif(Dev)), domain = NA)
#        }
#    }
#    if (!is.null(attr(th, "warn"))) 
#        fit$th.warn <- attr(th, "warn")
#    if (iter > control$maxit) {
#        warning("alternation limit reached")
#        fit$th.warn <- gettext("alternation limit reached")
#    }
#    if (length(offset) && attr(Terms, "intercept")) {
#        null.deviance <- if (length(Terms)) 
#            glm.fitter(X[, "(Intercept)", drop = FALSE], Y, w, 
#                offset = offset, family = fam, control = list(maxit = control$maxit, 
#                  epsilon = control$epsilon, trace = control$trace > 
#                    1), intercept = TRUE)$deviance
#        else fit$deviance
#        fit$null.deviance <- null.deviance
#    }
#    class(fit) <- c("negbin", "glm", "lm")
#    fit$terms <- Terms
#    fit$formula <- as.vector(attr(Terms, "formula"))
#    Call$init.theta <- signif(as.vector(th), 10)
#    Call$link <- link
#    fit$call <- Call
#    if (model) 
#        fit$model <- mf
#    fit$na.action <- attr(mf, "na.action")
#    if (x) 
#        fit$x <- X
#    if (!y) 
#        fit$y <- NULL
#    fit$theta <- as.vector(th)
#    fit$SE.theta <- attr(th, "SE")
#    fit$twologlik <- as.vector(2 * Lm)
#    fit$aic <- -fit$twologlik + 2 * fit$rank + 2
#    fit$contrasts <- attr(X, "contrasts")
#    fit$xlevels <- .getXlevels(Terms, mf)
#    fit$method <- method
#    fit$control <- control
#    fit$offset <- offset
#    fit



binomial_link_cache = binomial()

fast_logistic_regression = function(Xmm, y){
	mod = glm.fit(Xmm, y, family = binomial_link_cache)
	list(b = mod$coefficients, w = mod$weights)
}

fast_logistic_regression_with_var = function(Xmm, y){
	mod = fast_logistic_regression(Xmm, y)
	XtWX = eigen_Xt_times_diag_w_times_X_cpp(Xmm, mod$w)	
	mod$ssq_b_2 = eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(XtWX, 2)
	mod
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


#' Robust Beta Regression
#'
#' Performs Beta regression that if it fails, keeps dropping one column from the model matrix until it works
#'
#' @param form_obj	The formula
#' @param data_obj  The data frame to run beta regression on 
#' @return			The Beta regression model object
#' @export
robust_betareg = function(form_obj, data_obj){
	repeat {
		tryCatch({
			mod = suppressWarnings(betareg::betareg(form_obj, data = data_obj))
			return(mod)
		}, error = function(e){})
		data_obj = data_obj[, 1 : (ncol(data_obj) - 1), drop = FALSE] #chop off one column at a time until it works
		if (ncol(data_obj) == 1){
			break
		}
	}
	NULL
}

##' Robust Negative Binomial Regression
##'
##' Performs Negative Binomial regression that if it fails, keeps dropping one column from the model matrix until it works
##'
##' @param form_obj	The formula
##' @param data_obj  The data frame to run Negative Binomial regression on 
##' @return			The Negative Binomial regression model object
##' @export
#robust_negbinreg = function(form_obj, data_obj){
#	repeat {
#		tryCatch({
#			mod = suppressWarnings(MASS::glm.nb(form_obj, data = data_obj))
#			return(mod)
#		}, error = function(e){})
#		data_obj = data_obj[, 1 : (ncol(data_obj) - 1), drop = FALSE] #chop off one column at a time until it works
#		if (ncol(data_obj) == 0){
#			break
#		}
#	}
#	NA
#}

#' Sample Mode from an array
#'
#' Compute sample mode from an array of numbers
#'
#' @param data	The array of numbers
#' @return		The sample mode
#' @export
sample_mode = function(data){
	as.numeric(names(sort(-table(data)))[1])
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