# Package installation cache
package_cache = new.env(parent = emptyenv())

#' Check if a package is installed and cache the result
#'
#' @param package_name Character scalar. The name of the package.
#' @return Logical scalar. TRUE if installed, FALSE otherwise.
#' @keywords internal
#' @export
check_package_installed = function(package_name) {
	if (!exists(package_name, envir = package_cache, inherits = FALSE)) {
		# We use the real requireNamespace here once
		res = requireNamespace(package_name, quietly = TRUE)
		assign(package_name, res, envir = package_cache)
	}
	get(package_name, envir = package_cache)
}

# Helper for parallel progress bars
print_progress = function(pb, i, total) {
	if (!is.null(pb)) {
		utils::setTxtProgressBar(pb, i)
	}
}

assert_greedy_experimental_design_installed = function(caller) {
	if (!check_package_installed("GreedyExperimentalDesign")) {
		stop("Package 'GreedyExperimentalDesign' is required for ", caller, ".")
	}
}

assert_optimal_blocks_libraries_installed = function(caller) {
	required_pkgs = c("ompr", "ompr.roi", "ROI.plugin.glpk", "randomizr")
	missing_pkgs = required_pkgs[!vapply(required_pkgs, check_package_installed, logical(1))]
	if (length(missing_pkgs) > 0L) {
		stop("Packages ", paste(missing_pkgs, collapse = ", "), " are required for ", caller, ".")
	}
}

assert_blocktools_installed = function(caller) {
	if (!check_package_installed("blockTools"))
		stop("Package 'blockTools' is required for ", caller, ".")
}

assert_anticlust_installed = function(caller) {
	if (!check_package_installed("anticlust"))
		stop("Package 'anticlust' is required for ", caller, ".")
}

#' Logit
#'
#' Calculates the logit i.e., log(p / (1 - p))
#'
#' @param	p		The value between 0 and 1 non inclusive
#' @param	zero_one_logit_clamp	The clamping amount for exact 0 and 1 values
#' @return	Its corresponding logit value as a real number
#' @examples
#' logit(0.25)
#' @export
logit = function(p, zero_one_logit_clamp = .Machine$double.eps){
	p = pmax(zero_one_logit_clamp, pmin(1 - zero_one_logit_clamp, p))
	log(p / (1 - p))
}

#' Inverse Logit Function
#'
#' Computes the inverse logit of a real number or vector.
#'
#' @param    x               Any real number
#' @param	 zero_one_logit_clamp	The clamping amount
#' @return   Its corresponding inverse logit value between 0 and 1 non inclusive
#' @export
#' @examples
#' inv_logit(c(-1, 0, 1))
inv_logit = function(x, zero_one_logit_clamp = .Machine$double.eps){
	p = 1 / (1 + exp(-x))
	pmax(zero_one_logit_clamp, pmin(1 - zero_one_logit_clamp, p))
}


drop_linearly_dependent_cols = function(M){
	M = as.matrix(M)
	js = seq_len(ncol(M))
	if (ncol(M) > 0){
		# Use a standard tolerance for rank detection
		tol = 1e-7
		rank = matrix_rank_cpp(M, tol = tol)
		if (rank != ncol(M)){
			# Use the same tolerance for R's qr() to be consistent
			qrX = qr(M, tol = tol)
			# QR pivot contains indices of columns in order of their 'independence'
			# but we should trust our matrix_rank_cpp's rank estimate.
			actual_rank = min(rank, qrX$rank)
			js = qrX$pivot[seq_len(actual_rank)]
			M = M[, js, drop = FALSE]
		}
	}
	list(M = M, js = js)
}

drop_highly_correlated_cols = function(M, threshold = 0.99){
	M = as.matrix(M)
	js = seq_len(ncol(M))
	if (ncol(M) <= 1) return(list(M = M, js = js))
	# Drop zero-variance (constant) columns first so they don't produce NA
	# correlations that could accidentally flag the treatment column for removal
	const_cols = which(apply(M, 2, var) == 0)
	if (length(const_cols) > 0) {
		M = M[, -const_cols, drop = FALSE]
		js = js[-const_cols]
		if (ncol(M) <= 1) return(list(M = M, js = js))
	}
	repeat {
		R = suppressWarnings(stats::cor(M))
		pairs = which(abs(R) > threshold, arr.ind = TRUE)
		pairs = pairs[pairs[, 1] < pairs[, 2], , drop = FALSE]
		if (nrow(pairs) == 0) break
		j_kill = unique(pairs[, 2])
		M = M[, -j_kill, drop = FALSE]
		js = js[-j_kill]
		if (ncol(M) <= 1) break
	}
	list(M = M, js = js)
}

assertResponseType = function(response_type, needed_response_type){
	if (!(response_type %in% needed_response_type)){
		stop("This type of inference is only available for ", paste(needed_response_type, collapse = "/"), " responses.")
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
#' that if it fails, keeps trying with a different initialization point until a maximum number
#' of iterations
#'
#' @param	y						The response vector
#' @param	dead					The censoring vector (1 if dead/uncensored and 0 if censored)
#' @param	cov_matrix_or_vector  The model matrix
#' @param	dist  				The parametric distribution form (default is Weibull)
#' @param	num_max_iter			Maximum # of iterations to repeat (default is 50)
#' @return	The Survival regression model object
#' @export
#' @examples
#' y <- c(1.2, 2.4, 1.8, 3.1, 2.7, 4.0)
#' dead <- c(1, 1, 0, 1, 0, 1)
#' x <- data.frame(x1 = c(-1, 0, 1, 0, 1, 2))
#' robust_survreg(y, dead, x)
robust_survreg = function(y, dead, cov_matrix_or_vector, dist = "weibull", num_max_iter = 50){
	robust_survreg_with_surv_object(survival::Surv(y, dead), cov_matrix_or_vector, dist = dist, num_max_iter = num_max_iter)
}

#' Robust Survival Regression
#'
#' Performs Survival regression from a survival::Surv object
#' that if it fails, keeps trying with a different initialization point until a maximum number
#' of iterations
#'
#' @param surv_object                     The survival object (built from the response vector
#'   and censoring vector)
#' @param	cov_matrix_or_vector  The model matrix
#' @param	dist  				The parametric distribution form (default is Weibull)
#' @param	num_max_iter			Maximum # of iterations to repeat (default is 50)
#' @return	The Survival regression model object
#' @examples
#' surv_obj <- survival::Surv(c(1.2, 2.4, 1.8, 3.1, 2.7, 4.0), c(1, 1, 0, 1, 0, 1))
#' x <- data.frame(x1 = c(-1, 0, 1, 0, 1, 2))
#' robust_survreg_with_surv_object(surv_obj, x)
#' @export
robust_survreg_with_surv_object = function(surv_object, cov_matrix_or_vector, dist = "weibull", num_max_iter = 50){
	surv_reg_formula = surv_object ~ .
	X_mat = as.matrix(cov_matrix_or_vector)

	# Eliminate columns that may be causing multicollinearity before attempting model fit
	X_mat = drop_highly_correlated_cols(X_mat)$M
	X_mat = drop_linearly_dependent_cols(X_mat)$M

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
#' Performs Negative Binomial regression that if it fails, keeps dropping one column from the
#' model matrix until it works
#'
#' @param	form_obj	The formula
#' @param	data_obj  The data frame to run Negative Binomial regression on
#' @return	The Negative Binomial regression model object
#' @export
#' @examples
#' dat <- data.frame(y = c(0, 1, 1, 2, 3, 4), x1 = c(-1, 0, 1, 0, 1, 2))
#' robust_negbinreg(y ~ x1, dat)
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
#' @param	data	The array of numbers
#' @return	The sample mode
#' @examples
#' sample_mode(c(1, 2, 2, 3, 3, 3))
#' @export
sample_mode = function(data){
	sample_mode_cpp(data)
}

#' Lean GLM Summary
#'
#' Same as summary.glm except it doesn't calculate the residuals as this takes a long time
#'
#' @param	object		The GLM object
#' @param	dispersion	See GLM documentation
#' @param	correlation	See GLM documentation
#' @param	symbolic.cor	See GLM documentation
#' @param	...			Other parameters (currently unused)
#' @return	The summary of the GLM
#' @examples
#' mod <- glm(c(0, 1, 0, 1, 1, 0) ~ c(-1, 0, 1, 0, 1, 2), family = binomial())
#' summary_glm_lean(mod)$
#'   coefficients
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
#' @param	x A numeric vector.
#' @return	The mean of the vector.
#' @name mean_cpp
#' @rdname mean_cpp
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
#' @param	x A numeric vector.
#' @return	The variance of the vector.
#' @name var_cpp
#' @rdname var_cpp
#' @examples
#' \dontrun{
#' x <- rnorm(100)
#' var_cpp(x)
#' }
NULL

.compute_kk_basic_match_data = function(X, n, y, w, m_vec){
	if (is.null(m_vec)){
		m_vec = rep(NA_integer_, n)
	}
	m_vec[is.na(m_vec)] = 0
	compute_zhang_match_data_cpp(w, m_vec, y, X)
}

# Cached variant: the X/m structural part (X_matched_diffs, X_matched_diffs_full,
# X_reservoir, m) depends only on m_vec and X, never on y or w.
#
# Cache hierarchy:
#   1. des_priv (design object's private env) — shared across ALL inference objects
#      on the same design.  Used when the current m_vec equals the design's own m_vec,
#      i.e., for the original inference and every randomization iteration (m_vec fixed,
#      only y/w permuted).
#   2. private_env (inference object's private env) — local fallback for bootstrap
#      resamples, which have a different m_vec_b and must not corrupt the design cache.
#
# Data is NEVER written to the global environment.
.compute_kk_basic_match_data_cached = function(private_env, des_priv, X, n, y, w, m_vec){
	if (is.null(m_vec)) m_vec = rep(NA_integer_, n)
	m_vec[is.na(m_vec)] = 0L

	# --- Fast path 1: design-level structural cache ---
	if (!is.null(des_priv) &&
	    !is.null(des_priv$kk_xm_structural) &&
	    identical(m_vec, des_priv$kk_xm_m_vec)){
		wy = compute_kk_wy_stats_cpp(as.numeric(y), as.integer(w), as.integer(m_vec))
		return(c(des_priv$kk_xm_structural, wy))
	}

	# --- Fast path 2: inference-level structural cache (bootstrap case) ---
	if (!is.null(private_env$kk_xm_structural) &&
	    identical(m_vec, private_env$kk_xm_m_vec)){
		wy = compute_kk_wy_stats_cpp(as.numeric(y), as.integer(w), as.integer(m_vec))
		return(c(private_env$kk_xm_structural, wy))
	}

	# --- Full computation ---
	full = compute_zhang_match_data_cpp(as.integer(w), as.integer(m_vec), as.numeric(y), as.matrix(X))
	structural = full[c("m", "X_matched_diffs", "X_matched_diffs_full", "X_reservoir")]

	# Store in the design when the current m_vec is the design's own m_vec (so all
	# inference objects on this design share it); otherwise store locally (bootstrap).
	if (!is.null(des_priv)){
		des_m = des_priv$m
		if (is.null(des_m)) des_m = rep(NA_integer_, n)
		des_m[is.na(des_m)] = 0L
		if (identical(m_vec, des_m)){
			des_priv$kk_xm_structural = structural
			des_priv$kk_xm_m_vec      = m_vec
		} else {
			private_env$kk_xm_structural = structural
			private_env$kk_xm_m_vec      = m_vec
		}
	} else {
		private_env$kk_xm_structural = structural
		private_env$kk_xm_m_vec      = m_vec
	}
	full
}

.compute_kk_lin_basic_match_data = function(X, n, y, w, m_vec){
	if (is.null(m_vec)){
		m_vec = rep(NA_integer_, n)
	}
	m_vec[is.na(m_vec)] = 0
	compute_kk_lin_match_data_cpp(w, m_vec, y, X)
}

# Cached variant for the lin (means + diffs) C++ path. Same hierarchy as above.
.compute_kk_lin_basic_match_data_cached = function(private_env, des_priv, X, n, y, w, m_vec){
	if (is.null(m_vec)) m_vec = rep(NA_integer_, n)
	m_vec[is.na(m_vec)] = 0L

	if (!is.null(des_priv) &&
	    !is.null(des_priv$kk_lin_xm_structural) &&
	    identical(m_vec, des_priv$kk_lin_xm_m_vec)){
		wy = compute_kk_lin_wy_stats_cpp(as.numeric(y), as.integer(w), as.integer(m_vec))
		return(c(des_priv$kk_lin_xm_structural, wy))
	}

	if (!is.null(private_env$kk_lin_xm_structural) &&
	    identical(m_vec, private_env$kk_lin_xm_m_vec)){
		wy = compute_kk_lin_wy_stats_cpp(as.numeric(y), as.integer(w), as.integer(m_vec))
		return(c(private_env$kk_lin_xm_structural, wy))
	}

	full = compute_kk_lin_match_data_cpp(as.integer(w), as.integer(m_vec), as.numeric(y), as.matrix(X))
	structural = full[c("m", "X_matched_diffs_full", "X_matched_means_full", "X_reservoir")]

	if (!is.null(des_priv)){
		des_m = des_priv$m
		if (is.null(des_m)) des_m = rep(NA_integer_, n)
		des_m[is.na(des_m)] = 0L
		if (identical(m_vec, des_m)){
			des_priv$kk_lin_xm_structural = structural
			des_priv$kk_lin_xm_m_vec      = m_vec
		} else {
			private_env$kk_lin_xm_structural = structural
			private_env$kk_lin_xm_m_vec      = m_vec
		}
	} else {
		private_env$kk_lin_xm_structural = structural
		private_env$kk_lin_xm_m_vec      = m_vec
	}
	full
}

# Computes and caches the structural bootstrap components (i_reservoir, pair_rows, n_reservoir)
# for designs with a match vector (KK14, FixedBinaryMatch). Idempotent: no-op if already cached.
.init_kk_bootstrap_structure = function(des_priv){
	if (!is.null(des_priv$kk_boot_pair_rows)) return(invisible(NULL))
	m_vec = des_priv$m
	n = des_priv$n
	if (is.null(m_vec)){
		des_priv$kk_boot_i_reservoir  = seq_len(n)
		des_priv$kk_boot_n_reservoir  = n
		des_priv$kk_boot_pair_rows    = matrix(integer(0), nrow = 0L, ncol = 2L)
		return(invisible(NULL))
	}
	m_vec_int = as.integer(m_vec)
	m_vec_int[is.na(m_vec_int)] = 0L
	i_reservoir = which(m_vec_int == 0L)
	m_max = max(m_vec_int)
	pair_rows = if (m_max > 0L) {
		pr = matrix(integer(0), nrow = m_max, ncol = 2L)
		for (pid in seq_len(m_max)) pr[pid, ] = which(m_vec_int == pid)
		pr
	} else {
		matrix(integer(0), nrow = 0L, ncol = 2L)
	}
	des_priv$kk_boot_i_reservoir = i_reservoir
	des_priv$kk_boot_n_reservoir = length(i_reservoir)
	des_priv$kk_boot_pair_rows   = pair_rows
	invisible(NULL)
}

# Draws a KK-aware bootstrap sample: resamples reservoir subjects iid and matched pairs as units.
# Returns list(i_b, m_vec_b) compatible with bootstrap_subset_inference.
.draw_kk_bootstrap_indices = function(des_priv){
	.init_kk_bootstrap_structure(des_priv)
	draw_kk_bootstrap_sample_cpp(
		i_reservoir  = des_priv$kk_boot_i_reservoir,
		pair_rows    = des_priv$kk_boot_pair_rows,
		n_reservoir  = des_priv$kk_boot_n_reservoir
	)
}

.extract_se_from_rq_fit = function(fit, coef_name){
	is_bad_se = function(x) !is.finite(x) || x <= 0 || x > 1e6

	se = tryCatch({
		s_fit = suppressWarnings(summary(fit, se = "nid", covariance = TRUE))
		ct = s_fit$coefficients
		if (coef_name %in% rownames(ct)) ct[coef_name, "Std. Error"] else NA_real_
	}, error = function(e) NA_real_)

	if (is_bad_se(se)){
		se = tryCatch({
			s_fit = suppressWarnings(summary(fit, se = "iid", covariance = TRUE))
			ct = s_fit$coefficients
			if (coef_name %in% rownames(ct)) ct[coef_name, "Std. Error"] else NA_real_
		}, error = function(e) NA_real_)
	}

	if (is_bad_se(se)) NA_real_ else se
}

.complete_pair_index_matrix = function(pair_id){
	pair_id = as.integer(pair_id)
	valid = !is.na(pair_id) & pair_id > 0L
	if (!any(valid)) return(matrix(integer(0), ncol = 2))
	pair_rows = split(which(valid), pair_id[valid])
	pair_rows = pair_rows[lengths(pair_rows) == 2L]
	if (length(pair_rows) == 0L) return(matrix(integer(0), ncol = 2))
	pair_mat = do.call(rbind, lapply(pair_rows, function(idx) sort(as.integer(idx))))
	storage.mode(pair_mat) = "integer"
	pair_mat
}

.weibull_aft_margin_terms = function(y, eta, sigma){
	y = pmax(as.numeric(y), .Machine$double.xmin)
	log_t = log(y)
	log_H = (log_t - eta) / sigma
	H = exp(pmin(log_H, 700))
	log_f = log_H - log(sigma) - log_t - H
	list(H = H, log_f = log_f)
}

.clayton_copula_logA = function(H1, H2, theta){
	h1 = theta * H1
	h2 = theta * H2
	m = pmax(h1, h2, 0)
	inner = exp(h1 - m) + exp(h2 - m) - exp(-m)
	inner = pmax(inner, .Machine$double.xmin)
	m + log(inner)
}

.extract_survreg_start = function(y, dead, Xmm){
	full_names = c("(Intercept)", colnames(Xmm))
	start_beta = stats::setNames(rep(0, length(full_names)), full_names)
	start_log_sigma = 0

	mod_fast = tryCatch(fast_weibull_regression(y, dead, Xmm), error = function(e) NULL)
	if (!is.null(mod_fast) && !is.null(mod_fast$coefficients)){
		common = intersect(names(start_beta), names(mod_fast$coefficients))
		start_beta[common] = mod_fast$coefficients[common]
		if (!is.null(mod_fast$log_sigma) && is.finite(mod_fast$log_sigma)){
			start_log_sigma = mod_fast$log_sigma
		}
		return(list(beta = start_beta, log_sigma = start_log_sigma))
	}

	mod = robust_survreg_with_surv_object(survival::Surv(y, dead), Xmm)
	if (is.null(mod)) return(list(beta = start_beta, log_sigma = start_log_sigma))

	mod_coef = c(mod$coefficients, "log(scale)" = log(mod$scale))
	common = intersect(names(start_beta), names(mod_coef))
	start_beta[common] = mod_coef[common]
	if (is.finite(mod_coef["log(scale)"])){
		start_log_sigma = mod_coef["log(scale)"]
	}
	list(beta = start_beta, log_sigma = start_log_sigma)
}

.fit_standard_weibull_aft_from_matrix = function(y, dead, Xmm, estimate_only = FALSE){
	if (length(y) == 0L || sum(dead) == 0L) return(NULL)
	mod_fast = tryCatch(fast_weibull_regression(y, dead, Xmm), error = function(e) NULL)
	if (!is.null(mod_fast) &&
	    !is.null(mod_fast$coefficients) &&
	    (isTRUE(estimate_only) || !is.null(mod_fast$vcov)) &&
	    "w" %in% names(mod_fast$coefficients) &&
	    (isTRUE(estimate_only) || "w" %in% rownames(mod_fast$vcov))){
		beta = as.numeric(mod_fast$coefficients["w"])
		ssq = if (isTRUE(estimate_only)) NA_real_ else as.numeric(mod_fast$vcov["w", "w"])
		if (is.finite(beta) && (isTRUE(estimate_only) || (is.finite(ssq) && ssq > 0))){
			return(list(beta = beta, ssq = ssq))
		}
	}

	mod = robust_survreg_with_surv_object(survival::Surv(y, dead), Xmm)
	if (is.null(mod) || is.null(mod$coefficients) || is.null(mod$var)) return(NULL)

	mod_coef = c(mod$coefficients, "log(scale)" = log(mod$scale))
	mod_vcov = mod$var
	coef_names = c(names(mod$coefficients), "log(scale)")
	colnames(mod_vcov) = rownames(mod_vcov) = coef_names
	if (!("w" %in% names(mod_coef)) || !("w" %in% rownames(mod_vcov))) return(NULL)

	beta = as.numeric(mod_coef["w"])
	ssq = as.numeric(mod_vcov["w", "w"])
	if (!is.finite(beta) || !is.finite(ssq) || ssq <= 0) return(NULL)
	list(beta = beta, ssq = ssq)
}

.extract_lognormal_start = function(y, dead, Xmm, event_indicator){
	full_names = c("(Intercept)", colnames(Xmm))
	start_beta = stats::setNames(rep(0, length(full_names)), full_names)
	start_log_sigma = 0

	mod = robust_survreg_with_surv_object(
		survival::Surv(y, event_indicator),
		Xmm,
		dist = "lognormal"
	)
	if (is.null(mod)) return(list(beta = start_beta, log_sigma = start_log_sigma))

	mod_coef = c(mod$coefficients, "log(scale)" = log(mod$scale))
	common = intersect(names(start_beta), names(mod_coef))
	start_beta[common] = mod_coef[common]
	if (is.finite(mod_coef["log(scale)"])){
		start_log_sigma = mod_coef["log(scale)"]
	}
	list(beta = start_beta, log_sigma = start_log_sigma)
}

.fit_dep_cens_transform_model = function(y, dead, Xmm, estimate_only = FALSE){
	y = pmax(as.numeric(y), .Machine$double.xmin)
	dead = as.integer(dead > 0)
	Xmm = as.matrix(Xmm)
	if (length(y) != nrow(Xmm) || length(dead) != nrow(Xmm)){
		stop("Dependent censoring transformation fit inputs must have matching row counts.")
	}
	if (sum(dead) == 0L || sum(1L - dead) == 0L) return(NULL)

	if (is.null(colnames(Xmm))){
		full_names = c("treatment", paste0("x", seq_len(max(ncol(Xmm) - 1L, 0L))))
		colnames(Xmm) = full_names[seq_len(ncol(Xmm))]
	}

	Xfull = cbind("(Intercept)" = 1, Xmm)
	num_beta = ncol(Xfull)
	log_y = log(y)

	start_event = .extract_lognormal_start(y, dead, Xmm, dead)
	start_cens = .extract_lognormal_start(y, dead, Xmm, 1L - dead)
	base_start = c(
		unname(start_event$beta),
		unname(start_cens$beta),
		start_event$log_sigma,
		start_cens$log_sigma
	)

	neg_loglik = function(par){
		beta_event = par[seq_len(num_beta)]
		beta_cens = par[num_beta + seq_len(num_beta)]
		log_sigma_event = par[2L * num_beta + 1L]
		log_sigma_cens = par[2L * num_beta + 2L]
		atanh_rho = par[2L * num_beta + 3L]

		if (!is.finite(log_sigma_event) || !is.finite(log_sigma_cens) ||
		    log_sigma_event < -8 || log_sigma_event > 8 ||
		    log_sigma_cens < -8 || log_sigma_cens > 8 ||
		    !is.finite(atanh_rho) || abs(atanh_rho) > 8){
			return(1e100)
		}

		sigma_event = exp(log_sigma_event)
		sigma_cens = exp(log_sigma_cens)
		rho = tanh(atanh_rho)
		one_minus_rho_sq = pmax(1 - rho^2, .Machine$double.eps)
		sd_cond = sqrt(one_minus_rho_sq)

		mu_event = as.vector(Xfull %*% beta_event)
		mu_cens = as.vector(Xfull %*% beta_cens)
		z_event = (log_y - mu_event) / sigma_event
		z_cens = (log_y - mu_cens) / sigma_cens

		log_f_event = stats::dnorm(z_event, log = TRUE) - log_sigma_event - log_y
		log_f_cens = stats::dnorm(z_cens, log = TRUE) - log_sigma_cens - log_y
		log_surv_cens_cond = stats::pnorm((rho * z_event - z_cens) / sd_cond, log.p = TRUE)
		log_surv_event_cond = stats::pnorm((rho * z_cens - z_event) / sd_cond, log.p = TRUE)

		loglik = dead * (log_f_event + log_surv_cens_cond) +
			(1 - dead) * (log_f_cens + log_surv_event_cond)
		if (any(!is.finite(loglik))) return(1e100)
		-sum(loglik)
	}

	starts = if (isTRUE(estimate_only)) {
		# For bootstrap iterations, use only one start (no correlation) for speed.
		list(c(base_start, 0))
	} else {
		list(
			c(base_start, 0),
			c(base_start, atanh(0.25)),
			c(base_start, atanh(-0.25))
		)
	}
	best = NULL
	for (start_par in starts){
		fit = tryCatch(
			stats::optim(
				par = start_par,
				fn = neg_loglik,
				method = "BFGS",
				hessian = !isTRUE(estimate_only),
				control = list(
					maxit = 2000, 
					# Tight tolerance for main estimate, looser for bootstrap iterations
					reltol = if (isTRUE(estimate_only)) 1e-7 else 1e-9
				)
			),
			error = function(e) NULL
		)
		if (is.null(fit) || !is.finite(fit$value)) next
		if (is.null(best) || fit$value < best$value) best = fit
	}
	if (is.null(best)) return(NULL)

	if (isTRUE(estimate_only)) {
		coefficients = best$par
		event_names = colnames(Xfull)
		cens_names = paste0("censoring_", event_names)
		param_names = c(event_names, cens_names, "log_scale_event", "log_scale_censoring", "atanh_rho")
		names(coefficients) = param_names
		return(list(coefficients = coefficients, vcov = NULL))
	}

	hess = best$hessian
	vcov_full = tryCatch(solve(hess), error = function(e) NULL)
	if (is.null(vcov_full) || any(!is.finite(diag(vcov_full)))) return(NULL)

	event_names = colnames(Xfull)
	cens_names = paste0("censoring_", event_names)
	param_names = c(event_names, cens_names, "log_scale_event", "log_scale_censoring", "atanh_rho")
	rownames(vcov_full) = colnames(vcov_full) = param_names

	coefficients = best$par
	names(coefficients) = param_names
	list(coefficients = coefficients, vcov = vcov_full)
}



.sanitize_proportion_response = function(y, interior = FALSE){
	assertNumeric(y, any.missing = FALSE)
	y = as.numeric(y)
	if (length(y) == 0L) return(y)
	if (isTRUE(interior)) {
		eps = .Machine$double.eps
		return(pmin(1 - eps, pmax(eps, y)))
	}
	pmin(1, pmax(0, y))
}
.softmax_three_from_logits = function(alpha0, alpha1){
	m = max(0, alpha0, alpha1)
	e0 = exp(alpha0 - m)
	e1 = exp(alpha1 - m)
	e2 = exp(-m)
	den = e0 + e1 + e2
	c(pi0 = e0 / den, pi1 = e1 / den, pib = e2 / den)
}

.build_zoib_start = function(y, Xmm){
	y = as.numeric(y)
	eps = .Machine$double.eps
	y_clip = pmin(pmax(y, eps), 1 - eps)
	beta_start = rep(0, ncol(Xmm) + 1L)
	names(beta_start) = c("(Intercept)", colnames(Xmm))

	glm_start = tryCatch(
		fast_logistic_regression_cpp(
			X = cbind(1, Xmm),
			y = y_clip
		),
		error = function(e) NULL
	)
	if (!is.null(glm_start) && length(glm_start$b) == length(beta_start)){
		if (all(is.finite(glm_start$b))){
			beta_start = as.numeric(glm_start$b)
			names(beta_start) = c("(Intercept)", colnames(Xmm))
		}
	}

	pi0 = mean(y == 0)
	pi1 = mean(y == 1)
	pib = max(1 - pi0 - pi1, 1e-4)
	pi0 = min(max(pi0, 1e-4), 1 - 2e-4)
	pi1 = min(max(pi1, 1e-4), 1 - pi0 - 1e-4)
	pib = max(1 - pi0 - pi1, 1e-4)

	c(
		unname(beta_start),
		log_phi = log(10),
		alpha0 = log(pi0 / pib),
		alpha1 = log(pi1 / pib)
	)
}

.fit_zero_one_inflated_beta = function(y, Xmm, estimate_only = FALSE){
	y = as.numeric(y)
	Xmm = as.matrix(Xmm)
	if (length(y) != nrow(Xmm)){
		stop("Zero/one-inflated beta fit inputs must have matching row counts.")
	}
	if (!all(is.finite(y)) || any(y < 0 | y > 1)){
		stop("Zero/one-inflated beta requires y in [0, 1].")
	}
	if (sum(y > 0 & y < 1) == 0L) return(NULL)

	if (is.null(colnames(Xmm))){
		full_names = c("treatment", paste0("x", seq_len(max(ncol(Xmm) - 1L, 0L))))
		colnames(Xmm) = full_names[seq_len(ncol(Xmm))]
	}

	Xfull = cbind("(Intercept)" = 1, Xmm)
	p = ncol(Xfull)
	is_zero = y == 0
	is_one = y == 1
	is_beta = !(is_zero | is_one)
	y_beta = y[is_beta]
	X_beta = Xfull[is_beta, , drop = FALSE]

	neg_loglik = function(par){
		beta = par[seq_len(p)]
		log_phi = par[p + 1L]
		alpha0 = par[p + 2L]
		alpha1 = par[p + 3L]

		if (!is.finite(log_phi) || log_phi < -12 || log_phi > 12 ||
		    !is.finite(alpha0) || abs(alpha0) > 25 ||
		    !is.finite(alpha1) || abs(alpha1) > 25){
			return(1e100)
		}

		phi = exp(log_phi)
		eta = as.vector(Xfull %*% beta)
		mu = pmin(pmax(inv_logit(eta), 1e-8), 1 - 1e-8)
		pis = .softmax_three_from_logits(alpha0, alpha1)

		ll = numeric(length(y))
		ll[is_zero] = log(pis["pi0"])
		ll[is_one] = log(pis["pi1"])
		if (any(is_beta)){
			mu_beta = mu[is_beta]
			shape1 = pmax(mu_beta * phi, 1e-8)
			shape2 = pmax((1 - mu_beta) * phi, 1e-8)
			ll[is_beta] = log(pis["pib"]) + stats::dbeta(y_beta, shape1 = shape1, shape2 = shape2, log = TRUE)
		}
		if (any(!is.finite(ll))) return(1e100)
		-sum(ll)
	}

	start0 = .build_zoib_start(y, Xmm)
	starts = list(start0)

	best = NULL
	best_val = Inf
	for (start_par in starts){
		fit = tryCatch(
			fast_zero_one_inflated_beta_cpp(Xfull, y, start_par),
			error = function(e) NULL
		)
		if (is.null(fit) || !is.finite(fit$neg_loglik)) next
		if (is.null(best) || fit$neg_loglik < best_val){
			best = fit
			best_val = fit$neg_loglik
		}
	}
	if (is.null(best)) return(NULL)

	best_params = as.numeric(best$coefficients)
	param_names = c(colnames(Xfull), "log_phi", "alpha0", "alpha1")
	coef_full = best_params
	names(coef_full) = param_names

	if (estimate_only) {
		return(list(
			coefficients = coef_full,
			vcov = NULL
		))
	}

	vcov_full = best$vcov
	if (!is.matrix(vcov_full) || any(dim(vcov_full) != length(param_names))){
		vcov_full = NULL
	}
	if (is.null(vcov_full) || any(!is.finite(diag(vcov_full)))){
		vcov_full = tryCatch(numDeriv::hessian(neg_loglik, best_params), error = function(e) NULL)
		vcov_full = tryCatch(solve(vcov_full), error = function(e) NULL)
	}
	if (is.null(vcov_full) || any(!is.finite(diag(vcov_full)))){
		hess_alt = tryCatch(numDeriv::hessian(neg_loglik, best_params), error = function(e) NULL)
		if (!is.null(hess_alt)){
			vcov_full = tryCatch(MASS::ginv(hess_alt), error = function(e) NULL)
		}
	}
	if (is.null(vcov_full) || any(!is.finite(diag(vcov_full)))){
		# Keep the MLE when curvature is too unstable for a usable covariance matrix.
		vcov_full = matrix(NA_real_, nrow = length(param_names), ncol = length(param_names))
	}
	rownames(vcov_full) = colnames(vcov_full) = param_names

	list(
		coefficients = coef_full,
		vcov = vcov_full
	)
}

.fit_clayton_weibull_aft = function(y, dead, Xmm, pair_id, include_singletons = FALSE, starts = NULL, estimate_only = FALSE){
	y = as.numeric(y)
	dead = as.integer(dead > 0)
	Xmm = as.matrix(Xmm)
	if (is.null(colnames(Xmm))){
		full_names = c("w", paste0("x", seq_len(max(ncol(Xmm) - 1L, 0L))))
		colnames(Xmm) = full_names[seq_len(ncol(Xmm))]
	}
	if (length(y) != nrow(Xmm) || length(dead) != nrow(Xmm) || length(pair_id) != nrow(Xmm)){
		stop("Clayton copula fit inputs must have matching row counts.")
	}

	pair_idx = .complete_pair_index_matrix(pair_id)
	if (nrow(pair_idx) == 0L && !include_singletons) return(NULL)

	pair_rows = if (nrow(pair_idx) > 0L) sort(unique(as.vector(pair_idx))) else integer(0)
	singleton_rows = if (include_singletons) setdiff(seq_len(nrow(Xmm)), pair_rows) else integer(0)
	rows_used = sort(unique(c(pair_rows, singleton_rows)))
	if (length(rows_used) == 0L || sum(dead[rows_used]) == 0L) return(NULL)

	Xfull = cbind("(Intercept)" = 1, Xmm)
	num_beta = ncol(Xfull)

	# Pre-extract constant indices and values for the likelihood function to avoid overhead
	has_pairs = nrow(pair_idx) > 0L
	if (has_pairs){
		i1 = pair_idx[, 1]
		i2 = pair_idx[, 2]
		d1 = dead[i1]
		d2 = dead[i2]
		mask00 = d1 == 0L & d2 == 0L
		mask10 = d1 == 1L & d2 == 0L
		mask01 = d1 == 0L & d2 == 1L
		mask11 = d1 == 1L & d2 == 1L
	}
	has_singletons = length(singleton_rows) > 0L
	if (has_singletons){
		d_sg = dead[singleton_rows]
		d_sg_comp = 1 - d_sg
	}

	neg_loglik = function(par){
		log_sigma = par[num_beta + 1L]
		log_theta = par[num_beta + 2L]
		if (!is.finite(log_sigma) || !is.finite(log_theta) ||
		    log_sigma < -8 || log_sigma > 8 || log_theta < -12 || log_theta > 6){
			return(1e100)
		}
		sigma = exp(log_sigma)
		theta = exp(log_theta)
		eta = as.vector(Xfull %*% par[seq_len(num_beta)])
		margin_terms = .weibull_aft_margin_terms(y, eta, sigma)
		H = margin_terms$H
		log_f = margin_terms$log_f

		loglik = 0
		if (has_pairs){
			H1 = H[i1]
			H2 = H[i2]
			logf1 = log_f[i1]
			logf2 = log_f[i2]
			logA = .clayton_copula_logA(H1, H2, theta)
			pair_ll = numeric(length(i1))

			pair_ll[mask00] = -(1 / theta) * logA[mask00]
			pair_ll[mask10] = logf1[mask10] + (-1 / theta - 1) * logA[mask10] + (theta + 1) * H1[mask10]
			pair_ll[mask01] = logf2[mask01] + (-1 / theta - 1) * logA[mask01] + (theta + 1) * H2[mask01]
			pair_ll[mask11] = log(theta + 1) + logf1[mask11] + logf2[mask11] +
				(-1 / theta - 2) * logA[mask11] + (theta + 1) * (H1[mask11] + H2[mask11])

			if (any(!is.finite(pair_ll))) return(1e100)
			loglik = loglik + sum(pair_ll)
		}

		if (has_singletons){
			sg_ll = d_sg * log_f[singleton_rows] - d_sg_comp * H[singleton_rows]
			if (any(!is.finite(sg_ll))) return(1e100)
			loglik = loglik + sum(sg_ll)
		}

		if (!is.finite(loglik)) return(1e100)
		-loglik
	}

	if (is.null(starts)){
		start = .extract_survreg_start(y[rows_used], dead[rows_used], Xmm[rows_used, , drop = FALSE])
		start_par_base = c(unname(start$beta), start$log_sigma)
		starts = list(
			c(start_par_base, log(0.10)),
			c(start_par_base, log(0.50)),
			c(start_par_base, log(1.50))
		)
	}

	best = NULL
	# Use a slightly coarser tolerance for randomization draws if nsim is high
	control_list = list(maxit = 2000, reltol = 1e-9)

	for (start_par in starts){
		fit = tryCatch(
			stats::optim(
				par = start_par,
				fn = neg_loglik,
				method = "BFGS",
				hessian = !isTRUE(estimate_only),
				control = control_list
			),
			error = function(e) NULL
		)
		if (is.null(fit) || !is.finite(fit$value)) next
		if (is.null(best) || fit$value < best$value){
			best = fit
		}
	}
	if (is.null(best)) return(NULL)

	beta_hat = best$par[seq_len(num_beta)]
	names(beta_hat) = colnames(Xfull)
	
	if (isTRUE(estimate_only)){
		return(list(
			beta = as.numeric(beta_hat["w"]),
			ssq = NA_real_,
			theta = exp(best$par[num_beta + 2L]),
			log_sigma = best$par[num_beta + 1L],
			best_par = best$par
		))
	}

	hess = best$hessian
	vcov_full = tryCatch(solve(hess), error = function(e) NULL)
	if (is.null(vcov_full) || any(!is.finite(diag(vcov_full)))){
		return(list(
			beta = as.numeric(beta_hat["w"]),
			ssq = NA_real_,
			theta = exp(best$par[num_beta + 2L]),
			log_sigma = best$par[num_beta + 1L],
			best_par = best$par
		))
	}

	rownames(vcov_full) = colnames(vcov_full) = c(colnames(Xfull), "log_sigma", "log_theta")
	ssq = as.numeric(vcov_full["w", "w"])
	if (!is.finite(beta_hat["w"]) || !is.finite(ssq) || ssq <= 0) {
		# If w is not finite or ssq is not valid, we still return the best_par for potential reuse
		return(list(
			beta = as.numeric(beta_hat["w"]),
			ssq = NA_real_,
			theta = exp(best$par[num_beta + 2L]),
			log_sigma = best$par[num_beta + 1L],
			best_par = best$par
		))
	}

	list(
		beta = as.numeric(beta_hat["w"]),
		ssq = ssq,
		theta = exp(best$par[num_beta + 2L]),
		log_sigma = best$par[num_beta + 1L],
		best_par = best$par
	)
}
