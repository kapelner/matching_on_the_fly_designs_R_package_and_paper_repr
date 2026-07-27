library(EDI)

# test-brt-smoothed-weibull-kernel.R only checks the fast kernel against itself (noise
# ordering, NULL vs zero noise), so it would still pass if the kernel fit the wrong model
# or read off the wrong coefficient. These tests pin the kernel to the reference paths:
# survreg, and the generic per-draw R fallback it short-circuits.

make_weibull_kk_design = function(seed, n = 60){
	set.seed(seed)
	des = DesignSeqOneByOneKK14$new(n = n, response_type = "survival", verbose = FALSE)
	for (i in 1:n) des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
	w01 = as.numeric(des$get_w() == 1)
	y_lat = rexp(n) * exp(0.7 * w01 + rnorm(n, sd = 0.3))
	cens = as.numeric(quantile(y_lat, 0.85))
	add_all_subject_responses_seq(des, pmin(y_lat, cens), deads = as.numeric(y_lat <= cens))
	des
}

test_that("compute_weibull_rand_bootstrap_parallel_cpp reproduces survreg on the same resamples", {
	skip_if_not_installed("survival")
	n = 60
	des = make_weibull_kk_design(20260726, n)
	inf = InferenceSurvivalKKWeibullMarginal$new(des, verbose = FALSE)
	inf$compute_estimate()
	ipriv = inf$.__enclos_env__$private

	y0 = as.numeric(ipriv$y)
	dead0 = as.integer(ipriv$dead)
	X_data = ipriv$get_X()
	best = ipriv$best_X_colnames
	Xc = if (!is.null(best) && length(intersect(best, colnames(X_data))) > 0) {
		as.matrix(X_data[, intersect(best, colnames(X_data)), drop = FALSE])
	} else as.matrix(X_data)

	set.seed(4242)
	B = 25
	i_mat = matrix(NA_integer_, n, B)
	w_mat = matrix(NA_integer_, n, B)
	for (b in seq_len(B)) {
		i_mat[, b] = sample.int(n, n, replace = TRUE)
		wb = rbinom(n, 1, 0.5); wb[1] = 1L; wb[2] = 0L
		w_mat[, b] = as.integer(wb)
	}
	delta = 0.3

	fast = EDI:::compute_weibull_rand_bootstrap_parallel_cpp(y0, dead0, Xc, i_mat, w_mat, delta, NULL, 1L)

	ref = rep(NA_real_, B)
	for (b in seq_len(B)) {
		idx = i_mat[, b]; wb = w_mat[, b]
		yb = y0[idx]
		yb[wb == 1] = yb[wb == 1] * exp(delta)   # AFT sharp null: multiplicative on the time scale
		Xb = Xc[idx, , drop = FALSE]
		dat = data.frame(.y__ = yb, .dead__ = dead0[idx], treatment = wb, Xb, check.names = FALSE)
		fit = tryCatch(suppressWarnings(survival::survreg(
			stats::as.formula(paste0("survival::Surv(.y__, .dead__) ~ ",
				paste(c("treatment", colnames(Xb)), collapse = " + "))),
			data = dat, dist = "weibull")), error = function(e) NULL)
		if (!is.null(fit)) ref[b] = as.numeric(stats::coef(fit)[["treatment"]])
	}

	both = is.finite(fast) & is.finite(ref)
	expect_gt(sum(both), B / 2)
	# lbfgs (tol 1e-8) vs survreg's Newton solver: agreement is optimizer-tolerance limited
	expect_equal(fast[both], ref[both], tolerance = 1e-3)
})

test_that("the fast Weibull BRT kernel matches the generic per-draw fallback it replaces", {
	n = 60
	des = make_weibull_kk_design(20260726, n)
	inf = InferenceSurvivalKKWeibullMarginal$new(des, verbose = FALSE)
	inf$compute_estimate()
	ipriv = inf$.__enclos_env__$private
	y0 = as.numeric(ipriv$y)

	set.seed(7)
	draws = ipriv$generate_rand_bootstrap_draws(30, materialize_w = TRUE)

	for (delta in c(0, 0.3)) {
		fast = ipriv$compute_fast_rand_bootstrap_distr(y0, draws, delta, "log")
		expect_false(is.null(fast))
		slow = vapply(seq_along(draws), function(b) {
			tryCatch(suppressWarnings(ipriv$run_rand_bootstrap_iteration(draws[[b]], delta, "log", y0)),
				error = function(e) NA_real_)
		}, numeric(1))

		# the two paths must agree on WHICH draws are estimable, not merely on the values
		expect_equal(is.finite(fast), is.finite(slow), info = paste("delta =", delta))
		both = is.finite(fast) & is.finite(slow)
		expect_gt(sum(both), 0)
		expect_equal(fast[both], slow[both], tolerance = 1e-3, info = paste("delta =", delta))
	}
})
