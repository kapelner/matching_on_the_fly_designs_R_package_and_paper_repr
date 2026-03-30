test_that("KK Wilcox rank-regression fast bootstrap matches the generic KK bootstrap", {
	skip_if_not_installed("Rfit")

	SlowInferenceAllKKWilcoxRegrMultiIVWC = R6::R6Class(
		"SlowInferenceAllKKWilcoxRegrMultiIVWC",
		inherit = InferenceAllKKWilcoxRegrMultiIVWC,
		private = list(
			compute_fast_bootstrap_distr = function(B, i_reservoir, n_reservoir, m, y, w, m_vec){
				NULL
			}
		)
	)

	set.seed(20260329)
	n = 30
	p = 5
	X = as.data.frame(matrix(rnorm(n * p), nrow = n, ncol = p))
	colnames(X) = paste0("x", seq_len(p))
	y = as.numeric(rnorm(n))

	des = DesignSeqOneByOneKK14$new(n = n, response_type = "continuous", verbose = FALSE)
	for (i in seq_len(n)) {
		w_i = des$add_subject_to_experiment_and_assign(X[i, , drop = FALSE])
		des$add_subject_response(i, y[i] + 0.2 * w_i)
	}

	fast_inf = InferenceAllKKWilcoxRegrMultiIVWC$new(des, num_cores = 1, verbose = FALSE)
	slow_inf = SlowInferenceAllKKWilcoxRegrMultiIVWC$new(des, num_cores = 1, verbose = FALSE)

	set.seed(44)
	fast_boot = suppressWarnings(
		fast_inf$approximate_bootstrap_distribution_beta_hat_T(B = 9, show_progress = FALSE)
	)
	set.seed(44)
	slow_boot = suppressWarnings(
		slow_inf$approximate_bootstrap_distribution_beta_hat_T(B = 9, show_progress = FALSE)
	)

	expect_equal(fast_boot, slow_boot, tolerance = 1e-10)
})

test_that("KK Wilcox rank-regression low-level exact solver matches Rfit formula path", {
	skip_if_not_installed("Rfit")

	old_fit_matched = function(dy, dX){
		dat = as.data.frame(as.matrix(dX))
		colnames(dat) = paste0("x", seq_len(ncol(dat)))
		dat$dy = dy
		mod = suppressWarnings(Rfit::rfit(dy ~ ., data = dat))
		summ = suppressWarnings(summary(mod))
		list(
			beta = as.numeric(summ$coefficients["(Intercept)", "Estimate"]),
			ssq = as.numeric(summ$coefficients["(Intercept)", "Std. Error"])^2
		)
	}

	old_fit_reservoir = function(y_r, w_r, X_r){
		dat = data.frame(y = y_r, w = w_r)
		X_covs = as.data.frame(as.matrix(X_r))
		colnames(X_covs) = paste0("x", seq_len(ncol(X_covs)))
		dat = cbind(dat, X_covs)
		mod = suppressWarnings(Rfit::rfit(y ~ ., data = dat))
		summ = suppressWarnings(summary(mod))
		list(
			beta = as.numeric(summ$coefficients["w", "Estimate"]),
			ssq = as.numeric(summ$coefficients["w", "Std. Error"])^2
		)
	}

	set.seed(20260330)
	n = 12
	X = as.data.frame(matrix(rnorm(n * 2), nrow = n, ncol = 2))
	colnames(X) = c("x1", "x2")
	y = as.numeric(rnorm(n))
	des = DesignSeqOneByOneKK14$new(n = n, response_type = "continuous", verbose = FALSE)
	for (i in seq_len(n)) {
		w_i = des$add_subject_to_experiment_and_assign(X[i, , drop = FALSE])
		des$add_subject_response(i, y[i] + 0.1 * w_i)
	}
	inf = InferenceAllKKWilcoxRegrMultiIVWC$new(des, num_cores = 1, verbose = FALSE)
	priv = inf$.__enclos_env__$private

	set.seed(20260331)
	dy = rnorm(40)
	dX = matrix(rnorm(40 * 6), nrow = 40, ncol = 6)
	colnames(dX) = paste0("x", seq_len(ncol(dX)))
	fast_matched = suppressWarnings(priv$fit_matched_component(dy, dX))
	old_matched = suppressWarnings(old_fit_matched(dy, dX))
	expect_equal(fast_matched, old_matched, tolerance = 1e-8)

	set.seed(20260401)
	y_r = rnorm(80)
	w_r = rep(0:1, each = 40)
	X_r = matrix(rnorm(80 * 6), nrow = 80, ncol = 6)
	colnames(X_r) = paste0("x", seq_len(ncol(X_r)))
	fast_reservoir = suppressWarnings(priv$fit_reservoir_component(y_r, w_r, X_r))
	old_reservoir = suppressWarnings(old_fit_reservoir(y_r, w_r, X_r))
	expect_equal(fast_reservoir, old_reservoir, tolerance = 1e-8)
})
