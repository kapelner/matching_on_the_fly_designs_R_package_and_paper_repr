test_that("OLS smart-start helpers match the expected transformed regressions", {
	set.seed(42)
	X_cov <- matrix(rnorm(80), ncol = 2)
	X <- cbind(1, X_cov)
	beta <- c(0.5, -0.75, 1.25)
	y_cont <- as.numeric(X %*% beta + rnorm(nrow(X), sd = 0.1))
	y_count <- rpois(nrow(X), lambda = exp(X %*% c(0.2, 0.1, -0.15)))

	expect_equal(
		as.numeric(EDI:::test_ols_warm_start_beta_cpp(X, y_cont)),
		as.numeric(coef(lm(y_cont ~ X - 1))),
		tolerance = 1e-8
	)
	expect_equal(
		as.numeric(EDI:::test_ols_warm_start_beta_on_log1p_cpp(X, y_count)),
		as.numeric(coef(lm(log1p(y_count) ~ X - 1))),
		tolerance = 1e-8
	)
})

test_that("start finalization falls back cleanly and preserves fixed values", {
	smart_bad <- c(NA_real_, 2, 3)
	legacy <- c(0, 0, 0)
	out <- EDI:::test_finalize_warm_start_beta_cpp(
		smart_bad,
		legacy,
		use_smart = TRUE,
		fixed_idx = 2L,
		fixed_values = 9
	)
	expect_equal(as.numeric(out), c(0, 9, 0))
})

test_that("Weibull smart starts use uncensored rows and ordinal thresholds are ordered", {
	set.seed(7)
	X_cov <- matrix(rnorm(60), ncol = 2)
	X <- cbind(1, X_cov)
	y <- exp(0.8 + X_cov[, 1] - 0.5 * X_cov[, 2] + rnorm(nrow(X_cov), sd = 0.05))
	dead <- rep(1, nrow(X_cov))
	dead[seq(1, nrow(X_cov), by = 3)] <- 0
	y[dead == 0] <- y[dead == 0] * 20

	weibull_start <- EDI:::test_weibull_aft_start_cpp(X, y, dead)
	uncens_coef <- as.numeric(coef(lm(log(y[dead == 1]) ~ X[dead == 1, ] - 1)))
	expect_equal(as.numeric(weibull_start$beta), uncens_coef, tolerance = 1e-6)

	y_ord <- as.numeric(cut(X_cov[, 1] + 0.4 * X_cov[, 2] + rnorm(nrow(X_cov)),
		breaks = c(-Inf, -0.5, 0.4, Inf), labels = FALSE))
	ord_start <- EDI:::test_ordinal_start_cpp(X_cov, y_ord, link = "logit")
	expect_equal(length(ord_start$params), ncol(X_cov) + length(unique(y_ord)) - 1L)
	expect_true(all(diff(as.numeric(ord_start$alpha)) > 0))
})
