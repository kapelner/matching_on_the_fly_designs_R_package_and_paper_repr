library(testthat)
library(EDI)

test_that("fast_zero_one_inflated_beta_cpp factors into betareg and multinom submodels", {
	skip_if_not_installed("betareg")
	skip_if_not_installed("nnet")
	set.seed(2001)
	n = 500L
	x1 = rnorm(n)
	x2 = rnorm(n)
	X = cbind(1, x1, x2)
	X_zero_one = cbind(1, x1)
	beta = c(-0.2, 0.5, -0.35)
	phi = 10
	gamma0 = c(-1.0, 0.6)
	gamma1 = c(-1.3, -0.4)
	mu = plogis(drop(X %*% beta))
	eta0 = drop(X_zero_one %*% gamma0)
	eta1 = drop(X_zero_one %*% gamma1)
	denom = 1 + exp(eta0) + exp(eta1)
	p0 = exp(eta0) / denom
	p1 = exp(eta1) / denom
	u = runif(n)
	y = ifelse(u < p0, 0, ifelse(u < p0 + p1, 1, rbeta(n, mu * phi, (1 - mu) * phi)))

	fit_cpp = EDI:::fast_zero_one_inflated_beta_cpp(X, X_zero_one, y)

	is_beta = y > 0 & y < 1
	fit_beta = betareg::betareg(y[is_beta] ~ x1[is_beta] + x2[is_beta])
	expect_equal(as.numeric(fit_cpp$b), as.numeric(stats::coef(fit_beta)[1:3]), tolerance = 1e-3)
	expect_equal(as.numeric(fit_cpp$log_phi), log(as.numeric(stats::coef(fit_beta)[4])), tolerance = 1e-3)

	cls = factor(ifelse(y <= 0, "zero", ifelse(y >= 1, "one", "beta")),
		levels = c("beta", "zero", "one"))
	fit_multinom = nnet::multinom(cls ~ x1, trace = FALSE, maxit = 500)
	multinom_coef = stats::coef(fit_multinom)[c("zero", "one"), , drop = FALSE]
	expect_equal(as.numeric(fit_cpp$zero_one_b0), as.numeric(multinom_coef["zero", ]), tolerance = 5e-3)
	expect_equal(as.numeric(fit_cpp$zero_one_b1), as.numeric(multinom_coef["one", ]), tolerance = 5e-3)
})

test_that("zero-one-inflated beta likelihood matches gamlss.dist::dBEINF", {
	skip_if_not_installed("gamlss.dist")
	set.seed(2006)
	n = 32L
	X = cbind(1, rnorm(n))
	X_zero_one = cbind(1, rnorm(n))
	y = c(0, 1, 0, 1, seq(0.04, 0.96, length.out = n - 4L))
	params = c(0.15, -0.2, log(8), -1.1, 0.25, -1.3, -0.15)
	p = ncol(X)
	p_zero_one = ncol(X_zero_one)
	mu = plogis(drop(X %*% params[seq_len(p)]))
	phi = exp(params[p + 1L])
	eta0 = drop(X_zero_one %*% params[p + 1L + seq_len(p_zero_one)])
	eta1 = drop(X_zero_one %*% params[p + 1L + p_zero_one + seq_len(p_zero_one)])
	loglik_ref = sum(gamlss.dist::dBEINF(
		y,
		mu = mu,
		sigma = sqrt(1 / (phi + 1)),
		nu = exp(eta0),
		tau = exp(eta1),
		log = TRUE
	))

	nll = function(par) {
		mu_par = plogis(drop(X %*% par[seq_len(p)]))
		phi_par = exp(par[p + 1L])
		eta0_par = drop(X_zero_one %*% par[p + 1L + seq_len(p_zero_one)])
		eta1_par = drop(X_zero_one %*% par[p + 1L + p_zero_one + seq_len(p_zero_one)])
		-sum(gamlss.dist::dBEINF(
			y,
			mu = mu_par,
			sigma = sqrt(1 / (phi_par + 1)),
			nu = exp(eta0_par),
			tau = exp(eta1_par),
			log = TRUE
		))
	}
	score_cpp = as.numeric(EDI:::get_zero_one_inflated_beta_score_cpp(X, X_zero_one, y, params))
	score_ref = -as.numeric(numDeriv::grad(nll, params))
	expect_equal(score_cpp, score_ref, tolerance = 1e-5)

	loglik_cpp = -nll(params)
	expect_equal(loglik_cpp, loglik_ref, tolerance = 1e-12)
})

test_that("fast_truncated_negbin_count_cpp matches glmmTMB truncated_nbinom2", {
	skip_if_not_installed("glmmTMB")
	set.seed(2002)
	n = 300L
	x1 = rnorm(n)
	x2 = rnorm(n)
	X = cbind(1, x1, x2)
	mu = exp(0.2 + 0.4 * x1 - 0.25 * x2)
	theta = 2.3
	y = numeric(n)
	for (i in seq_len(n)) {
		repeat {
			draw = rnbinom(1L, mu = mu[i], size = theta)
			if (draw > 0L) {
				y[i] = draw
				break
			}
		}
	}

	fit_cpp = EDI:::fast_truncated_negbin_count_cpp(X, y)
	fit_ref = glmmTMB::glmmTMB(
		y ~ x1 + x2,
		family = glmmTMB::truncated_nbinom2,
		data = data.frame(y = y, x1 = x1, x2 = x2)
	)

	expect_true(fit_cpp$converged)
	expect_equal(as.numeric(fit_cpp$b), as.numeric(glmmTMB::fixef(fit_ref)$cond), tolerance = 1e-5)
	expect_equal(as.numeric(tail(fit_cpp$params, 1L)), log(as.numeric(glmmTMB::sigma(fit_ref))), tolerance = 1e-5)
})

test_that("fast_cpoisson_combined_with_var_cpp reduces to canonical GLMs in single-component cases", {
	set.seed(2003)
	nd = 200L
	p = 2L
	X_diff = matrix(rnorm(nd * p), nd, p)
	n_k = sample(2:8, nd, replace = TRUE)
	yT = rbinom(nd, n_k, plogis(0.2 + X_diff %*% c(0.4, -0.3)))
	fit_pair_cpp = EDI:::fast_cpoisson_combined_with_var_cpp(
		yT, n_k, X_diff,
		numeric(0), numeric(0), matrix(numeric(0), 0L, p),
		fixed_idx = 1L, fixed_values = 0
	)
	fit_pair_ref = stats::glm(cbind(yT, n_k - yT) ~ cbind(1, X_diff) - 1, family = binomial())
	expect_true(fit_pair_cpp$converged)
	expect_equal(as.numeric(fit_pair_cpp$params[-1L]), as.numeric(stats::coef(fit_pair_ref)), tolerance = 1e-8)

	nR = 250L
	X_r = matrix(rnorm(nR * p), nR, p)
	w_r = rbinom(nR, 1L, 0.5)
	y_r = rpois(nR, exp(0.1 + 0.5 * w_r + X_r %*% c(0.2, -0.1)))
	fit_res_cpp = EDI:::fast_cpoisson_combined_with_var_cpp(
		numeric(0), numeric(0), matrix(numeric(0), 0L, p),
		y_r, w_r, X_r
	)
	fit_res_ref = stats::glm(y_r ~ w_r + X_r, family = poisson())
	expect_true(fit_res_cpp$converged)
	expect_equal(as.numeric(fit_res_cpp$params), as.numeric(stats::coef(fit_res_ref)), tolerance = 1e-8)
})

test_that("fast_dep_cens_transform_optim_cpp rho=0 score matches two lognormal survreg fits", {
	skip_if_not_installed("survival")
	set.seed(2004)
	n = 400L
	x1 = rnorm(n)
	x2 = rnorm(n)
	X = cbind(1, x1, x2)
	p = ncol(X)
	t_event = exp(drop(X %*% c(0.4, 0.3, -0.2)) + 0.7 * rnorm(n))
	t_cens = exp(drop(X %*% c(0.2, -0.15, 0.25)) + 0.9 * rnorm(n))
	y = pmin(t_event, t_cens)
	dead = as.integer(t_event <= t_cens)

	fit_event = survival::survreg(survival::Surv(y, dead) ~ x1 + x2, dist = "lognormal")
	fit_cens = survival::survreg(survival::Surv(y, 1L - dead) ~ x1 + x2, dist = "lognormal")
	params = c(
		as.numeric(stats::coef(fit_event)),
		as.numeric(stats::coef(fit_cens)),
		log(fit_event$scale),
		log(fit_cens$scale),
		0
	)
	score = as.numeric(EDI:::get_dep_cens_transform_score_cpp(X, y, dead, params))
	free_score = score[-length(score)]
	expect_lt(max(abs(free_score)), 1e-5)
})

test_that("fast_clayton_weibull_aft_optim_cpp singleton-only model matches survreg Weibull", {
	skip_if_not_installed("survival")
	set.seed(2005)
	n = 250L
	x1 = rnorm(n)
	x2 = rnorm(n)
	X = cbind(1, x1, x2)
	p = ncol(X)
	y = exp(drop(X %*% c(0.3, 0.4, -0.2)) + 0.8 * rnorm(n))
	dead = rbinom(n, 1L, 0.75)

	fit_ref = survival::survreg(survival::Surv(y, dead) ~ x1 + x2, dist = "weibull")
	start = c(as.numeric(stats::coef(fit_ref)), log(fit_ref$scale), -2)
	score = as.numeric(EDI:::get_clayton_weibull_aft_score_cpp(
		X, y, dead,
		matrix(0L, 0L, 2L),
		as.integer(seq_len(n) - 1L),
		start
	))
	hessian = EDI:::get_clayton_weibull_aft_hessian_cpp(
		X, y, dead,
		matrix(0L, 0L, 2L),
		as.integer(seq_len(n) - 1L),
		start
	)

	expect_lt(max(abs(score)), 1e-6)
	expect_true(all(is.finite(hessian)))
	expect_equal(hessian, t(hessian), tolerance = 1e-8)
})

test_that("Clayton Weibull AFT pair score matches copula package reference likelihood", {
	skip_if_not_installed("copula")
	set.seed(2007)
	n = 8L
	X = cbind(1, rnorm(n))
	y = exp(rnorm(n))
	dead = c(0, 1, 1, 0, 1, 1, 0, 0)
	pair_idx = matrix(c(0L, 1L, 2L, 3L, 4L, 5L), 3L, 2L, byrow = TRUE)
	singleton_rows = as.integer(c(6L, 7L))
	params = c(0.2, -0.1, log(0.8), log(1.4))

	copula_ref_nll = function(par) {
		p = ncol(X)
		beta = par[seq_len(p)]
		sigma = exp(par[p + 1L])
		theta = exp(min(par[p + 2L], 10))
		eta = drop(X %*% beta)
		log_y = log(y)
		H = exp((log_y - eta) / sigma)
		log_f = (log_y - eta) / sigma - log(sigma) - log_y - H
		S = exp(-H)
		cop = copula::claytonCopula(theta, dim = 2L)
		loglik = 0
		for (k in seq_len(nrow(pair_idx))) {
			i1 = pair_idx[k, 1L] + 1L
			i2 = pair_idx[k, 2L] + 1L
			u = c(S[i1], S[i2])
			d = c(dead[i1], dead[i2])
			if (all(d == 0)) {
				loglik = loglik + log(copula::pCopula(matrix(u, nrow = 1L), cop))
			} else if (all(d == 1)) {
				loglik = loglik + log_f[i1] + log_f[i2] +
					copula::dCopula(matrix(u, nrow = 1L), cop, log = TRUE)
			} else if (d[1L] == 1) {
				dC_du1 = copula::cCopula(matrix(u, nrow = 1L), cop)[1L, 2L]
				loglik = loglik + log_f[i1] + log(dC_du1)
			} else {
				dC_du2 = copula::cCopula(matrix(rev(u), nrow = 1L), cop)[1L, 2L]
				loglik = loglik + log_f[i2] + log(dC_du2)
			}
		}
		for (i0 in singleton_rows) {
			i = i0 + 1L
			loglik = loglik + if (dead[i] > 0.5) log_f[i] else log(S[i])
		}
		-loglik
	}

	score_cpp = as.numeric(EDI:::get_clayton_weibull_aft_score_cpp(
		X, y, dead, pair_idx, singleton_rows, params
	))
	score_ref = -as.numeric(numDeriv::grad(copula_ref_nll, params))
	expect_equal(score_cpp, score_ref, tolerance = 1e-5)
})
