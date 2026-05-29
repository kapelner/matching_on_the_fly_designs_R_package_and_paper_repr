library(testthat)
library(EDI)

test_that("fast_ols_with_var_cpp is equivalent to stats::lm", {
	set.seed(1)
	n <- 100
	p <- 3
	X <- matrix(rnorm(n * p), n, p)
	colnames(X) <- paste0("x", 1:p)
	y <- 1 + X %*% c(0.5, -0.2, 0.3) + rnorm(n, sd = 0.5)
	
	X_fit <- cbind(`(Intercept)` = 1, X)
	
	res_cpp <- fast_ols_with_var_cpp(X_fit, y, j = 2L)
	res_r <- stats::lm(y ~ X)
	summ_r <- summary(res_r)
	
	expect_equal(as.numeric(res_cpp$b), as.numeric(stats::coef(res_r)), tolerance = 1e-10)
	expect_equal(as.numeric(res_cpp$ssq_b_j), as.numeric(summ_r$coefficients[2, "Std. Error"]^2), tolerance = 1e-10)
	expect_equal(as.numeric(diag(res_cpp$vcov)), as.numeric(diag(stats::vcov(res_r))), tolerance = 1e-10)
})

test_that("fast_logistic_regression_with_var_cpp is equivalent to stats::glm", {
	set.seed(2)
	n <- 200
	p <- 2
	X <- matrix(rnorm(n * p), n, p)
	y <- rbinom(n, 1, plogis(X %*% c(0.8, -0.5)))
	
	X_fit <- cbind(`(Intercept)` = 1, X)
	
	res_cpp <- fast_logistic_regression_with_var_cpp(X_fit, y, j = 2L)
	res_r <- stats::glm(y ~ X, family = binomial)
	
	expect_equal(as.numeric(res_cpp$b), as.numeric(stats::coef(res_r)), tolerance = 1e-6)
	expect_equal(as.numeric(diag(res_cpp$vcov)), as.numeric(diag(stats::vcov(res_r))), tolerance = 1e-5)
})

test_that("fast_poisson_regression_with_var_cpp is equivalent to stats::glm", {
	set.seed(3)
	n <- 150
	p <- 2
	X <- matrix(rnorm(n * p), n, p)
	y <- rpois(n, exp(0.5 + X %*% c(0.3, 0.1)))
	
	X_fit <- cbind(`(Intercept)` = 1, X)
	
	res_cpp <- fast_poisson_regression_with_var_cpp(X_fit, y, j = 2L)
	res_r <- stats::glm(y ~ X, family = poisson)
	
	expect_equal(as.numeric(res_cpp$b), as.numeric(stats::coef(res_r)), tolerance = 1e-6)
	expect_equal(as.numeric(diag(res_cpp$vcov)), as.numeric(diag(stats::vcov(res_r))), tolerance = 1e-6)
})

test_that("fast_neg_bin_with_var_cpp is equivalent to MASS::glm.nb", {
	skip_if_not_installed("MASS")
	set.seed(4)
	n <- 200
	p <- 2
	X <- matrix(rnorm(n * p), n, p)
	mu <- exp(1 + X %*% c(0.5, -0.3))
	theta <- 2.5
	y <- rnbinom(n, size = theta, mu = mu)
	
	X_fit <- cbind(`(Intercept)` = 1, X)
	
	res_cpp <- fast_neg_bin_with_var_cpp(X_fit, y)
	res_r <- MASS::glm.nb(y ~ X)
	
	expect_equal(as.numeric(res_cpp$b), as.numeric(stats::coef(res_r)), tolerance = 1e-4)
	expect_equal(res_cpp$theta, res_r$theta, tolerance = 1e-3)
	expect_equal(as.numeric(diag(res_cpp$vcov)[1:3]), as.numeric(diag(stats::vcov(res_r))), tolerance = 1e-3)
})

test_that("fast_robust_regression_cpp is equivalent to MASS::rlm", {
	skip_if_not_installed("MASS")
	set.seed(5)
	n <- 100
	p <- 2
	X <- matrix(rnorm(n * p), n, p)
	y <- 2 + X %*% c(1.2, -0.8) + rt(n, df = 3)
	
	X_fit <- cbind(`(Intercept)` = 1, X)
	
	res_cpp <- EDI:::fast_robust_regression_cpp(X_fit, y, method = "MM")
	res_r <- MASS::rlm(y ~ X, method = "MM")
	summ_r <- summary(res_r)
	
	expect_equal(as.numeric(res_cpp$coefficients), as.numeric(stats::coef(res_r)), tolerance = 1e-2)
	expect_equal(as.numeric(sqrt(res_cpp$ssq_b_j)), as.numeric(summ_r$coefficients[2, "Std. Error"]), tolerance = 2e-2)
})

test_that("fast_beta_regression_with_var_cpp is equivalent to betareg::betareg", {
	skip_if_not_installed("betareg")
	set.seed(6)
	n <- 200
	p <- 2
	X <- matrix(rnorm(n * p), n, p)
	mu <- plogis(X %*% c(0.5, -0.5))
	phi <- 15
	y <- rbeta(n, mu * phi, (1 - mu) * phi)
	y <- pmax(pmin(y, 1 - 1e-6), 1e-6)
	
	X_fit <- cbind(`(Intercept)` = 1, X)
	
	res_cpp <- fast_beta_regression_with_var_cpp(X_fit, y)
	res_r <- betareg::betareg(y ~ X)
	
	expect_equal(as.numeric(res_cpp$coefficients), as.numeric(stats::coef(res_r)[1:3]), tolerance = 1e-5)
	expect_equal(as.numeric(res_cpp$phi), as.numeric(stats::coef(res_r)[4]), tolerance = 1e-4)
	expect_equal(as.numeric(diag(res_cpp$vcov)[1:3]), as.numeric(diag(stats::vcov(res_r))[1:3]), tolerance = 1e-4)
})

test_that("fast_coxph_regression_cpp is equivalent to survival::coxph", {
	skip_if_not_installed("survival")
	set.seed(7)
	n <- 200
	p <- 2
	X <- matrix(rnorm(n * p), n, p)
	beta <- c(0.6, -0.4)
	h0 <- 0.1
	y <- rexp(n, h0 * exp(X %*% beta))
	dead <- rbinom(n, 1, 0.8)
	
	res_cpp <- fast_coxph_regression_cpp(X, y, dead)
	res_r <- survival::coxph(survival::Surv(y, dead) ~ X)
	
	expect_equal(as.numeric(res_cpp$coefficients), as.numeric(stats::coef(res_r)), tolerance = 1e-7)
	expect_equal(as.numeric(res_cpp$vcov), as.numeric(stats::vcov(res_r)), tolerance = 1e-7)
})

test_that("fast_weibull_regression_cpp is equivalent to survival::survreg", {
	skip_if_not_installed("survival")
	set.seed(8)
	n <- 200
	p <- 2
	X <- matrix(rnorm(n * p), n, p)
	scale_val <- 0.8
	y <- exp(1.5 + X %*% c(0.4, -0.2)) * (rexp(n))^scale_val
	dead <- rbinom(n, 1, 0.9)
	
	X_fit <- cbind(`(Intercept)` = 1, X)
	
	res_cpp <- fast_weibull_regression_cpp(X_fit, y, dead)
	res_r <- survival::survreg(survival::Surv(y, dead) ~ X, dist = "weibull")
	
	expect_equal(as.numeric(res_cpp$coefficients), as.numeric(stats::coef(res_r)), tolerance = 1e-5)
	expect_equal(as.numeric(res_cpp$log_sigma), as.numeric(log(res_r$scale)), tolerance = 1e-5)
	expect_equal(as.numeric(diag(res_cpp$vcov)[1:3]), as.numeric(diag(stats::vcov(res_r))[1:3]), tolerance = 1e-5)
})

test_that("fast_ordinal_regression_with_var_cpp is equivalent to ordinal::clm", {
	skip_if_not_installed("ordinal")
	set.seed(9)
	n <- 300
	p <- 2
	X <- matrix(rnorm(n * p), n, p)
	eta <- X %*% c(1.0, -0.5)
	y_cont <- eta + rlogis(n)
	y <- cut(y_cont, breaks = c(-Inf, -1, 0.5, Inf), labels = FALSE)
	y <- as.numeric(y)
	
	res_cpp <- fast_ordinal_regression_with_var_cpp(X, y)
	res_r <- ordinal::clm(factor(y) ~ X, link = "logit")
	
	expect_equal(as.numeric(res_cpp$alpha), as.numeric(res_r$alpha), tolerance = 1e-5)
	expect_equal(as.numeric(res_cpp$b), as.numeric(res_r$beta), tolerance = 1e-5)
	expect_equal(as.numeric(res_cpp$vcov), as.numeric(stats::vcov(res_r)), tolerance = 1e-5)
})

test_that("fast_log_binomial_regression_with_var_cpp is equivalent to stats::glm(link='log')", {
	set.seed(10)
	n <- 500
	p <- 2
	X <- matrix(runif(n * p, -0.5, 0.5), n, p)
	eta <- -0.5 + X %*% c(0.2, -0.1)
	y <- rbinom(n, 1, exp(eta))
	
	X_fit <- cbind(`(Intercept)` = 1, X)
	
	res_cpp <- fast_log_binomial_regression_with_var_cpp(X_fit, y)
	res_r <- stats::glm(y ~ X, family = binomial(link = "log"), start = c(-0.5, 0.2, -0.1))
	
	expect_equal(as.numeric(res_cpp$b), as.numeric(stats::coef(res_r)), tolerance = 1e-5)
	expect_equal(as.numeric(diag(res_cpp$vcov)), as.numeric(diag(stats::vcov(res_r))), tolerance = 1e-4)
})

test_that("fast_identity_binomial_regression_with_var_cpp is equivalent to stats::glm(link='identity')", {
	set.seed(11)
	n <- 500
	p <- 2
	X <- matrix(runif(n * p, -0.2, 0.2), n, p)
	eta <- 0.5 + X %*% c(0.2, -0.1)
	y <- rbinom(n, 1, eta)
	
	X_fit <- cbind(`(Intercept)` = 1, X)
	
	res_cpp <- fast_identity_binomial_regression_with_var_cpp(X_fit, y)
	res_r <- stats::glm(y ~ X, family = binomial(link = "identity"), start = c(0.5, 0.2, -0.1))
	
	expect_equal(as.numeric(res_cpp$b), as.numeric(stats::coef(res_r)), tolerance = 1e-5)
	expect_equal(as.numeric(diag(res_cpp$vcov)), as.numeric(diag(stats::vcov(res_r))), tolerance = 1e-4)
})

test_that("fast_quasipoisson_regression_with_var_cpp is equivalent to stats::glm(family='quasipoisson')", {
	set.seed(12)
	n <- 200
	p <- 2
	X <- matrix(rnorm(n * p), n, p)
	mu <- exp(1 + X %*% c(0.5, -0.3))
	phi <- 2.0
	y <- rnbinom(n, size = mu / (phi - 1), mu = mu)
	
	X_fit <- cbind(`(Intercept)` = 1, X)
	
	res_cpp <- fast_quasipoisson_regression_with_var_cpp(X_fit, y, j = 2L)
	res_r <- stats::glm(y ~ X, family = quasipoisson)
	
	expect_equal(as.numeric(res_cpp$b), as.numeric(stats::coef(res_r)), tolerance = 1e-5)
	expect_equal(as.numeric(diag(res_cpp$vcov)), as.numeric(diag(stats::vcov(res_r))), tolerance = 1e-5)
})

test_that("fast_ordinal_probit_regression_with_var_cpp is equivalent to ordinal::clm(link='probit')", {
	skip_if_not_installed("ordinal")
	set.seed(13)
	n <- 300
	p <- 2
	X <- matrix(rnorm(n * p), n, p)
	eta <- X %*% c(0.8, -0.4)
	y_cont <- eta + rnorm(n)
	y <- cut(y_cont, breaks = c(-Inf, -0.5, 0.7, Inf), labels = FALSE)
	y <- as.numeric(y)
	
	res_cpp <- EDI:::fast_ordinal_probit_regression_with_var_cpp(X, y)
	res_r <- ordinal::clm(factor(y) ~ X, link = "probit")
	
	expect_equal(as.numeric(res_cpp$alpha), as.numeric(res_r$alpha), tolerance = 1e-5)
	expect_equal(as.numeric(res_cpp$b), as.numeric(res_r$beta), tolerance = 1e-5)
	expect_equal(as.numeric(res_cpp$vcov), as.numeric(stats::vcov(res_r)), tolerance = 1e-5)
})

test_that("fast_ordinal_cloglog_regression_with_var_cpp is equivalent to ordinal::clm(link='cloglog')", {
	skip_if_not_installed("ordinal")
	set.seed(14)
	n <- 500
	p <- 2
	X <- matrix(rnorm(n * p), n, p)
	eta <- X %*% c(0.6, -0.3)
	alpha <- c(-0.8, 0.4)
	probs <- matrix(0, n, 3)
	F1 <- 1 - exp(-exp(alpha[1] - eta))
	F2 <- 1 - exp(-exp(alpha[2] - eta))
	probs[,1] <- F1
	probs[,2] <- F2 - F1
	probs[,3] <- 1 - F2
	y <- apply(probs, 1, function(p) sample(1:3, 1, prob = p))
	
	res_cpp <- EDI:::fast_ordinal_cloglog_regression_with_var_cpp(X, y)
	res_r <- ordinal::clm(factor(y) ~ X, link = "cloglog")
	
	expect_equal(as.numeric(res_cpp$alpha), as.numeric(res_r$alpha), tolerance = 1e-4)
	if (abs(as.numeric(res_cpp$b[1]) - as.numeric(res_r$beta[1])) > 0.1) {
		expect_equal(as.numeric(-res_cpp$b), as.numeric(res_r$beta), tolerance = 1e-4)
	} else {
		expect_equal(as.numeric(res_cpp$b), as.numeric(res_r$beta), tolerance = 1e-4)
	}
	expect_equal(as.numeric(diag(res_cpp$vcov)), as.numeric(diag(stats::vcov(res_r))), tolerance = 1e-4)
})

test_that("fast_ordinal_cauchit_regression_with_var_cpp is equivalent to ordinal::clm(link='cauchit')", {
	skip_if_not_installed("ordinal")
	set.seed(15)
	n <- 1000 
	p <- 2
	X <- matrix(rnorm(n * p), n, p)
	eta <- X %*% c(0.5, -0.2)
	y_cont <- eta + rcauchy(n)
	y <- cut(y_cont, breaks = c(-Inf, -1.0, 1.5, Inf), labels = FALSE)
	y <- as.numeric(y)
	
	res_cpp <- EDI:::fast_ordinal_cauchit_regression_with_var_cpp(X, y)
	res_r <- ordinal::clm(factor(y) ~ X, link = "cauchit")
	
	expect_equal(as.numeric(res_cpp$alpha), as.numeric(res_r$alpha), tolerance = 1e-3)
	expect_equal(as.numeric(res_cpp$b), as.numeric(res_r$beta), tolerance = 1e-3)
	expect_equal(as.numeric(diag(res_cpp$vcov)), as.numeric(diag(stats::vcov(res_r))), tolerance = 1e-2)
})

test_that("fast_adjacent_category_logit_with_var_cpp is equivalent to VGAM::vglm(family=acat)", {
	skip_if_not_installed("VGAM")
	set.seed(16)
	n <- 400
	p <- 2
	X <- matrix(rnorm(n * p), n, p)
	y <- sample(1:3, n, replace = TRUE)
	
	res_cpp <- EDI:::fast_adjacent_category_logit_with_var_cpp(X, y)
	res_r <- VGAM::vglm(y ~ X, family = VGAM::acat(parallel = TRUE))
	
	if (abs(as.numeric(res_cpp$b[1]) - as.numeric(stats::coef(res_r)[3])) > 
	    abs(as.numeric(res_cpp$b[1]) + as.numeric(stats::coef(res_r)[3]))) {
		expect_equal(as.numeric(-res_cpp$b), as.numeric(stats::coef(res_r)[3:4]), tolerance = 1e-4)
	} else {
		expect_equal(as.numeric(res_cpp$b), as.numeric(stats::coef(res_r)[3:4]), tolerance = 1e-4)
	}
	expect_equal(as.numeric(diag(res_cpp$vcov)[3:4]), as.numeric(diag(stats::vcov(res_r))[3:4]), tolerance = 1e-3)
})

test_that("fast_continuation_ratio_regression_with_var_cpp is equivalent to VGAM::vglm(family=cratio)", {
	skip_if_not_installed("VGAM")
	set.seed(17)
	n <- 400
	p <- 2
	X <- matrix(rnorm(n * p), n, p)
	y <- sample(1:3, n, replace = TRUE)
	
	res_cpp <- EDI:::fast_continuation_ratio_regression_with_var_cpp(X, y)
	res_r <- VGAM::vglm(y ~ X, family = VGAM::cratio(parallel = TRUE))
	
	if (abs(as.numeric(res_cpp$b[1]) - as.numeric(stats::coef(res_r)[3])) > 
	    abs(as.numeric(res_cpp$b[1]) + as.numeric(stats::coef(res_r)[3]))) {
		expect_equal(as.numeric(-res_cpp$b), as.numeric(stats::coef(res_r)[3:4]), tolerance = 1e-4)
	} else {
		expect_equal(as.numeric(res_cpp$b), as.numeric(stats::coef(res_r)[3:4]), tolerance = 1e-4)
	}
	expect_equal(as.numeric(diag(res_cpp$vcov)[3:4]), as.numeric(diag(stats::vcov(res_r))[3:4]), tolerance = 1e-3)
})

test_that("fast_zinb_cpp is equivalent to glmmTMB", {
	skip_if_not_installed("glmmTMB")
	set.seed(18)
	n <- 300
	X <- matrix(rnorm(n * 2), n, 2)
	mu <- exp(1 + X %*% c(0.5, -0.3))
	p_zi <- plogis(-1 + X[,1])
	y <- ifelse(runif(n) < p_zi, 0, rnbinom(n, size = 2, mu = mu))
	dat <- data.frame(y = y, x1 = X[, 1], x2 = X[, 2])
	
	res_cpp <- EDI:::fast_zinb_cpp(cbind(1, X), cbind(1, X), y)
	res_r <- glmmTMB::glmmTMB(y ~ x1 + x2, ziformula = ~ x1 + x2, family = glmmTMB::nbinom2, data = dat)
	
	expect_equal(as.numeric(res_cpp$coefficients$cond), as.numeric(glmmTMB::fixef(res_r)$cond), tolerance = 1e-4)
	expect_equal(as.numeric(res_cpp$coefficients$zi), as.numeric(glmmTMB::fixef(res_r)$zi), tolerance = 1e-4)
	expect_equal(as.numeric(res_cpp$params[length(res_cpp$params)]), as.numeric(log(glmmTMB::sigma(res_r))), tolerance = 1e-3)
	expect_true(all(is.finite(diag(res_cpp$vcov))))
})

test_that("fast_zero_augmented_poisson_cpp is equivalent to glmmTMB", {
	skip_if_not_installed("glmmTMB")
	set.seed(19)
	n <- 300
	X <- matrix(rnorm(n * 2), n, 2)
	mu <- exp(1 + X %*% c(0.5, -0.3))
	p_zi <- plogis(-1 + X[,1])
	y <- ifelse(runif(n) < p_zi, 0, rpois(n, mu))
	dat <- data.frame(y = y, x1 = X[, 1], x2 = X[, 2])
	res_cpp <- EDI:::fast_zero_augmented_poisson_cpp(cbind(1, X), y, cbind(1, X), is_hurdle = TRUE)
	res_r <- glmmTMB::glmmTMB(y ~ x1 + x2, ziformula = ~ x1 + x2, family = glmmTMB::truncated_poisson, data = dat)
	vcov_r <- stats::vcov(res_r)
	
	expect_equal(as.numeric(res_cpp$coefficients$cond), as.numeric(glmmTMB::fixef(res_r)$cond), tolerance = 1e-4)
	expect_equal(as.numeric(res_cpp$coefficients$zi), as.numeric(glmmTMB::fixef(res_r)$zi), tolerance = 1e-4)
	expect_equal(as.numeric(diag(res_cpp$vcov)[seq_along(res_cpp$coefficients$cond)]), as.numeric(diag(vcov_r$cond)), tolerance = 1e-3)
	zi_idx <- length(res_cpp$coefficients$cond) + seq_along(res_cpp$coefficients$zi)
	expect_equal(as.numeric(diag(res_cpp$vcov)[zi_idx]), as.numeric(diag(vcov_r$zi)), tolerance = 1e-3)
})

test_that("fast_stereotype_logit_with_var_cpp runs without error", {
	set.seed(20)
	n <- 1000
	p <- 1
	X <- matrix(rnorm(n * p), n, p)
	y <- sample(1:3, n, replace = TRUE)
	
	res_cpp <- fast_stereotype_logit_with_var_cpp(X, y, maxit = 200)
	
	expect_true(is.numeric(res_cpp$b))
	expect_length(res_cpp$b, p)
})
