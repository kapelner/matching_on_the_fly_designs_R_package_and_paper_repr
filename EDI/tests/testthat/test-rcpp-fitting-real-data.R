library(testthat)
library(EDI)

test_that("Continuous: fast_ols_with_var_cpp on MASS::Boston", {
	skip_if_not_installed("MASS")
	data("Boston", package = "MASS")
	
	# Clean and prepare
	y <- Boston$medv
	X_data <- Boston[, setdiff(names(Boston), "medv")]
	X_fit <- cbind(`(Intercept)` = 1, as.matrix(X_data))
	
	res_cpp <- fast_ols_with_var_cpp(X_fit, y)
	res_r <- stats::lm(medv ~ ., data = Boston)
	
	expect_equal(as.numeric(res_cpp$b), as.numeric(stats::coef(res_r)), tolerance = 1e-10)
	expect_equal(as.numeric(diag(res_cpp$vcov)), as.numeric(diag(stats::vcov(res_r))), tolerance = 1e-10)
})

test_that("Incidence: fast_logistic_regression_with_var_cpp on MASS::birthwt", {
	skip_if_not_installed("MASS")
	data("birthwt", package = "MASS")
	
	# binary outcome: low birthweight
	y <- birthwt$low
	# Use subset of predictors
	X_data <- birthwt[, c("age", "lwt", "smoke", "ht", "ui")]
	X_fit <- cbind(`(Intercept)` = 1, as.matrix(X_data))
	
	res_cpp <- fast_logistic_regression_with_var_cpp(X_fit, y)
	res_r <- stats::glm(low ~ age + lwt + smoke + ht + ui, data = birthwt, family = binomial)
	
	expect_equal(as.numeric(res_cpp$b), as.numeric(stats::coef(res_r)), tolerance = 1e-6)
	expect_equal(as.numeric(diag(res_cpp$vcov)), as.numeric(diag(stats::vcov(res_r))), tolerance = 1e-5)
})

test_that("Count: fast_poisson_regression_with_var_cpp on MASS::quine", {
	skip_if_not_installed("MASS")
	data("quine", package = "MASS")
	
	y <- quine$Days
	# Convert factors to dummy variables manually for Rcpp
	X_data <- model.matrix(Days ~ Eth + Sex + Age, data = quine)
	
	res_cpp <- fast_poisson_regression_with_var_cpp(X_data, y)
	res_r <- stats::glm(Days ~ Eth + Sex + Age, data = quine, family = poisson)
	
	expect_equal(as.numeric(res_cpp$b), as.numeric(stats::coef(res_r)), tolerance = 1e-6)
	expect_equal(as.numeric(diag(res_cpp$vcov)), as.numeric(diag(stats::vcov(res_r))), tolerance = 1e-6)
})

test_that("Count: fast_neg_bin_with_var_cpp on MASS::quine", {
	skip_if_not_installed("MASS")
	data("quine", package = "MASS")
	
	y <- quine$Days
	X_data <- model.matrix(Days ~ Eth + Sex + Age, data = quine)
	
	res_cpp <- fast_neg_bin_with_var_cpp(X_data, y)
	res_r <- MASS::glm.nb(Days ~ Eth + Sex + Age, data = quine)
	
	expect_equal(as.numeric(res_cpp$b), as.numeric(stats::coef(res_r)), tolerance = 1e-4)
	expect_equal(res_cpp$theta, res_r$theta, tolerance = 1e-3)
	# Looser tolerance for real data vcov (observed vs expected info)
	expect_equal(as.numeric(diag(res_cpp$vcov)[1:ncol(X_data)]), as.numeric(diag(stats::vcov(res_r))), tolerance = 1e-2)
})

test_that("Proportion: fast_beta_regression_with_var_cpp on betareg::ReadingSkills", {
	skip_if_not_installed("betareg")
	data("ReadingSkills", package = "betareg")
	
	y <- ReadingSkills$accuracy
	# model accuracy as function of dyslexia and iq
	X_data <- model.matrix(accuracy ~ dyslexia + iq, data = ReadingSkills)
	
	res_cpp <- fast_beta_regression_with_var_cpp(X_data, y)
	res_r <- betareg::betareg(accuracy ~ dyslexia + iq, data = ReadingSkills)
	
	expect_equal(as.numeric(res_cpp$coefficients), as.numeric(stats::coef(res_r)[1:3]), tolerance = 1e-4)
	expect_equal(as.numeric(res_cpp$phi), as.numeric(stats::coef(res_r)[4]), tolerance = 1e-3)
})

test_that("Survival: fast_coxph_regression_cpp on survival::lung", {
	skip_if_not_installed("survival")
	lung_full <- survival::lung
	# Only omit NAs for variables used to maximize censored coverage
	vars <- c("time", "status", "age", "sex", "ph.ecog")
	lung <- na.omit(lung_full[, vars])
	
	y <- lung$time
	dead <- lung$status - 1
	X_data <- lung[, c("age", "sex", "ph.ecog")]
	X_mat <- as.matrix(X_data)
	
	res_cpp <- fast_coxph_regression_cpp(X_mat, y, dead)
	res_r <- survival::coxph(survival::Surv(time, status) ~ age + sex + ph.ecog, data = lung, ties = "breslow")
	
	expect_equal(as.numeric(res_cpp$coefficients), as.numeric(stats::coef(res_r)), tolerance = 1e-7)
	expect_equal(as.numeric(res_cpp$vcov), as.numeric(stats::vcov(res_r)), tolerance = 1e-7)
})

test_that("Survival: fast_weibull_regression_cpp on survival::lung", {
	skip_if_not_installed("survival")
	lung_full <- survival::lung
	vars <- c("time", "status", "age", "sex", "ph.ecog")
	lung <- na.omit(lung_full[, vars])
	
	y <- lung$time
	dead <- lung$status - 1
	X_mat <- model.matrix(~ age + sex + ph.ecog, data = lung)
	
	res_cpp <- fast_weibull_regression_cpp(X_mat, y, dead)
	res_r <- survival::survreg(survival::Surv(time, status) ~ age + sex + ph.ecog, data = lung, dist = "weibull")
	
	expect_equal(as.numeric(res_cpp$coefficients), as.numeric(stats::coef(res_r)), tolerance = 1e-5)
	expect_equal(as.numeric(res_cpp$log_sigma), as.numeric(log(res_r$scale)), tolerance = 1e-5)
	expect_equal(as.numeric(diag(res_cpp$vcov)[1:4]), as.numeric(diag(stats::vcov(res_r))[1:4]), tolerance = 1e-4)
})

test_that("Ordinal: fast_ordinal_regression_with_var_cpp on ordinal::wine", {
	skip_if_not_installed("ordinal")
	data("wine", package = "ordinal")
	
	y <- as.numeric(wine$rating)
	X_mat <- model.matrix(~ temp + contact, data = wine)[, -1] # No intercept
	
	res_cpp <- fast_ordinal_regression_with_var_cpp(X_mat, y)
	res_r <- ordinal::clm(rating ~ temp + contact, data = wine)
	
	expect_equal(as.numeric(res_cpp$alpha), as.numeric(res_r$alpha), tolerance = 1e-5)
	expect_equal(as.numeric(res_cpp$b), as.numeric(res_r$beta), tolerance = 1e-5)
	expect_equal(as.numeric(res_cpp$vcov), as.numeric(stats::vcov(res_r)), tolerance = 1e-5)
})

test_that("Diverse (n, p) Configurations: OLS large n, medium p", {
	set.seed(100)
	n <- 10000
	p <- 20
	X <- matrix(rnorm(n * p), n, p)
	beta <- rnorm(p)
	y <- X %*% beta + rnorm(n)
	
	res_cpp <- fast_ols_with_var_cpp(X, y, j = 1L)
	res_r <- stats::lm.fit(X, y)
	
	expect_equal(as.numeric(res_cpp$b), as.numeric(stats::coef(res_r)), tolerance = 1e-10)
})

test_that("Diverse (n, p) Configurations: Logistic small n, large p", {
	set.seed(101)
	n <- 50
	p <- 5
	X <- matrix(rnorm(n * p), n, p)
	y <- rbinom(n, 1, plogis(X %*% rnorm(p)))
	
	X_fit <- cbind(1, X)
	res_cpp <- fast_logistic_regression_with_var_cpp(X_fit, y)
	res_r <- stats::glm(y ~ X, family = binomial)
	
	expect_equal(as.numeric(res_cpp$b), as.numeric(stats::coef(res_r)), tolerance = 1e-5)
})

test_that("Diverse (n, p) Configurations: Poisson medium n, large p", {
	set.seed(102)
	n <- 500
	p <- 30
	X <- matrix(rnorm(n * p), n, p)
	beta <- rnorm(p) * 0.1
	y <- rpois(n, exp(0.5 + X %*% beta))
	
	X_fit <- cbind(1, X)
	res_cpp <- fast_poisson_regression_with_var_cpp(X_fit, y)
	res_r <- stats::glm(y ~ X, family = poisson)
	
	expect_equal(as.numeric(res_cpp$b), as.numeric(stats::coef(res_r)), tolerance = 1e-6)
})

test_that("Diverse (n, p) Configurations: Cox large n, medium p", {
	skip_if_not_installed("survival")
	set.seed(103)
	n <- 2000
	p <- 15
	X <- matrix(rnorm(n * p), n, p)
	beta <- rnorm(p) * 0.2
	y <- rexp(n, 0.1 * exp(X %*% beta))
	dead <- rbinom(n, 1, 0.7)
	
	res_cpp <- fast_coxph_regression_cpp(X, y, dead)
	res_r <- survival::coxph(survival::Surv(y, dead) ~ X, ties = "breslow")
	
	expect_equal(as.numeric(res_cpp$coefficients), as.numeric(stats::coef(res_r)), tolerance = 1e-7)
})

test_that("Diverse (n, p) Configurations: Ordinal medium n, medium p", {
	skip_if_not_installed("ordinal")
	set.seed(104)
	n <- 800
	p <- 10
	X <- matrix(rnorm(n * p), n, p)
	beta <- rnorm(p) * 0.5
	eta <- X %*% beta
	y_cont <- eta + rlogis(n)
	y <- cut(y_cont, breaks = c(-Inf, -1, 0, 1, Inf), labels = FALSE)
	y <- as.numeric(y)
	
	res_cpp <- fast_ordinal_regression_with_var_cpp(X, y)
	res_r <- ordinal::clm(factor(y) ~ X, link = "logit")
	
	expect_equal(as.numeric(res_cpp$b), as.numeric(res_r$beta), tolerance = 1e-5)
	expect_equal(as.numeric(res_cpp$alpha), as.numeric(res_r$alpha), tolerance = 1e-5)
})

test_that("Robust: fast_robust_regression_cpp on mtcars", {
	skip_if_not_installed("MASS")
	data("mtcars")
	
	y <- mtcars$mpg
	X_fit <- model.matrix(mpg ~ wt + hp + qsec, data = mtcars)
	
	res_cpp <- EDI:::fast_robust_regression_cpp(X_fit, y, method = "MM")
	res_r <- MASS::rlm(mpg ~ wt + hp + qsec, data = mtcars, method = "MM")
	summ_r <- summary(res_r)
	
	# Robust regression is an iterative process and can have multiple local minima
	# or different convergence paths. Using a looser tolerance.
	expect_equal(as.numeric(res_cpp$coefficients), as.numeric(stats::coef(res_r)), tolerance = 1e-1)
	expect_equal(as.numeric(sqrt(res_cpp$ssq_b_j)), as.numeric(summ_r$coefficients[2, "Std. Error"]), tolerance = 1e-1)
})

test_that("Count: fast_zinb_cpp on glmmTMB::Salamanders", {
	skip_if_not_installed("glmmTMB")
	data("Salamanders", package = "glmmTMB")

	X_cond <- model.matrix(~ mined + spp, data = glmmTMB::Salamanders)
	X_zi <- model.matrix(~ mined, data = glmmTMB::Salamanders)
	y <- glmmTMB::Salamanders$count

	res_cpp <- EDI:::fast_zinb_cpp(X = X_cond, Xzi = X_zi, y = y)
	res_r <- glmmTMB::glmmTMB(count ~ mined + spp, ziformula = ~ mined, family = glmmTMB::nbinom2, data = glmmTMB::Salamanders)

	expect_equal(as.numeric(res_cpp$coefficients$cond), as.numeric(glmmTMB::fixef(res_r)$cond), tolerance = 5e-3)
	expect_equal(as.numeric(res_cpp$coefficients$zi), as.numeric(glmmTMB::fixef(res_r)$zi), tolerance = 5e-3)
	expect_equal(as.numeric(res_cpp$params[length(res_cpp$params)]), as.numeric(log(glmmTMB::sigma(res_r))), tolerance = 5e-3)
	expect_true(all(is.finite(diag(res_cpp$vcov))))
})

test_that("Count: fast_zero_augmented_poisson_cpp on glmmTMB::Salamanders", {
	skip_if_not_installed("glmmTMB")
	data("Salamanders", package = "glmmTMB")

	X_cond <- model.matrix(~ mined + spp, data = glmmTMB::Salamanders)
	X_zi <- model.matrix(~ mined, data = glmmTMB::Salamanders)
	y <- glmmTMB::Salamanders$count

	res_cpp <- EDI:::fast_zero_augmented_poisson_cpp(X_cond, y, X_zi, is_hurdle = TRUE)
	res_r <- glmmTMB::glmmTMB(count ~ mined + spp, ziformula = ~ mined, family = glmmTMB::truncated_poisson, data = glmmTMB::Salamanders)
	vcov_r <- stats::vcov(res_r)

	expect_equal(as.numeric(res_cpp$coefficients$cond), as.numeric(glmmTMB::fixef(res_r)$cond), tolerance = 5e-3)
	expect_equal(as.numeric(res_cpp$coefficients$zi), as.numeric(glmmTMB::fixef(res_r)$zi), tolerance = 5e-3)
	expect_equal(as.numeric(diag(res_cpp$vcov)[seq_along(res_cpp$coefficients$cond)]), as.numeric(diag(vcov_r$cond)), tolerance = 5e-3)
	zi_idx <- length(res_cpp$coefficients$cond) + seq_along(res_cpp$coefficients$zi)
	expect_equal(as.numeric(diag(res_cpp$vcov)[zi_idx]), as.numeric(diag(vcov_r$zi)), tolerance = 5e-3)
})
