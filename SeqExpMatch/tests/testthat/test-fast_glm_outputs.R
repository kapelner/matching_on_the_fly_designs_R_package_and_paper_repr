context("Test fast GLM implementations against canonical R packages")

test_that("fast_logistic_regression_with_var_cpp matches glm", {
  set.seed(1)
  n <- 100
  p <- 2
  X <- matrix(rnorm(n * p), n, p)
  X <- cbind(1, X)
  beta_true <- c(0.5, 1.5, -1)
  linear_predictor <- X %*% beta_true
  prob <- 1 / (1 + exp(-linear_predictor))
  y <- rbinom(n, 1, prob)
  
  # fast version
  fast_mod <- SeqExpMatch:::fast_logistic_regression_with_var_cpp(X, y)
  
  # canonical version
  canon_mod <- glm(y ~ X - 1, family = binomial(link = "logit"))
  
  # Check coefficients
  expect_equal(as.numeric(fast_mod$b), as.numeric(coef(canon_mod)), tolerance = 1e-4)
  
  # Check standard error for treatment (second column)
  canon_ses = summary(canon_mod)$coefficients[, "Std. Error"]
  expect_equal(as.numeric(sqrt(fast_mod$ssq_b_2)), as.numeric(canon_ses[2]), tolerance = 1e-4)
})

test_that("fast_weibull_regression matches survreg", {
  set.seed(1)
  n <- 300 # increased for better convergence
  p <- 2
  X <- matrix(rnorm(n * p), n, p)
  
  true_intercept = 1
  true_betas = c(0.5, -0.5)
  eta <- true_intercept + X %*% true_betas
  true_log_sigma = 0.2
  
  y = rweibull(n, shape = 1/exp(true_log_sigma), scale = exp(eta))
  dead <- rbinom(n, 1, 0.9)
  
  # fast version (expects X without intercept)
  fast_mod <- SeqExpMatch:::fast_weibull_regression(y, dead, X)
  
  # canonical version
  canon_mod <- survival::survreg(survival::Surv(y, dead) ~ X, dist = "weibull")
  
  # Check coefficients
  expect_equal(as.numeric(fast_mod$coefficients), as.numeric(coef(canon_mod)), tolerance = 0.1)
  
  # Check log(sigma)
  expect_equal(fast_mod$log_sigma, log(canon_mod$scale), tolerance = 0.1)
  
  # Check standard errors
  expect_equal(as.numeric(fast_mod$std_errs), as.numeric(summary(canon_mod)$table[, "Std. Error"]), tolerance = 0.1)
})

test_that("fast_beta_regression_mle matches betareg", {
  if (!requireNamespace("betareg", quietly = TRUE)) {
    skip("betareg not installed")
  }
  set.seed(123)
  n <- 200
  p <- 2
  X <- matrix(rnorm(n * p), n, p)
  X_int <- cbind(1, X)
  beta_true <- c(-0.5, 0.8, -0.8)
  phi_true <- 15
  
  linear_predictor <- X_int %*% beta_true
  mu <- 1 / (1 + exp(-linear_predictor))
  
  y <- rbeta(n, mu * phi_true, (1 - mu) * phi_true)
  y = pmax(pmin(y, 1 - 1e-6), 1e-6)
  
  # fast version
  fast_mod <- SeqExpMatch:::fast_beta_regression_mle(y, X_int, compute_std_errs = TRUE)
  
  # canonical version
  canon_mod <- betareg::betareg(y ~ X)
  
  # Check coefficients
  expect_equal(as.numeric(fast_mod$coefficients), as.numeric(coef(canon_mod)[1:3]), tolerance = 0.05)
  
  # Check phi
  expect_equal(fast_mod$phi, as.numeric(coef(canon_mod)[4]), tolerance = 0.2)
  
  # Check standard errors (only for mean model)
  expect_equal(as.numeric(fast_mod$std_errs[1:3]), as.numeric(summary(canon_mod)$coefficients$mean[, "Std. Error"]), tolerance = 0.1)
})

test_that("fast_neg_bin_with_censoring_with_sd_cpp matches MASS::glm.nb", {
  set.seed(123)
  n <- 200
  p <- 2
  X <- matrix(rnorm(n * p), n, p)
  X_int <- cbind(1, X)
  beta_true <- c(1, 0.5, -0.5)
  mu <- exp(X_int %*% beta_true)
  theta_true <- 2.0
  
  y <- MASS::rnegbin(n, mu = mu, theta = theta_true)
  dead <- rep(1L, n)
  
  fast_mod <- SeqExpMatch:::fast_neg_bin_with_censoring_with_sd_cpp(X_int, y, dead)
  
  # canonical version
  canon_mod <- MASS::glm.nb(y ~ X)
  
  # Check coefficients
  expect_equal(as.numeric(fast_mod$b), as.numeric(coef(canon_mod)), tolerance = 0.05)
  
  # Check standard error for treatment (second column)
  canon_se_w = summary(canon_mod)$coefficients[2, "Std. Error"]
  # Extract ssq_b_2 safely
  ssq_b_2 = if ("ssq_b_2" %in% names(fast_mod)) fast_mod$ssq_b_2 else diag(solve(fast_mod$hess_fisher_info_matrix))[2]
  expect_equal(as.numeric(sqrt(ssq_b_2)), as.numeric(canon_se_w), tolerance = 0.05)
})