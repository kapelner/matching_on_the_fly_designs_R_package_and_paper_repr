pacman::p_load(Rcpp, RcppEigen, microbenchmark, fastLogisticRegressionWrap)


sourceCpp(file = "../SeqExpMatch/src/fast_matrix_rank.cpp")

X = matrix(rnorm(100^2), nrow = 100, ncol = 100)

microbenchmark(
  rcpp = matrix_rank_cpp(X),
  R = Matrix::rankMatrix(X),
  times = 100
)
rm(list = ls())


X = matrix(rnorm(10000*100), nrow = 100, ncol = 10000)

Rcpp::sourceCpp("../SeqExpMatch/src/fast_scale_cols.cpp")
microbenchmark(
  Rcpp = scale_columns_cpp(X),
  R = apply(X, 2, scale),
  times=10
)
rm(list = ls())

X = as.matrix(MASS::Boston[, 1:13])
y = MASS::Boston$medv
Rcpp::sourceCpp("../SeqExpMatch/src/fast_ols.cpp")

mod = fast_ols_with_var_cpp(cbind(1, X), y); abs(mod$b[2] / sqrt(mod$ssq_b_j))
mod$b
mod = lm(y ~ X); abs(coef(summary(mod))[2, 3])
mod$coefficients
microbenchmark(
  Rcpp_opt = {mod = fast_ols_cpp(cbind(1, X), y); abs(mod$b[2] / mod$s_b[2])},
  Rcpp = {mod = fast_ols_with_var_cpp(cbind(1, X), y); abs(mod$b[2] / sqrt(mod$ssq_b_j))},
  Rcppeigen = {mod = fastLm(cbind(1, X), y); mod$coefficients[2] / mod$stderr[2]},
  R = {mod = lm(y ~ X); abs(coef(summary(mod))[2, 3])},
  times = 100
)
rm(list = ls())

Bsub = MASS::Boston[1:200, ]
X = as.matrix(Bsub[, 1:13])
y = as.numeric(Bsub$medv > median(Bsub$medv))
Rcpp::sourceCpp("../SeqExpMatch/src/_helper_functions.cpp")
# Rcpp::sourceCpp("../SeqExpMatch/src/fast_logistic_regression.cpp")

# mod = fast_logistic_regression_cpp(cbind(1, X), y, start = rep(0, ncol(X) + 1))
# mod$b
mod = glm(y ~ X, family="binomial")
mod$coefficients
mod = glm.fit(cbind(1, X), y, family=binomial())
mod$coefficients
# W <- diag(mod$weights)
# XtWX <- eigen_Xt_times_diag_w_times_X_cpp(cbind(1, X), mod$weights)
# XtWX_inv <- solve(t(cbind(1, X)) %*% W %*% cbind(1, X))
# sqrt(diag(XtWX_inv))[2]
# sqrt(eigen_compute_single_entry_of_diagonal_matrix(XtWX, 2))

pacman::p_load(glmnet, speedglm)
mod <- speedglm(y ~ X, family = binomial())
mod$coefficients

mod = glmnet(
  x = X,
  y = y,
  family = "binomial",
  alpha = 0,      # ridge penalty (alpha=1 would be lasso)
  lambda = 0,     # no penalty → close to unpenalized MLE
  standardize = FALSE,
  intercept = TRUE
)
mod$beta

microbenchmark::microbenchmark(
  # Rcpp = {mod = fast_logistic_regression_cpp(cbind(1, X), y, start = rep(0, ncol(X) + 1))},
  # fastLogisticRegressionWrap = {mod = fastLogisticRegressionWrap::fast_logistic_regression(cbind(1, X), y)},
  Rglm = {mod = glm(y ~ X, family="binomial")},
  Rglmfit = {mod = glm.fit(cbind(1, X), y, family=binomial())},
  Rspeedglm = {mod = speedglm(y ~ X, family = binomial())},
  Rglmnet = {mod = glmnet(
    x = X,
    y = y,
    family = "binomial",
    alpha = 0,      # ridge penalty (alpha=1 would be lasso)
    lambda = 0,     # no penalty → close to unpenalized MLE
    standardize = FALSE,
    intercept = TRUE
  )
  mod$beta},
  Rglmnetopt = {mod = glmnet(
    x = X,
    y = y,
    family = "binomial",
    alpha = 0,      # ridge penalty (alpha=1 would be lasso)
    lambda = 0,     # no penalty → close to unpenalized MLE
    standardize = FALSE,
    intercept = TRUE,
    type.logistic = "modified.Newton"
  )
  mod$beta},
  times = 500
)

microbenchmark::microbenchmark(
  # Rcpp = {mod = fast_logistic_regression_cpp(cbind(1, X), y, start = rep(0, ncol(X) + 1)); abs(mod$b[1] / mod$s_b[1])},
  fastLogisticRegressionWrap = {mod = fastLogisticRegressionWrap::fast_logistic_regression(cbind(1, X), y, do_inference_on_var = 2); mod$se[2]},
  Rglm = {mod = glm(y ~ X, family="binomial"); abs(coef(summary(mod))[2, 3])},
  Rglmfit = {mod = glm.fit(cbind(1, X), y, family=binomial()); XtWX <- eigen_Xt_times_diag_w_times_X_cpp(cbind(1, X), mod$weights); sqrt(eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(XtWX, 2))},
  Rspeedglm = {mod = speedglm(y ~ X, family = binomial()); coef(summary(mod))[2,2]},
  times = 100
)
rm(list = ls())


pacman::p_load(Rcpp, RcppEigen, RcppNumerical, RcppParallel, StanHeaders, microbenchmark, glmmTMB)

X = as.matrix(MASS::Boston[MASS::Boston$medv < 50, 1:13])
y = as.integer(round(MASS::Boston$medv[MASS::Boston$medv < 50])) #count response

# Rcpp::sourceCpp("../SeqExpMatch/src/log_lik_nb.cpp")
# Rcpp::sourceCpp("../SeqExpMatch/src/fast_negbin_regression_stan.cpp")
mod = MASS::glm.nb(y ~ X)
# source("../SeqExpMatch/R/model_fit_helpers.R")
modfast = suppressWarnings(glmmTMB(y ~ X, family = nbinom2))
coef(summary(modfast))$cond[2, 2]

microbenchmark::microbenchmark(
  R_fast = {mod = glmmTMB(y ~ X, family = nbinom2); coef(summary(modfast))$cond[2, 1]},
  R = {mod = MASS::glm.nb(y ~ X); coef(mod)[2]},
  times = 100
)

microbenchmark::microbenchmark(
  R_fast = {mod = fast_glm_nb(cbind(1, X), y); mod$b[2]},
  R = {mod = MASS::glm.nb(y ~ X); abs(coef(summary(mod))[2, 3])},
  times = 100
)

pacman::p_load(betareg, microbenchmark)


source(file = "../SeqExpMatch/R/model_fit_helpers.R")
Rcpp::sourceCpp("../SeqExpMatch/src/_helper_functions.cpp")
Rcpp::sourceCpp("../SeqExpMatch/src/beta_regression_helpers.cpp")

X = as.matrix(MASS::Boston[, 1:13])
y = (MASS::Boston$medv - min(MASS::Boston$medv) + 0.01) / (max(MASS::Boston$medv) - min(MASS::Boston$medv) + 0.02)

mod = betareg(y ~ ., data.frame(X))
coef(mod)
coef(summary(mod))$mean[2,2]

mod = fast_beta_regression(Xmm = cbind(1, X), y = y)
mod$b
mod = fast_beta_regression_with_var(Xmm = cbind(1, X), y = y)
sqrt(mod$ssq_b_2)

profvis::profvis({
  for (i in 1 : 100){
    fast_beta_regression(Xmm = cbind(1, X), y = y)$b[1]
  }
})
microbenchmark::microbenchmark(
  R = {mod = betareg(y ~ ., data.frame(X)); coef(mod)[2]},
  R_fast = {fast_beta_regression(Xmm = cbind(1, X), y = y)$b[1]},
  times = 50
)

microbenchmark::microbenchmark(
  R = {mod = betareg(y ~ ., data.frame(X)); coef(summary(mod))$mean[2,2]},
  R_fast = {mod = fast_beta_regression_with_var(Xmm = cbind(1, X), y = y); sqrt(mod$ssq_b_2)},
  times = 50
)
rm(list = ls())



