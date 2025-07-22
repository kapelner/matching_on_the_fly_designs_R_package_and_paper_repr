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

mod = fast_ols_with_sd_cpp(cbind(1, X), y); abs(mod$b[2] / mod$s_b_2)
mod$b
mod = lm(y ~ X); abs(coef(summary(mod))[2, 3])
mod$coefficients
microbenchmark(
  Rcpp_opt = {mod = fast_ols_cpp(cbind(1, X), y); abs(mod$b[2] / mod$s_b[2])},
  Rcpp = {mod = fast_ols_with_sd_cpp(cbind(1, X), y); abs(mod$b[2] / mod$s_b_2)},
  Rcppeigen = {mod = fastLm(cbind(1, X), y); mod$coefficients[2] / mod$stderr[2]},
  R = {mod = lm(y ~ X); abs(coef(summary(mod))[2, 3])},
  times = 100
)
rm(list = ls())


X = as.matrix(MASS::Boston[, 1:13])
y = as.numeric(MASS::Boston$medv > median(MASS::Boston$medv))
Rcpp::sourceCpp("../SeqExpMatch/src/fast_logistic_regression.cpp")

# mod = fast_logistic_regression_cpp(cbind(1, X), y, start = rep(0, ncol(X) + 1))
# mod$b
mod = glm(y ~ X, family="binomial")
mod$coefficients
mod = glm.fit(cbind(1, X), y, family=binomial())
mod$coefficients
W <- diag(mod$weights)
XtWX <- eigen_Xt_times_diag_w_times_X_cpp(cbind(1, X), mod$weights)
XtWX_inv <- solve(t(cbind(1, X)) %*% W %*% cbind(1, X))
sqrt(diag(XtWX_inv))[2]
sqrt(eigen_compute_single_entry_of_diagonal_matrix(XtWX, 2))


microbenchmark(
  Rcpp = {mod = fast_logistic_regression_cpp(cbind(1, X), y, start = rep(0, ncol(X) + 1)); abs(mod$b[1] / mod$s_b[1])},
  fastLogisticRegressionWrap = {mod = fastLogisticRegressionWrap::fast_logistic_regression(cbind(1, X), y)},
  R = {mod = glm(y ~ X, family="binomial"); abs(coef(summary(mod))[2, 3])},
  Ropt = {mod = glm.fit(cbind(1, X), y, family=binomial()); XtWX <- eigen_Xt_times_diag_w_times_X_cpp(cbind(1, X), mod$weights); sqrt(eigen_compute_single_entry_of_diagonal_matrix(XtWX, 2))},
  times = 100
)
rm(list = ls())





