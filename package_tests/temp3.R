X = as.matrix(MASS::Boston[, 1:13])
y = MASS::Boston$medv

pacman::p_load(Rcpp)

Rcpp::sourceCpp("../SeqExpMatch/src/fast_ols.cpp")

microbenchmark::microbenchmark(
  Rcpp = {mod = fast_ols_cpp(X, y); abs(mod$b[1] / mod$s_b[1])},
  R = {mod = lm(y ~ X); abs(coef(summary(mod))[2, 3])},
  times = 100
)
