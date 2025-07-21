X = matrix(rnorm(10000*100), nrow = 100, ncol = 10000)

Rcpp::sourceCpp("../SeqExpMatch/src/fast_scale_cols.cpp")
microbenchmark::microbenchmark(
  Rcpp = scale_columns_cpp(X),
  R = apply(X, 2, scale),
  times=10
)
