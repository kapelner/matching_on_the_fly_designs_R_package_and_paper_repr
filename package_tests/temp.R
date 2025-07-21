X = matrix(rnorm(100^2), nrow = 100, ncol = 100)

pacman::p_load(Rcpp, RcppEigen, microbenchmark)
sourceCpp(file = "../SeqExpMatch/src/fast_matrix_rank.cpp")

microbenchmark::microbenchmark(
  rcpp = matrix_rank_cpp(X),
  R = Matrix::rankMatrix(X),
  times = 100
)