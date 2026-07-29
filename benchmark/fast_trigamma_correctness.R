library(Rcpp)

Sys.setenv(PKG_CPPFLAGS = "-w")
Rcpp::sourceCpp("benchmark/fast_trigamma_correctness.cpp")

grid <- c(
  seq(1e-6, 0.5, length.out = 200L),
  seq(0.5, 8, length.out = 400L),
  seq(8, 1e4, length.out = 400L)
)

res <- fast_trigamma_vs_boost_cpp(grid)

rel_err_boost <- abs(res$fast - res$boost) / pmax(abs(res$boost), 1e-300)
rel_err_r     <- abs(res$fast - res$r)     / pmax(abs(res$r), 1e-300)

cat("max |fast - boost| relative error:", max(rel_err_boost), "\n")
cat("max |fast - R::trigamma| relative error:", max(rel_err_r), "\n")

stopifnot(max(rel_err_boost) < 1e-10)
stopifnot(max(rel_err_r) < 1e-10)

cat("PASS: fast_trigamma matches Boost.Math and R::trigamma to within 1e-10 relative error\n")

vec_grid <- seq(0.5, 500, length.out = 5000L)
vec_out <- fast_trigamma_vec_vs_scalar_cpp(vec_grid)
scalar_out <- vapply(vec_grid, function(v) fast_trigamma_vs_boost_cpp(v)$fast, numeric(1L))
rel_err_vec <- abs(vec_out - scalar_out) / pmax(abs(scalar_out), 1e-300)
cat("max |fast_trigamma_vec - fast_trigamma| relative error:", max(rel_err_vec), "\n")
stopifnot(max(rel_err_vec) < 1e-14)
cat("PASS: fast_trigamma_vec matches fast_trigamma exactly across the grid\n")
