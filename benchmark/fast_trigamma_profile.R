library(Rcpp)

Sys.setenv(PKG_CPPFLAGS = "-w")
Rcpp::sourceCpp("benchmark/fast_trigamma_profile.cpp", cacheDir = "/tmp/fast_trigamma_profile_cache")

set.seed(1)
n <- 2000L
data_kind <- commandArgs(trailingOnly = TRUE)
data_kind <- if (length(data_kind) >= 2) data_kind[2] else "wide"
if (identical(data_kind, "realistic")) {
  # mimics y + theta at real negbin/beta-regression call sites: y ~ Poisson(8), theta in [2, 12]
  x <- rpois(n, 8) + runif(n, 2, 12)
} else {
  x <- runif(n, 0.01, 50)
}

which_fn <- commandArgs(trailingOnly = TRUE)
which_fn <- if (length(which_fn) >= 1) which_fn[1] else "scalar"
reps <- 200000L

if (which_fn == "scalar") {
  invisible(bench_fast_trigamma_scalar_loop_cpp(x, reps))
} else {
  invisible(bench_fast_trigamma_vec_loop_cpp(x, reps))
}
