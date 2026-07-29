library(Rcpp)

Sys.setenv(PKG_CPPFLAGS = "-w")
Rcpp::sourceCpp("benchmark/fast_pchisq_correctness.cpp")

# Grid covering the real call-site range (df = 1, 2, 3, ... small integers
# from LRT/score tests) plus a wider df range for general robustness, and
# statistic spanning near-zero (large p-value) through the extreme tail
# (tiny p-value, where relative error matters most since these p-values
# feed directly into significance decisions).
set.seed(2026)
n <- 4000L
df_grid <- c(
  sample(1:10, n * 0.5, replace = TRUE),                # realistic small-integer df
  exp(runif(n * 0.5, log(0.5), log(100)))                # wider general df
)
stat_grid <- c(
  runif(length(df_grid) * 0.3, 0, 5),                    # central region
  exp(runif(length(df_grid) * 0.5, log(1e-6), log(50))), # log-spaced, covers tiny stats
  exp(runif(length(df_grid) * 0.2, log(50), log(1e4)))   # extreme tail (tiny p-values)
)
df_grid <- rep(df_grid, length.out = length(stat_grid))

res <- fast_pchisq_vs_boost_cpp(stat_grid, df_grid)

rel_err_boost <- abs(res$fast - res$boost) / pmax(abs(res$boost), 1e-300)
rel_err_r     <- abs(res$fast - res$r)     / pmax(abs(res$r), 1e-300)

cat("fast_pchisq_upper: max |fast - boost| relative error:", max(rel_err_boost), "\n")
cat("fast_pchisq_upper: max |fast - R::pchisq| relative error:", max(rel_err_r), "\n")

worst_idx <- order(rel_err_r, decreasing = TRUE)[1:5]
cat("Worst 5 (vs R::pchisq):\n")
print(data.frame(statistic = res$statistic[worst_idx], df = res$df[worst_idx],
                  fast = res$fast[worst_idx], r = res$r[worst_idx], rel_err = rel_err_r[worst_idx]))

stopifnot(max(rel_err_boost) < 1e-9)
stopifnot(max(rel_err_r) < 1e-9)
cat("PASS: fast_pchisq_upper matches Boost.Math and R::pchisq to within 1e-9 relative error\n")

# Specific check on the exact real call-site pattern: df = 1 (both score-test
# sites) and small integer df (LRT site).
for (d in c(1, 2, 3, 5)) {
  stats <- exp(seq(log(1e-8), log(1e3), length.out = 500))
  r <- fast_pchisq_vs_boost_cpp(stats, rep(d, length(stats)))
  err <- max(abs(r$fast - r$r) / pmax(abs(r$r), 1e-300))
  cat(sprintf("df=%d: max relative error vs R::pchisq = %.3e\n", d, err))
  stopifnot(err < 1e-9)
}
cat("PASS: real call-site df values (1, 2, 3, 5) all within tolerance\n")
