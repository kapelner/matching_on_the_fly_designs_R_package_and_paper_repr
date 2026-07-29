library(Rcpp)

Sys.setenv(PKG_CPPFLAGS = "-w")
Rcpp::sourceCpp("benchmark/fast_pchisq_speed_compare.cpp", cacheDir = "/tmp/fast_pchisq_speed_cache")

set.seed(2026)
n <- 500L
# Realistic LRT/score-test statistic range (chi-squared draws with small df,
# the actual call-site distribution) plus df matching the 3 real call sites
# (df=1 for both score-test sites; small positive integer df for the LRT site).
df <- sample(1:5, n, replace = TRUE)
statistic <- rchisq(n, df = df)
reps_per_call <- 50L
n_rounds <- 30L

fast_ms <- numeric(n_rounds)
r_ms <- numeric(n_rounds)

invisible(time_fast_pchisq_cpp(statistic, df, 5L))
invisible(time_r_pchisq_cpp(statistic, df, 5L))

order_flags <- rep(c(TRUE, FALSE), length.out = n_rounds)
for (round in seq_len(n_rounds)) {
  if (order_flags[round]) {
    fast_ms[round] <- time_fast_pchisq_cpp(statistic, df, reps_per_call)
    r_ms[round]    <- time_r_pchisq_cpp(statistic, df, reps_per_call)
  } else {
    r_ms[round]    <- time_r_pchisq_cpp(statistic, df, reps_per_call)
    fast_ms[round] <- time_fast_pchisq_cpp(statistic, df, reps_per_call)
  }
}

speedup <- r_ms / fast_ms
median_fast <- median(fast_ms)
median_r <- median(r_ms)
median_speedup <- median(speedup)

set.seed(1)
boot_medians <- replicate(5000L, median(sample(speedup, replace = TRUE)))
ci <- quantile(boot_medians, c(0.025, 0.975))

wtest <- wilcox.test(r_ms, fast_ms, paired = TRUE, alternative = "greater")

cat(sprintf("N = %d, reps per timed call = %d, rounds = %d\n", n, reps_per_call, n_rounds))
cat(sprintf("fast_pchisq_upper median: %.4f ms (IQR %.4f-%.4f)\n",
            median_fast, quantile(fast_ms, 0.25), quantile(fast_ms, 0.75)))
cat(sprintf("R::pchisq median:         %.4f ms (IQR %.4f-%.4f)\n",
            median_r, quantile(r_ms, 0.25), quantile(r_ms, 0.75)))
cat(sprintf("Median paired speedup: %.3fx (bootstrap 95%% CI %.3fx-%.3fx)\n",
            median_speedup, ci[1], ci[2]))
cat(sprintf("Paired one-sided Wilcoxon p-value (R::pchisq slower): %.3g\n", wtest$p.value))

stopifnot(median_speedup > 1)
stopifnot(wtest$p.value < 0.05)
cat("PASS: fast_pchisq_upper is significantly faster than R::pchisq\n")
