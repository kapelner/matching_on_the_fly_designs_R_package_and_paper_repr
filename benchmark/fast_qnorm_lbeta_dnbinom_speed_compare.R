library(Rcpp)

Sys.setenv(PKG_CPPFLAGS = "-w")
Rcpp::sourceCpp("benchmark/fast_qnorm_lbeta_dnbinom_speed_compare.cpp", cacheDir = "/tmp/fast_qnorm_lbeta_dnbinom_speed_cache")

run_paired_bench <- function(name, fast_fn, r_fn, n_rounds = 30L) {
  set.seed(2026)
  order_flags <- rep(c(TRUE, FALSE), length.out = n_rounds)
  fast_ms <- numeric(n_rounds)
  r_ms <- numeric(n_rounds)
  invisible(fast_fn()); invisible(r_fn())  # warm-up, discarded
  for (round in seq_len(n_rounds)) {
    if (order_flags[round]) {
      fast_ms[round] <- fast_fn()
      r_ms[round]    <- r_fn()
    } else {
      r_ms[round]    <- r_fn()
      fast_ms[round] <- fast_fn()
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

  cat(sprintf("== %s ==\n", name))
  cat(sprintf("fast median: %.4f ms (IQR %.4f-%.4f)\n", median_fast, quantile(fast_ms, 0.25), quantile(fast_ms, 0.75)))
  cat(sprintf("R::   median: %.4f ms (IQR %.4f-%.4f)\n", median_r, quantile(r_ms, 0.25), quantile(r_ms, 0.75)))
  cat(sprintf("Median paired speedup: %.3fx (bootstrap 95%% CI %.3fx-%.3fx)\n", median_speedup, ci[1], ci[2]))
  cat(sprintf("Paired one-sided Wilcoxon p-value (R:: slower): %.3g\n", wtest$p.value))
  stopifnot(median_speedup > 1)
  stopifnot(wtest$p.value < 0.05)
  cat(sprintf("PASS: %s is significantly faster than its R:: counterpart\n\n", name))
  invisible(list(median_fast = median_fast, median_r = median_r, speedup = median_speedup))
}

n <- 500L
reps_per_call <- 50L

# qnorm: draw p from a realistic mix (mostly central, per the two real call sites)
p_x <- runif(n, 0.02, 0.98)
run_paired_bench("fast_qnorm vs R::qnorm5",
  function() time_fast_qnorm_cpp(p_x, reps_per_call),
  function() time_r_qnorm_cpp(p_x, reps_per_call))

# log_pnorm: draw x spanning the central + moderate-tail region
x_x <- rnorm(n, 0, 3)
run_paired_bench("fast_log_pnorm vs R::pnorm5(log=TRUE)",
  function() time_fast_log_pnorm_cpp(x_x, reps_per_call),
  function() time_r_pnorm_cpp(x_x, reps_per_call))

# lbeta: draw mu*phi / (1-mu)*phi style arguments
a_x <- exp(runif(n, log(0.5), log(200)))
b_x <- exp(runif(n, log(0.5), log(200)))
run_paired_bench("fast_lbeta vs R::lbeta",
  function() time_fast_lbeta_cpp(a_x, b_x, reps_per_call),
  function() time_r_lbeta_cpp(a_x, b_x, reps_per_call))

# dnbinom_mu: draw x/size/mu in a realistic hurdle-negbin truncation-loop range
x_nb <- rpois(n, 6)
size_nb <- exp(runif(n, log(0.5), log(20)))
mu_nb <- exp(runif(n, log(0.5), log(20)))
run_paired_bench("fast_dnbinom_mu vs R::dnbinom_mu",
  function() time_fast_dnbinom_mu_cpp(x_nb, size_nb, mu_nb, reps_per_call),
  function() time_r_dnbinom_mu_cpp(x_nb, size_nb, mu_nb, reps_per_call))

cat("ALL PASS: fast_qnorm, fast_log_pnorm, fast_lbeta, fast_dnbinom_mu are all significantly faster than their R:: counterparts\n")
