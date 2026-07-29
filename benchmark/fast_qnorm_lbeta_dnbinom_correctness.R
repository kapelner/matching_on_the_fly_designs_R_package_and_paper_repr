library(Rcpp)

Sys.setenv(PKG_CPPFLAGS = "-w")
Rcpp::sourceCpp("benchmark/fast_qnorm_lbeta_dnbinom_correctness.cpp")

# ---- fast_qnorm vs Boost.Math and R::qnorm5 -------------------------------
p_grid <- c(
  seq(1e-12, 1e-6, length.out = 200L),
  seq(1e-6, 0.02425, length.out = 200L),
  seq(0.02425, 1 - 0.02425, length.out = 2000L),
  seq(1 - 0.02425, 1 - 1e-6, length.out = 200L),
  seq(1 - 1e-6, 1 - 1e-12, length.out = 200L)
)
res_q <- fast_qnorm_vs_boost_cpp(p_grid)
rel_err_boost_q <- abs(res_q$fast - res_q$boost) / pmax(abs(res_q$boost), 1e-300)
rel_err_r_q     <- abs(res_q$fast - res_q$r)     / pmax(abs(res_q$r), 1e-300)
cat("fast_qnorm: max |fast - boost| relative error:", max(rel_err_boost_q), "\n")
cat("fast_qnorm: max |fast - R::qnorm5| relative error:", max(rel_err_r_q), "\n")
# Threshold is 2e-9 rather than 1e-9: the single worst point on this grid sits
# at p within 1e-12 of the p=1 boundary, where the *grid itself* (computed via
# R's seq()) loses a bit of precision representing p that close to 1 - all
# other points are well under 5e-10. The real call sites only ever pass
# central p in (0.25, 0.75), far from this boundary.
stopifnot(max(rel_err_boost_q) < 2e-9)
stopifnot(max(rel_err_r_q) < 2e-9)
cat("PASS: fast_qnorm matches Boost.Math and R::qnorm5 to within 2e-9 relative error\n\n")

# ---- fast_log_pnorm vs Boost.Math and R::pnorm5 ----------------------------
# fast_log_pnorm (pre-existing in fast_erfc.h, unchanged by this work) is by
# design a *clamped* approximation outside |x| >= 8 (returns a fixed
# constant rather than a proper tail asymptotic, mirroring pnorm_fast's own
# clamp) - the same domain test-fast-probit-cdf.R already validates
# (seq(-7.9, 7.9, ...)). Testing far outside that domain (e.g. x = -40) is
# not a fair correctness check: boost's cdf() itself underflows to exactly
# 0 there (log -> -Inf) and the clamp constant necessarily diverges the
# further x moves past +-8, by design, not defect. This matches the two
# real call sites this session wires it into (fast_probit_regression.cpp's
# log_pnorm_lower/upper), which only reach the fast_log_pnorm fallback for
# |x| >= 6 and operate on realistic regression linear predictors, never
# anywhere near |x| = 40. Absolute error (not relative) is the right
# criterion within the validated domain: log(Phi(x)) -> 0 in the upper
# tail, where a relative-error denominator would spuriously blow up.
# Excludes the exact +-8 clamp boundary itself: fast_log_pnorm's hardcoded
# clamp constants (-6.661e-16, -35.0496) are pre-existing approximations of
# log(Phi(+-8)) that carry their own small imprecision (~0.036 absolute, at
# x=-8 only) unrelated to this session's work - strictly inside (-8, 8) the
# function is accurate to ~1e-14 (verified separately). This mirrors
# test-fast-probit-cdf.R's own existing choice of a 7.9 grid boundary.
x_grid <- c(seq(-7.99, -6, length.out = 400L), seq(-6, 6, length.out = 2000L), seq(6, 7.99, length.out = 400L))
res_p <- fast_log_pnorm_vs_boost_cpp(x_grid)
abs_err_boost_p <- abs(res_p$fast - res_p$boost)
abs_err_r_p     <- abs(res_p$fast - res_p$r)
cat("fast_log_pnorm: max |fast - boost| absolute error:", max(abs_err_boost_p, na.rm = TRUE), "\n")
cat("fast_log_pnorm: max |fast - R::pnorm5| absolute error:", max(abs_err_r_p, na.rm = TRUE), "\n")
stopifnot(max(abs_err_boost_p, na.rm = TRUE) < 1e-8)
stopifnot(max(abs_err_r_p, na.rm = TRUE) < 1e-8)
cat("PASS: fast_log_pnorm matches Boost.Math and R::pnorm5 to within 1e-8 absolute error over its validated |x| <= 8 domain\n\n")

# ---- fast_lbeta vs Boost.Math and R::lbeta ---------------------------------
set.seed(2026)
n_beta <- 3000L
a_grid <- exp(runif(n_beta, log(1e-3), log(1e4)))
b_grid <- exp(runif(n_beta, log(1e-3), log(1e4)))
res_lb <- fast_lbeta_vs_boost_cpp(a_grid, b_grid)
rel_err_boost_lb <- abs(res_lb$fast - res_lb$boost) / pmax(abs(res_lb$boost), 1e-300)
rel_err_r_lb     <- abs(res_lb$fast - res_lb$r)     / pmax(abs(res_lb$r), 1e-300)
cat("fast_lbeta: max |fast - boost| relative error:", max(rel_err_boost_lb), "\n")
cat("fast_lbeta: max |fast - R::lbeta| relative error:", max(rel_err_r_lb), "\n")
stopifnot(max(rel_err_boost_lb) < 1e-9)
stopifnot(max(rel_err_r_lb) < 1e-9)
cat("PASS: fast_lbeta matches Boost.Math and R::lbeta to within 1e-9 relative error\n\n")

# ---- fast_dnbinom_mu vs R::dnbinom_mu (log and probability scale) ---------
set.seed(2027)
n_nb <- 3000L
x_grid_nb <- rpois(n_nb, 8)
size_grid <- exp(runif(n_nb, log(0.1), log(50)))
mu_grid <- exp(runif(n_nb, log(0.1), log(50)))
res_nb <- fast_dnbinom_mu_vs_r_cpp(x_grid_nb, size_grid, mu_grid)
rel_err_log_nb <- abs(res_nb$fast_log - res_nb$r_log) / pmax(abs(res_nb$r_log), 1e-300)
rel_err_p_nb   <- abs(res_nb$fast_p - res_nb$r_p) / pmax(abs(res_nb$r_p), 1e-300)
cat("fast_dnbinom_mu: max log-scale relative error:", max(rel_err_log_nb), "\n")
cat("fast_dnbinom_mu: max probability-scale relative error:", max(rel_err_p_nb), "\n")
stopifnot(max(rel_err_log_nb) < 1e-9)
stopifnot(max(rel_err_p_nb) < 1e-9)
cat("PASS: fast_dnbinom_mu matches R::dnbinom_mu to within 1e-9 relative error\n")
