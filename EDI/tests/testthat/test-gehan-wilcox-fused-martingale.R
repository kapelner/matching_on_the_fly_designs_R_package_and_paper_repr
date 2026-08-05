library(testthat)
library(EDI)

# Verify fast_gehan_wilcox_stats_cpp, the fused Peto-Prentice (rho=1) analog of
# fast_logrank_stats_cpp (TODO: speed up InferenceSurvivalGehanWilcox, which was
# 0.84x vs. canonical R because it called coxph(~1) + residuals + survfit(~1) +
# survdiff(rho=1) as four separate survival:: calls).
# Reference: survival::survdiff(rho=1) (score/var_score) and the canonical
# Peto-Prentice weighted martingale residual mean difference (beta_hat/se_beta_hat).

skip_if_not_installed("survival")

make_surv_data <- function(n, seed) {
    set.seed(seed)
    w    <- as.integer(sample(0:1, n, replace = TRUE))
    y    <- rexp(n, rate = exp(0.3 * w))
    dead <- as.integer(rbinom(n, 1, 0.8))
    list(w = w, y = y, dead = dead)
}

canonical_gehan_wilcox <- function(d) {
    surv_obj  <- survival::Surv(d$y, d$dead)
    cox_null  <- survival::coxph(surv_obj ~ 1)
    M         <- as.numeric(residuals(cox_null, type = "martingale"))
    km_all    <- survival::survfit(surv_obj ~ 1)
    idx       <- findInterval(d$y, km_all$time, left.open = TRUE)
    peto_w    <- c(1.0, km_all$surv)[idx + 1L]
    M_w       <- peto_w * M
    n1 <- sum(d$w == 1); n0 <- sum(d$w == 0)
    beta_ref  <- mean(M_w[d$w == 1]) - mean(M_w[d$w == 0])
    se_ref    <- sqrt(var(M_w[d$w == 1]) / n1 + var(M_w[d$w == 0]) / n0)
    df   <- data.frame(time = d$y, status = d$dead, w = factor(d$w))
    sd   <- survival::survdiff(survival::Surv(time, status) ~ w, data = df, rho = 1)
    list(beta_hat = beta_ref, se_beta_hat = se_ref,
         score = sd$obs[2] - sd$exp[2], var_score = sd$var[2, 2])
}

test_that("gehan-wilcox score and var_score match survival::survdiff(rho=1)", {
    d   <- make_surv_data(120L, seed = 1L)
    res <- EDI:::fast_gehan_wilcox_stats_cpp(d$w, d$y, d$dead)
    ref <- canonical_gehan_wilcox(d)
    expect_equal(res$score,     ref$score,     tolerance = 1e-8)
    expect_equal(res$var_score, ref$var_score, tolerance = 1e-8)
})

test_that("gehan-wilcox beta_hat and se_beta_hat match canonical Peto-Prentice weighted martingale residuals", {
    d   <- make_surv_data(150L, seed = 2L)
    res <- EDI:::fast_gehan_wilcox_stats_cpp(d$w, d$y, d$dead)
    ref <- canonical_gehan_wilcox(d)
    expect_equal(res$beta_hat,    ref$beta_hat,    tolerance = 1e-8)
    expect_equal(res$se_beta_hat, ref$se_beta_hat, tolerance = 1e-8)
})

test_that("gehan-wilcox handles all-dead and no-ties data correctly", {
    d <- make_surv_data(80L, seed = 3L)
    d$dead <- rep(1L, 80L)  # all events
    res <- EDI:::fast_gehan_wilcox_stats_cpp(d$w, d$y, d$dead)
    expect_true(is.finite(res$score))
    expect_true(is.finite(res$beta_hat))
    expect_true(is.finite(res$se_beta_hat))
})

test_that("gehan-wilcox output matches across multiple seeds", {
    for (seed in c(10L, 20L, 30L, 40L)) {
        d   <- make_surv_data(200L, seed = seed)
        res <- EDI:::fast_gehan_wilcox_stats_cpp(d$w, d$y, d$dead)
        ref <- canonical_gehan_wilcox(d)
        expect_equal(res$score,       ref$score,       tolerance = 1e-8, label = paste("score seed", seed))
        expect_equal(res$var_score,   ref$var_score,   tolerance = 1e-8, label = paste("var_score seed", seed))
        expect_equal(res$beta_hat,    ref$beta_hat,    tolerance = 1e-8, label = paste("beta_hat seed", seed))
        expect_equal(res$se_beta_hat, ref$se_beta_hat, tolerance = 1e-8, label = paste("se_beta_hat seed", seed))
    }
})

test_that("InferenceSurvivalGehanWilcox estimate/pval match canonical survival:: calls end-to-end", {
    set.seed(7L)
    n <- 60L
    seq_des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "survival", verbose = FALSE)
    for (i in seq_len(n)) {
        seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
    }
    w    <- as.integer((seq_des$get_w() + 1L) / 2L) # get_w() is {-1,+1}; Inference classes use {0,1}
    y    <- rexp(n, rate = exp(0.3 * w))
    dead <- as.integer(rbinom(n, 1, 0.8))
    seq_des$add_all_subject_responses(y, dead)

    inf <- InferenceSurvivalGehanWilcox$new(seq_des)
    est <- inf$compute_estimate()
    pv  <- inf$compute_asymp_two_sided_pval()

    ref <- canonical_gehan_wilcox(list(w = w, y = y, dead = dead))
    ref_chisq <- ref$score^2 / ref$var_score
    ref_pval  <- stats::pchisq(ref_chisq, df = 1, lower.tail = FALSE)

    expect_equal(est, ref$beta_hat, tolerance = 1e-8)
    expect_equal(pv,  ref_pval,     tolerance = 1e-8)
})
