options(warn = -1)
.libPaths(c(file.path(Sys.getenv("HOME"), "R", paste0(R.version$platform, "-library"), paste(R.version$major, sub("\\..*$", "", R.version$minor), sep = ".")), .libPaths()))
compiler::enableJIT(0) # Disable JIT to prevent R6 on-the-fly compilation overhead
suppressPackageStartupMessages(library(microbenchmark))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(VGAM))
suppressPackageStartupMessages(library(DescTools))
suppressPackageStartupMessages(library(ordinal))
suppressPackageStartupMessages(library(pscl))
suppressPackageStartupMessages(library(betareg))
suppressPackageStartupMessages(library(quantreg))
library(EDI)

set.seed(42)

STYLE_BLOCK = c(
  "<style>",
  "    body, .markdown-body, .container {",
  "        max-width: 1200px !important;",
  "        width: 100% !important;",
  "        margin: 0 auto !important;",
  "    }",
  "</style>"
)

# Global Config for Wald Tests (Full Inference)
N_WALD = 200 
B_TIME = 10 
TARGET_BATCH_MS = 200
MIN_RESOLVED_BATCH_MS = 10
MAX_INNER_REPS = 100000L
FAST_PATH_MICROBENCH_REPS = 2000L
FAST_PATH_THRESHOLD_MS = 0.01

# --- Data Generation Helper ---
generate_data = function(n = 200, p = 5, family = "logistic") {
    X = matrix(rnorm(n * p), n, p); X[, 1] = 1 
    beta = rnorm(p) * 0.5
    n_treat = floor(n / 2)
    w = sample(c(rep(1, n_treat), rep(0, n - n_treat)))
    eta = X %*% beta + 0.5 * w
    res = list(X = X, w = w)
    
    if (family %in% c("logistic", "log-binomial", "identity-binomial", "glmm_logistic", "robust", "probit")) {
        prob = if (family == "log-binomial") pmin(0.5, exp(eta - 2)) else plogis(eta)
        res$y = rbinom(n, 1, prob)
    } else if (family %in% c("poisson", "glmm_poisson", "zap", "quasipoisson", "hurdle_poisson_glmm", "modified_poisson")) {
        res$y = rpois(n, exp(eta))
    } else if (family %in% c("negbin", "zinb", "truncated_negbin", "hurdle_negbin")) {
        res$y = rnbinom(n, size = 2, mu = exp(eta))
    } else if (family %in% c("beta", "zoib", "fractional")) {
        mu = plogis(eta); phi = 10; res$y = pmax(pmin(rbeta(n, mu * phi, (1 - mu) * phi), 1 - 1e-6), 1e-6)
    } else if (family %in% c("weibull", "cox", "strat_cox", "dep_cens")) {
        res$y = rexp(n, exp(eta)); res$dead = rbinom(n, 1, 0.8)
    } else if (family %in% c("ordinal", "adj_cat", "continuation", "stereotype", "glmm_ordinal")) {
        p1 = 1 / (1 + exp(eta - 1))
        p2 = 1 / (1 + exp(eta + 1)) - p1
        p3 = 1 - p1 - p2
        probs = cbind(p1, p2, p3)
        probs = t(apply(probs, 1, function(x) pmax(x, 1e-6)))
        probs = probs / rowSums(probs)
        res$y = apply(probs, 1, function(p) sample(1:3, 1, prob = p))
    } else {
        res$y = as.numeric(eta + rnorm(n, 0, 0.5))
    }
    res
}

make_true_stratified_survival_data = function(d) {
    n = nrow(d$X)
    # Force low-cardinality covariates so the stratified Cox row benchmarks a
    # genuinely stratified fit rather than the unstratified fallback.
    strata_grid = as.matrix(expand.grid(x1 = 0:1, x2 = 0:2))
    idx = sample(rep(seq_len(nrow(strata_grid)), length.out = n))
    d$X[, 2] = strata_grid[idx, 1]
    d$X[, 3] = strata_grid[idx, 2]
    beta_cov = c(0.45, -0.35, 0.20, -0.15)
    eta = 0.5 * d$w + drop(d$X[, 2:5, drop = FALSE] %*% beta_cov)
    d$y = rexp(n, exp(eta))
    d$dead = rbinom(n, 1, 0.8)
    d
}

collect_timing_ms = function(expr, times = B_TIME, env = parent.frame(), target_batch_ms = TARGET_BATCH_MS, max_inner_reps = MAX_INNER_REPS, fast_path_microbenchmark_reps = FAST_PATH_MICROBENCH_REPS) {
    gctorture(FALSE)
    gc(verbose = FALSE)
    
    micro_time_ms = function() {
        gctorture(FALSE)
        gc(verbose = FALSE)
        bm_res = tryCatch({
            bench::mark(
                total = {
                    eval(expr, envir = env)
                },
                iterations = fast_path_microbenchmark_reps,
                check = FALSE,
                filter_gc = TRUE,
                memory = FALSE
            )
        }, error = function(e) NULL)
        if (is.null(bm_res)) {
            return(list(median_ms = NA_real_, samples_ms = numeric(0), method = "bench::mark"))
        }
        t_total = bm_res$time[[1]] * 1000
        median_total = as.numeric(bm_res$median[1]) * 1000
        
        list(
            median_ms = median_total,
            samples_ms = t_total,
            method = "bench::mark"
        )
    }
    inner_reps = 1L
    batch_ms = 0
    while (inner_reps < max_inner_reps) {
        batch_ms = system.time(for (j in seq_len(inner_reps)) {
            eval(expr, envir = env)
        })[["elapsed"]] * 1000
        if (is.finite(batch_ms) && batch_ms >= min(target_batch_ms, MIN_RESOLVED_BATCH_MS)) break
        inner_reps = min(max_inner_reps, inner_reps * 2L)
    }
    if (!is.finite(batch_ms) || batch_ms < MIN_RESOLVED_BATCH_MS) {
        return(micro_time_ms())
    }
    if (is.finite(batch_ms) && batch_ms > 0) {
        inner_reps = min(max_inner_reps, max(1L, ceiling(inner_reps * target_batch_ms / batch_ms)))
    }
    vals_ms = numeric(times)
    for (i in seq_len(times)) {
        gctorture(FALSE)
        gc(verbose = FALSE)
        t_total = system.time(for (j in seq_len(inner_reps)) {
            eval(expr, envir = env)
        })[["elapsed"]] * 1000
        
        vals_ms[i] = t_total / inner_reps
    }
    vals_ms = vals_ms[is.finite(vals_ms) & vals_ms > 0]
    median_ms = if (length(vals_ms) == 0L) NA_real_ else median(vals_ms)
    if (!is.finite(median_ms) || median_ms < FAST_PATH_THRESHOLD_MS) {
        return(micro_time_ms())
    }
    list(
        median_ms = median_ms,
        samples_ms = vals_ms,
        method = "system.time"
    )
}

timing_ttest_pval = function(x, y) {
    x = x[is.finite(x)]
    y = y[is.finite(y)]
    if (length(x) < 2L || length(y) < 2L) return(NA_real_)
    tryCatch(stats::t.test(x, y, var.equal = FALSE)$p.value, error = function(e) NA_real_)
}

make_edi_wald_bm = function(cls_name, d) {
    X_cov = d$X[, -1, drop = FALSE]
    colnames(X_cov) = paste0("x", seq_len(ncol(X_cov)))
    X_bm  = cbind(`(Intercept)` = 1, treatment = d$w, X_cov)
    X_ord = cbind(treatment = d$w, X_cov)
    y_bm  = as.numeric(d$y)
    w_bm  = as.integer(d$w)
    dead_bm = if (!is.null(d$dead)) as.integer(d$dead) else NULL

    e = new.env(parent = globalenv())
    e$d      = d
    e$X_bm   = X_bm
    e$X_ord  = X_ord
    e$y_bm   = y_bm
    e$w_bm   = w_bm
    e$dead_bm = dead_bm
    n = length(d$y)
    resp_type = if (!is.null(dead_bm)) "survival" else if (all(y_bm %in% c(0, 1))) "incidence" else "continuous"
    if (grepl("Count|Poisson|NegBin|ZINB|ZAP|Hurdle", cls_name)) resp_type = "count"
    if (grepl("Prop|Beta|ZOIB|Fractional", cls_name)) resp_type = "proportion"
    if (grepl("Ordinal|AdjCat|ContRatio|Stereotype|Ridit|Sign", cls_name)) resp_type = "ordinal"
    
    des = DesignFixediBCRD$new(n = n, response_type = resp_type)
    des$add_all_subjects_to_experiment(as.data.frame(X_cov))
    des$overwrite_all_subject_assignments(d$w)
    des$add_all_subject_responses(y_bm, deads = dead_bm)
    e$des = des

    e$compute_md_gradient = function(X_fit, theta, n_alpha, j_treat, base_step = 1e-6){
        n_params = length(theta)
        grad = numeric(n_params)
        for (j in seq_len(n_params)){
            step = max(base_step, base_step * (1 + abs(theta[j])))
            theta_plus = theta
            theta_minus = theta
            theta_plus[j] = theta[j] + step
            theta_minus[j] = theta[j] - step
            md_plus = gcomp_ordinal_proportional_odds_post_fit_cpp(X_fit, theta_plus[(n_alpha + 1):n_params], theta_plus[1:n_alpha], j_treat)$md
            md_minus = gcomp_ordinal_proportional_odds_post_fit_cpp(X_fit, theta_minus[(n_alpha + 1):n_params], theta_minus[1:n_alpha], j_treat)$md
            grad[j] = (md_plus - md_minus) / (2 * step)
        }
        grad
    }

    expr = switch(cls_name,
        # --- Continuous/OLS/Lin ---
        InferenceContinOLS = quote({
            res = fast_ols_with_var_cpp(X_bm, y_bm, j = 2L)
            est = res$b[2]
            se = sqrt(res$ssq_b_j)
            crit = stats::qt(0.975, df = length(y_bm) - ncol(X_bm))
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pt(-abs(t_stat), df = length(y_bm) - ncol(X_bm))
        }),

        InferenceContinLin = quote({
            Xc = scale(X_bm[, -1, drop = FALSE], scale = FALSE)
            X_int = Xc * w_bm
            X_lin = cbind(1, w_bm, Xc, X_int)
            colnames(X_lin)[1:2] = c("(Intercept)", "treatment")
            res_fit = fast_ols_cpp(X_lin, y_bm)
            post_fit = ols_hc2_post_fit_cpp(X_lin, y_bm, as.numeric(res_fit$b), 2L)
            est = res_fit$b[2]
            se = sqrt(post_fit$ssq_hat)
            crit = stats::qt(0.975, df = length(y_bm) - ncol(X_lin))
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pt(-abs(t_stat), df = length(y_bm) - ncol(X_lin))
        }),

        InferenceContinRobustRegr = quote({
            res = fast_robust_regression_cpp(X_bm, y_bm, j = 2L)
            est = res$b[2]
            se = sqrt(res$ssq_b_j)
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        InferenceContinQuantileRegr = quote({
            res = quantreg::rq(y_bm ~ X_bm[, -1])
            summary(res, se = "nid")$coefficients[2, 4]
        }),

        # --- Incidence/Logit/GLM ---
        InferenceIncidLogRegr = quote({
            res = fast_logistic_regression_with_var_cpp(X_bm, y_bm, j = 2L)
            est = res$b[2]
            se = sqrt(res$ssq_b_j)
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        InferenceIncidLogBinomial = quote({
            res = fast_log_binomial_regression_with_var_cpp(X_bm, y_bm, j = 2L)
            est = res$b[2]
            se = sqrt(res$ssq_b_j)
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        InferenceIncidProbitRegr = quote({
            res = fast_probit_regression_with_var_cpp(X_bm, y_bm, j = 2L)
            est = res$b[2]
            se = sqrt(res$ssq_b_j)
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        InferenceIncidRiskDiff = quote({
            res = fast_ols_with_var_cpp(X_bm, y_bm, j = 2L)
            est = res$b[2]
            se = sqrt(res$ssq_b_j)
            crit = stats::qt(0.975, df = length(y_bm) - ncol(X_bm))
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pt(-abs(t_stat), df = length(y_bm) - ncol(X_bm))
        }),

        InferenceIncidExactFisher = quote({
            tab = table(factor(w_bm, levels = c(1, 0)), factor(y_bm, levels = c(1, 0)))
            stats::fisher.test(tab)$p.value
        }),


        InferenceIncidNewcombeRiskDiff = quote({
            x_t = sum(y_bm[w_bm == 1])
            n_t = sum(w_bm == 1)
            x_c = sum(y_bm[w_bm == 0])
            n_c = sum(w_bm == 0)
            est = x_t / n_t - x_c / n_c
            p_fn = function(a) {
                ci = newcombe_independent_ci_cpp(x_t, n_t, x_c, n_c, a)
                if (0 < est) ci[1] else ci[2]
            }
            res = tryCatch(stats::uniroot(p_fn, interval = c(1e-10, 1 - 1e-10))$root, error = function(e) NA_real_)
            if (!is.finite(res)) 1.0 else res
        }),

        InferenceIncidMiettinenNurminenRiskDiff = quote({
            x_t = sum(y_bm[w_bm == 1])
            n_t = sum(w_bm == 1)
            x_c = sum(y_bm[w_bm == 0])
            n_c = sum(w_bm == 0)
            p_t = x_t / n_t
            p_c = x_c / n_c
            mn_pvalue_cpp(x_t, n_t, x_c, n_c, 0, p_t, p_c)
        }),

        # --- Count ---
        InferenceCountPoisson = quote({
            res = fast_poisson_regression_with_var_cpp(X_bm, y_bm, j = 2L)
            est = res$b[2]
            se = sqrt(res$ssq_b_j)
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        InferenceCountNegBin = quote({
            res = fast_neg_bin_with_var_cpp(X_bm, as.integer(y_bm))
            est = res$b[2]
            hess = res$hess_fisher_info_matrix
            vcov = solve(hess)
            se = sqrt(vcov[2, 2])
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        InferenceCountHurdlePoisson = quote({
            res = fast_zero_augmented_poisson_cpp(X_bm, y_bm, X_bm, is_hurdle = TRUE, estimate_only = FALSE)
            est = res$coefficients$cond[2]
            se = sqrt(res$vcov[2, 2])
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        InferenceCountZeroInflatedPoisson = quote({
            res = fast_zero_augmented_poisson_cpp(X_bm, y_bm, X_bm, is_hurdle = FALSE, estimate_only = FALSE)
            est = res$coefficients$cond[2]
            se = sqrt(res$vcov[2, 2])
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        InferenceCountZeroInflatedNegBin = quote({
            res = fast_zinb_cpp(X_bm, X_bm, y_bm, estimate_only = FALSE)
            est = res$params[2]
            se = sqrt(res$vcov[2, 2])
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        InferenceCountHurdleNegBin = quote({
            res = fast_hurdle_negbin_with_var_cpp(X_bm, as.integer(y_bm), X_bm, j = 2L)
            est = res$b[2]
            se = sqrt(res$ssq_b_j)
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        InferenceCountQuasiPoisson = quote({
            res = fast_quasipoisson_regression_with_var_cpp(X_bm, y_bm, j = 2L)
            est = res$b[2]
            se = sqrt(res$ssq_b_j)
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        InferenceCountRobustPoisson = quote({
            res = fast_poisson_regression_cpp(X_bm, y_bm, estimate_only = FALSE)
            est = res$b[2]
            bread = solve(res$XtWX)
            resid = y_bm - as.numeric(res$mu)
            meat = crossprod(X_bm, X_bm * (resid^2))
            vcov_robust = bread %*% meat %*% bread
            se = sqrt(vcov_robust[2, 2])
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        # --- Proportion ---
        InferencePropBetaRegr = quote({
            res = fast_beta_regression_with_var_cpp(X_bm, y_bm)
            est = res$coefficients[2]
            se = sqrt(res$vcov[2, 2])
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        # --- Ordinal ---
        InferenceOrdinalPropOddsRegr = quote({
            res = fast_ordinal_regression_with_var_cpp(X_ord, y_bm)
            est = res$b[1]
            se = sqrt(res$ssq_b_j)
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        InferenceOrdinalAdjCatLogitRegr = quote({
            res = fast_adjacent_category_logit_with_var_cpp(X_ord, y_bm)
            est = res$b[1]
            se = sqrt(res$ssq_b_j)
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        InferenceOrdinalContRatioRegr = quote({
            res = fast_continuation_ratio_regression_with_var_cpp(X_ord, y_bm)
            est = res$b[1]
            se = sqrt(res$ssq_b_j)
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        InferenceOrdinalOrderedProbitRegr = quote({
            res = fast_ordinal_probit_regression_with_var_cpp(X_ord, y_bm)
            est = res$b[1]
            se = sqrt(res$ssq_b_j)
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        InferenceOrdinalCloglogRegr = quote({
            res = fast_ordinal_cloglog_regression_with_var_cpp(X_ord, y_bm)
            est = res$b[1]
            se = sqrt(res$ssq_b_j)
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        InferenceOrdinalCauchitRegr = quote({
            res = fast_ordinal_cauchit_regression_with_var_cpp(X_ord, y_bm)
            est = res$b[1]
            se = sqrt(res$ssq_b_j)
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        InferenceOrdinalRidit = quote({
            res = EDI:::fast_ridit_analysis_cpp(w_bm, as.integer(y_bm), reference = "control")
            est = res$estimate
            se = res$se
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        # --- Survival ---
        InferenceSurvivalCoxPHRegr = quote({
            res = fast_coxph_regression_cpp(X_ord, y_bm, dead_bm, estimate_only = FALSE)
            est = res$coefficients[1]
            se = sqrt(res$vcov[1, 1])
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        InferenceSurvivalStratCoxPHRegr = quote({
            strat_inputs = build_strat_cox_canonical_inputs(d)
            if (!is.null(strat_inputs$strata)) {
                cache = build_stratified_cox_data_cache_cpp(strat_inputs$X, strat_inputs$y, strat_inputs$dead, strat_inputs$strata)
            } else {
                cache = build_cox_data_cache_cpp(strat_inputs$X, strat_inputs$y, strat_inputs$dead)
            }
            res = fast_coxph_regression_prebuilt_cpp(cache, estimate_only = FALSE)
            est = res$coefficients[1]
            se = sqrt(res$vcov[1, 1])
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        InferenceSurvivalWeibullRegr = quote({
            res = fast_weibull_regression_cpp(X_ord, y_bm, dead_bm, estimate_only = FALSE)
            est = res$params[2]
            se = sqrt(res$vcov[2, 2])
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        InferenceSurvivalLogRank = quote({
            res = EDI:::fast_logrank_stats_cpp(w_bm, y_bm, as.integer(dead_bm))
            chisq_stat = res$score ^ 2 / res$var_score
            stats::pchisq(chisq_stat, df = 1, lower.tail = FALSE)
        }),

        InferenceSurvivalGehanWilcox = quote({
            survival::survdiff(survival::Surv(y_bm, dead_bm) ~ w_bm, rho = 1)$pvalue
        }),

        InferenceSurvivalKMDiff = quote({
            fit = survival::survfit(survival::Surv(y_bm, dead_bm) ~ w_bm, conf.int = 0.95)
            q = stats::quantile(fit, 0.5)
            strata_names = rownames(q$lower)
            idx_T = grep("1", strata_names)
            idx_C = grep("0", strata_names)
            lo_T = q$lower[idx_T, 1]
            hi_T = q$upper[idx_T, 1]
            lo_C = q$lower[idx_C, 1]
            hi_C = q$upper[idx_C, 1]
            z = stats::qnorm(0.975)
            se = sqrt(((hi_T - lo_T) / (2 * z))^2 + ((hi_C - lo_C) / (2 * z))^2)
            est = q$quantile[idx_T, 1] - q$quantile[idx_C, 1]
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),

        # --- Simple/Other ---
        InferenceAllSimpleMeanDiffPooledVar = quote({
            y_t = y_bm[w_bm == 1]
            y_c = y_bm[w_bm == 0]
            n_t = length(y_t)
            n_c = length(y_c)
            s2_t = stats::var(y_t)
            s2_c = stats::var(y_c)
            df = n_t + n_c - 2L
            s2_pooled = ((n_t - 1L) * s2_t + (n_c - 1L) * s2_c) / df
            se = sqrt(s2_pooled * (1 / n_t + 1 / n_c))
            est = mean(y_t) - mean(y_c)
            crit = stats::qt(0.975, df = df)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pt(-abs(t_stat), df = df)
        }),

        InferenceAllSimpleWilcox = quote({
            EDI:::wilcox_hl_point_estimate_cpp(w_bm, y_bm)
        }),

        # --- GComp ---
        InferencePropGCompMeanDiff = quote({
            fit = fast_logistic_regression_cpp(X_bm, y_bm, estimate_only = TRUE)
            coef_hat = as.numeric(fit$b)
            mu_hat = stats::plogis(as.numeric(X_bm %*% coef_hat))
            res = gcomp_logistic_post_fit_cpp(X_bm, y_bm, coef_hat, mu_hat, 2L)
            est = res$rd
            se = res$se_rd
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),
 
        InferenceIncidGCompRiskDiff = quote({
            fit = fast_logistic_regression_cpp(X_bm, y_bm, estimate_only = TRUE)
            coef_hat = as.numeric(fit$b)
            mu_hat = stats::plogis(as.numeric(X_bm %*% coef_hat))
            res = gcomp_logistic_post_fit_cpp(X_bm, y_bm, coef_hat, mu_hat, 2L)
            est = res$rd
            se = res$se_rd
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),
 
        InferenceIncidGCompRiskRatio = quote({
            fit = fast_logistic_regression_cpp(X_bm, y_bm, estimate_only = TRUE)
            coef_hat = as.numeric(fit$b)
            mu_hat = stats::plogis(as.numeric(X_bm %*% coef_hat))
            res = gcomp_logistic_post_fit_cpp(X_bm, y_bm, coef_hat, mu_hat, 2L)
            t_stat = res$log_rr / res$se_log_rr
            2 * stats::pnorm(-abs(t_stat))
        }),
 
        InferenceOrdinalGCompMeanDiff = quote({
            fit = fast_ordinal_regression_with_var_cpp(X_ord, y_bm)
            coef_hat = as.numeric(fit$b)
            alpha_hat = as.numeric(fit$alpha)
            res = gcomp_ordinal_proportional_odds_post_fit_cpp(X_ord, coef_hat, alpha_hat, 1L)
            est = res$md
            theta = c(alpha_hat, coef_hat)
            grad = compute_md_gradient(X_ord, theta, length(alpha_hat), 1L)
            var_md = as.numeric(crossprod(grad, as.matrix(fit$vcov) %*% grad))
            se = sqrt(var_md)
            crit = stats::qnorm(0.975)
            ci = c(est - crit * se, est + crit * se)
            t_stat = est / se
            2 * stats::pnorm(-abs(t_stat))
        }),
 
        InferenceOrdinalJonckheereTerpstraTest = quote({
            exact_jonckheere_terpstra_pval_cpp(as.integer(y_bm), w_bm)$p_exact
        }),

        NULL
    )
    list(env = e, expr = expr)
}

build_strat_cox_canonical_inputs = function(d) {
    X_cov = d$X[, -1, drop = FALSE]
    strata_info = EDI:::compute_survival_strata_ids_cpp(as.matrix(X_cov))
    X_linear = matrix(numeric(0), nrow = nrow(X_cov), ncol = 0)
    if (!is.null(X_cov) && ncol(X_cov) > 0) {
        keep_cols = setdiff(seq_len(ncol(X_cov)), as.integer(strata_info$selected_cols))
        if (length(keep_cols) > 0) {
            full_design = cbind(w = d$w, X_cov[, keep_cols, drop = FALSE])
            reduced = EDI:::drop_linearly_dependent_cols(full_design)$M
            if ("w" %in% colnames(reduced)) {
                X_linear = reduced[, colnames(reduced) != "w", drop = FALSE]
            }
        }
    }
    strata_id = as.integer(strata_info$strata_id)
    informative_rows = integer(0)
    for (s in unique(strata_id)) {
        i_s = which(strata_id == s)
        if (length(i_s) < 2L) next
        if (length(unique(d$w[i_s])) < 2L) next
        if (!any(d$dead[i_s] == 1, na.rm = TRUE)) next
        informative_rows = c(informative_rows, i_s)
    }
    informative_rows = sort(unique(informative_rows))
    if (length(informative_rows) >= 4L) {
        x = if (ncol(X_linear) > 0) cbind(w = d$w[informative_rows], X_linear[informative_rows, , drop = FALSE]) else matrix(d$w[informative_rows], ncol = 1L, dimnames = list(NULL, "w"))
        return(list(
            X = as.matrix(x),
            y = as.numeric(d$y[informative_rows]),
            dead = as.numeric(d$dead[informative_rows]),
            strata = as.integer(strata_id[informative_rows])
        ))
    }
    x = if (ncol(X_linear) > 0) cbind(w = d$w, X_linear) else matrix(d$w, ncol = 1L, dimnames = list(NULL, "w"))
    list(
        X = as.matrix(x),
        y = as.numeric(d$y),
        dead = as.numeric(d$dead),
        strata = NULL
    )
}

compute_binary_gcomp_effect = function(beta, X_base, effect = c("RD", "RR")) {
    effect = match.arg(effect)
    X1 = X_base
    X0 = X_base
    X1[, "treatment"] = 1
    X0[, "treatment"] = 0
    risk1_i = stats::plogis(drop(X1 %*% beta))
    risk0_i = stats::plogis(drop(X0 %*% beta))
    risk1 = mean(risk1_i)
    risk0 = mean(risk0_i)
    if (effect == "RD") {
        grad1 = drop(crossprod(X1, risk1_i * (1 - risk1_i))) / nrow(X1)
        grad0 = drop(crossprod(X0, risk0_i * (1 - risk0_i))) / nrow(X0)
        list(est = risk1 - risk0, grad = grad1 - grad0)
    } else {
        grad1 = drop(crossprod(X1, risk1_i * (1 - risk1_i))) / nrow(X1)
        grad0 = drop(crossprod(X0, risk0_i * (1 - risk0_i))) / nrow(X0)
        log_rr = log(risk1) - log(risk0)
        list(est = exp(log_rr), log_est = log_rr, grad = grad1 / risk1 - grad0 / risk0)
    }
}

compute_ordinal_gcomp_md = function(theta, X_base, n_alpha) {
    alpha = theta[seq_len(n_alpha)]
    beta = theta[-seq_len(n_alpha)]
    X1 = X_base
    X0 = X_base
    X1[, "treatment"] = 1
    X0[, "treatment"] = 0
    mean_score = function(Xmat) {
        eta = drop(Xmat %*% beta)
        cumprob = vapply(alpha, function(a) stats::plogis(a - eta), numeric(nrow(Xmat)))
        if (is.vector(cumprob)) cumprob = matrix(cumprob, ncol = 1)
        probs = matrix(0, nrow = nrow(Xmat), ncol = n_alpha + 1L)
        probs[, 1L] = cumprob[, 1L]
        if (n_alpha > 1L) {
            for (j in 2:n_alpha) probs[, j] = cumprob[, j] - cumprob[, j - 1L]
        }
        probs[, n_alpha + 1L] = 1 - cumprob[, n_alpha]
        drop(probs %*% seq_len(n_alpha + 1L))
    }
    mean(mean_score(X1)) - mean(mean_score(X0))
}

finite_diff_grad = function(fn, theta, base_step = 1e-6) {
    grad = numeric(length(theta))
    for (j in seq_along(theta)) {
        step = max(base_step, base_step * (1 + abs(theta[j])))
        theta_plus = theta
        theta_minus = theta
        theta_plus[j] = theta[j] + step
        theta_minus[j] = theta[j] - step
        grad[j] = (fn(theta_plus) - fn(theta_minus)) / (2 * step)
    }
    grad
}

# --- Mapping Table (Full Wald Inference) ---
bench_specs = list(
    list(cls = "InferenceIncidLogRegr", pkg = "stats", func = "glm.fit+Wald", expr = quote({
        res = glm.fit(x = X_can, y = df$y, family = binomial())
        p = ncol(X_can); R = res$qr$qr[1:p, 1:p, drop=FALSE]; R[lower.tri(R)] = 0
        v = chol2inv(R); 2 * pnorm(-abs(res$coefficients[2] / sqrt(v[2,2])))
    })),
    list(cls = "InferenceContinOLS", pkg = "stats", func = "lm.fit+Wald", expr = quote({
        res = lm.fit(x = X_can, y = df$y)
        p = ncol(X_can); n = nrow(X_can); rss = sum(res$residuals^2); df_res = n - p; sig2 = rss / df_res
        v = chol2inv(res$qr$qr[1:p, 1:p, drop=FALSE]) * sig2; 2 * pt(-abs(res$coefficients[2] / sqrt(v[2,2])), df_res)
    })),
    list(cls = "InferenceCountPoisson", pkg = "stats", func = "glm.fit+Wald", expr = quote({
        res = glm.fit(x = X_can, y = df$y, family = poisson())
        p = ncol(X_can); R = res$qr$qr[1:p, 1:p, drop=FALSE]; R[lower.tri(R)] = 0
        v = chol2inv(R); 2 * pnorm(-abs(res$coefficients[2] / sqrt(v[2,2])))
    })),
    list(cls = "InferenceSurvivalCoxPHRegr", pkg = "survival", func = "coxph.fit(breslow)+Wald", expr = quote({
        x_vars = c("treatment", grep("^x", names(df), value=T))
        res = survival::coxph.fit(x = as.matrix(df[, x_vars, drop=F]), y = survival::Surv(df$y, df$dead), strata=NULL, offset=NULL, init=NULL, control=survival::coxph.control(), weights=NULL, method="breslow", rownames=as.character(1:nrow(df)))
        2 * pnorm(-abs(res$coefficients[1] / sqrt(res$var[1,1])))
    })),
    list(cls = "InferenceCountNegBin", pkg = "MASS", func = "glm.nb+summary", expr = quote({
        res = MASS::glm.nb(y ~ treatment + x1 + x2 + x3 + x4, data = df)
        summary(res)$coefficients[2, 4]
    })),
    list(cls = "InferencePropBetaRegr", pkg = "betareg", func = "betareg+summary", expr = quote({
        res = betareg::betareg(y ~ treatment + x1 + x2 + x3 + x4, data = df)
        summary(res)$coefficients$mean[2, 4]
    })),
    list(cls = "InferenceOrdinalPropOddsRegr", pkg = "ordinal", func = "clm+summary", expr = quote({
        res = ordinal::clm(factor(y, ordered=T) ~ treatment + x1 + x2 + x3 + x4, data = df)
        summary(res)$coefficients["treatment", 4]
    })),
    list(cls = "InferenceCountHurdlePoisson", pkg = "pscl", func = "hurdle+summary", expr = quote({
        res = pscl::hurdle(y ~ treatment + x1 + x2 + x3 + x4, data = df)
        summary(res)$coefficients$count[2, 4]
    })),
    list(cls = "InferenceCountZeroInflatedPoisson", pkg = "pscl", func = "zeroinfl+summary", expr = quote({
        res = pscl::zeroinfl(y ~ treatment + x1 + x2 + x3 + x4, data = df)
        summary(res)$coefficients$count[2, 4]
    })),
    list(cls = "InferenceCountZeroInflatedNegBin", pkg = "pscl", func = "zeroinfl(nb)+summary", expr = quote({
        res = pscl::zeroinfl(y ~ treatment + x1 + x2 + x3 + x4, data = df, dist="negbin")
        summary(res)$coefficients$count[2, 4]
    })),
    list(cls = "InferenceCountHurdleNegBin", pkg = "pscl", func = "hurdle(nb)+summary", expr = quote({
        res = pscl::hurdle(y ~ treatment + x1 + x2 + x3 + x4, data = df, dist="negbin")
        summary(res)$coefficients$count[2, 4]
    })),
    list(cls = "InferenceCountQuasiPoisson", pkg = "stats", func = "glm.fit+Wald(quasi)", expr = quote({
        res = glm.fit(x = X_can, y = df$y, family = quasipoisson())
        p = ncol(X_can); n = nrow(X_can); R = res$qr$qr[1:p, 1:p, drop=FALSE]; R[lower.tri(R)] = 0
        disp = sum(res$weights * res$residuals^2) / res$df.residual
        v = chol2inv(R) * disp; 2 * pnorm(-abs(res$coefficients[2] / sqrt(v[2,2])))
    })),
    list(cls = "InferenceSurvivalWeibullRegr", pkg = "survival", func = "survreg+summary", expr = quote({
        res = survival::survreg(survival::Surv(y, dead) ~ treatment + x1 + x2 + x3 + x4, data = df, dist="weibull")
        summary(res)$table[2, 4]
    })),
    list(cls = "InferenceContinRobustRegr", pkg = "MASS", func = "rlm+summary", expr = quote({
        res = MASS::rlm(x = X_can, y = df$y)
        summary(res)$coefficients[2, 3] 
    })),
    list(cls = "InferenceContinQuantileRegr", pkg = "quantreg", func = "rq+summary", expr = quote({
        res = quantreg::rq(df$y ~ X_can[,-1])
        summary(res, se="nid")$coefficients[2, 4]
    })),
    list(cls = "InferenceIncidLogBinomial", pkg = "stats", func = "glm.fit+Wald(log)", expr = quote({
        res = glm.fit(x = X_can, y = df$y, family=binomial(link="log"), start=c(-2, rep(0, ncol(X_can)-1)))
        p = ncol(X_can); R = res$qr$qr[1:p, 1:p, drop=FALSE]; R[lower.tri(R)] = 0
        v = chol2inv(R); 2 * pnorm(-abs(res$coefficients[2] / sqrt(v[2,2])))
    })),
    list(cls = "InferenceIncidProbitRegr", pkg = "stats", func = "glm.fit(probit)+Wald", expr = quote({
        res = glm.fit(x = X_can, y = df$y, family = binomial(link = "probit"))
        p = ncol(X_can); R = res$qr$qr[1:p, 1:p, drop = FALSE]; R[lower.tri(R)] = 0
        v = chol2inv(R); 2 * pnorm(-abs(res$coefficients[2] / sqrt(v[2, 2])))
    })),
    list(cls = "InferenceOrdinalAdjCatLogitRegr", pkg = "VGAM", func = "vglm+summary", expr = quote({
        res = VGAM::vglm(factor(y, ordered=T) ~ treatment + x1 + x2 + x3 + x4, VGAM::acat(), data=df)
        summary(res)@coef3[1, 4]
    })),
    list(cls = "InferenceOrdinalContRatioRegr", pkg = "VGAM", func = "vglm+summary", expr = quote({
        res = VGAM::vglm(factor(y, ordered=T) ~ treatment + x1 + x2 + x3 + x4, VGAM::cratio(), data=df)
        summary(res)@coef3[1, 4]
    })),
    list(cls = "InferenceSurvivalLogRank", pkg = "survival", func = "survdiff", expr = quote(survival::survdiff(survival::Surv(y, dead) ~ treatment, data = df)$pvalue)),
    list(cls = "InferenceSurvivalGehanWilcox", pkg = "survival", func = "coxph(null)+KM weighted residual mean diff + survdiff(rho=1)", expr = quote({
        surv_obj = survival::Surv(df$y, df$dead)
        cox_null = survival::coxph(surv_obj ~ 1)
        M = as.numeric(stats::residuals(cox_null, type = "martingale"))
        km_all = survival::survfit(surv_obj ~ 1)
        idx = findInterval(df$y, km_all$time, left.open = TRUE)
        peto_weights = c(1.0, km_all$surv)[idx + 1L]
        M_w = peto_weights * M
        est = mean(M_w[df$treatment == 1]) - mean(M_w[df$treatment == 0])
        p = survival::survdiff(surv_obj ~ treatment, data = df, rho = 1)$pvalue
        if (!is.finite(est) || !is.finite(p)) stop("Non-finite canonical Gehan-Wilcoxon estimate or p-value.")
        p
    })),
    list(cls = "InferenceAllSimpleMeanDiffPooledVar", pkg = "stats", func = "t.test(pool)", expr = quote(t.test(df$y[df$treatment==1], df$y[df$treatment==0], var.equal = TRUE)$p.value)),
    list(cls = "InferenceAllSimpleWilcox", pkg = "stats", func = "wilcox.test", expr = quote(wilcox.test(df$y[df$treatment==1], df$y[df$treatment==0])$p.value)),
    list(cls = "InferenceIncidExactFisher", pkg = "stats", func = "fisher.test", expr = quote(fisher.test(table(df$treatment, df$y))$p.value)),
    list(cls = "InferenceOrdinalJonckheereTerpstraTest", pkg = "clinfun", func = "jonckheere", expr = quote(clinfun::jonckheere.test(df$y, df$treatment)$p.value)),
    list(cls = "InferenceIncidMiettinenNurminenRiskDiff", pkg = "DescTools", func = "BinomDiffCI(mn)", expr = quote(DescTools::BinomDiffCI(sum(df$y[df$treatment==1]), sum(df$treatment==1), sum(df$y[df$treatment==0]), sum(df$treatment==0), method="mn"))),
    list(cls = "InferenceSurvivalStratCoxPHRegr", pkg = "survival", func = "coxph.fit(strat,breslow)+Wald", expr = quote({
        strat_inputs = build_strat_cox_canonical_inputs(d)
        res = survival::coxph.fit(
            x = strat_inputs$X,
            y = survival::Surv(strat_inputs$y, strat_inputs$dead),
            strata = strat_inputs$strata,
            offset = NULL,
            init = NULL,
            control = survival::coxph.control(),
            weights = NULL,
            method = "breslow",
            rownames = as.character(seq_along(strat_inputs$y))
        )
        2 * pnorm(-abs(res$coefficients[1] / sqrt(res$var[1,1])))
    })),
    list(cls = "InferenceContinLin", pkg = "stats", func = "lm.fit(interact)+Wald", expr = quote({
        X_int = model.matrix(~ treatment * (x1 + x2 + x3 + x4), data = df)
        res = lm.fit(x = X_int, y = df$y)
        p = ncol(X_int); n = nrow(X_int); rss = sum(res$residuals^2); df_res = n - p; sig2 = rss / df_res
        v = chol2inv(res$qr$qr[1:p, 1:p, drop=FALSE]) * sig2; 2 * pt(-abs(res$coefficients[2] / sqrt(v[2,2])), df_res)
    })),
    list(cls = "InferenceIncidRiskDiff", pkg = "stats", func = "prop.test", expr = quote(prop.test(table(df$treatment, df$y))$p.value)),
    list(cls = "InferenceSurvivalKMDiff", pkg = "survival", func = "survfit(median)+CI", expr = quote({
        fit = survival::survfit(survival::Surv(df$y, df$dead) ~ df$treatment)
        stats::quantile(fit, 0.5)
    })),
    list(cls = "InferenceCountRobustPoisson", pkg = "sandwich", func = "glm+vcovHC", expr = quote({
        mod = glm(y ~ treatment + x1 + x2 + x3 + x4, family=poisson, data=df)
        v = sandwich::vcovHC(mod, type="HC0")
        2 * pnorm(-abs(coef(mod)[2] / sqrt(v[2,2])))
    })),
    list(cls = "InferenceOrdinalRidit", pkg = "stats", func = "mean(ridit)", expr = quote({
        y = df$y; tab = table(y); cum = cumsum(tab); prev = c(0, cum[-length(cum)])
        ridit_map = (prev + 0.5 * tab) / length(y); r = as.numeric(ridit_map[as.character(y)])
        mean(r[df$treatment==1]) - mean(r[df$treatment==0])
    })),
    list(cls = "InferenceIncidNewcombeRiskDiff", pkg = "DescTools", func = "BinomDiffCI(score)", expr = quote(DescTools::BinomDiffCI(sum(df$y[df$treatment==1]), sum(df$treatment==1), sum(df$y[df$treatment==0]), sum(df$treatment==0), method="score")))
)

# Paths with no canonical mapping (None)
no_can_specs = list(
    list(cls = "InferencePropGCompMeanDiff", pkg = "stats", func = "glm(quasi)+gcomp+Wald", expr = quote({
        mod = glm(y ~ treatment + x1 + x2 + x3 + x4, family = quasibinomial(), data = df)
        beta = stats::coef(mod)
        vc = stats::vcov(mod)
        eff = compute_binary_gcomp_effect(beta, X_can, effect = "RD")
        se = sqrt(drop(t(eff$grad) %*% vc %*% eff$grad))
        2 * stats::pnorm(-abs(eff$est / se))
    })),
    list(cls = "InferenceIncidGCompRiskRatio", pkg = "stats", func = "glm+gcomp(RR)+Wald", expr = quote({
        mod = glm(y ~ treatment + x1 + x2 + x3 + x4, family = binomial(), data = df)
        beta = stats::coef(mod)
        vc = stats::vcov(mod)
        eff = compute_binary_gcomp_effect(beta, X_can, effect = "RR")
        se = sqrt(drop(t(eff$grad) %*% vc %*% eff$grad))
        2 * stats::pnorm(-abs(eff$log_est / se))
    })),
    list(cls = "InferenceIncidGCompRiskDiff", pkg = "stats", func = "glm+gcomp(RD)+Wald", expr = quote({
        mod = glm(y ~ treatment + x1 + x2 + x3 + x4, family = binomial(), data = df)
        beta = stats::coef(mod)
        vc = stats::vcov(mod)
        eff = compute_binary_gcomp_effect(beta, X_can, effect = "RD")
        se = sqrt(drop(t(eff$grad) %*% vc %*% eff$grad))
        2 * stats::pnorm(-abs(eff$est / se))
    })),
    list(cls = "InferenceOrdinalGCompMeanDiff", pkg = "ordinal", func = "clm+gcomp+Wald", expr = quote({
        mod = ordinal::clm(factor(y, ordered = TRUE) ~ treatment + x1 + x2 + x3 + x4, data = df)
        theta = stats::coef(mod)
        vc = stats::vcov(mod)
        n_alpha = length(mod$alpha)
        fn = function(th) compute_ordinal_gcomp_md(th, X_can[, -1, drop = FALSE], n_alpha)
        est = fn(theta)
        grad = finite_diff_grad(fn, theta)
        se = sqrt(drop(t(grad) %*% vc %*% grad))
        2 * stats::pnorm(-abs(est / se))
    }))
)
bench_specs = c(bench_specs, no_can_specs)

# --- Benchmark Runner ---
results = list()

safe_instantiate = function(cls_name, des) {
	inf_cls = get(cls_name)
	init_args = list(des)
    f = formals(inf_cls$public_methods$initialize)
    if ("smart_cold_start_default" %in% names(f)) init_args$smart_cold_start_default = TRUE
    if ("optimization_alg" %in% names(f)) {
        if (grepl("LogRegr|Poisson|ModifiedPoisson|QuasiPoisson", cls_name) && !grepl("Hurdle|ZeroInflated", cls_name)) {
             init_args$optimization_alg = "irls"
        }
    }
    do.call(inf_cls$new, init_args)
}

run_one = function(spec) {
    cls_name = spec$cls
    cat(sprintf("Benchmarking Wald Test %s...\n", cls_name))
    
    resp_type = "continuous"
    family = "continuous"
    if (grepl("Ordinal|AdjCat|ContRatio|Stereotype|Ridit|Sign", cls_name)) { resp_type = "ordinal"; family = "ordinal" }
    else if (grepl("Incid|Binomial|Wald|CMH|Fisher|Robins|Newcombe|Nurminen", cls_name)) { resp_type = "incidence"; family = "logistic" }
    else if (grepl("Count|Poisson|NegBin|ZINB|ZAP|Hurdle", cls_name)) { resp_type = "count"; family = "poisson" }
    else if (grepl("Prop|Beta|ZOIB|Fractional", cls_name)) { resp_type = "proportion"; family = "beta" }
    else if (grepl("Survival|Cox|Weibull|KM|Rank|LogRank|Gehan|RMST|RMDiff|LWACox|Clayton", cls_name)) { resp_type = "survival"; family = "cox" }
    
    if (cls_name == "InferenceIncidLogBinomial") family = "log-binomial"
    if (cls_name == "InferenceIncidProbitRegr") family = "probit"
    if (cls_name %in% c("InferenceCountNegBin", "InferenceCountZeroInflatedNegBin", "InferenceCountHurdleNegBin")) family = "negbin"
    
    n = N_WALD 
    use_stable_count_draw = cls_name %in% c(
        "InferenceCountHurdlePoisson",
        "InferenceCountHurdleNegBin",
        "InferenceCountQuasiPoisson",
        "InferenceCountRobustPoisson"
    )
    if (use_stable_count_draw) {
        old_seed = if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) .Random.seed else NULL
        set.seed(42)
        d = generate_data(n = n, family = family)
        if (is.null(old_seed)) {
            if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
        } else {
            .Random.seed <<- old_seed
        }
    } else {
        d = generate_data(n = n, family = family)
    }
    if (identical(cls_name, "InferenceSurvivalStratCoxPHRegr")) {
        d = make_true_stratified_survival_data(d)
    }
    pval_delta = if (cls_name == "InferenceIncidGCompRiskRatio") 1 else 0
    
    # Timing EDI (Model Fit + Wald Test)
    timing_edi = tryCatch({
        bm = make_edi_wald_bm(cls_name, d)
        if (is.null(bm)) stop("no bare metal mapping for this class")
        eval(bm$expr, envir = bm$env)  # validation run
        collect_timing_ms(bm$expr, env = bm$env)
    }, error = function(e) { cat("  EDI Error:", e$message, "\n"); list(median_ms = NA_real_, samples_ms = numeric(0)) })
    
    timing_can = tryCatch({
        if (!is.null(spec$expr) && spec$pkg != "None") {
            library(spec$pkg, character.only = TRUE)
            X_cols = d$X[,-1,drop=F]; colnames(X_cols) = paste0("x", 1:ncol(X_cols))
            df = data.frame(y = d$y, treatment = d$w, g = factor(rep(1:10, length.out = n)), dead = if(exists("dead", d)) d$dead else 1)
            df = cbind(df, X_cols); X_can = cbind(`(Intercept)` = 1, treatment = d$w, as.matrix(X_cols))
            if (identical(cls_name, "InferenceSurvivalStratCoxPHRegr")) {
                strat_inputs = build_strat_cox_canonical_inputs(d)
                X_can_strat = strat_inputs$X
                y_can_strat = strat_inputs$y
                dead_can_strat = strat_inputs$dead
                strata_can = strat_inputs$strata
            }
            
            eval(spec$expr)
            collect_timing_ms(spec$expr)
        } else NA_real_
    }, error = function(e) { cat("  Canonical Error:", e$message, "\n"); list(median_ms = NA_real_, samples_ms = numeric(0)) })
    if (is.atomic(timing_can) && length(timing_can) == 1L && is.na(timing_can)) timing_can = list(median_ms = NA_real_, samples_ms = numeric(0))
    
    results[[length(results) + 1]] <<- data.table(
        Class = cls_name, Response = resp_type, EDI_Time_ms = timing_edi$median_ms,
        Canonical_Pkg = spec$pkg, Canonical_Func = spec$func, Canonical_Time_ms = timing_can$median_ms,
        Speedup = if (!is.na(timing_can$median_ms) && !is.na(timing_edi$median_ms) && timing_edi$median_ms > 0) timing_can$median_ms / timing_edi$median_ms else NA_real_,
        Timing_Pval = timing_ttest_pval(timing_edi$samples_ms, timing_can$samples_ms)
    )
}

unique_specs = bench_specs[!duplicated(sapply(bench_specs, `[[`, "cls"))]
for (i in seq_along(unique_specs)) {
    cat(sprintf("[%d/%d] ", i, length(unique_specs)))
    run_one(unique_specs[[i]])
}

# Finalize
dt = rbindlist(results)
saveRDS(dt, "benchmark/wald_test_results.rds")
if (identical(Sys.getenv("WALD_SAVE_ONLY"), "1")) {
    cat("Wald Test Benchmark timing results saved to benchmark/wald_test_results.rds.\n")
    quit(save = "no", status = 0)
}
response_levels = c("all", "continuous", "incidence", "count", "proportion", "survival", "ordinal")
dt[, Response := factor(Response, levels = response_levels)]
setorder(dt, Response, Class)
dt[, Speedup_Num := Speedup]
dt[, Speedup := ifelse(!is.na(Speedup), paste0(round(Speedup, 2), "x"), "NA")]
format_pval = function(x) {
  ifelse(is.na(x), "NA", vapply(x, function(v) sprintf("%.3g", v), character(1)))
}
format_pval_stars = function(x) {
  ifelse(
    is.na(x),
    "",
    ifelse(x < 0.001, "***", ifelse(x < 0.01, "**", ifelse(x < 0.05, "*", "")))
  )
}
row_bg_color = function(speedup, pval) {
  if (!is.finite(speedup) || is.na(pval)) return("#eceff1")
  if (pval < 0.05 && speedup > 1) return("#d9fdd3")
  ""
}
format_ms = function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(x < 0.01, format(x, scientific = TRUE, digits = 3, trim = TRUE), format(round(x, 2), nsmall = 2, trim = TRUE))
  )
}
dt[, EDI_Time_ms := format_ms(EDI_Time_ms)]
dt[, Canonical_Time_ms := format_ms(Canonical_Time_ms)]
dt[, Timing_Row_Color := mapply(row_bg_color, Speedup_Num, Timing_Pval, USE.NAMES = FALSE)]
dt[, Timing_Pval_Stars := format_pval_stars(Timing_Pval)]
dt[, Timing_Pval := format_pval(Timing_Pval)]
dt[, Response := as.character(Response)]

table_lines = c(
  "<table>",
  "  <thead>",
  "    <tr><th>Class</th><th>Response</th><th>EDI Time (ms)</th><th>Canonical Pkg</th><th>Canonical Func</th><th>Canonical Time (ms)</th><th>Speedup</th><th>Timing Pval</th><th></th></tr>",
  "  </thead>",
  "  <tbody>"
)
row_style = ifelse(
  nzchar(dt$Timing_Row_Color),
  paste0(" style=\"background-color: ", dt$Timing_Row_Color, ";\""),
  ""
)
table_rows = paste0(
  "    <tr", row_style,
  "><td>", dt$Class,
  "</td><td>", dt$Response,
  "</td><td>", dt$EDI_Time_ms,
  "</td><td>", dt$Canonical_Pkg,
  "</td><td>", dt$Canonical_Func,
  "</td><td>", dt$Canonical_Time_ms,
  "</td><td>", dt$Speedup,
  "</td><td>", dt$Timing_Pval,
  "</td><td>", dt$Timing_Pval_Stars,
  "</td></tr>"
)
table_lines = c(table_lines, table_rows, "  </tbody>", "</table>")

# Overwrite the previous Wald table
report_lines = readLines("package_metadata/benchmark_model_fits.md")
style_start = grep("^<style>$", report_lines)
if (length(style_start)) {
    report_lines = report_lines[seq_len(style_start[1] - 1L)]
}
wald_start = grep("## Wald Test Performance", report_lines)
if (length(wald_start)) {
    report_lines = report_lines[1:(wald_start-1)]
}
# Refresh the generated timestamp
ts_line = grep("^_Generated:", report_lines)
if (length(ts_line)) {
    report_lines[ts_line[1]] = paste0("_Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "_")
}

header = c(
  "## Wald Test Performance (Full Inference)",
  "",
  "This table compares the performance of **Full Inference** (Model Fit + Standard Error calculation + P-value derivation).",
  "Unlike the point-estimation table above, these results include the computational cost of the variance-covariance matrix (Hessian or Fisher Information) and the Wald test statistic calculation.",
  "All paths (EDI and Canonical) use a reduced sample size ($N=200$) for this full-inference benchmark to ensure iterative stability.",
  "EDI timings in this table correspond to fixed `iBCRD` design objects.",
  "**Stratified Cox Exception**: For `InferenceSurvivalStratCoxPHRegr`, the benchmark injects low-cardinality covariates before outcome generation so the row exercises a genuinely stratified Cox fit rather than the unstratified fallback.",
  "EDI regression models (Logistic, Poisson) are benchmarked using the **IRLS** optimizer for these Wald tests.",
  "**Note on Accessors**: EDI count-regression timings in this table explicitly call both `compute_wald_confidence_interval()` and `compute_wald_two_sided_pval()` so those rows remain Wald-only and cannot dispatch to bootstrap or likelihood-based fallback paths.",
  "**Solver-Only Prebuilds**: Benchmark setup prebuilds exposed observed-data design matrices, reduced design matrices, strata IDs, and other fixed working inputs outside the timed region when the implementation exposes those hooks. The timed region then measures the full-inference kernel on those fixed inputs.",
  "**Limitation**: Some canonical comparators only expose formula-based APIs rather than comparable low-level fit kernels. Those rows remain included, but their canonical timings may still contain formula/model-frame overhead beyond the numerical solver, variance, and p-value work itself.",
  paste0("**Timing Note**: All timings are medians over ", B_TIME, " warmed runs measured with adaptive batched `system.time`; paths below ", FAST_PATH_THRESHOLD_MS, " ms use `microbenchmark(times = ", FAST_PATH_MICROBENCH_REPS, ")` instead."),
  "**Timing P-Value**: `Timing Pval` reports a Welch two-sample t-test comparing the EDI and canonical timing replicate distributions for each row. The unlabeled final column marks thresholds with `***` for p < 0.001, `**` for p < 0.01, and `*` for p < 0.05.",
  "**Row Highlighting**: Light green rows indicate `Speedup > 1` and `Timing Pval < 0.05`; light grey rows indicate `NA` timing comparisons.",
  "",
  table_lines
)

METHODOLOGY_BLOCK = c(
  "## Garbage Collection and Cache Management",
  "",
  "To ensure that the benchmark results are highly precise, reproducible, and represent the actual computation speed of the numerical solvers, the benchmarking harness uses the following garbage collection and cache management strategies:",
  "",
  "### 1. Garbage Collection (GC) Filtering",
  "Garbage collection cycles run automatically by the R interpreter and can introduce significant, arbitrary pauses that skew timing measurements. To isolate the execution time of the code from R's GC overhead:",
  "* **GC Disabling**: We disable R's memory stress-testing mode using `gctorture(FALSE)` before running timing loops.",
  "* **Proactive Compaction**: In the `system.time()` path, we invoke `gc(verbose = FALSE)` immediately before timing each replicate. This starts the timer on a clean, compacted heap, minimizing the likelihood of triggering an automatic garbage collection cycle mid-replicate.",
  "* **Automatic Filtering**: In the microbenchmarking path, we utilize the `bench::mark()` engine with the `filter_gc = TRUE` parameter, which automatically tracks and discards timing iterations during which a garbage collection event occurred.",
  "",
  "### 2. Cold-Start Guarantee for EDI and Symmetric Warm-Up for Both Sides",
  "Both EDI and canonical timing expressions receive a single **validation/warm-up call** executed once before the calibration loop begins. This puts the machine code and working data into the instruction and data caches in the same warmed state for both sides, so the official timed replicates start on equal footing.",
  "",
  "EDI timings call exported C++ functions directly — no R6 objects are instantiated during benchmarking. As a result, **no R6 result caches exist to manage**. Each call to the C++ solver (e.g. `fast_logistic_regression_with_var_cpp`, `fast_ordinal_regression_cpp`) starts from a freshly zero-initialized parameter vector (or a model-specific data-driven initialization when `smart_cold_start = TRUE`). No prior-fit results are carried across timing repetitions, so every replication is a genuine cold start for the numerical optimizer."
)

writeLines(c(report_lines, "", header, "", METHODOLOGY_BLOCK, "", STYLE_BLOCK), "package_metadata/benchmark_model_fits.md")
cat("Wald Test Benchmark complete.\n")
