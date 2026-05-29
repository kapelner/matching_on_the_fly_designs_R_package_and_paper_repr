.libPaths(c(file.path(Sys.getenv("HOME"), "R", paste0(R.version$platform, "-library"), paste(R.version$major, sub("\\..*$", "", R.version$minor), sep = ".")), .libPaths()))
library(EDI)
library(microbenchmark)
library(data.table)
library(survival)

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

# Global Config
N_GLM = 1000
N_SURV = 500
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
    micro_time_ms = function() {
        mb_ms = microbenchmark(eval(expr, envir = env), times = fast_path_microbenchmark_reps)$time / 1e6
        mb_ms = mb_ms[is.finite(mb_ms) & mb_ms > 0]
        list(
            median_ms = if (length(mb_ms) == 0L) 0 else median(mb_ms),
            samples_ms = mb_ms,
            method = "microbenchmark"
        )
    }
    inner_reps = 1L
    batch_ms = 0
    while (inner_reps < max_inner_reps) {
        batch_ms = system.time(for (j in seq_len(inner_reps)) eval(expr, envir = env))[["elapsed"]] * 1000
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
        vals_ms[i] = system.time(for (j in seq_len(inner_reps)) eval(expr, envir = env))[["elapsed"]] * 1000 / inner_reps
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

# Builds a bare-metal timing environment and expression for EDI.
# All design matrices and fixed inputs are prebuilt here (outside the timed region).
# The returned expr calls the public (or internal) C++ function directly.
make_edi_bm = function(cls_name, d) {
    X_cov = d$X[, -1, drop = FALSE]
    colnames(X_cov) = paste0("x", seq_len(ncol(X_cov)))
    # GLM design (intercept + treatment + covariates)
    X_bm  = cbind(`(Intercept)` = 1, treatment = d$w, X_cov)
    # No-intercept design for ordinal and survival models
    X_ord = cbind(treatment = d$w, X_cov)
    y_bm  = as.numeric(d$y)
    w_bm  = as.integer(d$w)
    dead_bm = if (!is.null(d$dead)) as.numeric(d$dead) else NULL

    e = new.env(parent = globalenv())
    e$X_bm   = X_bm
    e$X_ord  = X_ord
    e$y_bm   = y_bm
    e$w_bm   = w_bm
    e$dead_bm = dead_bm

    expr = switch(cls_name,

        # --- GLM classes (intercept design) ---
        InferenceIncidLogRegr           = ,
        InferencePropFractionalLogit    = quote(fast_logistic_regression_cpp(X_bm, y_bm, estimate_only = TRUE)),

        InferenceContinOLS              = ,
        InferenceIncidRiskDiff          = quote(fast_ols_cpp(X_bm, y_bm)),

        InferenceCountPoisson           = ,
        InferenceCountQuasiPoisson      = ,
        InferenceCountRobustPoisson     = ,
        InferenceIncidModifiedPoisson   = quote(fast_poisson_regression_cpp(X_bm, y_bm, estimate_only = TRUE, optimization_alg = "irls")),

        InferenceCountNegBin            = quote(fast_neg_bin_cpp(X_bm, as.integer(y_bm), estimate_only = TRUE)),
        InferencePropBetaRegr           = quote(fast_beta_regression_cpp(X_bm, y_bm, estimate_only = TRUE)),
        InferenceContinRobustRegr       = quote(fast_robust_regression_cpp(X_bm, y_bm, estimate_only = TRUE)),
        InferenceIncidLogBinomial       = quote(fast_log_binomial_regression_cpp(X_bm, y_bm)),
        InferenceIncidProbitRegr        = quote(fast_probit_regression_cpp(X_bm, y_bm, estimate_only = TRUE)),
        InferenceIncidBinomialIdentityRiskDiff = quote(fast_identity_binomial_regression_cpp(X_bm, y_bm)),

        InferenceCountHurdlePoisson     = quote(fast_zero_augmented_poisson_cpp(X_bm, y_bm, X_bm, is_hurdle = TRUE,  estimate_only = TRUE)),
        InferenceCountZeroInflatedPoisson = quote(fast_zero_augmented_poisson_cpp(X_bm, y_bm, X_bm, is_hurdle = FALSE, estimate_only = TRUE)),
        InferenceCountZeroInflatedNegBin  = quote(fast_zinb_cpp(X_bm, X_bm, y_bm, estimate_only = TRUE)),
        InferenceCountHurdleNegBin      = quote(EDI:::fast_hurdle_negbin_cpp(X_bm, as.integer(y_bm), X_bm, estimate_only = TRUE)),

        # No EDI C++ kernel — pre-build R6 object outside timed region, time compute_estimate() only
        InferenceContinQuantileRegr = {
            des = DesignFixediBCRD$new(n = nrow(X_bm), response_type = "continuous")
            des$add_all_subjects_to_experiment(as.data.frame(X_cov))
            des$overwrite_all_subject_assignments(d$w)
            des$add_all_subject_responses(y_bm)
            e$inf_obj_quantile = tryCatch(InferenceContinQuantileRegr$new(des, smart_cold_start_default = FALSE), error = function(err) InferenceContinQuantileRegr$new(des))
            quote({
                inf_obj_quantile$.__enclos_env__$private$cached_mod = NULL
                inf_obj_quantile$.__enclos_env__$private$cached_values = list()
                inf_obj_quantile$compute_estimate(estimate_only = TRUE)
            })
        },

        # --- Ordinal classes (no-intercept design) ---
        InferenceOrdinalPropOddsRegr    = quote(fast_ordinal_regression_cpp(X_ord, y_bm, estimate_only = TRUE)),
        InferenceOrdinalAdjCatLogitRegr = quote(fast_adjacent_category_logit_cpp(X_ord, y_bm)),
        InferenceOrdinalContRatioRegr   = quote(fast_continuation_ratio_regression_cpp(X_ord, y_bm)),
        InferenceOrdinalOrderedProbitRegr = quote(fast_ordinal_probit_regression_cpp(X_ord, y_bm, estimate_only = TRUE)),
        InferenceOrdinalCloglogRegr     = quote(fast_ordinal_cloglog_regression_cpp(X_ord, y_bm, estimate_only = TRUE)),
        InferenceOrdinalCauchitRegr     = quote(fast_ordinal_cauchit_regression_cpp(X_ord, y_bm, estimate_only = TRUE)),

        # --- Survival classes (no-intercept design) ---
        # CoxPH and StratCox use survival::coxph.fit internally — no EDI C++ kernel; use R6 approach
        InferenceSurvivalCoxPHRegr = {
            des = DesignFixediBCRD$new(n = nrow(X_ord), response_type = "survival")
            des$add_all_subjects_to_experiment(as.data.frame(X_cov))
            des$overwrite_all_subject_assignments(d$w)
            des$add_all_subject_responses(y_bm, deads = dead_bm)
            e$inf_obj_coxph = tryCatch(InferenceSurvivalCoxPHRegr$new(des, smart_cold_start_default = FALSE), error = function(err) InferenceSurvivalCoxPHRegr$new(des))
            quote({
                inf_obj_coxph$.__enclos_env__$private$cached_mod = NULL
                inf_obj_coxph$.__enclos_env__$private$cached_values = list()
                inf_obj_coxph$compute_estimate(estimate_only = TRUE)
            })
        },
        InferenceSurvivalStratCoxPHRegr = {
            des = DesignFixediBCRD$new(n = nrow(X_ord), response_type = "survival")
            des$add_all_subjects_to_experiment(as.data.frame(X_cov))
            des$overwrite_all_subject_assignments(d$w)
            des$add_all_subject_responses(y_bm, deads = dead_bm)
            e$inf_obj_strat_cox = tryCatch(InferenceSurvivalStratCoxPHRegr$new(des, smart_cold_start_default = FALSE), error = function(err) InferenceSurvivalStratCoxPHRegr$new(des))
            quote({
                inf_obj_strat_cox$.__enclos_env__$private$cached_mod = NULL
                inf_obj_strat_cox$.__enclos_env__$private$cached_values = list()
                inf_obj_strat_cox$compute_estimate(estimate_only = TRUE)
            })
        },
        InferenceSurvivalWeibullRegr    = quote(fast_weibull_regression_cpp(X_ord, y_bm, dead_bm, estimate_only = TRUE)),
        InferenceSurvivalLogRank        = quote(EDI:::fast_logrank_stats_cpp(w_bm, y_bm, dead_bm)),
        InferenceSurvivalKMDiff         = quote(EDI:::get_survival_stat_diff(y_bm, dead_bm, w_bm, "median")),
        InferenceSurvivalRestrictedMeanDiff = quote(EDI:::get_survival_stat_diff(y_bm, dead_bm, w_bm, "restricted_mean")),

        InferenceAllSimpleWilcox = quote(EDI:::wilcox_hl_point_estimate_cpp(w_bm, y_bm)),

        # Lin regression: pre-build centered interaction design; call lm.fit directly


        # --- GComp classes: fast_logistic_regression_cpp + C++ marginalisation ---
        InferencePropGCompMeanDiff = quote({
            b = fast_logistic_regression_cpp(X_bm, y_bm, estimate_only = TRUE)$b
            EDI:::gcomp_fractional_logit_point_estimate_cpp(X_bm, b, 2L)$md
        }),
        InferenceIncidGCompRiskDiff = quote({
            b = fast_logistic_regression_cpp(X_bm, y_bm, estimate_only = TRUE)$b
            EDI:::gcomp_logistic_point_estimate_cpp(X_bm, b, 2L)$md
        }),
        InferenceIncidGCompRiskRatio = quote({
            b = fast_logistic_regression_cpp(X_bm, y_bm, estimate_only = TRUE)$b
            res = EDI:::gcomp_logistic_point_estimate_cpp(X_bm, b, 2L)
            res$mean1 / res$mean0
        }),
        InferenceOrdinalGCompMeanDiff = quote({
            fit = fast_ordinal_regression_cpp(X_ord, y_bm, estimate_only = TRUE)
            EDI:::gcomp_ordinal_proportional_odds_post_fit_cpp(
                X_fit     = X_ord,
                coef_hat  = as.numeric(fit$b),
                alpha_hat = as.numeric(fit$alpha),
                j_treat   = 1L
            )$md
        }),

        NULL  # unknown class
    )

    if (is.null(expr)) return(NULL)
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

# --- Mapping Table (Targeted subset requested by user) ---
# library() calls are used to ensure it crashes if pkgs are missing.
bench_specs = list(
    list(cls = "InferenceIncidLogRegr", pkg = "stats", func = "glm.fit", expr = quote(glm.fit(x = X_can, y = df$y, family = binomial()))),
    list(cls = "InferenceContinOLS", pkg = "stats", func = "lm.fit", expr = quote(lm.fit(x = X_can, y = df$y))),
    list(cls = "InferenceCountPoisson", pkg = "stats", func = "glm.fit", expr = quote(glm.fit(x = X_can, y = df$y, family = poisson()))),
    list(cls = "InferenceSurvivalCoxPHRegr", pkg = "survival", func = "coxph.fit(breslow)", expr = quote({
        x_vars = c("treatment", grep("^x", names(df), value=T))
        survival::coxph.fit(x = as.matrix(df[, x_vars, drop=F]), y = survival::Surv(df$y, df$dead), strata=NULL, offset=NULL, init=NULL, control=survival::coxph.control(), weights=NULL, method="breslow", rownames=as.character(1:nrow(df)))
    })),
    list(cls = "InferenceCountNegBin", pkg = "MASS", func = "glm.nb", expr = quote(MASS::glm.nb(y ~ treatment + x1 + x2 + x3 + x4, data = df))),
    list(cls = "InferencePropBetaRegr", pkg = "betareg", func = "betareg.fit", expr = quote(betareg::betareg.fit(x = X_can, y = df$y))),
    list(cls = "InferenceOrdinalPropOddsRegr", pkg = "ordinal", func = "clm", expr = quote(ordinal::clm(factor(y, ordered=T) ~ treatment + x1 + x2 + x3 + x4, data = df))),
    list(cls = "InferenceCountHurdlePoisson", pkg = "pscl", func = "hurdle", expr = quote(pscl::hurdle(y ~ treatment + x1 + x2 + x3 + x4, data = df))),
    list(cls = "InferenceCountZeroInflatedPoisson", pkg = "pscl", func = "zeroinfl", expr = quote(pscl::zeroinfl(y ~ treatment + x1 + x2 + x3 + x4, data = df))),
    list(cls = "InferenceCountZeroInflatedNegBin", pkg = "pscl", func = "zeroinfl(nb)", expr = quote(pscl::zeroinfl(y ~ treatment + x1 + x2 + x3 + x4, data = df, dist="negbin"))),
    list(cls = "InferenceCountHurdleNegBin", pkg = "pscl", func = "hurdle(nb)", expr = quote(pscl::hurdle(y ~ treatment + x1 + x2 + x3 + x4, data = df, dist="negbin"))),
    list(cls = "InferenceCountQuasiPoisson", pkg = "stats", func = "glm.fit(quasi)", expr = quote(glm.fit(x = X_can, y = df$y, family = quasipoisson()))),
    list(cls = "InferenceSurvivalWeibullRegr", pkg = "survival", func = "survreg", expr = quote(survival::survreg(survival::Surv(y, dead) ~ treatment + x1 + x2 + x3 + x4, data = df, dist="weibull"))),
    list(cls = "InferenceContinRobustRegr", pkg = "MASS", func = "rlm(MM)", expr = quote(MASS::rlm(x = X_can, y = df$y, method = "MM"))),
    list(cls = "InferenceContinQuantileRegr", pkg = "quantreg", func = "rq.fit", expr = quote(quantreg::rq.fit(x = X_can, y = df$y))),
    list(cls = "InferencePropFractionalLogit", pkg = "stats", func = "glm.fit(quasi)", expr = quote(glm.fit(x = X_can, y = df$y, family=quasibinomial()))),
    list(cls = "InferenceIncidLogBinomial", pkg = "stats", func = "glm.fit(log)", expr = quote({
        glm.fit(x = X_can, y = df$y, family=binomial(link="log"), start=c(-2, rep(0, ncol(X_can)-1)))
    })),
    list(cls = "InferenceIncidProbitRegr", pkg = "stats", func = "glm.fit(probit)", expr = quote(glm.fit(x = X_can, y = df$y, family=binomial(link="probit")))),
    list(cls = "InferenceIncidBinomialIdentityRiskDiff", pkg = "stats", func = "glm.fit(ident)", expr = quote(glm.fit(x = X_can, y = df$y, family=binomial(link="identity"), start=c(0.5, rep(0, ncol(X_can)-1))))),
    list(cls = "InferenceOrdinalAdjCatLogitRegr", pkg = "VGAM", func = "vglm(acat)", expr = quote(VGAM::vglm(factor(y, ordered=T) ~ treatment + x1 + x2 + x3 + x4, VGAM::acat(), data=df))),
    list(cls = "InferenceOrdinalContRatioRegr", pkg = "VGAM", func = "vglm(cratio)", expr = quote(VGAM::vglm(factor(y, ordered=T) ~ treatment + x1 + x2 + x3 + x4, VGAM::cratio(), data=df))),
    list(cls = "InferenceOrdinalOrderedProbitRegr", pkg = "ordinal", func = "clm(probit)", expr = quote(ordinal::clm(factor(y, ordered=T) ~ treatment + x1 + x2 + x3 + x4, data = df, link = "probit"))),
    list(cls = "InferenceOrdinalCloglogRegr", pkg = "ordinal", func = "clm(cloglog)", expr = quote(ordinal::clm(factor(y, ordered=T) ~ treatment + x1 + x2 + x3 + x4, data = df, link = "cloglog"))),
    list(cls = "InferenceOrdinalCauchitRegr", pkg = "ordinal", func = "clm(cauchit)", expr = quote(ordinal::clm(factor(y, ordered=T) ~ treatment + x1 + x2 + x3 + x4, data = df, link = "cauchit"))),
    list(cls = "InferenceSurvivalLogRank", pkg = "survival", func = "survdiff", expr = quote(survival::survdiff(survival::Surv(y, dead) ~ treatment, data = df))),
    list(cls = "InferenceAllSimpleWilcox", pkg = "stats", func = "HL median pairwise diff", expr = quote({
        y_t = df$y[df$treatment == 1]
        y_c = df$y[df$treatment == 0]
        stats::median(as.numeric(outer(y_t, y_c, "-")))
    }), scale = 0.5),
    list(cls = "InferenceIncidModifiedPoisson", pkg = "stats", func = "glm.fit(modified)", expr = quote(glm.fit(x = X_can, y = df$y, family = poisson()))),
    list(cls = "InferenceSurvivalStratCoxPHRegr", pkg = "survival", func = "coxph.fit(strat)", expr = quote({
        survival::coxph.fit(
            x = X_can_strat,
            y = survival::Surv(y_can_strat, dead_can_strat),
            strata = strata_can,
            offset = NULL,
            init = NULL,
            control = survival::coxph.control(),
            weights = NULL,
            method = "breslow",
            rownames = as.character(seq_along(y_can_strat))
        )
    })),
    list(cls = "InferenceIncidRiskDiff", pkg = "stats", func = "lm.fit(LPM)", expr = quote({
        fit = lm.fit(x = X_can, y = df$y)
        as.numeric(stats::coef(fit)[2L])
    })),
    list(cls = "InferenceSurvivalKMDiff", pkg = "survival", func = "survfit(median)", expr = quote({
        fit = survival::survfit(survival::Surv(df$y, df$dead) ~ df$treatment)
        q = stats::quantile(fit, 0.5)
        as.numeric(q$quantile[2] - q$quantile[1])
    })),
    list(cls = "InferenceSurvivalRestrictedMeanDiff", pkg = "survival", func = "survfit(rmean)", expr = quote({
        fit = survival::survfit(survival::Surv(df$y, df$dead) ~ df$treatment)
        res = survival:::survmean(fit, rmean="common")
        # The column name can be "*rmean" or just "rmean" depending on survival version
        col_idx = grep("rmean", colnames(res$matrix))
        rmst = res$matrix[, col_idx]
        as.numeric(rmst[2] - rmst[1])
    })),
    list(cls = "InferenceCountRobustPoisson", pkg = "stats", func = "glm.fit", expr = quote(glm.fit(x = X_can, y = df$y, family = poisson())))
)


# Paths with no canonical mapping (None)
no_can_specs = list(
    list(cls = "InferencePropGCompMeanDiff", pkg = "stats", func = "glm.fit(quasi)+gcomp", expr = quote({
        mod = glm.fit(x = X_can, y = df$y, family = quasibinomial())
        X1 = X_can; X0 = X_can
        X1[, "treatment"] = 1
        X0[, "treatment"] = 0
        mean(stats::plogis(drop(X1 %*% mod$coefficients))) - mean(stats::plogis(drop(X0 %*% mod$coefficients)))
    })),
    list(cls = "InferenceIncidGCompRiskRatio", pkg = "stats", func = "glm.fit+gcomp(RR)", expr = quote({
        mod = glm.fit(x = X_can, y = df$y, family = binomial())
        X1 = X_can; X0 = X_can
        X1[, "treatment"] = 1
        X0[, "treatment"] = 0
        risk1 = mean(stats::plogis(drop(X1 %*% mod$coefficients)))
        risk0 = mean(stats::plogis(drop(X0 %*% mod$coefficients)))
        risk1 / risk0
    })),
    list(cls = "InferenceIncidGCompRiskDiff", pkg = "stats", func = "glm.fit+gcomp(RD)", expr = quote({
        mod = glm.fit(x = X_can, y = df$y, family = binomial())
        X1 = X_can; X0 = X_can
        X1[, "treatment"] = 1
        X0[, "treatment"] = 0
        mean(stats::plogis(drop(X1 %*% mod$coefficients))) - mean(stats::plogis(drop(X0 %*% mod$coefficients)))
    })),
    list(cls = "InferenceOrdinalGCompMeanDiff", pkg = "ordinal", func = "clm+gcomp", expr = quote({
        mod = ordinal::clm(factor(y, ordered = TRUE) ~ treatment + x1 + x2 + x3 + x4, data = df)
        df1 = df; df0 = df
        df1$treatment = 1
        df0$treatment = 0
        p1 = predict(mod, newdata = df1, type = "prob")$fit
        p0 = predict(mod, newdata = df0, type = "prob")$fit
        score = seq_len(ncol(p1))
        mean(drop(p1 %*% score)) - mean(drop(p0 %*% score))
    }))
)
bench_specs = c(bench_specs, no_can_specs)

# --- Benchmark Runner ---
results = list()

run_one = function(spec) {
    cls_name = spec$cls
    cat(sprintf("Benchmarking %s...\n", cls_name))
    
    # Improved heuristic mapping: Ordinal MUST come before Prop/PropOdds
    resp_type = "continuous"
    family = "continuous"
    if (grepl("Ordinal|AdjCat|ContRatio|Stereotype|Jonckheere|Ridit|Sign", cls_name)) { resp_type = "ordinal"; family = "ordinal" }
    else if (grepl("Incid|Binomial|Wald|CMH|Fisher|Zhang|Robins|Newcombe|Nurminen", cls_name)) { resp_type = "incidence"; family = "logistic" }
    else if (grepl("Count|Poisson|NegBin|ZINB|ZAP|Hurdle", cls_name)) { resp_type = "count"; family = "poisson" }
    else if (grepl("Prop|Beta|ZOIB|Fractional", cls_name)) { resp_type = "proportion"; family = "beta" }
    else if (grepl("Survival|Cox|Weibull|KM|Rank|LogRank|Gehan|RMST|RMDiff|LWACox|Clayton", cls_name)) { resp_type = "survival"; family = "cox" }
    else if (grepl("Wilcox|MeanDiff|OLS|Lin|Robust|Quantile|Bai", cls_name)) { resp_type = "continuous"; family = "continuous" }
    
    if (cls_name == "InferenceIncidLogBinomial") family = "log-binomial"
    
    scale = if (!is.null(spec$scale)) spec$scale else 1.0
    fast_path_microbenchmark_reps = if (!is.null(spec$fast_path_microbenchmark_reps)) spec$fast_path_microbenchmark_reps else FAST_PATH_MICROBENCH_REPS
    b_time = if (!is.null(spec$b_time_override)) as.integer(spec$b_time_override) else B_TIME
    n = round(N_GLM * scale)
    if (grepl("Survival", cls_name)) n = round(N_SURV * scale)
    
    d = generate_data(n = n, family = family)
    if (identical(cls_name, "InferenceSurvivalStratCoxPHRegr")) {
        d = make_true_stratified_survival_data(d)
    }
    
    # Timing EDI (bare metal: call exported C++ functions directly with pre-built inputs)
    timing_edi = tryCatch({
        bm = make_edi_bm(cls_name, d)
        if (is.null(bm)) stop("no bare metal mapping for this class")
        eval(bm$expr, envir = bm$env)  # validation run
        collect_timing_ms(bm$expr, times = b_time, env = bm$env, fast_path_microbenchmark_reps = fast_path_microbenchmark_reps)
    }, error = function(e) {
        cat("  EDI Error:", e$message, "\n")
        list(median_ms = NA_real_, samples_ms = numeric(0))
    })

    # Timing Canonical (Crashes if pkg missing)
    timing_can = list(median_ms = NA_real_, samples_ms = numeric(0))
    if (!is.null(spec$expr) && spec$pkg != "None") {
        library(spec$pkg, character.only = TRUE)
        X_cols = d$X[,-1,drop=F]
        colnames(X_cols) = paste0("x", 1:ncol(X_cols))
        df = data.frame(y = d$y, treatment = d$w, g = factor(rep(1:10, length.out = n)), dead = if(exists("dead", d)) d$dead else 1)
        df = cbind(df, X_cols)
        X_can = cbind(`(Intercept)` = 1, treatment = d$w, as.matrix(X_cols))
        if (identical(cls_name, "InferenceSurvivalStratCoxPHRegr")) {
            strat_inputs = build_strat_cox_canonical_inputs(d)
            X_can_strat = strat_inputs$X
            y_can_strat = strat_inputs$y
            dead_can_strat = strat_inputs$dead
            strata_can = strat_inputs$strata
        }

        timing_can = tryCatch(
            collect_timing_ms(spec$expr, times = b_time, fast_path_microbenchmark_reps = fast_path_microbenchmark_reps),
            error = function(e) {
                cat("  Canonical Error:", e$message, "\n")
                list(median_ms = NA_real_, samples_ms = numeric(0))
            }
        )
    }
    
    results[[length(results) + 1]] <<- data.table(
        Class = cls_name, Response = resp_type, EDI_Time_ms = timing_edi$median_ms,
        Canonical_Pkg = spec$pkg, Canonical_Func = spec$func, Canonical_Time_ms = timing_can$median_ms,
        Speedup = if (!is.na(timing_can$median_ms) && !is.na(timing_edi$median_ms) && timing_edi$median_ms > 0) timing_can$median_ms / timing_edi$median_ms else NA_real_,
        Timing_Pval = timing_ttest_pval(timing_edi$samples_ms, timing_can$samples_ms)
    )
}

# Run the requested list
# First unique classes to avoid duplicates
unique_specs = bench_specs[!duplicated(sapply(bench_specs, `[[`, "cls"))]

for (i in seq_along(unique_specs)) {
    cat(sprintf("[%d/%d] ", i, length(unique_specs)))
    run_one(unique_specs[[i]])
}

# Finalize
dt = rbindlist(results)
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
table_rows = mapply(function(cls, resp, edi, pkg, func, can, speed, pval, stars, bg) {
  style = if (nzchar(bg)) paste0(" style=\"background-color: ", bg, ";\"") else ""
  paste0(
    "    <tr", style, "><td>", cls, "</td><td>", resp, "</td><td>", edi, "</td><td>", pkg,
    "</td><td>", func, "</td><td>", can, "</td><td>", speed, "</td><td>", pval, "</td><td>", stars, "</td></tr>"
  )
}, dt$Class, dt$Response, dt$EDI_Time_ms, dt$Canonical_Pkg, dt$Canonical_Func,
dt$Canonical_Time_ms, dt$Speedup, dt$Timing_Pval, dt$Timing_Pval_Stars, dt$Timing_Row_Color,
SIMPLIFY = TRUE, USE.NAMES = FALSE)
table_lines = c(table_lines, table_rows, "  </tbody>", "</table>")

cmd_config = function(var) {
  out = tryCatch(
    system2(file.path(R.home("bin"), "R"), c("CMD", "config", var), stdout = TRUE, stderr = TRUE),
    error = function(e) character(0)
  )
  out = out[nzchar(out)]
  if (length(out) == 0L || any(grepl("^ERROR:", out))) "unavailable" else paste(out, collapse = " ")
}

env_flag = function(var, default) {
  val = Sys.getenv(var, unset = NA_character_)
  if (is.na(val) || !nzchar(val)) default else val
}

compile_context_lines = function() {
  if (exists("edi_build_info_cpp", mode = "function")) {
    info = tryCatch(edi_build_info_cpp(), error = function(e) NULL)
    if (!is.null(info)) {
      so_path = file.path(system.file("libs", package = "EDI"), paste0("EDI", .Platform$dynlib.ext))
      so_info = if (file.exists(so_path)) file.info(so_path) else NULL
      return(c(
        "## Compilation Context",
        "",
        "These rows are read from build metadata compiled into the loaded `EDI` shared object via `edi_build_info_cpp()`.",
        "",
        "**Compilation warning:** EDI model-fit timings are sensitive to the compiler flags used to build the loaded `EDI.so`. If EDI is compiled without the proper optimized flags, or with flags that are known to degrade these kernels such as problematic LTO builds, the benchmark can show substantial performance regressions that reflect the binary build rather than the modeling algorithms.",
        "",
        paste0("*   **EDI shared object:** `", so_path, "`"),
        paste0("*   **EDI shared object mtime:** `", if (!is.null(so_info)) format(so_info$mtime) else "unknown", "`"),
        paste0("*   **Capture method:** `", info$capture_method, "`"),
        paste0("*   **Build timestamp:** `", info$build_timestamp, "`"),
        paste0("*   **Build host:** `", info$build_host, "`"),
        paste0("*   **R version at build:** `", info$r_version, "`"),
        paste0("*   **R `CXX20` at build:** `", info$r_cxx20, "`"),
        paste0("*   **R `CXX20STD` at build:** `", info$r_cxx20std, "`"),
        paste0("*   **R `CXX20FLAGS` at build:** `", info$r_cxx20flags, "`"),
        paste0("*   **R `SHLIB_OPENMP_CXXFLAGS` at build:** `", info$r_shlib_openmp_cxxflags, "`"),
        paste0("*   **Build env at build:** `EDI_PORTABLE=", info$env_edi_portable,
               "`, `EDI_DISABLE_VECTORIZATION=", info$env_edi_disable_vectorization,
               "`, `EDI_NATIVE_SPEED=", info$env_edi_native_speed,
               "`, `EDI_NATIVE_LTO=", info$env_edi_native_lto, "`"),
        paste0("*   **Package `PKG_CPPFLAGS` at build:** `", info$pkg_cppflags, "`"),
        paste0("*   **Package `PKG_CXXFLAGS` at build:** `", info$pkg_cxxflags, "`"),
        paste0("*   **Package `PKG_LIBS` at build:** `", info$pkg_libs, "`"),
        paste0("*   **Compiler reported by binary:** `", info$compiler, "`"),
        paste0("*   **Compiler optimization macro enabled:** `", info$compiler_optimize_macro, "`"),
        paste0("*   **Compiler fast-math macro enabled:** `", info$compiler_fast_math_macro, "`"),
        paste0("*   **Eigen vectorization disabled macro enabled:** `", info$eigen_dont_vectorize_macro, "`"),
        ""
      ))
    }
  }

  edi_portable = env_flag("EDI_PORTABLE", "0")
  edi_disable_vectorization = env_flag("EDI_DISABLE_VECTORIZATION", "0")
  edi_native_speed = env_flag("EDI_NATIVE_SPEED", "1")
  edi_native_lto = env_flag("EDI_NATIVE_LTO", "0")

  pkg_cppflags = "-I../inst/include"
  pkg_cxxflags = "$(SHLIB_OPENMP_CXXFLAGS) -DNDEBUG -DEIGEN_NO_DEBUG -Wno-ignored-attributes"
  pkg_libs = "$(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -ltbb12 -fstack-protector"

  if (identical(edi_portable, "0")) {
    pkg_cxxflags = paste(pkg_cxxflags, "-march=native -mtune=native")
  }
  if (identical(edi_disable_vectorization, "1")) {
    pkg_cppflags = paste(pkg_cppflags, "-DEIGEN_DONT_VECTORIZE -DEIGEN_UNALIGNED_VECTORIZE=0")
    pkg_cxxflags = paste(pkg_cxxflags, "-fno-tree-vectorize")
  }
  if (identical(edi_native_speed, "1")) {
    pkg_cxxflags = paste(pkg_cxxflags, "override CXXFLAGS+=-O3")
  }
  if (identical(edi_native_lto, "1")) {
    pkg_cxxflags = paste(pkg_cxxflags, "-flto")
    pkg_libs = paste(pkg_libs, "-flto")
  } else {
    pkg_cxxflags = paste(pkg_cxxflags, "-fno-lto")
  }

  c(
    "## Compilation Context",
    "",
    "The installed EDI package does not expose compiled-in build metadata via `edi_build_info_cpp()`. These rows therefore record only the current toolchain and `EDI/src/Makevars`-derived context visible to the benchmark process; they are not proof of the flags used to build the loaded `EDI.so`.",
    "",
    "**Compilation warning:** EDI model-fit timings are sensitive to the compiler flags used to build the loaded `EDI.so`. If EDI is compiled without the proper optimized flags, or with flags that are known to degrade these kernels such as problematic LTO builds, the benchmark can show substantial performance regressions that reflect the binary build rather than the modeling algorithms.",
    "",
    paste0("*   **EDI library path:** `", system.file(package = "EDI"), "`"),
    paste0("*   **R version:** `", R.version.string, "`"),
    paste0("*   **R `CXX20`:** `", cmd_config("CXX20"), "`"),
    paste0("*   **R `CXX20STD`:** `", cmd_config("CXX20STD"), "`"),
    paste0("*   **R `CXX20FLAGS`:** `", cmd_config("CXX20FLAGS"), "`"),
    paste0("*   **R `SHLIB_OPENMP_CXXFLAGS`:** `", cmd_config("SHLIB_OPENMP_CXXFLAGS"), "`"),
    paste0("*   **Build env:** `EDI_PORTABLE=", edi_portable,
           "`, `EDI_DISABLE_VECTORIZATION=", edi_disable_vectorization,
           "`, `EDI_NATIVE_SPEED=", edi_native_speed,
           "`, `EDI_NATIVE_LTO=", edi_native_lto, "`"),
    paste0("*   **Effective package `PKG_CPPFLAGS`:** `", pkg_cppflags, "`"),
    paste0("*   **Effective package `PKG_CXXFLAGS`:** `", pkg_cxxflags, "`"),
    paste0("*   **Effective package `PKG_LIBS`:** `", pkg_libs, "`"),
    ""
  )
}

report = c(
  "# EDI Exhaustive C++ Model Fit Benchmarks",
  "",
  paste0("_Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "_"),
  "",
  "This report compares the performance of EDI's Rcpp-optimized model fitting paths against **low-level** canonical R implementations (e.g., `glm.fit`, `lm.fit`, `coxph.fit`) where possible.",
  "",
  compile_context_lines(),
  "## Benchmark Dataset Specification",
  "",
  "All benchmarks were performed on a synthetic clinical-trial-scale dataset generated for each response type. The data generation process ensures numerical stability and fair solver comparison by using the following parameters:",
  "",
  "*   **Sample Size ($N$):** 1,000 subjects for most models; 500 subjects for survival models. Exact and trend tests may use smaller scaled samples (N=100-500) as noted in the results.",
  "*   **Predictors ($p$):** 5 total predictors, including a global intercept, a balanced binary treatment assignment from fixed `iBCRD`, and 4 continuous covariates ($X \\sim \\text{Normal}(0, 1)$).",
  "*   **Effect Sizes:** Covariate coefficients are sampled from $\\text{Normal}(0, 0.5)$. The treatment coefficient is set to 0.5 in the linear predictor so the benchmarked treatment effect is meaningfully separated from zero.",
  "*   **EDI Design Template:** EDI benchmark objects are instantiated on a fixed `iBCRD` design.",
  "*   **Response Generation:**",
  "    *   **Continuous:** Linear model with additive $\\text{Normal}(0, 0.5)$ noise.",
  "    *   **Incidence:** Binary outcomes via a Logistic link.",
  "    *   **Count:** Integer outcomes via Poisson or Negative Binomial distributions with an exponential link.",
  "    *   **Proportion:** Continuous outcomes in $(0, 1)$ via a Beta distribution with a logit link.",
  "    *   **Survival:** Exponentially distributed event times with approximately 20% random censoring.",
  "    *   **Ordinal:** 3-level categorical outcomes generated from the same ordinal construction used elsewhere in the benchmark suite.",
  "*   **Stratified Cox Exception:** For `InferenceSurvivalStratCoxPHRegr`, the benchmark injects low-cardinality covariates before outcome generation so the row exercises a genuinely stratified Cox fit rather than the unstratified fallback.",
  "",
  "## Methodology",
  "",
  "*   **Bare Metal EDI Timing:** EDI rows call the exported C++ functions directly (e.g., `fast_logistic_regression_cpp`, `fast_ordinal_regression_cpp`) with all design matrices and fixed inputs pre-built outside the timed region. There is no R6 object instantiation, no cached state management, no warm start storage, and no standard error computation in the timed region — only the raw numerical solver.",
  "*   **Apples-to-Apples Canonical Timing:** Canonical R timings likewise call the lowest-level publicly exposed interfaces (e.g., `glm.fit`, `lm.fit`, `coxph.fit`) with pre-built design matrices. If a canonical package exposes no low-level function, the formula-based API is used instead.",
  "*   **Low-Level Comparison:** Both EDI and canonical timings are measured on pre-built numeric matrices, removing formula parsing, model-frame construction, and R6/S3/S4 dispatch overhead from the timed region wherever the API permits.",
  "*   **Limitation:** Some canonical comparators only expose formula-based APIs. Those rows remain included but their canonical timings carry formula/model-frame overhead not present in the EDI bare-metal timing.",
  paste0("*   **Averaging:** All timings are medians over ", B_TIME, " cold estimate-only timing samples measured with adaptive batched `system.time`; paths below ", FAST_PATH_THRESHOLD_MS, " ms use `microbenchmark(times = ", FAST_PATH_MICROBENCH_REPS, ")` instead."),
  "*   **Timing P-Value:** `Timing Pval` reports a Welch two-sample t-test comparing the EDI and canonical timing replicate distributions for each row. The unlabeled final column marks thresholds with `***` for p < 0.001, `**` for p < 0.01, and `*` for p < 0.05.",
  "*   **Row Highlighting:** Light green rows indicate `Speedup > 1` and `Timing Pval < 0.05`; light grey rows indicate `NA` timing comparisons.",
  "*   **Constraints**: Matched-pair/KK and highly custom paths are excluded as per user request.",
  "",
  "## Results",
  "",
  table_lines
)
writeLines(c(report, "", STYLE_BLOCK), "package_metadata/benchmark_model_fits.md")
cat("Benchmark complete.\n")
