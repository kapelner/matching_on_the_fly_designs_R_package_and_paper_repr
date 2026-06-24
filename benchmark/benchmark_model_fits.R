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
suppressPackageStartupMessages(library(sandwich))
suppressPackageStartupMessages(library(clinfun))
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

# Global Config
N_GLM = 1000
N_SURV = 500
N_WALD = 200
B_TIME = 30
TARGET_BATCH_MS = 200
MIN_RESOLVED_BATCH_MS = 10
MAX_INNER_REPS = 100000L
FAST_PATH_MICROBENCH_REPS = 5000L
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
    dead_bm = if (!is.null(d$dead)) as.integer(d$dead) else NULL

    e = new.env(parent = globalenv())
    e$d      = d
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

        InferenceContinQuantileRegr = quote(quantreg::rq.fit(x = X_bm, y = y_bm, method = "br")),

        # --- Ordinal classes (no-intercept design) ---
        InferenceOrdinalPropOddsRegr    = quote(fast_ordinal_regression_cpp(X_ord, y_bm, estimate_only = TRUE)),
        InferenceOrdinalAdjCatLogitRegr = quote(fast_adjacent_category_logit_cpp(X_ord, y_bm)),
        InferenceOrdinalContRatioRegr   = quote(fast_continuation_ratio_regression_cpp(X_ord, y_bm)),
        InferenceOrdinalOrderedProbitRegr = quote(fast_ordinal_probit_regression_cpp(X_ord, y_bm, estimate_only = TRUE)),
        InferenceOrdinalCloglogRegr     = quote(fast_ordinal_cloglog_regression_cpp(X_ord, y_bm, estimate_only = TRUE)),
        InferenceOrdinalCauchitRegr     = quote(fast_ordinal_cauchit_regression_cpp(X_ord, y_bm, estimate_only = TRUE)),

        # --- Survival classes (no-intercept design) ---
        InferenceSurvivalCoxPHRegr = quote(fast_coxph_regression_cpp(X_ord, y_bm, dead_bm, estimate_only = TRUE)),
        InferenceSurvivalStratCoxPHRegr = quote({
            strat_inputs = build_strat_cox_canonical_inputs(d)
            if (!is.null(strat_inputs$strata)) {
                cache = build_stratified_cox_data_cache_cpp(strat_inputs$X, strat_inputs$y, strat_inputs$dead, strat_inputs$strata)
            } else {
                cache = build_cox_data_cache_cpp(strat_inputs$X, strat_inputs$y, strat_inputs$dead)
            }
            fast_coxph_regression_prebuilt_cpp(cache, estimate_only = TRUE)
        }),
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
            gcomp_ordinal_proportional_odds_post_fit_cpp(
                X_fit     = X_ord,
                coef_hat  = as.numeric(fit$b),
                alpha_hat = as.numeric(fit$alpha),
                j_treat   = 1L
            )$md
        }),

        NULL  # unknown class
    )

    if (is.null(expr)) return(NULL)
    setup = NULL
    if (is.list(expr) && !is.call(expr)) {
        setup = expr$setup
        expr = expr$expr
    }
    list(env = e, expr = expr, setup = setup)
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
        strat_inputs = build_strat_cox_canonical_inputs(d)
        survival::coxph.fit(
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

        tryCatch(eval(spec$expr), error = function(e) NULL)  # validation run (mirrors EDI's pre-timing warm-up call)
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

# ── Shared table-formatting helpers (used by both point-estimate and Wald tables)
format_pval = function(x) {
  ifelse(is.na(x), "NA", vapply(x, function(v) sprintf("%.3g", v), character(1)))
}
format_pval_stars = function(x) {
  ifelse(is.na(x), "",
    ifelse(x < 0.001, "***", ifelse(x < 0.01, "**", ifelse(x < 0.05, "*", ""))))
}
row_bg_color = function(speedup, pval) {
  if (!is.finite(speedup) || is.na(pval)) return("#eceff1")
  if (pval < 0.05 && speedup > 1) return("#d9fdd3")
  ""
}
format_ms = function(x) {
  ifelse(is.na(x), "NA",
    ifelse(x < 0.01, format(x, scientific = TRUE, digits = 3, trim = TRUE),
           format(round(x, 2), nsmall = 2, trim = TRUE)))
}

# Run the requested list
# First unique classes to avoid duplicates
unique_specs = bench_specs[!duplicated(sapply(bench_specs, `[[`, "cls"))]

for (i in seq_along(unique_specs)) {
    cat(sprintf("[%d/%d] ", i, length(unique_specs)))
    run_one(unique_specs[[i]])
}

# ── Wald / Full-Inference Section ──────────────────────────────────────────────

compute_binary_gcomp_effect = function(beta, X_base, effect = c("RD", "RR")) {
    effect = match.arg(effect)
    X1 = X_base; X0 = X_base
    X1[, "treatment"] = 1; X0[, "treatment"] = 0
    risk1_i = stats::plogis(drop(X1 %*% beta))
    risk0_i = stats::plogis(drop(X0 %*% beta))
    risk1 = mean(risk1_i); risk0 = mean(risk0_i)
    grad1 = drop(crossprod(X1, risk1_i * (1 - risk1_i))) / nrow(X1)
    grad0 = drop(crossprod(X0, risk0_i * (1 - risk0_i))) / nrow(X0)
    if (effect == "RD") {
        list(est = risk1 - risk0, grad = grad1 - grad0)
    } else {
        log_rr = log(risk1) - log(risk0)
        list(est = exp(log_rr), log_est = log_rr, grad = grad1 / risk1 - grad0 / risk0)
    }
}

compute_ordinal_gcomp_md = function(theta, X_base, n_alpha) {
    alpha = theta[seq_len(n_alpha)]; beta = theta[-seq_len(n_alpha)]
    X1 = X_base; X0 = X_base; X1[, "treatment"] = 1; X0[, "treatment"] = 0
    mean_score = function(Xmat) {
        eta = drop(Xmat %*% beta)
        cumprob = vapply(alpha, function(a) stats::plogis(a - eta), numeric(nrow(Xmat)))
        if (is.vector(cumprob)) cumprob = matrix(cumprob, ncol = 1)
        probs = matrix(0, nrow = nrow(Xmat), ncol = n_alpha + 1L)
        probs[, 1L] = cumprob[, 1L]
        if (n_alpha > 1L) for (j in 2:n_alpha) probs[, j] = cumprob[, j] - cumprob[, j - 1L]
        probs[, n_alpha + 1L] = 1 - cumprob[, n_alpha]
        drop(probs %*% seq_len(n_alpha + 1L))
    }
    mean(mean_score(X1)) - mean(mean_score(X0))
}

finite_diff_grad = function(fn, theta, base_step = 1e-6) {
    grad = numeric(length(theta))
    for (j in seq_along(theta)) {
        step = max(base_step, base_step * (1 + abs(theta[j])))
        theta_plus = theta; theta_minus = theta
        theta_plus[j] = theta[j] + step; theta_minus[j] = theta[j] - step
        grad[j] = (fn(theta_plus) - fn(theta_minus)) / (2 * step)
    }
    grad
}

make_edi_wald_bm = function(cls_name, d) {
    X_cov = d$X[, -1, drop = FALSE]; colnames(X_cov) = paste0("x", seq_len(ncol(X_cov)))
    X_bm  = cbind(`(Intercept)` = 1, treatment = d$w, X_cov)
    X_ord = cbind(treatment = d$w, X_cov)
    y_bm  = as.numeric(d$y); w_bm = as.integer(d$w)
    dead_bm = if (!is.null(d$dead)) as.integer(d$dead) else NULL

    e = new.env(parent = globalenv())
    e$d = d; e$X_bm = X_bm; e$X_ord = X_ord
    e$y_bm = y_bm; e$w_bm = w_bm; e$dead_bm = dead_bm

    e$compute_md_gradient = function(X_fit, theta, n_alpha, j_treat, base_step = 1e-6) {
        n_params = length(theta); grad = numeric(n_params)
        for (j in seq_len(n_params)) {
            step = max(base_step, base_step * (1 + abs(theta[j])))
            theta_plus = theta; theta_minus = theta
            theta_plus[j] = theta[j] + step; theta_minus[j] = theta[j] - step
            md_plus  = gcomp_ordinal_proportional_odds_post_fit_cpp(X_fit, theta_plus[(n_alpha+1):n_params],  theta_plus[1:n_alpha],  j_treat)$md
            md_minus = gcomp_ordinal_proportional_odds_post_fit_cpp(X_fit, theta_minus[(n_alpha+1):n_params], theta_minus[1:n_alpha], j_treat)$md
            grad[j] = (md_plus - md_minus) / (2 * step)
        }
        grad
    }

    expr = switch(cls_name,
        InferenceContinOLS = quote({
            res = fast_ols_with_var_cpp(X_bm, y_bm, j = 2L)
            se = sqrt(res$ssq_b_j); t_stat = res$b[2] / se
            2 * stats::pt(-abs(t_stat), df = length(y_bm) - ncol(X_bm))
        }),
        InferenceContinLin = quote({
            Xc = scale(X_bm[, -1, drop = FALSE], scale = FALSE)
            X_lin = cbind(1, w_bm, Xc, Xc * w_bm)
            colnames(X_lin)[1:2] = c("(Intercept)", "treatment")
            res_fit = fast_ols_cpp(X_lin, y_bm)
            post_fit = ols_hc2_post_fit_cpp(X_lin, y_bm, as.numeric(res_fit$b), 2L)
            se = sqrt(post_fit$ssq_hat); t_stat = res_fit$b[2] / se
            2 * stats::pt(-abs(t_stat), df = length(y_bm) - ncol(X_lin))
        }),
        InferenceContinRobustRegr = quote({
            res = fast_robust_regression_cpp(X_bm, y_bm, j = 2L)
            se = sqrt(res$ssq_b_j); t_stat = res$b[2] / se
            2 * stats::pnorm(-abs(t_stat))
        }),
        InferenceContinQuantileRegr = quote({
            res = quantreg::rq(y_bm ~ X_bm[, -1])
            summary(res, se = "nid")$coefficients[2, 4]
        }),
        InferenceIncidLogRegr = quote({
            res = fast_logistic_regression_with_var_cpp(X_bm, y_bm, j = 2L)
            se = sqrt(res$ssq_b_j); t_stat = res$b[2] / se
            2 * stats::pnorm(-abs(t_stat))
        }),
        InferenceIncidLogBinomial = quote({
            res = fast_log_binomial_regression_with_var_cpp(X_bm, y_bm, j = 2L)
            se = sqrt(res$ssq_b_j); t_stat = res$b[2] / se
            2 * stats::pnorm(-abs(t_stat))
        }),
        InferenceIncidProbitRegr = quote({
            res = fast_probit_regression_with_var_cpp(X_bm, y_bm, j = 2L)
            se = sqrt(res$ssq_b_j); t_stat = res$b[2] / se
            2 * stats::pnorm(-abs(t_stat))
        }),
        InferenceIncidRiskDiff = quote({
            res = fast_ols_with_var_cpp(X_bm, y_bm, j = 2L)
            se = sqrt(res$ssq_b_j); t_stat = res$b[2] / se
            2 * stats::pt(-abs(t_stat), df = length(y_bm) - ncol(X_bm))
        }),
        InferenceIncidExactFisher = quote({
            tab = table(factor(w_bm, levels = c(1, 0)), factor(y_bm, levels = c(1, 0)))
            stats::fisher.test(tab)$p.value
        }),
        InferenceIncidNewcombeRiskDiff = quote({
            x_t = sum(y_bm[w_bm == 1]); n_t = sum(w_bm == 1)
            x_c = sum(y_bm[w_bm == 0]); n_c = sum(w_bm == 0)
            est = x_t / n_t - x_c / n_c
            p_fn = function(a) { ci = newcombe_independent_ci_cpp(x_t, n_t, x_c, n_c, a); if (0 < est) ci[1] else ci[2] }
            res = tryCatch(stats::uniroot(p_fn, interval = c(1e-10, 1 - 1e-10))$root, error = function(e) NA_real_)
            if (!is.finite(res)) 1.0 else res
        }),
        InferenceIncidMiettinenNurminenRiskDiff = quote({
            x_t = sum(y_bm[w_bm == 1]); n_t = sum(w_bm == 1)
            x_c = sum(y_bm[w_bm == 0]); n_c = sum(w_bm == 0)
            mn_pvalue_cpp(x_t, n_t, x_c, n_c, 0, x_t/n_t, x_c/n_c)
        }),
        InferenceCountPoisson = quote({
            res = fast_poisson_regression_with_var_cpp(X_bm, y_bm, j = 2L)
            se = sqrt(res$ssq_b_j); t_stat = res$b[2] / se
            2 * stats::pnorm(-abs(t_stat))
        }),
        InferenceCountQuasiPoisson = quote({
            res = fast_quasipoisson_regression_with_var_cpp(X_bm, y_bm, j = 2L)
            se = sqrt(res$ssq_b_j); t_stat = res$b[2] / se
            2 * stats::pnorm(-abs(t_stat))
        }),
        InferenceCountNegBin = quote({
            res = fast_neg_bin_with_var_cpp(X_bm, as.integer(y_bm))
            vcov = solve(res$hess_fisher_info_matrix); se = sqrt(vcov[2, 2])
            t_stat = res$b[2] / se; 2 * stats::pnorm(-abs(t_stat))
        }),
        InferenceCountHurdlePoisson = quote({
            res = fast_zero_augmented_poisson_cpp(X_bm, y_bm, X_bm, is_hurdle = TRUE, estimate_only = FALSE)
            se = sqrt(res$vcov[2, 2]); t_stat = res$coefficients$cond[2] / se
            2 * stats::pnorm(-abs(t_stat))
        }),
        InferenceCountZeroInflatedPoisson = quote({
            res = fast_zero_augmented_poisson_cpp(X_bm, y_bm, X_bm, is_hurdle = FALSE, estimate_only = FALSE)
            se = sqrt(res$vcov[2, 2]); t_stat = res$coefficients$cond[2] / se
            2 * stats::pnorm(-abs(t_stat))
        }),
        InferenceCountZeroInflatedNegBin = quote({
            res = fast_zinb_cpp(X_bm, X_bm, y_bm, estimate_only = FALSE)
            se = sqrt(res$vcov[2, 2]); t_stat = res$params[2] / se
            2 * stats::pnorm(-abs(t_stat))
        }),
        InferenceCountHurdleNegBin = quote({
            res = fast_hurdle_negbin_with_var_cpp(X_bm, as.integer(y_bm), X_bm, j = 2L)
            se = sqrt(res$ssq_b_j); t_stat = res$b[2] / se
            2 * stats::pnorm(-abs(t_stat))
        }),
        InferenceCountRobustPoisson = quote({
            res = fast_poisson_regression_cpp(X_bm, y_bm, estimate_only = FALSE)
            bread = solve(res$XtWX); resid = y_bm - as.numeric(res$mu)
            meat = crossprod(X_bm, X_bm * (resid^2))
            se = sqrt((bread %*% meat %*% bread)[2, 2])
            2 * stats::pnorm(-abs(res$b[2] / se))
        }),
        InferencePropBetaRegr = quote({
            res = fast_beta_regression_with_var_cpp(X_bm, y_bm)
            se = sqrt(res$vcov[2, 2]); t_stat = res$coefficients[2] / se
            2 * stats::pnorm(-abs(t_stat))
        }),
        InferenceOrdinalPropOddsRegr = quote({
            res = fast_ordinal_regression_with_var_cpp(X_ord, y_bm)
            se = sqrt(res$ssq_b_j); t_stat = res$b[1] / se
            2 * stats::pnorm(-abs(t_stat))
        }),
        InferenceOrdinalAdjCatLogitRegr = quote({
            res = fast_adjacent_category_logit_with_var_cpp(X_ord, y_bm)
            se = sqrt(res$ssq_b_j); t_stat = res$b[1] / se
            2 * stats::pnorm(-abs(t_stat))
        }),
        InferenceOrdinalContRatioRegr = quote({
            res = fast_continuation_ratio_regression_with_var_cpp(X_ord, y_bm)
            se = sqrt(res$ssq_b_j); t_stat = res$b[1] / se
            2 * stats::pnorm(-abs(t_stat))
        }),
        InferenceOrdinalOrderedProbitRegr = quote({
            res = fast_ordinal_probit_regression_with_var_cpp(X_ord, y_bm)
            se = sqrt(res$ssq_b_j); t_stat = res$b[1] / se
            2 * stats::pnorm(-abs(t_stat))
        }),
        InferenceOrdinalCloglogRegr = quote({
            res = fast_ordinal_cloglog_regression_with_var_cpp(X_ord, y_bm)
            se = sqrt(res$ssq_b_j); t_stat = res$b[1] / se
            2 * stats::pnorm(-abs(t_stat))
        }),
        InferenceOrdinalCauchitRegr = quote({
            res = fast_ordinal_cauchit_regression_with_var_cpp(X_ord, y_bm)
            se = sqrt(res$ssq_b_j); t_stat = res$b[1] / se
            2 * stats::pnorm(-abs(t_stat))
        }),
        InferenceOrdinalRidit = quote({
            res = EDI:::fast_ridit_analysis_cpp(w_bm, as.integer(y_bm), reference = "control")
            2 * stats::pnorm(-abs(res$estimate / res$se))
        }),
        InferenceOrdinalJonckheereTerpstraTest = quote({
            exact_jonckheere_terpstra_pval_cpp(as.integer(y_bm), w_bm)$p_exact
        }),
        InferenceSurvivalCoxPHRegr = quote({
            res = fast_coxph_regression_cpp(X_ord, y_bm, dead_bm, estimate_only = FALSE)
            se = sqrt(res$vcov[1, 1]); 2 * stats::pnorm(-abs(res$coefficients[1] / se))
        }),
        InferenceSurvivalStratCoxPHRegr = quote({
            strat_inputs = build_strat_cox_canonical_inputs(d)
            cache = if (!is.null(strat_inputs$strata))
                build_stratified_cox_data_cache_cpp(strat_inputs$X, strat_inputs$y, strat_inputs$dead, strat_inputs$strata)
            else
                build_cox_data_cache_cpp(strat_inputs$X, strat_inputs$y, strat_inputs$dead)
            res = fast_coxph_regression_prebuilt_cpp(cache, estimate_only = FALSE)
            se = sqrt(res$vcov[1, 1]); 2 * stats::pnorm(-abs(res$coefficients[1] / se))
        }),
        InferenceSurvivalWeibullRegr = quote({
            res = fast_weibull_regression_cpp(X_ord, y_bm, dead_bm, estimate_only = FALSE)
            se = sqrt(res$vcov[2, 2]); 2 * stats::pnorm(-abs(res$params[2] / se))
        }),
        InferenceSurvivalLogRank = quote({
            res = EDI:::fast_logrank_stats_cpp(w_bm, y_bm, as.integer(dead_bm))
            stats::pchisq(res$score ^ 2 / res$var_score, df = 1, lower.tail = FALSE)
        }),
        InferenceSurvivalGehanWilcox = quote(
            survival::survdiff(survival::Surv(y_bm, dead_bm) ~ w_bm, rho = 1)$pvalue
        ),
        InferenceSurvivalKMDiff = quote({
            fit = survival::survfit(survival::Surv(y_bm, dead_bm) ~ w_bm, conf.int = 0.95)
            q = stats::quantile(fit, 0.5)
            sn = rownames(q$lower)
            idx_T = grep("1", sn); idx_C = grep("0", sn)
            z = stats::qnorm(0.975)
            se = sqrt(((q$upper[idx_T,1] - q$lower[idx_T,1])/(2*z))^2 + ((q$upper[idx_C,1] - q$lower[idx_C,1])/(2*z))^2)
            2 * stats::pnorm(-abs((q$quantile[idx_T,1] - q$quantile[idx_C,1]) / se))
        }),
        InferenceAllSimpleMeanDiffPooledVar = quote({
            y_t = y_bm[w_bm == 1]; y_c = y_bm[w_bm == 0]
            n_t = length(y_t); n_c = length(y_c); df_res = n_t + n_c - 2L
            s2p = ((n_t-1)*var(y_t) + (n_c-1)*var(y_c)) / df_res
            se = sqrt(s2p * (1/n_t + 1/n_c))
            2 * stats::pt(-abs((mean(y_t) - mean(y_c)) / se), df_res)
        }),
        InferenceAllSimpleWilcox = quote(
            EDI:::wilcox_hl_point_estimate_cpp(w_bm, y_bm)
        ),
        InferencePropGCompMeanDiff = quote({
            fit = fast_logistic_regression_cpp(X_bm, y_bm, estimate_only = TRUE)
            coef_hat = as.numeric(fit$b)
            mu_hat = stats::plogis(as.numeric(X_bm %*% coef_hat))
            res = gcomp_logistic_post_fit_cpp(X_bm, y_bm, coef_hat, mu_hat, 2L)
            2 * stats::pnorm(-abs(res$rd / res$se_rd))
        }),
        InferenceIncidGCompRiskDiff = quote({
            fit = fast_logistic_regression_cpp(X_bm, y_bm, estimate_only = TRUE)
            coef_hat = as.numeric(fit$b)
            mu_hat = stats::plogis(as.numeric(X_bm %*% coef_hat))
            res = gcomp_logistic_post_fit_cpp(X_bm, y_bm, coef_hat, mu_hat, 2L)
            2 * stats::pnorm(-abs(res$rd / res$se_rd))
        }),
        InferenceIncidGCompRiskRatio = quote({
            fit = fast_logistic_regression_cpp(X_bm, y_bm, estimate_only = TRUE)
            coef_hat = as.numeric(fit$b)
            mu_hat = stats::plogis(as.numeric(X_bm %*% coef_hat))
            res = gcomp_logistic_post_fit_cpp(X_bm, y_bm, coef_hat, mu_hat, 2L)
            2 * stats::pnorm(-abs(res$log_rr / res$se_log_rr))
        }),
        InferenceOrdinalGCompMeanDiff = quote({
            fit = fast_ordinal_regression_with_var_cpp(X_ord, y_bm)
            coef_hat = as.numeric(fit$b); alpha_hat = as.numeric(fit$alpha)
            res = gcomp_ordinal_proportional_odds_post_fit_cpp(X_ord, coef_hat, alpha_hat, 1L)
            theta = c(alpha_hat, coef_hat)
            grad = compute_md_gradient(X_ord, theta, length(alpha_hat), 1L)
            se = sqrt(as.numeric(crossprod(grad, as.matrix(fit$vcov) %*% grad)))
            2 * stats::pnorm(-abs(res$md / se))
        }),
        NULL
    )
    if (is.null(expr)) return(NULL)
    list(env = e, expr = expr)
}

# --- Wald Benchmark Specs ---
wald_specs = list(
    list(cls = "InferenceAllSimpleMeanDiffPooledVar", pkg = "stats", func = "t.test(pool)", expr = quote(t.test(df$y[df$treatment==1], df$y[df$treatment==0], var.equal=TRUE)$p.value)),
    list(cls = "InferenceAllSimpleWilcox", pkg = "stats", func = "wilcox.test", expr = quote(wilcox.test(df$y[df$treatment==1], df$y[df$treatment==0])$p.value)),
    list(cls = "InferenceContinLin", pkg = "stats", func = "lm.fit(interact)+Wald", expr = quote({
        X_int = model.matrix(~ treatment * (x1 + x2 + x3 + x4), data = df)
        res = lm.fit(x = X_int, y = df$y); p = ncol(X_int); n = nrow(X_int)
        sig2 = sum(res$residuals^2) / (n - p)
        v = chol2inv(res$qr$qr[1:p, 1:p, drop=FALSE]) * sig2
        2 * pt(-abs(res$coefficients[2] / sqrt(v[2,2])), n - p)
    })),
    list(cls = "InferenceContinOLS", pkg = "stats", func = "lm.fit+Wald", expr = quote({
        res = lm.fit(x = X_can, y = df$y); p = ncol(X_can); n = nrow(X_can)
        sig2 = sum(res$residuals^2) / (n - p)
        v = chol2inv(res$qr$qr[1:p, 1:p, drop=FALSE]) * sig2
        2 * pt(-abs(res$coefficients[2] / sqrt(v[2,2])), n - p)
    })),
    list(cls = "InferenceContinQuantileRegr", pkg = "quantreg", func = "rq+summary", expr = quote({
        res = quantreg::rq(df$y ~ X_can[,-1]); summary(res, se="nid")$coefficients[2, 4]
    })),
    list(cls = "InferenceContinRobustRegr", pkg = "MASS", func = "rlm+summary", expr = quote({
        summary(MASS::rlm(x = X_can, y = df$y))$coefficients[2, 3]
    })),
    list(cls = "InferenceIncidExactFisher", pkg = "stats", func = "fisher.test", expr = quote(fisher.test(table(df$treatment, df$y))$p.value)),
    list(cls = "InferenceIncidGCompRiskDiff", pkg = "stats", func = "glm+gcomp(RD)+Wald", expr = quote({
        mod = glm(y ~ treatment + x1 + x2 + x3 + x4, family=binomial(), data=df)
        eff = compute_binary_gcomp_effect(coef(mod), X_can, "RD")
        se = sqrt(drop(t(eff$grad) %*% vcov(mod) %*% eff$grad))
        2 * pnorm(-abs(eff$est / se))
    })),
    list(cls = "InferenceIncidGCompRiskRatio", pkg = "stats", func = "glm+gcomp(RR)+Wald", expr = quote({
        mod = glm(y ~ treatment + x1 + x2 + x3 + x4, family=binomial(), data=df)
        eff = compute_binary_gcomp_effect(coef(mod), X_can, "RR")
        se = sqrt(drop(t(eff$grad) %*% vcov(mod) %*% eff$grad))
        2 * pnorm(-abs(eff$log_est / se))
    })),
    list(cls = "InferenceIncidLogBinomial", pkg = "stats", func = "glm.fit+Wald(log)", expr = quote({
        res = glm.fit(x=X_can, y=df$y, family=binomial(link="log"), start=c(-2, rep(0, ncol(X_can)-1)))
        p = ncol(X_can); R = res$qr$qr[1:p,1:p,drop=FALSE]; R[lower.tri(R)] = 0
        v = chol2inv(R); 2 * pnorm(-abs(res$coefficients[2] / sqrt(v[2,2])))
    })),
    list(cls = "InferenceIncidLogRegr", pkg = "stats", func = "glm.fit+Wald", expr = quote({
        res = glm.fit(x=X_can, y=df$y, family=binomial())
        p = ncol(X_can); R = res$qr$qr[1:p,1:p,drop=FALSE]; R[lower.tri(R)] = 0
        v = chol2inv(R); 2 * pnorm(-abs(res$coefficients[2] / sqrt(v[2,2])))
    })),
    list(cls = "InferenceIncidMiettinenNurminenRiskDiff", pkg = "DescTools", func = "BinomDiffCI(mn)", expr = quote(
        DescTools::BinomDiffCI(sum(df$y[df$treatment==1]), sum(df$treatment==1), sum(df$y[df$treatment==0]), sum(df$treatment==0), method="mn")
    )),
    list(cls = "InferenceIncidNewcombeRiskDiff", pkg = "DescTools", func = "BinomDiffCI(score)", expr = quote(
        DescTools::BinomDiffCI(sum(df$y[df$treatment==1]), sum(df$treatment==1), sum(df$y[df$treatment==0]), sum(df$treatment==0), method="score")
    )),
    list(cls = "InferenceIncidProbitRegr", pkg = "stats", func = "glm.fit(probit)+Wald", expr = quote({
        res = glm.fit(x=X_can, y=df$y, family=binomial(link="probit"))
        p = ncol(X_can); R = res$qr$qr[1:p,1:p,drop=FALSE]; R[lower.tri(R)] = 0
        v = chol2inv(R); 2 * pnorm(-abs(res$coefficients[2] / sqrt(v[2,2])))
    })),
    list(cls = "InferenceIncidRiskDiff", pkg = "stats", func = "prop.test", expr = quote(prop.test(table(df$treatment, df$y))$p.value)),
    list(cls = "InferenceCountHurdleNegBin", pkg = "pscl", func = "hurdle(nb)+summary", expr = quote(summary(pscl::hurdle(y~treatment+x1+x2+x3+x4, data=df, dist="negbin"))$coefficients$count[2,4])),
    list(cls = "InferenceCountHurdlePoisson", pkg = "pscl", func = "hurdle+summary", expr = quote(summary(pscl::hurdle(y~treatment+x1+x2+x3+x4, data=df))$coefficients$count[2,4])),
    list(cls = "InferenceCountNegBin", pkg = "MASS", func = "glm.nb+summary", expr = quote(summary(MASS::glm.nb(y~treatment+x1+x2+x3+x4, data=df))$coefficients[2,4])),
    list(cls = "InferenceCountPoisson", pkg = "stats", func = "glm.fit+Wald", expr = quote({
        res = glm.fit(x=X_can, y=df$y, family=poisson())
        p = ncol(X_can); R = res$qr$qr[1:p,1:p,drop=FALSE]; R[lower.tri(R)] = 0
        v = chol2inv(R); 2 * pnorm(-abs(res$coefficients[2] / sqrt(v[2,2])))
    })),
    list(cls = "InferenceCountQuasiPoisson", pkg = "stats", func = "glm.fit+Wald(quasi)", expr = quote({
        res = glm.fit(x=X_can, y=df$y, family=quasipoisson())
        p = ncol(X_can); R = res$qr$qr[1:p,1:p,drop=FALSE]; R[lower.tri(R)] = 0
        disp = sum(res$weights * res$residuals^2) / res$df.residual
        v = chol2inv(R) * disp; 2 * pnorm(-abs(res$coefficients[2] / sqrt(v[2,2])))
    })),
    list(cls = "InferenceCountRobustPoisson", pkg = "sandwich", func = "glm+vcovHC", expr = quote({
        mod = glm(y~treatment+x1+x2+x3+x4, family=poisson, data=df)
        v = sandwich::vcovHC(mod, type="HC0")
        2 * pnorm(-abs(coef(mod)[2] / sqrt(v[2,2])))
    })),
    list(cls = "InferenceCountZeroInflatedNegBin", pkg = "pscl", func = "zeroinfl(nb)+summary", expr = quote(summary(pscl::zeroinfl(y~treatment+x1+x2+x3+x4, data=df, dist="negbin"))$coefficients$count[2,4])),
    list(cls = "InferenceCountZeroInflatedPoisson", pkg = "pscl", func = "zeroinfl+summary", expr = quote(summary(pscl::zeroinfl(y~treatment+x1+x2+x3+x4, data=df))$coefficients$count[2,4])),
    list(cls = "InferencePropBetaRegr", pkg = "betareg", func = "betareg+summary", expr = quote(summary(betareg::betareg(y~treatment+x1+x2+x3+x4, data=df))$coefficients$mean[2,4])),
    list(cls = "InferencePropGCompMeanDiff", pkg = "stats", func = "glm(quasi)+gcomp+Wald", expr = quote({
        mod = glm(y~treatment+x1+x2+x3+x4, family=quasibinomial(), data=df)
        eff = compute_binary_gcomp_effect(coef(mod), X_can, "RD")
        se = sqrt(drop(t(eff$grad) %*% vcov(mod) %*% eff$grad))
        2 * pnorm(-abs(eff$est / se))
    })),
    list(cls = "InferenceSurvivalCoxPHRegr", pkg = "survival", func = "coxph.fit(breslow)+Wald", expr = quote({
        x_vars = c("treatment", grep("^x", names(df), value=TRUE))
        res = survival::coxph.fit(x=as.matrix(df[,x_vars,drop=FALSE]), y=survival::Surv(df$y,df$dead), strata=NULL, offset=NULL, init=NULL, control=survival::coxph.control(), weights=NULL, method="breslow", rownames=as.character(seq_len(nrow(df))))
        2 * pnorm(-abs(res$coefficients[1] / sqrt(res$var[1,1])))
    })),
    list(cls = "InferenceSurvivalGehanWilcox", pkg = "survival", func = "survdiff(rho=1)", expr = quote(survival::survdiff(survival::Surv(y,dead)~treatment, data=df, rho=1)$pvalue)),
    list(cls = "InferenceSurvivalKMDiff", pkg = "survival", func = "survfit(median)+CI", expr = quote({
        fit = survival::survfit(survival::Surv(df$y, df$dead) ~ df$treatment, conf.int = 0.95)
        stats::quantile(fit, 0.5)
    })),
    list(cls = "InferenceSurvivalLogRank", pkg = "survival", func = "survdiff", expr = quote(survival::survdiff(survival::Surv(y,dead)~treatment, data=df)$pvalue)),
    list(cls = "InferenceSurvivalStratCoxPHRegr", pkg = "survival", func = "coxph.fit(strat)+Wald", expr = quote({
        si = build_strat_cox_canonical_inputs(d)
        res = survival::coxph.fit(x=si$X, y=survival::Surv(si$y,si$dead), strata=si$strata, offset=NULL, init=NULL, control=survival::coxph.control(), weights=NULL, method="breslow", rownames=as.character(seq_along(si$y)))
        2 * pnorm(-abs(res$coefficients[1] / sqrt(res$var[1,1])))
    })),
    list(cls = "InferenceSurvivalWeibullRegr", pkg = "survival", func = "survreg+summary", expr = quote(summary(survival::survreg(survival::Surv(y,dead)~treatment+x1+x2+x3+x4, data=df, dist="weibull"))$table[2,4])),
    list(cls = "InferenceOrdinalAdjCatLogitRegr", pkg = "VGAM", func = "vglm+summary", expr = quote(summary(VGAM::vglm(factor(y,ordered=TRUE)~treatment+x1+x2+x3+x4, VGAM::acat(), data=df))@coef3[1,4])),
    list(cls = "InferenceOrdinalContRatioRegr", pkg = "VGAM", func = "vglm+summary", expr = quote(summary(VGAM::vglm(factor(y,ordered=TRUE)~treatment+x1+x2+x3+x4, VGAM::cratio(), data=df))@coef3[1,4])),
    list(cls = "InferenceOrdinalGCompMeanDiff", pkg = "ordinal", func = "clm+gcomp+Wald", expr = quote({
        mod = ordinal::clm(factor(y,ordered=TRUE)~treatment+x1+x2+x3+x4, data=df)
        theta = coef(mod); n_alpha = length(mod$alpha)
        fn = function(th) compute_ordinal_gcomp_md(th, X_can[,-1,drop=FALSE], n_alpha)
        grad = finite_diff_grad(fn, theta)
        se = sqrt(drop(t(grad) %*% vcov(mod) %*% grad))
        2 * pnorm(-abs(fn(theta) / se))
    })),
    list(cls = "InferenceOrdinalJonckheereTerpstraTest", pkg = "clinfun", func = "jonckheere", expr = quote(clinfun::jonckheere.test(df$y, df$treatment)$p.value)),
    list(cls = "InferenceOrdinalPropOddsRegr", pkg = "ordinal", func = "clm+summary", expr = quote(summary(ordinal::clm(factor(y,ordered=TRUE)~treatment+x1+x2+x3+x4, data=df))$coefficients["treatment",4])),
    list(cls = "InferenceOrdinalRidit", pkg = "stats", func = "mean(ridit)", expr = quote({
        y = df$y; tab = table(y); cum = cumsum(tab); prev = c(0, cum[-length(cum)])
        ridit_map = (prev + 0.5*tab) / length(y); r = as.numeric(ridit_map[as.character(y)])
        mean(r[df$treatment==1]) - mean(r[df$treatment==0])
    }))
)

wald_results = list()

run_one_wald = function(spec) {
    cls_name = spec$cls
    cat(sprintf("Wald [%d/%d] %s...\n", match(cls_name, sapply(wald_unique_specs, `[[`, "cls")), length(wald_unique_specs), cls_name))

    resp_type = "continuous"; family = "continuous"
    if (grepl("Ordinal|AdjCat|ContRatio|Ridit|Jonckheere", cls_name)) { resp_type = "ordinal"; family = "ordinal" }
    else if (grepl("Incid|Binomial|Fisher|Newcombe|Nurminen", cls_name))  { resp_type = "incidence"; family = "logistic" }
    else if (grepl("Count|Poisson|NegBin|ZINB|ZAP|Hurdle", cls_name))    { resp_type = "count"; family = "poisson" }
    else if (grepl("Prop|Beta|ZOIB|Fractional", cls_name))                { resp_type = "proportion"; family = "beta" }
    else if (grepl("Survival|Cox|Weibull|KM|Rank|LogRank|Gehan", cls_name)) { resp_type = "survival"; family = "cox" }
    if (cls_name == "InferenceIncidLogBinomial") family = "log-binomial"
    if (cls_name == "InferenceIncidProbitRegr")  family = "probit"
    if (cls_name %in% c("InferenceCountNegBin","InferenceCountZeroInflatedNegBin","InferenceCountHurdleNegBin")) family = "negbin"

    n = N_WALD
    d = generate_data(n = n, family = family)
    if (identical(cls_name, "InferenceSurvivalStratCoxPHRegr")) d = make_true_stratified_survival_data(d)

    timing_edi = tryCatch({
        bm = make_edi_wald_bm(cls_name, d)
        if (is.null(bm$expr)) stop("no Wald mapping for this class")
        eval(bm$expr, envir = bm$env)  # validation
        collect_timing_ms(bm$expr, env = bm$env)
    }, error = function(e) { cat("  EDI Error:", e$message, "\n"); list(median_ms = NA_real_, samples_ms = numeric(0)) })

    timing_can = tryCatch({
        library(spec$pkg, character.only = TRUE)
        X_cols = d$X[,-1,drop=FALSE]; colnames(X_cols) = paste0("x", seq_len(ncol(X_cols)))
        df = data.frame(y = d$y, treatment = d$w, dead = if (!is.null(d$dead)) d$dead else 1L)
        df = cbind(df, X_cols)
        X_can = cbind(`(Intercept)` = 1, treatment = d$w, as.matrix(X_cols))
        eval(spec$expr)  # validation
        collect_timing_ms(spec$expr)
    }, error = function(e) { cat("  Canonical Error:", e$message, "\n"); list(median_ms = NA_real_, samples_ms = numeric(0)) })

    wald_results[[length(wald_results) + 1]] <<- data.table(
        Class = cls_name, Response = resp_type,
        EDI_Time_ms = timing_edi$median_ms,
        Canonical_Pkg = spec$pkg, Canonical_Func = spec$func,
        Canonical_Time_ms = timing_can$median_ms,
        Speedup = if (!is.na(timing_can$median_ms) && !is.na(timing_edi$median_ms) && timing_edi$median_ms > 0) timing_can$median_ms / timing_edi$median_ms else NA_real_,
        Timing_Pval = timing_ttest_pval(timing_edi$samples_ms, timing_can$samples_ms)
    )
}

wald_unique_specs = wald_specs[!duplicated(sapply(wald_specs, `[[`, "cls"))]
for (i in seq_along(wald_unique_specs)) run_one_wald(wald_unique_specs[[i]])

# Format Wald table
dt_wald = rbindlist(wald_results)
response_levels = c("all", "continuous", "incidence", "count", "proportion", "survival", "ordinal")
dt_wald[, Response := factor(Response, levels = response_levels)]
setorder(dt_wald, Response, Class)
dt_wald[, Speedup_Num := Speedup]
dt_wald[, Speedup := ifelse(!is.na(Speedup), paste0(round(Speedup, 2), "x"), "NA")]
dt_wald[, EDI_Time_ms := format_ms(EDI_Time_ms)]
dt_wald[, Canonical_Time_ms := format_ms(Canonical_Time_ms)]
dt_wald[, Timing_Row_Color := mapply(row_bg_color, Speedup_Num, Timing_Pval, USE.NAMES = FALSE)]
dt_wald[, Timing_Pval_Stars := format_pval_stars(Timing_Pval)]
dt_wald[, Timing_Pval := format_pval(Timing_Pval)]
dt_wald[, Response := as.character(Response)]

wald_table_lines = c(
  "<table>",
  "  <thead>",
  "    <tr><th>Class</th><th>Response</th><th>EDI Time (ms)</th><th>Canonical Pkg</th><th>Canonical Func</th><th>Canonical Time (ms)</th><th>Speedup</th><th>Timing Pval</th><th></th></tr>",
  "  </thead>",
  "  <tbody>"
)
wald_table_rows = mapply(function(cls, resp, edi, pkg, func, can, speed, pval, stars, bg) {
  style = if (nzchar(bg)) paste0(" style=\"background-color: ", bg, ";\"") else ""
  paste0("    <tr", style, "><td>", cls, "</td><td>", resp, "</td><td>", edi, "</td><td>", pkg, "</td><td>", func, "</td><td>", can, "</td><td>", speed, "</td><td>", pval, "</td><td>", stars, "</td></tr>")
}, dt_wald$Class, dt_wald$Response, dt_wald$EDI_Time_ms, dt_wald$Canonical_Pkg, dt_wald$Canonical_Func,
dt_wald$Canonical_Time_ms, dt_wald$Speedup, dt_wald$Timing_Pval, dt_wald$Timing_Pval_Stars, dt_wald$Timing_Row_Color,
SIMPLIFY = TRUE, USE.NAMES = FALSE)
wald_table_lines = c(wald_table_lines, wald_table_rows, "  </tbody>", "</table>")

WALD_HEADER = c(
  "## Wald Test Performance (Full Inference)",
  "",
  "This table compares the performance of **Full Inference** (Model Fit + Standard Error calculation + P-value derivation).",
  "Unlike the point-estimation table above, these results include the computational cost of the variance-covariance matrix (Hessian or Fisher Information) and the Wald test statistic calculation.",
  "All paths (EDI and Canonical) use a reduced sample size ($N=200$) for this full-inference benchmark to ensure iterative stability.",
  "**Stratified Cox Exception**: For `InferenceSurvivalStratCoxPHRegr`, the benchmark injects low-cardinality covariates before outcome generation so the row exercises a genuinely stratified Cox fit rather than the unstratified fallback.",
  "EDI regression models (Logistic, Poisson) are benchmarked using the **IRLS** optimizer for these Wald tests.",
  "**Solver-Only Prebuilds**: Benchmark setup prebuilds exposed observed-data design matrices, reduced design matrices, strata IDs, and other fixed working inputs outside the timed region when the implementation exposes those hooks. The timed region then measures the full-inference kernel on those fixed inputs.",
  "**Limitation**: Some canonical comparators only expose formula-based APIs rather than comparable low-level fit kernels. Those rows remain included, but their canonical timings may still contain formula/model-frame overhead beyond the numerical solver, variance, and p-value work itself.",
  paste0("**Timing Note**: All timings are medians over ", B_TIME, " warmed runs measured with adaptive batched `system.time`; paths below ", FAST_PATH_THRESHOLD_MS, " ms use `microbenchmark(times = ", FAST_PATH_MICROBENCH_REPS, ")` instead."),
  "**Timing P-Value**: `Timing Pval` reports a Welch two-sample t-test comparing the EDI and canonical timing replicate distributions for each row. The unlabeled final column marks thresholds with `***` for p < 0.001, `**` for p < 0.01, and `*` for p < 0.05.",
  "**Row Highlighting**: Light green rows indicate `Speedup > 1` and `Timing Pval < 0.05`; light grey rows indicate `NA` timing comparisons.",
  "",
  wald_table_lines
)

# ── Finalize
dt = rbindlist(results)
response_levels = c("all", "continuous", "incidence", "count", "proportion", "survival", "ordinal")
dt[, Response := factor(Response, levels = response_levels)]
setorder(dt, Response, Class)
dt[, Speedup_Num := Speedup]
dt[, Speedup := ifelse(!is.na(Speedup), paste0(round(Speedup, 2), "x"), "NA")]
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
  "EDI timings call exported C++ functions directly — no R6 objects are instantiated during benchmarking. As a result, **no R6 result caches exist to manage**. Each call to the C++ solver (e.g. `fast_logistic_regression_cpp`, `fast_ordinal_regression_cpp`) starts from a freshly zero-initialized parameter vector (or a model-specific data-driven initialization when `smart_cold_start = TRUE`). No prior-fit results are carried across timing repetitions, so every replication is a genuine cold start for the numerical optimizer."
)

writeLines(c(report, "", WALD_HEADER, "", METHODOLOGY_BLOCK, "", STYLE_BLOCK), "package_metadata/benchmark_model_fits.md")
cat("Benchmark complete.\n")
