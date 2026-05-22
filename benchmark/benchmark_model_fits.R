
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

prepare_solver_only_edi = function(inf_obj, cls_name) {
    priv = inf_obj$.__enclos_env__$private
    replace_binding = function(name, value) {
        if (!exists(name, envir = priv, inherits = FALSE)) return(invisible(FALSE))
        was_locked = bindingIsLocked(name, priv)
        if (was_locked) unlockBinding(name, priv)
        assign(name, value, envir = priv)
        if (was_locked) lockBinding(name, priv)
        invisible(TRUE)
    }

    if (exists("build_design_matrix", envir = priv, inherits = FALSE)) {
        X_built = tryCatch(priv$build_design_matrix(), error = function(e) NULL)
        if (!is.null(X_built)) {
            replace_binding("build_design_matrix", function() X_built)
        }
    }
    if (exists("create_design_matrix", envir = priv, inherits = FALSE)) {
        X_created = tryCatch(priv$create_design_matrix(), error = function(e) NULL)
        if (!is.null(X_created)) {
            replace_binding("create_design_matrix", function() X_created)
        }
    }

    if (identical(cls_name, "InferenceContinQuantileRegr") &&
        exists("build_design_matrix", envir = priv, inherits = FALSE) &&
        exists("reduce_design_matrix_for_quantile", envir = priv, inherits = FALSE)) {
        X_full = tryCatch(priv$build_design_matrix(), error = function(e) NULL)
        reduced = tryCatch(priv$reduce_design_matrix_for_quantile(X_full, reuse_factorizations = FALSE), error = function(e) NULL)
        if (!is.null(reduced)) {
            replace_binding("reduce_design_matrix_for_quantile", function(X_full, reuse_factorizations = FALSE) reduced)
        }
    }

    if (identical(cls_name, "InferenceContinLin") &&
        exists("build_lin_design_matrix", envir = priv, inherits = FALSE) &&
        exists("reduce_design_matrix_preserving_treatment", envir = priv, inherits = FALSE)) {
        X_lin = tryCatch(priv$build_lin_design_matrix(), error = function(e) NULL)
        reduced = tryCatch(priv$reduce_design_matrix_preserving_treatment(X_lin), error = function(e) NULL)
        if (!is.null(X_lin)) {
            replace_binding("build_lin_design_matrix", function() X_lin)
        }
        if (!is.null(reduced)) {
            replace_binding("reduce_design_matrix_preserving_treatment", function(X_full) reduced)
        }
    }

    if (identical(cls_name, "InferenceSurvivalStratCoxPHRegr") &&
        exists("compute_strata_info", envir = priv, inherits = FALSE)) {
        X_full = tryCatch(priv$X, error = function(e) NULL)
        strata_info = tryCatch(priv$compute_strata_info(X_full), error = function(e) NULL)
        if (!is.null(strata_info)) {
            replace_binding("compute_strata_info", function(X_full) strata_info)
            informative_rows = tryCatch(priv$get_informative_rows(strata_info$strata_id), error = function(e) integer(0))
            replace_binding("get_informative_rows", function(strata_id) informative_rows)
        }
    }

    invisible(inf_obj)
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
    list(cls = "InferenceSurvivalGehanWilcox", pkg = "survival", func = "coxph(null)+KM weighted residual mean diff", expr = quote({
        surv_obj = survival::Surv(df$y, df$dead)
        cox_null = survival::coxph(surv_obj ~ 1)
        M = as.numeric(stats::residuals(cox_null, type = "martingale"))
        km_all = survival::survfit(surv_obj ~ 1)
        idx = findInterval(df$y, km_all$time, left.open = TRUE)
        peto_weights = c(1.0, km_all$surv)[idx + 1L]
        M_w = peto_weights * M
        mean(M_w[df$treatment == 1]) - mean(M_w[df$treatment == 0])
    })),
    list(cls = "InferenceAllSimpleMeanDiffPooledVar", pkg = "stats", func = "t.test(pool)", expr = quote(t.test(df$y[df$treatment==1], df$y[df$treatment==0], var.equal = TRUE))),
    list(cls = "InferenceAllSimpleWilcox", pkg = "stats", func = "HL median pairwise diff", expr = quote({
        y_t = df$y[df$treatment == 1]
        y_c = df$y[df$treatment == 0]
        stats::median(as.numeric(outer(y_t, y_c, "-")))
    }), scale = 0.5),
    list(cls = "InferenceIncidExactFisher", pkg = "stats", func = "fisher.test", expr = quote(fisher.test(table(df$treatment, df$y))), scale = 0.1),
    list(cls = "InferenceIncidMiettinenNurminenRiskDiff", pkg = "DescTools", func = "BinomDiffCI(mn)", expr = quote(DescTools::BinomDiffCI(sum(df$y[df$treatment==1]), sum(df$treatment==1), sum(df$y[df$treatment==0]), sum(df$treatment==0), method="mn"))),
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
    list(cls = "InferenceContinLin", pkg = "stats", func = "lm.fit(interact)", expr = quote({
        X_int = model.matrix(~ treatment * (x1 + x2 + x3 + x4), data = df)
        lm.fit(x = X_int, y = df$y)
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
    list(cls = "InferenceCountRobustPoisson", pkg = "stats", func = "glm.fit", expr = quote(glm.fit(x = X_can, y = df$y, family = poisson()))),
    list(cls = "InferenceOrdinalRidit", pkg = "stats", func = "mean(ridit; ref=control)", expr = quote({
        y = df$y
        tab = table(y[df$treatment == 0])
        cum = cumsum(tab)
        prev = c(0, cum[-length(cum)])
        ridit_map = (prev + 0.5 * tab) / sum(tab)
        r = as.numeric(ridit_map[as.character(y)])
        mean(r[df$treatment == 1]) - 0.5
    })),
    list(cls = "InferenceIncidNewcombeRiskDiff", pkg = "DescTools", func = "BinomDiffCI(score)", expr = quote(DescTools::BinomDiffCI(sum(df$y[df$treatment==1]), sum(df$treatment==1), sum(df$y[df$treatment==0]), sum(df$treatment==0), method="score")))
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
    
    # Timing EDI
    timing_edi = tryCatch({
        des = DesignFixediBCRD$new(n = n, response_type = resp_type)
        des$add_all_subjects_to_experiment(as.data.frame(d$X[,-1,drop=F]))
        des$overwrite_all_subject_assignments(d$w)
        des$add_all_subject_responses(d$y, deads = if(exists("dead", d)) d$dead else NULL)
        
        inf_cls = get(cls_name)
        # Robustly check for smart_cold_start_default
        inf_obj = tryCatch(inf_cls$new(des, smart_cold_start_default = TRUE), error = function(e) inf_cls$new(des))
        prepare_solver_only_edi(inf_obj, cls_name)
        
        # Determine the benchmark expression
        bench_expr = if (cls_name == "InferenceAllSimpleWilcox") {
            quote({
                inf_obj$.__enclos_env__$private$hl_point_estimate(d$y, d$w)
            })
        } else {
            quote({
                # Clear cache to ensure pure solver execution is timed each iteration
                inf_obj$.__enclos_env__$private$cached_mod = NULL
                inf_obj$.__enclos_env__$private$cached_values = list()
                inf_obj$compute_estimate(estimate_only = TRUE)
            })
        }
        
        # Warmup once
        eval(bench_expr)
        
        collect_timing_ms(bench_expr, times = b_time, fast_path_microbenchmark_reps = fast_path_microbenchmark_reps)
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
  if (pval < 0.05 && speedup < 1) return("#ffd9d9")
  if (pval >= 0.05) return("#fff4bf")
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

report = c(
  "# EDI Exhaustive C++ Model Fit Benchmarks",
  "",
  "This report compares the performance of EDI's Rcpp-optimized model fitting paths against **low-level** canonical R implementations (e.g., `glm.fit`, `lm.fit`, `coxph.fit`) where possible.",
  "",
  "## Benchmark Dataset Specification",
  "",
  "All benchmarks were performed on a synthetic clinical-trial-scale dataset generated for each response type. The data generation process ensures numerical stability and fair solver comparison by using the following parameters:",
  "",
  "*   **Sample Size ($N$):** 1,000 subjects for most models; 500 subjects for survival models. Exact and trend tests may use smaller scaled samples (N=100-500) as noted in the results.",
  "*   **Predictors ($p$):** 5 total predictors, including a global intercept, a balanced binary treatment assignment from fixed `iBCRD`, and 4 continuous covariates ($X \\sim \\text{Normal}(0, 1)$).",
  "*   **Effect Sizes:** Covariate coefficients are sampled from $\\text{Normal}(0, 0.5)$, matching the warm-start benchmark data model. The treatment coefficient is set to 0.5 in the linear predictor so the benchmarked treatment effect is meaningfully separated from zero.",
  "*   **EDI Design Template:** EDI benchmark objects are instantiated on a fixed `iBCRD` design.",
  "*   **Response Generation:**",
  "    *   **Continuous:** Linear model with additive $\\text{Normal}(0, 0.5)$ noise.",
  "    *   **Incidence:** Binary outcomes via a Logistic link.",
  "    *   **Count:** Integer outcomes via Poisson or Negative Binomial distributions with an exponential link.",
  "    *   **Proportion:** Continuous outcomes in $(0, 1)$ via a Beta distribution with a logit link.",
  "    *   **Survival:** Exponentially distributed event times with approximately 20% random censoring.",
  "    *   **Ordinal:** 3-level categorical outcomes generated from the same ordinal construction used in the warm-start benchmark.",
  "",
  "## Methodology",
  "",
  "*   **Pure Solver Timing:** Results reflect the time taken for the core numerical optimization. We exclude R6 object instantiation, design matrix construction, and standard error estimation (which often uses different R-side matrix inversion logic) to isolate the efficiency of the underlying C++ backends.",
  "*   **Solver-Only Prebuilds:** For EDI rows, benchmark setup prebuilds exposed observed-data design matrices, reduced design matrices, and other fixed working inputs outside the timed region when the implementation exposes those hooks. Canonical rows using low-level matrix interfaces are likewise timed on prebuilt inputs.",
  "*   **Smart Cold Starts:** EDI models were initialized with `smart_cold_start = TRUE`, utilizing package-optimized heuristic starting values.",
  "*   **Randomization Design:** EDI timings in this table correspond to `iBCRD` design objects.",
  "*   **Low-Level Comparison:** Canonical R timings use the fastest available internal interfaces (e.g., `.fit` functions) to remove R's formula parsing and environment management overhead.",
  "*   **Limitation:** Some canonical comparators only expose formula-based APIs rather than comparable low-level fit kernels. Those rows remain included, but their canonical timings may still contain formula/model-frame overhead beyond the numerical solver itself.",
  paste0("*   **Averaging:** All timings are medians over ", B_TIME, " warmed runs measured with adaptive batched `system.time`; paths below ", FAST_PATH_THRESHOLD_MS, " ms use `microbenchmark(times = ", FAST_PATH_MICROBENCH_REPS, ")` instead."),
  "*   **Timing P-Value:** `Timing Pval` reports a Welch two-sample t-test comparing the EDI and canonical timing replicate distributions for each row. The unlabeled final column marks thresholds with `***` for p < 0.001, `**` for p < 0.01, and `*` for p < 0.05.",
  "*   **Row Highlighting:** Light green rows indicate `Speedup > 1` and `Timing Pval < 0.05`; light red rows indicate `Speedup < 1` and `Timing Pval < 0.05`; light yellow rows indicate `Timing Pval >= 0.05`; light grey rows indicate `NA` timing comparisons.",
  "*   **Constraints**: Matched-pair/KK and highly custom paths are excluded as per user request.",
  "",
  "## Results",
  "",
  table_lines
)
writeLines(c(report, "", STYLE_BLOCK), "package_metadata/benchmark_model_fits.md")
cat("Benchmark complete.\n")
