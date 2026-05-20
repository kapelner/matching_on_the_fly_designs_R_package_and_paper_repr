
if (!requireNamespace("pkgload", quietly = TRUE)) stop("The 'pkgload' package is required.")
pkgload::load_all("EDI", quiet = TRUE)
library(EDI)
library(microbenchmark)
library(data.table)
library(survival)

set.seed(42)

# Global Config
N_GLM = 1000
N_SURV = 500
B_TIME = 3 
TARGET_BATCH_MS = 50
MIN_RESOLVED_BATCH_MS = 10
MAX_INNER_REPS = 100000L
FAST_PATH_MICROBENCH_REPS = 500L
FAST_PATH_THRESHOLD_MS = 0.01

# --- Data Generation Helper ---
generate_data = function(n = 200, p = 5, family = "logistic") {
    X = matrix(rnorm(n * p), n, p); X[, 1] = 1 
    beta = rnorm(p) * 0.1 
    eta = X %*% beta
    n_treat = floor(n / 2)
    w = sample(c(rep(1, n_treat), rep(0, n - n_treat)))
    eta = eta + 0.1 * w 
    res = list(X = X, w = w)
    
    if (family %in% c("logistic", "log-binomial", "identity-binomial", "glmm_logistic", "robust", "probit")) {
        prob = if (family == "log-binomial") pmin(0.5, exp(eta - 2)) else plogis(eta)
        res$y = rbinom(n, 1, prob)
    } else if (family %in% c("poisson", "glmm_poisson", "zap", "quasipoisson", "hurdle_poisson_glmm", "modified_poisson")) {
        res$y = rpois(n, exp(pmin(eta, 2)))
    } else if (family %in% c("negbin", "zinb", "truncated_negbin", "hurdle_negbin")) {
        res$y = rnbinom(n, size = 2, mu = exp(pmin(eta, 2)))
    } else if (family %in% c("beta", "zoib", "fractional")) {
        mu = plogis(eta); phi = 10; res$y = pmax(pmin(rbeta(n, mu * phi, (1 - mu) * phi), 1 - 1e-6), 1e-6)
    } else if (family %in% c("weibull", "cox", "strat_cox", "dep_cens")) {
        res$y = rexp(n, 1/exp(pmin(eta, 2))); res$dead = rbinom(n, 1, 0.8)
    } else if (family %in% c("ordinal", "adj_cat", "continuation", "stereotype", "glmm_ordinal")) {
        p1 = plogis(-1.5 - eta); p_le_2 = plogis(0 - eta); p_le_3 = plogis(1.5 - eta)
        p2 = pmax(0, p_le_2 - p1); p3 = pmax(0, p_le_3 - p_le_2); p4 = pmax(0, 1 - p_le_3)
        probs = cbind(p1, p2, p3, p4); probs = probs / rowSums(probs)
        res$y = apply(probs, 1, function(p) sample(1:4, 1, prob = p))
    } else {
        res$y = as.numeric(eta + rnorm(n, 0, 0.5))
    }
    res
}

time_expr_ms = function(expr, times = B_TIME, env = parent.frame(), target_batch_ms = TARGET_BATCH_MS, max_inner_reps = MAX_INNER_REPS, fast_path_microbenchmark_reps = FAST_PATH_MICROBENCH_REPS) {
    micro_time_ms = function() {
        mb_ms = microbenchmark(eval(expr, envir = env), times = fast_path_microbenchmark_reps)$time / 1e6
        mb_ms = mb_ms[is.finite(mb_ms) & mb_ms > 0]
        if (length(mb_ms) == 0L) 0 else median(mb_ms)
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
    median_ms = median(vals_ms[is.finite(vals_ms)])
    if (!is.finite(median_ms) || median_ms < FAST_PATH_THRESHOLD_MS) {
        return(micro_time_ms())
    }
    median_ms
}

# --- Mapping Table (Targeted subset requested by user) ---
# library() calls are used to ensure it crashes if pkgs are missing.
bench_specs = list(
    list(cls = "InferenceIncidLogRegr", pkg = "stats", func = "glm.fit", expr = quote(glm.fit(x = X_can, y = df$y, family = binomial()))),
    list(cls = "InferenceContinOLS", pkg = "stats", func = "lm.fit", expr = quote(lm.fit(x = X_can, y = df$y))),
    list(cls = "InferenceCountPoisson", pkg = "stats", func = "glm.fit", expr = quote(glm.fit(x = X_can, y = df$y, family = poisson()))),
    list(cls = "InferenceSurvivalCoxPHRegr", pkg = "survival", func = "coxph.fit", expr = quote({
        x_vars = c("treatment", grep("^x", names(df), value=T))
        survival::coxph.fit(x = as.matrix(df[, x_vars, drop=F]), y = survival::Surv(df$y, df$dead), strata=NULL, offset=NULL, init=NULL, control=survival::coxph.control(), weights=NULL, method="efron", rownames=as.character(1:nrow(df)))
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
    list(cls = "InferenceContinRobustRegr", pkg = "MASS", func = "rlm", expr = quote(MASS::rlm(x = X_can, y = df$y))),
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
    list(cls = "InferenceOrdinalCloglogRegr", pkg = "ordinal", func = "clm(cll)", expr = quote(ordinal::clm(factor(y, ordered=T) ~ treatment + x1 + x2 + x3 + x4, data = df, link = "cloglog"))),
    list(cls = "InferenceOrdinalCauchitRegr", pkg = "ordinal", func = "clm(cauchit)", expr = quote(ordinal::clm(factor(y, ordered=T) ~ treatment + x1 + x2 + x3 + x4, data = df, link = "cauchit"))),
    list(cls = "InferenceSurvivalLogRank", pkg = "survival", func = "survdiff", expr = quote(survival::survdiff(survival::Surv(y, dead) ~ treatment, data = df))),
    list(cls = "InferenceSurvivalGehanWilcox", pkg = "survival", func = "survdiff(rho=1)", expr = quote(survival::survdiff(survival::Surv(y, dead) ~ treatment, data = df, rho = 1))),
    list(cls = "InferenceAllSimpleMeanDiffPooledVar", pkg = "stats", func = "t.test(pool)", expr = quote(t.test(df$y[df$treatment==1], df$y[df$treatment==0], var.equal = TRUE))),
    list(cls = "InferenceAllSimpleWilcox", pkg = "stats", func = "wilcox.test", expr = quote(wilcox.test(df$y[df$treatment==1], df$y[df$treatment==0])), scale = 0.5),
    list(cls = "InferenceIncidExactFisher", pkg = "stats", func = "fisher.test", expr = quote(fisher.test(table(df$treatment, df$y))), scale = 0.1),
    list(cls = "InferenceOrdinalJonckheereTerpstraTest", pkg = "clinfun", func = "jonckheere", expr = quote(clinfun::jonckheere.test(df$y, df$treatment))),
    list(cls = "InferenceIncidMiettinenNurminenRiskDiff", pkg = "DescTools", func = "BinomDiffCI(mn)", expr = quote(DescTools::BinomDiffCI(sum(df$y[df$treatment==1]), sum(df$treatment==1), sum(df$y[df$treatment==0]), sum(df$treatment==0), method="mn"))),
    list(cls = "InferenceIncidModifiedPoisson", pkg = "stats", func = "glm.fit(modified)", expr = quote(glm.fit(x = X_can, y = df$y, family = poisson()))),
    list(cls = "InferenceSurvivalStratCoxPHRegr", pkg = "survival", func = "coxph.fit(strat)", expr = quote(survival::coxph.fit(x = as.matrix(df[, grep("^x", names(df)), drop=F]), y = survival::Surv(df$y, df$dead), strata=as.integer(df$g), offset=NULL, init=NULL, control=survival::coxph.control(), weights=NULL, method="efron", rownames=as.character(1:nrow(df))))),
    list(cls = "InferenceContinLin", pkg = "stats", func = "lm.fit(interact)", expr = quote({
        X_int = model.matrix(~ treatment * (x1 + x2 + x3 + x4), data = df)
        lm.fit(x = X_int, y = df$y)
    })),
    list(cls = "InferenceIncidRiskDiff", pkg = "stats", func = "prop.test", expr = quote({
        tab = table(df$treatment, df$y)
        prop.test(tab)$estimate[2] - prop.test(tab)$estimate[1]
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
    list(cls = "InferenceOrdinalRidit", pkg = "stats", func = "mean(ridit)", expr = quote({
        y = df$y
        tab = table(y)
        cum = cumsum(tab)
        prev = c(0, cum[-length(cum)])
        ridit_map = (prev + 0.5 * tab) / length(y)
        r = as.numeric(ridit_map[as.character(y)])
        mean(r[df$treatment==1]) - mean(r[df$treatment==0])
    })),
    list(cls = "InferenceIncidExactZhang", pkg = "Exact", func = "exact.test(z)", expr = quote({
        Exact::exact.test(table(df$treatment, df$y), method="z-pooled", tsmethod="central")
    }), fast_path_microbenchmark_reps = 5000L, b_time_override = 20L),
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
    list(cls = "InferenceSurvivalDepCensTransformRegr", pkg = "None", func = "None", expr = NULL),
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
    })),
    list(cls = "InferencePropZeroOneInflatedBetaRegr", pkg = "None", func = "None", expr = NULL)
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
    t_edi = tryCatch({
        des = DesignFixediBCRD$new(n = n, response_type = resp_type)
        des$add_all_subjects_to_experiment(as.data.frame(d$X[,-1,drop=F]))
        des$overwrite_all_subject_assignments(d$w)
        des$add_all_subject_responses(d$y, deads = if(exists("dead", d)) d$dead else NULL)
        
        inf_cls = get(cls_name)
        # Robustly check for smart_cold_start_default
        inf_obj = tryCatch(inf_cls$new(des, smart_cold_start_default = TRUE), error = function(e) inf_cls$new(des))
        
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
        
        time_expr_ms(bench_expr, times = b_time, fast_path_microbenchmark_reps = fast_path_microbenchmark_reps)
    }, error = function(e) {
        cat("  EDI Error:", e$message, "\n")
        NA_real_
    })

    # Timing Canonical (Crashes if pkg missing)
    t_can = NA_real_
    if (!is.null(spec$expr) && spec$pkg != "None") {
        library(spec$pkg, character.only = TRUE)
        X_cols = d$X[,-1,drop=F]
        colnames(X_cols) = paste0("x", 1:ncol(X_cols))
        df = data.frame(y = d$y, treatment = d$w, g = factor(rep(1:10, length.out = n)), dead = if(exists("dead", d)) d$dead else 1)
        df = cbind(df, X_cols)
        X_can = cbind(`(Intercept)` = 1, treatment = d$w, as.matrix(X_cols))

        t_can = time_expr_ms(spec$expr, times = b_time, fast_path_microbenchmark_reps = fast_path_microbenchmark_reps)
    }
    
    results[[length(results) + 1]] <<- data.table(
        Class = cls_name, Response = resp_type, EDI_Time_ms = t_edi,
        Canonical_Pkg = spec$pkg, Canonical_Func = spec$func, Canonical_Time_ms = t_can,
        Speedup = if (!is.na(t_can) && !is.na(t_edi) && t_edi > 0) t_can / t_edi else NA_real_
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
dt[, Speedup := ifelse(!is.na(Speedup), paste0(round(Speedup, 2), "x"), "NA")]
format_ms = function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(x < 0.01, format(x, scientific = TRUE, digits = 3, trim = TRUE), format(round(x, 2), nsmall = 2, trim = TRUE))
  )
}
dt[, EDI_Time_ms := format_ms(EDI_Time_ms)]
dt[, Canonical_Time_ms := format_ms(Canonical_Time_ms)]
dt[, Response := as.character(Response)]

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
  "*   **Effect Sizes:** Coefficients are sampled from $\\text{Normal}(0, 0.1)$ to ensure reasonable event rates and avoid separation issues in logistic/ordinal models.",
  "*   **EDI Design Template:** EDI benchmark objects are instantiated on a fixed `iBCRD` design.",
  "*   **Response Generation:**",
  "    *   **Continuous:** Linear model with additive $\\text{Normal}(0, 0.5)$ noise.",
  "    *   **Incidence:** Binary outcomes via a Logistic link.",
  "    *   **Count:** Integer outcomes via Poisson or Negative Binomial distributions with an exponential link.",
  "    *   **Proportion:** Continuous outcomes in $(0, 1)$ via a Beta distribution with a logit link.",
  "    *   **Survival:** Exponentially distributed event times with approximately 20% random censoring.",
  "    *   **Ordinal:** 4-level categorical outcomes generated from a Proportional Odds model.",
  "",
  "## Methodology",
  "",
  "*   **Pure Solver Timing:** Results reflect the time taken for the core numerical optimization. We exclude R6 object instantiation, design matrix construction, and standard error estimation (which often uses different R-side matrix inversion logic) to isolate the efficiency of the underlying C++ backends.",
  "*   **Smart Cold Starts:** EDI models were initialized with `smart_cold_start = TRUE`, utilizing package-optimized heuristic starting values.",
  "*   **Randomization Design:** EDI timings in this table correspond to `iBCRD` design objects.",
  "*   **Low-Level Comparison:** Canonical R timings use the fastest available internal interfaces (e.g., `.fit` functions) to remove R's formula parsing and environment management overhead.",
  paste0("*   **Averaging:** All timings are medians over ", B_TIME, " warmed runs measured with adaptive batched `system.time`; paths below ", FAST_PATH_THRESHOLD_MS, " ms use `microbenchmark(times = ", FAST_PATH_MICROBENCH_REPS, ")` instead."),
  "*   **Constraints**: Matched-pair/KK and highly custom paths are excluded as per user request.",
  "",
  "## Results",
  "",
  "| Class | Response | EDI Time (ms) | Canonical Pkg | Canonical Func | Canonical Time (ms) | Speedup |",
  "| :--- | :--- | :--- | :--- | :--- | :--- | :--- |",
  paste0("| ", dt$Class, " | ", dt$Response, " | ", dt$EDI_Time_ms, " | ", 
         dt$Canonical_Pkg, " | ", dt$Canonical_Func, " | ", dt$Canonical_Time_ms, " | ", 
         dt$Speedup, " |")
)
writeLines(report, "package_metadata/benchmark_model_fits.md")
cat("Benchmark complete.\n")
