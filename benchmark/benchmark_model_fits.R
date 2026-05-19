
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

# --- Data Generation Helper ---
generate_data = function(n = 200, p = 5, family = "logistic", n_groups = 10) {
    X = matrix(rnorm(n * p), n, p); X[, 1] = 1 
    beta = rnorm(p) * 0.1 
    eta = X %*% beta
    w = rbinom(n, 1, 0.5); eta = eta + 0.1 * w 
    res = list(X = X, group_id = rep(1:n_groups, length.out = n), strata = rep(1:5, length.out = n), w = w)
    
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

# --- Mapping Table ---
bench_specs = list(
    list(cls = "InferenceIncidLogRegr", pkg = "stats", func = "glm.fit", expr = quote(glm.fit(x = X_can, y = df$y, family = binomial()))),
    list(cls = "InferenceContinOLS", pkg = "stats", func = "lm.fit", expr = quote(lm.fit(x = X_can, y = df$y))),
    list(cls = "InferenceCountPoisson", pkg = "stats", func = "glm.fit", expr = quote(glm.fit(x = X_can, y = df$y, family = poisson()))),
    list(cls = "InferenceSurvivalCoxPHRegr", pkg = "survival", func = "coxph.fit", expr = quote({
        survival::coxph.fit(x = X_can[,-1,drop=F], y = survival::Surv(df$y, df$dead), strata=NULL, offset=NULL, init=NULL, control=survival::coxph.control(), weights=NULL, method="efron", rownames=as.character(1:nrow(df)))
    })),
    list(cls = "InferenceCountNegBin", pkg = "MASS", func = "glm.nb", expr = quote(MASS::glm.nb(y ~ w + x1 + x2 + x3 + x4, data = df))),
    list(cls = "InferencePropBetaRegr", pkg = "betareg", func = "betareg.fit", expr = quote(betareg::betareg.fit(x = X_can, y = df$y))),
    list(cls = "InferenceOrdinalPropOddsRegr", pkg = "ordinal", func = "clm", expr = quote(ordinal::clm(factor(y, ordered=T) ~ w + x1 + x2 + x3 + x4, data = df))),
    list(cls = "InferenceCountHurdlePoisson", pkg = "pscl", func = "hurdle", expr = quote(pscl::hurdle(y ~ w + x1 + x2 + x3 + x4, data = df))),
    list(cls = "InferenceCountZeroInflatedPoisson", pkg = "pscl", func = "zeroinfl", expr = quote(pscl::zeroinfl(y ~ w + x1 + x2 + x3 + x4, data = df))),
    list(cls = "InferenceCountZeroInflatedNegBin", pkg = "pscl", func = "zeroinfl(nb)", expr = quote(pscl::zeroinfl(y ~ w + x1 + x2 + x3 + x4, data = df, dist="negbin"))),
    list(cls = "InferenceCountHurdleNegBin", pkg = "pscl", func = "hurdle(nb)", expr = quote(pscl::hurdle(y ~ w + x1 + x2 + x3 + x4, data = df, dist="negbin"))),
    list(cls = "InferenceCountQuasiPoisson", pkg = "stats", func = "glm.fit(quasi)", expr = quote(glm.fit(x = X_can, y = df$y, family = quasipoisson()))),
    list(cls = "InferenceSurvivalWeibullRegr", pkg = "survival", func = "survreg", expr = quote(survival::survreg(survival::Surv(y, dead) ~ w + x1 + x2 + x3 + x4, data = df, dist="weibull"))),
    list(cls = "InferenceContinRobustRegr", pkg = "MASS", func = "rlm", expr = quote(MASS::rlm(x = X_can, y = df$y))),
    list(cls = "InferenceContinQuantileRegr", pkg = "quantreg", func = "rq.fit", expr = quote(quantreg::rq.fit(x = X_can, y = df$y))),
    list(cls = "InferencePropFractionalLogit", pkg = "stats", func = "glm.fit(quasi)", expr = quote(glm.fit(x = X_can, y = df$y, family=quasibinomial()))),
    list(cls = "InferenceIncidLogBinomial", pkg = "stats", func = "glm.fit(log)", expr = quote(glm.fit(x = X_can, y = df$y, family=binomial(link="log"), start=c(-2, rep(0, ncol(X_can)-1))))),
    list(cls = "InferenceIncidProbitRegr", pkg = "stats", func = "glm.fit(probit)", expr = quote(glm.fit(x = X_can, y = df$y, family=binomial(link="probit")))),
    list(cls = "InferenceIncidBinomialIdentityRiskDiff", pkg = "stats", func = "glm.fit(ident)", expr = quote(glm.fit(x = X_can, y = df$y, family=binomial(link="identity"), start=c(0.5, rep(0, ncol(X_can)-1))))),
    list(cls = "InferenceOrdinalAdjCatLogitRegr", pkg = "VGAM", func = "vglm(acat)", expr = quote(VGAM::vglm(factor(y, ordered=T) ~ w + x1 + x2 + x3 + x4, VGAM::acat(), data=df))),
    list(cls = "InferenceOrdinalContRatioRegr", pkg = "VGAM", func = "vglm(cratio)", expr = quote(VGAM::vglm(factor(y, ordered=T) ~ w + x1 + x2 + x3 + x4, VGAM::cratio(), data=df))),
    list(cls = "InferenceOrdinalOrderedProbitRegr", pkg = "ordinal", func = "clm(probit)", expr = quote(ordinal::clm(factor(y, ordered=T) ~ w + x1 + x2 + x3 + x4, data = df, link = "probit"))),
    list(cls = "InferenceOrdinalCloglogRegr", pkg = "ordinal", func = "clm(cll)", expr = quote(ordinal::clm(factor(y, ordered=T) ~ w + x1 + x2 + x3 + x4, data = df, link = "cloglog"))),
    list(cls = "InferenceOrdinalCauchitRegr", pkg = "ordinal", func = "clm(cauchit)", expr = quote(ordinal::clm(factor(y, ordered=T) ~ w + x1 + x2 + x3 + x4, data = df, link = "cauchit"))),
    list(cls = "InferenceSurvivalLogRank", pkg = "survival", func = "survdiff", expr = quote(survival::survdiff(survival::Surv(y, dead) ~ w, data = df))),
    list(cls = "InferenceSurvivalGehanWilcox", pkg = "survival", func = "survdiff(rho=1)", expr = quote(survival::survdiff(survival::Surv(y, dead) ~ w, data = df, rho = 1))),
    list(cls = "InferenceAllSimpleMeanDiff", pkg = "base", func = "mean diff", expr = quote(mean(df$y[df$w==1]) - mean(df$y[df$w==0]))),
    list(cls = "InferenceAllSimpleMeanDiffPooledVar", pkg = "stats", func = "t.test(pool)", expr = quote(t.test(df$y[df$w==1], df$y[df$w==0], var.equal = TRUE))),
    list(cls = "InferenceAllSimpleWilcox", pkg = "stats", func = "wilcox.test", expr = quote(wilcox.test(df$y[df$w==1], df$y[df$w==0])), scale = 0.5),
    list(cls = "InferenceIncidExactFisher", pkg = "stats", func = "fisher.test", expr = quote(fisher.test(table(df$w, df$y))), scale = 0.1),
    list(cls = "InferenceIncidCMH", pkg = "stats", func = "mantelhaen", expr = quote(mantelhaen.test(table(df$w, df$y, df$g)))),
    list(cls = "InferenceOrdinalPairedSignTest", pkg = "stats", func = "binom.test", expr = quote(binom.test(sum(df$y[df$w==1] > df$y[df$w==0]), length(df$y)/2)))
)

# --- Benchmark Runner ---
results = list()

run_one = function(spec) {
    cls_name = spec$cls
    cat(sprintf("Benchmarking %s...\n", cls_name))
    
    # Improved heuristic mapping
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
    n = round(N_GLM * scale)
    if (grepl("Survival", cls_name)) n = round(N_SURV * scale)
    
    d = generate_data(n = n, family = family)
    
    # Timing EDI
    t_edi = tryCatch({
        des = DesignFixedBernoulli$new(n = n, response_type = resp_type)
        des$add_all_subjects_to_experiment(as.data.frame(d$X[,-1,drop=F]))
        des$overwrite_all_subject_assignments(d$w)
        des$add_all_subject_responses(d$y, deads = d$dead)
        
        inf_cls = get(cls_name)
        inf_obj = if ("smart_cold_start_default" %in% names(formals(inf_cls$public_methods$initialize))) {
            inf_cls$new(des, smart_cold_start_default = TRUE)
        } else {
            inf_cls$new(des)
        }
        
        bench_expr = if (cls_name == "InferenceAllSimpleWilcox") {
            quote(inf_obj$.__enclos_env__$private$hl_point_estimate(d$y, d$w))
        } else {
            quote(inf_obj$compute_estimate(estimate_only = TRUE))
        }
        
        median(microbenchmark(eval(bench_expr), times = B_TIME)$time) / 1e6
    }, error = function(e) NA_real_)
    
    # Timing Canonical
    t_can = NA_real_
    if (!is.null(spec$expr) && requireNamespace(spec$pkg, quietly = TRUE)) {
        X_cols = d$X[,-1,drop=F]
        colnames(X_cols) = paste0("x", 1:ncol(X_cols))
        df = data.frame(y = d$y, w = d$w, g = factor(d$group_id), dead = if(!is.null(d$dead)) d$dead else 1)
        df = cbind(df, X_cols)
        X_can = cbind(`(Intercept)` = 1, treatment = d$w, as.matrix(X_cols))
        
        t_can = tryCatch({
            median(microbenchmark(eval(spec$expr), times = B_TIME)$time) / 1e6
        }, error = function(e) NA_real_)
    }
    
    results[[length(results) + 1]] <<- data.table(
        Class = cls_name, Response = resp_type, EDI_Time_ms = t_edi,
        Canonical_Pkg = if(is.null(spec$pkg)) "None" else spec$pkg, 
        Canonical_Func = if(is.null(spec$func)) "None" else spec$func, 
        Canonical_Time_ms = t_can,
        Speedup = if (!is.na(t_can) && !is.na(t_edi) && t_edi > 0) t_can / t_edi else NA_real_
    )
}

add_res = function(cls, resp, t_edi, cpkg, cfunc, t_can) {
    results[[length(results) + 1]] <<- data.table(
        Class = cls, Response = resp, EDI_Time_ms = as.numeric(t_edi),
        Canonical_Pkg = as.character(cpkg), Canonical_Func = as.character(cfunc), 
        Canonical_Time_ms = as.numeric(t_can),
        Speedup = if (!is.na(t_can) && !is.na(t_edi) && t_edi > 0) t_can / t_edi else NA_real_
    )
}

# 1. Run defined specs
covered_initial = sapply(bench_specs, `[[`, "cls")
all_inf = grep("^Inference", getNamespaceExports("EDI"), value = TRUE)
remaining = setdiff(all_inf[!grepl("Abstract|Mixin|Suite|IVWC|OneLik|KK", all_inf)], covered_initial)

for (i in seq_along(bench_specs)) {
    cat(sprintf("[%d/%d] ", i, length(bench_specs) + length(remaining)))
    run_one(bench_specs[[i]])
}

# 2. Sweep remaining paths
cat("Sweeping remaining paths...\n")
for (i in seq_along(remaining)) {
    cat(sprintf("[%d/%d] ", length(bench_specs) + i, length(bench_specs) + length(remaining)))
    run_one(list(cls = remaining[i], pkg = "None", func = "None", expr = NULL))
}

# Finalize
dt = rbindlist(results)
dt[, Speedup := ifelse(!is.na(Speedup), paste0(round(Speedup, 2), "x"), "NA")]
dt[, EDI_Time_ms := round(EDI_Time_ms, 2)]
dt[, Canonical_Time_ms := round(Canonical_Time_ms, 2)]

report = c(
  "# EDI Exhaustive Model Fit Benchmarks",
  "",
  "This report compares the performance of EDI's Rcpp-optimized model fitting paths against **low-level** canonical R implementations (e.g., `glm.fit`, `lm.fit`, `coxph.fit`) where possible.",
  "Timings represent pure solver execution (excluding R6 object instantiation overhead) with `smart_cold_start = TRUE` enabled for EDI.",
  "",
  "| Class | Response | EDI Time (ms) | Canonical Pkg | Canonical Func | Canonical Time (ms) | Speedup |",
  "| :--- | :--- | :--- | :--- | :--- | :--- | :--- |",
  paste0("| ", dt$Class, " | ", dt$Response, " | ", dt$EDI_Time_ms, " | ", 
         dt$Canonical_Pkg, " | ", dt$Canonical_Func, " | ", dt$Canonical_Time_ms, " | ", 
         dt$Speedup, " |")
)
writeLines(report, "package_metadata/benchmark_model_fits.md")
cat("Benchmark complete.\n")
