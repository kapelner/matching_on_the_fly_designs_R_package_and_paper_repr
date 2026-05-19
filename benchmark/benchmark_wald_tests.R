
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
B_TIME = 1 

# --- Data Generation Helper (Identical to primary benchmark) ---
generate_data = function(n = 200, p = 5, family = "logistic") {
    # coefficients from Normal(0, 0.05) to ensure stability in Hessian calculation
    X = matrix(rnorm(n * p), n, p); X[, 1] = 1 
    beta = rnorm(p) * 0.05 
    eta = X %*% beta
    w = rbinom(n, 1, 0.5); eta = eta + 0.1 * w 
    res = list(X = X, w = w)
    
    if (family %in% c("logistic", "log-binomial", "identity-binomial", "glmm_logistic", "robust", "probit")) {
        prob = if (family == "log-binomial") pmin(0.5, exp(eta - 2)) else plogis(eta)
        res$y = rbinom(n, 1, prob)
    } else if (family %in% c("poisson", "glmm_poisson", "zap", "quasipoisson", "hurdle_poisson_glmm", "modified_poisson")) {
        res$y = rpois(n, exp(pmin(eta, 1)))
    } else if (family %in% c("negbin", "zinb", "truncated_negbin", "hurdle_negbin")) {
        res$y = rnbinom(n, size = 2, mu = exp(pmin(eta, 1)))
    } else if (family %in% c("beta", "zoib", "fractional")) {
        mu = plogis(eta); phi = 10; res$y = pmax(pmin(rbeta(n, mu * phi, (1 - mu) * phi), 1 - 1e-6), 1e-6)
    } else if (family %in% c("weibull", "cox", "strat_cox", "dep_cens")) {
        res$y = rexp(n, 1/exp(pmin(eta, 1))); res$dead = rbinom(n, 1, 0.8)
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

# --- Mapping Table (Full Wald Inference: Fit + SE + P-value) ---
# Note: expressions are modified to include full inference.
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
    list(cls = "InferenceSurvivalCoxPHRegr", pkg = "survival", func = "coxph.fit+Wald", expr = quote({
        x_vars = c("treatment", grep("^x", names(df), value=T))
        res = survival::coxph.fit(x = as.matrix(df[, x_vars, drop=F]), y = survival::Surv(df$y, df$dead), strata=NULL, offset=NULL, init=NULL, control=survival::coxph.control(), weights=NULL, method="efron", rownames=as.character(1:nrow(df)))
        2 * pnorm(-abs(res$coefficients[1] / sqrt(res$var[1,1])))
    })),
    list(cls = "InferenceCountNegBin", pkg = "MASS", func = "glm.nb+summary", expr = quote({
        res = MASS::glm.nb(y ~ treatment + x1 + x2 + x3 + x4, data = df)
        summary(res)$coefficients[2, 4]
    })),
    list(cls = "InferencePropBetaRegr", pkg = "betareg", func = "betareg.fit+Wald", expr = quote({
        res = betareg::betareg.fit(x = X_can, y = df$y)
        # betareg.fit doesn't give SEs easily. Use summary if needed or manual Hessian.
        # summary(betareg(...)) is most common.
        summary(betareg::betareg(df$y ~ X_can[,-1]))$coefficients$mean[2, 4]
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
        summary(res)$coefficients[2, 3] # rlm gives t-stat, p-value needs derivation
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
    list(cls = "InferenceOrdinalAdjCatLogitRegr", pkg = "VGAM", func = "vglm+summary", expr = quote({
        res = VGAM::vglm(factor(y, ordered=T) ~ treatment + x1 + x2 + x3 + x4, VGAM::acat(), data=df)
        summary(res)@coef3[1, 4]
    })),
    list(cls = "InferenceOrdinalContRatioRegr", pkg = "VGAM", func = "vglm+summary", expr = quote({
        res = VGAM::vglm(factor(y, ordered=T) ~ treatment + x1 + x2 + x3 + x4, VGAM::cratio(), data=df)
        summary(res)@coef3[1, 4]
    })),
    list(cls = "InferenceSurvivalLogRank", pkg = "survival", func = "survdiff", expr = quote(survival::survdiff(survival::Surv(y, dead) ~ treatment, data = df))),
    list(cls = "InferenceSurvivalGehanWilcox", pkg = "survival", func = "survdiff(rho=1)", expr = quote(survival::survdiff(survival::Surv(y, dead) ~ treatment, data = df, rho = 1))),
    list(cls = "InferenceAllSimpleMeanDiffPooledVar", pkg = "stats", func = "t.test(pool)", expr = quote(t.test(df$y[df$treatment==1], df$y[df$treatment==0], var.equal = TRUE))),
    list(cls = "InferenceAllSimpleWilcox", pkg = "stats", func = "wilcox.test", expr = quote(wilcox.test(df$y[df$treatment==1], df$y[df$treatment==0]))),
    list(cls = "InferenceIncidExactFisher", pkg = "stats", func = "fisher.test", expr = quote(fisher.test(table(df$treatment, df$y)))),
    list(cls = "InferenceIncidCMH", pkg = "stats", func = "mantelhaen", expr = quote(mantelhaen.test(table(df$treatment, df$y, df$g)))),
    list(cls = "InferenceOrdinalJonckheereTerpstraTest", pkg = "clinfun", func = "jonckheere", expr = quote(clinfun::jonckheere.test(df$y, df$treatment))),
    list(cls = "InferenceIncidMiettinenNurminenRiskDiff", pkg = "DescTools", func = "BinomDiffCI(mn)", expr = quote(DescTools::BinomDiffCI(sum(df$y[df$treatment==1]), sum(df$treatment==1), sum(df$y[df$treatment==0]), sum(df$treatment==0), method="mn"))),
    list(cls = "InferenceSurvivalStratCoxPHRegr", pkg = "survival", func = "coxph.fit(strat)+Wald", expr = quote({
        res = survival::coxph.fit(x = as.matrix(df[, grep("^x", names(df)), drop=F]), y = survival::Surv(df$y, df$dead), strata=as.integer(df$g), offset=NULL, init=NULL, control=survival::coxph.control(), weights=NULL, method="efron", rownames=as.character(1:nrow(df)))
        2 * pnorm(-abs(res$coefficients[1] / sqrt(res$var[1,1])))
    })),
    list(cls = "InferenceContinLin", pkg = "stats", func = "lm.fit(interact)+Wald", expr = quote({
        X_int = model.matrix(~ treatment * (x1 + x2 + x3 + x4), data = df)
        res = lm.fit(x = X_int, y = df$y)
        p = ncol(X_int); n = nrow(X_int); rss = sum(res$residuals^2); df_res = n - p; sig2 = rss / df_res
        v = chol2inv(res$qr$qr[1:p, 1:p, drop=FALSE]) * sig2; 2 * pt(-abs(res$coefficients[2] / sqrt(v[2,2])), df_res)
    })),
    list(cls = "InferenceIncidRiskDiff", pkg = "stats", func = "prop.test", expr = quote(prop.test(table(df$treatment, df$y)))),
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
        # SE for Ridit is usually bootstrap, just time the point estimate for consistency with primary table
        y = df$y; tab = table(y); cum = cumsum(tab); prev = c(0, cum[-length(cum)])
        ridit_map = (prev + 0.5 * tab) / length(y); r = as.numeric(ridit_map[as.character(y)])
        mean(r[df$treatment==1]) - mean(r[df$treatment==0])
    })),
    list(cls = "InferenceIncidExactZhang", pkg = "Exact", func = "exact.test(z)", expr = quote(Exact::exact.test(table(df$treatment, df$y), method="z-pooled", tsmethod="central"))),
    list(cls = "InferenceIncidNewcombeRiskDiff", pkg = "DescTools", func = "BinomDiffCI(score)", expr = quote(DescTools::BinomDiffCI(sum(df$y[df$treatment==1]), sum(df$treatment==1), sum(df$y[df$treatment==0]), sum(df$treatment==0), method="score")))
)

# --- Benchmark Runner ---
results = list()

run_one = function(spec) {
    cls_name = spec$cls
    cat(sprintf("Benchmarking Wald Test %s...\n", cls_name))
    
    resp_type = "continuous"
    family = "continuous"
    if (grepl("Ordinal|AdjCat|ContRatio|Stereotype|Ridit|Sign", cls_name)) { resp_type = "ordinal"; family = "ordinal" }
    else if (grepl("Incid|Binomial|Wald|CMH|Fisher|Zhang|Robins|Newcombe|Nurminen", cls_name)) { resp_type = "incidence"; family = "logistic" }
    else if (grepl("Count|Poisson|NegBin|ZINB|ZAP|Hurdle", cls_name)) { resp_type = "count"; family = "poisson" }
    else if (grepl("Prop|Beta|ZOIB|Fractional", cls_name)) { resp_type = "proportion"; family = "beta" }
    else if (grepl("Survival|Cox|Weibull|KM|Rank|LogRank|Gehan|RMST|RMDiff|LWACox|Clayton", cls_name)) { resp_type = "survival"; family = "cox" }
    
    if (cls_name == "InferenceIncidLogBinomial") family = "log-binomial"
    
    scale = if (!is.null(spec$scale)) spec$scale else 1.0
    n = round(N_GLM * scale)
    if (grepl("Survival", cls_name)) n = round(N_SURV * scale)
    d = generate_data(n = n, family = family)
    
    # Timing EDI (Model Fit + Wald Test)
    t_edi = tryCatch({
        des = DesignFixedBernoulli$new(n = n, response_type = resp_type)
        des$add_all_subjects_to_experiment(as.data.frame(d$X[,-1,drop=F]))
        des$overwrite_all_subject_assignments(d$w)
        des$add_all_subject_responses(d$y, deads = if(exists("dead", d)) d$dead else NULL)
        
        inf_cls = get(cls_name)
        inf_obj = tryCatch(inf_cls$new(des, smart_cold_start_default = TRUE), error = function(e) inf_cls$new(des))
        
        bench_expr = quote({
            # Clear cache
            inf_obj$.__enclos_env__$private$cached_mod = NULL
            inf_obj$.__enclos_env__$private$cached_values = list()
            # Compute full estimate including variance/SE
            inf_obj$compute_estimate(estimate_only = FALSE)
            # Derive p-value
            inf_obj$compute_asymp_two_sided_pval(delta = 0)
        })
        
        # Warmup
        eval(bench_expr)
        median(microbenchmark(eval(bench_expr), times = B_TIME)$time) / 1e6
    }, error = function(e) { cat("  EDI Error:", e$message, "\n"); NA_real_ })
    
    # Timing Canonical
    t_can = tryCatch({
        if (!is.null(spec$expr) && spec$pkg != "None") {
            library(spec$pkg, character.only = TRUE)
            X_cols = d$X[,-1,drop=F]; colnames(X_cols) = paste0("x", 1:ncol(X_cols))
            df = data.frame(y = d$y, treatment = d$w, g = factor(rep(1:10, length.out = n)), dead = if(exists("dead", d)) d$dead else 1)
            df = cbind(df, X_cols); X_can = cbind(`(Intercept)` = 1, treatment = d$w, as.matrix(X_cols))
            
            median(microbenchmark(eval(spec$expr), times = B_TIME)$time) / 1e6
        } else NA_real_
    }, error = function(e) { cat("  Canonical Error:", e$message, "\n"); NA_real_ })
    
    results[[length(results) + 1]] <<- data.table(
        Class = cls_name, Response = resp_type, EDI_Time_ms = t_edi,
        Canonical_Pkg = spec$pkg, Canonical_Func = spec$func, Canonical_Time_ms = t_can,
        Speedup = if (!is.na(t_can) && !is.na(t_edi) && t_edi > 0) t_can / t_edi else NA_real_
    )
}

unique_specs = bench_specs[!duplicated(sapply(bench_specs, `[[`, "cls"))]
for (i in seq_along(unique_specs)) {
    cat(sprintf("[%d/%d] ", i, length(unique_specs)))
    run_one(unique_specs[[i]])
}

# Finalize
dt = rbindlist(results)
dt[, Speedup := ifelse(!is.na(Speedup), paste0(round(Speedup, 2), "x"), "NA")]
dt[, EDI_Time_ms := round(EDI_Time_ms, 2)]
dt[, Canonical_Time_ms := round(Canonical_Time_ms, 2)]

header = c(
  "## Wald Test Performance (Full Inference)",
  "",
  "This table compares the performance of **Full Inference** (Model Fit + Standard Error calculation + P-value derivation).",
  "Unlike the point-estimation table above, these results include the computational cost of the variance-covariance matrix (Hessian or Fisher Information) and the Wald test statistic calculation.",
  "",
  "| Class | Response | EDI Time (ms) | Canonical Pkg | Canonical Func | Canonical Time (ms) | Speedup |",
  "| :--- | :--- | :--- | :--- | :--- | :--- | :--- |"
)
rows = paste0("| ", dt$Class, " | ", dt$Response, " | ", dt$EDI_Time_ms, " | ", 
              dt$Canonical_Pkg, " | ", dt$Canonical_Func, " | ", dt$Canonical_Time_ms, " | ", 
              dt$Speedup, " |")

# Append to file
write(c("", header, rows), file = "package_metadata/benchmark_model_fits.md", append = TRUE)
cat("Wald Test Benchmark complete.\n")
