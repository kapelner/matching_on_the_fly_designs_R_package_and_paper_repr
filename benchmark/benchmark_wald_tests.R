
if (!requireNamespace("pkgload", quietly = TRUE)) stop("The 'pkgload' package is required.")
compiler::enableJIT(0) # Disable JIT to prevent R6 on-the-fly compilation overhead
pkgload::load_all("EDI", quiet = TRUE)
library(EDI)
library(microbenchmark)
library(data.table)
library(survival)

set.seed(42)

# Global Config for Wald Tests (Full Inference)
N_WALD = 200 
B_TIME = 3 
TARGET_BATCH_MS = 50
MIN_RESOLVED_BATCH_MS = 10
MAX_INNER_REPS = 100000L
FAST_PATH_MICROBENCH_REPS = 500L
FAST_PATH_THRESHOLD_MS = 0.01

# --- Data Generation Helper ---
generate_data = function(n = 200, p = 5, family = "logistic") {
    X = matrix(rnorm(n * p), n, p); X[, 1] = 1 
    beta = rnorm(p) * 0.05 
    eta = X %*% beta
    n_treat = floor(n / 2)
    w = sample(c(rep(1, n_treat), rep(0, n - n_treat)))
    eta = eta + 0.1 * w 
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

time_expr_ms = function(expr, times = B_TIME, env = parent.frame(), target_batch_ms = TARGET_BATCH_MS, max_inner_reps = MAX_INNER_REPS) {
    micro_time_ms = function() {
        mb_ms = microbenchmark(eval(expr, envir = env), times = FAST_PATH_MICROBENCH_REPS)$time / 1e6
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
    list(cls = "InferenceOrdinalAdjCatLogitRegr", pkg = "VGAM", func = "vglm+summary", expr = quote({
        res = VGAM::vglm(factor(y, ordered=T) ~ treatment + x1 + x2 + x3 + x4, VGAM::acat(), data=df)
        summary(res)@coef3[1, 4]
    })),
    list(cls = "InferenceOrdinalContRatioRegr", pkg = "VGAM", func = "vglm+summary", expr = quote({
        res = VGAM::vglm(factor(y, ordered=T) ~ treatment + x1 + x2 + x3 + x4, VGAM::cratio(), data=df)
        summary(res)@coef3[1, 4]
    })),
    list(cls = "InferenceSurvivalLogRank", pkg = "survival", func = "survdiff", expr = quote(survival::survdiff(survival::Surv(y, dead) ~ treatment, data = df)$pvalue)),
    list(cls = "InferenceSurvivalGehanWilcox", pkg = "survival", func = "survdiff(rho=1)", expr = quote(survival::survdiff(survival::Surv(y, dead) ~ treatment, data = df, rho = 1)$pvalue)),
    list(cls = "InferenceAllSimpleMeanDiffPooledVar", pkg = "stats", func = "t.test(pool)", expr = quote(t.test(df$y[df$treatment==1], df$y[df$treatment==0], var.equal = TRUE)$p.value)),
    list(cls = "InferenceAllSimpleWilcox", pkg = "stats", func = "wilcox.test", expr = quote(wilcox.test(df$y[df$treatment==1], df$y[df$treatment==0])$p.value)),
    list(cls = "InferenceIncidExactFisher", pkg = "stats", func = "fisher.test", expr = quote(fisher.test(table(df$treatment, df$y))$p.value)),
    list(cls = "InferenceIncidCMH", pkg = "stats", func = "mantelhaen", expr = quote(mantelhaen.test(table(df$treatment, df$y, df$g))$p.value)),
    list(cls = "InferenceOrdinalJonckheereTerpstraTest", pkg = "clinfun", func = "jonckheere", expr = quote(clinfun::jonckheere.test(df$y, df$treatment)$p.value)),
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
    list(cls = "InferenceSurvivalDepCensTransformRegr", pkg = "None", func = "None", expr = NULL),
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
    })),
    list(cls = "InferencePropZeroOneInflatedBetaRegr", pkg = "None", func = "None", expr = NULL)
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
    else if (grepl("Incid|Binomial|Wald|CMH|Fisher|Zhang|Robins|Newcombe|Nurminen", cls_name)) { resp_type = "incidence"; family = "logistic" }
    else if (grepl("Count|Poisson|NegBin|ZINB|ZAP|Hurdle", cls_name)) { resp_type = "count"; family = "poisson" }
    else if (grepl("Prop|Beta|ZOIB|Fractional", cls_name)) { resp_type = "proportion"; family = "beta" }
    else if (grepl("Survival|Cox|Weibull|KM|Rank|LogRank|Gehan|RMST|RMDiff|LWACox|Clayton", cls_name)) { resp_type = "survival"; family = "cox" }
    
    if (cls_name == "InferenceIncidLogBinomial") family = "log-binomial"
    
    n = N_WALD 
    d = generate_data(n = n, family = family)
    pval_delta = if (cls_name == "InferenceIncidGCompRiskRatio") 1 else 0
    
    # Timing EDI (Model Fit + Wald Test)
    t_edi = tryCatch({
        des = DesignFixediBCRD$new(n = n, response_type = resp_type)
        X_df = as.data.frame(d$X[,-1,drop=F])
        colnames(X_df) = paste0("x", 1:ncol(X_df))
        if (cls_name == "InferenceSurvivalStratCoxPHRegr") X_df$g = factor(rep(1:10, length.out = n))
        
        des$add_all_subjects_to_experiment(X_df)
        des$overwrite_all_subject_assignments(d$w)
        des$add_all_subject_responses(d$y, deads = if(exists("dead", d)) d$dead else NULL)
        
        inf_obj = safe_instantiate(cls_name, des)
        
        # Pre-cache method existence
        has_asymp = "compute_asymp_two_sided_pval" %in% names(inf_obj)
        has_exact = "compute_exact_two_sided_pval_for_treatment_effect" %in% names(inf_obj)
        priv_env = inf_obj$.__enclos_env__$private
        
        # Determine the benchmark expression
        bench_expr = if (cls_name == "InferenceAllSimpleWilcox") {
            quote({
                priv_env$hl_point_estimate(d$y, d$w)
            })
        } else {
            if (has_asymp) {
                quote({
                    priv_env$cached_mod = NULL
                    priv_env$cached_values = list()
                    inf_obj$compute_estimate(estimate_only = FALSE)
                    p = inf_obj$compute_asymp_two_sided_pval(delta = pval_delta)
                    if (!is.finite(p)) stop("Non-finite EDI p-value.")
                    p
                })
            } else if (has_exact) {
                quote({
                    priv_env$cached_mod = NULL
                    priv_env$cached_values = list()
                    inf_obj$compute_estimate(estimate_only = FALSE)
                    p = inf_obj$compute_exact_two_sided_pval_for_treatment_effect()
                    if (!is.finite(p)) stop("Non-finite EDI p-value.")
                    p
                })
            } else {
                quote({
                    priv_env$cached_mod = NULL
                    priv_env$cached_values = list()
                    inf_obj$compute_estimate(estimate_only = FALSE)
                })
            }
        }
        
        # Warmup
        eval(bench_expr)
        time_expr_ms(bench_expr)
    }, error = function(e) { cat("  EDI Error:", e$message, "\n"); NA_real_ })
    
    # Timing Canonical
    t_can = tryCatch({
        if (!is.null(spec$expr) && spec$pkg != "None") {
            library(spec$pkg, character.only = TRUE)
            X_cols = d$X[,-1,drop=F]; colnames(X_cols) = paste0("x", 1:ncol(X_cols))
            df = data.frame(y = d$y, treatment = d$w, g = factor(rep(1:10, length.out = n)), dead = if(exists("dead", d)) d$dead else 1)
            df = cbind(df, X_cols); X_can = cbind(`(Intercept)` = 1, treatment = d$w, as.matrix(X_cols))
            
            eval(spec$expr)
            time_expr_ms(spec$expr)
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

# Overwrite the previous Wald table
report_lines = readLines("package_metadata/benchmark_model_fits.md")
wald_start = grep("## Wald Test Performance", report_lines)
if (length(wald_start)) {
    report_lines = report_lines[1:(wald_start-1)]
}

header = c(
  "## Wald Test Performance (Full Inference)",
  "",
  "This table compares the performance of **Full Inference** (Model Fit + Standard Error calculation + P-value derivation).",
  "Unlike the point-estimation table above, these results include the computational cost of the variance-covariance matrix (Hessian or Fisher Information) and the Wald test statistic calculation.",
  "All paths (EDI and Canonical) use a reduced sample size ($N=200$) for this full-inference benchmark to ensure iterative stability.",
  "EDI timings in this table correspond to fixed `iBCRD` design objects.",
  "EDI regression models (Logistic, Poisson) are benchmarked using the **IRLS** optimizer for these Wald tests.",
  "**Note on Coverage**: `InferenceIncidExactZhang` does not appear in this Wald table because it is not applicable under `iBCRD`.",
  "**Note on Coverage**: `InferenceIncidCMH` is retained for coverage, but under `iBCRD` its EDI asymptotic p-value may be non-finite, in which case the row is reported as `NA`.",
  "**Note on Slowdowns**: For some non-parametric tests (e.g. Jonckheere-Terpstra), EDI computes an **exact** p-value while R's counterpart uses a normal approximation for $N=200$, leading to a speedup < 1x.",
  paste0("**Timing Note**: All timings are medians over ", B_TIME, " warmed runs measured with adaptive batched `system.time`; paths below ", FAST_PATH_THRESHOLD_MS, " ms use `microbenchmark(times = ", FAST_PATH_MICROBENCH_REPS, ")` instead."),
  "",
  "| Class | Response | EDI Time (ms) | Canonical Pkg | Canonical Func | Canonical Time (ms) | Speedup |",
  "| :--- | :--- | :--- | :--- | :--- | :--- | :--- |"
)
rows = paste0("| ", dt$Class, " | ", dt$Response, " | ", dt$EDI_Time_ms, " | ", 
              dt$Canonical_Pkg, " | ", dt$Canonical_Func, " | ", dt$Canonical_Time_ms, " | ", 
              dt$Speedup, " |")

writeLines(c(report_lines, "", header, rows), "package_metadata/benchmark_model_fits.md")
cat("Wald Test Benchmark complete.\n")
