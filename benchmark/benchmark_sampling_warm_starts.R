
if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The 'pkgload' package is required.")
}
pkgload::load_all("EDI", TRUE)
library(EDI)
library(data.table)

set.seed(42)

# LARGE SCALE Parameters for meaningful work (>1ms per fit)
N = 1000; P = 80
R = 1; B = 1; J_count = 1; PB = 1

all_classes = c(
    "InferenceAllSimpleMeanDiff", "InferenceAllSimpleMeanDiffPooledVar", "InferenceAllSimpleWilcox",
    "InferenceContinKKGLMM", "InferenceContinKKOLSIVWC", "InferenceContinKKOLSOneLik", "InferenceContinKKQuantileRegrIVWC", "InferenceContinKKQuantileRegrOneLik",
    "InferenceContinKKRobustRegrIVWC", "InferenceContinKKRobustRegrOneLik", "InferenceContinLin", "InferenceContinOLS", "InferenceContinQuantileRegr", "InferenceContinRobustRegr",
    "InferenceCountCompositeLikelihood", "InferenceCountHurdleNegBin", "InferenceCountHurdlePoisson", "InferenceCountKKCondPoissonOneLik",
    "InferenceCountKKGLMM", "InferenceCountKKHurdlePoissonIVWC", "InferenceCountKKHurdlePoissonOneLik", "InferenceCountNegBin", "InferenceCountPoisson", "InferenceCountPoissonKKGEE",
    "InferenceCountQuasiPoisson", "InferenceCountRobustPoisson", "InferenceCountZeroInflatedNegBin", "InferenceCountZeroInflatedPoisson",
    "InferenceIncidBinomialIdentityRiskDiff", "InferenceIncidExactFisher", "InferenceIncidExactZhang", 
    "InferenceIncidGCompRiskDiff", "InferenceIncidGCompRiskRatio", "InferenceIncidKKCondLogitIVWC", "InferenceIncidKKCondLogitOneLik", 
    "InferenceIncidKKGCompRiskDiff", "InferenceIncidKKGCompRiskRatio", "InferenceIncidKKGEE", "InferenceIncidKKModifiedPoisson",
    "InferenceIncidKKNewcombeRiskDiff", "InferenceIncidLogBinomial", "InferenceIncidLogRegr", "InferenceIncidModifiedPoisson",
    "InferenceIncidNewcombeRiskDiff", "InferenceIncidProbitRegr", "InferenceIncidRiskDiff", "InferenceIncidenceWald",
    "InferenceOrdinalAdjCatLogitRegr", "InferenceOrdinalCauchitRegr", "InferenceOrdinalCloglogRegr", "InferenceOrdinalContRatioRegr", "InferenceOrdinalGCompMeanDiff",
    "InferenceOrdinalJonckheereTerpstraTest", "InferenceOrdinalKKCLMM", "InferenceOrdinalKKCondAdjCatLogitRegr", "InferenceOrdinalKKGEE", "InferenceOrdinalKKGLMM",
    "InferenceOrdinalMultiStereotypeLogitRegr", "InferenceOrdinalOrderedProbitRegr", "InferenceOrdinalRidit", "InferenceOrdinalStereotypeLogitRegr",
    "InferencePropBetaRegr", "InferencePropFractionalLogit", "InferencePropGCompMeanDiff", "InferencePropKKGEE", "InferencePropKKQuantileRegrIVWC",
    "InferencePropKKQuantileRegrOneLik", "InferencePropZeroOneInflatedBetaRegr",
    "InferenceSurvivalCoxPHRegr", "InferenceSurvivalGehanWilcox", "InferenceSurvivalKKClaytonCopulaIVWC",
    "InferenceSurvivalKKClaytonCopulaOneLik", "InferenceSurvivalKKLWACoxIVWC", "InferenceSurvivalKKLWACoxOneLik", "InferenceSurvivalKKRankRegrIVWC",
    "InferenceSurvivalKKStratCoxOneLik", "InferenceSurvivalKKWeibullFrailtyIVWC", "InferenceSurvivalKKWeibullFrailtyOneLik",
    "InferenceSurvivalKMDiff", "InferenceSurvivalLogRank", "InferenceSurvivalRestrictedMeanDiff", "InferenceSurvivalStratCoxPHRegr", "InferenceSurvivalWeibullRegr"
)

infer_metadata = function(class_name) {
    rt = "continuous"
    if (grepl("Incid|Incidence", class_name)) rt = "incidence"
    else if (grepl("Count", class_name)) rt = "count"
    else if (grepl("Prop|Beta", class_name)) rt = "proportion"
    else if (grepl("Survival|KM|Cox|Weibull|LogRank|RMST|RMDiff", class_name)) rt = "survival"
    else if (grepl("Ordinal|Ridit|Jonckheere", class_name)) rt = "ordinal"
    design = "bernoulli"; if (grepl("KK|CMH|Robins|Zhang", class_name)) design = "kk"
    list(rt = rt, design = design)
}

generate_data = function(n, p, rt) {
    X = matrix(rnorm(n * p), n, p); X[, 1] = 1; beta = rnorm(p) * 0.1; eta = X %*% beta
    res = list(X = X, rt = rt, dead = rep(1, n))
    if (rt == "incidence") res$y = rbinom(n, 1, 1 / (1 + exp(-eta)))
    else if (rt == "count") res$y = rpois(n, exp(pmin(pmax(eta, -2), 2)))
    else if (rt == "proportion") res$y = pmax(pmin(rbeta(n, (1 / (1 + exp(-eta))) * 10, (1 - (1 / (1 + exp(-eta)))) * 10), 1 - 1e-6), 1e-6)
    else if (rt == "survival") { res$y = rexp(n, exp(pmin(pmax(eta, -2), 2))); res$dead = rbinom(n, 1, 0.8) }
    else if (rt == "ordinal") {
        p1 = 1 / (1 + exp(-(-1 - eta))); p_le_2 = 1 / (1 + exp(-(0 - eta))); p_le_3 = 1 / (1 + exp(-(1 - eta)))
        probs = cbind(p1, pmax(0, p_le_2 - p1), pmax(0, p_le_3 - p_le_2), pmax(0, 1 - p_le_3))
        res$y = apply(probs / rowSums(probs), 1, function(p) sample(1:4, 1, prob = p))
    } else { res$y = as.numeric(eta + rnorm(n)) }
    res
}

# Split paths into Normal and Heavy
heavy_patterns = "GLMM|GEE|Clayton|Frailty|CompositeLikelihood"
all_classes_sorted = c(
    all_classes[!grepl(heavy_patterns, all_classes)],
    all_classes[grepl(heavy_patterns, all_classes)]
)

results_csv = "truly_exhaustive_no_na_results.csv"
if (file.exists(results_csv)) {
    done_paths = fread(results_csv)$Path
} else {
    done_paths = c()
}

for (cls_name in all_classes_sorted) {
    if (cls_name %in% done_paths) next
    
    cat(sprintf("\n>>> %s... ", cls_name)); flush.console()
    tryCatch({
        meta = infer_metadata(cls_name)
        raw_data = generate_data(N, P, meta$rt)
        formula_str = paste("~", paste(paste0("x", 1:(P-1)), collapse = " + "))
        
        ord_levels = if (meta$rt == "ordinal") as.character(1:4) else NULL
        
        des_obj = if (meta$design == "kk") DesignSeqOneByOneKK14$new(response_type = meta$rt, n = N, model_formula = as.formula(formula_str))
                  else DesignFixedBernoulli$new(response_type = meta$rt, n = N, model_formula = as.formula(formula_str))
        if (!is.null(ord_levels)) des_obj$.__enclos_env__$private$ordinal_levels = ord_levels
        
        X_df = as.data.frame(raw_data$X[, -1, drop = FALSE]); colnames(X_df) = paste0("x", 1:(P-1))
        if (meta$design == "kk") {
            tryCatch({ for (i in 1:N) des_obj$add_one_subject_to_experiment_and_assign(X_df[i, , drop=FALSE]) }, 
                     error = function(e) { des_obj$overwrite_all_subject_assignments(rbinom(N, 1, 0.5)) })
        } else {
            des_obj$add_all_subjects_to_experiment(X_df); des_obj$assign_w_to_all_subjects()
        }
        des_obj$add_all_subject_responses(raw_data$y, deads = raw_data$dead)
        
        inf_class = tryCatch(get(cls_name, envir = as.environment("package:EDI")), error = function(e) tryCatch(get(cls_name), error = function(e2) NULL))
        if (is.null(inf_class)) { cat("Class not found "); next }
        inf_obj = tryCatch(inf_class$new(des_obj), error = function(e) { cat(sprintf("Init failed: %s ", conditionMessage(e))); NULL })
        if (is.null(inf_obj)) next
        
        inf_cold = inf_obj$clone(deep = TRUE); inf_cold$.__enclos_env__$private$fit_warm_start_enabled = FALSE; inf_cold$.__enclos_env__$private$smart_cold_start_default = FALSE
        inf_warm = inf_obj$clone(deep = TRUE); inf_warm$.__enclos_env__$private$fit_warm_start_enabled = TRUE; inf_warm$.__enclos_env__$private$smart_cold_start_default = TRUE

        # 1. Rand
        t_rc = system.time({ tryCatch(inf_cold$compute_rand_two_sided_pval(r = R, show_progress = FALSE), error = function(e) NULL) })["elapsed"]
        t_rw = system.time({ tryCatch(inf_warm$compute_rand_two_sided_pval(r = R, show_progress = FALSE), error = function(e) NULL) })["elapsed"]
        
        # 2. NP Boot
        t_bc = system.time({ tryCatch(inf_cold$compute_bootstrap_confidence_interval(B = B, show_progress = FALSE), error = function(e) NULL) })["elapsed"]
        t_bw = system.time({ tryCatch(inf_warm$compute_bootstrap_confidence_interval(B = B, show_progress = FALSE), error = function(e) NULL) })["elapsed"]

        # 3. JK
        t_jc = system.time({ tryCatch(for (i in 1:J_count) { w=rep(1,N); w[i]=0; inf_cold$compute_estimate_with_bootstrap_weights(w, TRUE) }, error = function(e) NULL) })["elapsed"]
        t_jw = system.time({ tryCatch({ inf_warm$compute_estimate(); for (i in 1:J_count) { w=rep(1,N); w[i]=0; inf_warm$compute_estimate_with_bootstrap_weights(w, TRUE) } }, error = function(e) NULL) })["elapsed"]

        # 4. Param Boot
        t_pc = 0; t_pw = 0
        if ("supports_lik_ratio_param_bootstrap" %in% names(inf_obj) && inf_obj$supports_lik_ratio_param_bootstrap()) {
            t_pc = system.time({ tryCatch(inf_cold$compute_lik_ratio_bootstrap_confidence_interval(B = PB, show_progress = FALSE), error = function(e) NULL) })["elapsed"]
            t_pw = system.time({ tryCatch(inf_warm$compute_lik_ratio_bootstrap_confidence_interval(B = PB, show_progress = FALSE), error = function(e) NULL) })["elapsed"]
        }

        res = data.table(
            Path = cls_name,
            Rand_C = t_rc, Rand_W = t_rw,
            Boot_C = t_bc, Boot_W = t_bw,
            JK_C = t_jc, JK_W = t_jw,
            PB_C = t_pc, PB_W = t_pw
        )
        fwrite(res, results_csv, append = file.exists(results_csv))
        cat("Done.")
    }, error = function(e) { cat(sprintf("Failed: %s ", conditionMessage(e))) })
}

cat("\n\n### EXHAUSTIVE BENCHMARK BATCH FINISHED ###\n")
