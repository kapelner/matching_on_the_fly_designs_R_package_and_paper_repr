
if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The 'pkgload' package is required.")
}
pkgload::load_all("/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI", quiet = TRUE)
library(data.table)

set.seed(42)

# HIGH LOAD Parameters for noise-free results
N_VAL = 500; P_VAL = 20
R_VAL = 10; B_VAL = 10; J_VAL = 5; PB_VAL = 10

all_classes = c(
    "InferenceAllSimpleMeanDiff", "InferenceAllSimpleMeanDiffPooledVar", "InferenceAllSimpleWilcox",
    "InferenceContinLin", "InferenceContinOLS", "InferenceContinQuantileRegr", "InferenceContinRobustRegr",
    "InferenceCountCompositeLikelihood", "InferenceCountHurdleNegBin", "InferenceCountHurdlePoisson", 
    "InferenceCountNegBin", "InferenceCountPoisson", 
    "InferenceCountQuasiPoisson", "InferenceCountRobustPoisson", "InferenceCountZeroInflatedNegBin", "InferenceCountZeroInflatedPoisson",
    "InferenceIncidBinomialIdentityRiskDiff", "InferenceIncidExactFisher", "InferenceIncidExactZhang", 
    "InferenceIncidGCompRiskDiff", "InferenceIncidGCompRiskRatio", 
    "InferenceIncidLogBinomial", "InferenceIncidLogRegr", "InferenceIncidModifiedPoisson",
    "InferenceIncidNewcombeRiskDiff", "InferenceIncidProbitRegr", "InferenceIncidRiskDiff",
    "InferenceOrdinalAdjCatLogitRegr", "InferenceOrdinalCauchitRegr", "InferenceOrdinalCloglogRegr", "InferenceOrdinalContRatioRegr", "InferenceOrdinalGCompMeanDiff",
    "InferenceOrdinalOrderedProbitRegr", "InferenceOrdinalRidit", "InferenceOrdinalStereotypeLogitRegr",
    "InferencePropBetaRegr", "InferencePropFractionalLogit", "InferencePropGCompMeanDiff", 
    "InferencePropZeroOneInflatedBetaRegr",
    "InferenceSurvivalCoxPHRegr", "InferenceSurvivalGehanWilcox", 
    "InferenceSurvivalKMDiff", "InferenceSurvivalLogRank", "InferenceSurvivalRestrictedMeanDiff", "InferenceSurvivalStratCoxPHRegr", "InferenceSurvivalWeibullRegr",
    "InferenceContinKKOLSIVWC", "InferenceContinKKOLSOneLik", "InferenceContinKKQuantileRegrIVWC", "InferenceContinKKQuantileRegrOneLik",
    "InferenceContinKKRobustRegrIVWC", "InferenceContinKKRobustRegrOneLik", "InferenceCountKKCondPoissonOneLik",
    "InferenceCountKKHurdlePoissonOneLik", "InferenceCountPoissonKKGEE",
    "InferenceIncidKKCondLogitIVWC", "InferenceIncidKKCondLogitOneLik", "InferenceIncidKKGCompRiskDiff", "InferenceIncidKKGCompRiskRatio", "InferenceIncidKKGEE", "InferenceIncidKKModifiedPoisson",
    "InferenceIncidKKNewcombeRiskDiff", "InferenceOrdinalKKCondAdjCatLogitRegr", "InferenceOrdinalKKGEE",
    "InferencePropKKGEE", "InferencePropKKQuantileRegrIVWC", "InferencePropKKQuantileRegrOneLik",
    "InferenceSurvivalKKClaytonCopulaIVWC", "InferenceSurvivalKKClaytonCopulaOneLik", "InferenceSurvivalKKRankRegrIVWC"
)

infer_metadata = function(class_name) {
    rt = "continuous"
    if (grepl("Incid|Incidence", class_name)) rt = "incidence"
    else if (grepl("Count", class_name)) rt = "count"
    else if (grepl("Prop|Beta", class_name)) rt = "proportion"
    else if (grepl("Survival|KM|Cox|Weibull|LogRank|RMST|RMDiff", class_name)) rt = "survival"
    else if (grepl("Ordinal|Ridit|Jonckheere", class_name)) rt = "ordinal"
    list(rt = rt)
}

generate_data = function(n, p, rt) {
    X = matrix(rnorm(n * p), n, p); X[, 1] = 1; beta = rnorm(p) * 0.1; eta = X %*% beta
    res = list(X = X, rt = rt, dead = rep(1, n))
    if (rt == "incidence") res$y = rbinom(n, 1, 1 / (1 + exp(-eta)))
    else if (rt == "count") res$y = rpois(n, exp(pmin(pmax(eta, -2), 2)))
    else if (rt == "proportion") res$y = pmax(pmin(rbeta(n, (1 / (1 + exp(-eta))) * 10, (1 - (1 / (1 + exp(-eta)))) * 10), 1 - 1e-6), 1e-6)
    else if (rt == "survival") { res$y = rexp(n, exp(pmin(pmax(eta, -3), 3))); res$dead = rbinom(n, 1, 0.8) }
    else if (rt == "ordinal") {
        p1 = 1 / (1 + exp(-(-1 - eta))); p_le_2 = 1 / (1 + exp(-(0 - eta))); p_le_3 = 1 / (1 + exp(-(1 - eta)))
        probs = cbind(p1, pmax(0, p_le_2 - p1), pmax(0, p_le_3 - p_le_2), pmax(0, 1 - p_le_3))
        res$y = apply(probs / rowSums(probs), 1, function(p) sample(1:4, 1, prob = p))
    } else { res$y = as.numeric(eta + rnorm(n)) }
    res
}

results_csv = "/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/comprehensive_noise_free_results.csv"
if (file.exists(results_csv)) file.remove(results_csv)

for (cls_name in all_classes) {
    cat(sprintf("\n>>> %s... ", cls_name)); flush.console()
    tryCatch({
        meta = infer_metadata(cls_name)
        raw_data = generate_data(N_VAL, P_VAL, meta$rt)
        formula_str = paste("~", paste(paste0("x", 1:(P_VAL-1)), collapse = " + "))
        des_obj = DesignFixedBinaryMatch$new(response_type = meta$rt, n = N_VAL, model_formula = as.formula(formula_str))
        ord_levels = if (meta$rt == "ordinal") as.character(1:4) else NULL
        if (!is.null(ord_levels)) des_obj$.__enclos_env__$private$ordinal_levels = ord_levels
        
        X_df = as.data.frame(raw_data$X[, -1, drop = FALSE]); colnames(X_df) = paste0("x", 1:(P_VAL-1))
        des_obj$add_all_subjects_to_experiment(X_df); des_obj$assign_w_to_all_subjects()
        des_obj$add_all_subject_responses(raw_data$y, deads = raw_data$dead)
        
        inf_obj = get(cls_name)$new(des_obj)
        inf_cold = inf_obj$clone(deep = TRUE); inf_cold$.__enclos_env__$private$fit_warm_start_enabled = FALSE
        inf_warm = inf_obj$clone(deep = TRUE); inf_warm$.__enclos_env__$private$fit_warm_start_enabled = TRUE
        
        # MLE Anchor
        mle_fit = tryCatch(inf_warm$compute_estimate(), error = function(e) NULL)
        if (!is.null(mle_fit)) {
            pw = inf_warm$.__enclos_env__$private
            if (is.list(mle_fit)) { pw$fit_warm_start = mle_fit$b %||% mle_fit$params; pw$fit_warm_start_fisher = mle_fit$fisher_information }
            else { pw$fit_warm_start = as.numeric(mle_fit) }
        }
        
        mock_ctx = list(n_units = N_VAL, row_to_unit = 1:N_VAL)
        inf_cold$.__enclos_env__$private$current_bayesian_bootstrap_context = mock_ctx
        inf_warm$.__enclos_env__$private$current_bayesian_bootstrap_context = mock_ctx

        t_rc = tryCatch(system.time({ inf_cold$compute_rand_two_sided_pval(r = R_VAL, show_progress = FALSE) })["elapsed"], error = function(e) NA_real_)
        t_rw = tryCatch(system.time({ inf_warm$compute_rand_two_sided_pval(r = R_VAL, show_progress = FALSE) })["elapsed"], error = function(e) NA_real_)
        
        t_bc = tryCatch(system.time({ inf_cold$compute_bootstrap_confidence_interval(B = B_VAL, show_progress = FALSE) })["elapsed"], error = function(e) NA_real_)
        t_bw = tryCatch(system.time({ inf_warm$compute_bootstrap_confidence_interval(B = B_VAL, show_progress = FALSE) })["elapsed"], error = function(e) NA_real_)
        
        t_jc = tryCatch(system.time({ for (k in 1:J_VAL) { w_jk=rep(1,N_VAL); w_jk[k]=0; inf_cold$compute_estimate_with_bootstrap_weights(w_jk, TRUE) } })["elapsed"], error = function(e) NA_real_)
        t_jw = tryCatch(system.time({ for (k in 1:J_VAL) { w_jk=rep(1,N_VAL); w_jk[k]=0; inf_warm$compute_estimate_with_bootstrap_weights(w_jk, TRUE) } })["elapsed"], error = function(e) NA_real_)
        
        t_pc = NA_real_; t_pw = NA_real_
        # CRITICAL: Use correct method check to avoid "non-function" crash
        has_pb = "supports_lik_ratio_param_bootstrap" %in% names(inf_obj)
        if (has_pb && inf_obj$supports_lik_ratio_param_bootstrap()) {
            t_pc = tryCatch(system.time({ inf_cold$compute_lik_ratio_bootstrap_two_sided_pval(B = PB_VAL, show_progress = FALSE) })["elapsed"], error = function(e) NA_real_)
            t_pw = tryCatch(system.time({ inf_warm$compute_lik_ratio_bootstrap_two_sided_pval(B = PB_VAL, show_progress = FALSE) })["elapsed"], error = function(e) NA_real_)
        }
        res = data.table(Path = cls_name, Rand_C = t_rc, Rand_W = t_rw, Boot_C = t_bc, Boot_W = t_bw, JK_C = t_jc, JK_W = t_jw, PB_C = t_pc, PB_W = t_pw)
        fwrite(res, results_csv, append = file.exists(results_csv), col.names = !file.exists(results_csv))
        cat("Done."); flush.console()
    }, error = function(e) { cat(sprintf("Failed: %s ", conditionMessage(e))) })
}
