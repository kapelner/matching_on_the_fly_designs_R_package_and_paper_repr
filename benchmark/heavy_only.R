
if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The 'pkgload' package is required.")
}
pkgload::load_all("/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI", quiet = TRUE)
library(data.table)

set.seed(42)

# SMALL SCALE for heavy paths
N_SIZE_H = 100; P_SIZE_H = 5
RAND_H = 5; BOOT_H = 5; JK_H = 3; PB_H = 5

heavy_classes = c(
    "InferenceContinKKGLMM", "InferenceContinKKOLSIVWC", "InferenceContinKKOLSOneLik", "InferenceContinKKQuantileRegrIVWC", "InferenceContinKKQuantileRegrOneLik",
    "InferenceContinKKRobustRegrIVWC", "InferenceContinKKRobustRegrOneLik", "InferenceCountKKCondPoissonOneLik",
    "InferenceCountKKGLMM", "InferenceCountKKHurdlePoissonIVWC", "InferenceCountKKHurdlePoissonOneLik", "InferenceCountPoissonKKGEE",
    "InferenceIncidKKCondLogitIVWC", "InferenceIncidKKCondLogitOneLik", "InferenceIncidKKGCompRiskDiff", "InferenceIncidKKGCompRiskRatio", "InferenceIncidKKGEE", "InferenceIncidKKModifiedPoisson",
    "InferenceIncidKKNewcombeRiskDiff", "InferenceOrdinalKKCLMM", "InferenceOrdinalKKCondAdjCatLogitRegr", "InferenceOrdinalKKGEE", "InferenceOrdinalKKGLMM",
    "InferencePropKKGEE", "InferencePropKKQuantileRegrIVWC", "InferencePropKKQuantileRegrOneLik",
    "InferenceSurvivalKKClaytonCopulaIVWC", "InferenceSurvivalKKClaytonCopulaOneLik", "InferenceSurvivalKKLWACoxIVWC", "InferenceSurvivalKKLWACoxOneLik", "InferenceSurvivalKKRankRegrIVWC",
    "InferenceSurvivalKKStratCoxOneLik", "InferenceSurvivalKKWeibullFrailtyIVWC", "InferenceSurvivalKKWeibullFrailtyOneLik"
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

results_csv = "/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/heavy_paths_only.csv"
if (file.exists(results_csv)) file.remove(results_csv)

for (cls_name in heavy_classes) {
    cat(sprintf("\n>>> %s... ", cls_name)); flush.console()
    tryCatch({
        meta = infer_metadata(cls_name)
        raw_data = generate_data(N_SIZE_H, P_SIZE_H, meta$rt)
        formula_str = paste("~", paste(paste0("x", 1:(P_SIZE_H-1)), collapse = " + "))
        
        # Use DesignFixedBinaryMatch to satisfy KK path requirements
        des_obj = DesignFixedBinaryMatch$new(response_type = meta$rt, n = N_SIZE_H, model_formula = as.formula(formula_str))
        ord_levels = if (meta$rt == "ordinal") as.character(1:4) else NULL
        if (!is.null(ord_levels)) des_obj$.__enclos_env__$private$ordinal_levels = ord_levels
        
        X_df = as.data.frame(raw_data$X[, -1, drop = FALSE]); colnames(X_df) = paste0("x", 1:(P_SIZE_H-1))
        des_obj$add_all_subjects_to_experiment(X_df); des_obj$assign_w_to_all_subjects()
        des_obj$add_all_subject_responses(raw_data$y, deads = raw_data$dead)
        
        inf_obj = get(cls_name)$new(des_obj)
        inf_cold = inf_obj$clone(deep = TRUE); inf_cold$.__enclos_env__$private$fit_warm_start_enabled = FALSE
        inf_warm = inf_obj$clone(deep = TRUE); inf_warm$.__enclos_env__$private$fit_warm_start_enabled = TRUE
        mock_ctx = list(n_units = N_SIZE_H, row_to_unit = 1:N_SIZE_H)
        inf_cold$.__enclos_env__$private$current_bayesian_bootstrap_context = mock_ctx
        inf_warm$.__enclos_env__$private$current_bayesian_bootstrap_context = mock_ctx

        t_rc = system.time({ tryCatch(inf_cold$compute_rand_two_sided_pval(r = RAND_H, show_progress = FALSE), error = function(e) NULL) })["elapsed"]
        t_rw = system.time({ tryCatch(inf_warm$compute_rand_two_sided_pval(r = RAND_H, show_progress = FALSE), error = function(e) NULL) })["elapsed"]
        t_bc = system.time({ tryCatch(inf_cold$compute_bootstrap_confidence_interval(B = BOOT_H, show_progress = FALSE), error = function(e) NULL) })["elapsed"]
        t_bw = system.time({ tryCatch(inf_warm$compute_bootstrap_confidence_interval(B = BOOT_H, show_progress = FALSE), error = function(e) NULL) })["elapsed"]
        t_jc = system.time({ tryCatch(for (k in 1:JK_H) { w_jk=rep(1,N_SIZE_H); w_jk[k]=0; inf_cold$compute_estimate_with_bootstrap_weights(w_jk, TRUE) }, error = function(e) NULL) })["elapsed"]
        tryCatch(inf_warm$compute_estimate(), error = function(e) NULL)
        t_jw = system.time({ tryCatch(for (k in 1:JK_H) { w_jk=rep(1,N_SIZE_H); w_jk[k]=0; inf_warm$compute_estimate_with_bootstrap_weights(w_jk, TRUE) }, error = function(e) NULL) })["elapsed"]
        
        t_pc = NA_real_; t_pw = NA_real_
        if (inf_obj$supports_lik_ratio_param_bootstrap()) {
            t_pc = system.time({ tryCatch(inf_cold$compute_lik_ratio_bootstrap_confidence_interval(B = PB_H, show_progress = FALSE), error = function(e) NULL) })["elapsed"]
            t_pw = system.time({ tryCatch(inf_warm$compute_lik_ratio_bootstrap_confidence_interval(B = PB_H, show_progress = FALSE), error = function(e) NULL) })["elapsed"]
        }
        res = data.table(Path = cls_name, Rand_C = t_rc, Rand_W = t_rw, Boot_C = t_bc, Boot_W = t_bw, JK_C = t_jc, JK_W = t_jw, PB_C = t_pc, PB_W = t_pw)
        fwrite(res, results_csv, append = file.exists(results_csv), col.names = !file.exists(results_csv))
        cat("Done."); flush.console()
    }, error = function(e) { cat(sprintf("Failed: %s ", conditionMessage(e))) })
}
