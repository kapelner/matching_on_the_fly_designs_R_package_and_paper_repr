
if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The 'pkgload' package is required.")
}
pkgload::load_all("/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI", quiet = TRUE)
library(data.table)

set.seed(42)

# DEFINITIVE PARAMETERS
N_V = 1000; P_V = 30
R_V = 50; B_V = 50; J_V = 10; PB_V = 50

all_classes = c(
    "InferenceAllSimpleMeanDiff", "InferenceContinLin", "InferenceContinOLS",
    "InferenceCountPoisson", "InferenceCountNegBin", "InferencePropBetaRegr",
    "InferenceSurvivalWeibullRegr", "InferenceOrdinalPropOddsRegr",
    "InferenceCountHurdlePoisson", "InferenceCountHurdleNegBin", "InferenceCountRobustPoisson",
    "InferenceIncidLogRegr", "InferenceIncidProbitRegr", "InferencePropFractionalLogit",
    "InferenceSurvivalKMDiff", "InferenceSurvivalStratCoxPHRegr",
    "InferenceContinKKOLSOneLik", "InferenceCountKKCondPoissonOneLik", "InferenceIncidKKCondLogitOneLik"
)

results_csv = "final_high_precision_results.csv"
done_paths = if (file.exists(results_csv)) fread(results_csv)$Path else c()

for (cls_name in all_classes) {
    if (cls_name %in% done_paths) next
    cat(sprintf("\n>>> %s... ", cls_name)); flush.console()
    tryCatch({
        X = matrix(rnorm(N_V * P_V), N_V, P_V); X[, 1] = 1; beta = rnorm(P_V) * 0.1; eta = X %*% beta
        rt = if (grepl("Incid", cls_name)) "incidence" else if (grepl("Count", cls_name)) "count" else if (grepl("Prop", cls_name)) "proportion" else if (grepl("Survival", cls_name)) "survival" else if (grepl("Ordinal", cls_name)) "ordinal" else "continuous"
        y = if (rt == "incidence") rbinom(N_V, 1, plogis(eta)) 
            else if (rt == "count") rpois(N_V, exp(pmin(pmax(eta, -2), 2))) 
            else if (rt == "proportion") pmax(pmin(rbeta(N_V, 10, 10), 1-1e-6), 1e-6) 
            else if (rt == "survival") rexp(N_V, 1) 
            else if (rt == "ordinal") sample(1:4, N_V, replace=TRUE) 
            else as.numeric(eta + rnorm(N_V))
        
        formula_str = paste("~", paste(paste0("V", 2:P_V), collapse = " + "))
        des_obj = DesignFixedBinaryMatch$new(response_type = rt, n = N_V, model_formula = as.formula(formula_str))
        if (rt == "ordinal") des_obj$.__enclos_env__$private$ordinal_levels = as.character(1:4)
        X_df = as.data.frame(X[, -1, drop = FALSE]); colnames(X_df) = paste0("V", 2:P_V)
        des_obj$add_all_subjects_to_experiment(X_df); des_obj$assign_w_to_all_subjects()
        des_obj$add_all_subject_responses(y, dead = rep(1, N_V))
        
        inf_w = get(cls_name)$new(des_obj)
        inf_w$.__enclos_env__$private$fit_warm_start_enabled = TRUE
        mle = tryCatch(inf_w$compute_estimate(), error = function(e) NULL)
        if (!is.null(mle)) {
            pw = inf_w$.__enclos_env__$private
            if (is.list(mle)) {
                pw$fit_warm_start = mle$b %||% mle$params
                pw$fit_warm_start_fisher = mle$fisher_information
                pw$fit_warm_start_type = if (!is.null(mle$params)) "params" else "beta"
            } else { pw$fit_warm_start = as.numeric(mle) }
        }
        
        inf_c = inf_w$clone(deep = TRUE); inf_c$.__enclos_env__$private$fit_warm_start_enabled = FALSE
        mock_ctx = list(n_units = N_V, row_to_unit = 1:N_V)
        inf_c$.__enclos_env__$private$current_bayesian_bootstrap_context = mock_ctx
        inf_w$.__enclos_env__$private$current_bayesian_bootstrap_context = mock_ctx

        t_rc = system.time({ tryCatch(inf_c$compute_rand_two_sided_pval(r = R_V, show_progress = FALSE), error = function(e) NULL) })["elapsed"]
        t_rw = system.time({ tryCatch(inf_w$compute_rand_two_sided_pval(r = R_V, show_progress = FALSE), error = function(e) NULL) })["elapsed"]
        
        t_bc = system.time({ tryCatch(inf_c$compute_bootstrap_confidence_interval(B = B_V, show_progress = FALSE), error = function(e) NULL) })["elapsed"]
        t_bw = system.time({ tryCatch(inf_w$compute_bootstrap_confidence_interval(B = B_V, show_progress = FALSE), error = function(e) NULL) })["elapsed"]
        
        t_jc = system.time({ tryCatch(for (k in 1:J_V) { w_jk=rep(1,N_V); w_jk[k]=0; inf_c$compute_estimate_with_bootstrap_weights(w_jk, TRUE) }, error = function(e) NULL) })["elapsed"]
        t_jw = system.time({ tryCatch(for (k in 1:J_V) { w_jk=rep(1,N_V); w_jk[k]=0; inf_w$compute_estimate_with_bootstrap_weights(w_jk, TRUE) }, error = function(e) NULL) })["elapsed"]
        
        t_pc = NA_real_; t_pw = NA_real_
        is_pb_sup = tryCatch(inf_w$.__enclos_env__$private$supports_lik_ratio_param_bootstrap(), error = function(e) FALSE)
        if (isTRUE(is_pb_sup)) {
            t_pc = system.time({ tryCatch(inf_c$compute_lik_ratio_bootstrap_two_sided_pval(B = PB_V, show_progress = FALSE), error = function(e) NULL) })["elapsed"]
            t_pw = system.time({ tryCatch(inf_w$compute_lik_ratio_bootstrap_two_sided_pval(B = PB_V, show_progress = FALSE), error = function(e) NULL) })["elapsed"]
        }
        
        res = data.table(Path = cls_name, Rand_C = t_rc, Rand_W = t_rw, Boot_C = t_bc, Boot_W = t_bw, JK_C = t_jc, JK_W = t_jw, PB_C = t_pc, PB_W = t_pw)
        fwrite(res, results_csv, append = file.exists(results_csv), col.names = !file.exists(results_csv))
        cat("Done."); flush.console()
    }, error = function(e) { cat(sprintf("Failed: %s ", conditionMessage(e))) })
}
