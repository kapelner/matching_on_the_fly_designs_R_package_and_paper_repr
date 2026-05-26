
if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The 'pkgload' package is required.")
}
pkgload::load_all("/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI", quiet = TRUE)
library(data.table)

set.seed(42)

# HIGH SIGNAL AUDIT Parameters
# Large enough N to make fitting dominate, but small enough to be fast.
N_A = 1200; P_A = 20
R_A = 100; B_A = 100; J_A = 30; PB_A = 50

# Target list: only those that showed 0.0% or noise in previous runs
classes = c(
    "InferenceIncidLogRegr", "InferenceCountPoisson", "InferenceCountNegBin",
    "InferencePropBetaRegr", "InferenceSurvivalWeibullRegr", "InferenceContinOLS",
    "InferenceContinLin", "InferenceCountHurdlePoisson", "InferenceIncidProbitRegr",
    "InferenceSurvivalKMDiff", "InferenceContinKKOLSOneLik", "InferenceCountKKCondPoissonOneLik",
    "InferenceIncidKKCondLogitOneLik", "InferencePropGCompMeanDiff", "InferenceOrdinalRidit",
    "InferenceAllSimpleMeanDiff", "InferenceAllSimpleWilcox"
)

results_csv = "high_precision_audit_results.csv"
if (file.exists(results_csv)) file.remove(results_csv)

for (cls_name in classes) {
    cat(sprintf("\n>>> AUDITING %s... ", cls_name)); flush.console()
    
    tryCatch({
        X = matrix(rnorm(N_A * P_A), N_A, P_A); X[, 1] = 1; beta = rnorm(P_A) * 0.1; eta = X %*% beta
        rt = if (grepl("Incid", cls_name)) "incidence" else if (grepl("Count", cls_name)) "count" else if (grepl("Prop", cls_name)) "proportion" else if (grepl("Survival", cls_name)) "survival" else if (grepl("Ordinal", cls_name)) "ordinal" else "continuous"
        
        y = if (rt == "incidence") rbinom(N_A, 1, 1 / (1 + exp(-eta))) 
            else if (rt == "count") rpois(N_A, exp(pmin(pmax(eta, -2), 2))) 
            else if (rt == "proportion") pmax(pmin(rbeta(N_A, 10, 10), 1-1e-6), 1e-6) 
            else if (rt == "survival") rexp(N_A, 1) 
            else if (rt == "ordinal") sample(1:4, N_A, replace=TRUE) 
            else as.numeric(eta + rnorm(N_A))
        
        formula_str = paste("~", paste(paste0("V", 2:P_A), collapse = " + "))
        # Use BinaryMatch for all to ensure IVWC/KK paths work
        des_obj = DesignFixedBinaryMatch$new(response_type = rt, n = N_A, model_formula = as.formula(formula_str))
        if (rt == "ordinal") des_obj$.__enclos_env__$private$ordinal_levels = as.character(1:4)
        
        X_df = as.data.frame(X[, -1, drop = FALSE]); colnames(X_df) = paste0("V", 2:P_A)
        des_obj$add_all_subjects_to_experiment(X_df); des_obj$assign_w_to_all_subjects()
        des_obj$add_all_subject_responses(y, deads = rep(1, N_A))
        
        inf_w = get(cls_name)$new(des_obj)
        inf_w$.__enclos_env__$private$fit_warm_start_enabled = TRUE
        
        # Prime for solution anchor
        mle = tryCatch(inf_w$compute_estimate(), error = function(e) NULL)
        if (!is.null(mle)) {
            pw = inf_w$.__enclos_env__$private
            if (is.list(mle)) {
                pw$fit_warm_start = mle$b %||% mle$params
                pw$fit_warm_start_fisher = mle$fisher_information
                pw$fit_warm_start_type = if (!is.null(mle$params)) "params" else "beta"
            } else {
                pw$fit_warm_start = as.numeric(mle)
            }
        }
        
        inf_c = inf_w$clone(deep = TRUE)
        inf_c$.__enclos_env__$private$fit_warm_start_enabled = FALSE
        
        mock_ctx = list(n_units = N_A, row_to_unit = 1:N_A)
        inf_c$.__enclos_env__$private$current_bayesian_bootstrap_context = mock_ctx
        inf_w$.__enclos_env__$private$current_bayesian_bootstrap_context = mock_ctx

        cat("R "); flush.console()
        t_rc = system.time({ inf_c$compute_rand_two_sided_pval(r = R_A, show_progress = FALSE) })["elapsed"]
        t_rw = system.time({ inf_w$compute_rand_two_sided_pval(r = R_A, show_progress = FALSE) })["elapsed"]
        
        cat("B "); flush.console()
        t_bc = system.time({ inf_c$compute_bootstrap_confidence_interval(B = B_A, show_progress = FALSE) })["elapsed"]
        t_bw = system.time({ inf_w$compute_bootstrap_confidence_interval(B = B_A, show_progress = FALSE) })["elapsed"]
        
        cat("J "); flush.console()
        t_jc = system.time({ for (k in 1:J_A) { w_jk=rep(1,N_A); w_jk[k]=0; inf_c$compute_estimate_with_bootstrap_weights(w_jk, TRUE) } })["elapsed"]
        t_jw = system.time({ for (k in 1:J_A) { w_jk=rep(1,N_A); w_jk[k]=0; inf_w$compute_estimate_with_bootstrap_weights(w_jk, TRUE) } })["elapsed"]
        
        t_pc = NA_real_; t_pw = NA_real_
        is_pb_sup = tryCatch(inf_w$.__enclos_env__$private$supports_lik_ratio_param_bootstrap(), error = function(e) FALSE)
        if (isTRUE(is_pb_sup)) {
            cat("P "); flush.console()
            t_pc = system.time({ tryCatch(inf_c$compute_lik_ratio_bootstrap_two_sided_pval(B = PB_A, show_progress = FALSE), error = function(e) NULL) })["elapsed"]
            t_pw = system.time({ tryCatch(inf_w$compute_lik_ratio_bootstrap_two_sided_pval(B = PB_A, show_progress = FALSE), error = function(e) NULL) })["elapsed"]
        }
        
        res = data.table(Path = cls_name, Rand_C = t_rc, Rand_W = t_rw, Boot_C = t_bc, Boot_W = t_bw, JK_C = t_jc, JK_W = t_jw, PB_C = t_pc, PB_W = t_pw)
        fwrite(res, results_csv, append = file.exists(results_csv), col.names = !file.exists(results_csv))
        cat("Done."); flush.console()
    }, error = function(e) { cat(sprintf("Failed: %s ", conditionMessage(e))) })
}
