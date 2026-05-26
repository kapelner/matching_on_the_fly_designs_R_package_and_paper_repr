
library(data.table)
all_classes = c(
    "InferenceIncidLogRegr", "InferenceCountPoisson", "InferenceCountNegBin",
    "InferencePropBetaRegr", "InferenceSurvivalWeibullRegr", "InferenceContinOLS",
    "InferenceContinLin", "InferenceContinRobustRegr", "InferenceContinQuantileRegr",
    "InferenceCountHurdlePoisson", "InferenceCountHurdleNegBin", "InferenceCountRobustPoisson",
    "InferenceIncidProbitRegr", "InferenceIncidModifiedPoisson", "InferenceIncidRiskDiff",
    "InferencePropFractionalLogit", "InferenceSurvivalKMDiff", "InferenceSurvivalStratCoxPHRegr",
    "InferenceContinKKOLSOneLik", "InferenceCountKKCondPoissonOneLik", "InferenceIncidKKCondLogitOneLik"
)

results_csv = "comprehensive_high_precision_results.csv"

for (cls in all_classes) {
    cat(sprintf("Running %s...\n", cls))
    cmd = sprintf("Rscript -e 'pkgload::load_all(\"EDI\"); N=800; P=30; R=50; B=50; J=10; PB=40; res_file=\"%s\"; cls=\"%s\"; # actual measurement code here... '", results_csv, cls)
    # Actually, I'll just write the full command into a temporary script for each class
    script_content = sprintf("
        pkgload::load_all(\"EDI\", quiet = TRUE)
        library(data.table)
        set.seed(42)
        N=800; P=30; R=50; B=50; J=10; PB=40
        cls_name = \"%s\"
        results_csv = \"%s\"
        
        X = matrix(rnorm(N * P), N, P); X[, 1] = 1; beta = rnorm(P) * 0.1; eta = X %%*%% beta
        rt = if (grepl(\"Incid\", cls_name)) \"incidence\" else if (grepl(\"Count\", cls_name)) \"count\" else if (grepl(\"Prop\", cls_name)) \"proportion\" else if (grepl(\"Survival\", cls_name)) \"survival\" else if (grepl(\"Ordinal\", cls_name)) \"ordinal\" else \"continuous\"
        y = if (rt == \"incidence\") rbinom(N, 1, 1 / (1 + exp(-eta))) else if (rt == \"count\") rpois(N, exp(pmin(pmax(eta, -2), 2))) else if (rt == \"proportion\") pmax(pmin(rbeta(N, 1, 1), 1-1e-6), 1e-6) else if (rt == \"survival\") rexp(N, 1) else if (rt == \"ordinal\") sample(1:4, N, replace=TRUE) else as.numeric(eta + rnorm(N))
        
        formula_str = paste(\"~\", paste(paste0(\"V\", 2:P), collapse = \" + \"))
        des_obj = DesignFixedBinaryMatch$new(response_type = rt, n = N, model_formula = as.formula(formula_str))
        X_df = as.data.frame(X[, -1, drop = FALSE]); colnames(X_df) = paste0(\"V\", 2:P)
        des_obj$add_all_subjects_to_experiment(X_df); des_obj$assign_w_to_all_subjects()
        des_obj$add_all_subject_responses(y, deads = rep(1, N))
        
        inf_w = get(cls_name)$new(des_obj)
        inf_w$.__enclos_env__$private$fit_warm_start_enabled = TRUE
        mle = tryCatch(inf_w$compute_estimate(), error = function(e) NULL)
        if (!is.null(mle)) {
            pw = inf_w$.__enclos_env__$private
            if (is.list(mle)) {
                pw$fit_warm_start = mle$b %%||%% mle$params
                pw$fit_warm_start_fisher = mle$fisher_information
            } else { pw$fit_warm_start = as.numeric(mle) }
        }
        
        inf_c = inf_w$clone(deep = TRUE); inf_c$.__enclos_env__$private$fit_warm_start_enabled = FALSE
        mock_ctx = list(n_units = N, row_to_unit = 1:N)
        inf_c$.__enclos_env__$private$current_bayesian_bootstrap_context = mock_ctx
        inf_w$.__enclos_env__$private$current_bayesian_bootstrap_context = mock_ctx

        t_rc = system.time({ tryCatch(inf_c$compute_rand_two_sided_pval(r = R, show_progress = FALSE), error = function(e) NULL) })[\"elapsed\"]
        t_rw = system.time({ tryCatch(inf_w$compute_rand_two_sided_pval(r = R, show_progress = FALSE), error = function(e) NULL) })[\"elapsed\"]
        t_bc = system.time({ tryCatch(inf_c$compute_bootstrap_confidence_interval(B = B, show_progress = FALSE), error = function(e) NULL) })[\"elapsed\"]
        t_bw = system.time({ tryCatch(inf_w$compute_bootstrap_confidence_interval(B = B, show_progress = FALSE), error = function(e) NULL) })[\"elapsed\"]
        t_jc = system.time({ tryCatch(for (k in 1:J) { w_jk=rep(1,N); w_jk[k]=0; inf_c$compute_estimate_with_bootstrap_weights(w_jk, TRUE) }, error = function(e) NULL) })[\"elapsed\"]
        t_jw = system.time({ tryCatch(for (k in 1:J) { w_jk=rep(1,N); w_jk[k]=0; inf_w$compute_estimate_with_bootstrap_weights(w_jk, TRUE) }, error = function(e) NULL) })[\"elapsed\"]
        
        t_pc = NA_real_; t_pw = NA_real_
        is_pb_sup = tryCatch(inf_w$.__enclos_env__$private$supports_lik_ratio_param_bootstrap(), error = function(e) FALSE)
        if (isTRUE(is_pb_sup)) {
            t_pc = system.time({ tryCatch(inf_c$compute_lik_ratio_bootstrap_two_sided_pval(B = PB, show_progress = FALSE), error = function(e) NULL) })[\"elapsed\"]
            t_pw = system.time({ tryCatch(inf_w$compute_lik_ratio_bootstrap_two_sided_pval(B = PB, show_progress = FALSE), error = function(e) NULL) })[\"elapsed\"]
        }
        
        res = data.table(Path = cls_name, Rand_C = t_rc, Rand_W = t_rw, Boot_C = t_bc, Boot_W = t_bw, JK_C = t_jc, JK_W = t_jw, PB_C = t_pc, PB_W = t_pw)
        fwrite(res, results_csv, append = file.exists(results_csv), col.names = !file.exists(results_csv))
    ", cls, results_csv)
    writeLines(script_content, "worker.R")
    system("Rscript worker.R")
}
