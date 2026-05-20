
if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The 'pkgload' package is required.")
}
pkgload::load_all("/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI", quiet = TRUE)
library(data.table)

set.seed(42)

# BALANCED LOAD
N = 800; P = 30
R = 10; B = 10; J_count = 5; PB = 10

heavy_specs = list(
    list(name = "InferenceIncidLogRegr", rt = "incidence"),
    list(name = "InferenceCountPoisson", rt = "count"),
    list(name = "InferenceCountNegBin", rt = "count"),
    list(name = "InferencePropBetaRegr", rt = "proportion"),
    list(name = "InferenceSurvivalWeibullRegr", rt = "survival")
)

generate_data = function(n, p, rt) {
    X = matrix(rnorm(n * p), n, p); X[, 1] = 1; beta = rnorm(p) * 0.1; eta = X %*% beta
    res = list(X = X, rt = rt, dead = rep(1, n))
    if (rt == "incidence") res$y = rbinom(n, 1, 1 / (1 + exp(-eta)))
    else if (rt == "count") res$y = rpois(n, exp(pmin(pmax(eta, -2), 2)))
    else if (rt == "proportion") res$y = pmax(pmin(rbeta(n, (1 / (1 + exp(-eta))) * 10, (1 - (1 / (1 + exp(-eta)))) * 10), 1 - 1e-6), 1e-6)
    else if (rt == "survival") { res$y = rexp(n, exp(pmin(pmax(eta, -3), 3))); res$dead = rbinom(n, 1, 0.8) }
    res
}

results_csv = "/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/heavy_paths_results.csv"

for (spec in heavy_specs) {
    cat(sprintf("\n>>> %s... ", spec$name)); flush.console()
    tryCatch({
        raw_data = generate_data(N, P, spec$rt)
        formula_str = paste("~", paste(paste0("x", 1:(P-1)), collapse = " + "))
        des_obj = DesignFixedBernoulli$new(response_type = spec$rt, n = N, model_formula = as.formula(formula_str))
        X_df = as.data.frame(raw_data$X[, -1, drop = FALSE]); colnames(X_df) = paste0("x", 1:(P-1))
        des_obj$add_all_subjects_to_experiment(X_df); des_obj$assign_w_to_all_subjects()
        des_obj$add_all_subject_responses(raw_data$y, deads = raw_data$dead)
        
        inf_obj = get(spec$name)$new(des_obj)
        inf_cold = inf_obj$clone(deep = TRUE); inf_cold$.__enclos_env__$private$fit_warm_start_enabled = FALSE; inf_cold$.__enclos_env__$private$smart_cold_start_default = FALSE
        inf_warm = inf_obj$clone(deep = TRUE); inf_warm$.__enclos_env__$private$fit_warm_start_enabled = TRUE; inf_warm$.__enclos_env__$private$smart_cold_start_default = TRUE
        
        mock_ctx = list(n_units = N, row_to_unit = 1:N)
        inf_cold$.__enclos_env__$private$current_bayesian_bootstrap_context = mock_ctx
        inf_warm$.__enclos_env__$private$current_bayesian_bootstrap_context = mock_ctx

        cat("Rand... "); flush.console()
        t_rc = system.time({ inf_cold$compute_rand_two_sided_pval(r = R, show_progress = FALSE) })["elapsed"]
        t_rw = system.time({ inf_warm$compute_rand_two_sided_pval(r = R, show_progress = FALSE) })["elapsed"]
        
        cat("Boot... "); flush.console()
        t_bc = system.time({ inf_cold$compute_bootstrap_confidence_interval(B = B, show_progress = FALSE) })["elapsed"]
        t_bw = system.time({ inf_warm$compute_bootstrap_confidence_interval(B = B, show_progress = FALSE) })["elapsed"]
        
        cat("JK... "); flush.console()
        t_jc = system.time({ for (k in 1:J_count) { w_jk=rep(1,N); w_jk[k]=0; inf_cold$compute_estimate_with_bootstrap_weights(w_jk, TRUE) } })["elapsed"]
        t_jw = system.time({ inf_warm$compute_estimate(); for (k in 1:J_count) { w_jk=rep(1,N); w_jk[k]=0; inf_warm$compute_estimate_with_bootstrap_weights(w_jk, TRUE) } })["elapsed"]
        
        t_pc = NA_real_; t_pw = NA_real_
        if (inf_obj$supports_lik_ratio_param_bootstrap()) {
            cat("PB... "); flush.console()
            t_pc = system.time({ inf_cold$compute_lik_ratio_bootstrap_confidence_interval(B = PB, show_progress = FALSE) })["elapsed"]
            t_pw = system.time({ inf_warm$compute_lik_ratio_bootstrap_confidence_interval(B = PB, show_progress = FALSE) })["elapsed"]
        }
        res = data.table(Path = spec$name, Rand_C = t_rc, Rand_W = t_rw, Boot_C = t_bc, Boot_W = t_bw, JK_C = t_jc, JK_W = t_jw, PB_C = t_pc, PB_W = t_pw)
        fwrite(res, results_csv, append = file.exists(results_csv), col.names = !file.exists(results_csv))
        cat("Done."); flush.console()
    }, error = function(e) { cat(sprintf("Failed: %s ", conditionMessage(e))) })
}
