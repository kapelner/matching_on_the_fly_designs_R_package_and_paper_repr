
if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The 'pkgload' package is required to benchmark the current EDI source tree.")
}
pkgload::load_all("EDI", quiet = TRUE)
library(EDI)
library(microbenchmark)
library(data.table)

set.seed(42)

# Helper to generate data using the simulation framework's logic
generate_benchmark_data = function(n = 500, p = 5, response_type = "continuous") {
    # Reuse SimulationFramework's data generation
    # X from stats::rnorm, cond_exp_func_model="linear"
    # beta_x = seq(1, -1, length.out = p), rescaled to norm_sq_beta_vec=1
    
    # Actually, let's just use the exported helpers
    raw = generate_covariate_dataset(n = n, p = p, cond_exp_func_model = "linear")
    X = as.matrix(raw$X)
    y_cont = raw$y_cont
    
    # Apply treatment effect (betaT = 1) to treated (w)
    w = rep(0:1, length.out = n)
    # betaT = 1 is additive for most, log-mult for count/survival
    y_linear_model = y_cont
    if (response_type %in% c("count", "survival")) {
        y_linear_model = y_linear_model + w * 1.0 # log scale
    } else {
        y_linear_model = y_linear_model + w * 1.0 # linear/logit scale
    }
    
    y = transform_cont_y_based_on_response_type(y_linear_model, response_type)
    
    # Survival needs dead indicator
    dead = NULL
    if (response_type == "survival") {
        dead = rbinom(n, 1, 0.75)
    }
    
    # GLMMs need group_id
    group_id = rep(1:25, length.out = n)
    
    list(X = X, y = y, w = w, dead = dead, group_id = group_id)
}

run_comp = function(name, fit_fun, data, alg = NULL, extra_args = list()) {
    cat(sprintf("\n--- %s ---\n", name))
    
    args_base = list(X = data$X, y = data$y)
    if (!is.null(data$dead)) args_base$dead = data$dead
    if (!is.null(alg)) args_base$optimization_alg = alg
    if (grepl("GLMM", name) || name == "LMM") {
        args_base$group_id = data$group_id
        if (name != "LMM") args_base$j_T = 1L
    }
    if (name == "ZINB") args_base$Xzi = data$X
    if (name == "Zero-One Beta") args_base$X_zero_one = data$X
    args_base = c(args_base, extra_args)
    
    args_naive = c(args_base, list(smart_cold_start = FALSE))
    args_smart = c(args_base, list(smart_cold_start = TRUE))
    
    # Initial runs to capture iterations and verify convergence
    res_naive = tryCatch(do.call(fit_fun, args_naive), error = function(e) { cat("Naive failed:", e$message, "\n"); NULL })
    res_smart = tryCatch(do.call(fit_fun, args_smart), error = function(e) { cat("Smart failed:", e$message, "\n"); NULL })
    
    it_naive = if (!is.null(res_naive)) (res_naive$iterations %||% res_naive$niter %||% NA) else NA
    it_smart = if (!is.null(res_smart)) (res_smart$iterations %||% res_smart$niter %||% NA) else NA
    
    cat(sprintf("Iterations: Naive=%s, Smart=%s (Improvement: %.1f%%)\n", 
                it_naive, it_smart, (it_naive - it_smart)/it_naive * 100))

    bm = microbenchmark(
        Naive = do.call(fit_fun, args_naive),
        Smart = do.call(fit_fun, args_smart),
        times = 3
    )
    print(bm)
    
    # Extract median times in ms
    summary_bm = summary(bm)
    t_naive = summary_bm$median[summary_bm$expr == "Naive"]
    t_smart = summary_bm$median[summary_bm$expr == "Smart"]
    
    list(
        Model = name,
        It_Naive = it_naive,
        It_Smart = it_smart,
        It_Improv = sprintf("%.1f%%", (it_naive - it_smart)/it_naive * 100),
        Time_Naive_ms = round(t_naive, 2),
        Time_Smart_ms = round(t_smart, 2),
        Time_Improv = sprintf("%.1f%%", (t_naive - t_smart)/t_naive * 100)
    )
}

# Mapping benchmark items
configs = list(
    list(name = "Logistic", fun = fast_logistic_regression_with_var_cpp, rt = "incidence"),
    list(name = "Poisson", fun = fast_poisson_regression_with_var_cpp, rt = "count", alg = "irls"),
    list(name = "Weibull", fun = fast_weibull_regression_cpp, rt = "survival", alg = "newton_raphson"),
    list(name = "NegBin", fun = fast_neg_bin_with_var_cpp, rt = "count", alg = "newton_raphson"),
    list(name = "Beta", fun = fast_beta_regression_with_var_cpp, rt = "proportion", alg = "newton_raphson"),
    list(name = "Ordinal", fun = fast_ordinal_regression_with_var_cpp, rt = "ordinal", alg = "newton_raphson"),
    list(name = "Stereotype", fun = fast_stereotype_logit_with_var_cpp, rt = "ordinal"),
    list(name = "ZINB", fun = fast_zinb_cpp, rt = "count"),
    # list(name = "Zero-One Beta", fun = fast_zero_one_inflated_beta_cpp, rt = "proportion"),
    list(name = "LogisticGLMM", fun = fast_logistic_glmm_cpp, rt = "incidence"),
    list(name = "PoissonGLMM", fun = fast_poisson_glmm_cpp, rt = "count")
)

results_all = list()
for (cfg in configs) {
    cat(sprintf("Processing config: %s\n", cfg$name))
    data = generate_benchmark_data(n = 100, p = 7, response_type = cfg$rt)
    results_all[[cfg$name]] = run_comp(cfg$name, cfg$fun, data, alg = cfg$alg)
}

final_dt = rbindlist(results_all)
print(final_dt)

saveRDS(final_dt, "benchmark_cold_starts_results.RData")
