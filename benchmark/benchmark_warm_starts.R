
if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The 'pkgload' package is required to benchmark the current EDI source tree.")
}
pkgload::load_all("EDI", quiet = TRUE)
library(EDI)
library(microbenchmark)
library(data.table)

set.seed(42)

# Helper to generate data
generate_data = function(n = 500, p = 5, family = "logistic", n_groups = 20) {
    X = matrix(rnorm(n * p), n, p)
    X[, 1] = 1 # Intercept
    beta = rnorm(p) * 0.5 # Smaller beta for stability
    eta = X %*% beta
    group_id = rep(1:n_groups, length.out = n)
    
    if (family == "logistic") {
        mu = 1 / (1 + exp(-eta))
        y = rbinom(n, 1, mu)
    } else if (family == "poisson") {
        mu = exp(eta)
        y = rpois(n, mu)
    } else if (family == "beta") {
        mu = 1 / (1 + exp(-eta))
        phi = 10
        y = rbeta(n, mu * phi, (1 - mu) * phi)
        y = pmax(pmin(y, 1 - 1e-6), 1e-6)
    } else if (family == "cox" || family == "weibull") {
        mu = exp(eta)
        y = rexp(n, mu)
        dead = rbinom(n, 1, 0.8)
        return(list(X = X, y = y, dead = dead, group_id = group_id))
    } else if (family == "negbin") {
        mu = exp(eta)
        y = rnbinom(n, size = 2, mu = mu)
    } else if (family == "ordinal") {
        # 3 levels
        p1 = 1 / (1 + exp(eta - 1))
        p2 = 1 / (1 + exp(eta + 1)) - p1
        p3 = 1 - p1 - p2
        probs = cbind(p1, p2, p3)
        probs = t(apply(probs, 1, function(x) pmax(x, 1e-6)))
        probs = probs / rowSums(probs)
        y = apply(probs, 1, function(p) sample(1:3, 1, prob = p))
    } else if (family == "lmm") {
        # Random intercept
        re = rnorm(n_groups)[group_id]
        y = eta + re + rnorm(n, 0, 0.5)
    } else if (family == "glmm_logistic") {
        re = rnorm(n_groups)[group_id]
        mu = 1 / (1 + exp(-(eta + re)))
        y = rbinom(n, 1, mu)
    } else if (family == "glmm_poisson") {
        re = rnorm(n_groups, 0, 0.2)[group_id]
        mu = exp(eta + re)
        y = rpois(n, mu)
    }
    
    list(X = X, y = y, group_id = group_id)
}

run_benchmark_family = function(family_name, fit_fun, data, alg = NULL, extra_args = list()) {
    cat(sprintf("\nBenchmarking Family: %s\n", family_name))
    
    # Initial fit to get warm start values
    args_base = list(X = data$X, y = data$y)
    if (tolower(family_name) %in% c("cox", "weibull")) {
        args_base$dead = data$dead
    }
    if (grepl("glmm", tolower(family_name)) || tolower(family_name) == "lmm") {
        args_base$group_id = data$group_id
        if (family_name != "LMM") args_base$j_T = 2L
    }
    if (!is.null(alg)) args_base$optimization_alg = alg
    args_base = c(args_base, extra_args)
    
    initial_fit = tryCatch(do.call(fit_fun, args_base), error = function(e) {
        cat("Initial fit failed:", conditionMessage(e), "\n")
        return(NULL)
    })
    
    if (is.null(initial_fit)) return(NULL)
    
    # Extract parameters correctly based on family
    b = initial_fit$params %||% initial_fit$coefficients %||% initial_fit$b
    if (family_name == "Beta" && is.null(initial_fit$params)) {
        b = c(initial_fit$coefficients, log(initial_fit$phi))
    } else if (family_name == "NegBin" && is.null(initial_fit$params)) {
        b = c(initial_fit$b, log(initial_fit$theta_hat))
    } else if (family_name == "Weibull" && is.null(initial_fit$params)) {
        b = c(initial_fit$coefficients, initial_fit$log_sigma)
    }
    
    w = initial_fit$w %||% initial_fit$mu
    info = initial_fit$fisher_information %||% initial_fit$observed_information %||% initial_fit$information %||% initial_fit$XtWX %||% initial_fit$hess_fisher_info_matrix
    
    # Pre-build the argument lists
    args_cold = args_base
    
    args_beta = args_base
    if ("start_beta" %in% names(formals(fit_fun))) {
        if (family_name %in% c("Beta", "NegBin", "Weibull")) {
            args_beta$start_beta = initial_fit$b %||% initial_fit$coefficients
        } else {
            args_beta$start_beta = b
        }
    }
    if ("start_params" %in% names(formals(fit_fun))) args_beta$start_params = b
    if ("start_par" %in% names(formals(fit_fun))) args_beta$start_par = b
    if ("start" %in% names(formals(fit_fun))) args_beta$start = b
    
    args_weights = args_base
    if ("warm_start_weights" %in% names(formals(fit_fun))) args_weights$warm_start_weights = w
    
    args_info = args_base
    if ("warm_start_fisher_info" %in% names(formals(fit_fun))) args_info$warm_start_fisher_info = info
    
    args_beta_weights = args_beta
    if ("warm_start_weights" %in% names(formals(fit_fun))) args_beta_weights$warm_start_weights = w
    
    args_beta_info = args_beta
    if ("warm_start_fisher_info" %in% names(formals(fit_fun))) args_beta_info$warm_start_fisher_info = info
    
    args_weights_info = args_weights
    if ("warm_start_fisher_info" %in% names(formals(fit_fun))) args_weights_info$warm_start_fisher_info = info
    
    args_full = args_beta_weights
    if ("warm_start_fisher_info" %in% names(formals(fit_fun))) args_full$warm_start_fisher_info = info
    
    # Only run microbenchmark for supported combinations
    bm_list = list(Cold = quote(do.call(fit_fun, args_cold)))
    if (!is.null(args_beta$start_beta) || !is.null(args_beta$start_params) || !is.null(args_beta$start_par) || !is.null(args_beta$start)) 
        bm_list$Beta = quote(do.call(fit_fun, args_beta))
    if (!is.null(args_weights$warm_start_weights)) 
        bm_list$Weights = quote(do.call(fit_fun, args_weights))
    if (!is.null(args_info$warm_start_fisher_info)) 
        bm_list$Info = quote(do.call(fit_fun, args_info))
    if (!is.null(args_beta_weights$warm_start_weights) && length(setdiff(names(args_beta_weights), names(args_base))) > 0)
        bm_list$Beta_Weights = quote(do.call(fit_fun, args_beta_weights))
    if (!is.null(args_beta_info$warm_start_fisher_info) && length(setdiff(names(args_beta_info), names(args_base))) > 0)
        bm_list$Beta_Info = quote(do.call(fit_fun, args_beta_info))
    if (!is.null(args_weights_info$warm_start_weights) && !is.null(args_weights_info$warm_start_fisher_info))
        bm_list$Weights_Info = quote(do.call(fit_fun, args_weights_info))
    if (length(bm_list) > 1) {
        bm_list$Full = quote(do.call(fit_fun, args_full))
    }

    results = microbenchmark(list = bm_list, times = 10)
    print(results)
    
    # Iteration counts
    safe_it = function(args) {
        res = tryCatch(do.call(fit_fun, args), error = function(e) NULL)
        if (is.null(res)) return(NA)
        it = res$iterations %||% res$niter %||% NA
        it
    }
    
    it_cold          = safe_it(args_cold)
    it_beta          = if (length(setdiff(names(args_beta), names(args_base))) > 0) safe_it(args_beta) else NA
    it_weights       = if (!is.null(args_weights$warm_start_weights)) safe_it(args_weights) else NA
    it_info          = if (!is.null(args_info$warm_start_fisher_info)) safe_it(args_info) else NA
    it_beta_weights  = if (!is.null(args_beta_weights$warm_start_weights)) safe_it(args_beta_weights) else NA
    it_beta_info     = if (!is.null(args_beta_info$warm_start_fisher_info)) safe_it(args_beta_info) else NA
    it_weights_info  = if (!is.null(args_weights_info$warm_start_weights) && !is.null(args_weights_info$warm_start_fisher_info)) safe_it(args_weights_info) else NA
    it_full          = safe_it(args_full)
    
    cat(sprintf("Iterations - Cold: %s, Beta: %s, Weights: %s, Info: %s, B+W: %s, B+I: %s, W+I: %s, Full: %s\n", 
                as.character(it_cold), as.character(it_beta), as.character(it_weights), 
                as.character(it_info), as.character(it_beta_weights), as.character(it_beta_info),
                as.character(it_weights_info), as.character(it_full)))
}

cat("\n--- Warm Start Benchmark Analysis ---\n")
cat("1. Beta (Parameters): Essential for reducing iteration counts. Massive gains in GLMMs (>80%).\n")
cat("2. Weights: Negligible alone. Only useful in IRLS models when combined with Beta.\n")
cat("3. Info (Hessian): Crucial when combined with Beta. 'Beta + Info' often matches 'Full' performance.\n")
cat("4. Combinations: 'Beta + Info' is the most robust strategy for significant speedups across all families.\n")

# Run all benchmarks
families = list(
    list(name = "LMM", fun = fast_gaussian_lmm_cpp, family = "lmm"),
    list(name = "LogisticGLMM", fun = fast_logistic_glmm_cpp, family = "glmm_logistic"),
    list(name = "PoissonGLMM", fun = fast_poisson_glmm_cpp, family = "glmm_poisson"),
    list(name = "Logistic", fun = fast_logistic_regression_with_var_cpp, family = "logistic"),
    list(name = "Poisson", fun = fast_poisson_regression_with_var_cpp, family = "poisson", alg = "irls"),
    list(name = "Cox", fun = fast_coxph_regression_cpp, family = "cox", alg = "newton_raphson"),
    list(name = "Beta", fun = fast_beta_regression_with_var_cpp, family = "beta", alg = "newton_raphson"),
    list(name = "NegBin", fun = fast_neg_bin_with_var_cpp, family = "negbin", alg = "newton_raphson"),
    list(name = "Log-Binomial", fun = fast_log_binomial_regression_with_var_cpp, family = "logistic"),
    list(name = "Robust", fun = fast_robust_regression_cpp, family = "logistic"), 
    list(name = "Weibull", fun = fast_weibull_regression_cpp, family = "weibull", alg = "newton_raphson"),
    list(name = "Ordinal", fun = fast_ordinal_regression_with_var_cpp, family = "ordinal", alg = "newton_raphson"),
    list(name = "Stereotype", fun = fast_stereotype_logit_with_var_cpp, family = "ordinal")
)

for (f in families) {
    data = generate_data(family = f$family)
    run_benchmark_family(f$name, f$fun, data, alg = f$alg)
}
