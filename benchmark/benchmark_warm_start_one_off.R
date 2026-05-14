
if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The 'pkgload' package is required to benchmark the current EDI source tree.")
}
pkgload::load_all("EDI", quiet = TRUE)
library(EDI)
library(microbenchmark)
library(data.table)

set.seed(42)

# Helper to generate data
generate_data = function(n = 1000, p = 10, family = "logistic", n_groups = 50) {
    X = matrix(rnorm(n * p), n, p)
    X[, 1] = 1 # Intercept
    beta = rnorm(p) * 0.5
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

# Final results storage
all_results = list()

run_comprehensive_warm_benchmark = function(family_name, fit_fun, data, alg = NULL, extra_args = list()) {
    cat(sprintf("\n--- Benchmarking Family: %s ---\n", family_name))
    
    # Base arguments
    args_base = list(X = data$X, y = data$y)
    if (tolower(family_name) %in% c("cox", "weibull", "weibull (aft)")) {
        args_base$dead = data$dead
    }
    if (grepl("glmm", tolower(family_name)) || tolower(family_name) == "lmm") {
        args_base$group_id = data$group_id
        if (family_name != "LMM") args_base$j_T = 2L
    }
    if (!is.null(alg)) args_base$optimization_alg = alg
    args_base = c(args_base, extra_args)
    
    # 1. Obtain "Truth" (MLE) for this dataset
    cat("Fitting to obtain MLE...\n")
    mle_fit = tryCatch(do.call(fit_fun, args_base), error = function(e) {
        cat("MLE fit failed:", conditionMessage(e), "\n")
        return(NULL)
    })
    
    if (is.null(mle_fit)) return(NULL)
    
    # Extract warm start components
    b = mle_fit$params %||% mle_fit$coefficients %||% mle_fit$b %||% mle_fit$beta
    if (family_name == "Beta" && is.null(mle_fit$params)) b = c(mle_fit$coefficients, log(mle_fit$phi))
    if (family_name == "NegBin" && is.null(mle_fit$params)) b = c(mle_fit$b, log(mle_fit$theta_hat))
    if (family_name == "Weibull (AFT)" && is.null(mle_fit$params)) b = c(mle_fit$coefficients, mle_fit$log_sigma)
    
    w = mle_fit$w %||% mle_fit$mu
    info = mle_fit$fisher_information %||% mle_fit$observed_information %||% mle_fit$information %||% mle_fit$XtWX %||% mle_fit$hess_fisher_info_matrix
    
    # Function to prepare arguments for a specific combination
    prep_args = function(use_beta = FALSE, use_weights = FALSE, use_info = FALSE) {
        a = args_base
        if (use_beta) {
            if ("start_beta" %in% names(formals(fit_fun))) {
                # Some functions want just beta, others want full params
                if (family_name %in% c("Beta", "NegBin", "Weibull (AFT)")) {
                    a$start_beta = mle_fit$b %||% mle_fit$coefficients
                } else {
                    a$start_beta = b
                }
            }
            if ("start_params" %in% names(formals(fit_fun))) a$start_params = b
            if ("start_par" %in% names(formals(fit_fun))) a$start_par = b
            if ("start" %in% names(formals(fit_fun))) a$start = b
        }
        if (use_weights && "warm_start_weights" %in% names(formals(fit_fun))) {
            a$warm_start_weights = w
        }
        if (use_info && "warm_start_fisher_info" %in% names(formals(fit_fun))) {
            a$warm_start_fisher_info = info
        }
        a
    }
    
    # Define the 8 combinations
    combos = list(
        "Cold"          = prep_args(F, F, F),
        "Beta"          = prep_args(T, F, F),
        "Weights"       = prep_args(F, T, F),
        "Info"          = prep_args(F, F, T),
        "Beta+Weights"  = prep_args(T, T, F),
        "Beta+Info"     = prep_args(T, F, T),
        "Weights+Info"  = prep_args(F, T, T),
        "Full"          = prep_args(T, T, T)
    )
    
    # Filter combos that actually differ from Cold (e.g. if function doesn't support weights)
    valid_combos = list()
    for (n in names(combos)) {
        if (n == "Cold" || length(setdiff(names(combos[[n]]), names(args_base))) > 0) {
            valid_combos[[n]] = combos[[n]]
        }
    }
    
    # Iteration counts
    cat("Checking iteration counts...\n")
    its = sapply(valid_combos, function(a) {
        res = tryCatch(do.call(fit_fun, a), error = function(e) NULL)
        if (is.null(res)) return(NA)
        res$iterations %||% res$niter %||% 1 # If already converged, might be 0 or 1
    })
    
    # Microbenchmark
    cat("Microbenchmarking...\n")
    bm_calls = lapply(names(valid_combos), function(n) {
        substitute(do.call(fit_fun, valid_combos[[n]]), list(n = n))
    })
    names(bm_calls) = names(valid_combos)
    
    # Using microbenchmark on the calls
    bm_res = microbenchmark(list = bm_calls, times = 20)
    
    print(bm_res)
    
    # Summary for this family
    summ = as.data.table(bm_res)
    summ = summ[, .(time_ms = mean(time) / 1e6), by = expr]
    summ[, iterations := its[as.character(expr)]]
    summ[, family := family_name]
    
    all_results[[family_name]] <<- summ
}

# Run across families
families = list(
    list(name = "Logistic", fun = fast_logistic_regression_with_var_cpp, family = "logistic"),
    list(name = "Poisson", fun = fast_poisson_regression_with_var_cpp, family = "poisson", alg = "irls"),
    list(name = "Cox", fun = fast_coxph_regression_cpp, family = "cox", alg = "newton_raphson"),
    list(name = "Beta", fun = fast_beta_regression_with_var_cpp, family = "beta", alg = "newton_raphson"),
    list(name = "NegBin", fun = fast_neg_bin_with_var_cpp, family = "negbin", alg = "newton_raphson"),
    list(name = "LogisticGLMM", fun = fast_logistic_glmm_cpp, family = "glmm_logistic"),
    list(name = "PoissonGLMM", fun = fast_poisson_glmm_cpp, family = "glmm_poisson"),
    list(name = "Weibull (AFT)", fun = fast_weibull_regression_cpp, family = "weibull"),
    list(name = "Ordinal", fun = fast_ordinal_regression_with_var_cpp, family = "ordinal")
)

for (f in families) {
    data = generate_data(family = f$family)
    run_comprehensive_warm_benchmark(f$name, f$fun, data, alg = f$alg)
}

# Combine and show final table
final_tab = rbindlist(all_results)
cat("\n\n### COMPREHENSIVE ONE-OFF WARM START BENEFIT SUMMARY ###\n")
cat("Values are mean time in ms. Iterations in parentheses.\n\n")

# Pivot table for better viewing
pivoted = dcast(final_tab, family ~ expr, value.var = c("time_ms", "iterations"))

# Sort columns for logical flow: Cold -> Beta -> Weights -> Info -> ... -> Full
col_order = c("family", 
              "time_ms_Cold", "iterations_Cold",
              "time_ms_Beta", "iterations_Beta",
              "time_ms_Weights", "iterations_Weights",
              "time_ms_Info", "iterations_Info",
              "time_ms_Beta+Weights", "iterations_Beta+Weights",
              "time_ms_Beta+Info", "iterations_Beta+Info",
              "time_ms_Weights+Info", "iterations_Weights+Info",
              "time_ms_Full", "iterations_Full")

# Subset only columns that exist
existing_cols = intersect(col_order, names(pivoted))
pivoted = pivoted[, ..existing_cols]

print(pivoted)

# Pretty Markdown-ish table for console
cat("\n| Family | Cold | Beta | Weights | Info | Full |\n")
cat("| :--- | :---: | :---: | :---: | :---: | :---: |\n")
for (i in 1:nrow(pivoted)) {
    r = pivoted[i]
    fmt = function(prefix) {
        t_col = paste0("time_ms_", prefix)
        i_col = paste0("iterations_", prefix)
        if (t_col %in% names(r) && !is.na(r[[t_col]])) {
            sprintf("%.2f (%d)", r[[t_col]], r[[i_col]])
        } else {
            "-"
        }
    }
    cat(sprintf("| %s | %s | %s | %s | %s | %s |\n",
                r$family, fmt("Cold"), fmt("Beta"), fmt("Weights"), fmt("Info"), fmt("Full")))
}
