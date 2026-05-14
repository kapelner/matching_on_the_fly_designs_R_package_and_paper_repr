
if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The 'pkgload' package is required to benchmark the current EDI source tree.")
}
pkgload::load_all("EDI", quiet = TRUE)
library(EDI)
library(microbenchmark)
library(data.table)

set.seed(42)

# Helper to generate data (Smaller for demo)
generate_data = function(n = 250, p = 3, family = "logistic", n_groups = 10) {
    X = matrix(rnorm(n * p), n, p)
    X[, 1] = 1 
    beta = rnorm(p) * 0.5
    eta = X %*% beta
    group_id = rep(1:n_groups, length.out = n)
    
    if (family == "logistic") {
        mu = 1 / (1 + exp(-eta))
        y = rbinom(n, 1, mu)
    } else if (family == "poisson") {
        mu = exp(eta)
        y = rpois(n, mu)
    } else if (family == "glmm_logistic") {
        re = rnorm(n_groups)[group_id]
        mu = 1 / (1 + exp(-(eta + re)))
        y = rbinom(n, 1, mu)
    } else if (family == "negbin") {
        mu = exp(eta)
        y = rnbinom(n, size = 2, mu = mu)
    } else if (family == "weibull") {
        mu = exp(eta)
        y = rexp(n, 1/mu) 
        dead = rbinom(n, 1, 0.8)
        return(list(X = X, y = y, dead = dead, group_id = group_id))
    }
    
    list(X = X, y = y, group_id = group_id)
}

# Heuristic starts
get_reasonable_starts = function(data, family_name) {
    X = data$X
    y = data$y
    n = nrow(X)
    p = ncol(X)
    
    b_ols = NULL; w = NULL; info = NULL
    
    if (family_name == "Logistic (IRLS)" || family_name == "LogisticGLMM") {
        p_avg = mean(y)
        mu = rep(p_avg, n)
        w = mu * (1 - mu)
        mu_init = (y + 0.5) / 2
        eta_init = qlogis(mu_init)
        b_ols = fast_ols_cpp(X, eta_init)$coefficients
        info = t(X) %*% (X * w)
        if (family_name == "LogisticGLMM") b_ols = c(b_ols, 0) 
    } else if (family_name == "Poisson (IRLS)") {
        mu = y + 0.1
        eta = log(mu)
        b_ols = fast_ols_cpp(X, eta)$coefficients
        w = mu
        info = t(X) %*% (X * w)
    } else if (family_name == "NegBin") {
        mu = y + 0.1; eta = log(mu); theta = 1.0 
        b_ols = c(fast_ols_cpp(X, eta)$coefficients, log(theta))
        w = mu / (1 + mu / theta); info = t(X) %*% (X * w)
    } else if (family_name == "Weibull (AFT)") {
        idx = which(data$dead == 1)
        if (length(idx) > p) {
            X_f = X[idx, ]; y_f_log = log(y[idx] + 0.1)
            b_ols = c(fast_ols_cpp(X_f, y_f_log)$coefficients, log(sd(y_f_log) * sqrt(6)/pi))
        } else {
            b_ols = c(fast_ols_cpp(X, log(y + 0.1))$coefficients, 0)
        }
        info = diag(length(b_ols)) 
    }
    list(b = b_ols, w = w, info = info)
}

all_results = list()

run_comprehensive_warm_benchmark = function(family_name, fit_fun, data, alg = NULL, extra_args = list()) {
    cat(sprintf("\n--- Benchmarking Family: %s ---\n", family_name))
    args_base = list(X = data$X, y = data$y)
    if ("dead" %in% names(formals(fit_fun))) args_base$dead = data$dead
    if ("group_id" %in% names(formals(fit_fun))) args_base$group_id = data$group_id
    if ("j_T" %in% names(formals(fit_fun))) args_base$j_T = 2L
    if (!is.null(alg)) args_base$optimization_alg = alg
    args_base = c(args_base, extra_args)
    
    fit_warm = function() {
        starts = get_reasonable_starts(data, family_name)
        a = args_base
        if (family_name == "Logistic (IRLS)") {
            a$warm_start_weights = starts$w
        } else {
            if ("start_params" %in% names(formals(fit_fun))) a$start_params = starts$b
            else if ("warm_start_beta" %in% names(formals(fit_fun))) a$warm_start_beta = starts$b
            else if ("start_par" %in% names(formals(fit_fun))) a$start_par = starts$b
            else if ("start" %in% names(formals(fit_fun))) a$start = starts$b
            if (!is.null(starts$w) && "warm_start_weights" %in% names(formals(fit_fun))) a$warm_start_weights = starts$w
            if (!is.null(starts$info) && "warm_start_fisher_info" %in% names(formals(fit_fun))) a$warm_start_fisher_info = starts$info
        }
        do.call(fit_fun, a)
    }
    
    bm_res = microbenchmark(
        Cold   = do.call(fit_fun, args_base),
        Warm_Flow = fit_warm(),
        times = 3
    )
    print(bm_res)
    summ = as.data.table(bm_res)
    summ = summ[, .(time_ms = mean(time) / 1e6), by = expr]
    summ[, family := family_name]
    all_results[[family_name]] <<- summ

# Run across families
families = list(
    list(name = "Logistic (IRLS)", fun = fast_logistic_regression_with_var_cpp, family = "logistic", alg = "irls"),
    list(name = "Poisson (IRLS)", fun = fast_poisson_regression_with_var_cpp, family = "poisson", alg = "irls"),
    list(name = "NegBin", fun = fast_neg_bin_with_var_cpp, family = "negbin"),
    list(name = "Weibull (AFT)", fun = fast_weibull_regression_cpp, family = "weibull")
)

for (f in families) {
    cat(sprintf("\n--- Starting %s ---\n", f$name))
    data = generate_data(family = f$family)
    run_comprehensive_warm_benchmark(f$name, f$fun, data, alg = f$alg)
}


final_tab = rbindlist(all_results)
cat("\n\n### END-TO-END WARM START BENEFIT SUMMARY (Heuristic: fast_ols_cpp) ###\n")
pivoted = dcast(final_tab, family ~ expr, value.var = "time_ms")
pivoted[, speedup := (Cold - Warm_Flow) / Cold * 100]
print(pivoted[order(-speedup)])
