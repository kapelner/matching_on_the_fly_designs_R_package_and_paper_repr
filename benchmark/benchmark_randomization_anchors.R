
if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The 'pkgload' package is required.")
}
pkgload::load_all("EDI", quiet = TRUE)
library(EDI)
library(data.table)

set.seed(42)

# Tiny for verification
N = 10
P = 7
R = 5

# Helper to generate data with "Strong Signal"
generate_strong_signal_data = function(n = 100, p = 7, family = "logistic") {
    X = matrix(rnorm(n * p), n, p)
    X[, 1] = 1 
    beta = rnorm(p) * 2.5
    eta = X %*% beta
    if (family == "logistic") {
        y = rbinom(n, 1, 1 / (1 + exp(-eta)))
    } else if (family == "poisson") {
        y = rpois(n, exp(pmin(pmax(eta, -4), 4)))
    } else if (family == "negbin") {
        y = rnbinom(n, size = 2, mu = exp(pmin(pmax(eta, -4), 4)))
    }
    list(X = X, y = y)
}

all_benchmarks = list()

benchmark_sampling_v2 = function(family_name, fit_fun, data, extra_args = list()) {
    cat(sprintf("\n--- Benchmarking %s ---\n", family_name))
    args_obs = list(X = data$X, y = data$y)
    if ("j_T" %in% names(formals(fit_fun))) args_obs$j_T = 1L
    args_obs = c(args_obs, extra_args)
    obs_fit = do.call(fit_fun, args_obs)
    b_mle = obs_fit$params %||% obs_fit$b
    w_mle = obs_fit$w %||% obs_fit$mu
    i_mle = obs_fit$fisher_information %||% obs_fit$observed_information %||% obs_fit$XtWX
    perms = lapply(1:R, function(i) sample(data$X[, P]))
    
    cat("  Cold... "); flush.console()
    t_rand_cold = system.time({
        for (p_val in perms) {
            a = args_obs; a$X[, P] = p_val; a$smart_cold_start = FALSE
            do.call(fit_fun, a)
        }
    })["elapsed"]
    
    cat("Warm-MLE... "); flush.console()
    t_rand_warm_mle = system.time({
        for (p_val in perms) {
            a = args_obs; a$X[, P] = p_val
            a$start_beta = b_mle; a$warm_start_weights = w_mle; a$warm_start_fisher_info = i_mle
            do.call(fit_fun, a)
        }
    })["elapsed"]
    
    cat("Warm-Null... "); flush.console()
    t_rand_warm_null = system.time({
        a1 = args_obs; a1$X[, P] = perms[[1]]; a1$smart_cold_start = TRUE
        fit1 = do.call(fit_fun, a1)
        if (!is.null(fit1)) {
            b_null = fit1$params %||% fit1$b
            w_null = fit1$w %||% fit1$mu
            i_null = fit1$fisher_information %||% fit1$observed_information %||% fit1$XtWX
            for (i in 2:R) {
                a = args_obs; a$X[, P] = perms[[i]]
                a$start_beta = b_null; a$warm_start_weights = w_null; a$warm_start_fisher_info = i_null
                do.call(fit_fun, a)
            }
        }
    })["elapsed"]
    cat("Done.\n")
    all_benchmarks[[family_name]] <<- data.table(family = family_name, Cold = t_rand_cold, Warm_MLE = t_rand_warm_mle, Warm_Null = t_rand_warm_null)
}

# Just Logistic for speed test
f = list(name = "Logistic", fun = fast_logistic_regression_with_var_cpp, family = "logistic")
data = generate_strong_signal_data(n = N, p = P, family = f$family)
benchmark_sampling_v2(f$name, f$fun, data)

print(all_benchmarks[[1]])
