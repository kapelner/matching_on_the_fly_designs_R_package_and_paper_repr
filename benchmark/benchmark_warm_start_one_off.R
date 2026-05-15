
if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The 'pkgload' package is required.")
}
pkgload::load_all("EDI", quiet = TRUE)
library(EDI)
library(microbenchmark)
library(data.table)

set.seed(42)

# Helper to generate data (Very small for speed)
generate_data = function(n = 200, p = 2, family = "logistic", n_groups = 10) {
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
    } else if (family == "negbin" || family == "zinb") {
        mu = exp(eta)
        y = rnbinom(n, size = 2, mu = mu)
        if (family == "zinb") {
            z = rbinom(n, 1, 0.2)
            y[z == 1] = 0
        }
    } else if (family == "beta") {
        mu = 1 / (1 + exp(-eta))
        phi = 10
        y = rbeta(n, mu * phi, (1 - mu) * phi)
        y = pmax(pmin(y, 1 - 1e-6), 1e-6)
    } else if (family == "weibull" || family == "cox") {
        mu = exp(eta)
        y = rexp(n, 1/mu) 
        dead = rbinom(n, 1, 0.8)
        return(list(X = X, y = y, dead = dead, group_id = group_id))
    } else if (family == "ordinal") {
        p1 = 1 / (1 + exp(-(-1 - eta)))
        p_le_2 = 1 / (1 + exp(-(1 - eta)))
        p2 = pmax(0, p_le_2 - p1)
        p3 = pmax(0, 1 - p_le_2)
        probs = cbind(p1, p2, p3)
        probs = probs / rowSums(probs)
        y = apply(probs, 1, function(p) sample(1:3, 1, prob = p))
    }
    
    list(X = X, y = y, group_id = group_id)
}

# Heuristic starts using fast_ols_cpp
get_reasonable_starts = function(data, family_name) {
    X = data$X; y = data$y; n = nrow(X); p = ncol(X)
    b_ols = NULL; w = NULL; info = NULL
    
    if (family_name == "Logistic (IRLS)" || family_name == "LogisticGLMM") {
        p_avg = mean(y)
        mu = rep(p_avg, n)
        w = mu * (1 - mu)
        mu_init = (y + 0.5) / 2
        eta_init = qlogis(mu_init)
        b_ols = fast_ols_cpp(X, eta_init)$coefficients
        if (family_name == "LogisticGLMM") b_ols = c(b_ols, 0)
        info = t(X) %*% (X * w)
    } else if (family_name == "Poisson (IRLS)" || family_name == "ZAP") {
        eta = log(y + 0.1)
        b_ols = fast_ols_cpp(X, eta)$coefficients
        w = exp(eta)
        info = t(X) %*% (X * w)
        if (family_name == "ZAP") b_ols = c(b_ols, b_ols) 
    } else if (family_name == "NegBin" || family_name == "ZINB") {
        eta = log(y + 0.1)
        b_ols = c(fast_ols_cpp(X, eta)$coefficients, 0) 
        w = exp(eta) / 2 
        info = t(X) %*% (X * w)
        if (family_name == "ZINB") b_ols = c(b_ols, b_ols)
    } else if (family_name == "Weibull (AFT)") {
        idx = which(data$dead == 1)
        if (length(idx) > p) {
            X_f = X[idx, ]; y_f_log = log(y[idx] + 0.1)
            b_ols = c(fast_ols_cpp(X_f, y_f_log)$coefficients, log(sd(y_f_log)*sqrt(6)/pi))
        } else {
            b_ols = c(fast_ols_cpp(X, log(y + 0.1))$coefficients, 0)
        }
        info = diag(length(b_ols))
    } else if (family_name == "Beta") {
        eta = qlogis(pmax(pmin(y, 0.99), 0.01))
        b_ols = c(fast_ols_cpp(X, eta)$coefficients, log(10))
        info = diag(length(b_ols))
    } else if (family_name == "Ordinal") {
        K = length(unique(y))
        b_ols = c(seq(-1, 1, length.out = K-1), rep(0, p))
        info = diag(length(b_ols))
    }
    list(b = b_ols, w = w, info = info)
}

all_results = list()

run_comprehensive_warm_benchmark = function(family_name, fit_fun, data, alg = NULL, extra_args = list()) {
    cat(sprintf("Benchmarking %s... ", family_name))
    args_base = list(X = data$X, y = data$y)
    if ("dead" %in% names(formals(fit_fun))) args_base$dead = data$dead
    if ("group_id" %in% names(formals(fit_fun))) args_base$group_id = data$group_id
    if ("j_T" %in% names(formals(fit_fun))) args_base$j_T = 2L
    if (!is.null(alg)) args_base$optimization_alg = alg
    if (family_name %in% c("ZINB", "ZAP")) {
        args_base$Xzi = data$X
        if (family_name == "ZAP") args_base$is_hurdle = TRUE
    }
    args_base = c(args_base, extra_args)
    
    fit_warm = function() {
        starts = get_reasonable_starts(data, family_name)
        a = args_base
        if (family_name == "Logistic (IRLS)") {
            a$warm_start_weights = starts$w
        } else {
            if ("start_params" %in% names(formals(fit_fun))) a$start_params = starts$b
            else if ("start_beta" %in% names(formals(fit_fun))) a$start_beta = starts$b
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
    cat("Done.\n")
    summ = as.data.table(bm_res)
    summ = summ[, .(time_ms = mean(time) / 1e6), by = expr]
    summ[, family := family_name]
    all_results[[family_name]] <<- summ
}

families = list(
    list(name = "Logistic (IRLS)", fun = fast_logistic_regression_with_var_cpp, family = "logistic", alg = "irls"),
    list(name = "Poisson (IRLS)", fun = fast_poisson_regression_with_var_cpp, family = "poisson", alg = "irls"),
    list(name = "NegBin", fun = fast_neg_bin_with_var_cpp, family = "negbin"),
    list(name = "Beta", fun = fast_beta_regression_with_var_cpp, family = "beta"),
    list(name = "Weibull (AFT)", fun = fast_weibull_regression_cpp, family = "weibull"),
    list(name = "Ordinal", fun = fast_ordinal_regression_with_var_cpp, family = "ordinal"),
    list(name = "LogisticGLMM", fun = fast_logistic_glmm_cpp, family = "glmm_logistic"),
    list(name = "ZINB", fun = fast_zinb_cpp, family = "negbin"),
    list(name = "ZAP", fun = fast_zero_augmented_poisson_cpp, family = "poisson")
)

for (f in families) {
    data = generate_data(family = f$family)
    run_comprehensive_warm_benchmark(f$name, f$fun, data, alg = f$alg)
}

final_tab = rbindlist(all_results)
cat("\n\n### END-TO-END WARM START BENEFIT SUMMARY (Heuristic: fast_ols_cpp) ###\n")
pivoted = dcast(final_tab, family ~ expr, value.var = "time_ms")
pivoted[, speedup := (Cold - Warm_Flow) / Cold * 100]
print(pivoted[order(-speedup)])
