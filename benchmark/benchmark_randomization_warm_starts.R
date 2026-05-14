
if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The 'pkgload' package is required to benchmark the current EDI source tree.")
}
pkgload::load_all("EDI", quiet = TRUE)
library(EDI)
library(data.table)

set.seed(42)

# Helper to generate data
generate_data = function(n = 150, p = 15, family = "logistic", n_groups = 15) {
    X = matrix(rnorm(n * p), n, p)
    X[, 1] = 1 # Intercept
    X[, 2] = rbinom(n, 1, 0.5)
    beta = rnorm(p) * 0.2 # Small beta for stability in high P
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

final_results = list()

run_randomization_benchmark = function(family_name, fit_fun, data, alg = NULL, extra_args = list(), nsim = 100) {
    cat(sprintf("Benchmarking %s...\n", family_name))
    
    f_names = names(formals(fit_fun))
    args_obs = list(X = data$X, y = data$y)
    if ("dead" %in% f_names && !is.null(data$dead)) args_obs$dead = data$dead
    if ("group_id" %in% f_names && !is.null(data$group_id)) args_obs$group_id = data$group_id
    if ("j_T" %in% f_names) args_obs$j_T = 1L
    if ("Xzi" %in% f_names && !is.null(data$Xzi)) args_obs$Xzi = data$Xzi
    if ("X_hurdle" %in% f_names && !is.null(data$X_hurdle)) args_obs$X_hurdle = data$X_hurdle
    if ("optimization_alg" %in% f_names && !is.null(alg)) args_obs$optimization_alg = alg
    for (n in names(extra_args)) args_obs[[n]] = extra_args[[n]]
    
    obs_fit = tryCatch(do.call(fit_fun, args_obs), error = function(e) NULL)
    if (is.null(obs_fit)) { cat("Observed fit failed.\n"); return(NULL) }
    
    b = obs_fit$params %||% obs_fit$beta %||% obs_fit$b
    if (is.null(b) && !is.null(obs_fit$coefficients)) b = as.numeric(unlist(obs_fit$coefficients))
    
    # Extract warm starts with proper length detection
    if (family_name == "Beta" && length(b) < ncol(data$X) + 1) b = c(b, log(obs_fit$phi))
    if (family_name == "NegBin" && length(b) < ncol(data$X) + 1) b = c(b, log(obs_fit$theta_hat))
    if (family_name %in% c("Weibull (AFT)", "Weibull Frailty") && length(b) < ncol(data$X) + 1) b = c(b, obs_fit$log_sigma)
    
    w = obs_fit$w %||% obs_fit$mu
    info = obs_fit$fisher_information %||% obs_fit$observed_information %||% obs_fit$information %||% obs_fit$XtWX %||% obs_fit$hess_fisher_info_matrix
    
    treatment_idx = 2
    perms = replicate(nsim, sample(data$X[, treatment_idx]))
    
    run_loop = function(ws_type) {
        tm = system.time({
            for (i in 1:nsim) {
                args = args_obs
                args$X[, treatment_idx] = perms[, i]
                if (ws_type == "Beta" || ws_type == "Full") {
                    if ("warm_start_beta" %in% f_names) args$warm_start_beta = if (family_name %in% c("Beta", "NegBin", "Weibull (AFT)")) b[1:ncol(data$X)] else b
                    if ("start_params" %in% f_names) args$start_params = b
                    if ("start_par" %in% f_names) args$start_par = b
                    if ("start" %in% f_names) args$start = b
                }
                if (ws_type == "Weights" || ws_type == "Full") if ("warm_start_weights" %in% f_names) args$warm_start_weights = w
                if (ws_type == "Info" || ws_type == "Full") if ("warm_start_fisher_info" %in% f_names) args$warm_start_fisher_info = info
                tryCatch(do.call(fit_fun, args), error = function(e) NULL)
            }
        })
        as.numeric(tm["elapsed"])
    }
    
    res = list(Cold = run_loop("Cold"))
    test_args = args_obs
    if ("warm_start_beta" %in% f_names) test_args$warm_start_beta = b
    if ("start_params" %in% f_names) test_args$start_params = b
    if ("start_par" %in% f_names) test_args$start_par = b
    if ("start" %in% f_names) test_args$start = b
    if (length(setdiff(names(test_args), names(args_obs))) > 0) res$Beta = run_loop("Beta")
    if ("warm_start_weights" %in% f_names && !is.null(w)) res$Weights = run_loop("Weights")
    if ("warm_start_fisher_info" %in% f_names && !is.null(info)) res$Info = run_loop("Info")
    if (length(res) > 1) res$Full = run_loop("Full")
    
    final_results[[family_name]] <<- res
}

# Hit exactly 17 paths
families = list(
    list(name = "Logistic", fun = fast_logistic_regression_with_var_cpp, family = "logistic"),
    list(name = "Poisson", fun = fast_poisson_regression_with_var_cpp, family = "poisson", alg = "irls"),
    list(name = "Cox", fun = fast_coxph_regression_cpp, family = "cox", alg = "newton_raphson"),
    list(name = "Beta", fun = fast_beta_regression_with_var_cpp, family = "beta", alg = "newton_raphson"),
    list(name = "NegBin", fun = fast_neg_bin_with_var_cpp, family = "negbin", alg = "newton_raphson"),
    list(name = "LogisticGLMM", fun = fast_logistic_glmm_cpp, family = "glmm_logistic"),
    list(name = "PoissonGLMM", fun = fast_poisson_glmm_cpp, family = "glmm_poisson"),
    list(name = "OrdinalGLMM", fun = fast_ordinal_glmm_cpp, family = "ordinal", extra_args = list(K = 3L)),
    list(name = "Ordinal", fun = fast_ordinal_regression_with_var_cpp, family = "ordinal", alg = "newton_raphson"),
    list(name = "ZINB", fun = fast_zinb_cpp, family = "negbin", extra_data = list(Xzi = matrix(rnorm(150*15), 150, 15))),
    list(name = "ZAP", fun = fast_zero_augmented_poisson_cpp, family = "poisson", extra_data = list(Xzi = matrix(rnorm(150*15), 150, 15)), extra_args = list(is_hurdle = TRUE)),
    list(name = "Weibull (AFT)", fun = fast_weibull_regression_cpp, family = "weibull"),
    list(name = "AdjCatLogit", fun = fast_adjacent_category_logit_cpp, family = "ordinal"),
    list(name = "ContRatio", fun = fast_continuation_ratio_regression_cpp, family = "ordinal"),
    list(name = "Stereotype", fun = fast_stereotype_logit_cpp, family = "ordinal"),
    list(name = "GEE (Binomial)", fun = gee_pairs_singletons_cpp, family = "glmm_logistic", 
         extra_data = list(group_id = rep(1:75, each = 2)), extra_args = list(family_str = "binomial")),
    list(name = "LMM (Gaussian)", fun = fast_gaussian_lmm_cpp, family = "lmm")
)

for (f in families) {
    data = generate_data(n = 150, p = 15, family = f$family, n_groups = 15)
    if (!is.null(f$extra_data)) for (n in names(f$extra_data)) data[[n]] = f$extra_data[[n]]
    run_randomization_benchmark(f$name, f$fun, data, alg = f$alg, extra_args = f$extra_args %||% list(), nsim = 100)
}

# Final table assembly
cat("\n\n### FULL RANDOMIZATION BENCHMARK RESULTS (N=150, P=15, 100 permutations)\n\n")
cat("| Path | Cold | Beta-Only | Weights-Only | Info-Only | Full |\n")
cat("| :--- | :---: | :---: | :---: | :---: | :---: |\n")

for (name in names(final_results)) {
    r = final_results[[name]]
    cat(sprintf("| %s | %.3f | %s | %s | %s | %s |\n",
        name,
        r$Cold,
        if (!is.null(r$Beta)) sprintf("%.3f", r$Beta) else "-",
        if (!is.null(r$Weights)) sprintf("%.3f", r$Weights) else "-",
        if (!is.null(r$Info)) sprintf("%.3f", r$Info) else "-",
        if (!is.null(r$Full)) sprintf("%.3f", r$Full) else "-"
    ))
}
