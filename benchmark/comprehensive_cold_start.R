
if (!requireNamespace("pkgload", quietly = TRUE)) stop("The 'pkgload' package is required.")
pkgload::load_all("EDI", quiet = TRUE)
library(EDI)
library(microbenchmark)
library(data.table)

set.seed(42)

# Clinical trial scale data generation
generate_data = function(n = 100, p = 7, family = "logistic", n_groups = 10) {
    X = matrix(rnorm(n * p), n, p)
    X[, 1] = 1 
    beta = rnorm(p) * 2.5 # Strong signal
    eta = X %*% beta
    group_id = rep(1:n_groups, length.out = n)
    strata = rep(1:5, length.out = n)
    
    res = list(X = X, group_id = group_id, strata = strata)
    
    if (family %in% c("logistic", "log-binomial", "identity-binomial", "glmm_logistic", "robust")) {
        mu = 1 / (1 + exp(-eta))
        res$y = rbinom(n, 1, mu)
    } else if (family %in% c("poisson", "glmm_poisson", "zap", "quasipoisson", "hurdle_poisson_glmm")) {
        mu = exp(pmin(eta, 4))
        res$y = rpois(n, mu)
    } else if (family %in% c("negbin", "zinb", "truncated_negbin", "hurdle_negbin")) {
        mu = exp(pmin(eta, 4))
        res$y = rnbinom(n, size = 2, mu = mu)
        if (family == "zinb") {
            z = rbinom(n, 1, 0.2)
            res$y[z == 1] = 0
        }
        if (family == "truncated_negbin") {
            res$y[res$y == 0] = 1
        }
    } else if (family %in% c("beta", "zoib")) {
        mu = 1 / (1 + exp(-eta))
        phi = 10
        res$y = rbeta(n, mu * phi, (1 - mu) * phi)
        res$y = pmax(pmin(res$y, 1 - 1e-6), 1e-6)
    } else if (family %in% c("weibull", "cox", "strat_cox", "dep_cens")) {
        mu = exp(pmin(eta, 4))
        res$y = rexp(n, 1/mu) 
        res$dead = rbinom(n, 1, 0.8)
    } else if (family %in% c("ordinal", "adj_cat", "continuation", "stereotype", "glmm_ordinal")) {
        # 4 levels
        p1 = 1 / (1 + exp(-(-3 - eta)))
        p_le_2 = 1 / (1 + exp(-(0 - eta)))
        p_le_3 = 1 / (1 + exp(-(3 - eta)))
        p2 = pmax(0, p_le_2 - p1)
        p3 = pmax(0, p_le_3 - p_le_2)
        p4 = pmax(0, 1 - p_le_3)
        probs = cbind(p1, p2, p3, p4)
        probs = probs / rowSums(probs)
        res$y = apply(probs, 1, function(p) sample(1:4, 1, prob = p))
    }
    res
}

all_results = list()

run_bm = function(name, fit_fun, data, extra_args = list()) {
    cat(sprintf("Benchmarking %s... ", name))
    flush.console()
    
    args_base = list(X = data$X, y = data$y)
    f_names = names(formals(fit_fun))
    
    if ("dead" %in% f_names) args_base$dead = data$dead
    if ("group_id" %in% f_names) args_base$group_id = data$group_id
    if ("K" %in% f_names) args_base$K = length(unique(data$y))
    if ("j_T" %in% f_names) args_base$j_T = 2L
    if ("strata" %in% f_names) args_base$strata = data$strata
    if ("X_hurdle" %in% f_names) args_base$X_hurdle = data$X
    if ("Xzi" %in% f_names) args_base$Xzi = data$X
    if ("is_hurdle" %in% f_names) args_base$is_hurdle = TRUE
    if ("X_zero_one" %in% f_names) args_base$X_zero_one = data$X
    
    args_base = c(args_base, extra_args)
    
    # Verify both can run
    s_fit = tryCatch(do.call(fit_fun, c(args_base, list(smart_cold_start = TRUE))), error = function(e) { cat(sprintf("Smart failed: %s ", e$message)); NULL })
    d_fit = tryCatch(do.call(fit_fun, c(args_base, list(smart_cold_start = FALSE))), error = function(e) { cat(sprintf("Dumb failed: %s ", e$message)); NULL })
    
    if (is.null(s_fit) || is.null(d_fit)) {
        cat("Model failed. Skipping.\n")
        return(NULL)
    }

    s_iter = if (!is.null(s_fit$iterations)) s_fit$iterations else if (!is.null(s_fit$niter)) s_fit$niter else NA
    d_iter = if (!is.null(d_fit$iterations)) d_fit$iterations else if (!is.null(d_fit$niter)) d_fit$niter else NA
    
    bm = microbenchmark(
        Smart = do.call(fit_fun, c(args_base, list(smart_cold_start = TRUE))),
        Dumb  = do.call(fit_fun, c(args_base, list(smart_cold_start = FALSE))),
        times = 10
    )
    cat(sprintf("Iters S:%s D:%s. Done.\n", as.character(s_iter), as.character(d_iter)))
    summ = as.data.table(bm)
    summ = summ[, .(time_ms = mean(time) / 1e6), by = expr]
    summ[, family := name]
    summ[, iterations := ifelse(expr == "Smart", s_iter, d_iter)]
    all_results[[name]] <<- summ
}

# Expanded list of all paths with smart_cold_start
models = list(
    list(name = "Logistic", fun = fast_logistic_regression_with_var_cpp, family = "logistic", args = list(optimization_alg = "irls")),
    list(name = "Poisson", fun = fast_poisson_regression_with_var_cpp, family = "poisson", args = list(optimization_alg = "irls")),
    list(name = "NegBin", fun = fast_neg_bin_with_var_cpp, family = "negbin"),
    list(name = "Beta", fun = fast_beta_regression_with_var_cpp, family = "beta"),
    list(name = "Weibull", fun = fast_weibull_regression_cpp, family = "weibull"),
    list(name = "CoxPH", fun = fast_coxph_regression_cpp, family = "cox"),
    list(name = "Stratified Cox", fun = fast_stratified_coxph_regression_cpp, family = "strat_cox"),
    list(name = "Ordinal (Logit)", fun = fast_ordinal_regression_with_var_cpp, family = "ordinal"),
    list(name = "Ordinal (Cauchit)", fun = fast_ordinal_cauchit_regression_with_var_cpp, family = "ordinal"),
    list(name = "Ordinal (Probit)", fun = fast_ordinal_probit_regression_with_var_cpp, family = "ordinal"),
    list(name = "Ordinal (Cloglog)", fun = fast_ordinal_cloglog_regression_with_var_cpp, family = "ordinal"),
    list(name = "Adj. Cat. Logit", fun = fast_adjacent_category_logit_with_var_cpp, family = "adj_cat"),
    list(name = "Continuation Ratio", fun = fast_continuation_ratio_regression_with_var_cpp, family = "continuation"),
    list(name = "Stereotype Logit", fun = fast_stereotype_logit_with_var_cpp, family = "stereotype"),
    list(name = "Log-Binomial", fun = fast_log_binomial_regression_with_var_cpp, family = "log-binomial"),
    list(name = "Identity Binomial", fun = fast_identity_binomial_regression_with_var_cpp, family = "identity-binomial"),
    list(name = "ZINB", fun = fast_zinb_cpp, family = "zinb"),
    list(name = "ZAP", fun = fast_zero_augmented_poisson_cpp, family = "zap"),
    list(name = "ZOIB", fun = fast_zero_one_inflated_beta_cpp, family = "zoib"),
    list(name = "Logistic GLMM", fun = fast_logistic_glmm_cpp, family = "glmm_logistic"),
    list(name = "Poisson GLMM", fun = fast_poisson_glmm_cpp, family = "glmm_poisson"),
    list(name = "Ordinal GLMM", fun = fast_ordinal_glmm_cpp, family = "glmm_ordinal"),
    list(name = "Hurdle Poisson GLMM", fun = fast_hurdle_poisson_glmm_cpp, family = "hurdle_poisson_glmm"),
    list(name = "Hurdle NegBin", fun = fast_hurdle_negbin_with_var_cpp, family = "hurdle_negbin"),
    list(name = "Truncated NegBin", fun = fast_truncated_negbin_count_cpp, family = "truncated_negbin"),
    list(name = "Robust Regression", fun = fast_robust_regression_cpp, family = "robust"),
    list(name = "Dep Censoring Trans", fun = fast_dep_cens_transform_optim_cpp, family = "dep_cens"),
    list(name = "Quasipoisson", fun = fast_quasipoisson_regression_with_var_cpp, family = "quasipoisson")
)

for (m in models) {
    d = generate_data(family = m$family)
    run_bm(m$name, m$fun, d, extra_args = m$args)
}

final = rbindlist(all_results)
pivoted = dcast(final, family ~ expr, value.var = c("time_ms", "iterations"))
pivoted[, speedup := (time_ms_Dumb - time_ms_Smart) / time_ms_Dumb * 100]
print(pivoted[order(-speedup)])
