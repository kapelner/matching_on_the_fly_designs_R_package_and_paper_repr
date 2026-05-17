
if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The 'pkgload' package is required to benchmark the current EDI source tree.")
}
pkgload::load_all("EDI", quiet = TRUE)
library(EDI)

set.seed(42)

# Helper to generate a completed design
setup_benchmark_des = function(n = 150, p = 5, family = "logistic") {
    X = matrix(rnorm(n * p), n, p)
    X[, 1] = 1 # Intercept
    beta = rnorm(p) * 0.2
    
    if (family == "logistic") {
        y = rbinom(n, 1, plogis(X %*% beta))
        des = DesignFixedBernoulli$new(n = n, response_type = "incidence")
        des$add_all_subjects_to_experiment(as.data.frame(X[, -1]))
        des$assign_w_to_all_subjects()
        des$add_all_subject_responses(y)
    } else if (family == "negbin") {
        y = rnbinom(n, size = 2, mu = exp(X %*% beta))
        des = DesignFixedBernoulli$new(n = n, response_type = "count")
        des$add_all_subjects_to_experiment(as.data.frame(X[, -1]))
        des$assign_w_to_all_subjects()
        des$add_all_subject_responses(y)
    }
    des
}

run_full_resampling_benchmark = function(name, des, inf_class) {
    cat(sprintf("\n--- %s Resampling Benchmarks ---\n", name))
    
    # 1. Randomization (100 perms)
    inf_cold = inf_class$new(des, smart_cold_start_default = FALSE)
    tm_rand_cold = system.time(inf_cold$approximate_randomization_distribution_beta_hat_T(r = 100, show_progress = FALSE))["elapsed"]
    
    inf_warm = inf_class$new(des, smart_cold_start_default = TRUE)
    tm_rand_warm = system.time(inf_warm$approximate_randomization_distribution_beta_hat_T(r = 100, show_progress = FALSE))["elapsed"]
    
    # 2. Jackknife
    tm_jk_cold = system.time(inf_cold$compute_jackknife_standard_error())["elapsed"]
    tm_jk_warm = system.time(inf_warm$compute_jackknife_standard_error())["elapsed"]
    
    # 3. Non-parametric Bootstrap (100 draws)
    tm_boot_cold = system.time(inf_cold$approximate_bootstrap_distribution_beta_hat_T(B = 100, show_progress = FALSE))["elapsed"]
    tm_boot_warm = system.time(inf_warm$approximate_bootstrap_distribution_beta_hat_T(B = 100, show_progress = FALSE))["elapsed"]
    
    # 4. CI Inversion (Likelihood Ratio)
    inf_cold$set_testing_type("lik_ratio")
    tm_ci_cold = system.time(inf_cold$compute_asymp_confidence_interval())["elapsed"]
    
    inf_warm$set_testing_type("lik_ratio")
    tm_ci_warm = system.time(inf_warm$compute_asymp_confidence_interval())["elapsed"]
    
    # 5. Parametric Bootstrap (LRT calibration, 20 draws)
    # We use a smaller B because parametric boot is expensive
    tm_pboot_cold = system.time(inf_cold$compute_bootstrapped_lik_ratio_test_pval(B = 20, show_progress = FALSE))["elapsed"]
    tm_pboot_warm = system.time(inf_warm$compute_bootstrapped_lik_ratio_test_pval(B = 20, show_progress = FALSE))["elapsed"]
    
    res = data.frame(
        Procedure = c("Randomization (100)", "Jackknife (N)", "NP-Bootstrap (100)", "CI Inversion (LR)", "Param-Bootstrap (20)"),
        Cold_sec = as.numeric(c(tm_rand_cold, tm_jk_cold, tm_boot_cold, tm_ci_cold, tm_pboot_cold)),
        Warm_sec = as.numeric(c(tm_rand_warm, tm_jk_warm, tm_boot_warm, tm_ci_warm, tm_pboot_warm))
    )
    res$Speedup = sprintf("%.1f%%", (1 - res$Warm_sec / res$Cold_sec) * 100)
    print(res)
    res
}

# Logistic (Medium Tier)
des_log = setup_benchmark_des(family = "logistic")
res_log = run_full_resampling_benchmark("Logistic (Medium)", des_log, InferenceIncidLogRegr)

# NegBin (Heavy Tier)
des_nb = setup_benchmark_des(family = "negbin")
res_nb = run_full_resampling_benchmark("NegBin (Heavy)", des_nb, InferenceCountNegBin)

# Final formatted output for the .md file
cat("\n\n### Formatted for warm_starts.md\n\n")

print_table = function(name, df) {
    cat(sprintf("#### %s Results\n\n", name))
    cat("| Procedure | Cold Start (s) | Warm Start (s) | Speedup |\n")
    cat("| :--- | :---: | :---: | :---: |\n")
    for(i in 1:nrow(df)) {
        cat(sprintf("| %s | %.3f | %.3f | **%s** |\n", df$Procedure[i], df$Cold_sec[i], df$Warm_sec[i], df$Speedup[i]))
    }
    cat("\n")
}

print_table("Logistic (Medium Tier)", res_log)
print_table("Negative Binomial (Heavy Tier)", res_nb)
