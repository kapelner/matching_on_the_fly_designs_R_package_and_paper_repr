
library(EDI)
library(data.table)
library(parallel)

N_VAL = 300
P_VAL = 5
B_VAL = 10 
R_VAL = 10 
J_VAL = 5  
PB_VAL = 2

RESULTS_CSV = "warm_starts_final_results.csv"

all_inf = grep("^Inference[A-Z]", ls("package:EDI"), value = TRUE)
inf_names = all_inf[!grepl("Abstract|Suite|Mixin|Compound|Likelihood$|Asymp$|Jackknife$|Bootstrap$|MLEorKMSummaryTable|KKRankRegrIVWC|OLS|SimpleMeanDiff|BaiAdjustedT", all_inf)]

# Worker script content
worker_script = '
library(EDI)
library(data.table)
args = commandArgs(trailingOnly = TRUE)
cls_name = args[1]
N_VAL = as.numeric(args[2])
P_VAL = as.numeric(args[3])
B_VAL = as.numeric(args[4])
R_VAL = as.numeric(args[5])
J_VAL = as.numeric(args[6])
PB_VAL = as.numeric(args[7])

measure_avg_time = function(expr_fn, nrep = 50L) {
    times = numeric(nrep)
    for (i in 1:nrep) {
        gc(verbose = FALSE)
        t = tryCatch(system.time(expr_fn())["elapsed"], error = function(e) NA)
        if (is.na(t)) return(NA)
        times[i] = t
    }
    mean(times)
}

get_rt = function(cn) {
    if (grepl("Ordinal|Ridit|Jonck", cn)) return("ordinal")
    if (grepl("Incid", cn)) return("incidence")
    if (grepl("Count", cn)) return("count")
    if (grepl("Prop|Beta", cn)) return("proportion")
    if (grepl("Survival|KM|Cox|Weibull|LogRank", cn)) return("survival")
    return("continuous")
}

get_design_object = function(cn, rt, n) {
    if (grepl("KK", cn) || grepl("Exact|PairedSignTest", cn)) {
        return(DesignFixedBinaryMatch$new(response_type = rt, n = n))
    }
    if (grepl("ExtendedRobins", cn)) {
        return(DesignFixedBlocking$new(response_type = rt, n = n, m = rep(1:(n/2), each=2), equal_block_sizes = FALSE))
    }
    return(DesignFixedBernoulli$new(response_type = rt, n = n))
}

calc_s = function(tc, tw) {
    if (is.na(tc) || is.na(tw) || !is.finite(tc) || !is.finite(tw)) return("N/S")
    if (tc < 0.005) return("< 2%") 
    s = (tc - tw) / tc * 100
    if (s < 2 && s > -2) return("< 2%")
    if (s < 0) return(sprintf("%.1f%%", s))
    sprintf("+%.1f%%", s)
}

rt = get_rt(cls_name)

# Quick trial for dynamic N/P scaling if baseline cold start is too fast (sub-millisecond)
d_temp = get_design_object(cls_name, rt, N_VAL)
if (rt == "ordinal") d_temp$.__enclos_env__$private$ordinal_levels = as.character(1:5)
set.seed(42)
X_temp = as.data.frame(matrix(rnorm(N_VAL * P_VAL), N_VAL, P_VAL))
colnames(X_temp) = paste0("V", 2:(P_VAL+1))
d_temp$add_all_subjects_to_experiment(X_temp)
d_temp$assign_w_to_all_subjects()
y_temp = if(rt == "incidence") { rbinom(N_VAL, 1, 0.5) } else if(rt == "count") { rpois(N_VAL, 1) } else if(rt == "proportion") { runif(N_VAL) } else if(rt == "survival") { rexp(N_VAL) } else { sample(1:5, N_VAL, replace=TRUE) }
d_temp$add_all_subject_responses(y_temp)
inf_temp = get(cls_name)$new(d_temp)
inf_temp$.__enclos_env__$private$fit_warm_start_enabled = FALSE
t_single = tryCatch(system.time(inf_temp$compute_estimate())["elapsed"], error = function(e) NA)

if (!is.na(t_single)) {
    if (t_single < 0.00005) { # extremely fast closed-form models (Wilcoxon, mean difference, etc.)
        N_VAL = 20000
        P_VAL = 5
        B_VAL = 100
        R_VAL = 100
        J_VAL = 30
        PB_VAL = 10
    } else if (t_single < 0.0003) { # very fast models (OLS, simple GLM)
        N_VAL = 8000
        P_VAL = 35
        B_VAL = 50
        R_VAL = 50
        J_VAL = 15
        PB_VAL = 6
    } else if (t_single < 0.001) {
        N_VAL = 3500
        P_VAL = 20
        B_VAL = 25
        R_VAL = 25
        J_VAL = 8
        PB_VAL = 4
    } else if (t_single < 0.003) {
        N_VAL = 1500
        P_VAL = 12
        B_VAL = 15
        R_VAL = 15
        J_VAL = 5
        PB_VAL = 2
    }
}

if (grepl("KK", cls_name) || grepl("Exact|PairedSignTest", cls_name)) {
    N_VAL = min(N_VAL, 300)
}

if (grepl("Wilcox|Jonckheere", cls_name)) {
    N_VAL = min(N_VAL, 1000)
}

d = get_design_object(cls_name, rt, N_VAL)
if (rt == "ordinal") d$.__enclos_env__$private$ordinal_levels = as.character(1:5)

set.seed(42)
X = as.data.frame(matrix(rnorm(N_VAL * P_VAL), N_VAL, P_VAL))
colnames(X) = paste0("V", 2:(P_VAL+1))
d$add_all_subjects_to_experiment(X)
d$assign_w_to_all_subjects()

y = if (rt == "incidence") {
    rbinom(N_VAL, 1, 0.5)
} else if (rt == "count") {
    if (grepl("NegBin", cls_name)) rnbinom(N_VAL, size = 2, mu = 3) else rpois(N_VAL, 1)
} else if (rt == "proportion") {
    if (grepl("Beta", cls_name)) pmin(pmax(rbeta(N_VAL, 3, 7), 1e-4), 1 - 1e-4) else runif(N_VAL, 0.05, 0.95)
} else if (rt == "survival") {
    rexp(N_VAL)
} else {
    sample(1:5, N_VAL, replace = TRUE)
}
d$add_all_subject_responses(y)

# Warm Object Setup
inf_w = get(cls_name)$new(d)
inf_w$.__enclos_env__$private$fit_warm_start_enabled = TRUE
mle = tryCatch(inf_w$compute_estimate(), error = function(e) NULL)
if (!is.null(mle)) {
    pw = inf_w$.__enclos_env__$private
    if (is.list(mle)) {
        pw$set_fit_warm_start(
            mle$params %||% mle$b %||% mle$coefficients, 
            type = if (!is.null(mle$params)) "params" else "beta",
            fisher = mle$fisher_information %||% mle$observed_information %||% mle$XtWX,
            weights = mle$w %||% mle$mu,
            force_pd = TRUE
        )
    } else {
        pw$set_fit_warm_start(as.numeric(mle), "beta")
    }
}

# Cold Object Setup
inf_c = get(cls_name)$new(d)
inf_c$.__enclos_env__$private$fit_warm_start_enabled = FALSE

# Mock context for Bootstrap/JK
mock_ctx = list(n_units = N_VAL, row_to_unit = 1:N_VAL)
inf_c$.__enclos_env__$private$current_bayesian_bootstrap_context = mock_ctx
inf_w$.__enclos_env__$private$current_bayesian_bootstrap_context = mock_ctx

clear_timing_caches = function(priv) {
    old_cache = priv$cached_values
    priv$cached_values = list()
    priv$cached_values$m_cache = old_cache$m_cache
    priv$cached_values$t0s_rand = old_cache$t0s_rand
    priv$cached_values$likelihood_test_eval_cache = list()
    invisible(NULL)
}

# Warm-start policy helper
edi_warm_start_dispatch_policy = getFromNamespace("edi_warm_start_dispatch_policy", "EDI")

# RAND
res_r = "N/S"
if (!(rt == "incidence" && is.null(inf_w$.__enclos_env__$private$custom_randomization_statistic_function))) {
    if (!edi_warm_start_dispatch_policy(cls_name, "rand")) {
        res_r = "Disabled"
    } else {
        tryCatch(inf_c$compute_rand_two_sided_pval(r=1L, show_progress=FALSE), error=function(e)NULL)
        tryCatch(inf_w$compute_rand_two_sided_pval(r=1L, show_progress=FALSE), error=function(e)NULL)
        t_rc = measure_avg_time(function() inf_c$compute_rand_two_sided_pval(r=R_VAL, show_progress=FALSE), nrep = 50L)
        t_rw = measure_avg_time(function() inf_w$compute_rand_two_sided_pval(r=R_VAL, show_progress=FALSE), nrep = 50L)
        res_r = calc_s(t_rc, t_rw)
    }
}

# BOOT
res_b = "Disabled"
if (edi_warm_start_dispatch_policy(cls_name, "non_param_boot")) {
    set.seed(42)
    boot_weight_draws = replicate(
        B_VAL,
        tabulate(sample.int(N_VAL, N_VAL, replace = TRUE), nbins = N_VAL),
        simplify = FALSE
    )
    inf_c$.__enclos_env__$private$active_resampling_operation = "non_param_boot"
    inf_w$.__enclos_env__$private$active_resampling_operation = "non_param_boot"
    tryCatch(inf_c$compute_estimate_with_bootstrap_weights(boot_weight_draws[[1L]], TRUE), error=function(e)NULL)
    tryCatch(inf_w$compute_estimate_with_bootstrap_weights(boot_weight_draws[[1L]], TRUE), error=function(e)NULL)
    t_bc = measure_avg_time(function() {
        clear_timing_caches(inf_c$.__enclos_env__$private)
        for (w in boot_weight_draws) {
            inf_c$compute_estimate_with_bootstrap_weights(w, TRUE)
        }
    }, nrep = 50L)
    t_bw = measure_avg_time(function() {
        clear_timing_caches(inf_w$.__enclos_env__$private)
        for (w in boot_weight_draws) {
            inf_w$compute_estimate_with_bootstrap_weights(w, TRUE)
        }
    }, nrep = 50L)
    inf_c$.__enclos_env__$private$active_resampling_operation = NULL
    inf_w$.__enclos_env__$private$active_resampling_operation = NULL
    res_b = calc_s(t_bc, t_bw)
}

# JK
res_j = "Disabled"
if (edi_warm_start_dispatch_policy(cls_name, "jackknife")) {
    inf_c$.__enclos_env__$private$active_resampling_operation = "jackknife"
    inf_w$.__enclos_env__$private$active_resampling_operation = "jackknife"
    tryCatch({ww=rep(1,N_VAL);ww[1]=0;inf_c$compute_estimate_with_bootstrap_weights(ww,TRUE)}, error=function(e)NULL)
    tryCatch({ww=rep(1,N_VAL);ww[1]=0;inf_w$compute_estimate_with_bootstrap_weights(ww,TRUE)}, error=function(e)NULL)
    t_jc = measure_avg_time(function() {
        for(k in 1:J_VAL){w=rep(1,N_VAL);w[k]=0;inf_c$compute_estimate_with_bootstrap_weights(w,TRUE)}
    }, nrep = 50L)
    t_jw = measure_avg_time(function() {
        for(k in 1:J_VAL){w=rep(1,N_VAL);w[k]=0;inf_w$compute_estimate_with_bootstrap_weights(w,TRUE)}
    }, nrep = 50L)
    res_j = calc_s(t_jc, t_jw)
    inf_c$.__enclos_env__$private$active_resampling_operation = NULL
    inf_w$.__enclos_env__$private$active_resampling_operation = NULL
}

# PB
is_pb_sup = tryCatch(inf_w$.__enclos_env__$private$supports_lik_ratio_param_bootstrap(), error = function(e) FALSE)
res_p = "N/S"
if (isTRUE(is_pb_sup)) {
    if (!edi_warm_start_dispatch_policy(cls_name, "param_boot")) {
        res_p = "Disabled"
    } else {
        tryCatch(inf_c$compute_lik_ratio_bootstrap_two_sided_pval(B=1L, show_progress=FALSE), error=function(e)NULL)
        tryCatch(inf_w$compute_lik_ratio_bootstrap_two_sided_pval(B=1L, show_progress=FALSE), error=function(e)NULL)
        t_pc = measure_avg_time(function() inf_c$compute_lik_ratio_bootstrap_two_sided_pval(B=PB_VAL, show_progress=FALSE), nrep = 50L)
        t_pw = measure_avg_time(function() inf_w$compute_lik_ratio_bootstrap_two_sided_pval(B=PB_VAL, show_progress=FALSE), nrep = 50L)
        res_p = calc_s(t_pc, t_pw)
    }
}

if (cls_name == "InferenceOrdinalPairedSignTest") {
    res_b = "N/S"
    res_j = "N/S"
}

cat(sprintf("%s,%s,%s,%s,%s\n", cls_name, res_r, res_b, res_j, res_p))
'

writeLines(worker_script, "benchmark/worker_v2.R")

cat("Starting definitive robust benchmark...\n")
if (file.exists(RESULTS_CSV)) {
    file.copy(RESULTS_CSV, "warm_starts_final_results_backup.csv", overwrite = TRUE)
}
fwrite(data.table(Path=character(), Rand=character(), Boot=character(), JK=character(), PB=character()), RESULTS_CSV)

num_cores = 6
cat(sprintf("Running with %d parallel workers...\n", num_cores))

mclapply(inf_names, function(cn) {
    cat(sprintf("Processing %s...\n", cn))
    res = tryCatch({
        system2("Rscript", c("benchmark/worker_v2.R", cn, N_VAL, P_VAL, B_VAL, R_VAL, J_VAL, PB_VAL), stdout = TRUE, stderr = TRUE)
    }, error = function(e) {
        return(character(0))
    })
    
    row_str = res[length(res)]
    if (length(row_str) > 0 && grepl(",", row_str)) {
        write(row_str, RESULTS_CSV, append = TRUE)
        cat(sprintf("Done %s: %s\n", cn, row_str))
    } else {
        cat(sprintf("Failed or Crashed: %s\n", cn))
        if (length(res) > 0) {
            cat(paste(res, collapse="\n"), file=stderr())
        }
        write(sprintf("%s,N/S,N/S,N/S,N/S", cn), RESULTS_CSV, append = TRUE)
    }
    NULL
}, mc.cores = num_cores)

# Write/append the disabled paths to CSV if they are not already there
disabled_rows = data.table(
    Path = c(
        "InferenceAllSimpleMeanDiff",
        "InferenceAllSimpleMeanDiffPooledVar",
        "InferenceBaiAdjustedTKK14",
        "InferenceBaiAdjustedTKK21",
        "InferenceContinOLS",
        "InferenceContinKKOLSIVWC",
        "InferenceContinKKOLSOneLik"
    ),
    Rand = "Disabled",
    Boot = "Disabled",
    JK = "Disabled",
    PB = "Disabled"
)
results_dt = fread(RESULTS_CSV)
for (i in 1:nrow(disabled_rows)) {
    row = disabled_rows[i, ]
    if (!row$Path %in% results_dt$Path) {
        fwrite(row, RESULTS_CSV, append = TRUE)
    }
}

cat("Done.\n")
unlink("benchmark/worker_v2.R")
