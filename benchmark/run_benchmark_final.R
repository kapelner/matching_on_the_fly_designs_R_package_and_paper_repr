
library(EDI)
library(data.table)
library(parallel)

N_VAL = as.integer(Sys.getenv("WARM_START_BENCH_N", unset = "300"))
P_VAL = as.integer(Sys.getenv("WARM_START_BENCH_P", unset = "5"))
B_VAL = as.integer(Sys.getenv("WARM_START_BENCH_B", unset = "10"))
R_VAL = as.integer(Sys.getenv("WARM_START_BENCH_R", unset = "10"))
J_VAL = as.integer(Sys.getenv("WARM_START_BENCH_J", unset = "5"))
PB_VAL = as.integer(Sys.getenv("WARM_START_BENCH_PB", unset = "2"))
NREP = as.integer(Sys.getenv("WARM_START_BENCH_REPS", unset = "100"))
FIXED_N = as.integer(Sys.getenv("WARM_START_BENCH_FIXED_N", unset = "1"))

RESULTS_CSV = Sys.getenv("WARM_START_BENCH_RESULTS", unset = "warm_starts_final_results.csv")

md = readLines("package_metadata/warm_starts.md", warn = FALSE)
path_lines = grep("<td[^>]*>Inference", md, value = TRUE)
inf_names = unique(sub(".*<td[^>]*>(Inference[^<]+)</td>.*", "\\1", path_lines))
excluded_paths = c(
    "InferenceAllSimpleMeanDiff",
    "InferenceAllSimpleMeanDiffPooledVar"
)
inf_names = setdiff(inf_names, excluded_paths)
path_filter = Sys.getenv("WARM_START_BENCH_PATHS", unset = "")
if (nzchar(path_filter)) {
    requested_paths = trimws(strsplit(path_filter, ",", fixed = TRUE)[[1L]])
    inf_names = intersect(inf_names, requested_paths)
}

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
PB_VAL = as.numeric(args[7])\nNREP = as.integer(args[8])\nFIXED_N = as.logical(as.integer(args[9]))

`%||%` = function(a, b) if (!is.null(a)) a else b

measure_times = function(expr_fn, nrep = NREP) {
    times = numeric(nrep)
    for (i in 1:nrep) {
        gc(verbose = FALSE)
        t = tryCatch(system.time(expr_fn())["elapsed"], error = function(e) NA)
        if (is.na(t)) return(NA)
        times[i] = t
    }
    times
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

calc_sig_s = function(cold_times, warm_times, alpha = 0.01) {
    if (length(cold_times) == 1L && is.na(cold_times)) return("N/S")
    if (length(warm_times) == 1L && is.na(warm_times)) return("N/S")
    ok = is.finite(cold_times) & is.finite(warm_times)
    cold_times = cold_times[ok]
    warm_times = warm_times[ok]
    if (length(cold_times) < 2L) return("N/S")
    warm_wins = sum(warm_times < cold_times)
    cold_wins = sum(cold_times < warm_times)
    non_ties = warm_wins + cold_wins
    if (non_ties < 2L) return("same")
    pval = tryCatch(
        stats::prop.test(c(warm_wins, cold_wins), c(non_ties, non_ties), correct = FALSE)$p.value,
        error = function(e) NA_real_
    )
    if (is.na(pval) || !is.finite(pval) || pval >= alpha) return("same")
    calc_s(mean(cold_times), mean(warm_times))
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

if (!isTRUE(FIXED_N) && !is.na(t_single)) {
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

if (!isTRUE(FIXED_N) && (grepl("KK", cls_name) || grepl("Exact|PairedSignTest", cls_name))) {
    N_VAL = min(N_VAL, 300)
}

if (!isTRUE(FIXED_N) && grepl("Wilcox|Jonckheere", cls_name)) {
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

set_treatment_for_timing = function(priv, w_new) {
    priv$w = as.numeric(w_new)
    priv$cached_design_matrix = NULL
    priv$cached_w_for_design_matrix = NULL
    priv$cached_reduced_X = NULL
    priv$cached_j_treat_for_reduced = NULL
    if (!is.null(priv$des_obj_priv_int)) {
        priv$des_obj_priv_int$w = priv$w
        priv$des_obj_priv_int$all_subject_data_cache = list()
    }
    invisible(NULL)
}

restore_treatment_for_timing = function(priv, w_orig) {
    set_treatment_for_timing(priv, w_orig)
    invisible(NULL)
}

compute_rand_estimate_only = function(inf_obj) {
    priv = inf_obj$.__enclos_env__$private
    res = tryCatch(
        priv$compute_treatment_estimate_during_randomization_inference(estimate_only = TRUE),
        error = function(e) tryCatch(
            priv$compute_treatment_estimate_during_randomization_inference(),
            error = function(e2) tryCatch(inf_obj$compute_estimate(estimate_only = TRUE), error = function(e3) NA_real_)
        )
    )
    if (is.null(res)) res = priv$cached_values$beta_hat_T
    as.numeric(res)[1L]
}

load_param_boot_response_for_timing = function(priv, draw) {
    worker_data = draw$worker_data %||% draw
    y_new = worker_data$y
    if (is.null(y_new)) return(FALSE)
    priv$y = as.numeric(y_new)
    priv$y_temp = priv$y
    if (!is.null(worker_data$dead)) {
        priv$dead = as.numeric(worker_data$dead)
        priv$any_censoring = any(priv$dead == 0)
    }
    extra_names = setdiff(names(worker_data), c("y", "dead"))
    for (nm in extra_names) {
        priv[[nm]] = worker_data[[nm]]
    }
    if (!is.null(priv$des_obj_priv_int)) {
        priv$des_obj_priv_int$y = priv$y
        if (!is.null(worker_data$dead)) priv$des_obj_priv_int$dead = priv$dead
        priv$des_obj_priv_int$all_subject_data_cache = list()
    }
    invisible(TRUE)
}

restore_param_boot_response_for_timing = function(priv, y_orig, dead_orig) {
    priv$y = y_orig
    priv$y_temp = y_orig
    priv$dead = dead_orig
    priv$any_censoring = !is.null(dead_orig) && any(dead_orig == 0)
    if (!is.null(priv$des_obj_priv_int)) {
        priv$des_obj_priv_int$y = y_orig
        priv$des_obj_priv_int$dead = dead_orig
        priv$des_obj_priv_int$all_subject_data_cache = list()
    }
    invisible(NULL)
}

simulate_param_boot_draws_for_timing = function(priv, B, delta = 0) {
    if (!isTRUE(tryCatch(priv$supports_lik_ratio_param_bootstrap(), error = function(e) FALSE))) return(NULL)
    spec = tryCatch(priv$get_likelihood_test_spec(), error = function(e) NULL)
    if (is.null(spec)) return(NULL)
    eval_obs = tryCatch(
        priv$get_memoized_likelihood_test_eval(
            delta = delta,
            testing_type = "lik_ratio",
            spec = spec,
            include_full_negloglik = TRUE,
            include_null_negloglik = TRUE
        ),
        error = function(e) NULL
    )
    if (is.null(eval_obs) || isTRUE(eval_obs$invalid) || is.null(eval_obs$null_fit)) return(NULL)
    draws = vector("list", as.integer(B))
    for (b in seq_len(as.integer(B))) {
        draws[[b]] = tryCatch(priv$simulate_under_lik_null(spec, delta, eval_obs$null_fit), error = function(e) NULL)
        if (is.null(draws[[b]]) || is.null(draws[[b]]$worker_data$y)) return(NULL)
    }
    draws
}

# RAND
res_r = "N/S"
set.seed(42)
rand_w_draws = replicate(R_VAL, sample(inf_w$.__enclos_env__$private$w), simplify = FALSE)
orig_w_c = inf_c$.__enclos_env__$private$w
orig_w_w = inf_w$.__enclos_env__$private$w
inf_c$.__enclos_env__$private$active_resampling_operation = "rand"
inf_w$.__enclos_env__$private$active_resampling_operation = "rand"
tryCatch({ set_treatment_for_timing(inf_c$.__enclos_env__$private, rand_w_draws[[1L]]); clear_timing_caches(inf_c$.__enclos_env__$private); compute_rand_estimate_only(inf_c) }, error=function(e)NULL)
tryCatch({ set_treatment_for_timing(inf_w$.__enclos_env__$private, rand_w_draws[[1L]]); clear_timing_caches(inf_w$.__enclos_env__$private); compute_rand_estimate_only(inf_w) }, error=function(e)NULL)
t_rc = measure_times(function() {
    for (w_perm in rand_w_draws) {
        set_treatment_for_timing(inf_c$.__enclos_env__$private, w_perm)
        clear_timing_caches(inf_c$.__enclos_env__$private)
        compute_rand_estimate_only(inf_c)
    }
}, nrep = NREP)
t_rw = measure_times(function() {
    for (w_perm in rand_w_draws) {
        set_treatment_for_timing(inf_w$.__enclos_env__$private, w_perm)
        clear_timing_caches(inf_w$.__enclos_env__$private)
        compute_rand_estimate_only(inf_w)
    }
}, nrep = NREP)
restore_treatment_for_timing(inf_c$.__enclos_env__$private, orig_w_c)
restore_treatment_for_timing(inf_w$.__enclos_env__$private, orig_w_w)
inf_c$.__enclos_env__$private$active_resampling_operation = NULL
inf_w$.__enclos_env__$private$active_resampling_operation = NULL
res_r = calc_sig_s(t_rc, t_rw)

# BOOT
res_b = "N/S"
if (TRUE) {
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
    t_bc = measure_times(function() {
        for (w in boot_weight_draws) {
            clear_timing_caches(inf_c$.__enclos_env__$private)
            inf_c$compute_estimate_with_bootstrap_weights(w, TRUE)
        }
    }, nrep = NREP)
    t_bw = measure_times(function() {
        for (w in boot_weight_draws) {
            clear_timing_caches(inf_w$.__enclos_env__$private)
            inf_w$compute_estimate_with_bootstrap_weights(w, TRUE)
        }
    }, nrep = NREP)
    inf_c$.__enclos_env__$private$active_resampling_operation = NULL
    inf_w$.__enclos_env__$private$active_resampling_operation = NULL
    res_b = calc_sig_s(t_bc, t_bw)
}

# JK
res_j = "N/S"
if (TRUE) {
    inf_c$.__enclos_env__$private$active_resampling_operation = "jackknife"
    inf_w$.__enclos_env__$private$active_resampling_operation = "jackknife"
    tryCatch({ww=rep(1,N_VAL);ww[1]=0;inf_c$compute_estimate_with_bootstrap_weights(ww,TRUE)}, error=function(e)NULL)
    tryCatch({ww=rep(1,N_VAL);ww[1]=0;inf_w$compute_estimate_with_bootstrap_weights(ww,TRUE)}, error=function(e)NULL)
    t_jc = measure_times(function() {
        for(k in 1:J_VAL){clear_timing_caches(inf_c$.__enclos_env__$private);w=rep(1,N_VAL);w[k]=0;inf_c$compute_estimate_with_bootstrap_weights(w,TRUE)}
    }, nrep = NREP)
    t_jw = measure_times(function() {
        for(k in 1:J_VAL){clear_timing_caches(inf_w$.__enclos_env__$private);w=rep(1,N_VAL);w[k]=0;inf_w$compute_estimate_with_bootstrap_weights(w,TRUE)}
    }, nrep = NREP)
    res_j = calc_sig_s(t_jc, t_jw)
    inf_c$.__enclos_env__$private$active_resampling_operation = NULL
    inf_w$.__enclos_env__$private$active_resampling_operation = NULL
}

# PB
is_pb_sup = tryCatch(inf_w$.__enclos_env__$private$supports_lik_ratio_param_bootstrap(), error = function(e) FALSE)
res_p = "N/S"
if (isTRUE(is_pb_sup)) {
    set.seed(42)
    pb_draws = simulate_param_boot_draws_for_timing(inf_w$.__enclos_env__$private, PB_VAL, delta = 0)
    if (!is.null(pb_draws)) {
        orig_y_c = inf_c$.__enclos_env__$private$y
        orig_y_w = inf_w$.__enclos_env__$private$y
        orig_dead_c = inf_c$.__enclos_env__$private$dead
        orig_dead_w = inf_w$.__enclos_env__$private$dead
        inf_c$.__enclos_env__$private$active_resampling_operation = "param_boot"
        inf_w$.__enclos_env__$private$active_resampling_operation = "param_boot"
        tryCatch({ load_param_boot_response_for_timing(inf_c$.__enclos_env__$private, pb_draws[[1L]]); clear_timing_caches(inf_c$.__enclos_env__$private); inf_c$compute_estimate(estimate_only = TRUE) }, error=function(e)NULL)
        tryCatch({ load_param_boot_response_for_timing(inf_w$.__enclos_env__$private, pb_draws[[1L]]); clear_timing_caches(inf_w$.__enclos_env__$private); inf_w$compute_estimate(estimate_only = TRUE) }, error=function(e)NULL)
        t_pc = measure_times(function() {
            for (draw in pb_draws) {
                load_param_boot_response_for_timing(inf_c$.__enclos_env__$private, draw)
                clear_timing_caches(inf_c$.__enclos_env__$private)
                inf_c$compute_estimate(estimate_only = TRUE)
            }
        }, nrep = NREP)
        t_pw = measure_times(function() {
            for (draw in pb_draws) {
                load_param_boot_response_for_timing(inf_w$.__enclos_env__$private, draw)
                clear_timing_caches(inf_w$.__enclos_env__$private)
                inf_w$compute_estimate(estimate_only = TRUE)
            }
        }, nrep = NREP)
        restore_param_boot_response_for_timing(inf_c$.__enclos_env__$private, orig_y_c, orig_dead_c)
        restore_param_boot_response_for_timing(inf_w$.__enclos_env__$private, orig_y_w, orig_dead_w)
        inf_c$.__enclos_env__$private$active_resampling_operation = NULL
        inf_w$.__enclos_env__$private$active_resampling_operation = NULL
        res_p = calc_sig_s(t_pc, t_pw)
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

num_cores = 6
cat(sprintf("Running with %d parallel workers...\n", num_cores))

rows = mclapply(inf_names, function(cn) {
    cat(sprintf("Processing %s...\n", cn))
    res = tryCatch({
        system2("Rscript", c("benchmark/worker_v2.R", cn, N_VAL, P_VAL, B_VAL, R_VAL, J_VAL, PB_VAL, NREP, FIXED_N), stdout = TRUE, stderr = TRUE)
    }, error = function(e) {
        return(character(0))
    })
    
    row_str = res[length(res)]
    if (length(row_str) > 0 && grepl(",", row_str)) {
        cat(sprintf("Done %s: %s\n", cn, row_str))
        return(row_str)
    } else {
        cat(sprintf("Failed or Crashed: %s\n", cn))
        if (length(res) > 0) {
            cat(paste(res, collapse="\n"), file=stderr())
        }
        return(sprintf("%s,N/S,N/S,N/S,N/S", cn))
    }
}, mc.cores = num_cores)

rows = unlist(rows, use.names = FALSE)
parts = strsplit(rows, ",", fixed = TRUE)
dt = rbindlist(lapply(parts, function(x) {
    if (length(x) != 5L) x = c(x[1L], rep("N/S", 4L))
    data.table(Path = x[1L], Rand = x[2L], Boot = x[3L], JK = x[4L], PB = x[5L])
}))
fwrite(dt, RESULTS_CSV)

cat("Done.\n")
unlink("benchmark/worker_v2.R")
