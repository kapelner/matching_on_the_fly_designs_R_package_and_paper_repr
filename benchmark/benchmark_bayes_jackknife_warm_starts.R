library(EDI)
library(data.table)
library(parallel)

set.seed(42)

RESULTS_CSV = Sys.getenv("WARM_START_BENCH_RESULTS", unset = "warm_starts_bayes_jackknife_results.csv")
N_VAL = as.integer(Sys.getenv("WARM_START_BENCH_N", unset = "300"))
P_VAL = as.integer(Sys.getenv("WARM_START_BENCH_P", unset = "5"))
BB_VAL = as.integer(Sys.getenv("WARM_START_BENCH_BB", unset = "10"))
J_VAL = as.integer(Sys.getenv("WARM_START_BENCH_J", unset = "5"))
NREP = as.integer(Sys.getenv("WARM_START_BENCH_REPS", unset = "100"))
FIXED_N = as.integer(Sys.getenv("WARM_START_BENCH_FIXED_N", unset = "1"))

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

worker_script = '
library(EDI)

args = commandArgs(trailingOnly = TRUE)
cls_name = args[1]
N_VAL = as.integer(args[2])
P_VAL = as.integer(args[3])
BB_VAL = as.integer(args[4])
J_VAL = as.integer(args[5])
NREP = as.integer(args[6])
FIXED_N = as.logical(as.integer(args[7]))

`%||%` = function(a, b) if (!is.null(a)) a else b

measure_times = function(expr_fn, nrep = NREP) {
    times = numeric(nrep)
    for (i in seq_len(nrep)) {
        gc(verbose = FALSE)
        t = tryCatch(system.time(expr_fn())["elapsed"], error = function(e) NA_real_)
        if (is.na(t) || !is.finite(t)) return(NA_real_)
        times[i] = as.numeric(t)
    }
    times
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

get_rt = function(cn) {
    if (grepl("Ordinal|Ridit|Jonck", cn)) return("ordinal")
    if (grepl("Incid", cn)) return("incidence")
    if (grepl("Count", cn)) return("count")
    if (grepl("Prop|Beta", cn)) return("proportion")
    if (grepl("Survival|KM|Cox|Weibull|LogRank", cn)) return("survival")
    "continuous"
}

get_design_object = function(cn, rt, n) {
    if (grepl("KK", cn) || grepl("Exact|PairedSignTest", cn)) {
        return(DesignFixedBinaryMatch$new(response_type = rt, n = n))
    }
    if (grepl("ExtendedRobins", cn)) {
        return(DesignFixedBlocking$new(response_type = rt, n = n, m = rep(seq_len(n / 2), each = 2), equal_block_sizes = FALSE))
    }
    DesignFixedBernoulli$new(response_type = rt, n = n)
}

rt = get_rt(cls_name)
if (!isTRUE(FIXED_N) && (grepl("KK", cls_name) || grepl("Exact|PairedSignTest", cls_name))) {
    N_VAL = min(N_VAL, 300L)
}
if (!isTRUE(FIXED_N) && grepl("Wilcox|Jonckheere", cls_name)) {
    N_VAL = min(N_VAL, 1000L)
}

d = get_design_object(cls_name, rt, N_VAL)
if (rt == "ordinal") d$.__enclos_env__$private$ordinal_levels = as.character(1:5)

set.seed(42)
X = as.data.frame(matrix(rnorm(N_VAL * P_VAL), N_VAL, P_VAL))
colnames(X) = paste0("V", 2:(P_VAL + 1L))
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
} else if (rt == "continuous") {
    as.numeric(as.matrix(X) %*% rnorm(P_VAL, 0, 0.1) + rnorm(N_VAL))
} else {
    sample(1:5, N_VAL, replace = TRUE)
}
d$add_all_subject_responses(y, dead = rep(1, N_VAL))

res_bb = "N/S"
res_jk = "N/S"

inf_w = get(cls_name)$new(d)
inf_c = get(cls_name)$new(d)
priv_w = inf_w$.__enclos_env__$private
priv_c = inf_c$.__enclos_env__$private
priv_w$fit_warm_start_enabled = TRUE
priv_c$fit_warm_start_enabled = FALSE

mle = tryCatch(inf_w$compute_estimate(), error = function(e) NULL)
if (!is.null(mle)) {
    if (is.list(mle)) {
        priv_w$set_fit_warm_start(
            mle$params %||% mle$b %||% mle$coefficients,
            type = if (!is.null(mle$params)) "params" else "beta",
            fisher = mle$fisher_information %||% mle$observed_information %||% mle$XtWX,
            weights = mle$w %||% mle$mu,
            force_pd = TRUE
        )
    }
}

mock_ctx = list(n_units = N_VAL, row_to_unit = seq_len(N_VAL))
priv_c$current_bayesian_bootstrap_context = mock_ctx
priv_w$current_bayesian_bootstrap_context = mock_ctx

clear_timing_caches = function(priv) {
    old_cache = priv$cached_values
    priv$cached_values = list()
    priv$cached_values$m_cache = old_cache$m_cache
    priv$cached_values$t0s_rand = old_cache$t0s_rand
    priv$cached_values$likelihood_test_eval_cache = list()
    invisible(NULL)
}

if (grepl("BaiAdjusted|IVWC", cls_name) || cls_name == "InferenceOrdinalPairedSignTest") {
    res_bb = "N/S"
} else {
    set.seed(42)
    bb_draws = replicate(
        BB_VAL,
        priv_w$bayesian_bootstrap_sample_weights(weighting_unit_type = NULL),
        simplify = FALSE
    )
    priv_c$active_resampling_operation = "bayesian_boot"
    priv_w$active_resampling_operation = "bayesian_boot"
    tryCatch({
        draw = bb_draws[[1L]]
        priv_c$current_bayesian_bootstrap_context = draw$context
        inf_c$compute_estimate_with_bootstrap_weights(draw$subject_or_block_weights, TRUE)
    }, error = function(e) NULL)
    tryCatch({
        draw = bb_draws[[1L]]
        priv_w$current_bayesian_bootstrap_context = draw$context
        inf_w$compute_estimate_with_bootstrap_weights(draw$subject_or_block_weights, TRUE)
    }, error = function(e) NULL)
    t_bc = measure_times(function() {
        for (draw in bb_draws) {
            clear_timing_caches(priv_c)
            priv_c$current_bayesian_bootstrap_context = draw$context
            inf_c$compute_estimate_with_bootstrap_weights(draw$subject_or_block_weights, TRUE)
        }
    })
    t_bw = measure_times(function() {
        for (draw in bb_draws) {
            clear_timing_caches(priv_w)
            priv_w$current_bayesian_bootstrap_context = draw$context
            inf_w$compute_estimate_with_bootstrap_weights(draw$subject_or_block_weights, TRUE)
        }
    })
    priv_c$active_resampling_operation = NULL
    priv_w$active_resampling_operation = NULL
    res_bb = calc_sig_s(t_bc, t_bw)
}

if (cls_name == "InferenceOrdinalPairedSignTest") {
    res_jk = "N/S"
} else {
    priv_c$active_resampling_operation = "jackknife"
    priv_w$active_resampling_operation = "jackknife"
    tryCatch({ ww = rep(1, N_VAL); ww[1] = 0; inf_c$compute_estimate_with_bootstrap_weights(ww, TRUE) }, error = function(e) NULL)
    tryCatch({ ww = rep(1, N_VAL); ww[1] = 0; inf_w$compute_estimate_with_bootstrap_weights(ww, TRUE) }, error = function(e) NULL)
    t_jc = measure_times(function() {
        for (k in seq_len(J_VAL)) {
            clear_timing_caches(priv_c)
            w = rep(1, N_VAL)
            w[k] = 0
            inf_c$compute_estimate_with_bootstrap_weights(w, TRUE)
        }
    })
    t_jw = measure_times(function() {
        for (k in seq_len(J_VAL)) {
            clear_timing_caches(priv_w)
            w = rep(1, N_VAL)
            w[k] = 0
            inf_w$compute_estimate_with_bootstrap_weights(w, TRUE)
        }
    })
    res_jk = calc_sig_s(t_jc, t_jw)
}

cat(sprintf("%s,%s,%s\\n", cls_name, res_bb, res_jk))
'

writeLines(worker_script, "benchmark/worker_bayes_jackknife.R")

num_cores = min(6L, max(1L, detectCores() - 1L))
cat(sprintf("Running Bayesian-bootstrap and jackknife timings for %d paths with %d workers...\\n", length(inf_names), num_cores))

rows = mclapply(inf_names, function(cn) {
    cat(sprintf("Processing %s...\\n", cn))
    res = tryCatch({
        system2("Rscript", c("benchmark/worker_bayes_jackknife.R", cn, N_VAL, P_VAL, BB_VAL, J_VAL, NREP, FIXED_N), stdout = TRUE, stderr = TRUE)
    }, error = function(e) character(0))
    row_str = res[length(res)]
    if (length(row_str) > 0 && grepl(",", row_str)) {
        cat(sprintf("Done %s: %s\\n", cn, row_str))
        return(row_str)
    } else {
        cat(sprintf("Failed %s\\n", cn))
        if (length(res) > 0) cat(paste(res, collapse = "\\n"), file = stderr())
        return(sprintf("%s,N/S,N/S", cn))
    }
}, mc.cores = num_cores)

rows = unlist(rows, use.names = FALSE)
parts = strsplit(rows, ",", fixed = TRUE)
dt = rbindlist(lapply(parts, function(x) {
    if (length(x) != 3L) x = c(x[1L], rep("N/S", 2L))
    data.table(Path = x[1L], Bayesian = x[2L], JK = x[3L])
}))
fwrite(dt, RESULTS_CSV)

unlink("benchmark/worker_bayes_jackknife.R")
cat("Done. Results: ", RESULTS_CSV, "\\n", sep = "")
