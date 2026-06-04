library(EDI)
library(data.table)
library(parallel)

set.seed(42)

RESULTS_CSV = Sys.getenv("WARM_START_BENCH_RESULTS", unset = "warm_starts_bayes_jackknife_results.csv")
N_VAL = as.integer(Sys.getenv("WARM_START_BENCH_N", unset = "300"))
P_VAL = as.integer(Sys.getenv("WARM_START_BENCH_P", unset = "5"))
BB_VAL = as.integer(Sys.getenv("WARM_START_BENCH_BB", unset = "10"))
J_VAL = as.integer(Sys.getenv("WARM_START_BENCH_J", unset = "5"))
NREP = as.integer(Sys.getenv("WARM_START_BENCH_REPS", unset = "12"))

md = readLines("package_metadata/warm_starts.md", warn = FALSE)
path_lines = grep("<td .*<b>Inference", md, value = TRUE)
inf_names = unique(sub(".*<b>(Inference[^<]+)</b>.*", "\\1", path_lines))
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

`%||%` = function(a, b) if (!is.null(a)) a else b

measure_avg_time = function(expr_fn, nrep = NREP) {
    times = numeric(nrep)
    for (i in seq_len(nrep)) {
        gc(verbose = FALSE)
        t = tryCatch(system.time(expr_fn())["elapsed"], error = function(e) NA_real_)
        if (is.na(t) || !is.finite(t)) return(NA_real_)
        times[i] = as.numeric(t)
    }
    stats::median(times)
}

calc_s = function(tc, tw) {
    if (is.na(tc) || is.na(tw) || !is.finite(tc) || !is.finite(tw)) return("N/S")
    if (tc < 0.005) return("< 2%")
    s = (tc - tw) / tc * 100
    if (s < 2 && s > -2) return("< 2%")
    if (s < 0) return(sprintf("%.1f%%", s))
    sprintf("+%.1f%%", s)
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

edi_warm_start_dispatch_policy = getFromNamespace("edi_warm_start_dispatch_policy", "EDI")

rt = get_rt(cls_name)
if (grepl("KK", cls_name) || grepl("Exact|PairedSignTest", cls_name)) {
    N_VAL = min(N_VAL, 300L)
}
if (grepl("Wilcox|Jonckheere", cls_name)) {
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
    } else {
        priv_w$set_fit_warm_start(as.numeric(mle), "beta")
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

if (cls_name %in% c("InferenceAllSimpleMeanDiff", "InferenceAllSimpleMeanDiffPooledVar",
                    "InferenceContinOLS", "InferenceContinKKOLSOneLik")) {
    res_bb = "Disabled"
    res_jk = "Disabled"
} else if (grepl("BaiAdjusted|IVWC", cls_name) || cls_name == "InferenceOrdinalPairedSignTest") {
    res_bb = "N/S"
} else if (!edi_warm_start_dispatch_policy(cls_name, "bayesian_boot")) {
    res_bb = "Disabled"
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
    t_bc = measure_avg_time(function() {
        clear_timing_caches(priv_c)
        for (draw in bb_draws) {
            priv_c$current_bayesian_bootstrap_context = draw$context
            inf_c$compute_estimate_with_bootstrap_weights(draw$subject_or_block_weights, TRUE)
        }
    })
    t_bw = measure_avg_time(function() {
        clear_timing_caches(priv_w)
        for (draw in bb_draws) {
            priv_w$current_bayesian_bootstrap_context = draw$context
            inf_w$compute_estimate_with_bootstrap_weights(draw$subject_or_block_weights, TRUE)
        }
    })
    priv_c$active_resampling_operation = NULL
    priv_w$active_resampling_operation = NULL
    res_bb = calc_s(t_bc, t_bw)
}

if (cls_name %in% c("InferenceAllSimpleMeanDiff", "InferenceAllSimpleMeanDiffPooledVar",
                    "InferenceContinOLS", "InferenceContinKKOLSOneLik")) {
    res_jk = "Disabled"
} else if (cls_name == "InferenceOrdinalPairedSignTest") {
    res_jk = "N/S"
} else if (!edi_warm_start_dispatch_policy(cls_name, "jackknife")) {
    res_jk = "Disabled"
} else {
    priv_c$active_resampling_operation = "jackknife"
    priv_w$active_resampling_operation = "jackknife"
    tryCatch({ ww = rep(1, N_VAL); ww[1] = 0; inf_c$compute_estimate_with_bootstrap_weights(ww, TRUE) }, error = function(e) NULL)
    tryCatch({ ww = rep(1, N_VAL); ww[1] = 0; inf_w$compute_estimate_with_bootstrap_weights(ww, TRUE) }, error = function(e) NULL)
    t_jc = measure_avg_time(function() {
        for (k in seq_len(J_VAL)) {
            w = rep(1, N_VAL)
            w[k] = 0
            inf_c$compute_estimate_with_bootstrap_weights(w, TRUE)
        }
    })
    t_jw = measure_avg_time(function() {
        for (k in seq_len(J_VAL)) {
            w = rep(1, N_VAL)
            w[k] = 0
            inf_w$compute_estimate_with_bootstrap_weights(w, TRUE)
        }
    })
    res_jk = calc_s(t_jc, t_jw)
}

cat(sprintf("%s,%s,%s\\n", cls_name, res_bb, res_jk))
'

writeLines(worker_script, "benchmark/worker_bayes_jackknife.R")
fwrite(data.table(Path = character(), Bayesian = character(), JK = character()), RESULTS_CSV)

num_cores = min(6L, max(1L, detectCores() - 1L))
cat(sprintf("Running Bayesian-bootstrap and jackknife timings for %d paths with %d workers...\\n", length(inf_names), num_cores))

mclapply(inf_names, function(cn) {
    cat(sprintf("Processing %s...\\n", cn))
    res = tryCatch({
        system2("Rscript", c("benchmark/worker_bayes_jackknife.R", cn, N_VAL, P_VAL, BB_VAL, J_VAL, NREP), stdout = TRUE, stderr = TRUE)
    }, error = function(e) character(0))
    row_str = res[length(res)]
    if (length(row_str) > 0 && grepl(",", row_str)) {
        write(row_str, RESULTS_CSV, append = TRUE)
        cat(sprintf("Done %s: %s\\n", cn, row_str))
    } else {
        write(sprintf("%s,N/S,N/S", cn), RESULTS_CSV, append = TRUE)
        cat(sprintf("Failed %s\\n", cn))
        if (length(res) > 0) cat(paste(res, collapse = "\\n"), file = stderr())
    }
    NULL
}, mc.cores = num_cores)

unlink("benchmark/worker_bayes_jackknife.R")
cat("Done. Results: ", RESULTS_CSV, "\\n", sep = "")
