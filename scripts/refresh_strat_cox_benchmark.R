.libPaths(c(file.path(Sys.getenv("HOME"), "R", paste0(R.version$platform, "-library"), paste(R.version$major, sub("\\..*$", "", R.version$minor), sep = ".")), .libPaths()))
library(EDI)
library(bench)
library(data.table)
library(survival)

TARGET = "InferenceSurvivalStratCoxPHRegr"

# ---- Shared timing infrastructure ----
B_TIME = 10
TARGET_BATCH_MS = 200
MIN_RESOLVED_BATCH_MS = 10
MAX_INNER_REPS = 100000L
FAST_PATH_MICROBENCH_REPS = 2000L
FAST_PATH_THRESHOLD_MS = 0.01

collect_timing_ms = function(expr, setup = NULL, times = B_TIME, env = parent.frame(),
    target_batch_ms = TARGET_BATCH_MS, max_inner_reps = MAX_INNER_REPS,
    fast_path_microbenchmark_reps = FAST_PATH_MICROBENCH_REPS) {
  gctorture(FALSE)
  gc(verbose = FALSE)
  if (is.null(setup)) setup = quote({})
  
  micro_time_ms = function() {
    gctorture(FALSE)
    gc(verbose = FALSE)
    bm_res = tryCatch({
      bench::mark(
        total = {
          eval(setup, envir = env)
          eval(expr, envir = env)
        },
        setup_only = {
          eval(setup, envir = env)
        },
        iterations = fast_path_microbenchmark_reps,
        check = FALSE,
        filter_gc = TRUE,
        memory = FALSE
      )
    }, error = function(e) NULL)
    if (is.null(bm_res)) {
      return(list(median_ms = NA_real_, samples_ms = numeric(0), method = "bench::mark"))
    }
    t_total = bm_res$time[[which(bm_res$expression == "total")]] * 1000
    t_setup = bm_res$time[[which(bm_res$expression == "setup_only")]] * 1000
    
    median_total = as.numeric(bm_res$median[bm_res$expression == "total"]) * 1000
    median_setup = as.numeric(bm_res$median[bm_res$expression == "setup_only"]) * 1000
    median_diff = max(0, median_total - median_setup)
    
    samples_ms = t_total - median_setup
    samples_ms = samples_ms[is.finite(samples_ms) & samples_ms > 0]
    list(
      median_ms = median_diff,
      samples_ms = samples_ms,
      method = "bench::mark"
    )
  }
  inner_reps = 1L; batch_ms = 0
  while (inner_reps < max_inner_reps) {
    batch_ms = system.time(for (j in seq_len(inner_reps)) {
      eval(setup, envir = env)
      eval(expr, envir = env)
    })[["elapsed"]] * 1000
    if (is.finite(batch_ms) && batch_ms >= min(target_batch_ms, MIN_RESOLVED_BATCH_MS)) break
    inner_reps = min(max_inner_reps, inner_reps * 2L)
  }
  if (!is.finite(batch_ms) || batch_ms < MIN_RESOLVED_BATCH_MS) return(micro_time_ms())
  if (is.finite(batch_ms) && batch_ms > 0)
    inner_reps = min(max_inner_reps, max(1L, ceiling(inner_reps * target_batch_ms / batch_ms)))
  vals_ms = numeric(times)
  for (i in seq_len(times)) {
    gctorture(FALSE)
    gc(verbose = FALSE)
    t_total = system.time(for (j in seq_len(inner_reps)) {
      eval(setup, envir = env)
      eval(expr, envir = env)
    })[["elapsed"]] * 1000
    
    gctorture(FALSE)
    gc(verbose = FALSE)
    t_setup = system.time(for (j in seq_len(inner_reps)) {
      eval(setup, envir = env)
    })[["elapsed"]] * 1000
    
    vals_ms[i] = max(0, t_total - t_setup) / inner_reps
  }
  vals_ms = vals_ms[is.finite(vals_ms) & vals_ms > 0]
  median_ms = if (length(vals_ms) == 0L) NA_real_ else median(vals_ms)
  if (!is.finite(median_ms) || median_ms < FAST_PATH_THRESHOLD_MS) return(micro_time_ms())
  list(median_ms = median_ms, samples_ms = vals_ms, method = "system.time")
}

timing_ttest_pval = function(x, y) {
  x = x[is.finite(x)]; y = y[is.finite(y)]
  if (length(x) < 2L || length(y) < 2L) return(NA_real_)
  tryCatch(stats::t.test(x, y, var.equal = FALSE)$p.value, error = function(e) NA_real_)
}

format_ms = function(x) {
  ifelse(is.na(x), "NA",
    ifelse(x < 0.01, format(x, scientific = TRUE, digits = 3, trim = TRUE),
           format(round(x, 2), nsmall = 2, trim = TRUE)))
}
format_pval = function(x) {
  ifelse(is.na(x), "NA", vapply(x, function(v) sprintf("%.3g", v), character(1)))
}
format_pval_stars = function(x) {
  ifelse(is.na(x), "",
    ifelse(x < 0.001, "***", ifelse(x < 0.01, "**", ifelse(x < 0.05, "*", ""))))
}
row_bg_color = function(speedup, pval) {
  if (!is.finite(speedup) || is.na(pval)) return("#eceff1")
  if (pval < 0.05 && speedup > 1) return("#d9fdd3")
  if (pval < 0.05 && speedup < 1) return("#ffd9d9")
  if (pval >= 0.05) return("#fff4bf")
  ""
}

make_tr = function(cls, resp, edi_ms, pkg, func, can_ms, speedup_num, pval_num) {
  speedup_str = if (!is.na(speedup_num)) paste0(round(speedup_num, 2), "x") else "NA"
  bg = row_bg_color(speedup_num, pval_num)
  style = if (nzchar(bg)) paste0(" style=\"background-color: ", bg, ";\"") else ""
  paste0(
    "    <tr", style, "><td>", cls, "</td><td>", resp, "</td><td>", format_ms(edi_ms),
    "</td><td>", pkg, "</td><td>", func, "</td><td>", format_ms(can_ms),
    "</td><td>", speedup_str, "</td><td>", format_pval(pval_num),
    "</td><td>", format_pval_stars(pval_num), "</td></tr>"
  )
}

# ================================================================
# Data generation helpers (from benchmark_model_fits.R)
# ================================================================
N_SURV = 500

generate_data = function(n = 500, p = 5) {
  X = matrix(rnorm(n * p), n, p); X[, 1] = 1
  beta = rnorm(p) * 0.5
  n_treat = floor(n / 2)
  w = sample(c(rep(1, n_treat), rep(0, n - n_treat)))
  eta = X %*% beta + 0.5 * w
  y = rexp(n, exp(eta)); dead = rbinom(n, 1, 0.8)
  list(X = X, w = w, y = y, dead = dead)
}

make_true_stratified_survival_data = function(d) {
  n = nrow(d$X)
  strata_grid = as.matrix(expand.grid(x1 = 0:1, x2 = 0:2))
  idx = sample(rep(seq_len(nrow(strata_grid)), length.out = n))
  d$X[, 2] = strata_grid[idx, 1]; d$X[, 3] = strata_grid[idx, 2]
  beta_cov = c(0.45, -0.35, 0.20, -0.15)
  eta = 0.5 * d$w + drop(d$X[, 2:5, drop = FALSE] %*% beta_cov)
  d$y = rexp(n, exp(eta)); d$dead = rbinom(n, 1, 0.8)
  d
}

build_strat_cox_canonical_inputs = function(d) {
  X_cov = d$X[, -1, drop = FALSE]
  strata_info = EDI:::compute_survival_strata_ids_cpp(as.matrix(X_cov))
  X_linear = matrix(numeric(0), nrow = nrow(X_cov), ncol = 0)
  if (!is.null(X_cov) && ncol(X_cov) > 0) {
    keep_cols = setdiff(seq_len(ncol(X_cov)), as.integer(strata_info$selected_cols))
    if (length(keep_cols) > 0) {
      full_design = cbind(w = d$w, X_cov[, keep_cols, drop = FALSE])
      reduced = EDI:::drop_linearly_dependent_cols(full_design)$M
      if ("w" %in% colnames(reduced))
        X_linear = reduced[, colnames(reduced) != "w", drop = FALSE]
    }
  }
  strata_id = as.integer(strata_info$strata_id)
  informative_rows = integer(0)
  for (s in unique(strata_id)) {
    i_s = which(strata_id == s)
    if (length(i_s) < 2L || length(unique(d$w[i_s])) < 2L || !any(d$dead[i_s] == 1, na.rm = TRUE)) next
    informative_rows = c(informative_rows, i_s)
  }
  informative_rows = sort(unique(informative_rows))
  if (length(informative_rows) >= 4L) {
    x = if (ncol(X_linear) > 0)
      cbind(w = d$w[informative_rows], X_linear[informative_rows, , drop = FALSE])
    else matrix(d$w[informative_rows], ncol = 1L, dimnames = list(NULL, "w"))
    return(list(X = as.matrix(x), y = as.numeric(d$y[informative_rows]),
                dead = as.numeric(d$dead[informative_rows]),
                strata = as.integer(strata_id[informative_rows])))
  }
  x = if (ncol(X_linear) > 0) cbind(w = d$w, X_linear) else matrix(d$w, ncol = 1L, dimnames = list(NULL, "w"))
  list(X = as.matrix(x), y = as.numeric(d$y), dead = as.numeric(d$dead), strata = NULL)
}

# ================================================================
# prepare_solver_only_edi  (from benchmark_wald_tests.R)
# ================================================================
prepare_solver_only_edi = function(inf_obj, cls_name) {
  priv = inf_obj$.__enclos_env__$private
  replace_binding = function(name, value) {
    if (!exists(name, envir = priv, inherits = FALSE)) return(invisible(FALSE))
    was_locked = bindingIsLocked(name, priv)
    if (was_locked) unlockBinding(name, priv)
    assign(name, value, envir = priv)
    if (was_locked) lockBinding(name, priv)
    invisible(TRUE)
  }
  replace_binding("harden", FALSE)
  if (identical(cls_name, "InferenceSurvivalStratCoxPHRegr") &&
      exists("compute_strata_info", envir = priv, inherits = FALSE)) {
    X_full = tryCatch(priv$X, error = function(e) NULL)
    strata_info = tryCatch(priv$compute_strata_info(X_full), error = function(e) NULL)
    if (!is.null(strata_info)) {
      replace_binding("compute_strata_info", function(X_full) strata_info)
      informative_rows = tryCatch(priv$get_informative_rows(strata_info$strata_id), error = function(e) integer(0))
      replace_binding("get_informative_rows", function(strata_id) informative_rows)
      keep_cols = setdiff(seq_len(ncol(X_full)), as.integer(strata_info$selected_cols))
      X_linear = tryCatch(
        if (length(keep_cols) > 0L)
          priv$reduce_covariates_preserving_treatment(X_full[, keep_cols, drop = FALSE])
        else matrix(numeric(0), nrow = nrow(X_full), ncol = 0L),
        error = function(e) NULL)
      if (!is.null(X_linear)) {
        replace_binding("reduce_covariates_preserving_treatment", function(X) X_linear)
        if (length(informative_rows) >= 4L) {
          rcpp_inp = tryCatch(priv$build_rcpp_inputs(informative_rows, X_linear), error = function(e) NULL)
          if (!is.null(rcpp_inp))
            replace_binding("build_rcpp_inputs", function(rows, X_linear) rcpp_inp)
        }
      }
    }
  }
  invisible(inf_obj)
}

safe_instantiate = function(cls_name, des) {
  tryCatch(
    get(cls_name)$new(des),
    error = function(e) get(cls_name)$new(des)
  )
}

# ================================================================
# Run: estimate-only (model-fits benchmark)
# ================================================================
set.seed(42)
d_mf = generate_data(n = N_SURV)
d_mf = make_true_stratified_survival_data(d_mf)

X_cov_mf = d_mf$X[, -1, drop = FALSE]
colnames(X_cov_mf) = paste0("x", seq_len(ncol(X_cov_mf)))
X_ord_mf  = cbind(treatment = d_mf$w, X_cov_mf)
y_mf      = as.numeric(d_mf$y)
dead_mf   = as.numeric(d_mf$dead)

# EDI object (model-fits style)
des_mf = DesignFixediBCRD$new(n = nrow(X_ord_mf), response_type = "survival")
des_mf$add_all_subjects_to_experiment(as.data.frame(X_cov_mf))
des_mf$overwrite_all_subject_assignments(d_mf$w)
des_mf$add_all_subject_responses(y_mf, deads = dead_mf)
e_mf = new.env(parent = globalenv())
e_mf$inf_obj_strat_cox = tryCatch(
  InferenceSurvivalStratCoxPHRegr$new(des_mf, smart_cold_start_default = FALSE),
  error = function(err) InferenceSurvivalStratCoxPHRegr$new(des_mf)
)
expr_edi_mf = quote(inf_obj_strat_cox$compute_estimate(estimate_only = TRUE))
setup_edi_mf = quote({
  priv = inf_obj_strat_cox$.__enclos_env__$private
  priv$cached_mod = NULL
  priv$cached_values = list()
  priv$cox_data_cache = NULL
})
cat("Benchmarking model-fits EDI (estimate_only)...\n")
eval(expr_edi_mf, envir = e_mf)   # validation run
t_edi_mf = collect_timing_ms(expr_edi_mf, setup = setup_edi_mf, env = e_mf)

# Canonical (model-fits style)
strat_inputs_mf = build_strat_cox_canonical_inputs(d_mf)
X_can_strat_mf  = strat_inputs_mf$X
y_can_strat_mf  = strat_inputs_mf$y
dead_can_strat_mf = strat_inputs_mf$dead
strata_can_mf   = strat_inputs_mf$strata
expr_can_mf = quote(
  survival::coxph.fit(
    x = X_can_strat_mf, y = survival::Surv(y_can_strat_mf, dead_can_strat_mf),
    strata = strata_can_mf, offset = NULL, init = NULL,
    control = survival::coxph.control(), weights = NULL, method = "breslow",
    rownames = as.character(seq_along(y_can_strat_mf))
  )
)
cat("Benchmarking model-fits canonical...\n")
t_can_mf = collect_timing_ms(expr_can_mf)

speedup_mf  = t_can_mf$median_ms / t_edi_mf$median_ms
pval_mf     = timing_ttest_pval(t_edi_mf$samples_ms, t_can_mf$samples_ms)
cat(sprintf("  EDI=%.2f ms  Can=%.2f ms  Speedup=%.2fx  pval=%.3g\n",
            t_edi_mf$median_ms, t_can_mf$median_ms, speedup_mf, pval_mf))

# ================================================================
# Run: Wald test benchmark
# ================================================================
set.seed(42)
N_WALD = 200
d_wd = generate_data(n = N_WALD)
d_wd = make_true_stratified_survival_data(d_wd)

X_cov_wd = d_wd$X[, -1, drop = FALSE]
colnames(X_cov_wd) = paste0("x", seq_len(ncol(X_cov_wd)))
y_wd   = as.numeric(d_wd$y); dead_wd = as.numeric(d_wd$dead)

des_wd = DesignFixediBCRD$new(n = N_WALD, response_type = "survival")
des_wd$add_all_subjects_to_experiment(as.data.frame(X_cov_wd))
des_wd$overwrite_all_subject_assignments(d_wd$w)
des_wd$add_all_subject_responses(y_wd, deads = dead_wd)

inf_obj_wd = safe_instantiate(TARGET, des_wd)
prepare_solver_only_edi(inf_obj_wd, TARGET)
priv_env_wd = inf_obj_wd$.__enclos_env__$private

expr_edi_wd = quote({
  inf_obj_wd$compute_estimate(estimate_only = FALSE)
  p = inf_obj_wd$compute_asymp_two_sided_pval(delta = 0)
  if (!is.finite(p)) stop("Non-finite EDI p-value.")
  p
})
setup_edi_wd = quote({
  priv = inf_obj_wd$.__enclos_env__$private
  priv$cached_mod = NULL
  priv$cached_values = list()
  priv$cox_data_cache = NULL
  if (exists("clear_fit_warm_start", envir = priv, inherits = FALSE)) {
    priv$clear_fit_warm_start()
  }
})
cat("Benchmarking Wald EDI...\n")
eval(expr_edi_wd)   # warmup
t_edi_wd = collect_timing_ms(expr_edi_wd, setup = setup_edi_wd)

# Canonical Wald
strat_inputs_wd = build_strat_cox_canonical_inputs(d_wd)
X_can_strat_wd  = strat_inputs_wd$X
y_can_strat_wd  = strat_inputs_wd$y
dead_can_strat_wd = strat_inputs_wd$dead
strata_can_wd   = strat_inputs_wd$strata
expr_can_wd = quote({
  res = survival::coxph.fit(
    x = X_can_strat_wd, y = survival::Surv(y_can_strat_wd, dead_can_strat_wd),
    strata = strata_can_wd, offset = NULL, init = NULL,
    control = survival::coxph.control(), weights = NULL, method = "breslow",
    rownames = as.character(seq_along(y_can_strat_wd))
  )
  v = as.matrix(res$var)
  2 * pnorm(-abs(res$coefficients[1] / sqrt(v[1, 1])))
})
cat("Benchmarking Wald canonical...\n")
t_can_wd = collect_timing_ms(expr_can_wd)

speedup_wd = t_can_wd$median_ms / t_edi_wd$median_ms
pval_wd    = timing_ttest_pval(t_edi_wd$samples_ms, t_can_wd$samples_ms)
cat(sprintf("  EDI=%.2f ms  Can=%.2f ms  Speedup=%.2fx  pval=%.3g\n",
            t_edi_wd$median_ms, t_can_wd$median_ms, speedup_wd, pval_wd))

# ================================================================
# Patch benchmark_model_fits.md
# ================================================================
md_path = "package_metadata/benchmark_model_fits.md"
lines = readLines(md_path)

# Build replacement rows
new_tr_mf = make_tr(TARGET, "survival",
  t_edi_mf$median_ms, "survival", "coxph.fit(strat)",
  t_can_mf$median_ms, speedup_mf, pval_mf)

new_tr_wd = make_tr(TARGET, "survival",
  t_edi_wd$median_ms, "survival", "coxph.fit(strat,breslow)+Wald",
  t_can_wd$median_ms, speedup_wd, pval_wd)

# Replace the two existing StratCox <tr> rows
pat_mf = paste0(".*<td>", TARGET, "</td>.*coxph\\.fit\\(strat\\)<.*")
pat_wd = paste0(".*<td>", TARGET, "</td>.*coxph\\.fit\\(strat,breslow\\)\\+Wald<.*")

idx_mf = grep(pat_mf, lines)
idx_wd = grep(pat_wd, lines)

if (length(idx_mf) == 1L) {
  lines[idx_mf] = new_tr_mf
  cat("Patched model-fits row at line", idx_mf, "\n")
} else {
  cat("WARNING: found", length(idx_mf), "model-fits rows — not patching\n")
}
if (length(idx_wd) == 1L) {
  lines[idx_wd] = new_tr_wd
  cat("Patched Wald row at line", idx_wd, "\n")
} else {
  cat("WARNING: found", length(idx_wd), "Wald rows — not patching\n")
}

writeLines(lines, md_path)
cat("Done.\n")
