#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  pkgload::load_all("EDI", quiet = TRUE)
})

args = commandArgs(trailingOnly = TRUE)
get_arg = function(name, default, cast = identity) {
  prefix = paste0("--", name, "=")
  hit = args[startsWith(args, prefix)]
  if (length(hit) == 0L) return(default)
  cast(sub(prefix, "", hit[[1L]], fixed = TRUE))
}

seed = get_arg("seed", 20260417L, as.integer)
n = get_arg("n", 200L, as.integer)
p = get_arg("p", 5L, as.integer)
B_fixed = get_arg("B-fixed", 8L, as.integer)
B_kk21 = get_arg("B-kk21", 6L, as.integer)
diagnostics = get_arg("diagnostics", FALSE, function(x) tolower(x) %in% c("1", "true", "t", "yes", "y"))

set.seed(seed)

ci_to_se = function(ci, alpha = 0.05) {
  (unname(ci[2L]) - unname(ci[1L])) / (2 * qnorm(1 - alpha / 2))
}

percent_discordant = function(des) {
  m = des$.__enclos_env__$private$m
  y = des$get_y()
  m = ifelse(is.na(m), 0L, m)
  pair_ids = sort(unique(m[m > 0L]))
  if (length(pair_ids) == 0L) return(NA_real_)
  sums = vapply(pair_ids, function(g) sum(y[m == g]), numeric(1L))
  100 * mean(sums == 1)
}

as_diag_df = function(diag) {
  if (is.null(diag)) {
    diag = list()
  }
  fields = c(
    "matched_component_estimate",
    "matched_component_se",
    "matched_observed_info_beta_T",
    "reservoir_component_estimate",
    "reservoir_component_se",
    "reservoir_observed_info_beta_T",
    "combined_component_estimate",
    "combined_component_se",
    "combined_observed_info_beta_T",
    "combined_minus_matched_observed_info_beta_T",
    "combined_minus_matched_plus_reservoir_observed_info_beta_T"
  )
  out = as.data.frame(as.list(setNames(rep(NA_real_, length(fields)), fields)))
  for (nm in fields) {
    if (!is.null(diag[[nm]])) out[[nm]] = as.numeric(diag[[nm]])
  }
  out
}

bench_pair = function(cls_old, cls_new, des, diagnostics = FALSE) {
  pct_disc = percent_discordant(des)

  old = cls_old$new(des, verbose = FALSE)
  t_old = system.time(ci_old <- old$compute_asymp_confidence_interval())[["elapsed"]]

  new = cls_new$new(des, verbose = FALSE)
  t_new = system.time(ci_new <- new$compute_asymp_confidence_interval())[["elapsed"]]
  diag_df = NULL
  if (isTRUE(diagnostics)) {
    diag_time = system.time(
      diag <- tryCatch(
        new$.__enclos_env__$private$diagnose_clogit_plus_glmm_components(),
        error = function(e) NULL
      )
    )[["elapsed"]]
    diag_df = as_diag_df(diag)
    diag_df$diagnostics_runtime_sec = diag_time
  }

  se_old = ci_to_se(ci_old)
  se_new = ci_to_se(ci_new)
  out = data.frame(
    se_clogit = se_old,
    se_plus = se_new,
    se_diff = se_new - se_old,
    runtime_clogit_sec = t_old,
    runtime_plus_sec = t_new,
    runtime_diff_sec = t_new - t_old,
    pct_discordant_pairs = pct_disc,
    ok = all(is.finite(c(se_old, se_new, t_old, t_new, pct_disc)))
  )
  if (isTRUE(diagnostics)) out = cbind(out, diag_df)
  out
}

make_fixed_binary = function(n = 200L, p = 5L) {
  pair_id = rep(seq_len(n / 2L), each = 2L)
  X_pair = matrix(rnorm((n / 2L) * p), ncol = p)
  X = X_pair[pair_id, , drop = FALSE] + matrix(rnorm(n * p, sd = 0.08), ncol = p)
  colnames(X) = paste0("x", seq_len(p))
  Xdf = as.data.frame(X)
  beta_x = seq(0.25, -0.15, length.out = p)

  des = FixedDesignBinaryMatch$new(n = n, response_type = "incidence", verbose = FALSE)
  des$add_all_subjects_to_experiment(Xdf)
  des$overwrite_all_subject_assignments(rep(c(1, 0), n / 2L))
  des$.__enclos_env__$private$ensure_bms_computed()

  m = des$.__enclos_env__$private$m
  w = des$get_w()
  u = rnorm(max(m), sd = 0.8)
  eta = -0.35 + 0.6 * w + as.vector(as.matrix(Xdf) %*% beta_x) + u[m]
  y = rbinom(n, 1, plogis(eta))
  sums = vapply(split(y, m), sum, numeric(1L))
  if (!any(sums == 1) || !any(sums %in% c(0, 2))) return(NULL)

  des$add_all_subject_responses(y)
  des
}

make_kk21 = function(n = 200L, p = 5L) {
  des = DesignSeqOneByOneKK21stepwise$new(n = n, response_type = "incidence", verbose = FALSE)
  Xmat = matrix(rnorm(n * p), ncol = p)
  colnames(Xmat) = paste0("x", seq_len(p))
  X = as.data.frame(Xmat)
  beta_x = seq(0.25, -0.15, length.out = p)

  for (i in seq_len(n)) {
    des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
  }

  m = des$.__enclos_env__$private$m
  m0 = ifelse(is.na(m), 0L, m)
  w = des$get_w()
  u = rnorm(max(1L, max(m0)), sd = 0.7)
  eta = -0.4 + 0.55 * w + as.vector(as.matrix(X) %*% beta_x) + ifelse(m0 > 0L, u[m0], 0)
  y = rbinom(n, 1, plogis(eta))

  if (!any(m0 > 0L)) return(NULL)
  sums = vapply(split(y[m0 > 0L], m0[m0 > 0L]), sum, numeric(1L))
  if (!any(sums == 1) || !any(sums %in% c(0, 2))) return(NULL)

  des$add_all_subject_responses(y)
  inf = InferenceIncidMultiKKClogitCombinedLikelihood$new(des, verbose = FALSE)
  kk = inf$.__enclos_env__$private$cached_values$KKstats
  if (is.null(kk) || kk$nRT <= 0 || kk$nRC <= 0) return(NULL)
  des
}

run_bench = function(make_fun, old_cls, new_cls, B, label, n, p, diagnostics = FALSE) {
  rows = list()
  i = 1L
  attempts = 0L
  while (i <= B && attempts < B * 12L) {
    attempts = attempts + 1L
    des = make_fun(n = n, p = p)
    if (is.null(des)) next
    out = tryCatch(
      bench_pair(old_cls, new_cls, des, diagnostics = diagnostics),
      error = function(e) {
        out = data.frame(
          se_clogit = NA_real_,
          se_plus = NA_real_,
          se_diff = NA_real_,
          runtime_clogit_sec = NA_real_,
          runtime_plus_sec = NA_real_,
          runtime_diff_sec = NA_real_,
          pct_discordant_pairs = NA_real_,
          ok = FALSE
        )
        if (isTRUE(diagnostics)) out = cbind(out, as_diag_df(NULL), diagnostics_runtime_sec = NA_real_)
        out
      }
    )
    out$rep = i
    out$label = label
    rows[[length(rows) + 1L]] = out
    i = i + 1L
  }
  do.call(rbind, rows)
}

summary_one = function(d) {
  d = d[d$ok, , drop = FALSE]
  data.frame(
    n_ok = nrow(d),
    se_mean_diff = mean(d$se_diff),
    se_median_diff = median(d$se_diff),
    pct_discordant_mean = mean(d$pct_discordant_pairs),
    pct_discordant_median = median(d$pct_discordant_pairs),
    clogit_runtime_mean = mean(d$runtime_clogit_sec),
    clogit_runtime_median = median(d$runtime_clogit_sec),
    plus_runtime_mean = mean(d$runtime_plus_sec),
    plus_runtime_median = median(d$runtime_plus_sec),
    runtime_diff_mean = mean(d$runtime_diff_sec),
    runtime_diff_median = median(d$runtime_diff_sec),
    runtime_ratio_median = median(d$runtime_plus_sec / pmax(d$runtime_clogit_sec, 1e-9))
  )
}

fixed = run_bench(
  make_fixed_binary,
  InferenceIncidMultiKKClogitCombinedLikelihood,
  InferenceIncidMultiKKClogitPlusGLMMCombinedLikelihood,
  B_fixed,
  sprintf("fixed_binary_match_n%d_p%d", n, p),
  n,
  p,
  diagnostics = diagnostics
)

kk21 = run_bench(
  make_kk21,
  InferenceIncidMultiKKClogitCombinedLikelihood,
  InferenceIncidMultiKKClogitPlusGLMMCombinedLikelihood,
  B_kk21,
  sprintf("kk21stepwise_combined_n%d_p%d", n, p),
  n,
  p,
  diagnostics = diagnostics
)

cat("CONFIG\n")
print(data.frame(seed = seed, n = n, p = p, B_fixed = B_fixed, B_kk21 = B_kk21, diagnostics = diagnostics))
cat("\nDETAIL fixed_binary_match\n")
print(fixed)
cat("\nSUMMARY fixed_binary_match\n")
print(summary_one(fixed))
cat("\nDETAIL kk21stepwise_combined\n")
print(kk21)
cat("\nSUMMARY kk21stepwise_combined\n")
print(summary_one(kk21))
