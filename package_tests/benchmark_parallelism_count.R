rm(list = ls())
set.seed(1)

library(SeqExpMatch)
library(ggplot2)
library(data.table)
library(dplyr)

# ---- diamonds dataset (mirrors _dataset_load.R logic) ----
max_n_dataset = 150
set.seed(1)
diamonds_subset = ggplot2::diamonds |>
  na.omit() |>
  slice_sample(n = max_n_dataset, replace = TRUE) |>
  mutate_if(where(is.factor), as.character)

finagle = function(y_cont) {
  y_scaled = scale(y_cont)
  list(
    continuous = y_scaled,
    incidence  = ifelse(y_cont > median(y_cont), 1, 0),
    proportion = (y_cont - min(y_cont) + 1e-6) / max(y_cont - min(y_cont) + 2e-6),
    count      = round(y_cont - min(y_cont)),
    survival   = y_scaled - min(y_scaled) + 0.1
  )
}

X_mat = diamonds_subset |>
  model.matrix(price ~ 0 + ., data = _) |>
  apply(2, scale) |>
  (`/`)(ncol(model.matrix(price ~ 0 + ., data = diamonds_subset))) |>
  as.data.table() |>
  select(where(~ !any(is.na(.))))

D = list(
  X          = X_mat,
  y_original = finagle(log(diamonds_subset$price))
)

# ---- benchmark settings ----
response_type = "count"
nsim          = 1000
pval_epsilon  = 0.007
SD_NOISE      = 0.1
beta_T        = 0
design_types  = c("CRD", "iBCRD", "Efron", "KK14", "KK21", "KK21stepwise")
num_cores_vec = 1:6

n      = nrow(D$X)
y_base = D$y_original[[response_type]]

cat(sprintf("Benchmarking: dataset=diamonds  response=%s  n=%d  nsim=%d\n",
            response_type, n, nsim))
cat(sprintf("Inference class: SeqDesignInferenceCountMultiNegBinRegr\n\n"))

safe_time = function(expr) {
  tryCatch({
    t = system.time(result <- expr)["elapsed"]
    list(time = t, result = result, error = NULL)
  }, error = function(e) {
    list(time = NA_real_, result = NULL, error = conditionMessage(e))
  })
}

results = list()

for (design_type in design_types) {
  cat(sprintf("=== Design: %-14s ===\n", design_type))

  set.seed(42)
  seq_des_obj = switch(design_type,
    KK21         = SeqDesignKK21$new(        response_type = response_type, n = n),
    KK21stepwise = SeqDesignKK21stepwise$new(response_type = response_type, n = n),
    KK14         = SeqDesignKK14$new(        response_type = response_type, n = n),
    CRD          = SeqDesignCRD$new(         response_type = response_type, n = n),
    Efron        = SeqDesignEfron$new(       response_type = response_type, n = n),
    iBCRD        = SeqDesigniBCRD$new(       response_type = response_type, n = n),
    stop("unknown design: ", design_type)
  )
  set.seed(42)
  for (t in seq_len(n)) {
    w_t = seq_des_obj$add_subject_to_experiment_and_assign(D$X[t, ])
    lambda_t = pmax(.Machine$double.eps, y_base[t] * exp(beta_T * (w_t == "T") + rnorm(1, 0, SD_NOISE)))
    y_t = rpois(1, lambda = lambda_t)
    seq_des_obj$add_subject_response(t, y_t)
  }

  for (nc in num_cores_vec) {
    cat(sprintf("  num_cores = %d\n", nc))

    # 1. Bootstrap p-value
    cat("    bootstrap pval ... ")
    inf_b = SeqDesignInferenceCountMultiNegBinRegr$new(seq_des_obj, num_cores = nc)
    r_boot = safe_time(
      as.numeric(inf_b$compute_bootstrap_two_sided_pval(B = 1000, na.rm = TRUE))
    )
    if (is.null(r_boot$error)) {
      cat(sprintf("%.2fs\n", r_boot$time))
    } else {
      cat(sprintf("ERROR: %s\n", r_boot$error))
    }

    # 2. Randomization p-value (cold cache) — not in standard tests for count, try anyway
    cat("    rand pval       ... ")
    inf_r = SeqDesignInferenceCountMultiNegBinRegr$new(seq_des_obj, num_cores = nc)
    r_rand = safe_time(
      as.numeric(inf_r$compute_two_sided_pval_for_treatment_effect_rand(
        nsim_exact_test = nsim, show_progress = FALSE))
    )
    if (is.null(r_rand$error)) {
      cat(sprintf("%.2fs\n", r_rand$time))
    } else {
      cat(sprintf("ERROR: %s\n", r_rand$error))
    }

    # 3. Randomization CI (cold cache)
    cat("    rand CI         ... ")
    inf_ci = SeqDesignInferenceCountMultiNegBinRegr$new(seq_des_obj, num_cores = nc)
    r_ci = safe_time(
      inf_ci$compute_confidence_interval_rand(
        nsim_exact_test = nsim, pval_epsilon = pval_epsilon, show_progress = FALSE)
    )
    if (is.null(r_ci$error)) {
      cat(sprintf("%.2fs  CI=[%.3f, %.3f]\n", r_ci$time, r_ci$result[1], r_ci$result[2]))
    } else {
      cat(sprintf("ERROR: %s\n", r_ci$error))
    }

    results[[length(results) + 1L]] = list(
      design           = design_type,
      num_cores        = nc,
      t_bootstrap_pval = r_boot$time,
      t_rand_pval      = r_rand$time,
      t_rand_ci        = r_ci$time,
      err_boot         = r_boot$error,
      err_rand         = r_rand$error,
      err_ci           = r_ci$error
    )
  }
  cat("\n")
}

dt = rbindlist(lapply(results, as.data.table))

# Speedup table relative to num_cores = 1
base = dt[num_cores == 1, .(design,
  b1 = t_bootstrap_pval, r1 = t_rand_pval, c1 = t_rand_ci)]
dt2 = merge(dt, base, by = "design")
dt2[, `:=`(
  speedup_boot = round(b1 / t_bootstrap_pval, 2),
  speedup_rand = round(r1 / t_rand_pval,      2),
  speedup_ci   = round(c1 / t_rand_ci,        2)
)]

cat("\n============================================================\n")
cat("  Timing (s) and Speedup vs 1-core  [SeqDesignInferenceCountMultiNegBinRegr]\n")
cat("============================================================\n")
print(
  dt2[, .(design, num_cores,
    boot_s  = round(t_bootstrap_pval, 2), boot_x  = speedup_boot,
    rand_s  = round(t_rand_pval,      2), rand_x  = speedup_rand,
    ci_s    = round(t_rand_ci,        2), ci_x    = speedup_ci)],
  digits = 3
)

out_file = "package_tests/benchmark_parallelism_count_results.csv"
fwrite(dt, out_file)
cat(sprintf("\nFull results written to %s\n", out_file))
