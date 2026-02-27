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

D = list(
  X = diamonds_subset |>
    model.matrix(price ~ 0 + ., data = _) |>
    apply(2, scale) |> (`/`)(ncol(model.matrix(price ~ 0 + ., data = diamonds_subset))) |>
    as.data.table() |>
    select(where(~ !any(is.na(.)))),
  y_original = finagle(log(diamonds_subset$price))
)

# ---- benchmark settings ----
response_type = "continuous"
nsim          = 351
pval_epsilon  = 0.007
SD_NOISE      = 0.1
beta_T        = 0
design_types  = c("CRD", "iBCRD", "Efron", "KK14", "KK21", "KK21stepwise")
num_cores_vec = 1:6

n      = nrow(D$X)
y_base = D$y_original[[response_type]]

cat(sprintf("Benchmarking: dataset=diamonds  response=%s  n=%d  nsim=%d\n\n",
            response_type, n, nsim))

results = list()

for (design_type in design_types) {
  cat(sprintf("=== Design: %-14s ===\n", design_type))

  # Build the design once per design type (same data for all num_cores)
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
    y_t = y_base[t] + beta_T * (w_t == "T") + rnorm(1, 0, SD_NOISE)
    seq_des_obj$add_subject_response(t, y_t)
  }

  for (nc in num_cores_vec) {
    cat(sprintf("  num_cores = %d  ", nc))

    # 1. Bootstrap p-value (fresh object)
    inf_b  = SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj, num_cores = nc)
    t_boot = system.time(
      pval_boot <- as.numeric(
        inf_b$compute_bootstrap_two_sided_pval(B = nsim, na.rm = TRUE))
    )["elapsed"]

    # 2. Randomization p-value (fresh object, cold cache)
    inf_r  = SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj, num_cores = nc)
    t_rand = system.time(
      pval_rand <- as.numeric(
        inf_r$compute_two_sided_pval_for_treatment_effect_rand(
          nsim_exact_test = nsim, show_progress = FALSE))
    )["elapsed"]

    # 3. Randomization CI (fresh object, cold cache: measures null-dist build + bisection)
    inf_ci = SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj, num_cores = nc)
    t_ci   = system.time(
      ci_rand <- inf_ci$compute_confidence_interval_rand(
        nsim_exact_test = nsim, pval_epsilon = pval_epsilon, show_progress = FALSE)
    )["elapsed"]

    cat(sprintf("boot=%.2fs  rand_pval=%.2fs  rand_ci=%.2fs\n",
                t_boot, t_rand, t_ci))

    results[[length(results) + 1L]] = list(
      design           = design_type,
      num_cores        = nc,
      t_bootstrap_pval = t_boot,
      t_rand_pval      = t_rand,
      t_rand_ci        = t_ci,
      pval_boot        = pval_boot,
      pval_rand        = pval_rand,
      ci_lo            = ci_rand[1],
      ci_hi            = ci_rand[2]
    )
  }
  cat("\n")
}

dt = rbindlist(lapply(results, as.data.table))

# Speedup table relative to num_cores = 1
base = dt[num_cores == 1, .(design, b1 = t_bootstrap_pval, r1 = t_rand_pval, c1 = t_rand_ci)]
dt2  = merge(dt, base, by = "design")
dt2[, `:=`(
  speedup_boot = round(b1 / t_bootstrap_pval, 2),
  speedup_rand = round(r1 / t_rand_pval,      2),
  speedup_ci   = round(c1 / t_rand_ci,        2)
)]

cat("\n============================================================\n")
cat("           Timing (seconds) and Speedup vs 1-core\n")
cat("============================================================\n")
print(
  dt2[, .(design, num_cores,
          boot_s   = round(t_bootstrap_pval, 2), boot_x  = speedup_boot,
          rand_s   = round(t_rand_pval,      2), rand_x  = speedup_rand,
          ci_s     = round(t_rand_ci,        2), ci_x    = speedup_ci)],
  digits = 3
)

out_file = "package_tests/benchmark_parallelism_results.csv"
fwrite(dt, out_file)
cat(sprintf("\nFull results written to %s\n", out_file))
