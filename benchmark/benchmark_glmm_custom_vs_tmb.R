devtools::load_all('EDI')
library(glmmTMB)
library(microbenchmark)
library(data.table)

set.seed(42)

# --- 1. Logistic GLMM Benchmark ---
benchmark_logistic_glmm = function(n = 200, n_groups = 50) {
  cat("\n--- Logistic GLMM Benchmark (n =", n, ", groups =", n_groups, ") ---\n")
  
  X = cbind(1, treatment = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n))
  group_id = rep(1:n_groups, each = n/n_groups)
  u = rnorm(n_groups, 0, 0.5)
  eta = X %*% c(-1, 0.5, 0.3, -0.2) + u[group_id]
  y = rbinom(n, 1, plogis(eta))
  
  # glmmTMB
  dat = data.frame(y = y, w = X[,2], x1 = X[,3], x2 = X[,4], g = factor(group_id))
  t1 = system.time({
    mod_tmb = glmmTMB(y ~ w + x1 + x2 + (1|g), data = dat, family = binomial)
  })
  res_tmb = fixef(mod_tmb)$cond
  se_tmb = sqrt(diag(vcov(mod_tmb)$cond))
  
  # Custom
  t2 = system.time({
    mod_custom = fast_logistic_glmm_cpp(X, as.numeric(y), as.integer(group_id), j_T = 1)
  })
  res_custom = mod_custom$b
  se_custom = sqrt(mod_custom$ssq_b_T) # Only for treatment
  
  cat("Estimates (Treatment):\n")
  cat("  glmmTMB: ", res_tmb["w"], " (SE: ", se_tmb["w"], ")\n")
  cat("  Custom:  ", res_custom[2], " (SE: ", mod_custom$ssq_b_T^0.5, ")\n")
  cat("Difference:", abs(res_tmb["w"] - res_custom[2]), "\n")
  
  mb = microbenchmark(
    glmmTMB = glmmTMB(y ~ w + x1 + x2 + (1|g), data = dat, family = binomial),
    Custom  = fast_logistic_glmm_cpp(X, as.numeric(y), as.integer(group_id), j_T = 1),
    times = 20
  )
  print(mb)
}

# --- 2. Ordinal GLMM Benchmark ---
benchmark_ordinal_glmm = function(n = 200, n_groups = 50) {
  cat("\n--- Ordinal GLMM Benchmark (n =", n, ", groups =", n_groups, ") ---\n")
  
  K = 4
  X = cbind(treatment = rbinom(n, 1, 0.5), x1 = rnorm(n))
  group_id = rep(1:n_groups, each = n/n_groups)
  u = rnorm(n_groups, 0, 0.5)
  eta = X %*% c(0.8, -0.4) + u[group_id]
  
  # Sample ordinal responses
  y = integer(n)
  alphas = c(-1, 0, 1.5)
  for(i in 1:n) {
    p = plogis(alphas - eta[i])
    probs = diff(c(0, p, 1))
    y[i] = sample(1:K, 1, prob = probs)
  }
  
  # glmmTMB
  dat = data.frame(y = factor(y, ordered=TRUE), w = X[,1], x1 = X[,2], g = factor(group_id))
  t1 = system.time({
    mod_tmb = glmmTMB(y ~ w + x1 + (1|g), data = dat, family = ordqz())
  })
  # glmmTMB ordqz uses different parameterization, comparison might be tricky 
  # but beta_T should be comparable.
  res_tmb = fixef(mod_tmb)$cond
  
  # Custom
  t2 = system.time({
    mod_custom = fast_ordinal_glmm_cpp(X, as.integer(y), as.integer(group_id), K = K, j_T = 0)
  })
  res_custom = mod_custom$b
  
  cat("Estimates (Treatment):\n")
  cat("  glmmTMB: ", res_tmb["w"], "\n")
  cat("  Custom:  ", res_custom[1], "\n")
  cat("Difference:", abs(res_tmb["w"] - res_custom[1]), "\n")

  mb = microbenchmark(
    glmmTMB = glmmTMB(y ~ w + x1 + (1|g), data = dat, family = ordqz()),
    Custom  = fast_ordinal_glmm_cpp(X, as.integer(y), as.integer(group_id), K = K, j_T = 0),
    times = 10
  )
  print(mb)
}

# --- 3. clogit + GLMM Benchmark ---
benchmark_clogit_plus_glmm = function(n_pairs = 50, n_reservoir = 100) {
  cat("\n--- clogit+GLMM Benchmark (pairs =", n_pairs, ", reservoir =", n_reservoir, ") ---\n")
  
  # Discordant pairs (clogit part)
  X_disc = matrix(rnorm(n_pairs * 2), ncol = 2) # X_diff
  y_disc = rbinom(n_pairs, 1, 0.5)
  
  # Reservoir (GLMM part)
  n_res = n_reservoir
  n_groups_res = n_res / 2
  X_conc = cbind(1, treatment = rbinom(n_res, 1, 0.5), x1 = rnorm(n_res))
  group_conc = rep(1:n_groups_res, each = 2)
  u = rnorm(n_groups_res, 0, 0.5)
  eta = X_conc %*% c(-0.5, 0.7, 0.2) + u[group_conc]
  y_conc = rbinom(n_res, 1, plogis(eta))
  
  start_par = c(-0.5, 0.7, 0.2, -1.0) # b0, bT, b1, log_sigma
  
  # Custom
  t1 = system.time({
    mod_custom = fast_clogit_plus_glmm_cpp(
      X_disc, as.numeric(y_disc), X_conc, as.numeric(y_conc),
      as.integer(group_conc), start_par,
      has_discordant = TRUE, has_concordant = TRUE
    )
  })
  
  cat("Custom Estimate (Treatment):", mod_custom$b[2], " (SE:", sqrt(mod_custom$ssq_b_j), ")\n")
  
  mb = microbenchmark(
    Custom = fast_clogit_plus_glmm_cpp(
      X_disc, as.numeric(y_disc), X_conc, as.numeric(y_conc),
      as.integer(group_conc), start_par,
      has_discordant = TRUE, has_concordant = TRUE
    ),
    times = 20
  )
  print(mb)
}

# Run Benchmarks
benchmark_logistic_glmm(400, 100)
benchmark_ordinal_glmm(400, 100)
benchmark_clogit_plus_glmm(100, 200)
