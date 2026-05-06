args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list(lib = NULL, out = NULL, reps = 5L)
  for (arg in args) {
    if (startsWith(arg, "--lib=")) out$lib <- sub("^--lib=", "", arg)
    if (startsWith(arg, "--out=")) out$out <- sub("^--out=", "", arg)
    if (startsWith(arg, "--reps=")) out$reps <- as.integer(sub("^--reps=", "", arg))
  }
  if (is.null(out$lib) || is.null(out$out)) {
    stop("usage: Rscript scripts/benchmark_simd_matrix_worker.R --lib=/path/to/lib --out=/path/to/out.csv [--reps=5]")
  }
  out
}

cfg <- parse_args(args)
.libPaths(c(cfg$lib, .libPaths()))
suppressPackageStartupMessages(library(EDI, lib.loc = cfg$lib))

timed_eval <- function(fun, reps) {
  invisible(fun())
  times <- numeric(reps)
  for (i in seq_len(reps)) {
    gc(FALSE)
    times[i] <- system.time(invisible(fun()))[["elapsed"]]
  }
  times
}

make_design_matrix <- function(n, p_extra, seed) {
  set.seed(seed)
  X_cov <- matrix(rnorm(n * p_extra), nrow = n, ncol = p_extra)
  cbind(1, rbinom(n, 1, 0.5), X_cov)
}

scenario_poisson_glmm <- function() {
  n_groups <- 120L
  group_size <- 5L
  n <- n_groups * group_size
  group_id <- rep(seq_len(n_groups), each = group_size)
  X <- make_design_matrix(n, 4L, 1001L)
  beta <- c(0.2, 0.15, -0.1, 0.08, 0.05, -0.04)
  u <- rep(rnorm(n_groups, sd = 0.45), each = group_size)
  y <- rpois(n, lambda = exp(drop(X %*% beta) + u))
  function() EDI:::fast_poisson_glmm_cpp(X, y, group_id, j_T = 2L, estimate_only = TRUE, n_gh = 12L, maxit = 120L, eps_g = 1e-5, optimization_alg = "lbfgs")
}

scenario_logistic_glmm <- function() {
  n_groups <- 140L
  group_size <- 4L
  n <- n_groups * group_size
  group_id <- rep(seq_len(n_groups), each = group_size)
  X <- make_design_matrix(n, 4L, 1002L)
  beta <- c(-0.3, 0.4, -0.2, 0.1, 0.15, -0.05)
  u <- rep(rnorm(n_groups, sd = 0.7), each = group_size)
  p <- 1 / (1 + exp(-(drop(X %*% beta) + u)))
  y <- rbinom(n, 1, p)
  function() EDI:::fast_logistic_glmm_cpp(X, y, group_id, j_T = 2L, estimate_only = TRUE, n_gh = 12L, maxit = 120L, eps_g = 1e-5, optimization_alg = "lbfgs")
}

scenario_clogit_plus_glmm <- function() {
  set.seed(1003L)
  n_disc <- 320L
  X_disc <- cbind(rbinom(n_disc, 1, 0.5), matrix(rnorm(n_disc * 3L), nrow = n_disc))
  beta_disc <- c(0.35, -0.15, 0.1, 0.05)
  p_disc <- 1 / (1 + exp(-drop(X_disc %*% beta_disc)))
  y_disc <- rbinom(n_disc, 1, p_disc)

  n_groups <- 90L
  group_size <- 4L
  n_conc <- n_groups * group_size
  group_conc <- rep(seq_len(n_groups), each = group_size)
  X_conc <- make_design_matrix(n_conc, 3L, 1004L)
  beta_conc <- c(-0.1, 0.25, -0.08, 0.06, 0.04)
  u <- rep(rnorm(n_groups, sd = 0.6), each = group_size)
  p_conc <- 1 / (1 + exp(-(drop(X_conc %*% beta_conc) + u)))
  y_conc <- rbinom(n_conc, 1, p_conc)
  start <- rep(0, ncol(X_conc) + 1L)
  function() EDI:::fast_clogit_plus_glmm_cpp(X_disc, y_disc, X_conc, y_conc, group_conc, start, has_discordant = TRUE, has_concordant = TRUE, estimate_only = TRUE, maxit = 120L, eps_g = 1e-5, optimization_alg = "lbfgs")
}

scenario_adjacent_category_logit <- function() {
  set.seed(1005L)
  n <- 900L
  X <- matrix(rnorm(n * 4L), nrow = n, ncol = 4L)
  eta <- 0.6 * X[, 1] - 0.35 * X[, 2] + 0.15 * X[, 3]
  alpha <- c(0.7, 0.2, -0.4)
  scores <- cbind(alpha[1] + alpha[2] + alpha[3] - 3 * eta,
                  alpha[2] + alpha[3] - 2 * eta,
                  alpha[3] - eta,
                  0)
  probs <- exp(scores - apply(scores, 1L, max))
  probs <- probs / rowSums(probs)
  y <- apply(probs, 1L, function(pr) sample.int(4L, size = 1L, prob = pr))
  function() EDI:::fast_adjacent_category_logit_cpp(X, y, maxit = 80L, tol = 1e-6, optimization_alg = "newton_raphson")
}

scenario_ordinal_clmm <- function() {
  set.seed(1006L)
  n_groups <- 100L
  group_size <- 4L
  n <- n_groups * group_size
  K <- 4L
  group_id <- rep(seq_len(n_groups), each = group_size)
  X <- matrix(rnorm(n * 3L), nrow = n, ncol = 3L)
  beta <- c(0.4, -0.2, 0.1)
  u <- rep(rnorm(n_groups, sd = 0.5), each = group_size)
  latent <- drop(X %*% beta) + u + rlogis(n)
  y <- cut(latent, breaks = c(-Inf, -0.5, 0.5, 1.4, Inf), labels = FALSE)
  function() EDI:::fast_ordinal_clmm_cpp(X, y, group_id, K = K, j_T = 2L, link = "logit", estimate_only = TRUE, n_gh = 10L, maxit = 120L, eps_g = 1e-5, optimization_alg = "lbfgs")
}

scenario_weibull_frailty <- function() {
  set.seed(1007L)
  n_groups <- 120L
  group_size <- 4L
  n <- n_groups * group_size
  group_id <- rep(seq_len(n_groups), each = group_size)
  X <- make_design_matrix(n, 3L, 1008L)
  beta <- c(2.2, -0.15, 0.2, -0.05, 0.08)
  sigma_eps <- 0.7
  sigma_u <- 0.45
  u <- rep(rnorm(n_groups, sd = sigma_u), each = group_size)
  log_t <- drop(X %*% beta) + u + sigma_eps * rlogis(n)
  y <- pmax(exp(log_t), 1e-6)
  dead <- rbinom(n, 1, 0.7)
  function() EDI:::fast_weibull_frailty_cpp(y, dead, X, group_id, estimate_only = TRUE, n_gh = 10L, maxit = 120L, eps_g = 1e-5, optimization_alg = "lbfgs")
}

scenarios <- list(
  poisson_glmm = scenario_poisson_glmm(),
  logistic_glmm = scenario_logistic_glmm(),
  clogit_plus_glmm = scenario_clogit_plus_glmm(),
  adjacent_category_logit = scenario_adjacent_category_logit(),
  ordinal_clmm = scenario_ordinal_clmm(),
  weibull_frailty = scenario_weibull_frailty()
)

rows <- list()
for (nm in names(scenarios)) {
  times <- timed_eval(scenarios[[nm]], cfg$reps)
  rows[[length(rows) + 1L]] <- data.frame(
    kernel = nm,
    rep = seq_along(times),
    elapsed = times,
    stringsAsFactors = FALSE
  )
}

dir.create(dirname(cfg$out), recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, rows), cfg$out, row.names = FALSE)
