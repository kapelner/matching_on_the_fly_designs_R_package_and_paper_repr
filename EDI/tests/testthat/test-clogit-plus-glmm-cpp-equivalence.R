library(testthat)
library(EDI)

# Tests comparing fast_clogit_plus_glmm_cpp against canonical packages.
#
# fast_clogit_plus_glmm_cpp has two components controlled by has_discordant /
# has_concordant:
#
#   concordant-only  (has_concordant=TRUE,  has_discordant=FALSE)
#     -- random-intercept logistic GLMM via 20-node GH quadrature
#     -- canonical reference: lme4::glmer (nAGQ=20)
#
#   discordant-only  (has_discordant=TRUE,  has_concordant=FALSE)
#     -- conditional logistic regression on matched-pair differences
#     -- canonical reference: survival::clogit
#
#   combined         (both TRUE)
#     -- no single canonical package; test convergence + score=0 fingerprint
#
# Parameter layout for concordant path:
#   params = [intercept, beta_T, ..., log_sigma]   (length p_conc + 1)
#
# Parameter layout for discordant-only path:
#   params = [beta_T, beta_x_diff, ...]            (length p_disc = ncol(X_disc))
#
# Tolerance notes:
#   concordant vs lme4 (nAGQ=20): betas to ~1e-4, log_sigma to ~1e-3
#   discordant vs survival::clogit: betas to ~1e-5

# ---------------------------------------------------------------------------
# Helper: empty design-matrix placeholder for the unused component
# ---------------------------------------------------------------------------
empty_X <- function() matrix(numeric(0), nrow = 0L, ncol = 1L)
empty_y <- function() numeric(0L)
empty_g <- function() integer(0L)

# ===========================================================================
# 1. Concordant-only (GLMM component) versus lme4::glmer
# ===========================================================================

test_that("fast_clogit_plus_glmm_cpp concordant-only matches lme4::glmer (nAGQ=20)", {
	skip_if_not_installed("lme4")
	set.seed(123)

	n_pairs <- 35L
	group_id <- rep(seq_len(n_pairs), each = 2L)
	w  <- rep(c(0L, 1L), n_pairs)
	x1 <- rnorm(2L * n_pairs)
	re <- rnorm(n_pairs, sd = 0.8)[group_id]   # sd=0.8 keeps MLE away from the -Inf boundary
	y  <- rbinom(2L * n_pairs, 1L, plogis(-0.2 + 0.6 * w + 0.4 * x1 + re))

	X_conc <- cbind(1, w, x1)

	res_cpp <- EDI:::fast_clogit_plus_glmm_cpp(
		X_disc        = empty_X(),
		y_disc        = empty_y(),
		X_conc        = X_conc,
		y_conc        = as.numeric(y),
		group_conc    = as.integer(group_id),
		has_discordant = FALSE,
		has_concordant = TRUE
	)

	dat      <- data.frame(y = y, w = w, x1 = x1, grp = factor(group_id))
	res_lme4 <- lme4::glmer(y ~ w + x1 + (1 | grp), data = dat,
	                        family = stats::binomial, nAGQ = 20L)

	beta_cpp  <- as.numeric(res_cpp$params)[seq_len(3L)]   # intercept, w, x1
	beta_lme4 <- as.numeric(lme4::fixef(res_lme4))
	lsig_lme4 <- log(sqrt(as.numeric(lme4::VarCorr(res_lme4)$grp[1L])))
	lsig_cpp  <- as.numeric(res_cpp$params)[4L]

	expect_true(res_cpp$converged, label = "concordant GLMM converged")
	expect_equal(beta_cpp, beta_lme4,  tolerance = 1e-3, label = "concordant GLMM betas vs lme4")
	expect_equal(lsig_cpp, lsig_lme4,  tolerance = 2e-3, label = "concordant GLMM log_sigma vs lme4")
})

# ===========================================================================
# 2. Discordant-only (clogit component) versus survival::clogit
# ===========================================================================

test_that("fast_clogit_plus_glmm_cpp discordant-only matches survival::clogit (no covariates)", {
	skip_if_not_installed("survival")
	# survival::clogit internally calls bare coxph; attaching ensures it is found
	suppressPackageStartupMessages(library(survival, quietly = TRUE))
	set.seed(11)

	n_pairs <- 80L
	beta_T  <- 1.2
	p_T <- plogis(beta_T)
	p_C <- 0.5

	y_T <- rbinom(n_pairs, 1L, p_T)
	y_C <- rbinom(n_pairs, 1L, p_C)

	disc_idx <- which(y_T != y_C)
	n_disc   <- length(disc_idx)
	skip_if(n_disc < 10L, "too few discordant pairs; regenerate")

	y_disc_vec <- as.numeric((y_T[disc_idx] - y_C[disc_idx] + 1) / 2)  # 1 if T won, 0 if C won
	X_disc_mat <- matrix(1, nrow = n_disc, ncol = 1L)
	colnames(X_disc_mat) <- "treatment"

	res_cpp <- EDI:::fast_clogit_plus_glmm_cpp(
		X_disc        = X_disc_mat,
		y_disc        = y_disc_vec,
		X_conc        = matrix(0, nrow = 0L, ncol = 2L),
		y_conc        = empty_y(),
		group_conc    = empty_g(),
		has_discordant = TRUE,
		has_concordant = FALSE
	)

	# Build survival::clogit dataset: 2 rows per discordant pair
	strata_id  <- rep(seq_len(n_disc), each = 2L)
	trt_vec    <- rep(c(1L, 0L), n_disc)                 # T row then C row
	status_vec <- as.integer(c(rbind(y_disc_vec, 1 - y_disc_vec)))

	surv_dat   <- data.frame(status = status_vec, trt = trt_vec,
	                         pair = strata_id)
	res_surv   <- survival::clogit(status ~ trt + survival::strata(pair),
	                               data = surv_dat)

	beta_cpp  <- as.numeric(res_cpp$params)[1L]
	beta_surv <- as.numeric(stats::coef(res_surv)["trt"])

	expect_true(res_cpp$converged, label = "discordant clogit converged")
	expect_equal(beta_cpp, beta_surv, tolerance = 1e-5,
	             label = "discordant clogit beta_T vs survival::clogit")
})

test_that("fast_clogit_plus_glmm_cpp discordant-only matches survival::clogit (with covariate diffs)", {
	skip_if_not_installed("survival")
	suppressPackageStartupMessages(library(survival, quietly = TRUE))
	set.seed(13)

	n_pairs <- 100L
	beta_T  <- 0.9
	beta_x  <- 0.6
	x_T <- rnorm(n_pairs)
	x_C <- rnorm(n_pairs)
	y_T <- rbinom(n_pairs, 1L, plogis(beta_T + beta_x * x_T))
	y_C <- rbinom(n_pairs, 1L, plogis(         beta_x * x_C))

	disc_idx <- which(y_T != y_C)
	n_disc   <- length(disc_idx)
	skip_if(n_disc < 10L, "too few discordant pairs; regenerate")

	y_disc_vec  <- as.numeric((y_T[disc_idx] - y_C[disc_idx] + 1) / 2)
	delta_x     <- x_T[disc_idx] - x_C[disc_idx]          # T minus C difference
	X_disc_mat  <- cbind(treatment = 1, x_diff = delta_x)

	res_cpp <- EDI:::fast_clogit_plus_glmm_cpp(
		X_disc        = X_disc_mat,
		y_disc        = y_disc_vec,
		X_conc        = matrix(0, nrow = 0L, ncol = 2L),
		y_conc        = empty_y(),
		group_conc    = empty_g(),
		has_discordant = TRUE,
		has_concordant = FALSE
	)

	# Build survival::clogit dataset with actual covariate values per individual
	strata_id  <- rep(seq_len(n_disc), each = 2L)
	trt_vec    <- rep(c(1L, 0L), n_disc)
	x_vec      <- c(rbind(x_T[disc_idx], x_C[disc_idx]))
	status_vec <- as.integer(c(rbind(y_disc_vec, 1L - y_disc_vec)))

	surv_dat  <- data.frame(status = status_vec, trt = trt_vec,
	                        x = x_vec, pair = strata_id)
	res_surv  <- survival::clogit(status ~ trt + x + survival::strata(pair),
	                              data = surv_dat)

	beta_cpp  <- as.numeric(res_cpp$params)
	beta_surv <- as.numeric(stats::coef(res_surv)[c("trt", "x")])

	expect_true(res_cpp$converged, label = "discordant clogit (covariate) converged")
	expect_equal(beta_cpp, beta_surv, tolerance = 1e-5,
	             label = "discordant clogit betas vs survival::clogit (with covariate)")
})

# ===========================================================================
# 3. Combined model: convergence and score-at-zero fingerprint
# ===========================================================================

test_that("fast_clogit_plus_glmm_cpp combined model converges and score is near zero at optimum", {
	set.seed(17)

	# Concordant pairs (same outcome in both arms)
	n_conc <- 20L
	group_conc <- rep(seq_len(n_conc), each = 2L)
	w_conc  <- rep(c(0L, 1L), n_conc)
	x_conc  <- rnorm(2L * n_conc)
	re_conc <- rnorm(n_conc, sd = 0.4)[group_conc]
	y_conc_prob <- plogis(-0.1 + 0.5 * w_conc + 0.3 * x_conc + re_conc)
	y_conc_full <- rbinom(2L * n_conc, 1L, y_conc_prob)
	# Force concordant: keep only pairs where T and C had same outcome
	pair_y <- matrix(y_conc_full, nrow = 2L)
	conc_pairs <- which(pair_y[1L,] == pair_y[2L,])
	n_c <- length(conc_pairs)
	skip_if(n_c < 5L, "too few concordant pairs; regenerate")
	idx_c    <- as.vector(sapply(conc_pairs, function(p) c(2*p-1, 2*p)))
	X_conc   <- cbind(1, w_conc[idx_c], x_conc[idx_c])
	y_conc   <- as.numeric(y_conc_full[idx_c])
	g_conc   <- as.integer(group_conc[idx_c])

	# Discordant pairs (different outcome in T vs C)
	n_disc_gen <- 30L
	x_d  <- rnorm(n_disc_gen)
	y_T  <- rbinom(n_disc_gen, 1L, plogis(0.8 + 0.5 * x_d))
	y_C  <- rbinom(n_disc_gen, 1L, plogis(       0.5 * x_d))
	d_idx <- which(y_T != y_C)
	n_d   <- length(d_idx)
	skip_if(n_d < 5L, "too few discordant pairs; regenerate")
	y_disc   <- as.numeric((y_T[d_idx] - y_C[d_idx] + 1) / 2)
	X_disc   <- matrix(1, nrow = n_d, ncol = 1L)  # treatment indicator only

	res <- EDI:::fast_clogit_plus_glmm_cpp(
		X_disc        = X_disc,
		y_disc        = y_disc,
		X_conc        = X_conc,
		y_conc        = y_conc,
		group_conc    = g_conc,
		has_discordant = TRUE,
		has_concordant = TRUE
	)

	expect_true(res$converged, label = "combined model converged")
	expect_true(is.finite(res$beta_T), label = "combined model beta_T is finite")

	# Score at the optimum must be near zero for all free parameters
	score <- EDI:::get_clogit_plus_glmm_score_cpp(
		X_disc, y_disc, X_conc, y_conc, g_conc,
		as.numeric(res$params),
		has_discordant = TRUE, has_concordant = TRUE
	)
	expect_true(all(abs(score) < 1e-3),
	            label = "score at optimum is near zero (all components)")
})

# ===========================================================================
# 4. Structural / smoke tests
# ===========================================================================

test_that("fast_clogit_plus_glmm_cpp concordant-only returns correct parameter length", {
	set.seed(21)
	n  <- 20L
	g  <- rep(seq_len(n / 2L), each = 2L)
	X  <- cbind(1, rbinom(n, 1L, 0.5), rnorm(n))
	y  <- rbinom(n, 1L, 0.5)

	res <- EDI:::fast_clogit_plus_glmm_cpp(
		X_disc = empty_X(), y_disc = empty_y(),
		X_conc = X, y_conc = as.numeric(y),
		group_conc = as.integer(g),
		has_discordant = FALSE, has_concordant = TRUE
	)
	# p_conc = 3, so params = [intercept, beta_T, beta_x, log_sigma]
	expect_length(res$params, ncol(X) + 1L)
})

test_that("fast_clogit_plus_glmm_cpp discordant-only returns correct parameter length", {
	set.seed(23)
	n_d  <- 15L
	X_d  <- cbind(1, rnorm(n_d))
	y_d  <- rbinom(n_d, 1L, 0.5)

	res <- EDI:::fast_clogit_plus_glmm_cpp(
		X_disc = X_d, y_disc = as.numeric(y_d),
		X_conc = matrix(0, 0L, 2L), y_conc = empty_y(),
		group_conc = empty_g(),
		has_discordant = TRUE, has_concordant = FALSE
	)
	expect_length(res$params, ncol(X_d))
})
