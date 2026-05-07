args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1L]] else file.path(tempdir(), "smart_start_reports")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
	library(EDI)
	library(MASS)
})

n <- 100L
p_total <- 10L
n_cov <- p_total - 2L
alpha_ci <- 0.2
r_rand <- 60L
B_boot <- 60L

bench_reps <- function(expr_fn, reps) {
	times <- numeric(reps)
	values <- vector("list", reps)
	for (i in seq_len(reps)) {
		gc(FALSE)
		times[i] <- system.time(values[[i]] <- expr_fn())[["elapsed"]]
	}
	list(times = times, values = values)
}

median_or_na <- function(x) {
	x <- x[is.finite(x)]
	if (!length(x)) return(NA_real_)
	stats::median(x)
}

mean_or_na <- function(x) {
	x <- x[is.finite(x)]
	if (!length(x)) return(NA_real_)
	mean(x)
}

extract_iterations <- function(obj) {
	if (is.list(obj) && !is.null(obj$iterations)) return(as.numeric(obj$iterations)[1L])
	NA_real_
}

extract_converged <- function(obj) {
	if (is.list(obj) && !is.null(obj$converged)) return(isTRUE(obj$converged))
	if (is.numeric(obj)) return(all(is.finite(obj)))
	NA
}

extract_finite <- function(obj) {
	if (is.numeric(obj)) return(all(is.finite(obj)))
	if (is.list(obj) && !is.null(obj$converged)) {
		return(isTRUE(obj$converged))
	}
	TRUE
}

summarize_fit_policy <- function(model, policy, bench) {
	iters <- vapply(bench$values, extract_iterations, numeric(1))
	converged <- vapply(bench$values, extract_converged, logical(1))
	data.frame(
		model = model,
		policy = policy,
		median_elapsed = stats::median(bench$times),
		min_elapsed = min(bench$times),
		max_elapsed = max(bench$times),
		median_iterations = median_or_na(iters),
		mean_iterations = mean_or_na(iters),
		converged_rate = mean(converged, na.rm = TRUE),
		nonconvergence_rate = mean(!converged, na.rm = TRUE)
	)
}

summarize_reuse_policy <- function(model, reuse, bench) {
	finite_rate <- mean(vapply(bench$values, extract_finite, logical(1)), na.rm = TRUE)
	data.frame(
		model = model,
		reuse = reuse,
		median_elapsed = stats::median(bench$times),
		min_elapsed = min(bench$times),
		max_elapsed = max(bench$times),
		finite_rate = finite_rate
	)
}

add_gain_column <- function(df, group_col) {
	df$elapsed_gain_vs_off_or_legacy_pct <- NA_real_
	keys <- unique(df[[group_col]])
	state_col <- if ("policy" %in% names(df)) "policy" else "reuse"
	for (k in keys) {
		idx_on <- which(df[[group_col]] == k & df[[state_col]] %in% c("smart", "on"))
		idx_off <- which(df[[group_col]] == k & df[[state_col]] %in% c("legacy", "off"))
		if (length(idx_on) == 1L && length(idx_off) == 1L) {
			base <- df$median_elapsed[idx_off]
			if (is.finite(base) && base > 0) {
				df$elapsed_gain_vs_off_or_legacy_pct[idx_on] <- 100 * (base - df$median_elapsed[idx_on]) / base
				df$elapsed_gain_vs_off_or_legacy_pct[idx_off] <- 0
			}
		}
	}
	df
}

make_fixed_design <- function(response_type, X_cov, w, y, dead = NULL) {
	des <- FixedDesign$new(n = length(y), response_type = response_type, verbose = FALSE)
	des$add_all_subjects_to_experiment(as.data.frame(X_cov))
	des$overwrite_all_subject_assignments(w)
	if (is.null(dead)) des$add_all_subject_responses(y) else des$add_all_subject_responses(y, dead)
	des
}

make_seq_count_design <- function(X_cov, y) {
	des <- DesignSeqOneByOneBernoulli$new(n = nrow(X_cov), response_type = "count", verbose = FALSE)
	for (i in seq_len(nrow(X_cov))) {
		des$add_one_subject_to_experiment_and_assign(as.data.frame(X_cov[i, , drop = FALSE]))
	}
	for (i in seq_along(y)) {
		des$add_one_subject_response(i, y[i], 1)
	}
	des
}

with_reuse_settings <- function(inf_obj, fit_warm = TRUE, null_warm = TRUE, reusable_worker = TRUE) {
	priv <- inf_obj$.__enclos_env__$private
	priv$fit_warm_start_enabled <- fit_warm
	priv$null_fit_warm_start_enabled <- null_warm
	priv$reusable_bootstrap_worker_enabled <- reusable_worker
	priv$clear_fit_warm_start()
	priv$clear_likelihood_null_warm_cache()
	invisible(inf_obj)
}

set.seed(20260506)
X_cov <- replicate(n_cov, rnorm(n))
colnames(X_cov) <- paste0("x", seq_len(n_cov))
w <- rep(c(1, 0), length.out = n)
X_fit <- cbind(1, treatment = w, X_cov)

beta_cov <- seq(0.05, 0.4, length.out = n_cov)
linpred_logit <- -0.2 + 0.6 * w + drop(as.matrix(X_cov) %*% beta_cov)
linpred_pois <- 0.15 + 0.25 * w + drop(as.matrix(X_cov) %*% rev(beta_cov)) / 2

y_logit <- rbinom(n, 1, plogis(linpred_logit))
y_pois <- rpois(n, lambda = exp(linpred_pois))
y_nb <- MASS::rnegbin(n, mu = exp(0.1 + 0.25 * w + drop(as.matrix(X_cov) %*% beta_cov) / 3), theta = 2.5)
y_surv <- exp(0.8 + 0.3 * w + drop(as.matrix(X_cov) %*% beta_cov) / 5 + rnorm(n, sd = 0.25))
dead <- rbinom(n, 1, 0.85)
latent_ord <- 0.4 * w + drop(as.matrix(X_cov) %*% beta_cov) / 4 + rnorm(n)
y_ord <- as.numeric(cut(latent_ord, breaks = c(-Inf, -0.5, 0.4, Inf), labels = FALSE))

init_report <- do.call(rbind, list(
	summarize_fit_policy("logistic", "smart", bench_reps(function() EDI:::fast_logistic_regression_weighted_cpp(X_fit, y_logit, rep(1, n), smart_start = TRUE), 5L)),
	summarize_fit_policy("logistic", "legacy", bench_reps(function() EDI:::fast_logistic_regression_weighted_cpp(X_fit, y_logit, rep(1, n), smart_start = FALSE), 5L)),
	summarize_fit_policy("poisson", "smart", bench_reps(function() EDI:::fast_poisson_regression_cpp(X_fit, y_pois, smart_start = TRUE), 5L)),
	summarize_fit_policy("poisson", "legacy", bench_reps(function() EDI:::fast_poisson_regression_cpp(X_fit, y_pois, smart_start = FALSE), 5L)),
	summarize_fit_policy("negbin", "smart", bench_reps(function() EDI:::fast_neg_bin_cpp(X_fit, as.integer(y_nb), smart_start = TRUE), 5L)),
	summarize_fit_policy("negbin", "legacy", bench_reps(function() EDI:::fast_neg_bin_cpp(X_fit, as.integer(y_nb), smart_start = FALSE), 5L)),
	summarize_fit_policy("weibull", "smart", bench_reps(function() EDI:::fast_weibull_regression_cpp(X_fit, y_surv, dead, smart_start = TRUE, estimate_only = TRUE), 5L)),
	summarize_fit_policy("weibull", "legacy", bench_reps(function() EDI:::fast_weibull_regression_cpp(X_fit, y_surv, dead, smart_start = FALSE, estimate_only = TRUE), 5L)),
	summarize_fit_policy("ordinal_logit", "smart", bench_reps(function() EDI:::fast_ordinal_regression_cpp(X_fit[, -1, drop = FALSE], y_ord, smart_start = TRUE), 5L)),
	summarize_fit_policy("ordinal_logit", "legacy", bench_reps(function() EDI:::fast_ordinal_regression_cpp(X_fit[, -1, drop = FALSE], y_ord, smart_start = FALSE), 5L))
))
init_report <- add_gain_column(init_report, "model")

des_logit <- make_fixed_design("incidence", X_cov, w, y_logit)
des_pois <- make_fixed_design("count", X_cov, w, y_pois)
des_nb <- make_fixed_design("count", X_cov, w, y_nb)
des_weib <- make_fixed_design("survival", X_cov, w, y_surv, dead)
des_ord <- make_fixed_design("ordinal", X_cov, w, y_ord)

ci_models <- list(
	logistic = function() InferenceIncidLogRegr$new(des_logit, verbose = FALSE, smart_default = TRUE),
	poisson = function() InferenceCountPoisson$new(des_pois, verbose = FALSE, smart_default = TRUE),
	negbin = function() InferenceCountNegBin$new(des_nb, verbose = FALSE, smart_default = TRUE),
	weibull = function() InferenceSurvivalWeibullRegr$new(des_weib, verbose = FALSE, smart_default = TRUE),
	ordinal_logit = function() InferenceOrdinalPropOddsRegr$new(des_ord, verbose = FALSE, smart_default = TRUE)
)

bench_ci_type <- function(testing_type) {
	rows <- list()
	for (model in names(ci_models)) {
		for (reuse in c("on", "off")) {
			inf <- ci_models[[model]]()
			inf$set_testing_type(testing_type)
			if (reuse == "on") {
				with_reuse_settings(inf, fit_warm = TRUE, null_warm = TRUE, reusable_worker = TRUE)
			} else {
				with_reuse_settings(inf, fit_warm = FALSE, null_warm = FALSE, reusable_worker = TRUE)
			}
			bench <- bench_reps(function() inf$compute_asymp_confidence_interval(alpha = alpha_ci), 3L)
			rows[[length(rows) + 1L]] <- summarize_reuse_policy(model, reuse, bench)
		}
	}
	out <- do.call(rbind, rows)
	add_gain_column(out, "model")
}

score_ci_report <- bench_ci_type("score")
lik_ratio_ci_report <- bench_ci_type("lik_ratio")

X_cov_rand <- X_cov[1:n, , drop = FALSE]
y_pois_rand <- y_pois
des_rand <- make_seq_count_design(X_cov_rand, y_pois_rand)

rand_rows <- list()
for (reuse in c("on", "off")) {
	inf <- InferenceCountPoisson$new(des_rand, verbose = FALSE, smart_default = TRUE)
	if (reuse == "on") {
		with_reuse_settings(inf, fit_warm = TRUE, null_warm = TRUE, reusable_worker = TRUE)
		ctrl <- list(fit_warm_start_enable = TRUE, fit_reuse_factorizations = TRUE)
	} else {
		with_reuse_settings(inf, fit_warm = FALSE, null_warm = TRUE, reusable_worker = FALSE)
		ctrl <- list(fit_warm_start_enable = FALSE, fit_reuse_factorizations = FALSE)
	}
	bench <- bench_reps(function() {
		inf$compute_rand_confidence_interval(alpha = alpha_ci, r = r_rand, pval_epsilon = 0.05, show_progress = FALSE, ci_search_control = ctrl)
	}, 2L)
	rand_rows[[length(rand_rows) + 1L]] <- summarize_reuse_policy("poisson", reuse, bench)
}
randomization_report <- add_gain_column(do.call(rbind, rand_rows), "model")

boot_rows <- list()
for (reuse in c("on", "off")) {
	inf <- InferenceCountPoisson$new(des_rand, verbose = FALSE, smart_default = TRUE)
	if (reuse == "on") {
		with_reuse_settings(inf, fit_warm = TRUE, null_warm = TRUE, reusable_worker = TRUE)
	} else {
		with_reuse_settings(inf, fit_warm = FALSE, null_warm = TRUE, reusable_worker = FALSE)
	}
	bench <- bench_reps(function() {
		inf$compute_bootstrap_confidence_interval(alpha = alpha_ci, B = B_boot, show_progress = FALSE)
	}, 2L)
	boot_rows[[length(boot_rows) + 1L]] <- summarize_reuse_policy("poisson", reuse, bench)
}
bootstrap_report <- add_gain_column(do.call(rbind, boot_rows), "model")

utils::write.csv(init_report, file.path(out_dir, "initialization_report.csv"), row.names = FALSE)
utils::write.csv(score_ci_report, file.path(out_dir, "ci_inversion_score_report.csv"), row.names = FALSE)
utils::write.csv(lik_ratio_ci_report, file.path(out_dir, "ci_inversion_lik_ratio_report.csv"), row.names = FALSE)
utils::write.csv(randomization_report, file.path(out_dir, "randomization_ci_reuse_report.csv"), row.names = FALSE)
utils::write.csv(bootstrap_report, file.path(out_dir, "bootstrap_reuse_report.csv"), row.names = FALSE)

summary_lines <- c(
	"# Smart Start Benchmark Reports",
	"",
	sprintf("- n = %d", n),
	sprintf("- p = %d (intercept + treatment + %d covariates)", p_total, n_cov),
	sprintf("- direct-fit reps = %d", 5L),
	sprintf("- CI inversion reps = %d", 3L),
	sprintf("- randomization/bootstrap reps = %d", 2L),
	"",
	"Generated files:",
	"- initialization_report.csv",
	"- ci_inversion_score_report.csv",
	"- ci_inversion_lik_ratio_report.csv",
	"- randomization_ci_reuse_report.csv",
	"- bootstrap_reuse_report.csv"
)
writeLines(summary_lines, file.path(out_dir, "README.md"))
cat(out_dir, "\n")
