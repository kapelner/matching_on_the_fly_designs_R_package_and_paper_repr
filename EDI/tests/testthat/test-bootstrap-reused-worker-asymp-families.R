compare_bootstrap_fast_slow_asymp <- function(fast_inf, slow_inf, B = 5L, seed = 1L, tolerance = 1e-10){
	fast_inf$num_cores = 1L
	slow_inf$num_cores = 1L
	set.seed(seed)
	fast_boot = fast_inf$approximate_bootstrap_distribution_beta_hat_T(B = B, show_progress = FALSE)
	set.seed(seed)
	slow_boot = slow_inf$approximate_bootstrap_distribution_beta_hat_T(B = B, show_progress = FALSE)
	expect_equal(unname(fast_boot), unname(slow_boot), tolerance = tolerance)
}

make_fixed_design <- function(response_type, X, y_fun, dead_fun = NULL){
	des = FixedDesignBernoulli$new(n = nrow(X), response_type = response_type, verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	w = des$get_w()
	y = y_fun(w, X)
	if (is.null(dead_fun)) {
		des$add_all_subject_responses(y)
	} else {
		des$add_all_subject_responses(y, deads = dead_fun(w, X, y))
	}
	des
}

make_fixed_blocked_cluster_design <- function(X, y_fun, cluster_size = 2L){
	X_design = as.data.frame(X)
	strata_col = "block"
	cluster_col = ".assignment_only_cluster_id"
	cluster_ids = integer(nrow(X_design))
	next_cluster_id = 1L

	for (block in unique(X_design[[strata_col]])) {
		idx = which(X_design[[strata_col]] == block)
		cluster_ids[idx] = ((seq_along(idx) - 1L) %/% cluster_size) + next_cluster_id
		next_cluster_id = max(cluster_ids[idx]) + 1L
	}

	X_design[[cluster_col]] = factor(cluster_ids)
	des = FixedDesignBlockedCluster$new(
		n = nrow(X_design),
		strata_cols = strata_col,
		cluster_col = cluster_col,
		response_type = "continuous",
		verbose = FALSE
	)
	des$add_all_subjects_to_experiment(X_design)
	des$assign_w_to_all_subjects()
	w = des$get_w()
	des$add_all_subject_responses(y_fun(w, X_design))
	des
}

test_that("incidence and ordinal g-computation reusable workers match generic bootstrap", {
	SlowInferenceIncidUnivGCompRiskDiff = R6::R6Class(
		"SlowInferenceIncidUnivGCompRiskDiff",
		inherit = InferenceIncidUnivGCompRiskDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceIncidMultiGCompRiskDiff = R6::R6Class(
		"SlowInferenceIncidMultiGCompRiskDiff",
		inherit = InferenceIncidMultiGCompRiskDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceIncidUnivGCompRiskRatio = R6::R6Class(
		"SlowInferenceIncidUnivGCompRiskRatio",
		inherit = InferenceIncidUnivGCompRiskRatio,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceIncidMultiGCompRiskRatio = R6::R6Class(
		"SlowInferenceIncidMultiGCompRiskRatio",
		inherit = InferenceIncidMultiGCompRiskRatio,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceOrdinalUniGCompMeanDiff = R6::R6Class(
		"SlowInferenceOrdinalUniGCompMeanDiff",
		inherit = InferenceOrdinalUniGCompMeanDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceOrdinalMultiGCompMeanDiff = R6::R6Class(
		"SlowInferenceOrdinalMultiGCompMeanDiff",
		inherit = InferenceOrdinalMultiGCompMeanDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)

	set.seed(20260415)
	n = 42L
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))

	incid_des = make_fixed_design(
		"incidence",
		X,
		y_fun = function(w, X){
			stats::rbinom(nrow(X), 1, stats::plogis(-0.6 + 0.9 * w + 0.4 * X$x1 - 0.2 * X$x2))
		}
	)

	ordinal_des = make_fixed_design(
		"ordinal",
		X,
		y_fun = function(w, X){
			latent = -0.3 + 0.7 * w + 0.4 * X$x1 - 0.2 * X$x2 + stats::rlogis(nrow(X))
			as.integer(cut(latent, breaks = c(-Inf, -0.6, 0.2, 1.0, Inf), labels = FALSE))
		}
	)

	compare_bootstrap_fast_slow_asymp(
		InferenceIncidUnivGCompRiskDiff$new(incid_des, verbose = FALSE),
		SlowInferenceIncidUnivGCompRiskDiff$new(incid_des, verbose = FALSE),
		seed = 201
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceIncidMultiGCompRiskDiff$new(incid_des, verbose = FALSE),
		SlowInferenceIncidMultiGCompRiskDiff$new(incid_des, verbose = FALSE),
		seed = 202
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceIncidUnivGCompRiskRatio$new(incid_des, verbose = FALSE),
		SlowInferenceIncidUnivGCompRiskRatio$new(incid_des, verbose = FALSE),
		seed = 203
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceIncidMultiGCompRiskRatio$new(incid_des, verbose = FALSE),
		SlowInferenceIncidMultiGCompRiskRatio$new(incid_des, verbose = FALSE),
		seed = 204
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceOrdinalUniGCompMeanDiff$new(ordinal_des, verbose = FALSE),
		SlowInferenceOrdinalUniGCompMeanDiff$new(ordinal_des, verbose = FALSE),
		seed = 205
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceOrdinalMultiGCompMeanDiff$new(ordinal_des, verbose = FALSE),
		SlowInferenceOrdinalMultiGCompMeanDiff$new(ordinal_des, verbose = FALSE),
		seed = 206
	)
})

test_that("continuous lin, count negbin, and classical incidence estimators match generic bootstrap", {
	SlowInferenceContinMultLin = R6::R6Class(
		"SlowInferenceContinMultLin",
		inherit = InferenceContinMultLin,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountUnivNegBinRegr = R6::R6Class(
		"SlowInferenceCountUnivNegBinRegr",
		inherit = InferenceCountUnivNegBinRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountMultiNegBinRegr = R6::R6Class(
		"SlowInferenceCountMultiNegBinRegr",
		inherit = InferenceCountMultiNegBinRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountUnivHurdleNegBinRegr = R6::R6Class(
		"SlowInferenceCountUnivHurdleNegBinRegr",
		inherit = InferenceCountUnivHurdleNegBinRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountMultiHurdleNegBinRegr = R6::R6Class(
		"SlowInferenceCountMultiHurdleNegBinRegr",
		inherit = InferenceCountMultiHurdleNegBinRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceIncidUnivRiskDiff = R6::R6Class(
		"SlowInferenceIncidUnivRiskDiff",
		inherit = InferenceIncidUnivRiskDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceIncidUnivNewcombeRiskDiff = R6::R6Class(
		"SlowInferenceIncidUnivNewcombeRiskDiff",
		inherit = InferenceIncidUnivNewcombeRiskDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceIncidUnivMiettinenNurminenRiskDiff = R6::R6Class(
		"SlowInferenceIncidUnivMiettinenNurminenRiskDiff",
		inherit = InferenceIncidUnivMiettinenNurminenRiskDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)

	set.seed(20260416)
	n = 44L
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))

	continuous_des = make_fixed_design(
		"continuous",
		X,
		y_fun = function(w, X){
			0.5 * w + 0.4 * X$x1 - 0.3 * X$x2 + stats::rnorm(nrow(X), sd = 0.8)
		}
	)
	count_des = make_fixed_design(
		"count",
		X,
		y_fun = function(w, X){
			stats::rnbinom(nrow(X), mu = exp(0.2 + 0.45 * w + 0.25 * X$x1), size = 2.2)
		}
	)
	hurdle_des = make_fixed_design(
		"count",
		X,
		y_fun = function(w, X){
			is_zero = stats::rbinom(nrow(X), 1, stats::plogis(0.2 - 0.8 * w + 0.2 * X$x2))
			pos = 1 + stats::rnbinom(nrow(X), mu = exp(0.1 + 0.35 * w + 0.2 * X$x1), size = 1.8)
			ifelse(is_zero == 1, 0, pos)
		}
	)
	incid_des = make_fixed_design(
		"incidence",
		X,
		y_fun = function(w, X){
			stats::rbinom(nrow(X), 1, stats::plogis(-0.5 + 1.0 * w))
		}
	)

	compare_bootstrap_fast_slow_asymp(
		InferenceContinMultLin$new(continuous_des, verbose = FALSE),
		SlowInferenceContinMultLin$new(continuous_des, verbose = FALSE),
		seed = 207
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceCountUnivNegBinRegr$new(count_des, verbose = FALSE),
		SlowInferenceCountUnivNegBinRegr$new(count_des, verbose = FALSE),
		seed = 208
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceCountMultiNegBinRegr$new(count_des, verbose = FALSE),
		SlowInferenceCountMultiNegBinRegr$new(count_des, verbose = FALSE),
		seed = 209
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceCountUnivHurdleNegBinRegr$new(hurdle_des, verbose = FALSE),
		SlowInferenceCountUnivHurdleNegBinRegr$new(hurdle_des, verbose = FALSE),
		seed = 210
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceCountMultiHurdleNegBinRegr$new(hurdle_des, verbose = FALSE),
		SlowInferenceCountMultiHurdleNegBinRegr$new(hurdle_des, verbose = FALSE),
		seed = 211
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceIncidUnivRiskDiff$new(incid_des, verbose = FALSE),
		SlowInferenceIncidUnivRiskDiff$new(incid_des, verbose = FALSE),
		seed = 212
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceIncidUnivNewcombeRiskDiff$new(incid_des, verbose = FALSE),
		SlowInferenceIncidUnivNewcombeRiskDiff$new(incid_des, verbose = FALSE),
		seed = 213
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceIncidUnivMiettinenNurminenRiskDiff$new(incid_des, verbose = FALSE),
		SlowInferenceIncidUnivMiettinenNurminenRiskDiff$new(incid_des, verbose = FALSE),
		seed = 214
	)
})

test_that("continuous blocked-cluster bootstrap keeps multivariate workers finite", {
	set.seed(20260417)
	block = rep(c("a", "b", "c"), each = 5L)
	X = data.frame(
		block = factor(block),
		x1 = rnorm(length(block)),
		x2 = rnorm(length(block)),
		x3 = rnorm(length(block))
	)
	des = make_fixed_blocked_cluster_design(
		X,
		y_fun = function(w, X_design){
			0.6 * w + 0.4 * X_design$x1 - 0.25 * X_design$x2 + stats::rnorm(nrow(X_design), sd = 0.5)
		}
	)

	inf_ols = InferenceContinMultOLS$new(des, verbose = FALSE)
	inf_lin = InferenceContinMultLin$new(des, verbose = FALSE)
	inf_robust = InferenceContinMultiRobustRegr$new(des, method = "M", verbose = FALSE)

	set.seed(215)
	dbg_ols = inf_ols$approximate_bootstrap_distribution_beta_hat_T(B = 25L, show_progress = FALSE, debug = TRUE)
	set.seed(216)
	dbg_lin = inf_lin$approximate_bootstrap_distribution_beta_hat_T(B = 25L, show_progress = FALSE, debug = TRUE)
	set.seed(217)
	dbg_robust = inf_robust$approximate_bootstrap_distribution_beta_hat_T(B = 25L, show_progress = FALSE, debug = TRUE)

	expect_equal(sum(is.finite(dbg_ols$values)), 25L)
	expect_equal(sum(is.finite(dbg_lin$values)), 25L)
	expect_equal(sum(is.finite(dbg_robust$values)), 25L)
	expect_false(any(grepl("incompatible dimensions|number of rows of result", unlist(dbg_ols$errors))))
	expect_false(any(grepl("incompatible dimensions|number of rows of result", unlist(dbg_lin$errors))))
	expect_false(any(grepl("incompatible dimensions|number of rows of result", unlist(dbg_robust$errors))))

	set.seed(218)
	ci_ols = inf_ols$compute_bootstrap_confidence_interval(B = 51L, show_progress = FALSE)
	set.seed(219)
	ci_lin = inf_lin$compute_bootstrap_confidence_interval(B = 51L, show_progress = FALSE)
	set.seed(220)
	ci_robust = inf_robust$compute_bootstrap_confidence_interval(B = 51L, show_progress = FALSE)
	expect_true(all(is.finite(ci_ols)))
	expect_true(all(is.finite(ci_lin)))
	expect_true(all(is.finite(ci_robust)))
})

test_that("MLE and proportion families picked up through InferenceAsymp match generic bootstrap", {
	SlowInferenceOrdinalUniPropOddsRegr = R6::R6Class(
		"SlowInferenceOrdinalUniPropOddsRegr",
		inherit = InferenceOrdinalUniPropOddsRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferencePropUniBetaRegr = R6::R6Class(
		"SlowInferencePropUniBetaRegr",
		inherit = InferencePropUniBetaRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferencePropUniFractionalLogit = R6::R6Class(
		"SlowInferencePropUniFractionalLogit",
		inherit = InferencePropUniFractionalLogit,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferencePropMultiFractionalLogit = R6::R6Class(
		"SlowInferencePropMultiFractionalLogit",
		inherit = InferencePropMultiFractionalLogit,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferencePropUniZeroOneInflatedBetaRegr = R6::R6Class(
		"SlowInferencePropUniZeroOneInflatedBetaRegr",
		inherit = InferencePropUniZeroOneInflatedBetaRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferencePropMultiZeroOneInflatedBetaRegr = R6::R6Class(
		"SlowInferencePropMultiZeroOneInflatedBetaRegr",
		inherit = InferencePropMultiZeroOneInflatedBetaRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)

	set.seed(20260417)
	n = 40L
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))

	ordinal_des = make_fixed_design(
		"ordinal",
		X,
		y_fun = function(w, X){
			latent = 0.5 * w + 0.3 * X$x1 - 0.2 * X$x2 + stats::rlogis(nrow(X))
			as.integer(cut(latent, breaks = c(-Inf, -0.8, 0.0, 0.8, Inf), labels = FALSE))
		}
	)
	prop_des = make_fixed_design(
		"proportion",
		X,
		y_fun = function(w, X){
			stats::plogis(-0.3 + 0.8 * w + 0.2 * X$x1 + stats::rnorm(nrow(X), sd = 0.35))
		}
	)
	zoib_des = make_fixed_design(
		"proportion",
		X,
		y_fun = function(w, X){
			u = stats::runif(nrow(X))
			mu = stats::plogis(-0.4 + 0.9 * w + 0.2 * X$x1)
			out = stats::rbeta(nrow(X), shape1 = 6 * mu + 0.5, shape2 = 6 * (1 - mu) + 0.5)
			out[u < 0.15] = 0
			out[u > 0.85] = 1
			out
		}
	)

	compare_bootstrap_fast_slow_asymp(
		InferenceOrdinalUniPropOddsRegr$new(ordinal_des, verbose = FALSE),
		SlowInferenceOrdinalUniPropOddsRegr$new(ordinal_des, verbose = FALSE),
		seed = 215
	)
	compare_bootstrap_fast_slow_asymp(
		InferencePropUniBetaRegr$new(prop_des, verbose = FALSE),
		SlowInferencePropUniBetaRegr$new(prop_des, verbose = FALSE),
		seed = 216
	)
	compare_bootstrap_fast_slow_asymp(
		InferencePropUniFractionalLogit$new(prop_des, verbose = FALSE),
		SlowInferencePropUniFractionalLogit$new(prop_des, verbose = FALSE),
		seed = 217
	)
	compare_bootstrap_fast_slow_asymp(
		InferencePropMultiFractionalLogit$new(prop_des, verbose = FALSE),
		SlowInferencePropMultiFractionalLogit$new(prop_des, verbose = FALSE),
		seed = 218
	)
	compare_bootstrap_fast_slow_asymp(
		InferencePropUniZeroOneInflatedBetaRegr$new(zoib_des, verbose = FALSE),
		SlowInferencePropUniZeroOneInflatedBetaRegr$new(zoib_des, verbose = FALSE),
		seed = 219
	)
	compare_bootstrap_fast_slow_asymp(
		InferencePropMultiZeroOneInflatedBetaRegr$new(zoib_des, verbose = FALSE),
		SlowInferencePropMultiZeroOneInflatedBetaRegr$new(zoib_des, verbose = FALSE),
		seed = 220
	)
})

test_that("survival reusable workers match generic bootstrap", {
	SlowInferenceSurvivalUniCoxPHRegr = R6::R6Class(
		"SlowInferenceSurvivalUniCoxPHRegr",
		inherit = InferenceSurvivalUniCoxPHRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceSurvivalMultiCoxPHRegr = R6::R6Class(
		"SlowInferenceSurvivalMultiCoxPHRegr",
		inherit = InferenceSurvivalMultiCoxPHRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceSurvivalUniStratCoxPHRegr = R6::R6Class(
		"SlowInferenceSurvivalUniStratCoxPHRegr",
		inherit = InferenceSurvivalUniStratCoxPHRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceSurvivalMultiStratCoxPHRegr = R6::R6Class(
		"SlowInferenceSurvivalMultiStratCoxPHRegr",
		inherit = InferenceSurvivalMultiStratCoxPHRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceSurvivalKMDiff = R6::R6Class(
		"SlowInferenceSurvivalKMDiff",
		inherit = InferenceSurvivalKMDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceSurvivalLogRank = R6::R6Class(
		"SlowInferenceSurvivalLogRank",
		inherit = InferenceSurvivalLogRank,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceSurvivalRestrictedMeanDiff = R6::R6Class(
		"SlowInferenceSurvivalRestrictedMeanDiff",
		inherit = InferenceSurvivalRestrictedMeanDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceSurvivalGehanWilcox = R6::R6Class(
		"SlowInferenceSurvivalGehanWilcox",
		inherit = InferenceSurvivalGehanWilcox,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)

	set.seed(20260418)
	n = 46L
	X = data.frame(x1 = rnorm(n), x2 = sample(0:1, n, replace = TRUE), x3 = rnorm(n))
	surv_des = make_fixed_design(
		"survival",
		X,
		y_fun = function(w, X){
			base_rate = exp(-0.4 + 0.5 * w + 0.2 * X$x1 - 0.3 * X$x2)
			pmax(stats::rexp(nrow(X), rate = base_rate), 0.05)
		},
		dead_fun = function(w, X, y){
			stats::rbinom(length(y), 1, stats::plogis(1.0 - 0.12 * y + 0.2 * w))
		}
	)

	compare_bootstrap_fast_slow_asymp(
		InferenceSurvivalUniCoxPHRegr$new(surv_des, verbose = FALSE),
		SlowInferenceSurvivalUniCoxPHRegr$new(surv_des, verbose = FALSE),
		seed = 221
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceSurvivalMultiCoxPHRegr$new(surv_des, verbose = FALSE),
		SlowInferenceSurvivalMultiCoxPHRegr$new(surv_des, verbose = FALSE),
		seed = 222
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceSurvivalUniStratCoxPHRegr$new(surv_des, verbose = FALSE),
		SlowInferenceSurvivalUniStratCoxPHRegr$new(surv_des, verbose = FALSE),
		seed = 223
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceSurvivalMultiStratCoxPHRegr$new(surv_des, verbose = FALSE),
		SlowInferenceSurvivalMultiStratCoxPHRegr$new(surv_des, verbose = FALSE),
		seed = 224
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceSurvivalKMDiff$new(surv_des, verbose = FALSE),
		SlowInferenceSurvivalKMDiff$new(surv_des, verbose = FALSE),
		seed = 225
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceSurvivalLogRank$new(surv_des, verbose = FALSE),
		SlowInferenceSurvivalLogRank$new(surv_des, verbose = FALSE),
		seed = 226
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceSurvivalRestrictedMeanDiff$new(surv_des, verbose = FALSE),
		SlowInferenceSurvivalRestrictedMeanDiff$new(surv_des, verbose = FALSE),
		seed = 227
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceSurvivalGehanWilcox$new(surv_des, verbose = FALSE),
		SlowInferenceSurvivalGehanWilcox$new(surv_des, verbose = FALSE),
		seed = 228
	)
})

test_that("zero-augmented count reusable workers match generic bootstrap", {
	skip_if_not_installed("glmmTMB")

	SlowInferenceCountUnivZeroInflatedPoissonRegr = R6::R6Class(
		"SlowInferenceCountUnivZeroInflatedPoissonRegr",
		inherit = InferenceCountUnivZeroInflatedPoissonRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountMultiZeroInflatedPoissonRegr = R6::R6Class(
		"SlowInferenceCountMultiZeroInflatedPoissonRegr",
		inherit = InferenceCountMultiZeroInflatedPoissonRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountUnivZeroInflatedNegBinRegr = R6::R6Class(
		"SlowInferenceCountUnivZeroInflatedNegBinRegr",
		inherit = InferenceCountUnivZeroInflatedNegBinRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountMultiZeroInflatedNegBinRegr = R6::R6Class(
		"SlowInferenceCountMultiZeroInflatedNegBinRegr",
		inherit = InferenceCountMultiZeroInflatedNegBinRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountUnivHurdlePoissonRegr = R6::R6Class(
		"SlowInferenceCountUnivHurdlePoissonRegr",
		inherit = InferenceCountUnivHurdlePoissonRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountMultiHurdlePoissonRegr = R6::R6Class(
		"SlowInferenceCountMultiHurdlePoissonRegr",
		inherit = InferenceCountMultiHurdlePoissonRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)

	set.seed(20260419)
	n = 36L
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	des = make_fixed_design(
		"count",
		X,
		y_fun = function(w, X){
			is_zero = stats::rbinom(nrow(X), 1, stats::plogis(0.4 - 0.9 * w + 0.2 * X$x1))
			pos = stats::rpois(nrow(X), lambda = exp(0.3 + 0.4 * w + 0.2 * X$x2))
			ifelse(is_zero == 1, 0, pos)
		}
	)

	compare_bootstrap_fast_slow_asymp(
		InferenceCountUnivZeroInflatedPoissonRegr$new(des, verbose = FALSE),
		SlowInferenceCountUnivZeroInflatedPoissonRegr$new(des, verbose = FALSE),
		seed = 229,
		tolerance = 1e-8
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceCountMultiZeroInflatedPoissonRegr$new(des, verbose = FALSE),
		SlowInferenceCountMultiZeroInflatedPoissonRegr$new(des, verbose = FALSE),
		seed = 230,
		tolerance = 1e-8
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceCountUnivZeroInflatedNegBinRegr$new(des, verbose = FALSE),
		SlowInferenceCountUnivZeroInflatedNegBinRegr$new(des, verbose = FALSE),
		seed = 231,
		tolerance = 1e-8
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceCountMultiZeroInflatedNegBinRegr$new(des, verbose = FALSE),
		SlowInferenceCountMultiZeroInflatedNegBinRegr$new(des, verbose = FALSE),
		seed = 232,
		tolerance = 1e-8
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceCountUnivHurdlePoissonRegr$new(des, verbose = FALSE),
		SlowInferenceCountUnivHurdlePoissonRegr$new(des, verbose = FALSE),
		seed = 233,
		tolerance = 1e-8
	)
	compare_bootstrap_fast_slow_asymp(
		InferenceCountMultiHurdlePoissonRegr$new(des, verbose = FALSE),
		SlowInferenceCountMultiHurdlePoissonRegr$new(des, verbose = FALSE),
		seed = 234,
		tolerance = 1e-8
	)
})
