library(EDI)

# Canonical KM median per R's survival package: survival:::findq() finds the first
# time t where S(t) <= 0.5. When S(t) < 0.5 strictly, the median is t itself (a step
# function, no interpolation). Only when S(t) lands exactly on 0.5 does it average t
# with the *next* event time. If the curve never reaches 0.5, the result is NA.
canonical_median = function(y, dead, w = NULL) {
	fit = if (is.null(w)) {
		survival::survfit(survival::Surv(y, dead) ~ 1)
	} else {
		survival::survfit(survival::Surv(y, dead) ~ 1, weights = w)
	}
	as.numeric(stats::quantile(fit, 0.5)$quantile)
}

test_that("get_survival_stat_for_group matches canonical survfit median on a strict crossing (no interpolation)", {
	y = c(1, 2, 3, 4, 5)
	dead = c(1, 1, 1, 1, 1)
	expect_equal(EDI:::get_survival_stat_for_group(y, dead, "median"), canonical_median(y, dead))
	expect_equal(EDI:::get_survival_stat_for_group(y, dead, "median"), 3)
})

test_that("get_survival_stat_for_group matches canonical survfit median when S(t) lands exactly on 0.5", {
	y = as.numeric(1:10)
	dead = rep(1, 10)
	expect_equal(EDI:::get_survival_stat_for_group(y, dead, "median"), canonical_median(y, dead))
	expect_equal(EDI:::get_survival_stat_for_group(y, dead, "median"), 5.5)
})

test_that("get_survival_stat_for_group matches canonical survfit median with tied event times", {
	y = c(1, 1, 2, 3, 4, 5)
	dead = rep(1, 6)
	expect_equal(EDI:::get_survival_stat_for_group(y, dead, "median"), canonical_median(y, dead))
	expect_equal(EDI:::get_survival_stat_for_group(y, dead, "median"), 2.5)
})

test_that("get_survival_stat_for_group returns NA (not Inf) when the KM curve never reaches 0.5", {
	y = c(1, 2, 3, 4, 5)
	dead = c(1, 0, 0, 0, 0)
	expect_true(is.na(canonical_median(y, dead)))
	expect_true(is.na(EDI:::get_survival_stat_for_group(y, dead, "median")))
	expect_false(is.infinite(EDI:::get_survival_stat_for_group(y, dead, "median")))
})

test_that("get_survival_stat_diff matches the difference of canonical per-group medians", {
	y_t = c(1, 2, 3, 4)
	dead_t = c(1, 1, 1, 1)
	y_c = c(1, 2, 10, 20)
	dead_c = c(1, 1, 1, 1)

	y = c(y_t, y_c)
	dead = c(dead_t, dead_c)
	w = c(rep(1L, length(y_t)), rep(0L, length(y_c)))

	expected = canonical_median(y_t, dead_t) - canonical_median(y_c, dead_c)
	expect_equal(EDI:::get_survival_stat_diff(y, dead, w, "median"), expected)
	expect_equal(EDI:::get_survival_stat_diff(y, dead, w, "median"), -3.5)
})

test_that("get_survival_stat_diff returns NA when either group's median is non-estimable", {
	y_t = c(1, 2, 3, 4, 5)
	dead_t = c(1, 0, 0, 0, 0) # never crosses 0.5
	y_c = c(1, 2, 3, 4, 5)
	dead_c = c(1, 1, 1, 1, 1)

	y = c(y_t, y_c)
	dead = c(dead_t, dead_c)
	w = c(rep(1L, length(y_t)), rep(0L, length(y_c)))

	expect_true(is.na(EDI:::get_survival_stat_diff(y, dead, w, "median")))
})

test_that("compute_survival_stat_diff_rand_bootstrap_parallel_cpp's inline KM kernel matches canonical medians", {
	y_t = c(1, 2, 3, 4)
	dead_t = c(1, 1, 1, 1)
	y_c = c(1, 2, 10, 20)
	dead_c = c(1, 1, 1, 1)

	y0 = c(y_t, y_c)
	dead = as.integer(c(dead_t, dead_c))
	n = length(y0)
	i_mat = matrix(1:n, ncol = 1)
	w_mat = matrix(c(rep(1L, length(y_t)), rep(0L, length(y_c))), ncol = 1)

	expected = canonical_median(y_t, dead_t) - canonical_median(y_c, dead_c)
	res = EDI:::compute_survival_stat_diff_rand_bootstrap_parallel_cpp(
		as.numeric(y0), dead, i_mat, w_mat, 0, FALSE, NULL, 1L
	)
	expect_equal(res[1], expected)
})

test_that("weighted_km_median (R fallback) matches canonical survfit median, unweighted", {
	des = DesignSeqOneByOneBernoulli$new(n = 4, response_type = "survival")
	for (i in 1:4) des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	des$add_all_subject_responses(ys = c(1, 2, 3, 4), deads = c(1L, 1L, 1L, 1L))
	inf = InferenceSurvivalKMDiff$new(des, verbose = FALSE)
	priv = inf$.__enclos_env__$private

	y = as.numeric(1:10)
	dead = rep(1L, 10)
	expect_equal(priv$weighted_km_median(y, dead, rep(1, 10)), canonical_median(y, dead))
	expect_equal(priv$weighted_km_median(y, dead, rep(1, 10)), 5.5)
})

test_that("weighted_km_median (R fallback) matches canonical weighted survfit median", {
	des = DesignSeqOneByOneBernoulli$new(n = 4, response_type = "survival")
	for (i in 1:4) des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	des$add_all_subject_responses(ys = c(1, 2, 3, 4), deads = c(1L, 1L, 1L, 1L))
	inf = InferenceSurvivalKMDiff$new(des, verbose = FALSE)
	priv = inf$.__enclos_env__$private

	y = c(1, 2, 3, 4, 5, 6)
	dead = rep(1L, 6)
	wts = c(2, 1, 1, 3, 1, 2)
	expect_equal(priv$weighted_km_median(y, dead, wts), canonical_median(y, dead, wts))
	expect_equal(priv$weighted_km_median(y, dead, wts), 4)
})

test_that("weighted_km_median (R fallback) returns NA (not Inf) when the curve never reaches 0.5", {
	des = DesignSeqOneByOneBernoulli$new(n = 4, response_type = "survival")
	for (i in 1:4) des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	des$add_all_subject_responses(ys = c(1, 2, 3, 4), deads = c(1L, 1L, 1L, 1L))
	inf = InferenceSurvivalKMDiff$new(des, verbose = FALSE)
	priv = inf$.__enclos_env__$private

	y = c(1, 2, 3, 4, 5)
	dead = c(1, 0, 0, 0, 0)
	result = priv$weighted_km_median(y, dead, rep(1, 5))
	expect_true(is.na(result))
	expect_false(is.infinite(result))
})
