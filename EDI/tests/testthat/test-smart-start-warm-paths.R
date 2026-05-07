test_that("likelihood inversion paths populate and reuse null-fit warm-start caches", {
	set.seed(202)
	n <- 60
	x <- rnorm(n)
	w <- rep(c(1, 0), length.out = n)
	des <- FixedDesign$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = x))
	des$overwrite_all_subject_assignments(w)
	linpred <- -0.2 + 0.7 * w + 0.3 * x
	des$add_all_subject_responses(rbinom(n, 1, plogis(linpred)))

	inf <- InferenceIncidLogRegr$new(des, verbose = FALSE, smart_default = TRUE)
	inf$set_testing_type("lik_ratio")
	p1 <- inf$compute_asymp_two_sided_pval(0.05)
	p2 <- inf$compute_asymp_two_sided_pval(0.08)
	cache_lr <- inf$.__enclos_env__$private$likelihood_null_warm_cache[["likelihood_test:lik_ratio"]]
	expect_true(is.finite(p1) && is.finite(p2))
	expect_true(!is.null(cache_lr$start))
	expect_true(length(cache_lr$start) > 0L)

	inf$set_testing_type("score")
	ci <- inf$compute_asymp_confidence_interval(alpha = 0.2)
	cache_score <- inf$.__enclos_env__$private$likelihood_null_warm_cache[["likelihood_test:score"]]
	expect_true(all(is.finite(ci)))
	expect_true(!is.null(cache_score$start))
})

test_that("bootstrap reusable workers retain warm-start state across refits", {
	set.seed(303)
	n <- 50
	x <- rnorm(n)
	w <- rep(c(1, 0), length.out = n)
	des <- FixedDesign$new(n = n, response_type = "count", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = x))
	des$overwrite_all_subject_assignments(w)
	mu <- exp(0.25 + 0.4 * w + 0.2 * x)
	des$add_all_subject_responses(rpois(n, mu))

	inf <- InferenceCountPoisson$new(des, verbose = FALSE, smart_default = TRUE)
	worker_state <- inf$.__enclos_env__$private$create_bootstrap_worker_state()
	idx1 <- sample.int(n, n, replace = TRUE)
	idx2 <- sample.int(n, n, replace = TRUE)

	inf$.__enclos_env__$private$load_bootstrap_sample_into_worker(worker_state, idx1)
	val1 <- inf$.__enclos_env__$private$compute_bootstrap_worker_estimate(worker_state)
	start1 <- worker_state$worker_priv$fit_warm_start
	inf$.__enclos_env__$private$load_bootstrap_sample_into_worker(worker_state, idx2)
	val2 <- inf$.__enclos_env__$private$compute_bootstrap_worker_estimate(worker_state)
	start2 <- worker_state$worker_priv$fit_warm_start

	expect_true(is.finite(val1))
	expect_true(is.finite(val2))
	expect_true(!is.null(start1) && length(start1) > 0L)
	expect_true(!is.null(start2) && length(start2) > 0L)
})

test_that("randomization CI p-value search reuses object-level warm starts across nearby deltas", {
	set.seed(404)
	n <- 40
	x <- rnorm(n)
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "count", verbose = FALSE)
	for (i in seq_len(n)) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x = x[i]))
	}
	add_all_subject_responses_seq(des, rpois(n, lambda = exp(0.2 + 0.3 * des$.__enclos_env__$private$w)))

	inf <- InferenceCountPoisson$new(des, verbose = FALSE, smart_default = TRUE)
	priv <- inf$.__enclos_env__$private
	perms <- priv$generate_permutations(25L)
	cache_env <- new.env(parent = emptyenv())
	ctrl <- priv$normalize_randomization_ci_search_control(NULL, r = 25L, pval_epsilon = 0.05)

	p1 <- priv$compute_randomization_ci_pval_cached(inf, 25L, 0.05, "log", perms, ctrl, cache_env)
	start1 <- priv$fit_warm_start
	p2 <- priv$compute_randomization_ci_pval_cached(inf, 25L, 0.08, "log", perms, ctrl, cache_env)
	start2 <- priv$fit_warm_start

	expect_true(is.finite(p1) || is.na(p1))
	expect_true(is.finite(p2) || is.na(p2))
	expect_true(!is.null(start1) && length(start1) > 0L)
	expect_true(!is.null(start2) && length(start2) > 0L)
})
