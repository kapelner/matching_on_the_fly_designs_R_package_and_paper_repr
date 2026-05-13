run_likelihood_method_smoke_suite <- function(){
	old_seed = if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
		get(".Random.seed", envir = .GlobalEnv)
	} else {
		NULL
	}
	on.exit({
		if (!is.null(old_seed)) {
			assign(".Random.seed", old_seed, envir = .GlobalEnv)
		} else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
			rm(".Random.seed", envir = .GlobalEnv)
		}
	}, add = TRUE)
	set.seed(20240508)

	call_all_methods <- function(inf, label){
		is_unsupported_method_error <- function(err){
			grepl("does not expose a likelihood-test specification", conditionMessage(err), fixed = TRUE)
		}

		pval_methods <- c(
			"compute_wald_two_sided_pval",
			"compute_score_two_sided_pval",
			"compute_gradient_two_sided_pval",
			"compute_lik_ratio_two_sided_pval"
		)
		ci_methods <- c(
			"compute_wald_confidence_interval",
			"compute_score_confidence_interval",
			"compute_gradient_confidence_interval",
			"compute_lik_ratio_confidence_interval"
		)

		for (method_name in pval_methods){
			method_fn = inf[[method_name]]
			stopifnot(is.function(method_fn))
			cat(sprintf("  [%s] calling %s ...\n", label, method_name))
			val = tryCatch(
				method_fn(delta = 0),
				error = function(e) {
					if (is_unsupported_method_error(e)) {
						cat(sprintf("  [%s] skipping %s: %s\n", label, method_name, conditionMessage(e)))
						return(NULL)
					}
					stop(label, ": ", method_name, " failed: ", e$message)
				}
			)
			if (is.null(val)) next
			cat(sprintf("  [%s] %s = %s\n", label, method_name, paste(val, collapse=", ")))
			stopifnot(is.numeric(val), length(val) == 1L, is.finite(val), val >= 0, val <= 1)
		}

		for (method_name in ci_methods){
			method_fn = inf[[method_name]]
			stopifnot(is.function(method_fn))
			cat(sprintf("  [%s] calling %s ...\n", label, method_name))
			val = tryCatch(
				method_fn(alpha = 0.2),
				error = function(e) {
					if (is_unsupported_method_error(e)) {
						cat(sprintf("  [%s] skipping %s: %s\n", label, method_name, conditionMessage(e)))
						return(NULL)
					}
					stop(label, ": ", method_name, " failed: ", e$message)
				}
			)
			if (is.null(val)) next
			cat(sprintf("  [%s] %s = %s\n", label, method_name, paste(val, collapse=", ")))
			stopifnot(is.numeric(val), length(val) == 2L, all(is.finite(val)), val[1] <= val[2])
		}

		invisible(TRUE)
	}

	make_fixed_incidence_design <- function(n = 40L){
		x1 = rnorm(n)
		des = DesignFixedBernoulli$new(n = n, response_type = "incidence", verbose = FALSE)
		des$add_all_subjects_to_experiment(data.frame(x1 = x1))
		des$assign_w_to_all_subjects()
		w = des$get_w()
		y = rbinom(n, 1, plogis(-0.3 + 0.6 * w + 0.4 * x1))
		des$add_all_subject_responses(y)
		des
	}

	make_fixed_count_design <- function(n = 40L){
		x1 = rnorm(n)
		des = DesignFixedBernoulli$new(n = n, response_type = "count", verbose = FALSE)
		des$add_all_subjects_to_experiment(data.frame(x1 = x1))
		des$assign_w_to_all_subjects()
		w = des$get_w()
		y = rpois(n, exp(0.2 + 0.3 * w + 0.25 * x1))
		des$add_all_subject_responses(y)
		des
	}

	make_fixed_survival_design <- function(n = 40L){
		x1 = rnorm(n)
		des = DesignFixedBernoulli$new(n = n, response_type = "survival", verbose = FALSE)
		des$add_all_subjects_to_experiment(data.frame(x1 = x1))
		des$assign_w_to_all_subjects()
		w = des$get_w()
		y = rexp(n, rate = exp(-0.2 + 0.15 * w + 0.2 * x1))
		des$add_all_subject_responses(y, rep(1L, n))
		des
	}

	make_kk_incidence_design <- function(n = 16L){
		des = DesignSeqOneByOneKK14$new(n = n, response_type = "incidence", verbose = FALSE)
		x1 = rnorm(n)
		x2 = rnorm(n)
		for (i in seq_len(n)){
			w_i = des$add_one_subject_to_experiment_and_assign(data.frame(x1 = x1[i], x2 = x2[i]))
			y_i = rbinom(1, 1, plogis(-0.15 + 0.5 * w_i + 0.25 * x1[i] - 0.1 * x2[i]))
			des$add_one_subject_response(i, y_i, 1L)
		}
		des
	}

	make_kk_survival_design <- function(n = 16L){
		des = DesignSeqOneByOneKK14$new(n = n, response_type = "survival", verbose = FALSE)
		x1 = rnorm(n)
		x2 = rnorm(n)
		for (i in seq_len(n)){
			w_i = des$add_one_subject_to_experiment_and_assign(data.frame(x1 = x1[i], x2 = x2[i]))
			y_i = rexp(1, rate = exp(-0.1 + 0.2 * w_i + 0.15 * x1[i] - 0.05 * x2[i]))
			des$add_one_subject_response(i, y_i, 1L)
		}
		des
	}

	make_kk_ordinal_design <- function(n = 32L){
		des = DesignSeqOneByOneKK14$new(n = n, response_type = "ordinal", verbose = FALSE)
		x1 = rnorm(n)
		x2 = rnorm(n)
		for (i in seq_len(n)){
			w_i = des$add_one_subject_to_experiment_and_assign(data.frame(x1 = x1[i], x2 = x2[i]))
			y_i = sample.int(4L, 1L, prob = c(
				plogis(-0.8 + 0.4 * w_i + 0.2 * x1[i] - 0.1 * x2[i]),
				0.2,
				0.2,
				0.4
			))
			des$add_one_subject_response(i, y_i, 1L)
		}
		des
	}

	results = list()
	results$count_poisson = call_all_methods(
		InferenceCountPoisson$new(make_fixed_count_design(), model_formula = ~ x1, verbose = FALSE),
		"InferenceCountPoisson"
	)
	results$count_robust_poisson = call_all_methods(
		InferenceCountRobustPoisson$new(make_fixed_count_design(), model_formula = ~ x1, verbose = FALSE),
		"InferenceCountRobustPoisson"
	)
	results$count_quasi_poisson = call_all_methods(
		InferenceCountQuasiPoisson$new(make_fixed_count_design(), model_formula = ~ x1, verbose = FALSE),
		"InferenceCountQuasiPoisson"
	)
	results$count_zip = call_all_methods(
		InferenceCountZeroInflatedPoisson$new(make_fixed_count_design(), model_formula = ~ x1, use_rcpp = TRUE, verbose = FALSE, optimization_alg = "lbfgs"),
		"InferenceCountZeroInflatedPoisson"
	)
	results$count_hurdle_poisson = call_all_methods(
		InferenceCountHurdlePoisson$new(make_fixed_count_design(), model_formula = ~ x1, use_rcpp = TRUE, verbose = FALSE, optimization_alg = "lbfgs"),
		"InferenceCountHurdlePoisson"
	)
	results$count_hurdle_negbin = call_all_methods(
		InferenceCountHurdleNegBin$new(make_fixed_count_design(), model_formula = ~ x1, verbose = FALSE),
		"InferenceCountHurdleNegBin"
	)
	results$incidence_modified_poisson = call_all_methods(
		InferenceIncidModifiedPoisson$new(make_fixed_incidence_design(), model_formula = ~ x1, verbose = FALSE),
		"InferenceIncidModifiedPoisson"
	)
	results$incidence_kk_modified_poisson = call_all_methods(
		InferenceIncidKKModifiedPoisson$new(make_kk_incidence_design(), model_formula = ~ x1, verbose = FALSE),
		"InferenceIncidKKModifiedPoisson"
	)
	results$survival_cox = call_all_methods(
		InferenceSurvivalCoxPHRegr$new(make_fixed_survival_design(), model_formula = ~ x1, use_rcpp = TRUE, verbose = FALSE),
		"InferenceSurvivalCoxPHRegr"
	)
	results$survival_strat_cox = call_all_methods(
		InferenceSurvivalStratCoxPHRegr$new(make_fixed_survival_design(), model_formula = ~ x1, use_rcpp = TRUE, verbose = FALSE),
		"InferenceSurvivalStratCoxPHRegr"
	)
	results$kk_survival_strat_cox = call_all_methods(
		InferenceSurvivalKKStratCoxOneLik$new(make_kk_survival_design(), model_formula = ~ x1, verbose = FALSE),
		"InferenceSurvivalKKStratCoxOneLik"
	)
	results$kk_survival_lwa_cox = call_all_methods(
		InferenceSurvivalKKLWACoxOneLik$new(make_kk_survival_design(), model_formula = ~ x1, verbose = FALSE),
		"InferenceSurvivalKKLWACoxOneLik"
	)
	results$kk_survival_clayton = call_all_methods(
		InferenceSurvivalKKClaytonCopulaOneLik$new(make_kk_survival_design(n = 64L), model_formula = ~ x1, verbose = FALSE),
		"InferenceSurvivalKKClaytonCopulaOneLik"
	)
	results$ordinal_kk_glmm = call_all_methods(
		InferenceOrdinalKKGLMM$new(make_kk_ordinal_design(), model_formula = ~ x1, use_rcpp = TRUE, verbose = FALSE),
		"InferenceOrdinalKKGLMM"
	)
	invisible(results)
}
