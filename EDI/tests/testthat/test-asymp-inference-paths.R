library(testthat)
library(EDI)

# Helper to run all 4 asymp paths for a given inference object
test_all_asymp_paths <- function(inf, label = "") {
	supported <- inf$get_supported_testing_types()
	
	for (tt in c("wald", "score", "gradient", "lik_ratio")) {
		if (tt %in% supported) {
			test_that(sprintf("%s: %s paths work", label, tt), {
				inf$set_testing_type(tt)
				
				# Test p-value
				pval <- inf$compute_asymp_two_sided_pval(delta = 0)
				expect_true(is.numeric(pval))
				expect_true(pval >= 0 && pval <= 1)
				
				# Test CI
				ci <- inf$compute_asymp_confidence_interval(alpha = 0.05)
				expect_true(is.numeric(ci))
				expect_length(ci, 2)
				expect_true(all(is.finite(ci)))
				expect_true(ci[1] <= ci[2])
			})
		}
	}
}

test_that("Logistic Asymp paths", {
	set.seed(1)
	n <- 100
	X <- data.frame(x1 = rnorm(n))
	y <- rbinom(n, 1, plogis(0.5 + X$x1))
	des <- DesignFixedBernoulli$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	des$assign_w_to_all_subjects()
	des$add_all_subject_responses(y)
	
	inf <- InferenceIncidLogRegr$new(des)
	test_all_asymp_paths(inf, "Logistic")
})

test_that("Poisson Asymp paths", {
	set.seed(2)
	n <- 100
	X <- data.frame(x1 = rnorm(n))
	y <- rpois(n, exp(0.5 + X$x1))
	des <- DesignFixedBernoulli$new(n = n, response_type = "count", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	des$add_all_subject_responses(y)
	
	inf <- InferenceCountPoisson$new(des)
	test_all_asymp_paths(inf, "Poisson")
})

test_that("NegBin Asymp paths", {
	set.seed(3)
	n <- 100
	X <- data.frame(x1 = rnorm(n))
	y <- rnbinom(n, size = 2, mu = exp(0.5 + X$x1))
	des <- DesignFixedBernoulli$new(n = n, response_type = "count", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	des$add_all_subject_responses(y)
	
	inf <- InferenceCountNegBin$new(des)
	test_all_asymp_paths(inf, "NegBin")
})

test_that("Cox Asymp paths", {
	skip_if_not_installed("survival")
	set.seed(4)
	n <- 100
	X <- data.frame(x1 = rnorm(n))
	y <- rexp(n, 0.1 * exp(X$x1))
	dead <- rbinom(n, 1, 0.8)
	des <- DesignFixedBernoulli$new(n = n, response_type = "survival", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	des$add_all_subject_responses(y, dead)
	
	inf <- InferenceSurvivalCoxPHRegr$new(des)
	test_all_asymp_paths(inf, "Cox")
})

test_that("Ordinal Asymp paths", {
	set.seed(5)
	n <- 150
	X <- data.frame(x1 = rnorm(n))
	y_cont <- 0.5 * X$x1 + rnorm(n)
	y <- as.numeric(cut(y_cont, breaks = c(-Inf, -0.5, 0.5, Inf), labels = FALSE))
	des <- DesignFixedBernoulli$new(n = n, response_type = "ordinal", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	des$add_all_subject_responses(y)
	
	inf <- InferenceOrdinalPropOddsRegr$new(des)
	test_all_asymp_paths(inf, "Ordinal")
})

test_that("Beta Asymp paths", {
	set.seed(6)
	n <- 100
	X <- data.frame(x1 = rnorm(n))
	mu <- plogis(0.5 + X$x1)
	phi <- 10
	y <- rbeta(n, mu * phi, (1 - mu) * phi)
	y <- pmax(pmin(y, 1 - 1e-6), 1e-6)
	des <- DesignFixedBernoulli$new(n = n, response_type = "proportion", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	des$add_all_subject_responses(y)
	
	inf <- InferencePropBetaRegr$new(des)
	test_all_asymp_paths(inf, "Beta")
})

test_that("Ordinal (additional) Asymp paths", {
	set.seed(7)
	n <- 150
	X <- data.frame(x1 = rnorm(n))
	y_cont <- 0.5 * X$x1 + rnorm(n)
	y <- as.numeric(cut(y_cont, breaks = c(-Inf, -0.5, 0.5, Inf), labels = FALSE))
	des <- DesignFixedBernoulli$new(n = n, response_type = "ordinal", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	des$add_all_subject_responses(y)
	
	for (cls in list(InferenceOrdinalCloglogRegr, InferenceOrdinalOrderedProbitRegr, InferenceOrdinalCauchitRegr)) {
		inf <- cls$new(des)
		test_all_asymp_paths(inf, cls$classname)
	}
})

test_that("KK Design Asymp paths", {
	set.seed(8)
	n <- 60
	X <- data.frame(x1 = rnorm(n))
	y_cont <- 0.5 * X$x1 + rnorm(n)
	
	des <- DesignSeqOneByOneKK14$new(n = n, response_type = "continuous", verbose = FALSE)
	for (i in seq_len(n)) {
		w_i <- des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
		y_i <- 0.5 * w_i + y_cont[i]
		des$add_one_subject_response(i, y_i, 1)
	}
	
	inf <- InferenceContinKKOLSOneLik$new(des)
	test_all_asymp_paths(inf, "KK OLS OneLik")
})

test_that("KK Design continuous GLMM Asymp paths", {
	set.seed(10)
	n <- 60
	X <- data.frame(x1 = rnorm(n))
	
	des <- DesignSeqOneByOneKK14$new(n = n, response_type = "continuous", verbose = FALSE)
	for (i in seq_len(n)) {
		w_i <- des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
		y_i <- 0.5 * w_i + 0.25 * X$x1[i] + rnorm(1)
		des$add_one_subject_response(i, y_i, 1)
	}
	
	inf <- InferenceContinKKGLMM$new(des)
	test_all_asymp_paths(inf, "KK Continuous GLMM")
})

test_that("KK Design count GLMM Asymp paths", {
	set.seed(11)
	n <- 80
	X <- data.frame(x1 = rnorm(n))
	
	des <- DesignSeqOneByOneKK14$new(n = n, response_type = "count", verbose = FALSE)
	for (i in seq_len(n)) {
		w_i <- des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
		mu_i <- exp(0.25 + 0.35 * w_i + 0.2 * X$x1[i])
		y_i <- rpois(1, mu_i)
		des$add_one_subject_response(i, y_i, 1)
	}
	
	inf <- InferenceCountKKGLMM$new(des)
	test_all_asymp_paths(inf, "KK Count GLMM")
})

test_that("KK Design Incidence Asymp paths", {
	set.seed(9)
	n <- 100
	X <- data.frame(x1 = rnorm(n))
	
	des <- DesignSeqOneByOneKK14$new(n = n, response_type = "incidence", verbose = FALSE)
	for (i in seq_len(n)) {
		w_i <- des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
		p_i <- plogis(0.5 * w_i + 0.5 * X$x1[i])
		y_i <- rbinom(1, 1, p_i)
		des$add_one_subject_response(i, y_i, 1)
	}
	
	inf <- InferenceIncidKKClogitOneLik$new(des)
	test_all_asymp_paths(inf, "KK Clogit OneLik")
})
