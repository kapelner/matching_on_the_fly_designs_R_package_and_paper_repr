.libPaths(c("Rlib", .libPaths()))
library(testthat)
library(EDI)
library(data.table)

test_that("InferenceIncidenceExactBinomial matches binom.test for FixedDesignBinaryMatch", {
	skip_if_not_installed("GreedyExperimentalDesign")

	x_dat <- data.table(
		x1 = c(-2.00, -2.01, -1.00, -1.01, 1.00, 1.01, 2.00, 2.01),
		x2 = c(0, 0, 1, 1, 0, 0, 1, 1)
	)
	des <- FixedDesignBinaryMatch$new(n = nrow(x_dat), response_type = "incidence", verbose = FALSE)
	for (i in seq_len(nrow(x_dat))) {
		des$add_subject(x_dat[i, ])
	}
	des$randomize()
	des$.__enclos_env__$private$ensure_bms_computed()
	m <- as.integer(des$.__enclos_env__$private$m)
	w <- des$get_w()
	y <- integer(length(w))
	for (pair_id in sort(unique(m[m > 0L]))) {
		idx <- which(m == pair_id)
		if (pair_id <= 3L) {
			y[idx[w[idx] == 1L]] <- 1L
			y[idx[w[idx] == 0L]] <- 0L
		} else {
			y[idx[w[idx] == 1L]] <- 0L
			y[idx[w[idx] == 0L]] <- 1L
		}
	}
	des$add_all_subject_responses(y)

	inf <- InferenceIncidenceExactBinomial$new(des, verbose = FALSE)
	ref <- stats::binom.test(3L, 4L, p = 0.5)

	expect_equal(inf$compute_treatment_estimate(), log((3 + 0.5) / (1 + 0.5)), tolerance = 1e-12)
	expect_equal(inf$compute_exact_two_sided_pval_for_treatment_effect(delta = 0), ref$p.value, tolerance = 1e-12)
	expect_equal(unname(inf$compute_exact_confidence_interval(alpha = 0.05)), stats::qlogis(ref$conf.int), tolerance = 1e-12)
})

test_that("InferenceIncidenceExactBinomial ignores KK reservoir data", {
	x_dat <- data.table(x1 = c(-3, 3, -3.01, 3.01, 0, 0.01))
	build_inf <- function(reservoir_vals) {
		des <- DesignSeqOneByOneKK14$new(
			n = nrow(x_dat),
			response_type = "incidence",
			t_0_pct = 0.34,
			lambda = 0.99,
			verbose = FALSE
		)
		for (i in seq_len(nrow(x_dat))) {
			des$add_subject_to_experiment_and_assign(x_dat[i, ])
		}
		m <- as.integer(des$.__enclos_env__$private$m)
		w <- des$.__enclos_env__$private$w
		y <- integer(length(w))
		matched_ids <- sort(unique(m[m > 0L]))
		for (j in seq_along(matched_ids)) {
			idx <- which(m == matched_ids[j])
			if (j == 1L) {
				y[idx[w[idx] == 1L]] <- 1L
				y[idx[w[idx] == 0L]] <- 0L
			} else {
				y[idx[w[idx] == 1L]] <- 0L
				y[idx[w[idx] == 0L]] <- 1L
			}
		}
		reservoir_idx <- which(m == 0L)
		y[reservoir_idx] <- reservoir_vals
		des$add_all_subject_responses(y)
		InferenceIncidenceExactBinomial$new(des, verbose = FALSE)
	}

	inf_a <- build_inf(c(1L, 0L))
	inf_b <- build_inf(c(0L, 1L))

	expect_equal(inf_a$compute_treatment_estimate(), inf_b$compute_treatment_estimate(), tolerance = 1e-12)
	expect_equal(
		inf_a$compute_exact_two_sided_pval_for_treatment_effect(delta = 0),
		inf_b$compute_exact_two_sided_pval_for_treatment_effect(delta = 0),
		tolerance = 1e-12
	)
})

test_that("InferenceIncidenceExactBinomial rejects unsupported designs", {
	des <- DesignSeqOneByOneBernoulli$new(n = 10, response_type = "incidence", verbose = FALSE)
	for (i in seq_len(10)) {
		des$add_subject_to_experiment_and_assign(data.table(x1 = i))
	}
	des$add_all_subject_responses(rep(c(0L, 1L), length.out = 10))

	expect_error(
		InferenceIncidenceExactBinomial$new(des, verbose = FALSE),
		"Exact binomial incidence inference requires FixedDesignBinaryMatch or KK matching designs"
	)
})
