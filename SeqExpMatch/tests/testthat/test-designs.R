test_that("DesignSeqOneByOneBernoulli works", {
	# Fixed sample
	des <- DesignSeqOneByOneBernoulli$new(response_type = "continuous", n = 10, verbose = FALSE)
	expect_true(des$is_fixed_sample_size())
	expect_equal(des$get_n(), 10)

	x_new <- data.frame(x1 = rnorm(1))
	w <- des$add_one_subject_to_experiment_and_assign(x_new)
	expect_true(w %in% c(0, 1))
	expect_equal(des$get_t(), 1)

	des$add_one_subject_response(1, 0.5)
	expect_equal(des$get_y()[1], 0.5)

	# Not fixed sample
	des <- DesignSeqOneByOneBernoulli$new(response_type = "continuous", n = NULL, verbose = FALSE)
	expect_false(des$is_fixed_sample_size())
	des$add_one_subject_to_experiment_and_assign(x_new)
	expect_equal(des$get_t(), 1)
})

test_that("DesignSeqOneByOneEfron works", {
	des <- DesignSeqOneByOneEfron$new(response_type = "continuous", n = 10, prob_T = 0.5, verbose = FALSE)
	expect_true(des$is_fixed_sample_size())

	x_new <- data.frame(x1 = rnorm(1))
	w <- des$add_one_subject_to_experiment_and_assign(x_new)
	expect_true(w %in% c(0, 1))
})

test_that("DesignSeqOneByOneAtkinson works", {
	des <- DesignSeqOneByOneAtkinson$new(response_type = "continuous", n = 10, verbose = FALSE)
	x_new <- data.frame(x1 = rnorm(1), x2 = rnorm(1))
	w <- des$add_one_subject_to_experiment_and_assign(x_new)
	expect_true(w %in% c(0, 1))
})

test_that("DesignSeqOneByOneKK14 works", {
	des <- DesignSeqOneByOneKK14$new(response_type = "continuous", n = 10, verbose = FALSE)
	x_new <- data.frame(x1 = rnorm(1), x2 = rnorm(1))
	w <- des$add_one_subject_to_experiment_and_assign(x_new)
	expect_true(w %in% c(0, 1))
})

test_that("DesignSeqOneByOneKK21 works", {
	des <- DesignSeqOneByOneKK21$new(response_type = "continuous", n = 10, verbose = FALSE)
	x_new <- data.frame(x1 = rnorm(1), x2 = rnorm(1))
	w <- des$add_one_subject_to_experiment_and_assign(x_new)
	expect_true(w %in% c(0, 1))
})

test_that("DesignSeqOneByOneKK21stepwise works", {
	des <- DesignSeqOneByOneKK21stepwise$new(response_type = "continuous", n = 10, verbose = FALSE)
	x_new <- data.frame(x1 = rnorm(1), x2 = rnorm(1))
	w <- des$add_one_subject_to_experiment_and_assign(x_new)
	expect_true(w %in% c(0, 1))
})

test_that("DesignSeqOneByOneiBCRD works", {
	des <- DesignSeqOneByOneiBCRD$new(response_type = "continuous", n = 10, verbose = FALSE)
	x_new <- data.frame(x1 = rnorm(1))
	w <- des$add_one_subject_to_experiment_and_assign(x_new)
	expect_true(w %in% c(0, 1))
	expect_error(des$get_block_ids(), "undefined")
	for (i in 2:10) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
	}
	expect_identical(des$get_block_ids(), rep(1L, 10))
})

test_that("DesignSeqOneByOneUrn works", {
	des <- DesignSeqOneByOneUrn$new(response_type = "continuous", n = 12, alpha = 1, beta = 1, verbose = FALSE)
	for (i in 1:12) {
		w <- des$add_one_subject_to_experiment_and_assign(data.frame(x1 = i))
		expect_true(w %in% c(0, 1))
	}
	expect_equal(des$get_t(), 12)
	expect_true(all(des$get_w()[1:12] %in% c(0, 1)))

	W <- des$draw_ws_according_to_design(r = 3)
	expect_equal(dim(W), c(12, 3))
	expect_true(all(W %in% c(0, 1)))
})

test_that("DesignSeqOneByOneSPBR works", {
	des <- DesignSeqOneByOneSPBR$new(response_type = "continuous", strata_cols = "stratum", block_size = 4, n = 8, verbose = FALSE)
	X <- data.frame(
		stratum = rep(c("A", "B"), each = 4),
		x1 = rnorm(8)
	)
	for (i in 1:8) {
		w <- des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
		expect_true(w %in% c(0, 1))
	}

	ws <- des$get_w()[1:8]
	expect_equal(sum(ws[1:4]), 2)
	expect_equal(sum(ws[5:8]), 2)

	W <- des$draw_ws_according_to_design(r = 3)
	expect_equal(dim(W), c(8, 3))
	expect_true(all(colSums(W[1:4, , drop = FALSE]) == 2))
	expect_true(all(colSums(W[5:8, , drop = FALSE]) == 2))
})

test_that("DesignSeqOneByOneRandomBlockSize works", {
	des <- DesignSeqOneByOneRandomBlockSize$new(
		response_type = "continuous",
		strata_cols = "stratum",
		block_sizes = c(2, 4),
		n = 8,
		verbose = FALSE
	)
	X <- data.frame(
		stratum = rep(c("A", "B"), each = 4),
		x1 = rnorm(8)
	)
	for (i in 1:8) {
		w <- des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
		expect_true(w %in% c(0, 1))
	}

	ws <- des$get_w()[1:8]
	expect_equal(sum(ws[1:4]), 2)
	expect_equal(sum(ws[5:8]), 2)
})

test_that("Response types work", {
	types <- c("continuous", "incidence", "proportion", "count", "survival")
	for (rt in types) {
	des <- DesignSeqOneByOneBernoulli$new(n = 10, response_type = rt, verbose = FALSE)
	expect_equal(des$get_response_type(), rt)

	x_new <- data.frame(x1 = rnorm(1))
	des$add_one_subject_to_experiment_and_assign(x_new)

	val <- switch(rt,
					continuous = 1.5,
					incidence = 1,
					proportion = 0.5,
					count = 5,
					survival = 10)

	if (rt == "survival") {
		des$add_one_subject_response(1, val, dead = 1)
		expect_equal(des$get_dead()[1], 1)
	} else {
		des$add_one_subject_response(1, val)
	}
	expect_equal(des$get_y()[1], val)
	}
})

test_that("m is only populated for blocking and KK designs", {
	bernoulli_des <- DesignSeqOneByOneBernoulli$new(response_type = "continuous", n = 4, verbose = FALSE)
	bernoulli_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = 0))
	expect_null(bernoulli_des$get_m())

	blocking_des <- FixedDesignBlocking$new(response_type = "continuous", strata_cols = "stratum", n = 4, verbose = FALSE)
	blocking_des$add_all_subjects_to_experiment(data.frame(stratum = c("A", "A", "B", "B")))
	expect_equal(blocking_des$get_m(), c(1L, 1L, 2L, 2L))

	for (des in list(
		DesignSeqOneByOneKK14$new(response_type = "continuous", n = 6, verbose = FALSE),
		DesignSeqOneByOneKK21$new(response_type = "continuous", n = 6, verbose = FALSE),
		DesignSeqOneByOneKK21stepwise$new(response_type = "continuous", n = 6, verbose = FALSE)
	)) {
		expect_true(is.null(des$get_m()) || is.integer(des$get_m()) || is.numeric(des$get_m()))
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
		expect_true(is.null(des$get_m()) || is.integer(des$get_m()) || is.numeric(des$get_m()))
	}
})
