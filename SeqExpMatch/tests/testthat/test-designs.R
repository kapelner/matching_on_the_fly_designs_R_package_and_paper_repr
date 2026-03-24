test_that("SeqDesignBernoulli works", {
	# Fixed sample
	des <- SeqDesignBernoulli$new(n = 10, verbose = FALSE)
	expect_true(des$is_fixed_sample_size())
	expect_equal(des$get_n(), 10)

	x_new <- data.frame(x1 = rnorm(1))
	w <- des$add_subject_to_experiment_and_assign(x_new)
	expect_true(w %in% c(0, 1))
	expect_equal(des$get_t(), 1)

	des$add_subject_response(1, 0.5)
	expect_equal(des$get_y()[1], 0.5)

	# Not fixed sample
	des <- SeqDesignBernoulli$new(n = NULL, verbose = FALSE)
	expect_false(des$is_fixed_sample_size())
	des$add_subject_to_experiment_and_assign(x_new)
	expect_equal(des$get_t(), 1)
})

test_that("SeqDesignEfron works", {
	des <- SeqDesignEfron$new(n = 10, prob_T = 0.5, verbose = FALSE)
	expect_true(des$is_fixed_sample_size())

	x_new <- data.frame(x1 = rnorm(1))
	w <- des$add_subject_to_experiment_and_assign(x_new)
	expect_true(w %in% c(0, 1))
})

test_that("SeqDesignAtkinson works", {
	des <- SeqDesignAtkinson$new(n = 10, verbose = FALSE)
	x_new <- data.frame(x1 = rnorm(1), x2 = rnorm(1))
	w <- des$add_subject_to_experiment_and_assign(x_new)
	expect_true(w %in% c(0, 1))
})

test_that("SeqDesignKK14 works", {
	des <- SeqDesignKK14$new(n = 10, verbose = FALSE)
	x_new <- data.frame(x1 = rnorm(1), x2 = rnorm(1))
	w <- des$add_subject_to_experiment_and_assign(x_new)
	expect_true(w %in% c(0, 1))
})

test_that("SeqDesignKK21 works", {
	des <- SeqDesignKK21$new(n = 10, verbose = FALSE)
	x_new <- data.frame(x1 = rnorm(1), x2 = rnorm(1))
	w <- des$add_subject_to_experiment_and_assign(x_new)
	expect_true(w %in% c(0, 1))
})

test_that("SeqDesignKK21stepwise works", {
	des <- SeqDesignKK21stepwise$new(n = 10, verbose = FALSE)
	x_new <- data.frame(x1 = rnorm(1), x2 = rnorm(1))
	w <- des$add_subject_to_experiment_and_assign(x_new)
	expect_true(w %in% c(0, 1))
})

test_that("SeqDesigniBCRD works", {
	des <- SeqDesigniBCRD$new(n = 10, verbose = FALSE)
	x_new <- data.frame(x1 = rnorm(1))
	w <- des$add_subject_to_experiment_and_assign(x_new)
	expect_true(w %in% c(0, 1))
})

test_that("SeqDesignUrn works", {
	des <- SeqDesignUrn$new(n = 12, alpha = 1, beta = 1, verbose = FALSE)
	for (i in 1:12) {
		w <- des$add_subject_to_experiment_and_assign(data.frame(x1 = i))
		expect_true(w %in% c(0, 1))
	}
	expect_equal(des$get_t(), 12)
	expect_true(all(des$get_w()[1:12] %in% c(0, 1)))

	W <- des$draw_ws_according_to_design(r = 3)
	expect_equal(dim(W), c(12, 3))
	expect_true(all(W %in% c(0, 1)))
})

test_that("SeqDesignSPBR works", {
	des <- SeqDesignSPBR$new(strata_cols = "stratum", block_size = 4, n = 8, verbose = FALSE)
	X <- data.frame(
		stratum = rep(c("A", "B"), each = 4),
		x1 = rnorm(8)
	)
	for (i in 1:8) {
		w <- des$add_subject_to_experiment_and_assign(X[i, , drop = FALSE])
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

test_that("SeqDesignRandomBlockSize works", {
	des <- SeqDesignRandomBlockSize$new(
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
		w <- des$add_subject_to_experiment_and_assign(X[i, , drop = FALSE])
		expect_true(w %in% c(0, 1))
	}

	ws <- des$get_w()[1:8]
	expect_equal(sum(ws[1:4]), 2)
	expect_equal(sum(ws[5:8]), 2)
})

test_that("Response types work", {
	types <- c("continuous", "incidence", "proportion", "count", "survival")
	for (rt in types) {
	des <- SeqDesignBernoulli$new(n = 10, response_type = rt, verbose = FALSE)
	expect_equal(des$get_response_type(), rt)

	x_new <- data.frame(x1 = rnorm(1))
	des$add_subject_to_experiment_and_assign(x_new)

	val <- switch(rt,
					continuous = 1.5,
					incidence = 1,
					proportion = 0.5,
					count = 5,
					survival = 10)

	if (rt == "survival") {
		des$add_subject_response(1, val, dead = 1)
		expect_equal(des$get_dead()[1], 1)
	} else {
		des$add_subject_response(1, val)
	}
	expect_equal(des$get_y()[1], val)
	}
})

test_that("m is only populated for blocking and KK designs", {
	bernoulli_des <- SeqDesignBernoulli$new(n = 4, verbose = FALSE)
	bernoulli_des$add_subject_to_experiment_and_assign(data.frame(x1 = 0))
	expect_null(bernoulli_des$get_m())

	blocking_des <- FixedDesignBlocking$new(strata_cols = "stratum", n = 4, verbose = FALSE)
	blocking_des$add_subject(data.frame(stratum = "A"))
	blocking_des$add_subject(data.frame(stratum = "A"))
	blocking_des$add_subject(data.frame(stratum = "B"))
	blocking_des$add_subject(data.frame(stratum = "B"))
	expect_equal(blocking_des$get_m(), c(1L, 1L, 2L, 2L))

	for (des in list(
		SeqDesignKK14$new(n = 6, verbose = FALSE),
		SeqDesignKK21$new(n = 6, verbose = FALSE),
		SeqDesignKK21stepwise$new(n = 6, verbose = FALSE)
	)) {
		expect_true(is.null(des$get_m()) || is.integer(des$get_m()) || is.numeric(des$get_m()))
		des$add_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
		expect_true(is.null(des$get_m()) || is.integer(des$get_m()) || is.numeric(des$get_m()))
	}
})
