context("Pocock-Simon Minimization")

test_that("SeqDesignPocockSimon initialization and basic allocation", {
	n <- 20
	des <- SeqDesignPocockSimon$new(n = n, p_best = 0.9, verbose = FALSE)
	expect_true(des$is_fixed_sample_size())
	expect_equal(des$get_n(), n)

	# Add first subject
	x1 <- data.frame(c1 = "A", c2 = "X")
	w1 <- des$add_subject_to_experiment_and_assign(x1)
	expect_true(w1 %in% c(0, 1))
	expect_equal(des$get_t(), 1)
	
	# Add second subject with same covariates
	x2 <- data.frame(c1 = "A", c2 = "X")
	w2 <- des$add_subject_to_experiment_and_assign(x2)
	expect_true(w2 %in% c(0, 1))
})

test_that("SeqDesignPocockSimon balances covariates", {
	des <- SeqDesignPocockSimon$new(n = 20, p_best = 1.0, verbose = FALSE)
	
	set.seed(42)
	for (i in 1:10) {
		des$add_subject_to_experiment_and_assign(data.frame(cv = "L1"))
	}
	
	ws <- des$get_w()
	expect_equal(sum(ws[1:10] == 1), 5)
	expect_equal(sum(ws[1:10] == 0), 5)
	
	for (i in 1:10) {
		des$add_subject_to_experiment_and_assign(data.frame(cv = "L2"))
	}
	ws2 <- des$get_w()
	expect_equal(sum(ws2[11:20] == 1), 5)
	expect_equal(sum(ws2[11:20] == 0), 5)
})

test_that("SeqDesignPocockSimon redraw works", {
	n <- 10
	des <- SeqDesignPocockSimon$new(n = n, p_best = 0.8, verbose = FALSE)
	for (i in 1:n) {
		des$add_subject_to_experiment_and_assign(data.frame(v1 = sample(c("A", "B"), 1), v2 = sample(c("X", "Y"), 1)))
	}
	
	des$randomize()
	w_new <- des$get_w()
	
	expect_equal(length(w_new), n)
})

test_that("SeqDesignPocockSimon handles missing values", {
	des <- SeqDesignPocockSimon$new(n = 10, include_is_missing_as_a_new_feature = TRUE, verbose = FALSE)
	
	des$add_subject_to_experiment_and_assign(data.frame(v = "A"))
	des$add_subject_to_experiment_and_assign(data.frame(v = NA_character_))
	
	expect_equal(des$get_t(), 2)
})
