test_that("DesignSeqOneByOnePocockSimon works", {
	n = 20
	X = data.frame(
		gender = sample(c("M", "F"), n, replace = TRUE),
		age_cat = sample(c("Young", "Old"), n, replace = TRUE)
	)
	
	des = DesignSeqOneByOnePocockSimon$new(response_type = "continuous", strata_cols = c("gender", "age_cat"), n = n, verbose = FALSE)
	
	for (i in 1:n) {
		des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
	}
	
	w = des$get_w()
	expect_length(w, n)
	expect_true(all(w %in% c(-1, 1)))
	
	# Check balance (Pocock-Simon should be reasonably balanced)
	# For n=20 with {-1,+1} encoding, sum should be close to 0
	expect_lte(abs(sum(w)), 8)
	
	# Test draw_ws
	W2 = des$draw_ws_according_to_design(r = 1)
	expect_equal(dim(W2), c(n, 1L))
	expect_true(all(W2 %in% c(-1, 1)))
})
