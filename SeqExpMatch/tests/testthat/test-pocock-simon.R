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
	expect_true(all(w %in% c(0, 1)))
	
	# Check balance (Pocock-Simon should be reasonably balanced)
	# For n=20, it should be close to 10/10
	expect_lte(abs(sum(w) - n/2), 4) 
	
	# Test redraw
	des$assign_w_to_all_subjects()
	w2 = des$get_w()
	expect_length(w2, n)
})
