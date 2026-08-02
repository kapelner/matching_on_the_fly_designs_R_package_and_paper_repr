library(testthat)
library(EDI)

make_fork_seed_completed_fixed_design = function(seed = NULL, n = 10) {
	des = DesignFixedBernoulli$new(n = n, response_type = "continuous", seed = seed)
	set.seed(1)
	des$add_all_subjects_to_experiment(data.frame(x = rnorm(n)))
	des$assign_w_to_all_subjects()
	set.seed(2)
	des$add_all_subject_responses(rnorm(n))
	des
}

test_that("Inference seed: same seed + same num_cores gives same rand p-value (fork cluster)", {
	skip_on_os("windows")
	skip_if(
		isTRUE(EDI:::edi_env$mirai_has_been_used),
		"fork clusters cannot be started safely after mirai has been used in this R session"
	)
	des = make_fork_seed_completed_fixed_design(seed = NULL, n = 20)
	set_num_cores(2)
	on.exit(unset_num_cores(), add = TRUE)
	inf1 = InferenceAllSimpleMeanDiff$new(des); inf1$set_seed(42)
	inf2 = InferenceAllSimpleMeanDiff$new(des); inf2$set_seed(42)
	p1 = inf1$compute_rand_two_sided_pval(r = 99, show_progress = FALSE)
	p2 = inf2$compute_rand_two_sided_pval(r = 99, show_progress = FALSE)
	expect_identical(p1, p2)
})
