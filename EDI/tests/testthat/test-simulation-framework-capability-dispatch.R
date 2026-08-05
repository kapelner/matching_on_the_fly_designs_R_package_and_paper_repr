test_that("SimulationFramework discovers randomization inference by capability", {
	sf = SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(),
		inference_classes_and_params = list(),
		n = 4L,
		p = 1L,
		verbose = FALSE
	)
	priv = sf$.__enclos_env__$private
	priv$inf_types = c("rand_pval", "rand_ci")
	priv$current_response_type = "continuous"

	fake_inference = function(capabilities) {
		env = new.env(parent = emptyenv())
		env$supports = function(capability) {
			stats::setNames(capability %in% capabilities, capability)
		}
		env
	}

	expect_identical(
		priv$.valid_inference_types(fake_inference("randomization_test")),
		"rand_pval"
	)
	expect_identical(
		priv$.valid_inference_types(fake_inference(c("randomization_test", "randomization_ci"))),
		c("rand_pval", "rand_ci")
	)
	expect_identical(
		priv$.valid_inference_types(fake_inference(character())),
		character()
	)
})
