test_that("migration golden design builders are registered for all response types", {
	builders = inference_migration_golden_design_builders()
	expect_named(builders, c("continuous", "incidence", "count", "proportion", "ordinal", "survival"))
	expect_true(all(vapply(builders, is.function, logical(1L))))
})

test_that("migration output comparison helper compares supported methods", {
	ToyLegacy = R6::R6Class(
		"MigrationHarnessToyLegacy",
		lock_objects = FALSE,
		public = list(
			initialize = function(des_obj, verbose = FALSE) {
				private$des_obj = des_obj
			},
			compute_estimate = function(estimate_only = FALSE) 1.25,
			compute_wald_confidence_interval = function(alpha = 0.2) c(0.5, 2.0),
			compute_wald_two_sided_pval = function(delta = 0) 0.12
		),
		private = list(
			des_obj = NULL
		)
	)
	ToyMigrated = R6::R6Class(
		"MigrationHarnessToyMigrated",
		lock_objects = FALSE,
		public = list(
			initialize = function(des_obj, verbose = FALSE) {
				private$des_obj = des_obj
			},
			compute_estimate = function(estimate_only = FALSE) 1.25,
			compute_wald_confidence_interval = function(alpha = 0.2) c(0.5, 2.0),
			compute_wald_two_sided_pval = function(delta = 0) 0.12
		),
		private = list(
			des_obj = NULL
		)
	)

	expect_silent(expect_inference_migration_outputs_equal(
		legacy_class = ToyLegacy,
		migrated_class = ToyMigrated,
		design = list(response_type = "continuous"),
		method_calls = inference_migration_method_calls[c("estimate", "wald_ci", "wald_pval")]
	))
})

test_that("migration method snapshot helper checks capability-backed public methods", {
	before = list(
		class_name = "MigrationHarnessBefore",
		public_methods = c("compute_estimate"),
		capability_methods = character()
	)
	after = list(
		class_name = "MigrationHarnessAfter",
		public_methods = c("compute_estimate", "compute_wald_confidence_interval", "compute_wald_two_sided_pval"),
		capability_methods = c("compute_wald_confidence_interval", "compute_wald_two_sided_pval")
	)

	expect_silent(expect_inference_migration_method_snapshot_compatible(before, after))
})

test_that("migration private-state owner snapshot detects new duplicate owners", {
	Base = R6::R6Class(
		"MigrationHarnessPrivateBase",
		lock_objects = FALSE,
		private = list(
			cached_values = NULL
		)
	)
	ChildWithoutDuplicate = R6::R6Class(
		"MigrationHarnessPrivateChildWithoutDuplicate",
		lock_objects = FALSE,
		inherit = Base,
		private = list(
			model_cache = NULL
		)
	)
	ChildWithDuplicate = R6::R6Class(
		"MigrationHarnessPrivateChildWithDuplicate",
		lock_objects = FALSE,
		inherit = Base,
		private = list(
			cached_values = NULL
		)
	)

	before = inference_migration_duplicate_private_owners(ChildWithoutDuplicate)
	after = inference_migration_duplicate_private_owners(ChildWithDuplicate)
	expect_silent(expect_no_new_inference_migration_private_owner_duplicates(before, before))
	expect_identical(names(after), "cached_values")
	expect_silent(expect_no_new_inference_migration_private_owner_duplicates(before, after, allowed = "cached_values"))
})

test_that("migration dual-definition convention forbids package-level legacy aliases", {
	expect_silent(expect_inference_migration_dual_definition(
		legacy_class_name = inference_migration_legacy_name("Inference"),
		migrated_class_name = "Inference"
	))
})

test_that("migration family PR summary reports counts and newly migrated classes", {
	before = list(
		InferenceA = list(migration_status = "pending"),
		InferenceB = list(migration_status = "pending")
	)
	after = list(
		InferenceA = list(migration_status = "migrated"),
		InferenceB = list(migration_status = "pending")
	)
	summary = inference_migration_family_pr_summary(before, after)

	expect_named(summary, c("before_counts", "after_counts", "newly_migrated"))
	expect_identical(summary$newly_migrated, "InferenceA")
	expect_silent(expect_inference_migration_family_pr_summary_complete(summary))
})

test_that("custom randomization migration golden tests compare distribution and p-value", {
	custom_rand_public = list(
		fit = function(estimate_only = FALSE) {
			w = private$w
			y = private$y
			list(
				estimate = mean(y[w == 1L]) - mean(y[w == 0L]),
				model = NULL
			)
		}
	)
	CustomRandLegacy = R6::R6Class(
		"MigrationHarnessCustomRandLegacy",
		lock_objects = FALSE,
		inherit = EDI:::InferenceRand,
		public = c(
			custom_rand_public,
			list(compute_estimate = EDI:::InferenceCustomRand$public_methods$compute_estimate)
		)
	)
	CustomRandMigrated = R6::R6Class(
		"MigrationHarnessCustomRandMigrated",
		lock_objects = FALSE,
		inherit = EDI:::InferenceCustomRand,
		public = custom_rand_public
	)
	make_permutations = function(n, r) {
		list(
			w_mat = matrix(as.numeric((seq_len(n * r) + rep(seq_len(r), each = n)) %% 2L), nrow = n, ncol = r),
			m_mat = NULL
		)
	}
	method_calls_for = function(n, r = 9L) {
		permutations = make_permutations(n, r)
		list(
			randomization_distr = list(
				method = "approximate_randomization_distribution_beta_hat_T",
				args = list(r = r, show_progress = FALSE, permutations = permutations)
			),
			randomization_pval = list(
				method = "compute_rand_two_sided_pval",
				args = list(delta = 0, r = r, show_progress = FALSE, permutations = permutations)
			)
		)
	}

	builders = inference_migration_golden_design_builders()
	for (response_type in c("continuous", "incidence")) {
		design = builders[[response_type]](n = 12L, seed = 20260728L)
		expect_silent(expect_inference_migration_outputs_equal(
			legacy_class = CustomRandLegacy,
			migrated_class = CustomRandMigrated,
			design = design,
			method_calls = method_calls_for(n = 12L),
			tolerance = 1e-12
		))
	}
})

test_that("migrated custom randomization host exposes only randomization-test optional APIs", {
	EDI:::populate_inference_class_registry()
	methods = inference_migration_public_methods("InferenceCustomRand")
	randomization_methods = EDI:::inference_optional_method_names_for_capabilities("randomization_test")
	disallowed_methods = EDI:::inference_optional_method_names_for_capabilities(c(
		"randomization_ci",
		"randomization_bootstrap",
		"nonparametric_bootstrap",
		"bayesian_bootstrap",
		"jackknife"
	))

	expect_true(all(randomization_methods %in% methods))
	expect_equal(intersect(disallowed_methods, methods), character())
	expect_identical(EDI:::get_effective_components("InferenceCustomRand"), "RandomizationTest")
	expect_identical(EDI:::get_effective_capabilities("InferenceCustomRand"), "randomization_test")
	expect_silent(EDI:::mark_custom_randomization_classes_migrated("InferenceCustomRand"))
	expect_identical(
		EDI:::get_inference_hierarchy_migration_record("InferenceCustomRand")$migration_status,
		"migrated"
	)
})
