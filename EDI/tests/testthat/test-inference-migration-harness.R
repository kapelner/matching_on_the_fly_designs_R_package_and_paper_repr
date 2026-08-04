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
