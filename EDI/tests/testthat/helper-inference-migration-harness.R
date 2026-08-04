inference_migration_with_seed = function(seed, expr) {
	old_seed = if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
		get(".Random.seed", envir = .GlobalEnv)
	} else {
		NULL
	}
	on.exit({
		if (!is.null(old_seed)) {
			assign(".Random.seed", old_seed, envir = .GlobalEnv)
		} else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
			rm(".Random.seed", envir = .GlobalEnv)
		}
	}, add = TRUE)
	set.seed(seed)
	force(expr)
}

inference_migration_add_subjects = function(des, x) {
	for (i in seq_len(nrow(x))) {
		des$add_one_subject_to_experiment_and_assign(x[i, , drop = FALSE])
	}
	des
}

inference_migration_complete_design = function(response_type, n = 12L, seed = 20260728L) {
	inference_migration_with_seed(seed, {
		x = data.frame(
			x = seq(-1.1, 1.1, length.out = n),
			z = rep(c(0L, 1L, 1L, 0L), length.out = n)
		)
		des = DesignSeqOneByOneBernoulli$new(
			n = n,
			response_type = response_type,
			verbose = FALSE
		)
		inference_migration_add_subjects(des, x)

		w01 = as.integer(des$.__enclos_env__$private$w == 1L)
		linpred = -0.2 + 0.55 * w01 + 0.25 * x$x - 0.15 * x$z
		if (identical(response_type, "continuous")) {
			y = 0.4 + 0.8 * w01 + 0.35 * x$x - 0.2 * x$z + seq(-0.3, 0.3, length.out = n)
			add_all_subject_responses_seq(des, y)
		} else if (identical(response_type, "incidence")) {
			y = as.integer(stats::plogis(linpred) > stats::quantile(stats::plogis(linpred), 0.45))
			add_all_subject_responses_seq(des, y)
		} else if (identical(response_type, "count")) {
			y = as.integer(pmax(0L, round(exp(0.6 + 0.35 * w01 + 0.2 * x$x))))
			add_all_subject_responses_seq(des, y)
		} else if (identical(response_type, "proportion")) {
			y = pmin(0.95, pmax(0.05, stats::plogis(linpred)))
			add_all_subject_responses_seq(des, y)
		} else if (identical(response_type, "ordinal")) {
			score = linpred + seq(-0.4, 0.4, length.out = n)
			y = as.integer(cut(score, breaks = c(-Inf, -0.15, 0.25, 0.65, Inf), labels = FALSE))
			add_all_subject_responses_seq(des, y)
		} else if (identical(response_type, "survival")) {
			y_lat = exp(1.2 - 0.25 * w01 + 0.1 * x$x)
			cens = exp(1.35 + 0.05 * x$z)
			add_all_subject_responses_seq(des, pmin(y_lat, cens), deads = as.integer(y_lat <= cens))
		} else {
			stop(sprintf("No migration golden fixture for response type `%s`.", response_type), call. = FALSE)
		}
		des
	})
}

inference_migration_golden_design_builders = function() {
	response_types = c("continuous", "incidence", "count", "proportion", "ordinal", "survival")
	stats::setNames(lapply(response_types, function(response_type) {
		force(response_type)
		function(n = 12L, seed = 20260728L) {
			inference_migration_complete_design(response_type, n = n, seed = seed)
		}
	}), response_types)
}

inference_migration_public_methods = function(generator) {
	if (is.character(generator)) {
		generator = get(generator, envir = asNamespace("EDI"))
	}
	public_names = character()
	current = generator
	while (!is.null(current)) {
		public_names = c(public_names, names(current$public_methods) %||% character())
		current = current$get_inherit()
	}
	sort(unique(public_names))
}

inference_migration_private_owners = function(generator) {
	if (is.character(generator)) {
		generator = get(generator, envir = asNamespace("EDI"))
	}
	owners = list()
	current = generator
	while (!is.null(current)) {
		private_names = unique(c(
			names(current$private_methods) %||% character(),
			names(current$private_fields) %||% character()
		))
		for (private_name in private_names) {
			owners[[private_name]] = c(owners[[private_name]] %||% character(), current$classname)
		}
		current = current$get_inherit()
	}
	owners
}

inference_migration_duplicate_private_owners = function(generator) {
	owners = inference_migration_private_owners(generator)
	owners[vapply(owners, length, integer(1L)) > 1L]
}

inference_migration_expected_capability_methods = function(class_name) {
	capabilities = EDI:::get_effective_capabilities(class_name)
	sort(unique(as.character(unlist(
		EDI:::public_methods_for_capability[
			intersect(names(EDI:::public_methods_for_capability), capabilities)
		],
		use.names = FALSE
	))))
}

inference_migration_method_snapshot = function(class_name) {
	list(
		class_name = class_name,
		public_methods = inference_migration_public_methods(class_name),
		capability_methods = inference_migration_expected_capability_methods(class_name)
	)
}

expect_inference_migration_method_snapshot_compatible = function(before, after) {
	testthat::expect_true(all(before$public_methods %in% after$public_methods))
	testthat::expect_true(all(after$capability_methods %in% after$public_methods))
	invisible(TRUE)
}

expect_no_new_inference_migration_private_owner_duplicates = function(before, after, allowed = character()) {
	before_dupes = names(before)
	after_dupes = names(after)
	new_dupes = setdiff(after_dupes, c(before_dupes, allowed))
	testthat::expect_equal(new_dupes, character())
	invisible(TRUE)
}

inference_migration_legacy_name = function(class_name) {
	paste0(class_name, "Legacy")
}

expect_inference_migration_dual_definition = function(legacy_class_name, migrated_class_name) {
	testthat::expect_true(grepl("Legacy$", legacy_class_name))
	testthat::expect_false(identical(legacy_class_name, migrated_class_name))
	testthat::expect_false(exists(legacy_class_name, envir = asNamespace("EDI"), inherits = FALSE))
	testthat::expect_true(exists(migrated_class_name, envir = asNamespace("EDI"), inherits = FALSE))
	invisible(TRUE)
}

inference_migration_manifest_counts = function(manifest) {
	statuses = vapply(manifest, `[[`, character(1L), "migration_status")
	as.list(table(factor(statuses, levels = c("root", "infrastructure", "pending", "migrated"))))
}

inference_migration_family_pr_summary = function(before_manifest, after_manifest) {
	before_counts = inference_migration_manifest_counts(before_manifest)
	after_counts = inference_migration_manifest_counts(after_manifest)
	newly_migrated = sort(names(after_manifest)[vapply(names(after_manifest), function(name) {
		before_status = before_manifest[[name]]$migration_status %||% NA_character_
		after_status = after_manifest[[name]]$migration_status %||% NA_character_
		!identical(before_status, "migrated") && identical(after_status, "migrated")
	}, logical(1L))])
	list(
		before_counts = before_counts,
		after_counts = after_counts,
		newly_migrated = newly_migrated
	)
}

expect_inference_migration_family_pr_summary_complete = function(summary) {
	testthat::expect_named(summary, c("before_counts", "after_counts", "newly_migrated"))
	testthat::expect_true(length(summary$newly_migrated) > 0L)
	invisible(TRUE)
}

inference_migration_method_calls = list(
	estimate = list(method = "compute_estimate", args = list()),
	asymp_ci = list(method = "compute_asymp_confidence_interval", args = list(alpha = 0.2)),
	asymp_pval = list(method = "compute_asymp_two_sided_pval", args = list(delta = 0)),
	wald_ci = list(method = "compute_wald_confidence_interval", args = list(alpha = 0.2)),
	wald_pval = list(method = "compute_wald_two_sided_pval", args = list(delta = 0)),
	score_ci = list(method = "compute_score_confidence_interval", args = list(alpha = 0.2)),
	score_pval = list(method = "compute_score_two_sided_pval", args = list(delta = 0)),
	lik_ratio_ci = list(method = "compute_lik_ratio_confidence_interval", args = list(alpha = 0.2)),
	lik_ratio_pval = list(method = "compute_lik_ratio_two_sided_pval", args = list(delta = 0)),
	gradient_ci = list(method = "compute_gradient_confidence_interval", args = list(alpha = 0.2)),
	gradient_pval = list(method = "compute_gradient_two_sided_pval", args = list(delta = 0)),
	bootstrap_distr = list(method = "approximate_bootstrap_distribution_beta_hat_T", args = list(B = 9L, show_progress = FALSE)),
	bootstrap_ci = list(method = "compute_bootstrap_confidence_interval", args = list(alpha = 0.2, B = 9L, show_progress = FALSE)),
	bootstrap_pval = list(method = "compute_bootstrap_two_sided_pval", args = list(delta = 0, B = 9L, show_progress = FALSE)),
	randomization_distr = list(method = "approximate_randomization_distribution_beta_hat_T", args = list(r = 9L, show_progress = FALSE)),
	randomization_ci = list(method = "compute_rand_confidence_interval", args = list(alpha = 0.2, r = 9L, show_progress = FALSE)),
	randomization_pval = list(method = "compute_rand_two_sided_pval", args = list(delta = 0, r = 9L, show_progress = FALSE)),
	randomization_bootstrap_distr = list(method = "approximate_rand_bootstrap_distribution_beta_hat_T", args = list(B = 9L, show_progress = FALSE)),
	randomization_bootstrap_pval = list(method = "compute_rand_bootstrap_two_sided_pval", args = list(delta = 0, B = 9L, show_progress = FALSE)),
	jackknife_estimate = list(method = "compute_jackknife_estimate", args = list()),
	jackknife_se = list(method = "compute_jackknife_std_error", args = list()),
	jackknife_ci = list(method = "compute_jackknife_wald_confidence_interval", args = list(alpha = 0.2)),
	jackknife_pval = list(method = "compute_jackknife_wald_two_sided_pval", args = list(delta = 0)),
	param_boot_pval = list(method = "compute_lik_ratio_bootstrap_two_sided_pval", args = list(delta = 0, B = 9L, show_progress = FALSE)),
	param_boot_ci = list(method = "compute_lik_ratio_bootstrap_confidence_interval", args = list(alpha = 0.2, B = 9L, show_progress = FALSE))
)

inference_migration_unsupported_error = function(err) {
	msg = conditionMessage(err)
	grepl("not implemented|not supported|only supported|does not support|does not expose|Must be implemented", msg)
}

inference_migration_call_optional_method = function(obj, method_name, args) {
	method = obj[[method_name]]
	if (!is.function(method)) {
		return(list(status = "absent", value = NULL))
	}
	tryCatch(
		list(status = "ok", value = do.call(method, args)),
		error = function(e) {
			if (inference_migration_unsupported_error(e)) {
				return(list(status = "unsupported", value = NULL, message = conditionMessage(e)))
			}
			stop(e)
		}
	)
}

inference_migration_class_generator = function(class) {
	if (is.character(class)) {
		return(get(class, envir = asNamespace("EDI")))
	}
	class
}

expect_inference_migration_outputs_equal = function(
	legacy_class,
	migrated_class,
	design,
	legacy_args = list(),
	migrated_args = list(),
	method_calls = inference_migration_method_calls,
	tolerance = 1e-7
) {
	legacy_class = inference_migration_class_generator(legacy_class)
	migrated_class = inference_migration_class_generator(migrated_class)
	legacy = do.call(legacy_class$new, c(list(des_obj = design), legacy_args))
	migrated = do.call(migrated_class$new, c(list(des_obj = design), migrated_args))
	for (label in names(method_calls)) {
		spec = method_calls[[label]]
		legacy_result = inference_migration_call_optional_method(legacy, spec$method, spec$args)
		migrated_result = inference_migration_call_optional_method(migrated, spec$method, spec$args)
		if (legacy_result$status %in% c("absent", "unsupported") &&
				migrated_result$status %in% c("absent", "unsupported")) {
			next
		}
		testthat::expect_identical(
			legacy_result$status,
			migrated_result$status,
			info = sprintf("%s status mismatch for %s", label, spec$method)
		)
		if (identical(legacy_result$status, "ok")) {
			testthat::expect_equal(
				migrated_result$value,
				legacy_result$value,
				tolerance = tolerance,
				info = sprintf("%s output mismatch for %s", label, spec$method)
			)
		}
	}
	invisible(TRUE)
}
