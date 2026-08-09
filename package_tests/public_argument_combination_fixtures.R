#!/usr/bin/env Rscript

`%||%` = function(x, y) if (is.null(x)) y else x

high_priority_fixture_methods = c(
	"compute_asymp_two_sided_pval",
	"compute_asymp_confidence_interval",
	"compute_bootstrap_two_sided_pval",
	"compute_bootstrap_confidence_interval",
	"compute_bayesian_bootstrap_two_sided_pval",
	"compute_bayesian_bootstrap_confidence_interval",
	"compute_rand_two_sided_pval",
	"compute_rand_confidence_interval"
)

response_fixture_families = c(
	"continuous",
	"incidence",
	"proportion",
	"count",
	"survival",
	"ordinal"
)

fixture_tier_rank = c(smoke = 1L, ci = 2L, nightly = 3L, release = 4L)

load_edi_for_fixtures = function() {
	if (!requireNamespace("EDI", quietly = TRUE)) {
		stop("The EDI package must be installed or loadable before building public argument fixtures.", call. = FALSE)
	}
	invisible(TRUE)
}

edi_export = function(name) {
	load_edi_for_fixtures()
	getExportedValue("EDI", name)
}

make_fixture_covariates = function(n = 8L) {
	stopifnot(n >= 4L)
	data.frame(
		x1 = seq(-1, 1, length.out = n),
		x2 = rep(c(0, 1), length.out = n),
		stratum = factor(rep(c("low", "high"), each = ceiling(n / 2))[seq_len(n)]),
		cluster_id = factor(rep(seq_len(ceiling(n / 2)), each = 2L)[seq_len(n)]),
		stringsAsFactors = FALSE
	)
}

make_fixture_response = function(response_type, n = 8L, censored = FALSE) {
	y = switch(
		response_type,
		continuous = seq_len(n) / n,
		incidence = rep(c(0, 1), length.out = n),
		proportion = seq(0.2, 0.8, length.out = n),
		count = rep(0:3, length.out = n),
		survival = seq_len(n) + 0.5,
		ordinal = rep(1:3, length.out = n),
		stop("Unknown response_type: ", response_type, call. = FALSE)
	)
	dead = if (identical(response_type, "survival") && isTRUE(censored)) {
		rep(c(1, 0, 1, 1), length.out = n)
	} else {
		rep(1, n)
	}
	list(y = y, dead = dead)
}

deterministic_assignment = function(n, design_type) {
	if (design_type %in% c("clustered", "blocked_cluster")) {
		return(rep(c(0, 1), each = 2L, length.out = n))
	}
	rep(c(0, 1), length.out = n)
}

public_argument_fixture_specs = function() {
	response_specs = lapply(c("continuous", "incidence", "proportion", "count", "survival", "ordinal"), function(response_type) {
		list(
			fixture_id = paste("fixed_bernoulli", response_type, "smoke", sep = "_"),
			response_type = response_type,
			design_type = "fixed",
			class_name = "DesignFixedBernoulli",
			n = 8L,
			tier = "smoke"
		)
	})
	names(response_specs) = vapply(response_specs, function(x) x$fixture_id, character(1))
	extra_specs = list(
		survival_censored_smoke = list(
			fixture_id = "fixed_bernoulli_survival_censored_smoke",
			response_type = "survival",
			design_type = "fixed",
			class_name = "DesignFixedBernoulli",
			n = 8L,
			tier = "smoke",
			has_censoring = TRUE
		),
		sequential_bernoulli_continuous_smoke = list(
			fixture_id = "sequential_bernoulli_continuous_smoke",
			response_type = "continuous",
			design_type = "sequential",
			class_name = "DesignSeqOneByOneBernoulli",
			n = 8L,
			tier = "smoke"
		),
		fixed_blocking_continuous_smoke = list(
			fixture_id = "fixed_blocking_continuous_smoke",
			response_type = "continuous",
			design_type = "stratified",
			class_name = "DesignFixedBlocking",
			n = 8L,
			tier = "smoke",
			strata_cols = "stratum",
			B_target = 2L
		),
		fixed_cluster_continuous_smoke = list(
			fixture_id = "fixed_cluster_continuous_smoke",
			response_type = "continuous",
			design_type = "clustered",
			class_name = "DesignFixedCluster",
			n = 8L,
			tier = "smoke",
			cluster_col = "cluster_id"
		),
		fixed_blocked_cluster_continuous_smoke = list(
			fixture_id = "fixed_blocked_cluster_continuous_smoke",
			response_type = "continuous",
			design_type = "blocked_cluster",
			class_name = "DesignFixedBlockedCluster",
			n = 8L,
			tier = "smoke",
			strata_cols = "stratum",
			cluster_col = "cluster_id"
		),
		fixed_binary_match_continuous_smoke = list(
			fixture_id = "fixed_binary_match_continuous_smoke",
			response_type = "continuous",
			design_type = "matched",
			class_name = "DesignFixedBinaryMatch",
			n = 8L,
			tier = "smoke",
			m = rep(seq_len(4L), each = 2L)
		),
		fixed_greedy_continuous_smoke = list(
			fixture_id = "fixed_greedy_continuous_smoke",
			response_type = "continuous",
			design_type = "search",
			class_name = "DesignFixedGreedy",
			n = 8L,
			tier = "smoke",
			n_iter = 1L
		)
	)
	specs = c(response_specs, extra_specs)
	names(specs) = vapply(specs, function(spec) spec$fixture_id, character(1))
	specs
}

construct_public_design = function(spec, covariates) {
	constructor = edi_export(spec$class_name)
	common = list(
		response_type = spec$response_type,
		n = spec$n,
		seed = 20260809L,
		design_formula = ~ x1 + x2 + stratum
	)
	if (identical(spec$class_name, "DesignFixedBlocking")) {
		return(do.call(constructor$new, c(list(
			strata_cols = spec$strata_cols,
			B_target = spec$B_target %||% 2L,
			equal_block_sizes = TRUE
		), common)))
	}
	if (identical(spec$class_name, "DesignFixedCluster")) {
		return(do.call(constructor$new, c(list(cluster_col = spec$cluster_col), common)))
	}
	if (identical(spec$class_name, "DesignFixedBlockedCluster")) {
		return(do.call(constructor$new, c(list(strata_cols = spec$strata_cols, cluster_col = spec$cluster_col), common)))
	}
	if (identical(spec$class_name, "DesignFixedBinaryMatch")) {
		return(do.call(constructor$new, c(list(m = spec$m), common)))
	}
	if (identical(spec$class_name, "DesignFixedGreedy")) {
		return(do.call(constructor$new, c(list(n_iter = spec$n_iter %||% 1L), common)))
	}
	do.call(constructor$new, common)
}

populate_public_design = function(design, spec, covariates, response) {
	if (identical(spec$design_type, "sequential")) {
		for (i in seq_len(spec$n)) {
			design$add_one_subject_to_experiment_and_assign(covariates[i, , drop = FALSE])
			design$add_one_subject_response(i, response$y[i], response$dead[i])
		}
		return(invisible(design))
	}
	design$add_all_subjects_to_experiment(covariates)
	design$assign_w_to_all_subjects(w_precomputed = deterministic_assignment(spec$n, spec$design_type))
	design$add_all_subject_responses(response$y, deads = response$dead)
	invisible(design)
}

fixture_metadata = function(spec, covariates, design) {
	has_strata = !is.null(spec$strata_cols)
	has_cluster = !is.null(spec$cluster_col)
	has_matching = identical(spec$design_type, "matched")
	has_censoring = isTRUE(spec$has_censoring)
	list(
		fixture_id = spec$fixture_id,
		response_type = spec$response_type,
		design_type = spec$design_type,
		class_name = spec$class_name,
		n = spec$n,
		p = ncol(covariates),
		has_strata = has_strata,
		has_cluster = has_cluster,
		has_matching = has_matching,
		has_censoring = has_censoring,
		supports_jackknife = !has_cluster && spec$n >= 4L,
		supports_randomization = TRUE,
		supports_censoring = identical(spec$response_type, "survival"),
		has_finite_standard_errors = spec$n >= 4L,
		documented_nonestimable_behavior = FALSE,
		optional_backend_available = TRUE,
		capture_progress_output = FALSE,
		legal_response_transforms = legal_response_transforms(spec$response_type),
		available_columns = names(covariates),
		factor_columns = names(covariates)[vapply(covariates, is.factor, logical(1))],
		cluster_columns = if (has_cluster) spec$cluster_col else "cluster_id",
		strata_cols = spec$strata_cols %||% character(),
		cluster_col = spec$cluster_col %||% character(),
		strata_cols_must_be_factor = TRUE,
		allows_strata_cluster_overlap = FALSE,
		strata_level_count = length(unique(covariates$stratum)),
		n_exchangeable_units = if (has_cluster) length(unique(covariates[[spec$cluster_col]])) else spec$n,
		tier = spec$tier,
		compatible_method_families = high_priority_fixture_methods,
		notes = if (identical(spec$response_type, "ordinal")) "Ordinal fixture uses integer-coded ordered categories accepted by the public response API." else ""
	)
}

legal_response_transforms = function(response_type) {
	switch(
		response_type,
		continuous = c("identity", "rank", "scale"),
		incidence = c("identity", "risk"),
		proportion = c("identity", "logit"),
		count = c("identity", "log1p"),
		survival = c("identity", "log"),
		ordinal = c("identity", "rank"),
		"identity"
	)
}

build_public_argument_fixture = function(spec) {
	load_edi_for_fixtures()
	covariates = make_fixture_covariates(spec$n)
	response = make_fixture_response(spec$response_type, spec$n, censored = isTRUE(spec$has_censoring))
	design = construct_public_design(spec, covariates)
	populate_public_design(design, spec, covariates, response)
	fixture = list(
		fixture_id = spec$fixture_id,
		design = design,
		data = covariates,
		response = response$y,
		dead = response$dead,
		w = design$get_w(),
		metadata = fixture_metadata(spec, covariates, design)
	)
	validate_public_argument_fixture(fixture)
	fixture
}

build_public_argument_fixtures = function(tier = "smoke", fixture_ids = NULL) {
	specs = public_argument_fixture_specs()
	if (!is.null(fixture_ids)) {
		missing = setdiff(fixture_ids, names(specs))
		if (length(missing)) stop("Unknown fixture_ids: ", paste(missing, collapse = ", "), call. = FALSE)
		specs = specs[fixture_ids]
	}
	max_rank = fixture_tier_rank[[tier]]
	if (is.null(max_rank) || is.na(max_rank)) stop("Unknown tier: ", tier, call. = FALSE)
	specs = specs[vapply(specs, function(spec) fixture_tier_rank[[spec$tier]] <= max_rank, logical(1))]
	fixtures = lapply(specs, build_public_argument_fixture)
	names(fixtures) = vapply(fixtures, function(fixture) fixture$fixture_id, character(1))
	fixtures
}

validate_public_argument_fixture = function(fixture) {
	required_top = c("fixture_id", "design", "data", "response", "dead", "w", "metadata")
	missing_top = setdiff(required_top, names(fixture))
	if (length(missing_top)) stop("Fixture missing fields: ", paste(missing_top, collapse = ", "), call. = FALSE)
	meta = fixture$metadata
	required_meta = c(
		"fixture_id", "response_type", "design_type", "n", "p", "has_strata", "has_cluster",
		"has_matching", "has_censoring", "supports_jackknife", "supports_randomization",
		"available_columns", "tier"
	)
	missing_meta = setdiff(required_meta, names(meta))
	if (length(missing_meta)) stop("Fixture metadata missing fields: ", paste(missing_meta, collapse = ", "), call. = FALSE)
	if (!identical(fixture$fixture_id, meta$fixture_id)) stop("Fixture ID mismatch.", call. = FALSE)
	if (!identical(fixture$design$get_response_type(), meta$response_type)) stop("Response type mismatch.", call. = FALSE)
	if (!identical(as.integer(fixture$design$get_n()), as.integer(meta$n))) stop("n mismatch.", call. = FALSE)
	if (!identical(as.integer(fixture$design$get_t()), as.integer(meta$n))) stop("Fixture subjects were not fully added.", call. = FALSE)
	if (length(fixture$response) != meta$n) stop("Response length mismatch.", call. = FALSE)
	if (length(fixture$dead) != meta$n) stop("Censoring vector length mismatch.", call. = FALSE)
	if (length(fixture$w) != meta$n || any(!fixture$w %in% c(0, 1))) stop("Treatment assignment is not a legal binary vector.", call. = FALSE)
	if (!all(meta$available_columns %in% names(fixture$data))) stop("Metadata available_columns are not present in data.", call. = FALSE)
	if (isTRUE(meta$has_censoring) && !isTRUE(fixture$design$any_censoring())) stop("Censoring metadata mismatch.", call. = FALSE)
	invisible(TRUE)
}

fixture_inventory = function(fixtures) {
	rows = lapply(fixtures, function(fixture) {
		meta = fixture$metadata
		data.frame(
			fixture_id = meta$fixture_id,
			response_type = meta$response_type,
			design_type = meta$design_type,
			class_name = meta$class_name,
			n = meta$n,
			p = meta$p,
			has_strata = meta$has_strata,
			has_cluster = meta$has_cluster,
			has_matching = meta$has_matching,
			has_censoring = meta$has_censoring,
			supports_jackknife = meta$supports_jackknife,
			supports_randomization = meta$supports_randomization,
			available_columns = paste(meta$available_columns, collapse = ";"),
			tier = meta$tier,
			stringsAsFactors = FALSE
		)
	})
	out = do.call(rbind, rows)
	row.names(out) = NULL
	out
}
