manual_kk_bootstrap_reference <- function(inf, B = 11L, seed = 1L) {
	inf$num_cores = 1L
	priv = inf$.__enclos_env__$private
	if (is.null(priv$cached_values$KKstats)) {
		priv$compute_basic_match_data()
	}

	n = priv$n
	y = priv$y
	dead = priv$dead
	w = priv$w
	X = priv$get_X()
	m_vec = priv$m
	if (is.null(m_vec)) {
		m_vec = rep(NA_integer_, n)
	}
	m_vec[is.na(m_vec)] = 0L
	m = priv$cached_values$KKstats$m
	i_reservoir = which(m_vec == 0L)
	n_reservoir = length(i_reservoir)

	pair_indices = vector("list", m)
	if (m > 0L) {
		for (pair_id in seq_len(m)) {
			pair_indices[[pair_id]] = which(m_vec == pair_id)
		}
	}

	set.seed(seed)
	out = numeric(B)
	for (b in seq_len(B)) {
		i_reservoir_b = sample(i_reservoir, n_reservoir, replace = TRUE)
		if (m > 0L) {
			pairs_to_include = sample(seq_len(m), m, replace = TRUE)
			i_matched_b = integer(0)
			m_vec_b_matched = integer(0)
			for (new_pair_id in seq_len(m)) {
				orig_pid = pairs_to_include[new_pair_id]
				pair_idx = pair_indices[[orig_pid]]
				i_matched_b = c(i_matched_b, pair_idx)
				m_vec_b_matched = c(m_vec_b_matched, new_pair_id, new_pair_id)
			}
		} else {
			i_matched_b = integer(0)
			m_vec_b_matched = integer(0)
		}
		i_b = c(i_reservoir_b, i_matched_b)
		m_vec_b = c(rep.int(0L, n_reservoir), m_vec_b_matched)

		boot_inf = inf$duplicate(verbose = FALSE, make_fork_cluster = FALSE)
		boot_priv = boot_inf$.__enclos_env__$private
		boot_priv$y = y[i_b]
		boot_priv$y_temp = boot_priv$y
		boot_priv$dead = dead[i_b]
		boot_priv$w = w[i_b]
		boot_priv$X = if (is.null(X)) NULL else X[i_b, , drop = FALSE]
		boot_priv$cached_values = list()
		boot_priv$m = m_vec_b
		boot_priv$compute_basic_match_data()
		if (is.function(get0("compute_reservoir_and_match_statistics", envir = boot_priv, inherits = TRUE))) {
			boot_priv$compute_reservoir_and_match_statistics()
		}
		out[b] = as.numeric(boot_inf$compute_treatment_estimate(estimate_only = TRUE))[1L]
	}
	out
}

compare_kk_bootstrap_fast_reference <- function(inf, B = 11L, seed = 1L, tolerance = 1e-10) {
	inf$num_cores = 1L
	reference_boot = manual_kk_bootstrap_reference(inf, B = B, seed = seed)
	set.seed(seed)
	fast_boot = inf$approximate_bootstrap_distribution_beta_hat_T(B = B, show_progress = FALSE)
	expect_equal(unname(fast_boot), unname(reference_boot), tolerance = tolerance)
}

test_that("KK g-computation reusable bootstrap worker matches a direct manual reference", {
	set.seed(20260422)
	n = 48L
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
	des = DesignSeqOneByOneKK14$new(n = n, response_type = "incidence", verbose = FALSE)
	for (i in seq_len(n)) {
		des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
	}
	w = des$get_w()
	p = stats::plogis(-0.5 + 0.8 * w + 0.4 * X$x1 - 0.25 * X$x2)
	y = stats::rbinom(n, 1, p)
	des$add_all_subject_responses(y)

	compare_kk_bootstrap_fast_reference(
		InferenceIncidUnivKKGCompRiskDiff$new(des, verbose = FALSE),
		seed = 301
	)
	compare_kk_bootstrap_fast_reference(
		InferenceIncidMultiKKGCompRiskDiff$new(des, verbose = FALSE),
		seed = 302
	)
	compare_kk_bootstrap_fast_reference(
		InferenceIncidUnivKKGCompRiskRatio$new(des, verbose = FALSE),
		seed = 303
	)
	compare_kk_bootstrap_fast_reference(
		InferenceIncidMultiKKGCompRiskRatio$new(des, verbose = FALSE),
		seed = 304
	)
})
