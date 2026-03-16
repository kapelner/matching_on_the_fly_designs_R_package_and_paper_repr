ordinal_cond_clogit_initialize = function(self, private, seq_des_obj, num_cores = 1, verbose = FALSE){
	assertResponseType(seq_des_obj$get_response_type(), "ordinal")
	self$super$initialize(seq_des_obj, num_cores, verbose)
	assertNoCensoring(private$any_censoring)
}


ordinal_cond_clogit_compute_setup = function(private){
	match_indic = private$match_indic
	if (is.null(match_indic)) match_indic = rep(0L, private$n)
	match_indic[is.na(match_indic)] = 0L

	strata_ids = match_indic
	reservoir_idx = which(strata_ids == 0L)
	if (length(reservoir_idx) > 0L){
		strata_ids[reservoir_idx] = max(strata_ids) + seq_along(reservoir_idx)
	}

	y_ord = as.integer(factor(private$y, ordered = TRUE))
	K = max(y_ord)

	list(
		strata_ids = strata_ids,
		y_ord = y_ord,
		K = K,
		n_alpha = K - 1L
	)
}


ordinal_cond_clogit_assert_finite_se = function(private, model_label){
	if (!is.finite(private$cached_values$s_beta_hat_T)){
		stop(model_label, ": could not compute a finite standard error.")
	}
}


ordinal_cond_clogit_shared_univ = function(private, expand_fun){
	if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

	setup = ordinal_cond_clogit_compute_setup(private)
	if (setup$K < 2L){
		private$cached_values$beta_hat_T   = NA_real_
		private$cached_values$s_beta_hat_T = NA_real_
		private$cached_values$is_z         = TRUE
		return(invisible(NULL))
	}

	expanded = expand_fun(
		as.integer(setup$y_ord),
		as.integer(private$w),
		as.integer(setup$strata_ids),
		as.integer(setup$K)
	)

	mod = clogit_helper(expanded$y, data.frame(), expanded$w, expanded$strata)
	if (is.null(mod) || !is.finite(mod$b[1]) || !is.finite(mod$ssq_b_j) || mod$ssq_b_j <= 0){
		private$cached_values$beta_hat_T   = NA_real_
		private$cached_values$s_beta_hat_T = NA_real_
		private$cached_values$is_z         = TRUE
		return(invisible(NULL))
	}

	private$cached_values$beta_hat_T   = as.numeric(mod$b[1])
	private$cached_values$s_beta_hat_T = sqrt(as.numeric(mod$ssq_b_j))
	private$cached_values$is_z         = TRUE
}


ordinal_cond_clogit_shared_multi = function(private, expand_fun, trials_fun){
	if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

	setup = ordinal_cond_clogit_compute_setup(private)
	if (setup$K < 2L){
		private$cached_values$beta_hat_T   = NA_real_
		private$cached_values$s_beta_hat_T = NA_real_
		private$cached_values$is_z         = TRUE
		return(invisible(NULL))
	}

	expanded = expand_fun(
		as.integer(setup$y_ord),
		as.integer(private$w),
		as.integer(setup$strata_ids),
		as.integer(setup$K)
	)

	X = private$get_X()
	X_stack_list = list()
	for (i in seq_len(private$n)) {
		trials_i = trials_fun(setup$y_ord[i], setup$n_alpha)
		if (length(trials_i) > 0L) {
			X_stack_list[[i]] = matrix(rep(X[i, ], each = length(trials_i)), nrow = length(trials_i))
		}
	}
	X_stack = do.call(rbind, X_stack_list)

	mod = clogit_helper(expanded$y, as.data.frame(X_stack), expanded$w, expanded$strata)
	if (is.null(mod) || !is.finite(mod$b[1]) || !is.finite(mod$ssq_b_j) || mod$ssq_b_j <= 0){
		private$cached_values$beta_hat_T   = NA_real_
		private$cached_values$s_beta_hat_T = NA_real_
		private$cached_values$is_z         = TRUE
		return(invisible(NULL))
	}

	private$cached_values$beta_hat_T   = as.numeric(mod$b[1])
	private$cached_values$s_beta_hat_T = sqrt(as.numeric(mod$ssq_b_j))
	private$cached_values$is_z         = TRUE
}
