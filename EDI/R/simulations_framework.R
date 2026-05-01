# simulations.R
#
# Simulation framework for comparing experimental designs and inference methods
# in the EDI (Experimental Design with Inference) package.
#
# Usage:
#   devtools::load_all("EDI")
#   source("simulations.R")
#
#   designs = list(
#     DesignSeqOneByOneKK21 = list(lambda = 0.5),
#     DesignSeqOneByOneUrn  = list(alpha = 2, beta = 2)
#   )
#   inf_cls = list(
#     InferenceContinMultOLS = list(max_resample_attempts = 25L),
#     InferenceContinMultOLSKKIVWC
#   )
#   sim = SimulationFramework$new(
#     response_type    = "continuous",
#     design_classes_and_params = designs,
#     inference_classes_and_params = inf_cls,
#     n = 100, p = 5, cond_exp_func_model = "linear", Nrep = 50, betaT = 1,
#     inference_types_and_params = list(asymp_pval = list(delta = 0))
#   )
#   sim$run()
#   sim$summarize()

#' Generate Synthetic Simulation Covariates and Continuous Response
#'
#' @description
#' A helper function to generate synthetic covariates and a latent continuous response
#' identical to the logic used within \code{SimulationFramework}.  Covariates may be
#' supplied directly via \code{X_mat} or drawn randomly via \code{cov_draw_method};
#' exactly one of the two must be non-\code{NULL}.
#'
#' @param n Integer. Sample size (number of rows).
#' @param p Integer. Number of covariates (number of columns).
#' @param cond_exp_func_model Character scalar. Either \code{"linear"} (latent response is a
#'   weighted linear combination of covariates) or \code{"nonlinear"} (Friedman 1991
#'   function applied to the first five covariates; requires \code{p >= 5}).
#' @param norm_sq_beta_vec Positive numeric scalar. The desired squared Euclidean norm
#'   of the coefficient vector, i.e. \code{sum(beta^2)}.  The coefficient vector (or
#'   the overall Friedman scale) is rescaled so that this quantity equals
#'   \code{norm_sq_beta_vec}.  Default \code{1}.
#' @param X_mat Numeric matrix of dimensions \code{n x p}, or \code{NULL} (default).
#'   When supplied, this matrix is used directly as the covariate matrix and
#'   \code{cov_draw_method} must be \code{NULL}.
#' @param cov_draw_method A function used to draw \code{n * p} i.i.d. covariate
#'   values, or \code{NULL}.  The function must accept the number of draws as its
#'   first positional argument followed by any named arguments in
#'   \code{cov_draw_method_args}.  Default \code{stats::rnorm}.
#'   Must be \code{NULL} when \code{X_mat} is supplied.
#' @param cov_draw_method_args Named list of additional arguments forwarded to
#'   \code{cov_draw_method} beyond the sample-size first argument.
#'   Default is \code{list(mean = 0, sd = 1)}.
#'
#' @return A list with two elements: \code{X} (a data frame of covariates) and
#'   \code{y_cont} (a numeric vector of the latent continuous response).
#'
#' @export
generate_covariate_dataset = function(n, p,
                                      cond_exp_func_model  = c("linear", "nonlinear"),
                                      norm_sq_beta_vec     = 1,
                                      X_mat                = NULL,
                                      cov_draw_method      = stats::rnorm,
                                      cov_draw_method_args = list(mean = 0, sd = 1)) {
  cond_exp_func_model = match.arg(cond_exp_func_model)

  user_supplied_X   = !is.null(X_mat)
  user_supplied_cov = !is.null(cov_draw_method)
  if (user_supplied_X && user_supplied_cov) {
    stop("generate_covariate_dataset: supply exactly one of 'X_mat' or 'cov_draw_method', not both.")
  }
  if (!user_supplied_X && !user_supplied_cov) {
    stop("generate_covariate_dataset: one of 'X_mat' or 'cov_draw_method' must be non-NULL.")
  }

  if (user_supplied_X) {
    X_mat = as.matrix(X_mat)
  } else {
    X_mat = matrix(do.call(cov_draw_method, c(list(n * p), cov_draw_method_args)), nrow = n, ncol = p)
  }
  colnames(X_mat) = paste0("x", seq_len(p))
  X = as.data.frame(X_mat)

  if (cond_exp_func_model == "linear") {
    beta_x = seq(1, -1, length.out = p)
    beta_x = beta_x * sqrt(norm_sq_beta_vec / sum(beta_x^2))
    y_cont = as.numeric(X_mat %*% beta_x)
  } else {
    if (p < 5L) stop("Friedman nonlinear cond_exp_func_model requires p >= 5.")
    # Friedman (1991) function — first five covariates, x in [-1,1]
    beta_friedman = c(10, 20, 10, 5)
    scale_factor  = sqrt(norm_sq_beta_vec / sum(beta_friedman^2))
    y_cont = scale_factor * (
      10 * sin(pi * X_mat[, 1L] * X_mat[, 2L]) +
      20 * (X_mat[, 3L] - 0.5)^2 +
      10 *  X_mat[, 4L] +
       5 *  X_mat[, 5L]
    )
  }

  list(X = X, y_cont = y_cont)
}

#' Transform continuous latent signal to the response type scale
#'
#' @description
#' A helper function to transform a latent continuous signal to the scale
#' appropriate for a given \code{response_type}, identical to the logic used
#' within \code{SimulationFramework}.
#'
#' @param y_cont Numeric vector. The latent continuous response signal.
#' @param response_type Character scalar. One of \code{"continuous"}, \code{"incidence"},
#'   \code{"proportion"}, \code{"count"}, \code{"survival"}, \code{"ordinal"}.
#' @param n_ordinal_levels Positive integer. Number of ordinal categories when
#'   \code{response_type = "ordinal"}. Default \code{4L}.
#' @param proportion_epsilon Numeric scalar. Small value added to proportion to avoid 0 and 1. Default \code{1e-6}.
#' @param survival_min_time Numeric scalar. Minimum survival time and shift. Default \code{0.1}.
#' @param count_min_rate Integer scalar. Minimum baseline rate for count response. Default \code{0L}.
#' @param count_shift Numeric scalar. Constant added to counts after zero-centering. Default \code{0}.
#'
#' @return A numeric vector of transformed responses on the appropriate scale.
#'
#' @export
transform_cont_y_based_on_response_type = function(
  y_cont, 
  response_type, 
  n_ordinal_levels   = 4L,
  proportion_epsilon = 1e-6,
  survival_min_time  = 0.1,
  count_min_rate     = 0L,
  count_shift        = 0
) {
  y_sh = as.numeric(scale(y_cont))
  switch(response_type,
    continuous = y_sh,
    incidence  = stats::plogis(y_sh),
    proportion = {
      y_tmp = y_cont - min(y_cont) + proportion_epsilon
      y_tmp / (max(y_tmp) + proportion_epsilon)
    },
    count      = pmax(count_min_rate, round(y_cont - min(y_cont) + count_shift)),
    survival   = pmax(survival_min_time, y_sh - min(y_sh) + survival_min_time),
    ordinal    = as.integer(cut(
      y_cont,
      breaks = unique(stats::quantile(y_cont, probs = seq(0, 1, length.out = n_ordinal_levels + 1L))),
      include.lowest = TRUE
    )),
    stop("Unknown response_type: ", response_type)
  )
}

#' Simulation Framework for Experimental Designs and Inference Methods
#'
#' @description
#' An R6 class for benchmarking experimental designs and inference methods by
#' Monte Carlo simulation. Each replication generates synthetic covariates and
#' responses, runs every requested \code{(design, inference)} pair, and records
#' point estimates, confidence intervals, and p-values. Aggregated metrics are
#' available from \code{$summarize()}.
#'
#' @details
#' Covariates are drawn independently from \eqn{\mathrm{Uniform}(0, 1)}.
#' \itemize{
#'   \item \code{cond_exp_func_model = "linear"}: the base continuous signal is
#'     \eqn{y = X\beta} where \eqn{\beta} is evenly spaced from 1 to \eqn{-1}.
#'   \item \code{cond_exp_func_model = "nonlinear"}: the Friedman (1991) function
#'     \eqn{10\sin(\pi x_1 x_2) + 20(x_3-0.5)^2 + 10x_4 + 5x_5};
#'     requires \eqn{p \ge 5}.
#' }
#' The continuous base signal is transformed to the scale appropriate for
#' \code{response_type}.  Treatment effects are applied per-subject: additive
#' on the linear/logit/ordinal scale, log-multiplicative for count and survival.
#'
#' For each \code{(design, inference)} pair the framework runs whichever of the
#' following are supported by the inference class:
#' \itemize{
#'   \item \strong{asymptotic} (\code{InferenceAsymp} subclasses): Wald CI and
#'     p-value.
#'   \item \strong{bootstrap} (\code{InferenceBoot} subclasses): percentile CI
#'     and p-value.
#'   \item \strong{randomisation} (\code{InferenceRand} subclasses): p-value;
#'     additionally a test-inversion CI for \code{continuous},
#'     \code{proportion}, and \code{count} response types
#'     (\code{InferenceRandCI} subclasses).
#' }
#' Incompatible \code{(design, inference)} pairs (e.g.\ a KK-specific inference
#' class with a non-KK design) are silently skipped via \code{tryCatch}.
#'
#' Reported summary metrics include:
#' \itemize{
#'   \item \strong{MSE}: \eqn{\overline{(\hat\beta_T - \beta_T)^2}} over reps
#'     with a finite point estimate.
#'   \item \strong{coverage}: proportion of reps where \eqn{\beta_T} lies
#'     inside the CI (\code{NA} when no CI is available for that inference type).
#'   \item \strong{power}: proportion of p-values \eqn{< \alpha}; equals the
#'     empirical type-I error rate when \code{betaT = 0}.
#' }
#'
#' @export
SimulationFramework = R6::R6Class("SimulationFramework",
  lock_objects = FALSE,

  # ── public ─────────────────────────────────────────────────────────────────
  public = list(

    #' @description
    #' Create a new \code{SimulationFramework}.
    #'
    #' @param response_type \strong{(required)} Character scalar.  The type of
    #'   outcome variable.  One of \code{"continuous"}, \code{"incidence"},
    #'   \code{"proportion"}, \code{"count"}, \code{"survival"},
    #'   \code{"ordinal"}.
    #'
    #' @param design_classes_and_params \code{NULL} (default) or a list
    #'   describing design classes and optional constructor parameters.
    #'   Unnamed R6 class generators use default parameters, for example
    #'   \code{list(DesignSeqOneByOneKK21, FixedDesignBernoulli)}.  Named
    #'   entries use the entry name as the design class and the value as the
    #'   parameter list, for example
    #'   \code{list(DesignSeqOneByOneUrn = list(alpha = 2, beta = 2))}.
    #'   Duplicate named entries are allowed for repeated designs with
    #'   different parameters.
    #'   Each generator must be constructable with only \code{response_type} and
    #'   \code{n} plus any extra params supplied in this list.
    #'   \code{NULL} uses the package's standard design set.
    #'   Designs requiring \code{strata_cols}, \code{cluster_col}, or
    #'   \code{factors} have sensible defaults auto-injected (first covariate
    #'   column; second for \code{cluster_col}; \code{list(treatment=2)} for
    #'   \code{factors}) when not supplied in the parameter list.
    #'   Example:
    #'   \preformatted{design_classes_and_params = list(
    #'   DesignSeqOneByOneKK21 = list(lambda = 0.5, t_0_pct = 0.1),
    #'   DesignSeqOneByOneUrn  = list(alpha  = 2,   beta    = 2),
    #'   FixedDesignBernoulli  # default params
    #' )}
    #'   Commonly useful design constructor parameters:
    #'   \describe{
    #'     \item{\code{lambda}}{Matching-weight decay for
    #'       \code{KK14} / \code{KK21} / \code{KK21stepwise}.}
    #'     \item{\code{t_0_pct}}{Burn-in fraction for
    #'       \code{KK14} / \code{KK21} / \code{KK21stepwise}.}
    #'     \item{\code{morrison}}{Logical; Morrison correction for \code{KK14}.}
    #'     \item{\code{alpha}, \code{beta}}{Shape parameters for
    #'       \code{DesignSeqOneByOneUrn}.}
    #'     \item{\code{preferred_num_bins_for_continuous_covariate}}{Bin count for
    #'       \code{FixedDesignBlocking} and \code{FixedDesignBlockedCluster}.}
    #'     \item{\code{B_target}}{Target number of blocks for \code{FixedDesignBlocking}.}
   #'   }
    #'
    #' @param inference_classes_and_params \code{NULL} (default) or a list
    #'   describing inference classes and optional constructor parameters.
    #'   Unnamed R6 class generators use default parameters, for example
    #'   \code{list(InferenceContinMultOLS, InferenceContinMultOLSKKIVWC)}.
    #'   Named entries use the entry name as the inference class and the value
    #'   as the constructor parameter list, for example
    #'   \code{list(InferenceContinMultOLS = list(max_resample_attempts = 25L))}.
    #'   Duplicate named entries are allowed for repeated inference classes with
    #'   different parameters. Supplied parameters must be accepted by the
    #'   inference class constructor.
    #'   \code{NULL} selects a curated set for the given
    #'   \code{response_type}: several universal classes that work with any
    #'   design, plus representative KK-specific classes (silently skipped for
    #'   non-KK designs at runtime).
    #'
    #' @param n Integer.  Sample size per simulation replication.  Default
    #'   \code{100}.
    #'
    #' @param p Integer.  Number of covariates.  Must be \eqn{\ge 5} when
    #'   \code{cond_exp_func_model = "nonlinear"}.  Default \code{5}.
    #'
    #' @param cond_exp_func_model Character scalar or vector.  How the latent continuous signal is
    #'   constructed before transformation to the \code{response_type} scale.
    #'   \describe{
    #'     \item{\code{"linear"}}{Linear combination \eqn{X\beta} with
    #'       coefficients evenly spaced from 1 to \eqn{-1}.}
    #'     \item{\code{"nonlinear"}}{Friedman (1991) function
    #'       \eqn{10\sin(\pi x_1 x_2)+20(x_3-0.5)^2+10x_4+5x_5};
    #'       requires \eqn{p \ge 5}.}
    #'   }
    #'   Default \code{"linear"}.
    #'
    #' @param Nrep Positive integer.  Number of Monte Carlo replications.
    #'   Default \code{100}.
    #'
    #' @param betaT Numeric scalar.  True treatment effect added to treated
    #'   subjects' outcomes.  The scale is response-type specific: additive for
    #'   \code{continuous}, \code{proportion}, and \code{ordinal}; on the logit
    #'   scale for \code{incidence}; log-multiplicative for \code{count} and
    #'   \code{survival}.  Default \code{1}.  Set \code{betaT = 0} to check
    #'   type-I error.
    #'
    #' @param alpha Numeric in \eqn{(0,1)}.  Significance level used for all
    #'   confidence intervals and for computing power (\eqn{p < \alpha}).
    #'   Default \code{0.05}.
    #'
    #' @param B_boot Positive integer.  Bootstrap resamples per CI / p-value
    #'   call.  Default \code{201}.
    #'
    #' @param r_rand Positive integer.  Randomisation draws per rand p-value
    #'   call, and per bisection step of the rand CI.  Default \code{201}.
    #'
    #' @param pval_epsilon Numeric.  Bisection convergence tolerance for
    #'   randomisation-based CIs (\code{compute_rand_confidence_interval}).
    #'   Default \code{0.02}.
    #'
    #' @param sd_noise Numeric \eqn{> 0}. Standard deviation of independent
    #'   Gaussian noise added to each subject's outcome.  Default \code{1}.
    #'
    #' @param n_ordinal_levels Positive integer. Number of ordinal categories when
    #'   \code{response_type = "ordinal"}. Default \code{4L}.
    #'
	    #' @param proportion_epsilon Numeric scalar. Small value added to proportion
	    #'   base responses to avoid 0 and 1. Default \code{1e-6}.
	    #'
	    #' @param phi_proportion Positive numeric scalar. Precision parameter for
	    #'   beta-distributed observed proportion outcomes. The beta mean is
	    #'   \code{y_linear_model[i] + betaT * w[i]}. Default \code{100}.
	    #'
	    #' @param k_survival Positive numeric scalar. Scale parameter passed to
	    #'   the Weibull draw for observed survival outcomes. Default \code{2}.
	    #'
	    #' @param incidence_clamp Numeric scalar in \eqn{(0, 0.5)}. Clamp applied
	    #'   to the Bernoulli probability for observed incidence outcomes.
	    #'   Default \code{1e-9}.
	    #'
	    #' @param proportion_clamp Numeric scalar in \eqn{(0, 0.5)}. Clamp applied
	    #'   to the beta mean for observed proportion outcomes. Default \code{1e-9}.
	    #'
	    #' @param count_clamp Positive numeric scalar. Minimum Poisson mean for
	    #'   observed count outcomes. Default \code{1e-9}.
	    #'
	    #' @param survival_clamp Positive numeric scalar. Minimum Weibull shape for
	    #'   observed survival outcomes. Default \code{1e-9}.
	    #'
	    #' @param survival_min_time Numeric scalar. Minimum survival time and shift
    #'   for base responses. Default \code{0.1}.
    #'
    #' @param count_min_rate Integer scalar. Minimum baseline rate for count
    #'   responses. Default \code{0L}.
    #'
    #' @param count_shift Numeric scalar. Constant added to counts after
    #'   zero-centering for base responses. Default \code{0}.
    #'
    #' @param norm_sq_beta_vec Positive numeric scalar. The desired squared
    #'   Euclidean norm of the latent linear coefficient vector \eqn{\beta}.
    #'   The generated vector is scaled to match this norm. Default \code{1}.
    #'
    #' @param X_mat Numeric matrix of dimensions \code{n x p}, or \code{NULL} (default).
    #'   If provided, these fixed covariates are used for every replication.
    #'   In this case, \code{cov_draw_method} must be \code{NULL}.
    #'
    #' @param seed Integer or \code{NULL} (default).  Random seed for the
    #'   entire simulation run.
    #'
    #' @param cov_draw_method A function used to draw \code{n * p} i.i.d. covariate
    #'   values for every replication. The function must accept the total number
    #'   of values as its first argument, followed by arguments in
    #'   \code{cov_draw_method_args}.  Default \code{stats::rnorm}.
    #'   Must be \code{NULL} when \code{X_mat} is supplied.
    #'
    #' @param cov_draw_method_args Named list of additional arguments forwarded to
    #'   \code{cov_draw_method} beyond the sample-size first argument. Default is
    #'   \code{list(mean = 0, sd = 1)}.
    #'
    #' @param random_X_draws Logical. If \code{TRUE} (default), a new set of
    #'   covariates is drawn for every single replication. If \code{FALSE},
    #'   one set is drawn per \code{(n, p)} cell and shared across its
    #'   replications.
    #'
    #' @param prob_censoring Numeric in \eqn{[0,1]}.  Per-subject independent
    #'   censoring probability; applied only when
    #'   \code{response_type = "survival"}.  Default \code{0.2}.
    #'
    #' @param verbose Logical.  If \code{TRUE}, prints a message for every
    #'   replication and for every \code{(design, inference)} pair that is
    #'   skipped due to an error.  Default \code{FALSE}.
    #'
    #' @param keep_all_intermediate_data Logical. If \code{TRUE}, the framework
    #'   saves the instantiated design and inference objects for every replication.
    #'   These can be retrieved after the run using \code{$get_all_intermediate_data()}.
    #'   Warning: this can consume a lot of memory for many replications.
    #'   Default \code{FALSE}.
    #'
    #' @param turn_off_asserts_for_speed Logical. If \code{TRUE} (default),
    #'   all \pkg{checkmate} assertions across the package are globally
    #'   disabled during the simulation run to improve performance.
    #'
    #' @param results_filename Character scalar.  The filename for the CSV results
    #'   file.  Default \code{"simulation_framework_results.csv"}.
    #'
    #' @param continue_from_last_result_row Logical. If \code{TRUE} (default),
    #'   the framework loads existing results from \code{results_filename} and
    #'   skips previously completed replications.
    #'
    #' @param inference_types_and_params \code{NULL} (default) or a named list
    #'   from inference type to a named list of arguments for that type's function
    #'   invocation.  The list names control which inference outputs are computed.
    #'   Valid names are \code{"asymp_ci"}, \code{"asymp_pval"},
    #'   \code{"exact_ci"}, \code{"exact_pval"}, \code{"boot_ci"},
    #'   \code{"boot_pval"}, \code{"rand_ci"}, and \code{"rand_pval"}.
    #'   Each value must be a named list whose names are accepted by the
    #'   corresponding inference function.  \code{NULL} runs all eight types with
    #'   default invocation arguments.
    #'   Example:
    #'   \preformatted{inference_types_and_params = list(
    #'   asymp_pval = list(delta = 0),
    #'   boot_ci    = list(B = 99, type = "perc"),
    #'   rand_pval  = list(r = 999, transform_responses = TRUE)
    #' )}
    #'   When no \code{*_ci} type is requested, \code{coverage} is omitted from
    #'   \code{$summarize()}.  When no \code{*_pval} type is requested,
    #'   \code{power} is omitted.
    initialize = function(
      response_type,
      design_classes_and_params = NULL,
      inference_classes_and_params = NULL,
      n                     = 100L,
      p                     = 5L,
      cond_exp_func_model   = "linear",
      Nrep                  = 100L,
      betaT                 = 1,
      alpha                 = 0.05,
      B_boot                = 201L,
      r_rand                = 201L,
      pval_epsilon          = 0.02,
      sd_noise              = 1,
      n_ordinal_levels      = 4L,
      proportion_epsilon    = 1e-6,
      phi_proportion        = 100,
      k_survival            = 2,
      incidence_clamp       = 1e-9,
      proportion_clamp      = 1e-9,
      count_clamp           = 1e-9,
      survival_clamp        = 1e-9,
      survival_min_time     = 0.1,
      count_min_rate        = 0L,
      count_shift           = 0,
      norm_sq_beta_vec      = 1,
      X_mat                 = NULL,
      seed                  = NULL,
      cov_draw_method       = stats::rnorm,
      cov_draw_method_args  = list(mean = 0, sd = 1),
      random_X_draws        = TRUE,
      prob_censoring        = 0.25,
      verbose                    = TRUE,
      keep_all_intermediate_data = FALSE,
      turn_off_asserts_for_speed = TRUE,
      inference_types_and_params = NULL,
      results_filename      = "simulation_framework_results.csv",
      continue_from_last_result_row = TRUE
    ) {
      valid_rt = c("continuous", "incidence", "proportion",
                   "count", "survival", "ordinal")
      if (!response_type %in% valid_rt)
        stop("response_type must be one of: ", paste(valid_rt, collapse = ", "))

      n_values = unique(as.integer(n))
      p_values = unique(as.integer(p))
      betaT_values = unique(as.numeric(betaT))
      cond_exp_func_model_values = unique(as.character(cond_exp_func_model))

      if (length(n_values) == 0L || any(!is.finite(n_values)) || any(n_values <= 1L))
        stop("n must contain finite integers greater than 1")
      if (length(p_values) == 0L || any(!is.finite(p_values)) || any(p_values < 1L))
        stop("p must contain finite positive integers")
      if (length(betaT_values) == 0L || any(!is.finite(betaT_values)))
        stop("betaT must contain finite numeric values")
      if (length(cond_exp_func_model_values) == 0L ||
          any(!cond_exp_func_model_values %in% c("linear", "nonlinear")))
        stop("cond_exp_func_model must contain only 'linear' and/or 'nonlinear'")
      if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)))
        stop("seed must be NULL or one finite numeric value")
      if (!is.null(X_mat) && (length(n_values) > 1L || length(p_values) > 1L))
        stop("X_mat can only be used when n and p are scalar")
      if (!isTRUE(random_X_draws) && is.null(seed))
        stop("random_X_draws = FALSE requires seed to be non-NULL")
	      if (!is.finite(phi_proportion) || phi_proportion <= 0)
	        stop("phi_proportion must be finite and > 0")
	      if (!is.finite(k_survival) || k_survival <= 0)
	        stop("k_survival must be finite and > 0")
	      if (!is.finite(incidence_clamp) || incidence_clamp <= 0 || incidence_clamp >= 0.5)
	        stop("incidence_clamp must be finite and in (0, 0.5)")
	      if (!is.finite(proportion_clamp) || proportion_clamp <= 0 || proportion_clamp >= 0.5)
	        stop("proportion_clamp must be finite and in (0, 0.5)")
	      if (!is.finite(count_clamp) || count_clamp <= 0)
	        stop("count_clamp must be finite and > 0")
	      if (!is.finite(survival_clamp) || survival_clamp <= 0)
	        stop("survival_clamp must be finite and > 0")

      valid_inf_types = c("asymp_ci", "asymp_pval", "exact_ci", "exact_pval",
                          "boot_ci",  "boot_pval",  "rand_ci",  "rand_pval")
      inf_type_spec = private$.parse_inference_types_and_params(
        inference_types_and_params,
        valid_inf_types
      )
      inf_types = names(inf_type_spec)

      private$response_type    = response_type
      private$n_values         = n_values
      private$p_values         = p_values
      private$cond_exp_func_model_values = cond_exp_func_model_values
      private$Nrep             = as.integer(Nrep)
      private$betaT_values     = betaT_values
      private$alpha            = alpha
      private$B_boot           = as.integer(B_boot)
      private$r_rand           = as.integer(r_rand)
      private$pval_epsilon     = pval_epsilon
	      private$sd_noise             = sd_noise
	      private$n_ordinal_levels     = as.integer(n_ordinal_levels)
	      private$proportion_epsilon    = proportion_epsilon
	      private$phi_proportion        = phi_proportion
	      private$k_survival            = k_survival
	      private$incidence_clamp       = incidence_clamp
	      private$proportion_clamp      = proportion_clamp
	      private$count_clamp           = count_clamp
	      private$survival_clamp        = survival_clamp
	      private$survival_min_time     = survival_min_time
      private$count_min_rate        = as.integer(count_min_rate)
      private$count_shift           = count_shift
      private$norm_sq_beta_vec     = norm_sq_beta_vec
      private$X_mat                = X_mat
      private$seed                 = if (is.null(seed)) NULL else as.integer(seed)
      private$cov_draw_method      = cov_draw_method
      private$cov_draw_method_args = cov_draw_method_args
      private$random_X_draws       = random_X_draws
      private$prob_censoring       = prob_censoring
      private$verbose                    = verbose
      private$turn_off_asserts_for_speed = turn_off_asserts_for_speed
      private$keep_all_intermediate_data = keep_all_intermediate_data
      private$results_filename     = results_filename
      private$continue_from_last_result_row = continue_from_last_result_row
      private$inf_types        = inf_types
      private$inference_type_params = inf_type_spec
      private$param_grid       = private$.build_param_grid(
        n_values,
        p_values,
        betaT_values,
        cond_exp_func_model_values
      )
      private$current_n        = private$param_grid$n[[1L]]
      private$current_p        = private$param_grid$p[[1L]]
      private$current_betaT    = private$param_grid$betaT[[1L]]
      private$current_cond_exp_func_model = private$param_grid$cond_exp_func_model[[1L]]

      design_spec = private$.parse_design_classes_and_params(
        design_classes_and_params,
        parent.frame()
      )
      private$design_classes = design_spec$classes
      private$design_params  = design_spec$params

      inference_spec = private$.parse_inference_classes_and_params(
        inference_classes_and_params,
        parent.frame()
      )
      private$inference_classes = inference_spec$classes
      private$inference_constructor_params = inference_spec$params

      n_des = length(private$design_classes)
      n_inf = length(private$inference_classes)

      private$design_labels    = private$.compute_design_labels()
      private$inference_labels = private$.compute_inference_labels()

      private$raw_results = list()
      private$has_run     = FALSE
    },

    # ── run() ─────────────────────────────────────────────────────────────────
    #' @description
    #' Execute the simulation replications.
    #'
    #' @return The \code{SimulationFramework} object itself (invisibly).
    run = function() {
      if (!is.null(private$seed)) {
        had_seed = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        if (had_seed) {
          old_seed = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
          on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
        } else {
          on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
        }
        set.seed(private$seed)
      }

      # Disable assertions for the duration of the simulation for speed
      if (private$turn_off_asserts_for_speed){
        toggle_asserts(FALSE)
        on.exit(toggle_asserts(TRUE))        
      }

      n_des = length(private$design_classes)
      n_inf = length(private$inference_classes)
      n_met = length(private$inf_types)
      n_cells = nrow(private$param_grid)
      existing_results = private$.load_existing_results()
      n_existing = nrow(existing_results)
      
      # Pre-allocate results list to maximum possible size to avoid re-allocations
      private$raw_results = vector("list", n_existing + private$Nrep * n_des * n_inf * n_met * n_cells)
      private$results_idx = 0L
      private$seen_result_keys = character(0L)
      if (n_existing > 0L) {
        existing_rows = lapply(seq_len(n_existing), function(i) as.list(existing_results[i]))
        private$raw_results[seq_len(n_existing)] = existing_rows
        private$results_idx = n_existing
        private$seen_result_keys = vapply(existing_rows, private$.result_key_from_row, "", USE.NAMES = FALSE)
      }

      private$all_intermediate_data = vector("list", private$Nrep * n_cells)
      private$valid_combos          = list()
      private$seen_combo_keys      = character(0L)
      private$exact_warned_classes = character(0L)

      n_des = length(private$design_classes)
      n_inf = length(private$inference_classes)

      if (isTRUE(private$verbose)) {
        message(sprintf(
          "simulations: CEF_mod=%s  n=%s  p=%s  Nrep=%d  betaT=%s designs=%d inferences=%d",
          private$.format_values(private$cond_exp_func_model_values),
          private$.format_values(private$n_values),
          private$.format_values(private$p_values),
          private$Nrep,
          private$.format_values(private$betaT_values), 
          n_des, 
          n_inf
        ))
      }

      log_interval = max(1L, private$Nrep %/% 10L)
      shared_X_draws = list()
      if (!isTRUE(private$random_X_draws)) {
        np_grid = unique(private$param_grid[, .(n, p)])
        for (np_idx in seq_len(nrow(np_grid))) {
          n_i = np_grid$n[[np_idx]]
          p_i = np_grid$p[[np_idx]]
          X_i = if (is.null(private$X_mat)) {
            matrix(
              do.call(private$cov_draw_method, c(list(n_i * p_i), private$cov_draw_method_args)),
              nrow = n_i,
              ncol = p_i
            )
          } else {
            as.matrix(private$X_mat)
          }
          colnames(X_i) = paste0("x", seq_len(p_i))
          shared_X_draws[[paste(n_i, p_i, sep = "|||")]] = X_i
        }
      }

      had_seed_plan = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      if (had_seed_plan) {
        old_seed_plan = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      }
      planned_combos = list()
      for (cell_idx in seq_len(n_cells)) {
        private$current_n = private$param_grid$n[[cell_idx]]
        private$current_p = private$param_grid$p[[cell_idx]]
        private$current_betaT = private$param_grid$betaT[[cell_idx]]
        private$current_cond_exp_func_model = private$param_grid$cond_exp_func_model[[cell_idx]]
        rep_data = if (isTRUE(private$random_X_draws)) {
          private$.generate_data()
        } else {
          private$.generate_data_from_X(
            shared_X_draws[[paste(private$current_n, private$current_p, sep = "|||")]]
          )
        }
        cell_combos = private$.build_valid_combos_for_current_cell(rep_data)
        if (length(cell_combos) > 0L)
          planned_combos = c(planned_combos, cell_combos)
      }
      if (had_seed_plan) {
        assign(".Random.seed", old_seed_plan, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
      private$valid_combos = planned_combos
      private$seen_combo_keys = unique(vapply(
        planned_combos,
        function(combo) paste(combo$cond_exp_func_model, combo$n, combo$p, format(combo$betaT, digits = 17L), combo$design, combo$inference, sep = "|||"),
        ""
      ))
      private$progress_total = length(planned_combos) * private$Nrep
      planned_result_keys = unlist(lapply(seq_len(private$Nrep), function(rep_i) {
        vapply(planned_combos, function(combo) {
          private$.result_key_for_values(
            private$response_type,
            combo$cond_exp_func_model,
            combo$n,
            combo$p,
            combo$betaT,
            rep_i,
            combo$design,
            combo$inference,
            combo$inference_type
          )
        }, "")
      }), use.names = FALSE)
      private$progress_count = sum(private$seen_result_keys %in% planned_result_keys)
      private$progress_log_interval = max(1L, private$progress_total %/% 20L)
      private$progress_bar = NULL
      private$use_progress_bar = FALSE
      private$total_cells = n_cells
      if (isTRUE(private$verbose)) {
        message(sprintf("%d / %d runs skipped", private$progress_count, private$progress_total))
        private$use_progress_bar = TRUE
        if (private$progress_total > 0L) {
           # We don't use txtProgressBar anymore, just our dual bar
           # But we need a newline at the end
           on.exit(cat("\n"), add = TRUE)
        }
      }

      for (cell_idx in seq_len(n_cells)) {
        private$current_cell_idx = cell_idx
        private$current_n = private$param_grid$n[[cell_idx]]
        private$current_p = private$param_grid$p[[cell_idx]]
        private$current_betaT = private$param_grid$betaT[[cell_idx]]
        private$current_cond_exp_func_model = private$param_grid$cond_exp_func_model[[cell_idx]]

        # Pre-calculate tasks per rep for this cell
        private$tasks_per_rep = sum(vapply(private$valid_combos, function(c) {
          c$n == private$current_n && c$p == private$current_p && 
          format(c$betaT, digits = 17L) == format(private$current_betaT, digits = 17L) && 
          c$cond_exp_func_model == private$current_cond_exp_func_model
        }, logical(1L)))

        if (isTRUE(private$verbose) && !private$use_progress_bar)
          message(sprintf(
            "  Cell %d / %d: cond_exp_func_model=%s  n=%d  p=%d  betaT=%g",
            cell_idx, n_cells, private$current_cond_exp_func_model,
            private$current_n, private$current_p, private$current_betaT
          ))

        for (rep in seq_len(private$Nrep)) {
          private$current_rep_idx = rep
          private$current_task_in_rep_idx = 0L
          if (isTRUE(private$use_progress_bar)) {
            private$.draw_triple_progress_bar()
          } else if (isTRUE(private$verbose) && (rep %% log_interval == 0L || private$Nrep == 1L)) {
            message(sprintf("    Rep %d / %d", rep, private$Nrep))
          }

          rep_data = if (isTRUE(private$random_X_draws)) {
            private$.generate_data()
          } else {
            private$.generate_data_from_X(
              shared_X_draws[[paste(private$current_n, private$current_p, sep = "|||")]]
            )
          }
          X = rep_data$X
          y_linear_model = rep_data$y_linear_model
          true_mean_diff_ate = private$compute_true_mean_diff_ate(y_linear_model)
          rep_slot = (cell_idx - 1L) * private$Nrep + rep

          for (di in seq_along(private$design_classes)) {
            design_gen   = private$design_classes[[di]]
            design_name  = private$design_labels[[di]]
            design_extra = if (!is.null(private$design_params)) private$design_params[[di]] else list()

            des_obj = tryCatch(
              private$.build_design(design_gen, X, y_linear_model, design_extra),
              error = function(e) {
                NULL
              }
            )
            if (private$keep_all_intermediate_data) {
              private$all_intermediate_data[[rep_slot]]$n = private$current_n
              private$all_intermediate_data[[rep_slot]]$p = private$current_p
              private$all_intermediate_data[[rep_slot]]$betaT = private$current_betaT
              private$all_intermediate_data[[rep_slot]]$cond_exp_func_model = private$current_cond_exp_func_model
              private$all_intermediate_data[[rep_slot]]$y_linear_model = y_linear_model
              private$all_intermediate_data[[rep_slot]]$designs[[design_name]] = des_obj
            }
            if (is.null(des_obj)) next

            for (ii in seq_along(private$inference_classes)) {
              inf_gen  = private$inference_classes[[ii]]
              inf_name = private$inference_labels[[ii]]
              inf_ctor_extra = private$inference_constructor_params[[ii]]

              inf_obj = tryCatch({
                do.call(inf_gen$new, c(list(des_obj), inf_ctor_extra))
              }, error = function(e) {
                NULL
              })
              if (is.null(inf_obj)) next
              if (private$keep_all_intermediate_data)
                private$all_intermediate_data[[rep_slot]]$inferences[[design_name]][[inf_name]] = inf_obj

              te = if (is(inf_obj, "InferenceAllSimpleMeanDiff")) true_mean_diff_ate else private$current_betaT

              valid_inference_types = private$.valid_inference_types(inf_obj)

              combo_key = paste(
                private$current_cond_exp_func_model,
                private$current_n,
                private$current_p,
                format(private$current_betaT, digits = 17L),
                design_name,
                inf_name,
                sep = "|||"
              )
              pending_inference_types = valid_inference_types[!vapply(
                valid_inference_types,
                function(it) private$.result_key(rep, design_name, inf_name, it) %in% private$seen_result_keys,
                logical(1L)
              )]
              skipped_inference_types = setdiff(valid_inference_types, pending_inference_types)
              if (length(skipped_inference_types) > 0L) {
                for (it in skipped_inference_types) {
                  private$.log_skip(rep, design_name, inf_name, it)
                  private$current_task_in_rep_idx = private$current_task_in_rep_idx + 1L
                }
                if (isTRUE(private$use_progress_bar)) private$.draw_triple_progress_bar()
              }
              if (length(pending_inference_types) == 0L) next

              est = tryCatch({
                v = inf_obj$compute_estimate()
                if (is.null(v) || length(v) == 0L) NA_real_ else as.numeric(v)[1L]
              }, error = function(e) NA_real_)

              if (is(inf_obj, "InferenceAsymp") && private$.any_inf_type(c("asymp_ci", "asymp_pval"))) {
                if ("asymp_pval" %in% pending_inference_types) {
                  asymp_pval_args = private$.args_for_inf_type(inf_obj, "asymp_pval")
                  pval_a = tryCatch(
                    do.call(inf_obj$compute_asymp_two_sided_pval, asymp_pval_args),
                    error = function(e) NA_real_
                  )
                  private$.record(rep, design_name, inf_name, "asymp_pval", est, c(NA_real_, NA_real_), pval_a, te)
                }
                if ("asymp_ci" %in% pending_inference_types) {
                  asymp_ci_args = private$.args_for_inf_type(inf_obj, "asymp_ci", list(alpha = private$alpha))
                  ci_a = tryCatch(
                    do.call(inf_obj$compute_asymp_confidence_interval, asymp_ci_args),
                    error = function(e) c(NA_real_, NA_real_)
                  )
                  private$.record(rep, design_name, inf_name, "asymp_ci", est, ci_a, NA_real_, te)
                }
              }

              if (private$.any_inf_type(c("exact_ci", "exact_pval"))) {
                if (!is(inf_obj, "InferenceExact")) {
                  if (isTRUE(private$verbose) && !inf_name %in% private$exact_warned_classes) {
                    warning(sprintf(
                      "'%s' does not inherit InferenceExact; exact_* inference will be skipped for this class.",
                      inf_name))
                    private$exact_warned_classes = c(private$exact_warned_classes, inf_name)
                  }
                } else {
                  if ("exact_pval" %in% pending_inference_types) {
                    exact_pval_args = private$.args_for_inf_type(inf_obj, "exact_pval")
                    pval_e = tryCatch(
                      do.call(inf_obj$compute_exact_two_sided_pval_for_treatment_effect, exact_pval_args),
                      error = function(e) NA_real_
                    )
                    private$.record(rep, design_name, inf_name, "exact_pval", est, c(NA_real_, NA_real_), pval_e, te)
                  }
                  if ("exact_ci" %in% pending_inference_types) {
                    exact_ci_args = private$.args_for_inf_type(inf_obj, "exact_ci", list(alpha = private$alpha))
                    ci_e = tryCatch(
                      do.call(inf_obj$compute_exact_confidence_interval, exact_ci_args),
                      error = function(e) c(NA_real_, NA_real_)
                    )
                    private$.record(rep, design_name, inf_name, "exact_ci", est, ci_e, NA_real_, te)
                  }
                }
              }

              if (is(inf_obj, "InferenceBoot") && private$.any_inf_type(c("boot_ci", "boot_pval"))) {
                if ("boot_pval" %in% pending_inference_types) {
                  boot_pval_args = private$.args_for_inf_type(inf_obj, "boot_pval", list(B = private$B_boot, na.rm = TRUE))
                  pval_b = tryCatch(
                    do.call(inf_obj$compute_bootstrap_two_sided_pval, boot_pval_args),
                    error = function(e) NA_real_
                  )
                  private$.record(rep, design_name, inf_name, "boot_pval", est, c(NA_real_, NA_real_), pval_b, te)
                }
                if ("boot_ci" %in% pending_inference_types) {
                  boot_ci_args = private$.args_for_inf_type(
                    inf_obj,
                    "boot_ci",
                    list(B = private$B_boot, alpha = private$alpha, na.rm = TRUE, show_progress = FALSE)
                  )
                  ci_b = tryCatch(
                    do.call(inf_obj$compute_bootstrap_confidence_interval, boot_ci_args),
                    error = function(e) c(NA_real_, NA_real_)
                  )
                  private$.record(rep, design_name, inf_name, "boot_ci", est, ci_b, NA_real_, te)
                }
              }

              if (is(inf_obj, "InferenceRand") && private$.any_inf_type(c("rand_ci", "rand_pval"))) {
                if ("rand_pval" %in% pending_inference_types) {
                  rand_pval_args = private$.args_for_inf_type(
                    inf_obj,
                    "rand_pval",
                    list(r = private$r_rand, na.rm = TRUE, show_progress = FALSE)
                  )
                  pval_r = tryCatch(
                    do.call(inf_obj$compute_rand_two_sided_pval, rand_pval_args),
                    error = function(e) NA_real_
                  )
                  private$.record(rep, design_name, inf_name, "rand_pval", est, c(NA_real_, NA_real_), pval_r, te)
                }
                if (
                  "rand_ci" %in% pending_inference_types &&
                  is(inf_obj, "InferenceRandCI") &&
                  private$response_type %in% c("continuous", "proportion", "count")
                ) {
                  rand_ci_args = private$.args_for_inf_type(
                    inf_obj,
                    "rand_ci",
                    list(r = private$r_rand, alpha = private$alpha, pval_epsilon = private$pval_epsilon, show_progress = FALSE)
                  )
                  ci_r = tryCatch(
                    do.call(inf_obj$compute_rand_confidence_interval, rand_ci_args),
                    error = function(e) c(NA_real_, NA_real_)
                  )
                  private$.record(rep, design_name, inf_name, "rand_ci", est, ci_r, NA_real_, te)
                }
              }
            }
          }

        }
      }

      private$has_run = TRUE
      invisible(self)
    },

    # ── get_all_intermediate_data() ───────────────────────────────────────────
    #' @description
    #' Return all intermediate data saved during \code{$run()}, or \code{NULL}
    #' if \code{keep_all_intermediate_data} was \code{FALSE} (the default).
    #' Throws an error if \code{$run()} has not been called yet.
    #'
    #' @return A list of length \code{Nrep}.  Each element is a named list with
    #'   three entries:
    #'   \describe{
    #'     \item{\code{y_linear_model}}{Numeric vector of base responses for that replication.}
    #'     \item{\code{designs}}{Named list of instantiated design objects, one per
    #'       configured design class.}
    #'     \item{\code{inferences}}{Named list of lists of instantiated inference
    #'       objects, indexed first by design class name then by inference class name.}
    #'   }
    #'   Returns \code{NULL} when \code{keep_all_intermediate_data = FALSE}.
    get_all_intermediate_data = function() {
      if (!private$has_run) stop("Call $run() first.")
      if (!private$keep_all_intermediate_data) return(NULL)
      private$all_intermediate_data
    },

    # ── clear_all_intermediate_data_and_gc() ─────────────────────────────────
    #' @description
    #' Release all stored intermediate data and invoke the garbage collector.
    #' Useful after inspecting intermediate results to free memory before
    #' further processing.  Sets the internal store to \code{NULL} and calls
    #' \code{gc()}.
    #'
    #' @return The \code{SimulationFramework} object itself (invisibly).
    clear_all_intermediate_data_and_gc = function() {
      private$all_intermediate_data = NULL
      gc()
      invisible(self)
    },

    # ── get_results() ─────────────────────────────────────────────────────────
    #' @description
    #' Get the raw results of the simulation.
    #'
    #' @return A \code{data.table} containing one row per (replication, design,
    #'   inference class, and inference type.
    get_results = function() {
      if (!private$has_run) stop("Call $run() first.")
      if (private$results_idx == 0L)
        return(data.table::data.table(
          rep = integer(), cond_exp_func_model = character(), n = integer(),
          p = integer(), betaT = numeric(), design = character(), inference = character(),
          inference_type = character(), estimate = numeric(),
          ci_lo = numeric(), ci_hi = numeric(), pval = numeric()
        ))
      # Prune pre-allocated list to only include what was actually recorded
      data.table::rbindlist(private$raw_results[seq_len(private$results_idx)])
    },

    # ── summarize() ───────────────────────────────────────────────────────────
    #' @description
    #' Aggregate and summarize simulation results.
    #'
    #' @return A \code{data.table} with aggregated metrics (MSE, coverage,
    #'   power) grouped by design, inference class, and inference type.
    summarize = function() {
      if (!private$has_run) stop("Call $run() first.")

      # ── Reference grid of all valid (design, inference, inference_type) combos ────────
      if (length(private$valid_combos) == 0L) {
        message("No results."); return(invisible(NULL))
      }
      ref_grid = data.table::rbindlist(lapply(private$valid_combos, as.list))

      # ── Per-class params strings ───────────────────────────────────────────────
      inf_names = private$inference_labels
      design_params_map = stats::setNames(
        lapply(seq_along(private$design_classes), function(di)
          private$.params_to_str(if (!is.null(private$design_params)) private$design_params[[di]] else NULL)),
        private$design_labels)
      inf_params_map = stats::setNames(
        lapply(seq_along(private$inference_classes), function(ii)
          private$.params_to_str(private$inference_constructor_params[[ii]])),
        inf_names)

      ref_grid[, design_params    := unlist(design_params_map[design])]
      ref_grid[, inference_params := unlist(inf_params_map[inference])]
      ref_grid[, inference_type_params := vapply(inference_type, private$.params_for_inference_type_to_str, "")]

      # ── Aggregate raw results ─────────────────────────────────────────────────
      dt         = self$get_results()
      alpha      = private$alpha
      report_cov = any(grepl("_ci$",   private$inf_types))
      report_pow = any(grepl("_pval$", private$inf_types))

      if (nrow(dt) > 0L) {
        agg = dt[, {
          est_ok = .SD[is.finite(estimate)]
          ci_ok  = .SD[is.finite(ci_lo) & is.finite(ci_hi)]
          pv_ok  = pval[is.finite(pval)]
          row = list(
            MSE   = if (nrow(est_ok)) mean((est_ok$estimate - est_ok$true_estimand)^2) else NA_real_,
            n_est = nrow(est_ok)
          )
          if (report_cov) {
            row$coverage  = if (nrow(ci_ok)) mean(ci_ok$ci_lo <= ci_ok$true_estimand & ci_ok$true_estimand <= ci_ok$ci_hi) else NA_real_
            row$n_cov     = nrow(ci_ok)
            row$ci_length = if (nrow(ci_ok)) mean(ci_ok$ci_hi - ci_ok$ci_lo) else NA_real_
          }
          if (report_pow) {
            row$power = if (length(pv_ok)) mean(pv_ok < alpha) else NA_real_
            row$n_pow = length(pv_ok)
          }
          row
        }, by = .(cond_exp_func_model, n, p, betaT, design, inference, inference_type),
           .SDcols = c("estimate", "ci_lo", "ci_hi", "pval", "true_estimand")]
      } else {
        agg = data.table::data.table(cond_exp_func_model = character(), n = integer(), p = integer(),
                                     betaT = numeric(), design = character(), inference = character(),
                                     inference_type = character(), power = numeric(), MSE = numeric(),
                                     n_est = integer(), n_pow = integer())
      }

      # ── Right-join: every valid combo appears, NA for those with no data ──────
      result = agg[ref_grid, on = .(cond_exp_func_model, n, p, betaT, design, inference, inference_type)]
      
      # Ensure n_est and n_pow are present and replace NA with 0
      if (!"n_est" %in% names(result)) result[, n_est := 0L]
      if (!"n_pow" %in% names(result)) result[, n_pow := 0L]
      result[is.na(n_est), n_est := 0L]
      result[is.na(n_pow), n_pow := 0L]
      
      result[order(cond_exp_func_model, betaT, n, p, design, inference, inference_type)]
    },

    # ── print() ───────────────────────────────────────────────────────────────
    #' @description
    #' Print a summary of the \code{SimulationFramework} configuration and status.
    print = function() {
      cat("SimulationFramework\n")
      cat("  response_type :", private$response_type, "\n")
      cat("  cond_exp_func_model :", private$.format_values(private$cond_exp_func_model_values), "\n")
      cat("  n / p         :", private$.format_values(private$n_values), "/", private$.format_values(private$p_values), "\n")
      cat("  Nrep / betaT  :", private$Nrep, "/", private$.format_values(private$betaT_values), "\n")
      cat("  alpha / B_boot / r_rand :",
          private$alpha, "/", private$B_boot, "/", private$r_rand, "\n")
      cat("  inference_types:", paste(private$inf_types, collapse = ", "), "\n")
      if (any(vapply(private$design_params, length, integer(1L)) > 0L))
        cat("  design_classes_and_params: (", length(private$design_params), " designs)\n", sep = "")
      if (any(vapply(private$inference_constructor_params, length, integer(1L)) > 0L))
        cat("  inference_classes_and_params: (", length(private$inference_constructor_params), " inference classes)\n", sep = "")
      if (any(vapply(private$inference_type_params, length, integer(1L)) > 0L))
        cat("  inference_types_and_params: (", length(private$inference_type_params), " inference types)\n", sep = "")
      design_names = vapply(private$design_classes,  function(g) g$classname, "")
      inf_names    = private$inference_labels
      cat("  Designs (", length(design_names), "):",
          paste(design_names, collapse = ", "), "\n")
      cat("  Inference (", length(inf_names), "):",
          paste(inf_names,    collapse = ", "), "\n")
      if (private$has_run) {
        cat("  Status : completed\n")
        sm = self$summarize()
        if (!is.null(sm) && nrow(sm) > 0L) {
          cat(sprintf("\nSummary  (alpha = %g):\n", private$alpha))
          print(sm)
        }
      } else {
        cat("  Status : not yet run  (call $run())\n")
      }
      invisible(self)
    }
  ),

  # ── private ────────────────────────────────────────────────────────────────
  private = list(
    response_type    = NULL,
    n_values         = NULL,
    p_values         = NULL,
    cond_exp_func_model_values = NULL,
    Nrep             = NULL,
    betaT_values     = NULL,
    param_grid       = NULL,
    current_n        = NULL,
    current_p        = NULL,
    current_cond_exp_func_model = NULL,
    current_betaT    = NULL,
    alpha            = NULL,
    B_boot           = NULL,
    r_rand           = NULL,
    pval_epsilon     = NULL,
    sd_noise             = NULL,
	    n_ordinal_levels     = NULL,
	    proportion_epsilon   = NULL,
	    phi_proportion       = NULL,
	    k_survival           = NULL,
	    incidence_clamp      = NULL,
	    proportion_clamp     = NULL,
	    count_clamp          = NULL,
	    survival_clamp       = NULL,
	    survival_min_time    = NULL,
    count_min_rate       = NULL,
    count_shift          = NULL,
    norm_sq_beta_vec     = NULL,
    X_mat                = NULL,
    seed                 = NULL,
    cov_draw_method      = NULL,
    cov_draw_method_args = NULL,
    random_X_draws       = NULL,
    prob_censoring       = NULL,
    verbose          = NULL,
    results_filename = NULL,
    continue_from_last_result_row = NULL,
    design_params    = NULL,
    inference_constructor_params = NULL,
    inference_type_params = NULL,
    inf_types        = NULL,
    design_classes   = NULL,
    inference_classes = NULL,
    design_labels    = NULL,
    inference_labels = NULL,
    turn_off_asserts_for_speed = NULL,
    raw_results               = NULL,
    results_idx               = 0L,
    all_intermediate_data     = NULL,
    keep_all_intermediate_data = FALSE,
    has_run                   = FALSE,
    exact_warned_classes      = NULL,
    valid_combos              = NULL,
    seen_combo_keys  = NULL,
    seen_result_keys = NULL,
    total_cells              = 0L,
    current_cell_idx         = 0L,
    current_rep_idx          = 0L,
    current_task_in_rep_idx  = 0L,
    tasks_per_rep            = 0L,
    progress_total           = 0L,
    progress_count           = 0L,
    progress_bar             = NULL,
    use_progress_bar         = FALSE,
    progress_log_interval    = 0L,

    # ── Design spec parsing ───────────────────────────────────────────────────

    .parse_design_classes_and_params = function(spec, eval_env) {
      if (is.null(spec)) {
        classes = private$.default_design_classes()
        return(list(classes = classes, params = lapply(classes, function(...) list())))
      }
      if (!is.list(spec))
        stop("design_classes_and_params must be NULL or a list")

      nm = names(spec)
      if (is.null(nm)) nm = rep("", length(spec))

      classes = vector("list", length(spec))
      params  = vector("list", length(spec))

      for (i in seq_along(spec)) {
        entry_name = nm[[i]]
        entry      = spec[[i]]

        if (inherits(entry, "R6ClassGenerator")) {
          classes[[i]] = entry
          params[[i]]  = list()
          next
        }

        if (!nzchar(entry_name)) {
          stop(
            "design_classes_and_params[[", i, "]] must be an R6 class generator ",
            "or a named parameter list whose name is the design class"
          )
        }

        cls = private$.resolve_design_class(entry_name, eval_env)
        if (!inherits(cls, "R6ClassGenerator"))
          stop("design class '", entry_name, "' is not an R6 class generator")

        if (is.null(entry)) {
          entry = list()
        } else if (!is.list(entry)) {
          stop("design_classes_and_params[['", entry_name, "']] must be a list of parameters")
        }

        classes[[i]] = cls
        params[[i]]  = entry
      }

      list(classes = classes, params = params)
    },

    .resolve_design_class = function(class_name, eval_env) {
      if (exists(class_name, envir = eval_env, inherits = TRUE))
        return(get(class_name, envir = eval_env, inherits = TRUE))
      ns = asNamespace("EDI")
      if (exists(class_name, envir = ns, inherits = FALSE))
        return(get(class_name, envir = ns, inherits = FALSE))
      stop("could not find design class '", class_name, "'")
    },

    # ── Inference spec parsing ────────────────────────────────────────────────

    .parse_inference_classes_and_params = function(spec, eval_env) {
      if (is.null(spec)) {
        classes = private$.default_inference_classes()
        return(list(classes = classes, params = lapply(classes, function(...) list())))
      }
      if (!is.list(spec))
        stop("inference_classes_and_params must be NULL or a list")

      nm = names(spec)
      if (is.null(nm)) nm = rep("", length(spec))

      classes = vector("list", length(spec))
      params  = vector("list", length(spec))

      for (i in seq_along(spec)) {
        entry_name = nm[[i]]
        entry      = spec[[i]]

        if (inherits(entry, "R6ClassGenerator")) {
          classes[[i]] = entry
          params[[i]]  = list()
          next
        }

        if (!nzchar(entry_name)) {
          stop(
            "inference_classes_and_params[[", i, "]] must be an R6 class generator ",
            "or a named parameter list whose name is the inference class"
          )
        }

        cls = private$.resolve_inference_class(entry_name, eval_env)
        if (!inherits(cls, "R6ClassGenerator"))
          stop("inference class '", entry_name, "' is not an R6 class generator")

        if (is.null(entry)) {
          entry = list()
        } else if (!is.list(entry)) {
          stop("inference_classes_and_params[['", entry_name, "']] must be a list of parameters")
        }

        private$.validate_r6_init_args(cls, entry, "inference_classes_and_params")
        classes[[i]] = cls
        params[[i]]  = entry
      }

      list(classes = classes, params = params)
    },

    .resolve_inference_class = function(class_name, eval_env) {
      if (exists(class_name, envir = eval_env, inherits = TRUE))
        return(get(class_name, envir = eval_env, inherits = TRUE))
      ns = asNamespace("EDI")
      if (exists(class_name, envir = ns, inherits = FALSE))
        return(get(class_name, envir = ns, inherits = FALSE))
      stop("could not find inference class '", class_name, "'")
    },

    .validate_r6_init_args = function(r6gen, args, arg_name) {
      if (length(args) == 0L) return(invisible(TRUE))
      if (is.null(names(args)) || any(!nzchar(names(args)))) {
        stop(arg_name, " for ", r6gen$classname,
             " must be a named list of constructor arguments")
      }
      init_fn = private$.get_r6_init_fn(r6gen)
      if (is.null(init_fn)) {
        stop(r6gen$classname, " has no discoverable initialize() constructor; ",
             arg_name, " parameters cannot be validated")
      }
      fn_formals = names(formals(init_fn))
      if ("..." %in% fn_formals) return(invisible(TRUE))
      bad = setdiff(names(args), fn_formals)
      if (length(bad)) {
        stop(
          arg_name, " for ", r6gen$classname,
          " contains constructor argument(s) not accepted by initialize(): ",
          paste(bad, collapse = ", ")
        )
      }
      invisible(TRUE)
    },

    # ── Inference type parsing and invocation args ────────────────────────────

    .parse_inference_types_and_params = function(spec, valid_inf_types) {
      if (is.null(spec)) {
        return(stats::setNames(lapply(valid_inf_types, function(...) list()),
                               valid_inf_types))
      }
      if (!is.list(spec) || is.null(names(spec)) || any(!nzchar(names(spec)))) {
        stop("inference_types_and_params must be NULL or a named list")
      }
      bad = setdiff(names(spec), valid_inf_types)
      if (length(bad)) {
        stop("Invalid inference_types_and_params names: ", paste(bad, collapse = ", "),
             ".  Valid values: ", paste(valid_inf_types, collapse = ", "))
      }
      spec = spec[!duplicated(names(spec))]
      for (inf_type in names(spec)) {
        if (is.null(spec[[inf_type]])) {
          spec[[inf_type]] = list()
        } else if (!is.list(spec[[inf_type]])) {
          stop("inference_types_and_params[['", inf_type, "']] must be a named list")
        } else if (length(spec[[inf_type]]) > 0L &&
                   (is.null(names(spec[[inf_type]])) || any(!nzchar(names(spec[[inf_type]]))))) {
          stop("inference_types_and_params[['", inf_type, "']] must be a named list")
        }
      }
      spec
    },

    .has_inf_type = function(inf_type) {
      inf_type %in% private$inf_types
    },

    .any_inf_type = function(inf_types) {
      any(inf_types %in% private$inf_types)
    },

    .inf_type_method_name = function(inf_type) {
      switch(inf_type,
        asymp_ci   = "compute_asymp_confidence_interval",
        asymp_pval = "compute_asymp_two_sided_pval",
        exact_ci   = "compute_exact_confidence_interval",
        exact_pval = "compute_exact_two_sided_pval_for_treatment_effect",
        boot_ci    = "compute_bootstrap_confidence_interval",
        boot_pval  = "compute_bootstrap_two_sided_pval",
        rand_ci    = "compute_rand_confidence_interval",
        rand_pval  = "compute_rand_two_sided_pval",
        stop("Unknown inference type: ", inf_type)
      )
    },

    .args_for_inf_type = function(inf_obj, inf_type, defaults = list()) {
      user_args = private$inference_type_params[[inf_type]]
      if (is.null(user_args)) user_args = list()
      method_name = private$.inf_type_method_name(inf_type)
      private$.validate_method_args(inf_obj, method_name, user_args, inf_type)
      modifyList(defaults, user_args)
    },

    .validate_method_args = function(inf_obj, method_name, args, inf_type) {
      if (length(args) == 0L) return(invisible(TRUE))
      fn = tryCatch(inf_obj[[method_name]], error = function(e) NULL)
      if (!is.function(fn)) {
        stop("Cannot validate parameters for ", inf_type, ": function ",
             method_name, "() is not available on ", class(inf_obj)[1L])
      }
      fn_formals = names(formals(fn))
      if ("..." %in% fn_formals) return(invisible(TRUE))
      bad = setdiff(names(args), fn_formals)
      if (length(bad)) {
        stop("inference_types_and_params[['", inf_type, "']] contains argument(s) ",
             "not accepted by ", method_name, "(): ", paste(bad, collapse = ", "))
      }
      invisible(TRUE)
    },

    .params_for_inference_type_to_str = function(inference_type) {
      ps = private$.params_to_str(private$inference_type_params[[inference_type]])
      if (nchar(ps) > 0L) paste0(inference_type, "(", ps, ")") else ""
    },

    # ── R6 formals helpers ────────────────────────────────────────────────────

    # Walk an R6 class generator's inheritance chain to find the first class
    # that defines initialize(), and return that function.
    .get_r6_init_fn = function(r6gen) {
      gen = r6gen
      while (!is.null(gen)) {
        fn = tryCatch(gen$public_methods$initialize, error = function(e) NULL)
        if (!is.null(fn)) return(fn)
        gen = tryCatch(gen$get_inherit(), error = function(e) NULL)
      }
      NULL
    },

    # Filter an arg list to only those accepted by r6gen's initialize().
    # If initialize accepts '...' or cannot be found, all args pass through.
    .filter_by_formals = function(r6gen, args) {
      if (length(args) == 0L) return(args)
      init_fn = private$.get_r6_init_fn(r6gen)
      if (is.null(init_fn)) return(args)
      fn_formals = names(formals(init_fn))
      if ("..." %in% fn_formals) return(args)
      args[names(args) %in% fn_formals]
    },

    # Serialize a named params list to "k=v, ..." string.
    .params_to_str = function(p) {
      if (is.null(p) || length(p) == 0L) return("")
      kv = mapply(function(k, v) paste0(k, "=", paste(deparse(v), collapse = "")),
                  names(p), p, SIMPLIFY = TRUE)
      paste(kv, collapse = ", ")
    },

    .build_param_grid = function(n_values, p_values, betaT_values, cond_exp_func_model_values) {
      grid = data.table::as.data.table(expand.grid(
        cond_exp_func_model = cond_exp_func_model_values,
        n = n_values,
        p = p_values,
        betaT = betaT_values,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      ))
      grid = grid[!(cond_exp_func_model == "nonlinear" & p < 5L)]
      if (nrow(grid) == 0L)
        stop("No valid simulation cells remain after filtering cond_exp_func_model / p combinations")
      grid
    },

    .format_values = function(x) {
      if (length(x) == 1L) as.character(x) else paste0("c(", paste(x, collapse = ", "), ")")
    },

    .load_existing_results = function() {
      empty_dt = data.table::data.table(
        response_type = character(),
        cond_exp_func_model = character(),
        n = integer(),
        p = integer(),
        betaT = numeric(),
        rep = integer(),
        design = character(),
        inference = character(),
        inference_type = character(),
        estimate = numeric(),
        ci_lo = numeric(),
        ci_hi = numeric(),
        pval = numeric(),
        true_estimand = numeric()
      )
      if (!isTRUE(private$continue_from_last_result_row) || !file.exists(private$results_filename))
        return(empty_dt)

      dt = data.table::fread(private$results_filename)
      if (!"response_type" %in% names(dt))
        dt[, response_type := private$response_type]
      dt = dt[response_type == private$response_type]

      for (nm in names(empty_dt)) {
        if (!nm %in% names(dt))
          dt[, (nm) := empty_dt[[nm]]]
      }
      dt[, response_type := as.character(response_type)]
      dt[, cond_exp_func_model := as.character(cond_exp_func_model)]
      dt[, design := as.character(design)]
      dt[, inference := as.character(inference)]
      dt[, inference_type := as.character(inference_type)]
      dt[, n := as.integer(n)]
      dt[, p := as.integer(p)]
      dt[, rep := as.integer(rep)]
      dt[, betaT := as.numeric(betaT)]
      dt[, estimate := as.numeric(estimate)]
      dt[, ci_lo := as.numeric(ci_lo)]
      dt[, ci_hi := as.numeric(ci_hi)]
      dt[, pval := as.numeric(pval)]
      dt[, true_estimand := as.numeric(true_estimand)]
      dt[, names(empty_dt), with = FALSE]
    },

    .result_key_for_values = function(response_type, cond_exp_func_model, n, p, betaT, rep, design, inference, inference_type) {
      paste(
        response_type,
        cond_exp_func_model,
        n,
        p,
        format(betaT, digits = 17L),
        rep,
        design,
        inference,
        inference_type,
        sep = "|||"
      )
    },

    .result_key = function(rep, design, inference, inference_type) {
      private$.result_key_for_values(
        private$response_type,
        private$current_cond_exp_func_model,
        private$current_n,
        private$current_p,
        private$current_betaT,
        rep,
        design,
        inference,
        inference_type
      )
    },

    .result_key_from_row = function(row) {
      private$.result_key_for_values(
        row$response_type,
        row$cond_exp_func_model,
        row$n,
        row$p,
        row$betaT,
        row$rep,
        row$design,
        row$inference,
        row$inference_type
      )
    },

    .result_metadata_dt = function(rep, design, inference, inference_type) {
      data.table::data.table(
        response_type = private$response_type,
        cond_exp_func_model = private$current_cond_exp_func_model,
        n = as.integer(private$current_n),
        p = as.integer(private$current_p),
        betaT = as.numeric(private$current_betaT),
        rep = as.integer(rep),
        design = design,
        inference = inference,
        inference_type = inference_type
      )
    },

    .log_skip = function(rep, design, inference, inference_type) {
      invisible(NULL)
    },

    .valid_inference_types = function(inf_obj) {
      valid_inference_types = character(0L)
      if (is(inf_obj, "InferenceAsymp")) {
        valid_inference_types = c(
          valid_inference_types,
          intersect(private$inf_types, c("asymp_ci", "asymp_pval"))
        )
      }
      if (is(inf_obj, "InferenceExact")) {
        valid_inference_types = c(
          valid_inference_types,
          intersect(private$inf_types, c("exact_ci", "exact_pval"))
        )
      }
      if (is(inf_obj, "InferenceBoot")) {
        valid_inference_types = c(
          valid_inference_types,
          intersect(private$inf_types, c("boot_ci", "boot_pval"))
        )
      }
      if (is(inf_obj, "InferenceRand")) {
        valid_inference_types = c(
          valid_inference_types,
          intersect(private$inf_types, "rand_pval")
        )
        if (is(inf_obj, "InferenceRandCI") &&
            private$response_type %in% c("continuous", "proportion", "count")) {
          valid_inference_types = c(
            valid_inference_types,
            intersect(private$inf_types, "rand_ci")
          )
        }
      }
      valid_inference_types
    },

    .build_valid_combos_for_current_cell = function(rep_data) {
      X = rep_data$X
      y_linear_model = rep_data$y_linear_model
      combos = list()
      for (di in seq_along(private$design_classes)) {
        design_gen   = private$design_classes[[di]]
        design_name  = private$design_labels[[di]]
        design_extra = if (!is.null(private$design_params)) private$design_params[[di]] else list()
        des_obj = tryCatch(
          private$.build_design(design_gen, X, y_linear_model, design_extra),
          error = function(e) NULL
        )
        if (is.null(des_obj)) next
        for (ii in seq_along(private$inference_classes)) {
          inf_gen  = private$inference_classes[[ii]]
          inf_name = private$inference_labels[[ii]]
          inf_ctor_extra = private$inference_constructor_params[[ii]]
          inf_obj = tryCatch(
            do.call(inf_gen$new, c(list(des_obj), inf_ctor_extra)),
            error = function(e) NULL
          )
          if (is.null(inf_obj)) next
          valid_inference_types = private$.valid_inference_types(inf_obj)
          if (length(valid_inference_types) == 0L) next
          for (it in valid_inference_types) {
            combos[[length(combos) + 1L]] = list(
              cond_exp_func_model = private$current_cond_exp_func_model,
              n = private$current_n,
              p = private$current_p,
              betaT = private$current_betaT,
              design = design_name,
              inference = inf_name,
              inference_type = it
            )
          }
        }
      }
      combos
    },

    .advance_progress = function() {
      private$progress_count = private$progress_count + 1L
      private$current_task_in_rep_idx = private$current_task_in_rep_idx + 1L
      if (!isTRUE(private$verbose)) return(invisible(NULL))
      
      if (isTRUE(private$use_progress_bar)) {
        private$.draw_triple_progress_bar()
      } else if (private$progress_log_interval > 0L &&
                 (private$progress_count %% private$progress_log_interval == 0L ||
                  private$progress_count == private$progress_total)) {
        message(sprintf("Completed %d / %d runs", private$progress_count, private$progress_total))
      }
      invisible(NULL)
    },

    .draw_triple_progress_bar = function() {
      width = getOption("width", 80L)
      if (is.null(width) || width < 80L) width = 80L
      
      # Proportions
      cell_prop = if (private$total_cells > 0) (private$current_cell_idx - 1) / private$total_cells else 0
      rep_prop  = if (private$Nrep > 0) (private$current_rep_idx - 1) / private$Nrep else 0
      task_prop = if (private$tasks_per_rep > 0) private$current_task_in_rep_idx / private$tasks_per_rep else 0
      
      # Labels
      label_cell = sprintf("DGP:%d/%d", private$current_cell_idx, private$total_cells)
      label_rep  = sprintf("Rep:%d/%d", private$current_rep_idx, private$Nrep)
      label_task = sprintf("Des/Inf:%d/%d", private$current_task_in_rep_idx, private$tasks_per_rep)
      
      # Total label width
      total_label_width = nchar(label_cell) + nchar(label_rep) + nchar(label_task)
      
      # Available width for 3 bars and separators
      # Separators: 5 spaces + 6 brackets = 11 spaces. Plus \r = 12 total overhead.
      overhead = 11L
      bar_space = width - total_label_width - overhead
      if (bar_space < 15L) {
        bar_width1 = bar_width2 = bar_width3 = 5L
      } else {
        bar_width1 = bar_space %/% 3L
        bar_width2 = bar_space %/% 3L
        bar_width3 = bar_space - bar_width1 - bar_width2
      }
      
      make_bar = function(prop, b_width) {
        pct_str = sprintf(" %3d%% ", floor(prop * 100))
        n_pct = nchar(pct_str)
        fill = floor(prop * b_width)
        
        full_bar = paste0(strrep("=", fill), strrep(" ", b_width - fill))
        if (b_width >= n_pct) {
          start_pos = (b_width - n_pct) %/% 2 + 1
          substr(full_bar, start_pos, start_pos + n_pct - 1) = pct_str
        }
        sprintf("[%s]", full_bar)
      }
      
      msg = sprintf("\r%s %s %s %s %s %s", 
                    label_cell, make_bar(cell_prop, bar_width1),
                    label_rep,  make_bar(rep_prop,  bar_width2),
                    label_task, make_bar(task_prop, bar_width3))
      
      message(substr(msg, 1, width), appendLF = FALSE)
      if (exists("flush.console")) utils::flush.console()
    },

    .append_result_row_to_file = function(row) {
      file_exists = file.exists(private$results_filename)
      data.table::fwrite(row, private$results_filename, append = file_exists, col.names = !file_exists)
    },

    # Build unique per-instance design labels: "ClassName (params)" with
    # " [k]" suffix when two instances would otherwise share the same label.
    .compute_design_labels = function() {
      labels = vapply(seq_along(private$design_classes), function(di) {
        cls = private$design_classes[[di]]$classname
        ps  = private$.params_to_str(
          if (!is.null(private$design_params)) private$design_params[[di]] else NULL)
        if (nchar(ps) > 0L) paste0(cls, " (", ps, ")") else cls
      }, "")
      # Append [k] for any group of labels that are still identical
      for (lbl in unique(labels[duplicated(labels)])) {
        idx = which(labels == lbl)
        for (k in seq_along(idx)) labels[idx[k]] = paste0(lbl, " [", k, "]")
      }
      labels
    },

    .compute_inference_labels = function() {
      labels = vapply(seq_along(private$inference_classes), function(ii) {
        cls = private$inference_classes[[ii]]$classname
        ps  = private$.params_to_str(private$inference_constructor_params[[ii]])
        if (nchar(ps) > 0L) paste0(cls, " (", ps, ")") else cls
      }, "")
      for (lbl in unique(labels[duplicated(labels)])) {
        idx = which(labels == lbl)
        for (k in seq_along(idx)) labels[idx[k]] = paste0(lbl, " [", k, "]")
      }
      labels
    },

    # ── Data generation ───────────────────────────────────────────────────────

    .generate_data = function() {
      data = generate_covariate_dataset(
        n                    = private$current_n,
        p                    = private$current_p,
        cond_exp_func_model  = private$current_cond_exp_func_model,
        norm_sq_beta_vec     = private$norm_sq_beta_vec,
        X_mat                = private$X_mat,
        cov_draw_method      = private$cov_draw_method,
        cov_draw_method_args = private$cov_draw_method_args
      )
      data$y_linear_model = transform_cont_y_based_on_response_type(
        y_cont             = data$y_cont,
        response_type      = private$response_type,
        n_ordinal_levels   = private$n_ordinal_levels,
        proportion_epsilon = private$proportion_epsilon,
        survival_min_time  = private$survival_min_time,
        count_min_rate     = private$count_min_rate,
        count_shift        = private$count_shift
      )
      data$y_cont = NULL # SimulationFramework doesn't need the raw cont y anymore
      data
    },

    .generate_data_from_X = function(X_mat) {
      data = generate_covariate_dataset(
        n                    = private$current_n,
        p                    = private$current_p,
        cond_exp_func_model  = private$current_cond_exp_func_model,
        norm_sq_beta_vec     = private$norm_sq_beta_vec,
        X_mat                = X_mat,
        cov_draw_method      = NULL,
        cov_draw_method_args = private$cov_draw_method_args
      )
      data$y_linear_model = transform_cont_y_based_on_response_type(
        y_cont             = data$y_cont,
        response_type      = private$response_type,
        n_ordinal_levels   = private$n_ordinal_levels,
        proportion_epsilon = private$proportion_epsilon,
        survival_min_time  = private$survival_min_time,
        count_min_rate     = private$count_min_rate,
        count_shift        = private$count_shift
      )
      data$y_cont = NULL # SimulationFramework doesn't need the raw cont y anymore
      data
    },

    compute_true_mean_diff_ate = function(y_linear_model) {
      eta_c = y_linear_model
      eta_t = y_linear_model + private$current_betaT

      clamp = function(x, lo, hi) {
        pmin(hi, pmax(lo, x))
      }

      switch(private$response_type,
        continuous = private$current_betaT,
        incidence = {
          p_t = clamp(stats::plogis(eta_t), private$incidence_clamp, 1 - private$incidence_clamp)
          p_c = clamp(stats::plogis(eta_c), private$incidence_clamp, 1 - private$incidence_clamp)
          mean(p_t - p_c)
        },
        proportion = {
          mu_t = clamp(stats::plogis(eta_t), private$proportion_clamp, 1 - private$proportion_clamp)
          mu_c = clamp(stats::plogis(eta_c), private$proportion_clamp, 1 - private$proportion_clamp)
          mean(mu_t - mu_c)
        },
        count = {
          mu_t = pmax(private$count_clamp, exp(eta_t))
          mu_c = pmax(private$count_clamp, exp(eta_c))
          mean(mu_t - mu_c)
        },
        survival = {
          shape_t = pmax(private$survival_clamp, exp(eta_t))
          shape_c = pmax(private$survival_clamp, exp(eta_c))
          mean_t = private$k_survival * gamma(1 + 1 / shape_t)
          mean_c = private$k_survival * gamma(1 + 1 / shape_c)
          (1 - private$prob_censoring / 2) * mean(mean_t - mean_c)
        },
        ordinal = private$compute_true_ordinal_mean_diff(eta_c, eta_t),
        private$current_betaT
      )
    },

    compute_true_ordinal_mean_diff = function(eta_c, eta_t) {
      expected_ordinal = function(eta) {
        K = private$n_ordinal_levels
        if (private$sd_noise <= 0) {
          rounded_eta = sign(eta) * floor(abs(eta) + 0.5)
          return(pmin(K, pmax(1, rounded_eta)))
        }

        sigma = private$sd_noise
        probs = matrix(0, nrow = length(eta), ncol = K)
        probs[, 1L] = stats::pnorm((1.5 - eta) / sigma)
        if (K > 2L) {
          for (k in 2L:(K - 1L)) {
            lo = (k - 0.5 - eta) / sigma
            hi = (k + 0.5 - eta) / sigma
            probs[, k] = stats::pnorm(hi) - stats::pnorm(lo)
          }
        }
        if (K > 1L) {
          probs[, K] = 1 - stats::pnorm((K - 0.5 - eta) / sigma)
        }
        as.numeric(probs %*% seq_len(K))
      }

      mean(expected_ordinal(eta_t) - expected_ordinal(eta_c))
    },

    # Instantiate design and run the full experiment (assign + observe all n).
    .build_design = function(design_gen, X, y_linear_model, design_extra) {
      n       = private$current_n

      # Auto-inject required args that depend on the covariate matrix when the
      # user has not already supplied them via design_classes_and_params.
      init_fn = private$.get_r6_init_fn(design_gen)
      if (!is.null(init_fn)) {
        fn_formals = names(formals(init_fn))
        x_names    = names(X)
        if ("strata_cols" %in% fn_formals && !"strata_cols" %in% names(design_extra))
          design_extra$strata_cols = x_names[1L]
        if ("cluster_col" %in% fn_formals && !"cluster_col" %in% names(design_extra))
          design_extra$cluster_col = x_names[min(2L, length(x_names))]
        if ("factors"     %in% fn_formals && !"factors"     %in% names(design_extra))
          design_extra$factors = list(treatment = 2L)
      }

      des_obj = do.call(design_gen$new, c(
        list(response_type = private$response_type, n = n),
        design_extra
      ))

      if (inherits(des_obj, "DesignSeqOneByOne")) {
        # Sequential: assignment depends on prior responses so w is obtained
        # one subject at a time.  Call Rcpp with length-1 vectors per subject.
        for (t in seq_len(n)) {
          w_t = des_obj$add_one_subject_to_experiment_and_assign(
            X[t, , drop = FALSE])
          out = apply_treatment_and_noise_cpp(
            y_linear_model[t], w_t,
	            private$response_type, private$current_betaT,
	            private$sd_noise, private$prob_censoring,
	            private$n_ordinal_levels,
	            phi_proportion = private$phi_proportion,
	            k_survival = private$k_survival,
	            incidence_clamp = private$incidence_clamp,
	            proportion_clamp = private$proportion_clamp,
	            count_clamp = private$count_clamp,
	            survival_clamp = private$survival_clamp)
          des_obj$add_one_subject_response(t, out$y, out$dead)
        }
      } else {
        # Fixed: all assignments known upfront — vectorize across all n subjects.
        des_obj$add_all_subjects_to_experiment(X)
        des_obj$assign_w_to_all_subjects()
        w   = des_obj$get_w()
        out = apply_treatment_and_noise_cpp(
          y_linear_model, w,
	          private$response_type, private$current_betaT,
	          private$sd_noise, private$prob_censoring,
	          private$n_ordinal_levels,
	          phi_proportion = private$phi_proportion,
	          k_survival = private$k_survival,
	          incidence_clamp = private$incidence_clamp,
	          proportion_clamp = private$proportion_clamp,
	          count_clamp = private$count_clamp,
	          survival_clamp = private$survival_clamp)
        des_obj$add_all_subject_responses(out$y, out$dead)
      }
      des_obj
    },

    # Append one row to raw_results as a plain list (no data.table allocation).
    # rbindlist() in get_results() converts the whole list at once.
    .record = function(rep, design, inf_name, inference_type, est, ci, pval, true_estimand) {
      ci2 = if (length(ci) >= 2L) as.numeric(ci[1:2]) else c(NA_real_, NA_real_)
      if (all(is.finite(ci2)) && ci2[1L] > ci2[2L]) ci2 = rev(ci2)
      
      private$results_idx = private$results_idx + 1L
      row = list(
        response_type = private$response_type,
        rep           = as.integer(rep),
        cond_exp_func_model = private$current_cond_exp_func_model,
        n             = as.integer(private$current_n),
        p             = as.integer(private$current_p),
        betaT         = as.numeric(private$current_betaT),
        design        = design,
        inference     = inf_name,
        inference_type = inference_type,
        estimate      = if (is.null(est) || !is.finite(est)) NA_real_ else as.numeric(est),
        ci_lo         = ci2[1L],
        ci_hi         = ci2[2L],
        pval          = if (is.null(pval) || length(pval) == 0L || !is.finite(pval[1L]))
                          NA_real_ else as.numeric(pval[1L]),
        true_estimand = as.numeric(true_estimand)
      )
      private$raw_results[[private$results_idx]] = row
      private$seen_result_keys = c(
        private$seen_result_keys,
        private$.result_key(rep, design, inf_name, inference_type)
      )
      row_dt = data.table::as.data.table(row)
      private$.append_result_row_to_file(row_dt)
      private$.advance_progress()
    },

    # ── Defaults ──────────────────────────────────────────────────────────────

    .default_design_classes = function() {
      list(
        # ── Fixed ──────────────────────────────────────────────────────────────
        FixedDesignBernoulli,
        FixedDesigniBCRD,
        FixedDesignBinaryMatch,
        FixedDesignBlocking,                # strata_cols auto-injected if absent
        FixedDesignGreedy,
        FixedDesignMatchingGreedyPairSwitching,
        FixedDesignRerandomization,
        FixedDesignOptimalBlocks,
        FixedDesignCluster,                 # cluster_col auto-injected if absent
        FixedDesignBlockedCluster,          # strata_cols + cluster_col auto-injected if absent
        FixedDesignDOptimal,
        FixedDesignAOptimal,
        FixedDesignFactorial,               # factors auto-injected if absent

        # ── Sequential one-by-one ──────────────────────────────────────────────
        DesignSeqOneByOneBernoulli,
        DesignSeqOneByOneiBCRD,
        DesignSeqOneByOneEfron,
        DesignSeqOneByOneAtkinson,
        DesignSeqOneByOneUrn,
        DesignSeqOneByOneRandomBlockSize,   # strata_cols auto-injected if absent
        DesignSeqOneByOneSPBR,              # strata_cols auto-injected if absent
        DesignSeqOneByOnePocockSimon,       # strata_cols auto-injected if absent
        DesignSeqOneByOneKK21,
        DesignSeqOneByOneKK21stepwise,
        DesignSeqOneByOneKK14
      )
    },

    .default_inference_classes = function() {
      rt   = private$response_type
      univ = if (!(rt == "survival" && private$prob_censoring > 0)) list(InferenceAllSimpleMeanDiff) else list()

      type_specific = switch(rt,
        continuous = list(
          InferenceAllSimpleWilcox,
          InferenceContinMultOLS,
          InferenceContinMultLin,
          InferenceContinUnivRobustRegr,
          InferenceContinMultiRobustRegr,
          InferenceContinMultOLSKKIVWC,
          InferenceContinMultOLSKKOneLik
        ),
        incidence = list(
          InferenceIncidUnivLogRegr,
          InferenceIncidMultiLogRegr,
          InferenceIncidUnivModifiedPoisson,
          InferenceIncidMultiModifiedPoisson,
          InferenceIncidUnivKKClogitIVWC,
          InferenceIncidMultiKKClogitOneLik,
          InferenceIncidCMH,
          InferenceIncidExtendedRobins,
          InferenceIncidExactZhang,
          InferenceIncidExactFisher,
          InferenceIncidenceExactBinomial
        ),
        proportion = list(
          InferenceAllSimpleWilcox,
          InferencePropUniBetaRegr,
          InferencePropMultiBetaRegr,
          InferencePropUniFractionalLogit,
          InferencePropMultiFractionalLogit,
          InferencePropUnivKKGEE,
          InferencePropMultiKKQuantileRegrIVWC
        ),
        count = list(
          InferenceAllSimpleWilcox,
          InferenceCountUnivPoissonRegr,
          InferenceCountMultiPoissonRegr,
          InferenceCountUnivRobustPoissonRegr,
          InferenceCountMultiRobustPoissonRegr,
          InferenceCountPoissonUnivKKGEE,
          InferenceCountPoissonUnivKKCPoissonIVWC
        ),
        survival = list(
          InferenceSurvivalUniCoxPHRegr,
          InferenceSurvivalMultiCoxPHRegr,
          InferenceSurvivalLogRank,
          InferenceSurvivalRestrictedMeanDiff,
          InferenceSurvivalUnivKKStratCoxIVWC,
          InferenceSurvivalUnivKKLWACoxIVWC
        ),
        ordinal = list(
          InferenceOrdinalUniPropOddsRegr,
          InferenceOrdinalUniCumulProbitRegr,
          InferenceOrdinalUniCLLRegr,
          InferenceOrdinalUnivKKGEE
        ),
        stop("Unknown response_type: ", rt)
      )

      c(univ, type_specific)
    }
  )
)
