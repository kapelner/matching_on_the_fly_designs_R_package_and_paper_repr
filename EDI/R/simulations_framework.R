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
#     InferenceContinOLS = list(max_resample_attempts = 25L),
#     InferenceContinKKOLSIVWC
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
#' @examples
#' generate_covariate_dataset(n = 10, p = 5)
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
    if (!is.matrix(X_mat)) X_mat = as.matrix(X_mat)
  } else {
    X_mat = matrix(do.call(cov_draw_method, c(list(n * p), cov_draw_method_args)), nrow = n, ncol = p)
  }
  
  if (is.null(colnames(X_mat))) {
    colnames(X_mat) = paste0("x", seq_len(p))
  }
  X = data.table::as.data.table(X_mat)
  
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
#' @examples
#' transform_cont_y_based_on_response_type(rnorm(10), 'incidence')
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
      include.lowest = TRUE, labels = FALSE
    )),
    stop("Unknown response_type: ", response_type)
  )
}

# Internal helper for null-coalescing
`%||%` = function(a, b) if (!is.null(a)) a else b

# Walk an R6 class generator's inheritance chain to find the first class
# that defines initialize(), and return that function.
get_r6_init_fn = function(r6gen) {
  gen = r6gen
  while (!is.null(gen)) {
    fn = tryCatch(gen$public_methods$initialize, error = function(e) NULL)
    if (!is.null(fn)) return(fn)
    gen = tryCatch(gen$get_inherit(), error = function(e) NULL)
  }
  NULL
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
#' @examples
#' \donttest{
#' # Simple simulation with two designs and two inference methods
#' sim = SimulationFramework$new(
#'   response_type = "continuous",
#'   design_classes_and_params = list(
#'     DesignSeqOneByOneKK21 = list(lambda = 0.5),
#'     DesignSeqOneByOneBernoulli = list()
#'   ),
#'   inference_classes_and_params = list(
#'     InferenceContinOLS = list(),
#'     InferenceContinKKOLSIVWC = list()
#'   ),
#'   n = 100, p = 5, Nrep = 10, betaT = 1
#' )
#' sim$run()
#' sim$summarize()
#' }
#' @export
SimulationFramework = R6::R6Class("SimulationFramework",
  lock_objects = FALSE,

  # ── public ─────────────────────────────────────────────────────────────────
  public = list(

    #' @description
    #' Create a new \code{SimulationFramework}.
    #'
    #' @param response_type \strong{(required)} Character scalar or vector.  The type of
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
    #'   \code{list(InferenceContinOLS, InferenceContinKKOLSIVWC)}.
    #'   Named entries use the entry name as the inference class and the value
    #'   as the constructor parameter list, for example
    #'   \code{list(InferenceContinOLS = list(max_resample_attempts = 25L))}.
    #'   Duplicate named entries are allowed for repeated inference classes with
    #'   different parameters. Supplied parameters must be accepted by the
    #'   inference class constructor.
    #'   \code{NULL} selects a curated set for the given
    #'   \code{response_type}: several universal classes that work with any
    #'   design, plus representative KK-specific classes (silently skipped for
    #'   non-KK designs at runtime).
    #'
    #' @param n Integer scalar or vector.  Sample size per simulation replication.  Default
    #'   \code{100}.
    #'
    #' @param p Integer scalar or vector.  Number of covariates.  Must be \eqn{\ge 5} when
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
    #' @param betaT Numeric scalar or vector.  True treatment effect added to treated
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
    #'   \code{response_type = "survival"}.  Default \code{0.25}.
    #'
    #' @param verbose Logical.  If \code{TRUE}, prints a message for every
    #'   replication and for every \code{(design, inference)} pair that is
    #'   skipped due to an error.  Default \code{TRUE}.
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
    #' @param num_cores Positive integer.  Number of worker processes for
    #'   parallel execution of Monte Carlo replications.  Note that when
    #'   \code{num_cores > 1}, parallelization *within* individual inference
    #'   routines (e.g. bootstrap, randomization) is automatically disabled
    #'   to prevent thread oversubscription.  Default \code{1}.
    #'
    #' @param results_filename Character scalar. The filename for the results
    #'   file. Supported extensions are \code{.csv} and \code{.csv.bz2}.
    #'   Default \code{"simulation_framework_results.csv.bz2"}.
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
      num_cores             = 1L,
      seed                  = NULL,
      cov_draw_method       = stats::rnorm,
      cov_draw_method_args  = list(mean = 0, sd = 1),
      random_X_draws        = TRUE,
      prob_censoring        = 0.25,
      verbose                    = TRUE,
      keep_all_intermediate_data = FALSE,
      turn_off_asserts_for_speed = TRUE,
      inference_types_and_params = NULL,
      results_filename      = "simulation_framework_results.csv.bz2",
      continue_from_last_result_row = TRUE
    ) {
      valid_rt = c("continuous", "incidence", "proportion",
                   "count", "survival", "ordinal")
      if (any(!response_type %in% valid_rt))
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
      if (!is.character(results_filename) || length(results_filename) != 1L || is.na(results_filename))
        stop("results_filename must be a single non-missing character string")
      results_file_format = private$.results_file_format(results_filename)
      if (is.na(results_file_format)) {
        stop("results_filename must end in either '.csv' or '.csv.bz2'")
      }
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

      private$response_type_values = unique(as.character(response_type))
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
      private$num_cores            = as.integer(num_cores)
      private$seed                 = if (is.null(seed)) NULL else as.integer(seed)
      private$cov_draw_method      = if (is.null(X_mat)) cov_draw_method else NULL
      private$cov_draw_method_args = cov_draw_method_args
      private$random_X_draws       = random_X_draws
      private$prob_censoring       = prob_censoring
      private$verbose                    = verbose
      private$turn_off_asserts_for_speed = turn_off_asserts_for_speed
      
      if (as.integer(num_cores) > 1L && keep_all_intermediate_data) {
        stop("Multithreading (num_cores > 1) is incompatible with 'keep_all_intermediate_data = TRUE'.")
      }
      
      private$keep_all_intermediate_data = keep_all_intermediate_data
      private$results_filename     = results_filename
      private$continue_from_last_result_row = continue_from_last_result_row
      private$inf_types        = inf_types
      private$inference_type_params = inf_type_spec
      private$param_grid       = private$.build_param_grid(
        n_values,
        p_values,
        betaT_values,
        cond_exp_func_model_values,
        private$response_type_values
      )
      private$current_n        = private$param_grid$n[[1L]]
      private$current_p        = private$param_grid$p[[1L]]
      private$current_betaT    = private$param_grid$betaT[[1L]]
      private$current_cond_exp_func_model = private$param_grid$cond_exp_func_model[[1L]]
      private$current_response_type = private$param_grid$response_type[[1L]]

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
        on.exit(toggle_asserts(TRUE), add = TRUE)        
      }
      
      # ── Parallelism management ─────────────────────────────────────────────
      # Save state to restore on exit
      ns = asNamespace("EDI")
      prev_threads = getOption(".edi_last_set_threads")
      if (is.null(prev_threads)) prev_threads = 1L
      prev_global_cores = get_num_cores()
      prev_override = ns$edi_env$num_cores_override
      
      num_cores = private$num_cores
      set_package_threads(num_cores)
      
      # If the simulation is running its own parallel loop, we force 
      # nested core budget to 1 to prevent N*M explosion.
      if (num_cores > 1L && prev_global_cores > 1L) {
        set_num_cores(1L)
      }
      
      on.exit({
        set_package_threads(prev_threads)
        if (get_num_cores() != prev_global_cores) {
           set_num_cores(prev_global_cores)
        }
        assign("num_cores_override", prev_override, envir = ns$edi_env)
      }, add = TRUE)

      # Handle cleanup/setup
      if (!isTRUE(private$continue_from_last_result_row)) {
        if (file.exists(private$results_filename)) unlink(private$results_filename)
        private$.cleanup_results_staging_file()
      }
      
      private$simulation_start_time = as.numeric(Sys.time())

      n_des = length(private$design_classes)
      n_inf = length(private$inference_classes)
      n_met = length(private$inf_types)
      n_cells = nrow(private$param_grid)

      # Ensure staging file is ready if we are using bz2 and continuing
      if (isTRUE(private$continue_from_last_result_row)) {
        private$.ensure_staging_file_exists()
      }

      existing_results = private$.load_existing_results()
      n_existing = nrow(existing_results)
      
      if (isTRUE(private$verbose) && n_existing > 0L) {
        private$.message_stderr(sprintf("%d existing results loaded\n", n_existing))
      }

      # Pre-allocate results list for data.table batches (one per chunk or serial replication)
      private$raw_results = vector("list", 1L + private$Nrep * n_cells)
      private$results_idx = 0L

      # Use C++ ResultKeyStore for O(1) key lookups without R string interning overhead
      init_result_key_store_cpp(n_existing + 1000L) # Reserve extra space for new results
      if (n_existing > 0L) {
        # Optimization: Store the entire data.table as the first element instead of row-by-row lists
        private$raw_results[[1L]] = existing_results
        private$results_idx = 1L
        
        # Extract columns once as vectors (fast pointer-based extraction in data.table)
        rt_v = existing_results$response_type
        cm_v = existing_results$cond_exp_func_model
        n_v  = existing_results$n
        p_v  = existing_results$p
        bt_v = existing_results$betaT
        rp_v = existing_results$rep
        ds_v = existing_results$design
        if_v = existing_results$inference
        it_v = existing_results$inference_type

        # Chunking allows showing progress and avoids one massive allocation.
        chunk_size_indexing = 200000L
        indices = seq(1L, n_existing, by = chunk_size_indexing)
        
        for (i in seq_along(indices)) {
          if (isTRUE(private$verbose)) private$.draw_labeled_progress_bar("indexing existing results", (i - 1L) / length(indices))
          start = indices[i]
          end = min(n_existing, start + chunk_size_indexing - 1L)
          
          add_to_result_key_store_cpp(
            rt_v, cm_v, n_v, p_v, bt_v, rp_v, ds_v, if_v, it_v,
            start, end
          )
        }
        if (isTRUE(private$verbose)) {
          private$.draw_labeled_progress_bar("indexing existing results", 1)
          cat("\n", file = stderr())
        }
      }

      private$all_intermediate_data = vector("list", private$Nrep * n_cells)
      private$valid_combos          = list()
      private$seen_combo_keys      = character(0L)
      private$exact_warned_classes = character(0L)

      if (isTRUE(private$verbose)) {
        message(sprintf(
          "simulations: CEF_mod=%s  n=%s  p=%s  Nrep=%d  betaT=%s designs=%d inferences=%d num_cores=%d",
          private$.format_values(private$cond_exp_func_model_values),
          private$.format_values(private$n_values),
          private$.format_values(private$p_values),
          private$Nrep,
          private$.format_values(private$betaT_values), 
          n_des, 
          n_inf,
          num_cores
        ))
      }

      log_interval = max(1L, private$Nrep %/% 10L)
      shared_X_draws = list()
      if (!isTRUE(private$random_X_draws)) {
        np_grid = unique(private$param_grid[, .(n, p)])
        n_np = nrow(np_grid)
        
        # Optimization: convert X_mat once outside the loop
        X_mat_matrix = if (!is.null(private$X_mat)) as.matrix(private$X_mat) else NULL
        
        for (np_idx in seq_len(n_np)) {
          if (isTRUE(private$verbose)) private$.draw_labeled_progress_bar("generating covariate matrices", (np_idx - 1L) / n_np)
          n_i = np_grid$n[[np_idx]]
          p_i = np_grid$p[[np_idx]]
          X_i = if (is.null(X_mat_matrix)) {
            m = matrix(
              do.call(private$cov_draw_method, c(list(n_i * p_i), private$cov_draw_method_args)),
              nrow = n_i,
              ncol = p_i
            )
            colnames(m) = paste0("x", seq_len(p_i))
            m
          } else {
            X_mat_matrix # Names should already be there or handled by generate_covariate_dataset
          }
          shared_X_draws[[paste(n_i, p_i, sep = "|")]] = X_i
        }
        if (isTRUE(private$verbose)) {
          private$.draw_labeled_progress_bar("generating covariate matrices", 1)
          cat("\n", file = stderr())
        }
      }

      # Pre-calculate planned combos and cache cell metadata
      had_seed_plan = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      if (had_seed_plan) old_seed_plan = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      
      planned_combos_list = vector("list", n_cells)
      cell_reps_to_run    = vector("list", n_cells)
      validity_cache      = list()
      cell_tasks_count    = integer(n_cells)
      total_existing_planned_tasks = 0L
      
      for (cell_idx in seq_len(n_cells)) {
        n_i  = private$param_grid$n[[cell_idx]]
        p_i  = private$param_grid$p[[cell_idx]]
        dt_i = private$param_grid$cond_exp_func_model[[cell_idx]]
        bt_i = private$param_grid$betaT[[cell_idx]]
        rt_i = private$param_grid$response_type[[cell_idx]]

        # Show more granular progress during plan pre-calculation
        if (isTRUE(private$verbose)) {
          label = sprintf("pre-calculating plan (cell %d/%d: n=%d, p=%d, %s, %s)",
                          cell_idx, n_cells, n_i, p_i, dt_i, rt_i)
          private$.draw_labeled_progress_bar(label, (cell_idx - 1) / n_cells)
        }

        cache_key = paste(rt_i, n_i, p_i, dt_i, sep = "|")
        if (is.null(validity_cache[[cache_key]])) {
           private$current_n = n_i; private$current_p = p_i; private$current_betaT = bt_i; private$current_cond_exp_func_model = dt_i; private$current_response_type = rt_i
           rep_data = if (isTRUE(private$random_X_draws)) private$.generate_data() else private$.generate_data_from_X(shared_X_draws[[paste(n_i, p_i, sep = "|")]])
           validity_cache[[cache_key]] = private$.build_valid_combos_for_current_cell(rep_data)
        }
        
        cell_combos = lapply(validity_cache[[cache_key]], function(c) {
           c$betaT = bt_i; c$n = n_i; c$p = p_i; c$cond_exp_func_model = dt_i; c$response_type = rt_i; c
        })
        planned_combos_list[[cell_idx]] = cell_combos
        n_combos = length(cell_combos)
        cell_tasks_count[cell_idx] = n_combos
        
        # Identify missing replications for this specific cell
        reps_to_run = seq_len(private$Nrep)
        if (n_existing > 0L && n_combos > 0L) {
          # Use batches for safety with extremely large Nrep
          rep_batch_size = max(100L, 500000L %/% n_combos)
          rep_batches = split(seq_len(private$Nrep), (seq_len(private$Nrep) - 1) %/% rep_batch_size)
          
          c_des  = vapply(cell_combos, `[[`, "", "design")
          c_inf  = vapply(cell_combos, `[[`, "", "inference")
          c_type = vapply(cell_combos, `[[`, "", "inference_type")
          
          cell_req = rep(TRUE, private$Nrep)
          for (b_reps in rep_batches) {
            nb = length(b_reps)
            all_exists = check_in_result_key_store_cpp(
              rep(rt_i, nb * n_combos),
              rep(dt_i, nb * n_combos),
              rep(as.integer(n_i), nb * n_combos),
              rep(as.integer(p_i), nb * n_combos),
              rep(as.numeric(bt_i), nb * n_combos),
              rep(as.integer(b_reps), each = n_combos),
              rep(c_des, nb),
              rep(c_inf, nb),
              rep(c_type, nb)
            )
            total_existing_planned_tasks = total_existing_planned_tasks + sum(all_exists)
            exists_mat = matrix(all_exists, nrow = n_combos)
            cell_req[b_reps] = colSums(exists_mat) < n_combos
          }
          if (!private$keep_all_intermediate_data) {
             reps_to_run = which(cell_req)
          }
        }
        cell_reps_to_run[[cell_idx]] = reps_to_run
      }
      if (isTRUE(private$verbose)) {
        private$.draw_labeled_progress_bar("pre-calculating simulation plan", 1)
        cat("\n", file = stderr())
      }
      if (had_seed_plan) assign(".Random.seed", old_seed_plan, envir = .GlobalEnv)
      
      planned_combos = unlist(planned_combos_list, recursive = FALSE)
      private$valid_combos = planned_combos
      private$progress_total = length(planned_combos) * private$Nrep
      
      private$progress_count = total_existing_planned_tasks
      private$initial_progress_count = total_existing_planned_tasks
      private$progress_log_interval = max(1L, private$progress_total %/% 20L)
      private$progress_bar = NULL; private$use_progress_bar = FALSE
      private$total_cells = n_cells; private$last_progress_draw_time = 0
      private$current_task_label = "Des/Inf"
      
      # For ETA calculation, we need to know where we started in the cell hierarchy
      private$initial_cell_in_progress_prop = if (private$progress_total > 0) private$progress_count / private$progress_total else 0
      
      if (isTRUE(private$verbose)) {
        private$.print_plan_summary(planned_combos_list)
        private$use_progress_bar = TRUE

        # Precisely set initial progress bar state based on first missing task
        first_active_cell = which(vapply(cell_reps_to_run, length, 0L) > 0L)[1L]
        if (!is.na(first_active_cell)) {
          private$current_cell_idx = first_active_cell
          private$current_rep_idx  = cell_reps_to_run[[first_active_cell]][1L]
          private$tasks_per_rep = cell_tasks_count[first_active_cell]
          private$current_task_in_rep_idx = 0L
        } else {
          # Everything already done
          private$current_cell_idx = n_cells
          private$current_rep_idx  = private$Nrep
          private$tasks_per_rep = if (length(cell_tasks_count) > 0L) cell_tasks_count[[length(cell_tasks_count)]] else 0L
          private$current_task_in_rep_idx = private$tasks_per_rep
        }

        # If we are going to run and use the progress bar, ensure we clean up the 
        # multi-line block at the end of the run().
        on.exit({
           if (!is.null(private$progress_bar_drawn)) {
              cat("\n", file = stderr())
              private$progress_bar_drawn = NULL
           }
        }, add = TRUE)

        private$.draw_progress()
      }

      # Early exit if everything is already done
      if (private$progress_count >= private$progress_total) {
        if (isTRUE(private$verbose)) message("Simulation already complete.")
        private$.finish_run()
        return(invisible(self))
      }

      for (cell_idx in seq_len(n_cells)) {
        private$current_cell_idx = cell_idx
        private$current_n = private$param_grid$n[[cell_idx]]
        private$current_p = private$param_grid$p[[cell_idx]]
        private$current_betaT = private$param_grid$betaT[[cell_idx]]
        private$current_cond_exp_func_model = private$param_grid$cond_exp_func_model[[cell_idx]]
        private$current_response_type = private$param_grid$response_type[[cell_idx]]
        private$tasks_per_rep = cell_tasks_count[cell_idx]
        
        reps_to_run = cell_reps_to_run[[cell_idx]]
        run_required_v = rep(FALSE, private$Nrep)
        run_required_v[reps_to_run] = TRUE

        # Pre-format key prefix for this cell to avoid redundant formatting in workers
        cell_key_prefix = paste(private$current_response_type, private$current_cond_exp_func_model,
                                private$current_n, private$current_p,
                                private$current_betaT, sep = "|")

        cell_state = list(
          n = private$current_n, p = private$current_p, betaT = private$current_betaT,
          cond_exp_func_model = private$current_cond_exp_func_model, norm_sq_beta_vec = private$norm_sq_beta_vec,
          response_type = private$current_response_type, random_X_draws = private$random_X_draws,
          shared_X = if (!isTRUE(private$random_X_draws)) shared_X_draws[[paste(private$current_n, private$current_p, sep = "|")]] else NULL,
          X_mat = private$X_mat, cov_draw_method = private$cov_draw_method, cov_draw_method_args = private$cov_draw_method_args,
          design_classes = private$design_classes, design_labels = private$design_labels, design_params = private$design_params,
          inference_classes = private$inference_classes, inference_labels = private$inference_labels, inference_ctor_params = private$inference_constructor_params,
          inf_types = private$inf_types, alpha = private$alpha, B_boot = private$B_boot, r_rand = private$r_rand,
          pval_epsilon = private$pval_epsilon, sd_noise = private$sd_noise, n_ordinal_levels = private$n_ordinal_levels,
          phi_proportion = private$phi_proportion, k_survival = private$k_survival,
          incidence_clamp = private$incidence_clamp, proportion_clamp = private$proportion_clamp, count_clamp = private$count_clamp,
          survival_clamp = private$survival_clamp, survival_min_time = private$survival_min_time, count_min_rate = private$count_min_rate,
          count_shift = private$count_shift, prob_censoring = private$prob_censoring,
          inference_type_params = private$inference_type_params,
          cell_key_prefix = cell_key_prefix
        )

        num_cores_to_use = 1L
        private$current_task_label = "Des/Inf"
        if (num_cores > 1L) {
          has_java_designs = any(vapply(private$design_classes, function(cls) cls$classname %in% c("FixedDesignBinaryMatch", "FixedDesignGreedy", "FixedDesignMatchingGreedyPairSwitching", "FixedDesignRerandomization"), logical(1L)))
          if (has_java_designs) {
            if (isTRUE(private$verbose)) private$.message_stderr("Warning: Multithreading disabled for Java-based designs; forking is unstable. Serial execution used for this cell.\n")
          } else {
            num_cores_to_use = num_cores
            chunk_size = num_cores_to_use
            rep_chunks = split(seq_len(private$Nrep), (seq_len(private$Nrep) - 1) %/% chunk_size)
            for (chunk in rep_chunks) {
               # Identify which reps in this chunk actually need running
               active_reps = chunk[run_required_v[chunk]]
               
               if (length(active_reps) == 0L) {
                  private$current_rep_idx = chunk[[length(chunk)]]
                  private$tasks_per_rep = length(chunk)
                  private$current_task_in_rep_idx = length(chunk)
                  private$.draw_progress()
                  next
               }

               private$current_rep_idx = chunk[[1L]]
               private$tasks_per_rep = length(active_reps)
               private$current_task_in_rep_idx = 0L
               private$last_progress_draw_time = 0
               private$.draw_progress()
               RUN_REP_DETACHED = private$.run_single_replication_in_worker; environment(RUN_REP_DETACHED) = asNamespace("EDI")
               RUN_REP = function(rep_i) RUN_REP_DETACHED(rep_i, cell_state, is_forked = TRUE)
               jobs = lapply(active_reps, function(rep_i) {
                 parallel::mcparallel(
                   RUN_REP(rep_i),
                   name = as.character(rep_i),
                   mc.set.seed = TRUE,
                   silent = TRUE
                 )
               })
               names(jobs) = as.character(active_reps)
               chunk_results = vector("list", length(active_reps))
               names(chunk_results) = as.character(active_reps)
               completed_in_chunk = 0L
               while (completed_in_chunk < length(active_reps)) {
                 collected = parallel::mccollect(jobs, wait = FALSE, timeout = 0.1)
                 if (is.null(collected) || length(collected) == 0L) next
                 for (nm in names(collected)) {
                   if (is.null(chunk_results[[nm]])) {
                     chunk_results[[nm]] = collected[[nm]]
                     completed_in_chunk = completed_in_chunk + 1L
                   }
                 }
                 private$current_rep_idx = min(chunk[[length(chunk)]], chunk[[1L]] - 1L + completed_in_chunk)
                 private$current_task_in_rep_idx = completed_in_chunk
                 private$.draw_progress()
               }
               chunk_results = unname(chunk_results[as.character(active_reps)])
               chunk_results_ok = chunk_results[!vapply(chunk_results, function(x) is(x, "try-error") || is.null(x), logical(1L))]
               if (length(chunk_results_ok) > 0L) {
                 all_chunk_dt = data.table::rbindlist(lapply(chunk_results_ok, function(w) w$results_dt), use.names = TRUE, fill = TRUE)
                 all_chunk_skipped = sum(vapply(chunk_results_ok, function(w) w$skipped_count, 0L))
                 private$current_rep_idx = chunk[[length(chunk)]]
                 private$current_task_in_rep_idx = length(chunk)
                 private$.record_batch(all_chunk_dt, all_chunk_skipped)
               }
               failed_indices = which(vapply(chunk_results, function(x) is(x, "try-error") || is.null(x), logical(1L)))
               if (length(failed_indices) > 0L) {
                 # For failures, we only advance progress by the number of tasks
                 # that were actually pending for these replications (to avoid double-counting existing ones).
                 tasks_per_rep = cell_tasks_count[cell_idx]
                 tasks_to_add = sum(cell_req[failed_indices, , drop = FALSE])
                 private$progress_count = private$progress_count + tasks_to_add
                 private$.draw_progress()
               }
            }
          }
        }
        
        if (num_cores_to_use == 1L) {
          private$current_task_label = "Des/Inf"
          private$tasks_per_rep = cell_tasks_count[cell_idx]
          for (rep in seq_len(private$Nrep)) {
            if (!run_required_v[rep]) {
               # Fast-forward progress bar for fully skipped replications
               if (rep %% 100L == 0L || rep == private$Nrep) {
                  private$current_rep_idx = rep; private$current_task_in_rep_idx = 0L; private$.draw_progress()
               }
               next
            }
            private$current_rep_idx = rep; private$current_task_in_rep_idx = 0L
            private$last_progress_draw_time = 0
            private$.draw_progress()
            worker_out = private$.run_single_replication_in_worker(rep, cell_state, progress_cb = private$.advance_progress, is_forked = FALSE)
            if (private$keep_all_intermediate_data) {
               rep_data = if (isTRUE(private$random_X_draws)) private$.generate_data() else private$.generate_data_from_X(shared_X_draws[[paste(private$current_n, private$current_p, sep = "|")]])
               X = rep_data$X; y_linear_model = rep_data$y_linear_model; true_mean_diff_ate = private$compute_true_mean_diff_ate(y_linear_model); rep_slot = (cell_idx - 1L) * private$Nrep + rep
               for (di in seq_along(private$design_classes)) {
                 design_gen   = private$design_classes[[di]]; design_name  = private$design_labels[[di]]; design_extra = if (!is.null(private$design_params)) private$design_params[[di]] else list()
                 des_obj = private$.build_design(design_gen, X, y_linear_model, design_extra); private$all_intermediate_data[[rep_slot]]$designs[[design_name]] = des_obj
                 if (is.null(des_obj)) next
                 for (ii in seq_along(private$inference_classes)) {
                   inf_gen  = private$inference_classes[[ii]]; inf_name = private$inference_labels[[ii]]; inf_ctor_extra = private$inference_constructor_params[[ii]]
                   inf_obj = do.call(inf_gen$new, c(list(des_obj), inf_ctor_extra)); private$all_intermediate_data[[rep_slot]]$inferences[[design_name]][[inf_name]] = inf_obj
                 }
               }
               private$all_intermediate_data[[rep_slot]]$y_linear_model = y_linear_model
            }
            if (!is(worker_out, "try-error") && !is.null(worker_out)) {
               n_res = if (is.null(worker_out$results_dt)) 0L else nrow(worker_out$results_dt)
               private$current_task_in_rep_idx = worker_out$skipped_count + n_res
               private$.record_batch(worker_out$results_dt, worker_out$skipped_count)
            } else {
               # For failures, we only advance progress by the number of tasks
               # that were actually pending for this replication.
               tasks_to_add = sum(cell_req[rep, ])
               private$current_task_in_rep_idx = private$tasks_per_rep; private$progress_count = private$progress_count + tasks_to_add; private$.draw_progress()
            }
          }
        }
      }

      private$.finish_run()
      invisible(self)
    },
    #' @description
    #' Retrieve the stored intermediate data (design and inference objects)
    #' for every replication. Only available if \code{keep_all_intermediate_data = TRUE}
    #' was passed to the constructor.
    #'
    #' @return A nested list containing the intermediate data for each replication,
    #'   or \code{NULL} if not recorded.
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
          ci_lo = numeric(), ci_hi = numeric(), pval = numeric(),
          true_estimand = numeric()
        ))
      # Prune pre-allocated list to only include what was actually recorded
      data.table::rbindlist(private$raw_results[seq_len(private$results_idx)], use.names = TRUE, fill = TRUE)
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
      ref_grid = data.table::rbindlist(lapply(private$valid_combos, as.list), use.names = TRUE, fill = TRUE)
      
      # Collapse over betaT: remove it from the unique keys
      if ("betaT" %in% names(ref_grid)) ref_grid[, betaT := NULL]
      ref_grid = unique(ref_grid)

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
      
      # Grouping columns - EXCLUDE betaT to allow averaging across effect sizes
      by_cols = c("response_type", "cond_exp_func_model", "n", "p", "design", "inference", "inference_type")

      if (nrow(dt) > 0L) {
        # Optimization: Avoid .SD subsetting within groups. Use vectorized logic.
        agg = dt[, {
          est_fin = is.finite(estimate)
          m_row = list(
            MSE   = if (any(est_fin)) mean((estimate[est_fin] - true_estimand[est_fin])^2) else NA_real_,
            n_est = sum(est_fin)
          )
          if (report_cov) {
            ci_fin = is.finite(ci_lo) & is.finite(ci_hi)
            m_row$coverage  = if (any(ci_fin)) mean(ci_lo[ci_fin] <= true_estimand[ci_fin] & true_estimand[ci_fin] <= ci_hi[ci_fin]) else NA_real_
            m_row$n_cov     = sum(ci_fin)
            m_row$ci_length = if (any(ci_fin)) mean(ci_hi[ci_fin] - ci_lo[ci_fin]) else NA_real_
          }
          if (report_pow) {
            pv_fin = is.finite(pval)
            is_zero = (abs(betaT) < 1e-12) # Use tolerance for safety
            pv_fin_zero = pv_fin & is_zero
            pv_fin_nonzero = pv_fin & !is_zero

            m_row$power = if (any(pv_fin_nonzero)) mean(pval[pv_fin_nonzero] < alpha) else NA_real_
            m_row$n_pow = sum(pv_fin_nonzero)
            
            if (any(is_zero)) {
              m_row$size = if (any(pv_fin_zero)) mean(pval[pv_fin_zero] < alpha) else NA_real_
              m_row$n_size = sum(pv_fin_zero)
            }
          }
          m_row
        }, by = by_cols]
      } else {
        agg = data.table::data.table(response_type = character(), cond_exp_func_model = character(), n = integer(), p = integer(),
                                     design = character(), inference = character(),
                                     inference_type = character(), power = numeric(), MSE = numeric(),
                                     n_est = integer(), n_pow = integer())
      }

      # ── Right-join: every valid combo appears, NA for those with no data ──────
      data.table::setkeyv(agg, by_cols)
      data.table::setkeyv(ref_grid, by_cols)
      result = agg[ref_grid]
      
      # Ensure n_est, n_pow, n_size are present and replace NA with 0
      if (!"n_est" %in% names(result)) result[, n_est := 0L]
      if (!"n_pow" %in% names(result)) result[, n_pow := 0L]
      if (!"n_size" %in% names(result)) result[, n_size := 0L]
      result[is.na(n_est), n_est := 0L]
      result[is.na(n_pow), n_pow := 0L]
      result[is.na(n_size), n_size := 0L]
      
      result[order(cond_exp_func_model, n, p, design, inference, inference_type)]
    },

    # ── print() ───────────────────────────────────────────────────────────────
    #' @description
    #' Print a summary of the \code{SimulationFramework} configuration and status.
    print = function() {
      cat("SimulationFramework\n")
      cat("  response_type :", private$.format_values(private$response_type_values), "\n")
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
      design_names = gsub("Design", "", design_names)
      inf_names    = gsub("Inference", "", private$inference_labels)
      inf_names    = gsub("^(Contin|Count|Incid|Prop|Survival|Ordinal|All)", "", inf_names)
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
    .finish_run = function() {
      private$has_run = TRUE
      if (private$.results_file_format(private$results_filename) == "csv.bz2" &&
          file.exists(private$.results_staging_filename())) {
        private$.sync_results_bz2_from_staging()
      }
      private$.cleanup_results_staging_file()
      clear_result_key_store_cpp() # Free memory
      if (isTRUE(private$use_progress_bar)) {
        private$current_cell_idx = private$total_cells
        private$current_rep_idx = private$Nrep
        private$current_task_in_rep_idx = private$tasks_per_rep
        private$.draw_simulation_progress_bars()
        cat("\n", file = stderr())
        private$progress_bar_drawn = NULL
      }
    },
    response_type_values = NULL,
    current_response_type = NULL,
    n_values         = NULL,
    p_values         = NULL,
    cond_exp_func_model_values = NULL,
    Nrep             = NULL,
    betaT_values     = NULL,
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
    num_cores            = NULL,
    seed                 = NULL,
    cov_draw_method      = NULL,
    cov_draw_method_args = NULL,
    random_X_draws       = NULL,
    prob_censoring       = NULL,
    verbose          = NULL,
    results_filename = NULL,
    simulation_start_time = NULL,
    initial_progress_count = NULL,
    initial_cell_in_progress_prop = NULL,
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
    all_intermediate_data     = NULL,    keep_all_intermediate_data = FALSE,
    has_run                   = FALSE,
    exact_warned_classes      = NULL,
    valid_combos              = NULL,
    seen_combo_keys  = NULL,
    seen_result_keys = NULL,
    total_cells              = 0L,
    current_cell_idx         = 0L,
    last_progress_draw_time  = 0,
    current_task_label       = "Des/Inf",
    current_rep_idx          = 0L,
    current_task_in_rep_idx  = 0L,
    tasks_per_rep            = 0L,
    progress_total           = 0L,
    progress_count           = 0L,
    progress_bar             = NULL,
    use_progress_bar         = FALSE,
    progress_log_interval    = 0L,
    progress_bar_drawn       = NULL,

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
      init_fn = get_r6_init_fn(r6gen)
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

    # Filter an arg list to only those accepted by r6gen's initialize().
    # If initialize accepts '...' or cannot be found, all args pass through.
    .filter_by_formals = function(r6gen, args) {
      if (length(args) == 0L) return(args)
      init_fn = get_r6_init_fn(r6gen)
      if (is.null(init_fn)) return(args)
      fn_formals = names(formals(init_fn))
      if ("..." %in% fn_formals) return(args)
      args[names(args) %in% fn_formals]
    },

    .has_private_method_on_object = function(obj, method_name) {
      exists(method_name, envir = obj$.__enclos_env__$private, inherits = FALSE)
    },

    # Serialize a named params list to "k=v, ..." string.
    .params_to_str = function(p) {
      if (is.null(p) || length(p) == 0L) return("")
      kv = mapply(function(k, v) paste0(k, "=", paste(deparse(v), collapse = "")),
                  names(p), p, SIMPLIFY = TRUE)
      paste(kv, collapse = ", ")
    },

    .build_param_grid = function(n_values, p_values, betaT_values, cond_exp_func_model_values, response_type_values) {
      grid = data.table::as.data.table(expand.grid(
        response_type = response_type_values,
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
        return(private$.load_existing_results_from_staging_or_empty(empty_dt))

      if (isTRUE(private$verbose)) private$.draw_labeled_progress_bar("loading previous results", 0)
      dt = private$.read_results_file(show_progress = FALSE)
      if (isTRUE(private$verbose)) {
        private$.draw_labeled_progress_bar("loading previous results", 1)
        cat("\n", file = stderr())
      }
      if (!"response_type" %in% names(dt))
        dt[, response_type := private$response_type_values[[1L]]]
      dt = dt[response_type %in% private$response_type_values]

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

    .load_existing_results_from_staging_or_empty = function(empty_dt) {
      if (private$.results_file_format(private$results_filename) != "csv.bz2")
        return(empty_dt)
      staging_filename = private$.results_staging_filename()
      if (!file.exists(staging_filename))
        return(empty_dt)
      data.table::fread(staging_filename, showProgress = FALSE)
    },

    .results_file_format = function(filename) {
      if (grepl("\\.csv\\.bz2$", filename, ignore.case = TRUE)) return("csv.bz2")
      if (grepl("\\.csv$", filename, ignore.case = TRUE)) return("csv")
      NA_character_
    },

    .results_staging_filename = function() {
      path = private$.results_output_path()
      file.path(
        dirname(path),
        paste0(
          sub("\\.csv\\.bz2$", "", basename(path), ignore.case = TRUE),
          "__staging.csv"
        )
      )
    },
    .results_output_path = function() {
      if (grepl("^/", private$results_filename)) {
        private$results_filename
      } else {
        file.path(getwd(), private$results_filename)
      }
    },

    .copy_binary_stream = function(from, to, chunk_size = 1024L * 1024L) {
      repeat {
        bytes = readBin(from, what = "raw", n = chunk_size)
        if (length(bytes) == 0L) break
        writeBin(bytes, to)
      }
      invisible(NULL)
    },

    .read_results_file = function(show_progress = FALSE) {
      format = private$.results_file_format(private$results_filename)
      if (identical(format, "csv")) {
        return(data.table::fread(private$results_filename, showProgress = show_progress))
      }
      if (!identical(format, "csv.bz2")) {
        stop("Unsupported results file format: ", private$results_filename)
      }

      staging_filename = private$.results_staging_filename()
      if (file.exists(staging_filename)) {
        return(data.table::fread(staging_filename, showProgress = show_progress))
      }

      extracted_file = tempfile("simulation_framework_results_", fileext = ".csv")
      input_con = bzfile(private$.results_output_path(), open = "rb")
      output_con = file(extracted_file, open = "wb")
      on.exit(try(close(input_con), silent = TRUE), add = TRUE)
      on.exit(try(close(output_con), silent = TRUE), add = TRUE)
      on.exit(unlink(extracted_file, force = TRUE), add = TRUE)
      private$.copy_binary_stream(input_con, output_con)
      close(output_con)
      close(input_con)
      data.table::fread(extracted_file, showProgress = show_progress)
    },

    .result_key_for_values = function(response_type, cond_exp_func_model, n, p, betaT, rep, design, inference, inference_type) {
      paste(
        response_type,
        cond_exp_func_model,
        n,
        p,
        betaT,
        rep,
        design,
        inference,
        inference_type,
        sep = "|"
      )
    },

    .result_key = function(rep, design, inference, inference_type) {
      private$.result_key_for_values(
        private$current_response_type,
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
      if (!is.list(row)) {
        stop(".result_key_from_row expected a list, but got: ", typeof(row))
      }
      private$.result_key_for_values(
        row[["response_type"]],
        row[["cond_exp_func_model"]],
        row[["n"]],
        row[["p"]],
        row[["betaT"]],
        row[["rep"]],
        row[["design"]],
        row[["inference"]],
        row[["inference_type"]]
      )
    },

    .result_metadata_dt = function(rep, design, inference, inference_type) {
      data.table::data.table(
        response_type = private$current_response_type,
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
            private$current_response_type %in% c("continuous", "proportion", "count")) {
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
          private$.build_design(design_gen, X, y_linear_model, design_extra, skip_assignment = TRUE),
          error = function(e) NULL
        )
        if (is.null(des_obj)) next
        for (ii in seq_along(private$inference_classes)) {
          inf_gen  = private$inference_classes[[ii]]
          inf_name = private$inference_labels[[ii]]
          inf_ctor_extra = private$inference_constructor_params[[ii]]
          toggle_asserts(TRUE)
          inf_obj = tryCatch(
            do.call(inf_gen$new, c(list(des_obj), inf_ctor_extra)),
            error = function(e) NULL
          )
          toggle_asserts(FALSE)
          if (is.null(inf_obj)) next
          valid_inference_types = private$.valid_inference_types(inf_obj)
          if (length(valid_inference_types) == 0L) next
          for (it in valid_inference_types) {
            # Validate user-supplied params for this inference type early
            private$.args_for_inf_type(inf_obj, it)
            
            combos[[length(combos) + 1L]] = list(
              response_type = private$current_response_type,
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

    .run_single_replication_in_worker = function(rep_i, state, progress_cb = NULL, is_forked = FALSE) {
      # This runs in a worker process. It must be self-contained.
      # 1. Cap threads and nested parallelism to avoid N*M oversubscription.
      # We call these directly as the environment is set to the EDI namespace.
      set_package_threads(1L)
      if (is_forked) {
        unset_num_cores() # Clear any inherited clusters
      }
      
      # Pre-format key prefix for this replication
      rep_key_prefix = paste(state$cell_key_prefix, rep_i, sep = "|")

      # 1. Generate data
      data = generate_covariate_dataset(
        n                    = state$n,
        p                    = state$p,
        cond_exp_func_model  = state$cond_exp_func_model,
        norm_sq_beta_vec     = state$norm_sq_beta_vec,
        X_mat                = state$shared_X %||% state$X_mat,
        cov_draw_method      = if (is.null(state$shared_X) && is.null(state$X_mat)) state$cov_draw_method else NULL,
        cov_draw_method_args = state$cov_draw_method_args
      )
      y_linear_model = as.numeric(scale(data$y_cont))
      X = data$X
      
      # True ATE calculation logic (extracted from compute_true_mean_diff_ate)
      clamp = function(x, lo, hi) {
        pmin(hi, pmax(lo, x))
      }
      true_mean_diff_ate = switch(state$response_type,
        continuous = state$betaT,
        incidence = {
          p_t = clamp(stats::plogis(y_linear_model + state$betaT), state$incidence_clamp, 1 - state$incidence_clamp)
          p_c = clamp(stats::plogis(y_linear_model), state$incidence_clamp, 1 - state$incidence_clamp)
          mean(p_t - p_c)
        },
        proportion = {
          p_t = clamp(stats::plogis(y_linear_model + state$betaT), state$proportion_clamp, 1 - state$proportion_clamp)
          p_c = clamp(stats::plogis(y_linear_model), state$proportion_clamp, 1 - state$proportion_clamp)
          mean(p_t - p_c)
        },
        count = {
          r_t = clamp(exp(y_linear_model + state$betaT), state$count_clamp, Inf)
          r_c = clamp(exp(y_linear_model), state$count_clamp, Inf)
          mean(r_t - r_c)
        },
        NA_real_
      )

      results = list()
      result_keys = character()
      skipped_count = 0L

      # 2. Design and Inference loop
      for (di in seq_along(state$design_classes)) {
        design_gen   = state$design_classes[[di]]
        design_name  = state$design_labels[[di]]
        design_extra = if (!is.null(state$design_params)) state$design_params[[di]] else list()

        # Auto-inject covariate-dependent params (mirrors .build_design)
        init_fn_w = get_r6_init_fn(design_gen)
        if (!is.null(init_fn_w)) {
          fn_formals_w = names(formals(init_fn_w))
          x_names_w    = names(X)
        if ("strata_cols" %in% fn_formals_w &&
            !"strata_cols" %in% names(design_extra) &&
            !identical(design_gen$classname, "FixedDesignBlocking"))
          design_extra$strata_cols = x_names_w[1L]
          if ("cluster_col" %in% fn_formals_w && !"cluster_col" %in% names(design_extra))
            design_extra$cluster_col = x_names_w[min(2L, length(x_names_w))]
          if ("factors"     %in% fn_formals_w && !"factors"     %in% names(design_extra))
            design_extra$factors = list(treatment = 2L)
        }

        # Build design (extracted from .build_design)
        des_obj = tryCatch({
          d = do.call(design_gen$new, c(list(response_type = state$response_type, n = state$n), design_extra))
          if (inherits(d, "DesignSeqOneByOne")) {
            for (t in seq_len(state$n)) {
              w_t = d$add_one_subject_to_experiment_and_assign(X[t, , drop = FALSE])
              out = apply_treatment_and_noise_cpp(
                y_linear_model[t], w_t,
                state$response_type, state$betaT,
                state$sd_noise, state$prob_censoring,
                state$n_ordinal_levels,
                phi_proportion = state$phi_proportion,
                k_survival = state$k_survival,
                incidence_clamp = state$incidence_clamp,
                proportion_clamp = state$proportion_clamp,
                count_clamp = state$count_clamp,
                survival_clamp = state$survival_clamp)
              d$add_one_subject_response(t, out$y, out$dead)
            }
          } else {
            d$add_all_subjects_to_experiment(X)
            d$assign_w_to_all_subjects()
            w = d$get_w()
            out = apply_treatment_and_noise_cpp(
              y_linear_model, w,
              state$response_type, state$betaT,
              state$sd_noise, state$prob_censoring,
              state$n_ordinal_levels,
              phi_proportion = state$phi_proportion,
              k_survival = state$k_survival,
              incidence_clamp = state$incidence_clamp,
              proportion_clamp = state$proportion_clamp,
              count_clamp = state$count_clamp,
              survival_clamp = state$survival_clamp)
            d$add_all_subject_responses(out$y, out$dead)
          }
          d
        }, error = function(e) NULL)
        
        if (is.null(des_obj)) next

        for (ii in seq_along(state$inference_classes)) {
          inf_gen  = state$inference_classes[[ii]]
          inf_name = state$inference_labels[[ii]]
          inf_ctor_extra = state$inference_constructor_params[[ii]]

          inf_obj = tryCatch({
            do.call(inf_gen$new, c(list(des_obj), inf_ctor_extra))
          }, error = function(e) NULL)
          
          if (is.null(inf_obj)) next

          is_mean_diff = is(inf_obj, "InferenceAllSimpleMeanDiff") ||
                         is(inf_obj, "InferenceIncidenceWald") ||
                         is(inf_obj, "InferenceIncidCMH") ||
                         is(inf_obj, "InferenceIncidExtendedRobins") ||
                         is(inf_obj, "InferenceIncidRiskDiff") ||
                         is(inf_obj, "InferenceIncidMiettinenNurminenRiskDiff") ||
                         is(inf_obj, "InferenceIncidNewcombeRiskDiff") ||
                         is(inf_obj, "InferenceIncidKKNewcombeRiskDiff") ||
                         is(inf_obj, "InferenceIncidKKGCompRiskDiff") ||
                         is(inf_obj, "InferencePropGCompMeanDiff")
          te = if (is_mean_diff) true_mean_diff_ate else state$betaT

          # .valid_inference_types logic (extracted)
          valid_inference_types = character(0)
          if (is(inf_obj, "InferenceAsymp")) 
            valid_inference_types = c(valid_inference_types, intersect(state$inf_types, c("asymp_ci", "asymp_pval")))
          if (is(inf_obj, "InferenceExact"))
            valid_inference_types = c(valid_inference_types, intersect(state$inf_types, c("exact_ci", "exact_pval")))
          if (is(inf_obj, "InferenceBoot"))
            valid_inference_types = c(valid_inference_types, intersect(state$inf_types, c("boot_ci", "boot_pval")))
          if (is(inf_obj, "InferenceRand")) {
            valid_inference_types = c(valid_inference_types, intersect(state$inf_types, "rand_pval"))
            if (is(inf_obj, "InferenceRandCI") && state$response_type %in% c("continuous", "proportion", "count"))
              valid_inference_types = c(valid_inference_types, intersect(state$inf_types, "rand_ci"))
          }

          # Result key and pending check
          pending_inference_types = valid_inference_types[!check_in_result_key_store_cpp(
            rep(state$response_type, length(valid_inference_types)),
            rep(state$cond_exp_func_model, length(valid_inference_types)),
            rep(state$n, length(valid_inference_types)),
            rep(state$p, length(valid_inference_types)),
            rep(state$betaT, length(valid_inference_types)),
            rep(rep_i, length(valid_inference_types)),
            rep(design_name, length(valid_inference_types)),
            rep(inf_name, length(valid_inference_types)),
            valid_inference_types
          )]
          
          skipped_count = skipped_count + (length(valid_inference_types) - length(pending_inference_types))
          if (length(pending_inference_types) == 0L) {
            # Advance progress for skipped types
            if (!is.null(progress_cb)) {
               for (k in seq_along(valid_inference_types)) progress_cb()
            }
            next
          }

          est = tryCatch({
            v = inf_obj$compute_estimate()
            if (is.null(v) || length(v) == 0L) NA_real_ else as.numeric(v)[1L]
          }, error = function(e) NA_real_)

          # Helper to merge user-supplied params with defaults
          get_args = function(type, defaults = list()) {
            user_args = state$inference_type_params[[type]]
            if (is.null(user_args)) user_args = list()
            modifyList(defaults, user_args)
          }

          # helper for recording in worker
          local_record = function(type, ci, pval) {
            ci2 = if (length(ci) >= 2L) as.numeric(ci[1:2]) else c(NA_real_, NA_real_)
            if (all(is.finite(ci2)) && ci2[1L] > ci2[2L]) ci2 = rev(ci2)
            
            results[[length(results) + 1L]] <<- list(
              response_type = state$response_type,
              rep           = rep_i,
              cond_exp_func_model = state$cond_exp_func_model,
              n             = as.integer(state$n),
              p             = as.integer(state$p),
              betaT         = as.numeric(state$betaT),
              design        = design_name,
              inference     = inf_name,
              inference_type = type,
              estimate      = if (is.null(est) || !is.finite(est)) NA_real_ else as.numeric(est),
              ci_lo          = ci2[1L],
              ci_hi          = ci2[2L],
              pval          = if (is.null(pval) || length(pval) == 0L || !is.finite(pval[1L])) 
                                NA_real_ else as.numeric(pval[1L]),
              true_estimand = as.numeric(te)
            )
            if (!is.null(progress_cb)) progress_cb()
          }

          # Advance progress for already-present results (cached)
          if (!is.null(progress_cb) && length(valid_inference_types) > length(pending_inference_types)) {
             already_done = length(valid_inference_types) - length(pending_inference_types)
             for (k in seq_len(already_done)) progress_cb()
          }

          # Inference execution logic
          if (is(inf_obj, "InferenceAsymp") && any(c("asymp_ci", "asymp_pval") %in% pending_inference_types)) {
            if ("asymp_pval" %in% pending_inference_types) {
              args = get_args("asymp_pval")
              pval_a = tryCatch(do.call(inf_obj$compute_asymp_two_sided_pval, args), error = function(e) NA_real_)
              local_record("asymp_pval", c(NA_real_, NA_real_), pval_a)
            }
            if ("asymp_ci" %in% pending_inference_types) {
              args = get_args("asymp_ci", list(alpha = state$alpha))
              ci_a = tryCatch(do.call(inf_obj$compute_asymp_confidence_interval, args), error = function(e) c(NA_real_, NA_real_))
              local_record("asymp_ci", ci_a, NA_real_)
            }
          }
          
          if (is(inf_obj, "InferenceExact") && any(c("exact_ci", "exact_pval") %in% pending_inference_types)) {
            if ("exact_pval" %in% pending_inference_types) {
              args = get_args("exact_pval")
              pval_e = tryCatch(do.call(inf_obj$compute_exact_two_sided_pval_for_treatment_effect, args), error = function(e) NA_real_)
              local_record("exact_pval", c(NA_real_, NA_real_), pval_e)
            }
            if ("exact_ci" %in% pending_inference_types) {
              args = get_args("exact_ci", list(alpha = state$alpha))
              ci_e = tryCatch(do.call(inf_obj$compute_exact_confidence_interval, args), error = function(e) c(NA_real_, NA_real_))
              local_record("exact_ci", ci_e, NA_real_)
            }
          }
          
          if (is(inf_obj, "InferenceBoot") && any(c("boot_ci", "boot_pval") %in% pending_inference_types)) {
            if ("boot_pval" %in% pending_inference_types) {
              args = get_args("boot_pval", list(B = state$B_boot, na.rm = TRUE))
              pval_b = tryCatch(do.call(inf_obj$compute_bootstrap_two_sided_pval, args), error = function(e) NA_real_)
              local_record("boot_pval", c(NA_real_, NA_real_), pval_b)
            }
            if ("boot_ci" %in% pending_inference_types) {
              args = get_args("boot_ci", list(B = state$B_boot, alpha = state$alpha, na.rm = TRUE, show_progress = FALSE))
              ci_b = tryCatch(do.call(inf_obj$compute_bootstrap_confidence_interval, args), 
                              error = function(e) c(NA_real_, NA_real_))
              local_record("boot_ci", ci_b, NA_real_)
            }
          }
          
          if (is(inf_obj, "InferenceRand") && any(c("rand_ci", "rand_pval") %in% pending_inference_types)) {
            if ("rand_pval" %in% pending_inference_types) {
              args = get_args("rand_pval", list(r = state$r_rand, na.rm = TRUE, show_progress = FALSE))
              pval_r = tryCatch(do.call(inf_obj$compute_rand_two_sided_pval, args), error = function(e) NA_real_)
              local_record("rand_pval", c(NA_real_, NA_real_), pval_r)
            }
            if ("rand_ci" %in% pending_inference_types && is(inf_obj, "InferenceRandCI") && state$response_type %in% c("continuous", "proportion", "count")) {
              args = get_args("rand_ci", list(r = state$r_rand, alpha = state$alpha, pval_epsilon = state$pval_epsilon, show_progress = FALSE))
              ci_r = tryCatch(do.call(inf_obj$compute_rand_confidence_interval, args), 
                              error = function(e) c(NA_real_, NA_real_))
              local_record("rand_ci", ci_r, NA_real_)
            }
          }
        }
      }
      # Return results as a data.table for efficiency in master loop
      list(
        results_dt = if (length(results) > 0L) data.table::rbindlist(results) else NULL,
        skipped_count = skipped_count
      )
    },

    .advance_progress = function() {
      private$current_task_in_rep_idx = private$current_task_in_rep_idx + 1L
      if (!isTRUE(private$verbose)) return(invisible(NULL))
      
      if (isTRUE(private$use_progress_bar)) {
        private$.draw_simulation_progress_bars()
      }
      invisible(NULL)
    },

    .print_plan_summary = function(planned_combos_list) {
      n_cells = length(planned_combos_list)
      cat("Simulation Plan Summary:\n", file = stderr())
      
      # Group by response_type to keep it concise
      rt_summaries = list()
      for (cell_idx in seq_len(n_cells)) {
        rt = private$param_grid$response_type[[cell_idx]]
        combos = planned_combos_list[[cell_idx]]
        design_names = unique(vapply(combos, `[[`, "", "design"))
        design_names = gsub("Design", "", design_names)
        inference_names = unique(vapply(combos, `[[`, "", "inference"))
        inference_names = gsub("Inference", "", inference_names)
        inference_names = gsub("^(Contin|Count|Incid|Prop|Survival|Ordinal|All)", "", inference_names)

        if (is.null(rt_summaries[[rt]])) {
          rt_summaries[[rt]] = list(
            designs = design_names,
            inferences = inference_names,
            n_tasks = length(combos)
          )
        } else {
          rt_summaries[[rt]]$designs = unique(c(rt_summaries[[rt]]$designs, design_names))
          rt_summaries[[rt]]$inferences = unique(c(rt_summaries[[rt]]$inferences, inference_names))
          rt_summaries[[rt]]$n_tasks = rt_summaries[[rt]]$n_tasks + length(combos)
        }
      }

      response_types = names(rt_summaries)
      if (length(response_types) == 1L) {
        s = rt_summaries[[response_types[[1L]]]]
        cat(sprintf("  Designs (%d): %s\n", length(s$designs), paste(s$designs, collapse = ", ")), file = stderr())
        cat(sprintf("  Inferences (%d): %s\n", length(s$inferences), paste(s$inferences, collapse = ", ")), file = stderr())
      } else {
        for (rt in response_types) {
          s = rt_summaries[[rt]]
          cat(sprintf("  - Response Type: %s\n", rt), file = stderr())
          cat(sprintf("    Designs (%d): %s\n", length(s$designs), paste(s$designs, collapse = ", ")), file = stderr())
          cat(sprintf("    Inferences (%d): %s\n", length(s$inferences), paste(s$inferences, collapse = ", ")), file = stderr())
        }
      }
      cat("\n", file = stderr())
    },

    .draw_progress = function() {
      if (!isTRUE(private$verbose)) return(invisible(NULL))
      if (isTRUE(private$use_progress_bar)) {
        private$.draw_simulation_progress_bars()
      } else if (private$progress_log_interval > 0L &&
                 (private$progress_count %% private$progress_log_interval == 0L ||
                  private$progress_count == private$progress_total)) {
        message(sprintf("Completed %d / %d runs", private$progress_count, private$progress_total))
      }
    },

    .draw_labeled_progress_bar = function(label, prop) {
      width = getOption("width", 80L)
      if (is.null(width) || width < 80L) width = 80L
      
      bar_width = width - nchar(label) - 10L
      if (bar_width < 10L) bar_width = 10L
      
      make_bar = function(p, b_width) {
        pct_str = sprintf(" %3d%% ", floor(p * 100))
        n_pct = nchar(pct_str)
        fill = floor(p * b_width)
        full_bar = paste0(strrep("=", fill), strrep(" ", b_width - fill))
        if (b_width >= n_pct) {
          start_pos = (b_width - n_pct) %/% 2 + 1
          substr(full_bar, start_pos, start_pos + n_pct - 1) = pct_str
        }
        sprintf("[%s]", full_bar)
      }
      
      msg = sprintf("\r%s %s", label, make_bar(prop, bar_width))
      cat(substr(msg, 1, width), file = stderr())
      if (exists("flush.console")) utils::flush.console()
    },

    .message_stderr = function(msg) {
      if (isTRUE(private$progress_bar_drawn)) {
        # Move up 4 lines and clear them to end of screen
        cat("\033[4A\r\033[J", file = stderr())
        private$progress_bar_drawn = NULL
      }
      cat(msg, file = stderr())
      if (exists("flush.console")) utils::flush.console()
    },

    .draw_simulation_progress_bars = function() {
      now = as.numeric(Sys.time())
      # Throttle: only redraw if 100ms have passed OR we are at 100%
      is_done = private$current_cell_idx == private$total_cells &&
                private$current_rep_idx == private$Nrep &&
                private$current_task_in_rep_idx == private$tasks_per_rep

      if (!is_done && (now - private$last_progress_draw_time) < 0.1) return(invisible(NULL))
      private$last_progress_draw_time = now

      # Force a narrow width to prevent wrapping which breaks ANSI sequences.
      width = 60L

      task_total_display = max(private$tasks_per_rep, private$current_task_in_rep_idx)

      # Proportions (clamped to [0,1])
      task_in_progress_prop = max(0, min(1, if (task_total_display > 0) private$current_task_in_rep_idx / task_total_display else 0))
      rep_in_progress_prop  = max(0, min(1, if (private$Nrep > 0) (max(0, private$current_rep_idx - 1) + task_in_progress_prop) / private$Nrep else 0))
      cell_in_progress_prop = max(0, min(1, if (private$total_cells > 0) (max(0, private$current_cell_idx - 1) + rep_in_progress_prop) / private$total_cells else 0))
      
      eta_str = ""
      prop_done_this_run = (cell_in_progress_prop - private$initial_cell_in_progress_prop) / max(1e-6, 1 - private$initial_cell_in_progress_prop)
      prop_remaining     = 1 - cell_in_progress_prop
      
      if (prop_done_this_run > 0 && cell_in_progress_prop < 0.9999) {
        elapsed = now - private$simulation_start_time
        remaining = (elapsed / prop_done_this_run) * (prop_remaining / max(1e-6, 1 - private$initial_cell_in_progress_prop))
        
        d = floor(remaining / 86400)
        remaining = remaining %% 86400
        h = floor(remaining / 3600)
        remaining = remaining %% 3600
        m = floor(remaining / 60)
        s = round(remaining %% 60)
        
        parts = character()
        if (d > 0) parts = c(parts, paste0(d, "d"))
        if (h > 0) parts = c(parts, paste0(h, "h"))
        if (m > 0) parts = c(parts, paste0(m, "m"))
        parts = c(parts, paste0(s, "s"))
        
        eta_str = paste0("Time Left: ", paste(parts, collapse = " "))
      } else if (cell_in_progress_prop >= 0.9999) {
        eta_str = "Status: Completed."
      } else {
        eta_str = "Status: Estimating..."
      }

      make_bar_line = function(label, prop, b_width, digits = 0, label_width = 25) {
        padded_label = sprintf("%-*s", label_width, substr(label, 1, label_width))
        label_len = nchar(padded_label)
        bar_available = b_width - label_len - 2L
        if (bar_available < 10L) return(padded_label)

        fill = floor(prop * bar_available)
        fill = max(0, min(bar_available, fill))

        if (digits == 0) {
          pct_str = sprintf(" %d%% ", floor(prop * 100))
        } else {
          pct_str = sprintf(" %.1f%% ", prop * 100)
        }
        n_pct = nchar(pct_str)

        full_bar = paste0(strrep("=", fill), strrep(" ", bar_available - fill))
        if (bar_available >= n_pct) {
           start_pos = (bar_available - n_pct) %/% 2 + 1
           substr(full_bar, start_pos, start_pos + n_pct - 1) = pct_str
        }
        sprintf("%s[%s]", padded_label, full_bar)
      }

      line1 = make_bar_line(sprintf("DGP %d/%d Overall", private$current_cell_idx, private$total_cells), cell_in_progress_prop, width, digits = 1)
      line2 = make_bar_line(sprintf("Rep %d/%d", private$current_rep_idx, private$Nrep), rep_in_progress_prop, width)
      line3 = make_bar_line(sprintf("%s %d/%d", private$current_task_label, private$current_task_in_rep_idx, task_total_display), task_in_progress_prop, width)

      # Use \r and \033[2K on EACH line for maximum robustness.
      output_block = paste0(
        "\r\033[2K", eta_str, "\n",
        "\r\033[2K", line1, "\n",
        "\r\033[2K", line2, "\n",
        "\r\033[2K", line3, "\n"
      )

      if (is.null(private$progress_bar_drawn)) {
         cat(output_block, sep = "", file = stderr())
         private$progress_bar_drawn = TRUE
      } else {
         # Move up 4 lines, and print the block which clears each line as it goes.
         cat("\033[4A", output_block, sep = "", file = stderr())
      }

      if (exists("flush.console")) utils::flush.console()
    },
    .append_result_row_to_file = function(row) {
      format = private$.results_file_format(private$results_filename)
      if (identical(format, "csv")) {
        file_exists = file.exists(private$results_filename)
        data.table::fwrite(row, private$results_filename, append = file_exists, col.names = !file_exists)
        return(invisible(NULL))
      }
      if (!identical(format, "csv.bz2")) {
        stop("Unsupported results file format: ", private$results_filename)
      }

      staging_filename = private$.results_staging_filename()
      staging_exists = file.exists(staging_filename)
      data.table::fwrite(row, staging_filename, append = staging_exists, col.names = !staging_exists)
      invisible(NULL)
    },
    .sync_results_bz2_from_staging = function(staging_filename = private$.results_staging_filename()) {
      if (!file.exists(staging_filename)) {
        stop("Cannot update compressed results because staging CSV is missing: ", staging_filename)
      }

      results_path = private$.results_output_path()
      results_dir = dirname(results_path)
      if (!dir.exists(results_dir)) {
        dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
      }

      compressed_tmpfile = tempfile(
        pattern = paste0(sub("\\.csv\\.bz2$", "", basename(private$results_filename), ignore.case = TRUE), "_"),
        tmpdir = results_dir,
        fileext = ".csv.bz2"
      )
      on.exit(if (file.exists(compressed_tmpfile)) unlink(compressed_tmpfile), add = TRUE)

      input_con = file(staging_filename, open = "rb")
      output_con = bzfile(compressed_tmpfile, open = "wb")
      on.exit(try(close(input_con), silent = TRUE), add = TRUE)
      on.exit(try(close(output_con), silent = TRUE), add = TRUE)
      private$.copy_binary_stream(input_con, output_con)
      close(output_con)
      close(input_con)

      if (file.exists(results_path))
        unlink(results_path)
      if (!file.rename(compressed_tmpfile, results_path)) {
        stop("Failed to move temporary compressed results into place: ", private$results_filename)
      }
      invisible(NULL)
    },

    .cleanup_results_staging_file = function() {
      if (private$.results_file_format(private$results_filename) != "csv.bz2")
        return(invisible(NULL))
      staging_filename = private$.results_staging_filename()
      if (file.exists(staging_filename))
        unlink(staging_filename)
      invisible(NULL)
    },

    .ensure_staging_file_exists = function() {
      if (private$.results_file_format(private$results_filename) != "csv.bz2")
        return(invisible(NULL))
      
      staging_filename = private$.results_staging_filename()
      if (file.exists(staging_filename))
        return(invisible(NULL))
      
      results_path = private$.results_output_path()
      if (!file.exists(results_path))
        return(invisible(NULL))
      
      if (isTRUE(private$verbose)) message("Extracting existing results to staging file...")
      input_con = bzfile(results_path, open = "rb")
      output_con = file(staging_filename, open = "wb")
      on.exit(try(close(input_con), silent = TRUE), add = TRUE)
      on.exit(try(close(output_con), silent = TRUE), add = TRUE)
      private$.copy_binary_stream(input_con, output_con)
      invisible(NULL)
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
      data$y_linear_model = as.numeric(scale(data$y_cont))
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
      data$y_linear_model = as.numeric(scale(data$y_cont))
      data$y_cont = NULL # SimulationFramework doesn't need the raw cont y anymore
      data
    },

    compute_true_mean_diff_ate = function(y_linear_model) {
      eta_c = y_linear_model
      eta_t = y_linear_model + private$current_betaT

      clamp = function(x, lo, hi) {
        pmin(hi, pmax(lo, x))
      }

      switch(private$current_response_type,
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
    .build_design = function(design_gen, X, y_linear_model, design_extra, skip_assignment = FALSE) {
      n       = private$current_n

      # Auto-inject required args that depend on the covariate matrix when the
      # user has not already supplied them via design_classes_and_params.
      init_fn = get_r6_init_fn(design_gen)
      if (!is.null(init_fn)) {
        fn_formals = names(formals(init_fn))
        x_names    = names(X)
        if ("strata_cols" %in% fn_formals &&
            !"strata_cols" %in% names(design_extra) &&
            !identical(design_gen$classname, "FixedDesignBlocking"))
          design_extra$strata_cols = x_names[1L]
        if ("cluster_col" %in% fn_formals && !"cluster_col" %in% names(design_extra))
          design_extra$cluster_col = x_names[min(2L, length(x_names))]
        if ("factors"     %in% fn_formals && !"factors"     %in% names(design_extra))
          design_extra$factors = list(treatment = 2L)
      }

      des_obj = do.call(design_gen$new, c(
        list(response_type = private$current_response_type, n = n),
        design_extra
      ))

      if (skip_assignment) {
        # Bypass heavy assignment logic during validation phase.
        # Populate the minimum state needed for inference constructors that
        # validate against completed-design metadata such as block IDs or
        # binary-match structure.
        priv = des_obj$.__enclos_env__$private
        priv$Xraw = data.table::as.data.table(X)
        priv$Ximp = data.table::copy(priv$Xraw)
        priv$X = X
        priv$w = rep(c(0L, 1L), length.out = n)
        priv$y = rep(0, n)
        priv$y_original = priv$y
        priv$dead = rep(1L, n)  # 1 = uncensored; 0 would trigger "uncensored responses" asserts
        priv$t = n

        # Some fixed designs derive their blocking structure lazily from the
        # covariates. Build that structure here so validation-time assertions
        # (e.g. CMH / Extended Robins) see the same metadata as a fully
        # realized design.
        if (inherits(des_obj, "FixedDesignBinaryMatch") &&
            private$.has_private_method_on_object(des_obj, "ensure_bms_computed")) {
          priv$ensure_bms_computed()
        } else if (inherits(des_obj, "FixedDesignOptimalBlocks") &&
                   private$.has_private_method_on_object(des_obj, "get_or_compute_block_ids")) {
          priv$m = as.integer(priv$get_or_compute_block_ids())
        } else if (!is.null(priv$strata_cols) &&
                   length(priv$strata_cols) > 0L &&
                   private$.has_private_method_on_object(des_obj, "get_strata_keys")) {
          strata_keys = priv$get_strata_keys()
          if (length(strata_keys) == n) {
            priv$m = match(strata_keys, unique(strata_keys))
          }
        }

        return(des_obj)
      }

      if (inherits(des_obj, "DesignSeqOneByOne")) {
        # Sequential: assignment depends on prior responses so w is obtained
        # one subject at a time.  Call Rcpp with length-1 vectors per subject.
        for (t in seq_len(n)) {
          w_t = des_obj$add_one_subject_to_experiment_and_assign(
            X[t, , drop = FALSE])
          out = apply_treatment_and_noise_cpp(
            y_linear_model[t], w_t,
	            private$current_response_type, private$current_betaT,
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
	          private$current_response_type, private$current_betaT,
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

    # Append multiple rows to raw_results and write to disk in one go.
    .record_batch = function(rows, skipped_count = 0L, keys = NULL) {
      is_dt = data.table::is.data.table(rows)
      n_rows = if (is_dt) nrow(rows) else length(rows)
      if (n_rows == 0L && skipped_count == 0L) return(invisible(NULL))

      if (n_rows > 0L) {
        # 1. Update memory store
        # Optimization: Store data.tables directly in raw_results list
        private$results_idx = private$results_idx + 1L
        private$raw_results[[private$results_idx]] = rows

        # 2. Update seen keys
        if (is_dt) {
          add_to_result_key_store_cpp(
            rows$response_type, rows$cond_exp_func_model, rows$n, rows$p, rows$betaT,
            rows$rep, rows$design, rows$inference, rows$inference_type
          )
        } else {
          # rows is a list of lists - convert to vectors for C++
          add_to_result_key_store_cpp(
            vapply(rows, `[[`, "", "response_type"),
            vapply(rows, `[[`, "", "cond_exp_func_model"),
            vapply(rows, `[[`, 0L, "n"),
            vapply(rows, `[[`, 0L, "p"),
            vapply(rows, `[[`, 0, "betaT"),
            vapply(rows, `[[`, 0L, "rep"),
            vapply(rows, `[[`, "", "design"),
            vapply(rows, `[[`, "", "inference"),
            vapply(rows, `[[`, "", "inference_type")
          )
        }

        # 3. Batch write to disk
        private$.append_result_row_to_file(rows)
      }

      # 4. Advance progress count
      # ONLY add n_rows (the new ones). 
      # DO NOT add skipped_count because those tasks were already accounted for
      # in the initial progress_count calculation (the "scanning" phase).
      private$progress_count = private$progress_count + n_rows
      private$.draw_progress()

      invisible(NULL)
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
      rt   = private$response_type_values[[1L]]
      univ = if (!(rt == "survival" && private$prob_censoring > 0)) list(InferenceAllSimpleMeanDiff) else list()

      type_specific = switch(rt,
        continuous = list(
          InferenceAllSimpleWilcox,
          InferenceContinOLS,
          InferenceContinLin,
          InferenceContinRobustRegr,
          InferenceContinRobustRegr,
          InferenceContinKKOLSIVWC,
          InferenceContinKKOLSOneLik
        ),
        incidence = list(
          InferenceIncidLogRegr,
          InferenceIncidLogRegr,
          InferenceIncidModifiedPoisson,
          InferenceIncidModifiedPoisson,
          InferenceIncidKKClogitIVWC,
          InferenceIncidKKClogitOneLik,
          InferenceIncidCMH,
          InferenceIncidExtendedRobins,
          InferenceIncidExactZhang,
          InferenceIncidExactFisher,
          InferenceIncidenceExactBinomial
        ),
        proportion = list(
          InferenceAllSimpleWilcox,
          InferencePropBetaRegr,
          InferencePropBetaRegr,
          InferencePropFractionalLogit,
          InferencePropFractionalLogit,
          InferencePropKKGEE,
          InferencePropKKQuantileRegrIVWC
        ),
        count = list(
          InferenceAllSimpleWilcox,
          InferenceCountPoisson,
          InferenceCountPoisson,
          InferenceCountRobustPoisson,
          InferenceCountRobustPoisson,
          InferenceCountKKGEE,
          InferenceCountKKCPoissonIVWC
        ),
        survival = list(
          InferenceSurvivalCoxPHRegr,
          InferenceSurvivalCoxPHRegr,
          InferenceSurvivalLogRank,
          InferenceSurvivalRestrictedMeanDiff,
          InferenceSurvivalKKStratCoxIVWC,
          InferenceSurvivalKKLWACoxIVWC
        ),
        ordinal = list(
          InferenceOrdinalPropOddsRegr,
          InferenceOrdinalOrderedProbitRegr,
          InferenceOrdinalCloglogRegr,
          InferenceOrdinalKKGEE
        ),
        stop("Unknown response_type: ", rt)
      )

      c(univ, type_specific)
    }
  )
)
