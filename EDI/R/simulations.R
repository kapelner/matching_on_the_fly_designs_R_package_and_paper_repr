# simulations.R
#
# Simulation framework for comparing experimental designs and inference methods
# in the EDI (Experimental Design with Inference) package.
#
# Usage:
#   devtools::load_all("EDI")
#   source("simulations.R")
#
#   design_cls = list(DesignSeqOneByOneKK21, DesignSeqOneByOneUrn)
#   inf_cls    = list(InferenceContinMultOLS, InferenceContinMultOLSKKIVWC)
#   sim = SimulationFramework$new(
#     response_type    = "continuous",
#     design_classes   = design_cls,
#     inference_classes = inf_cls,
#     n = 100, p = 5, data_type = "linear", Nrep = 50, betaT = 1,
#     design_params    = list(list(lambda = 0.5), list(alpha = 2, beta = 2)),
#     inference_params = list(list(harden = FALSE), list())
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
#' @param data_type Character scalar. Either \code{"linear"} (latent response is a
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
#'   \code{cov_draw_method_args}.  Default \code{stats::runif}.
#'   Must be \code{NULL} when \code{X_mat} is supplied.
#' @param cov_draw_method_args Named list of additional arguments forwarded to
#'   \code{cov_draw_method} beyond the sample-size first argument.
#'   Default \code{list(min = -1, max = +1)}.
#'
#' @return A list with two elements: \code{X} (a data frame of covariates) and
#'   \code{y_cont} (a numeric vector of the latent continuous response).
#'
#' @export
generate_covariate_dataset = function(n, p,
                                      data_type            = c("linear", "nonlinear"),
                                      norm_sq_beta_vec     = 1,
                                      X_mat                = NULL,
                                      cov_draw_method      = stats::runif,
                                      cov_draw_method_args = list(min = -1, max = +1)) {
  data_type = match.arg(data_type)

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

  if (data_type == "linear") {
    beta_x = seq(1, -1, length.out = p)
    beta_x = beta_x * sqrt(norm_sq_beta_vec / sum(beta_x^2))
    y_cont = as.numeric(X_mat %*% beta_x)
  } else {
    if (p < 5L) stop("Friedman nonlinear data_type requires p >= 5.")
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

#' Transform Continuous Signal to Response Type Scale
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
#'   \code{response_type = "ordinal"}.  The continuous signal is cut into this many
#'   equal-frequency bins.  Default \code{4L}.  Ignored for all other response types.
#'
#' @return A numeric vector of base responses on the appropriate scale.
#'
#' @export
transform_cont_y_based_on_response_type = function(y_cont, response_type, n_ordinal_levels = 4L) {
  y_sh = as.numeric(scale(y_cont))
  switch(response_type,
    continuous = y_sh,
    incidence  = stats::plogis(y_sh),
    proportion = {
      y_tmp = y_cont - min(y_cont) + 1e-6
      y_tmp / (max(y_tmp) + 1e-6)
    },
    count      = pmax(1L, round(y_cont - min(y_cont) + 1)),
    survival   = pmax(0.1, y_sh - min(y_sh) + 0.1),
    ordinal    = as.integer(cut(
      y_cont,
      breaks = unique(stats::quantile(y_cont, probs = seq(0, 1, length.out = n_ordinal_levels + 1L))),      include.lowest = TRUE)),
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
#'   \item \code{data_type = "linear"}: the base continuous signal is
#'     \eqn{y = X\beta} where \eqn{\beta} is evenly spaced from 1 to \eqn{-1}.
#'   \item \code{data_type = "nonlinear"}: the Friedman (1991) function
#'     \eqn{10\sin(\pi x_1 x_2) + 20(x_3-0.5)^2 + 10x_4 + 5x_5};
#'     requires \eqn{p \ge 5}.
#' }
#' The continuous base signal is transformed to the scale appropriate for
#' \code{response_type} (logistic for incidence, exponentiated for
#' count/survival, etc.).  Treatment effects are applied per-subject: additive
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
#'     inside the CI (\code{NA} when no CI is available for that method).
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
    #' @param design_classes List of R6 class generators, for example
    #'   \code{list(DesignSeqOneByOneKK21, FixedDesignBernoulli)}).  Each
    #'   generator must be constructable with only \code{response_type} and
    #'   \code{n} (plus any extra params supplied via \code{design_params}).
    #'   \code{NULL} (default) uses the package's standard design set.
    #'   Designs requiring \code{strata_cols}, \code{cluster_col}, or
    #'   \code{factors} have sensible defaults auto-injected (first covariate
    #'   column; second for \code{cluster_col}; \code{list(treatment=2)} for
    #'   \code{factors}) when not supplied via \code{design_params}.
    #'
    #' @param inference_classes List of R6 class generators, for example
    #'   \code{list(InferenceContinMultOLS, InferenceContinMultOLSKKIVWC)}).
    #'   \code{NULL} (default) selects a curated set for the given
    #'   \code{response_type}: several universal classes that work with any
    #'   design, plus representative KK-specific classes (silently skipped for
    #'   non-KK designs at runtime).
    #'
    #' @param n Integer.  Sample size per simulation replication.  Default
    #'   \code{100}.
    #'
    #' @param p Integer.  Number of covariates.  Must be \eqn{\ge 5} when
    #'   \code{data_type = "nonlinear"}.  Default \code{5}.
    #'
    #' @param data_type Character scalar.  How the latent continuous signal is
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
    #'   randomisation-based CIs (\code{compute_confidence_interval_rand}).
    #'   Default \code{0.02}.
    #'
    #' @param sd_noise Numeric \eqn{> 0}. Standard deviation of independent
    #'   Gaussian noise added to each subject's outcome.  Default \code{1}.
    #'
    #' @param n_ordinal_levels Positive integer. Number of ordinal categories when
    #'   \code{response_type = "ordinal"}. Default \code{4}.
    #'
    #' @param norm_sq_beta_vec Positive numeric scalar. The desired squared
    #'   Euclidean norm of the latent linear coefficient vector \eqn{\beta}.
    #'   The generated vector is scaled to match this norm. Default \code{1}.
    #'
    #' @param X_mat Numeric matrix of dimensions \code{n x p}, or \code{NULL} (default).
    #'   If provided, these fixed covariates are used for every replication.
    #'   In this case, \code{cov_draw_method} must be \code{NULL}.
    #'
    #' @param cov_draw_method A function used to draw \code{n * p} i.i.d. covariate
    #'   values for every replication. The function must accept the total number
    #'   of values as its first argument, followed by arguments in
    #'   \code{cov_draw_method_args}.  Default \code{stats::runif}.
    #'   Must be \code{NULL} when \code{X_mat} is supplied.
    #'
    #' @param cov_draw_method_args Named list of additional arguments forwarded to
    #'   \code{cov_draw_method} beyond the sample-size first argument.
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
    #' @param design_params \code{NULL} (default) or a list of lists of the same
    #'   length as \code{design_classes}.  Each element is a named list of
    #'   additional constructor arguments passed to the corresponding design's
    #'   \code{$new()} call (matched positionally).  Use an empty \code{list()}
    #'   for a design that needs no extra arguments.  Example:
    #'   \preformatted{design_classes = list(DesignSeqOneByOneKK21,
    #'                         DesignSeqOneByOneUrn)
    #' design_params  = list(
    #'   list(lambda = 0.5, t_0_pct = 0.1),   # for KK21
    #'   list(alpha  = 2,   beta    = 2)        # for Urn
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
    #'     \item{\code{num_bins_for_continuous_covariate}}{Bin count for
    #'       \code{FixedDesignBlocking}.}
    #'   }
    #'
    #' @param inference_params \code{NULL} (default) or a list of lists of the
    #'   same length as \code{inference_classes}.  Each element is a named list
    #'   of arguments for the corresponding inference class (matched
    #'   positionally).  Use an empty \code{list()} for a class that needs no
    #'   extra arguments.
    #'
    #'   Arguments are routed automatically to the appropriate call site based
    #'   on argument name (see \code{.INF_ARG_ROUTES}):
    #'   \itemize{
    #'     \item constructor arguments such as \code{harden} and
    #'       \code{max_resample_attempts}
    #'     \item asymptotic p-value arguments such as \code{delta}
    #'     \item bootstrap arguments such as \code{na.rm} and \code{type}
    #'     \item randomization p-value arguments such as \code{delta},
    #'       \code{transform_responses}, and \code{permutations}
    #'     \item randomization CI arguments such as \code{ci_search_control}
    #'     \item shared randomization arguments such as \code{args_for_type}
    #'   }
    #'   Example:
    #'   \preformatted{inference_classes = list(InferenceContinMultOLS,
    #'                            InferenceContinMultOLSKKIVWC)
    #' inference_params = lapply(
    #'   list(
    #'     list(harden = FALSE), 
    #'     list()
    #'   ),
    #'   append, 
    #'   list(delta = 0)
    #' )}
    #'
    #' @param inf_types Character vector controlling which inference outputs are
    #'   computed and recorded.  \code{NULL} (default) runs all eight types.
    #'   Valid elements are \code{"asymp_ci"}, \code{"asymp_pval"},
    #'   \code{"exact_ci"}, \code{"exact_pval"}, \code{"boot_ci"},
    #'   \code{"boot_pval"}, \code{"rand_ci"}, and \code{"rand_pval"}.
    #'   When no \code{*_ci} type is requested, \code{coverage} is omitted from
    #'   \code{$summarize()}.  When no \code{*_pval} type is requested,
    #'   \code{power} is omitted.
    initialize = function(
      response_type,
      design_classes    = NULL,
      inference_classes = NULL,
      n                 = 100L,
      p                 = 5L,
      data_type         = "linear",
      Nrep              = 100L,
      betaT             = 1,
      alpha             = 0.05,
      B_boot            = 201L,
      r_rand            = 201L,
      pval_epsilon      = 0.02,
      sd_noise              = 1,
      n_ordinal_levels      = 4L,
      norm_sq_beta_vec      = 1,
      X_mat                 = NULL,
      cov_draw_method       = stats::runif,
      cov_draw_method_args  = list(min = -1, max = +1),
      prob_censoring        = 0.2,
      verbose                    = FALSE,
      keep_all_intermediate_data = FALSE,
      turn_off_asserts_for_speed = TRUE,
      design_params              = NULL,
      inference_params           = NULL,
      inf_types                  = NULL
    ) {
      valid_rt = c("continuous", "incidence", "proportion",
                   "count", "survival", "ordinal")
      if (!response_type %in% valid_rt)
        stop("response_type must be one of: ", paste(valid_rt, collapse = ", "))
      if (!data_type %in% c("linear", "nonlinear"))
        stop("data_type must be 'linear' or 'nonlinear'")
      if (data_type == "nonlinear" && p < 5L)
        stop("nonlinear Friedman function requires p >= 5")

      valid_inf_types = c("asymp_ci", "asymp_pval", "exact_ci", "exact_pval",
                          "boot_ci",  "boot_pval",  "rand_ci",  "rand_pval")
      if (is.null(inf_types)) {
        inf_types = valid_inf_types
      } else {
        bad = setdiff(inf_types, valid_inf_types)
        if (length(bad))
          stop("Invalid inf_types: ", paste(bad, collapse = ", "),
               ".  Valid values: ", paste(valid_inf_types, collapse = ", "))
        inf_types = unique(inf_types)
      }

      private$response_type    = response_type
      private$n                = as.integer(n)
      private$p                = as.integer(p)
      private$data_type        = data_type
      private$Nrep             = as.integer(Nrep)
      private$betaT            = betaT
      private$alpha            = alpha
      private$B_boot           = as.integer(B_boot)
      private$r_rand           = as.integer(r_rand)
      private$pval_epsilon     = pval_epsilon
      private$sd_noise             = sd_noise
      private$n_ordinal_levels     = as.integer(n_ordinal_levels)
      private$norm_sq_beta_vec     = norm_sq_beta_vec
      private$X_mat                = X_mat
      private$cov_draw_method      = cov_draw_method
      private$cov_draw_method_args = cov_draw_method_args
      private$prob_censoring       = prob_censoring
      private$verbose                    = verbose
      private$turn_off_asserts_for_speed = turn_off_asserts_for_speed
      private$keep_all_intermediate_data = keep_all_intermediate_data
      private$inf_types        = inf_types

      private$want_asymp_ci   = "asymp_ci"   %in% inf_types
      private$want_asymp_pval = "asymp_pval" %in% inf_types
      private$want_exact_ci   = "exact_ci"   %in% inf_types
      private$want_exact_pval = "exact_pval" %in% inf_types
      private$want_boot_ci    = "boot_ci"    %in% inf_types
      private$want_boot_pval  = "boot_pval"  %in% inf_types
      private$want_rand_ci    = "rand_ci"    %in% inf_types
      private$want_rand_pval  = "rand_pval"  %in% inf_types

      private$design_classes = if (is.null(design_classes)) {
        private$.default_design_classes()
      } else {
        design_classes
      }

      private$inference_classes = if (is.null(inference_classes)) {
        private$.default_inference_classes()
      } else {
        inference_classes
      }

      n_des = length(private$design_classes)
      n_inf = length(private$inference_classes)

      if (!is.null(design_params)) {
        if (!is.list(design_params) || length(design_params) != n_des)
          stop("design_params must be a list of lists of length ", n_des,
               " (one per design class)")
        for (i in seq_len(n_des))
          if (!is.list(design_params[[i]]))
            stop("design_params[[", i, "]] must be a named list")
      }

      if (!is.null(inference_params)) {
        if (!is.list(inference_params) || length(inference_params) != n_inf)
          stop("inference_params must be a list of lists of length ", n_inf,
               " (one per inference class)")
        for (i in seq_len(n_inf))
          if (!is.list(inference_params[[i]]))
            stop("inference_params[[", i, "]] must be a named list")
      }

      private$design_params    = design_params
      private$inference_params = inference_params
      private$design_labels    = private$.compute_design_labels()

      private$raw_results = list()
      private$has_run     = FALSE
    },

    # ── run() ─────────────────────────────────────────────────────────────────
    #' @description
    #' Execute the simulation replications.
    #'
    #' @return The \code{SimulationFramework} object itself (invisibly).
    run = function() {
      # Disable assertions for the duration of the simulation for speed
      if (private$turn_off_asserts_for_speed){
        toggle_asserts(FALSE)
        on.exit(toggle_asserts(TRUE))        
      }

      n_des = length(private$design_classes)
      n_inf = length(private$inference_classes)
      n_met = length(private$inf_types)
      
      # Pre-allocate results list to maximum possible size to avoid re-allocations
      private$raw_results = vector("list", private$Nrep * n_des * n_inf * n_met)
      private$results_idx = 0L

      private$all_intermediate_data = vector("list", private$Nrep)
      private$valid_combos          = list()
      private$seen_combo_keys      = character(0L)
      private$exact_warned_classes = character(0L)

      n_des = length(private$design_classes)
      n_inf = length(private$inference_classes)

      message(sprintf(
        "SimulationFramework: response_type=%s  data_type=%s  n=%d  p=%d  Nrep=%d  betaT=%g",
        private$response_type, private$data_type,
        private$n, private$p, private$Nrep, private$betaT
      ))
      message(sprintf("  Designs: %d   Inference classes: %d", n_des, n_inf))

      log_interval = max(1L, private$Nrep %/% 10L)

      for (rep in seq_len(private$Nrep)) {
        if (rep %% log_interval == 0L || private$verbose)
          message(sprintf("  Rep %d / %d", rep, private$Nrep))

        rep_data = private$.generate_data()
        X      = rep_data$X
        y_base = rep_data$y_base

        # Per-replication true estimand for InferenceAllSimpleMeanDiff.
        # The DGP applies betaT on a scale that may differ from the mean-difference
        # scale that InferenceAllSimpleMeanDiff estimates, so we derive the correct
        # comparison target analytically from y_base for each response type.
        #
        #   incidence : logit-shift  -> E[plogis(qlogis(p_b)+betaT) - p_b]
        #   count     : log-mult.    -> mean(y_base)*exp(sd^2/2)*(exp(betaT)-1)
        #   survival  : log-mult.    -> same as count, then scaled by (1-pc/2)
        #                              to account for uniform censoring of observed times
        #   otherwise : additive DGP -> betaT
        true_mean_diff_ate = switch(private$response_type,
          incidence = mean(stats::plogis(stats::qlogis(y_base) + private$betaT) - y_base),
          count     = mean(y_base) * exp(private$sd_noise^2 / 2) * (exp(private$betaT) - 1),
          survival  = mean(y_base) * exp(private$sd_noise^2 / 2) * (exp(private$betaT) - 1) *
                        (1 - private$prob_censoring / 2),
          private$betaT   # continuous, proportion, ordinal: additive DGP
        )

        for (di in seq_along(private$design_classes)) {
          design_gen   = private$design_classes[[di]]
          design_name  = private$design_labels[[di]]
          design_extra = if (!is.null(private$design_params)) private$design_params[[di]] else list()

          des_obj = tryCatch(
            private$.build_design(design_gen, X, y_base, design_extra),
            error = function(e) {
              if (private$verbose)
                message(sprintf("    [SKIP design %s]: %s", design_name, e$message))
              NULL
            }
          )
          if (private$keep_all_intermediate_data) {
            private$all_intermediate_data[[rep]]$y_base                    = y_base
            private$all_intermediate_data[[rep]]$designs[[design_name]]    = des_obj
          }
          if (is.null(des_obj)) next

          for (ii in seq_along(private$inference_classes)) {
            inf_gen  = private$inference_classes[[ii]]
            inf_name = inf_gen$classname

            # Route per-class params to their call-site destinations
            class_p = if (!is.null(private$inference_params)) private$inference_params[[ii]] else list()
            routed  = private$.route_inference_args(class_p)

            inf_obj = tryCatch({
              inf_extra = private$.filter_by_formals(inf_gen, routed$inf_init)
              do.call(inf_gen$new, c(list(des_obj), inf_extra))
            }, error = function(e) {
              if (private$verbose)
                message(sprintf("    [SKIP %s / %s]: %s",
                                inf_name, design_name, e$message))
              NULL
            })
            if (is.null(inf_obj)) next
            if (private$keep_all_intermediate_data)
              private$all_intermediate_data[[rep]]$inferences[[design_name]][[inf_name]] = inf_obj

            # True estimand on the scale this inference class estimates.
            # InferenceAllSimpleMeanDiff targets the mean-difference ATE, whose
            # true value must be derived analytically for non-additive DGPs
            # (incidence, count, survival).  All other classes target betaT directly.
            te = if (is(inf_obj, "InferenceAllSimpleMeanDiff"))
                   true_mean_diff_ate
                 else
                   private$betaT

            # Track valid (design, inference, method) combos — once per unique pair
            combo_key = paste0(design_name, "|||", inf_name)
            if (!combo_key %in% private$seen_combo_keys) {
              private$seen_combo_keys = c(private$seen_combo_keys, combo_key)
              valid_methods = character(0L)
              if (is(inf_obj, "InferenceAsymp") && (private$want_asymp_ci || private$want_asymp_pval))
                valid_methods = c(valid_methods, "asymp")
              if (is(inf_obj, "InferenceExact") && (private$want_exact_ci || private$want_exact_pval))
                valid_methods = c(valid_methods, "exact")
              if (is(inf_obj, "InferenceBoot") && (private$want_boot_ci || private$want_boot_pval))
                valid_methods = c(valid_methods, "boot")
              if (is(inf_obj, "InferenceRand") && (private$want_rand_ci || private$want_rand_pval))
                valid_methods = c(valid_methods, "rand")
              for (m in valid_methods)
                private$valid_combos[[length(private$valid_combos) + 1L]] =
                  list(design = design_name, inference = inf_name, method = m)
            }

            # Point estimate (shared across inference methods)
            est = tryCatch({
              v = inf_obj$compute_treatment_estimate()
              if (is.null(v) || length(v) == 0L) NA_real_ else as.numeric(v)[1L]
            }, error = function(e) NA_real_)

            # ── Asymptotic ────────────────────────────────────────────────────
            if (is(inf_obj, "InferenceAsymp") &&
                (private$want_asymp_ci || private$want_asymp_pval)) {
              pval_a = if (private$want_asymp_pval) {
                tryCatch(
                  do.call(inf_obj$compute_asymp_two_sided_pval_for_treatment_effect,
                          routed$asymp_pval),
                  error = function(e) NA_real_
                )
              } else NA_real_
              ci_a = if (private$want_asymp_ci) {
                tryCatch(
                  inf_obj$compute_asymp_confidence_interval(alpha = private$alpha),
                  error = function(e) c(NA_real_, NA_real_)
                )
              } else c(NA_real_, NA_real_)
              private$.record(rep, design_name, inf_name, "asymp", est, ci_a, pval_a, te)
            }

            # ── Exact ─────────────────────────────────────────────────────────
            if (private$want_exact_ci || private$want_exact_pval) {
              if (!is(inf_obj, "InferenceExact")) {
                if (!inf_name %in% private$exact_warned_classes) {
                  warning(sprintf(
                    "'%s' does not inherit InferenceExact; exact_* inference will be skipped for this class.",
                    inf_name))
                  private$exact_warned_classes =
                    c(private$exact_warned_classes, inf_name)
                }
              } else {
                pval_e = if (private$want_exact_pval) {
                  tryCatch(
                    inf_obj$compute_exact_two_sided_pval_for_treatment_effect(),
                    error = function(e) NA_real_
                  )
                } else NA_real_
                ci_e = if (private$want_exact_ci) {
                  tryCatch(
                    inf_obj$compute_exact_confidence_interval(alpha = private$alpha),
                    error = function(e) c(NA_real_, NA_real_)
                  )
                } else c(NA_real_, NA_real_)
                private$.record(rep, design_name, inf_name, "exact", est, ci_e, pval_e, te)
              }
            }

            # ── Bootstrap ─────────────────────────────────────────────────────
            if (is(inf_obj, "InferenceBoot") &&
                (private$want_boot_ci || private$want_boot_pval)) {
              pval_b = if (private$want_boot_pval) {
                tryCatch(
                  do.call(inf_obj$compute_bootstrap_two_sided_pval,
                          modifyList(list(B = private$B_boot, na.rm = TRUE),
                                     routed$boot_pval)),
                  error = function(e) NA_real_
                )
              } else NA_real_
              ci_b = if (private$want_boot_ci) {
                tryCatch(
                  do.call(inf_obj$compute_bootstrap_confidence_interval,
                          modifyList(list(B     = private$B_boot,
                                          alpha = private$alpha,
                                          na.rm = TRUE,
                                          show_progress = FALSE),
                                     routed$boot_ci)),
                  error = function(e) c(NA_real_, NA_real_)
                )
              } else c(NA_real_, NA_real_)
              private$.record(rep, design_name, inf_name, "boot", est, ci_b, pval_b, te)
            }

            # ── Randomisation ─────────────────────────────────────────────────
            if (is(inf_obj, "InferenceRand") &&
                (private$want_rand_ci || private$want_rand_pval)) {
              pval_r = if (private$want_rand_pval) {
                tryCatch(
                  do.call(inf_obj$compute_two_sided_pval_for_treatment_effect_rand,
                          modifyList(list(r             = private$r_rand,
                                          na.rm         = TRUE,
                                          show_progress = FALSE),
                                     routed$rand_pval)),
                  error = function(e) NA_real_
                )
              } else NA_real_
              ci_r = if (
                private$want_rand_ci &&
                is(inf_obj, "InferenceRandCI") &&
                private$response_type %in% c("continuous", "proportion", "count")
              ) {
                tryCatch(
                  do.call(inf_obj$compute_confidence_interval_rand,
                          modifyList(list(r             = private$r_rand,
                                          alpha         = private$alpha,
                                          pval_epsilon  = private$pval_epsilon,
                                          show_progress = FALSE),
                                     routed$rand_ci)),
                  error = function(e) c(NA_real_, NA_real_)
                )
              } else {
                c(NA_real_, NA_real_)
              }
              private$.record(rep, design_name, inf_name, "rand", est, ci_r, pval_r, te)
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
    #'     \item{\code{y_base}}{Numeric vector of base responses for that replication.}
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
    #'   inference class, inference method).
    get_results = function() {
      if (!private$has_run) stop("Call $run() first.")
      if (private$results_idx == 0L)
        return(data.table::data.table(
          rep = integer(), design = character(), inference = character(),
          method = character(), estimate = numeric(),
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
    #'   power) grouped by design, inference class, and inference method.
    summarize = function() {
      if (!private$has_run) stop("Call $run() first.")

      # ── Reference grid of all valid (design, inference, method) combos ────────
      if (length(private$valid_combos) == 0L) {
        message("No results."); return(invisible(NULL))
      }
      ref_grid = data.table::rbindlist(lapply(private$valid_combos, as.list))

      # ── Per-class params strings ───────────────────────────────────────────────
      inf_names = vapply(private$inference_classes, function(g) g$classname, "")
      design_params_map = stats::setNames(
        lapply(seq_along(private$design_classes), function(di)
          private$.params_to_str(if (!is.null(private$design_params)) private$design_params[[di]] else NULL)),
        private$design_labels)
      inf_params_map = stats::setNames(
        lapply(seq_along(private$inference_classes), function(ii)
          private$.params_to_str(if (!is.null(private$inference_params)) private$inference_params[[ii]] else NULL)),
        inf_names)

      ref_grid[, design_params    := unlist(design_params_map[design])]
      ref_grid[, inference_params := unlist(inf_params_map[inference])]

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
            MSE   = if (nrow(est_ok)) round(mean((est_ok$estimate - est_ok$true_estimand)^2), 5L) else NA_real_,
            n_est = nrow(est_ok)
          )
          if (report_cov) {
            row$coverage  = if (nrow(ci_ok)) round(mean(ci_ok$ci_lo <= ci_ok$true_estimand & ci_ok$true_estimand <= ci_ok$ci_hi), 3L) else NA_real_
            row$n_cov     = nrow(ci_ok)
            row$ci_length = if (nrow(ci_ok)) round(mean(ci_ok$ci_hi - ci_ok$ci_lo), 5L) else NA_real_
          }
          if (report_pow) {
            row$power = if (length(pv_ok)) round(mean(pv_ok < alpha), 3L) else NA_real_
            row$n_pow = length(pv_ok)
          }
          row
        }, by = .(design, inference, method),
           .SDcols = c("estimate", "ci_lo", "ci_hi", "pval", "true_estimand")]
      } else {
        agg = data.table::data.table(design = character(), inference = character(),
                                     method = character())
      }

      # ── Right-join: every valid combo appears, NA for those with no data ──────
      result = agg[ref_grid, on = .(design, inference, method)]
      result[order(design, inference, method)]
    },

    # ── print() ───────────────────────────────────────────────────────────────
    #' @description
    #' Print a summary of the \code{SimulationFramework} configuration and status.
    print = function() {
      cat("SimulationFramework\n")
      cat("  response_type :", private$response_type, "\n")
      cat("  data_type     :", private$data_type, "\n")
      cat("  n / p         :", private$n, "/", private$p, "\n")
      cat("  Nrep / betaT  :", private$Nrep, "/", private$betaT, "\n")
      cat("  alpha / B_boot / r_rand :",
          private$alpha, "/", private$B_boot, "/", private$r_rand, "\n")
      cat("  inf_types     :", paste(private$inf_types, collapse = ", "), "\n")
      if (!is.null(private$design_params))
        cat("  design_params : (positional list of length", length(private$design_params), ")\n")
      if (!is.null(private$inference_params))
        cat("  inference_params: (positional list of length", length(private$inference_params), ")\n")
      design_names = vapply(private$design_classes,  function(g) g$classname, "")
      inf_names    = vapply(private$inference_classes, function(g) g$classname, "")
      cat("  Designs (", length(design_names), "):",
          paste(design_names, collapse = ", "), "\n")
      cat("  Inference (", length(inf_names), "):",
          paste(inf_names,    collapse = ", "), "\n")
      if (private$has_run) {
        cat("  Status : completed\n")
        sm = self$summarize()
        if (!is.null(sm) && nrow(sm) > 0L) {
          cat(sprintf("\nSummary  (betaT = %g,  alpha = %g):\n",
                      private$betaT, private$alpha))
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
    n                = NULL,
    p                = NULL,
    data_type        = NULL,
    Nrep             = NULL,
    betaT            = NULL,
    alpha            = NULL,
    B_boot           = NULL,
    r_rand           = NULL,
    pval_epsilon     = NULL,
    sd_noise             = NULL,
    n_ordinal_levels     = NULL,
    norm_sq_beta_vec     = NULL,
    X_mat                = NULL,
    cov_draw_method      = NULL,
    cov_draw_method_args = NULL,
    prob_censoring       = NULL,
    verbose          = NULL,
    design_params    = NULL,
    inference_params = NULL,
    inf_types        = NULL,
    want_asymp_ci    = NULL,
    want_asymp_pval  = NULL,
    want_exact_ci    = NULL,
    want_exact_pval  = NULL,
    want_boot_ci     = NULL,
    want_boot_pval   = NULL,
    want_rand_ci     = NULL,
    want_rand_pval   = NULL,
    design_classes   = NULL,
    inference_classes = NULL,
    design_labels    = NULL,
    turn_off_asserts_for_speed = NULL,
    raw_results               = NULL,
    results_idx               = 0L,
    all_intermediate_data     = NULL,
    keep_all_intermediate_data = FALSE,
    has_run                   = FALSE,
    exact_warned_classes      = NULL,
    valid_combos              = NULL,
    seen_combo_keys  = NULL,

    # ── Inference arg routing ─────────────────────────────────────────────────
    #
    # Maps argument names found in inference_params to the call site(s) that
    # accept them.  Destinations:
    #   inf_init   → inf_gen$new()                          (filtered by formals)
    #   asymp_pval → compute_asymp_two_sided_pval_for_treatment_effect()
    #   boot_ci    → compute_bootstrap_confidence_interval()
    #   boot_pval  → compute_bootstrap_two_sided_pval()
    #   rand_pval  → compute_two_sided_pval_for_treatment_effect_rand()
    #   rand_ci    → compute_confidence_interval_rand()
    .INF_ARG_ROUTES = list(
      # constructor
      harden                = c("inf_init"),
      max_resample_attempts = c("inf_init"),
      # p-value null hypothesis — all three methods
      delta                 = c("asymp_pval", "boot_pval", "rand_pval"),
      # bootstrap
      na.rm                 = c("boot_ci", "boot_pval", "rand_pval"),
      type                  = c("boot_ci", "boot_pval", "rand_pval", "rand_ci"),
      # randomisation p-value
      transform_responses   = c("rand_pval"),
      permutations          = c("rand_pval"),
      zero_one_logit_clamp  = c("rand_pval"),
      # shared rand
      args_for_type         = c("rand_pval", "rand_ci"),
      # randomisation CI
      ci_search_control     = c("rand_ci")
    ),

    # Route a merged inference-params list to per-destination sub-lists.
    # Unknown arg names produce a warning; known args fan out to all their
    # destinations.
    .route_inference_args = function(params) {
      dest_keys = c("inf_init", "asymp_pval", "boot_ci",
                    "boot_pval", "rand_pval", "rand_ci")
      out = stats::setNames(
        lapply(seq_along(dest_keys), function(...) list()),
        dest_keys
      )
      if (length(params) == 0L) return(out)
      for (nm in names(params)) {
        dests = private$.INF_ARG_ROUTES[[nm]]
        if (is.null(dests)) {
          warning(sprintf(
            "SimulationFramework: unrecognised inference_params arg '%s'; ignored.", nm))
          next
        }
        for (d in dests) out[[d]][[nm]] = params[[nm]]
      }
      out
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

    # ── Data generation ───────────────────────────────────────────────────────

    .generate_data = function() {
      data = generate_covariate_dataset(
        n                    = private$n,
        p                    = private$p,
        data_type            = private$data_type,
        norm_sq_beta_vec     = private$norm_sq_beta_vec,
        X_mat                = private$X_mat,
        cov_draw_method      = private$cov_draw_method,
        cov_draw_method_args = private$cov_draw_method_args
      )
      data$y_base = transform_cont_y_based_on_response_type(
        y_cont           = data$y_cont,
        response_type    = private$response_type,
        n_ordinal_levels = private$n_ordinal_levels
      )
      data$y_cont = NULL # SimulationFramework doesn't need the raw cont y anymore
      data
    },

    # Instantiate design and run the full experiment (assign + observe all n).
    .build_design = function(design_gen, X, y_base, design_extra) {
      n       = private$n

      # Auto-inject required args that depend on the covariate matrix when the
      # user has not already supplied them via design_params.
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
            y_base[t], w_t,
            private$response_type, private$betaT,
            private$sd_noise, private$prob_censoring,
            private$n_ordinal_levels)
          des_obj$add_one_subject_response(t, out$y, out$dead)
        }
      } else {
        # Fixed: all assignments known upfront — vectorize across all n subjects.
        des_obj$add_all_subjects_to_experiment(X)
        des_obj$assign_w_to_all_subjects()
        w   = des_obj$get_w()
        out = apply_treatment_and_noise_cpp(
          y_base, w,
          private$response_type, private$betaT,
          private$sd_noise, private$prob_censoring,
          private$n_ordinal_levels)
        des_obj$add_all_subject_responses(out$y, out$dead)
      }
      des_obj
    },

    # Append one row to raw_results as a plain list (no data.table allocation).
    # rbindlist() in get_results() converts the whole list at once.
    .record = function(rep, design, inf_name, method, est, ci, pval, true_estimand) {
      ci2 = if (length(ci) >= 2L) as.numeric(ci[1:2]) else c(NA_real_, NA_real_)
      if (all(is.finite(ci2)) && ci2[1L] > ci2[2L]) ci2 = rev(ci2)
      
      private$results_idx = private$results_idx + 1L
      private$raw_results[[private$results_idx]] = list(
        rep           = as.integer(rep),
        design        = design,
        inference     = inf_name,
        method        = method,
        estimate      = if (is.null(est) || !is.finite(est)) NA_real_ else as.numeric(est),
        ci_lo         = ci2[1L],
        ci_hi         = ci2[2L],
        pval          = if (is.null(pval) || length(pval) == 0L || !is.finite(pval[1L]))
                          NA_real_ else as.numeric(pval[1L]),
        true_estimand = as.numeric(true_estimand)
      )
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
          InferenceContinMultOLSKKCombinedLikelihood
        ),
        incidence = list(
          InferenceIncidUnivLogRegr,
          InferenceIncidMultiLogRegr,
          InferenceIncidUnivModifiedPoisson,
          InferenceIncidMultiModifiedPoisson,
          InferenceIncidUnivKKClogitIVWC,
          InferenceIncidMultiKKClogitCombinedLikelihood,
          InferenceIncidAzriel,
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
