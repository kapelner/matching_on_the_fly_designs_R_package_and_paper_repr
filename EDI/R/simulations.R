# simulations.R
#
# Simulation framework for comparing experimental designs and inference methods
# in the EDI (Experimental Design with Inference) package.
#
# Usage:
#   devtools::load_all("EDI")
#   source("simulations.R")
#
#   sim = SimulationFramework$new(
#     response_type = "continuous",
#     n = 100, p = 5, data_type = "linear", Nrep = 50, betaT = 1,
#     design_params = list(
#       DesignSeqOneByOneKK21 = list(lambda = 0.5),
#       DesignSeqOneByOneUrn  = list(alpha = 2, beta = 2)
#     ),
#     inference_params = list(
#       all_inference          = list(delta = 0, na.rm = TRUE),
#       InferenceContinMultOLS = list(harden = FALSE)
#     )
#   )
#   sim$run()
#   sim$summarize()

#' Simulation Framework for Experimental Designs and Inference Methods
#'
#' @description
#' An R6 class that benchmarks experimental designs and inference methods via
#' Monte Carlo simulation.  For each replication it generates synthetic
#' covariate and response data, runs every requested \code{(design, inference)}
#' combination, and collects point estimates, confidence intervals, and
#' p-values.  After \code{Nrep} replications the aggregated metrics—MSE,
#' coverage, and power—are available through \code{$summarize()}.
#'
#' @details
#' \strong{Data generation}\cr
#' Covariates are drawn i.i.d. \eqn{\text{Uniform}(0,1)}.
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
#' \strong{Inference methods run}\cr
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
#' \strong{Reported metrics} (from \code{$summarize()})\cr
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
    #' @param design_classes List of R6 class generators (e.g.\
    #'   \code{list(DesignSeqOneByOneKK21, FixedDesignBernoulli)}).  Each
    #'   generator must be constructable with only \code{response_type} and
    #'   \code{n} (plus any extra params supplied via \code{design_params}).
    #'   \code{NULL} (default) uses all 16 standard designs: the five non-KK
    #'   sequential designs (\code{Bernoulli}, \code{iBCRD}, \code{Efron},
    #'   \code{Atkinson}, \code{Urn}), the two KK sequential designs
    #'   (\code{KK21}, \code{KK14}), and nine fixed designs (\code{Bernoulli},
    #'   \code{iBCRD}, \code{BinaryMatch}, \code{Greedy},
    #'   \code{MatchingGreedyPairSwitching}, \code{Rerandomization},
    #'   \code{DOptimal}, \code{AOptimal}, \code{OptimalBlocks}).
    #'
    #' @param inference_classes List of R6 class generators (e.g.\
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
    #' @param sd_noise Numeric \eqn{> 0}.  Standard deviation of i.i.d.\
    #'   Gaussian noise added to each subject's outcome.  Default \code{1}.
    #'
    #' @param prob_censoring Numeric in \eqn{[0,1]}.  Per-subject independent
    #'   censoring probability; applied only when
    #'   \code{response_type = "survival"}.  Default \code{0.2}.
    #'
    #' @param verbose Logical.  If \code{TRUE}, prints a message for every
    #'   replication and for every \code{(design, inference)} pair that is
    #'   skipped due to an error.  Default \code{FALSE}.
    #'
    #' @param design_params Named list of lists.  Each name must exactly match a
    #'   design class name present in \code{design_classes}, and the
    #'   corresponding value is a named list of additional constructor arguments
    #'   passed to that design's \code{$new()} call.  Design classes whose names
    #'   are absent receive no extra arguments.  Example:
    #'   \preformatted{design_params = list(
    #'     DesignSeqOneByOneKK21 = list(lambda = 0.5, t_0_pct = 0.1),
    #'     DesignSeqOneByOneKK14 = list(lambda = 0.3, morrison = TRUE),
    #'     DesignSeqOneByOneUrn  = list(alpha = 2, beta = 2)
    #'   )}
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
    #' @param inference_params Named list of lists.  Controls constructor
    #'   arguments and inference method call arguments for specific inference
    #'   classes, as well as a universal fallback entry.
    #'
    #'   Each non-\code{all_inference} name must exactly match an inference
    #'   class name present in \code{inference_classes}.  The corresponding
    #'   value is a named list of constructor arguments passed \emph{only} to
    #'   that class's \code{$new()} call, after filtering to args actually
    #'   present in its \code{initialize()} signature.
    #'
    #'   The special key \code{all_inference} provides arguments that apply to
    #'   \emph{every} inference class and \emph{every} inference method call.
    #'   These are routed automatically based on argument name:
    #'   \describe{
    #'     \item{Constructor args (to \code{$new()})}{\code{harden},
    #'       \code{max_resample_attempts} — filtered per class.}
    #'     \item{Asymptotic p-value}{\code{delta} — null treatment effect
    #'       (default \code{0}).}
    #'     \item{Bootstrap CI and p-value}{\code{na.rm} — drop \code{NA}
    #'       replicates (framework default \code{TRUE}); \code{type} — bootstrap
    #'       variant (\code{"perc"}, \code{"bca"}, …).}
    #'     \item{Randomisation p-value}{\code{delta}; \code{na.rm};
    #'       \code{transform_responses} — \code{"none"} | \code{"logit"} |
    #'       \code{"log"}; \code{permutations} — pre-computed permutation matrix;
    #'       \code{zero_one_logit_clamp}.}
    #'     \item{Randomisation CI}{\code{ci_search_control} — list controlling
    #'       the bisection search strategy.}
    #'     \item{Randomisation p-value and CI}{\code{type};
    #'       \code{args_for_type} — e.g.\
    #'       \code{list(Zhang = list(combination_method = "Fisher"))}.}
    #'   }
    #'   Class-specific entries take precedence over \code{all_inference} for
    #'   constructor arguments.  Example:
    #'   \preformatted{inference_params = list(
    #'     all_inference          = list(delta = 0, na.rm = TRUE,
    #'                                   harden = TRUE),
    #'     InferenceContinMultOLS = list(harden = FALSE)
    #'   )}
    #'   Here every class gets \code{harden = TRUE} from \code{all_inference},
    #'   except \code{InferenceContinMultOLS} which uses \code{harden = FALSE}.
    #'
    #' @param inf_types Character vector controlling which inference outputs are
    #'   computed and recorded.  \code{NULL} (default) runs all eight types.
    #'   Valid elements:
    #'   \describe{
    #'     \item{\code{"asymp_ci"}}{Asymptotic (Wald) confidence interval.}
    #'     \item{\code{"asymp_pval"}}{Asymptotic two-sided p-value.}
    #'     \item{\code{"exact_ci"}}{Exact confidence interval
    #'       (\code{InferenceExact} subclasses only; others warned and skipped).}
    #'     \item{\code{"exact_pval"}}{Exact two-sided p-value (same caveat).}
    #'     \item{\code{"boot_ci"}}{Bootstrap percentile confidence interval.}
    #'     \item{\code{"boot_pval"}}{Bootstrap two-sided p-value.}
    #'     \item{\code{"rand_ci"}}{Randomisation-based confidence interval.}
    #'     \item{\code{"rand_pval"}}{Randomisation two-sided p-value.}
    #'   }
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
      sd_noise          = 1,
      prob_censoring    = 0.2,
      verbose           = FALSE,
      design_params     = list(),
      inference_params  = list(),
      inf_types         = NULL
    ) {
      valid_rt = c("continuous", "incidence", "proportion",
                   "count", "survival", "ordinal")
      if (!response_type %in% valid_rt)
        stop("response_type must be one of: ", paste(valid_rt, collapse = ", "))
      if (!data_type %in% c("linear", "nonlinear"))
        stop("data_type must be 'linear' or 'nonlinear'")
      if (data_type == "nonlinear" && p < 5L)
        stop("nonlinear Friedman function requires p >= 5")
      if (!is.list(design_params))
        stop("design_params must be a named list")
      if (!is.list(inference_params))
        stop("inference_params must be a named list")

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
      private$sd_noise         = sd_noise
      private$prob_censoring   = prob_censoring
      private$verbose          = verbose
      private$design_params    = design_params
      private$inference_params = inference_params
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

      private$raw_results = list()
      private$has_run     = FALSE
    },

    # ── run() ─────────────────────────────────────────────────────────────────
    run = function() {
      private$raw_results       = list()
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

        for (design_gen in private$design_classes) {
          design_name  = design_gen$classname
          design_extra = private$design_params[[design_name]]
          if (is.null(design_extra)) design_extra = list()

          des_obj = tryCatch(
            private$.build_design(design_gen, X, y_base, design_extra),
            error = function(e) {
              if (private$verbose)
                message(sprintf("    [SKIP design %s]: %s", design_name, e$message))
              NULL
            }
          )
          if (is.null(des_obj)) next

          for (inf_gen in private$inference_classes) {
            inf_name = inf_gen$classname

            # Merge all_inference (lower priority) with class-specific (higher)
            all_p   = private$inference_params[["all_inference"]]
            if (is.null(all_p)) all_p = list()
            class_p = private$inference_params[[inf_name]]
            if (is.null(class_p)) class_p = list()
            merged  = modifyList(all_p, class_p)

            # Route merged params to their call-site destinations
            routed = private$.route_inference_args(merged)

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
              private$.record(rep, design_name, inf_name, "asymp", est, ci_a, pval_a)
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
                private$.record(rep, design_name, inf_name, "exact", est, ci_e, pval_e)
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
              private$.record(rep, design_name, inf_name, "boot", est, ci_b, pval_b)
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
              private$.record(rep, design_name, inf_name, "rand", est, ci_r, pval_r)
            }
          }
        }
      }

      private$has_run = TRUE
      invisible(self)
    },

    # ── get_results() ─────────────────────────────────────────────────────────
    # Returns a data.table with one row per (rep, design, inference, method).
    get_results = function() {
      if (!private$has_run) stop("Call $run() first.")
      if (length(private$raw_results) == 0L)
        return(data.table::data.table(
          rep = integer(), design = character(), inference = character(),
          method = character(), estimate = numeric(),
          ci_lo = numeric(), ci_hi = numeric(), pval = numeric()
        ))
      data.table::rbindlist(private$raw_results)
    },

    # ── summarize() ───────────────────────────────────────────────────────────
    # Returns a data.table with MSE, and conditionally coverage / power,
    # per (design, inference, method).
    summarize = function() {
      if (!private$has_run) stop("Call $run() first.")
      dt = self$get_results()
      if (nrow(dt) == 0L) { message("No results."); return(invisible(NULL)) }

      betaT      = private$betaT
      alpha      = private$alpha
      report_cov = any(grepl("_ci$",   private$inf_types))
      report_pow = any(grepl("_pval$", private$inf_types))

      result = dt[, {
        est_ok = estimate[is.finite(estimate)]
        ci_ok  = .SD[is.finite(ci_lo) & is.finite(ci_hi)]
        pv_ok  = pval[is.finite(pval)]
        row = list(
          MSE   = if (length(est_ok)) round(mean((est_ok - betaT)^2), 5L) else NA_real_,
          n_est = length(est_ok)
        )
        if (report_cov) {
          row$coverage = if (nrow(ci_ok)) round(mean(ci_ok$ci_lo <= betaT & betaT <= ci_ok$ci_hi), 3L) else NA_real_
          row$n_cov    = nrow(ci_ok)
        }
        if (report_pow) {
          row$power = if (length(pv_ok)) round(mean(pv_ok < alpha), 3L) else NA_real_
          row$n_pow = length(pv_ok)
        }
        row
      }, by = .(design, inference, method),
         .SDcols = c("estimate", "ci_lo", "ci_hi", "pval")]

      result[order(design, inference, method)]
    },

    # ── print() ───────────────────────────────────────────────────────────────
    print = function() {
      cat("SimulationFramework\n")
      cat("  response_type :", private$response_type, "\n")
      cat("  data_type     :", private$data_type, "\n")
      cat("  n / p         :", private$n, "/", private$p, "\n")
      cat("  Nrep / betaT  :", private$Nrep, "/", private$betaT, "\n")
      cat("  alpha / B_boot / r_rand :",
          private$alpha, "/", private$B_boot, "/", private$r_rand, "\n")
      cat("  inf_types     :", paste(private$inf_types, collapse = ", "), "\n")
      if (length(private$design_params))
        cat("  design_params :", paste(names(private$design_params), collapse = ", "), "\n")
      if (length(private$inference_params))
        cat("  inference_params:", paste(names(private$inference_params), collapse = ", "), "\n")
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
    sd_noise         = NULL,
    prob_censoring   = NULL,
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
    raw_results      = NULL,
    has_run          = FALSE,
    exact_warned_classes = NULL,

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

    # ── Data generation ───────────────────────────────────────────────────────

    .generate_data = function() {
      n = private$n; p = private$p
      X_mat = matrix(stats::runif(n * p), nrow = n, ncol = p)
      colnames(X_mat) = paste0("x", seq_len(p))
      X = as.data.frame(X_mat)

      if (private$data_type == "linear") {
        beta_x = seq(1, -1, length.out = p)
        y_cont = as.numeric(X_mat %*% beta_x)
      } else {
        # Friedman (1991) function — first five covariates, x in [0,1]
        y_cont = 10 * sin(pi * X_mat[, 1L] * X_mat[, 2L]) +
                 20 * (X_mat[, 3L] - 0.5)^2 +
                 10 *  X_mat[, 4L] +
                  5 *  X_mat[, 5L]
      }

      list(X = X, y_base = private$.y_cont_to_response(y_cont))
    },

    # Transform continuous signal to the appropriate response-type scale.
    .y_cont_to_response = function(y_cont) {
      rt  = private$response_type
      y_s = as.numeric(scale(y_cont))
      switch(rt,
        continuous = y_s,
        incidence  = stats::plogis(y_s),
        proportion = { y_sh = y_cont - min(y_cont) + 1e-6
                       y_sh / (max(y_sh) + 1e-6) },
        count      = pmax(1L, round(y_cont - min(y_cont) + 1)),
        survival   = pmax(0.1, y_s - min(y_s) + 0.1),
        ordinal    = as.integer(cut(
          y_cont,
          breaks = unique(quantile(y_cont, probs = seq(0, 1, length.out = 5))),
          include.lowest = TRUE)),
        stop("Unknown response_type: ", rt)
      )
    },

    # Apply treatment effect + noise to one subject; return list(y, dead).
    .apply_one = function(y_base_i, w_i) {
      rt  = private$response_type
      bt  = if (w_i == 1L) private$betaT else 0
      eps = stats::rnorm(1L, 0, private$sd_noise)
      y = switch(rt,
        continuous = y_base_i + bt + eps,
        incidence  = {
          p_b = if (is.finite(y_base_i) && y_base_i >= 0 && y_base_i <= 1)
                  y_base_i else stats::plogis(y_base_i)
          p_b = pmin(0.95, pmax(0.05, p_b))
          as.numeric(stats::rbinom(1L, 1L,
            stats::plogis(stats::qlogis(p_b) + bt + eps)))
        },
        proportion = pmin(1 - 1e-9, pmax(1e-9, y_base_i + bt + eps)),
        count = {
          lam = pmax(.Machine$double.eps, y_base_i * exp(bt + eps))
          as.numeric(stats::rpois(1L, lambda = lam))
        },
        survival = pmax(.Machine$double.eps, y_base_i * exp(bt + eps)),
        ordinal  = as.integer(max(1L, round(y_base_i + bt + eps))),
        stop("Unknown response_type: ", rt)
      )
      dead = 1L
      if (rt == "survival" && stats::runif(1L) < private$prob_censoring) {
        y    = stats::runif(1L, 0, y)
        dead = 0L
      }
      list(y = y, dead = dead)
    },

    # Instantiate design and run the full experiment (assign + observe all n).
    .build_design = function(design_gen, X, y_base, design_extra) {
      n       = private$n
      des_obj = do.call(design_gen$new, c(
        list(response_type = private$response_type, n = n),
        design_extra
      ))

      if (inherits(des_obj, "DesignSeqOneByOne")) {
        for (t in seq_len(n)) {
          w_t = des_obj$add_one_subject_to_experiment_and_assign(
            X[t, , drop = FALSE])
          out = private$.apply_one(y_base[t], w_t)
          des_obj$add_one_subject_response(t, out$y, out$dead)
        }
      } else {
        des_obj$add_all_subjects_to_experiment(X)
        des_obj$assign_w_to_all_subjects()
        w = des_obj$get_w()
        for (t in seq_len(n)) {
          out = private$.apply_one(y_base[t], w[t])
          des_obj$add_one_subject_response(t, out$y, out$dead)
        }
      }
      des_obj
    },

    # Append one data.table row to raw_results.
    .record = function(rep, design, inf_name, method, est, ci, pval) {
      ci2 = if (length(ci) >= 2L) as.numeric(ci[1:2]) else c(NA_real_, NA_real_)
      if (all(is.finite(ci2)) && ci2[1L] > ci2[2L]) ci2 = rev(ci2)
      private$raw_results[[length(private$raw_results) + 1L]] = data.table::data.table(
        rep       = as.integer(rep),
        design    = design,
        inference = inf_name,
        method    = method,
        estimate  = if (is.null(est) || !is.finite(est)) NA_real_ else est,
        ci_lo     = ci2[1L],
        ci_hi     = ci2[2L],
        pval      = if (is.null(pval) || length(pval) == 0L || !is.finite(pval[1L]))
                      NA_real_ else as.numeric(pval[1L])
      )
    },

    # ── Defaults ──────────────────────────────────────────────────────────────

    .default_design_classes = function() {
      list(
        DesignSeqOneByOneBernoulli,
        DesignSeqOneByOneiBCRD,
        DesignSeqOneByOneEfron,
        DesignSeqOneByOneAtkinson,
        DesignSeqOneByOneUrn,
        DesignSeqOneByOneKK21,
        DesignSeqOneByOneKK14,
        FixedDesignBernoulli,
        FixedDesigniBCRD,
        FixedDesignBinaryMatch,
        FixedDesignGreedy,
        FixedDesignMatchingGreedyPairSwitching,
        FixedDesignRerandomization,
        FixedDesignDOptimal,
        FixedDesignAOptimal,
        FixedDesignOptimalBlocks
      )
    },

    .default_inference_classes = function() {
      rt   = private$response_type
      univ = if (rt != "survival") list(InferenceAllSimpleMeanDiff) else list()

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
