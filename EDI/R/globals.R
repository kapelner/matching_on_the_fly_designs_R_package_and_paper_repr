# Internal environment for the EDI package to store global state
edi_env = new.env(parent = emptyenv())

# Closure to encapsulate the internal assertion override flag (used by SimulationFramework)
.assert_manager = (function() {
  internal_run_asserts = TRUE
  list(
    toggle = function(on = TRUE) {
      internal_run_asserts <<- isTRUE(on)
      invisible(internal_run_asserts)
    },
    should_run = function() {
      internal_run_asserts && isTRUE(getOption("edi.run_asserts", TRUE))
    }
  )
})()

#' Toggle the execution of assertions throughout the package
#' 
#' @description
#' This function enables or disables the internal input validation checks (assertions)
#' by setting the \code{options(edi.run_asserts = ...)} value.
#' Disabling assertions can provide a significant performance boost in heavy
#' simulations (often 10x-20x speedup), but it removes the safety rails that
#' prevent invalid data from reaching the internal algorithms.
#' 
#' \strong{Warning:} If assertions are disabled, passing malformed or invalid
#' data to package functions may result in cryptic R errors, incorrect
#' statistical results, or even hard system crashes (SEGFAULTs) at the C++ layer.
#' Only disable assertions if you are certain your data is pre-validated and
#' follows the package requirements exactly.
#' 
#' @param on Logical scalar. If TRUE (default), assertions are executed. If FALSE, they are skipped.
#' @keywords internal
#' @export
toggle_asserts = function(on = TRUE) {
  options(edi.run_asserts = isTRUE(on))
  invisible(isTRUE(on))
}

# private method
should_run_asserts = .assert_manager$should_run

#' Set the number of cores for parallelization
#' 
#' This function initializes a persistent parallel cluster (either a fork cluster
#' on Unix-like systems or a mirai cluster on others) to be used by all Design
#' and Inference objects. This avoids the overhead of creating clusters
#' repeatedly.
#'
#' @details
#' \code{set_num_cores()} sets a global upper bound for parallel work. It does
#' not guarantee that every inference routine will use all requested workers.
#' EDI's inference dispatcher applies a blocklist-first heuristic informed by
#' package benchmarks: workloads that have shown consistent multicore slowdowns
#' are forced to run serially, while the remaining workloads are allowed to use
#' their method-specific warmup heuristics and native thread caps.
#'
#' The default forced-serial blocklist covers incidence randomization confidence
#' intervals, bootstrap for non-regression KK Wilcoxon inference, bootstrap for
#' non-KK survival procedures, and bootstrap for incidence procedures. Do not
#' expect a universal "more cores is faster" rule.
#'
#' If you want to change the default policy, use
#' \code{set_parallel_dispatch_policy()}.
#' 
#' @param num_cores Integer number of worker processes to make available.
#' @param force_mirai If \code{TRUE}, forces the use of the \code{mirai} package
#'   even on systems where forking is available.
#' 
#' @return Invisible \code{NULL}.
#' 
#' @export
set_num_cores = function(num_cores, force_mirai = FALSE) {
  checkmate::assertCount(num_cores, positive = TRUE)
  checkmate::assertFlag(force_mirai)
  
  # Clear any existing clusters first
  unset_num_cores()

  # Set package/native thread budgets before workers are created so forked
  # children inherit the intended global thread settings.
  set_package_threads(num_cores)

  if (as.integer(num_cores) <= 1L) {
    return(invisible(NULL))
  }
  
  if (force_mirai || .Platform$OS.type != "unix") {
    if (!check_package_installed("mirai")) {
      stop("The 'mirai' package is required for parallelization on this system or when force_mirai = TRUE. Please install it.")
    }
    edi_env$mirai_has_been_used = TRUE
    # Initialize mirai daemons
    mirai::daemons(num_cores)
    edi_env$global_mirai_num_cores = num_cores
  } else {
    if (isTRUE(edi_env$mirai_has_been_used)) {
      stop(
        "Cannot switch from mirai-backed parallelism to fork-based parallelism in the same R session. ",
        "Restart R or keep using force_mirai = TRUE. This avoids the nng is not fork-reentrant safe panic."
      )
    }
    # Unix-like system, use forking
    edi_env$global_fork_cluster = parallel::makeForkCluster(num_cores)
  }
  
  invisible(NULL)
}


#' Unset the number of cores and stop parallel clusters
#' 
#' This function stops any global fork or mirai clusters stored in the package
#' environment and resets the core count to serial execution.
#' 
#' @return Invisible \code{NULL}.
#' 
#' @export
unset_num_cores = function() {
  # Handle fork cluster
  if (!is.null(edi_env$global_fork_cluster)) {
    cl = edi_env$global_fork_cluster
    edi_env$global_fork_cluster = NULL
    tryCatch(parallel::stopCluster(cl), error = function(e) invisible(NULL))
  }
  
  # Handle mirai cluster
  if (!is.null(edi_env$global_mirai_num_cores)) {
    if (check_package_installed("mirai")) {
      tryCatch(mirai::daemons(0), error = function(e) invisible(NULL)) # Stop daemons
    }
    edi_env$global_mirai_num_cores = NULL
  }
  
  # Reset package threads to 1
  set_package_threads(1L)
  
  invisible(NULL)
}


#' Get the default parallel dispatch policy
#'
#' Returns EDI's built-in blocklist-first dispatch policy. This is the policy
#' used when the user has not installed a custom dispatch override.
#'
#' Do not expect a universal "more cores is faster" rule.
#'
#' @return A named list describing the built-in dispatch policy.
#'
#' @export
get_parallel_dispatch_policy = function() {
  list(
    bootstrap = list(
      serial_inference_class_patterns = c(
        "^InferenceIncid",
        "^InferenceSurvival(?!.*KK)",
        "^InferenceAllKKWilcoxIVWC$"
      ),
      serial_response_types = c("incidence")
    ),
    rand_ci = list(
      serial_inference_class_patterns = c("^InferenceIncid"),
      serial_response_types = c("incidence")
    )
  )
}

edi_env$parallel_dispatch_policy_config = get_parallel_dispatch_policy()

#' Get the default bootstrap dispatch policy
#'
#' Returns EDI's built-in policy for choosing a default bootstrap type. Each
#' inference class can override the standard \code{"bca"} default via a regular
#' expression pattern.
#'
#' @details
#' The built-in overrides are empirical. Inference classes with repeated
#' bootstrap failures in the comprehensive test suite, especially classes whose
#' bootstrap-standard-error or related higher-order bootstrap paths are
#' numerically fragile, are dispatched to \code{"percentile"} by default. This
#' keeps the default bootstrap path estimate-only, faster, and less prone to
#' NA/NaN/Inf failures.
#'
#' @return A named list describing the default bootstrap type configuration.
#' @export
get_bootstrap_dispatch_policy = function() {
  list(
    default_type = "bca",
    inference_class_overrides = c(
      "^InferenceContinMultLin$" = "percentile",
      "^InferenceIncid(Univ|Multi)GCompRisk(Diff|Ratio)$" = "percentile",
      "^InferenceIncid(Univ|Multi)KKGCompRisk(Diff|Ratio)$" = "percentile",
      "^InferenceProp(Uni|Multi)GCompMeanDiff$" = "percentile",
      "^InferenceSurvival(Uni|Multi)DepCensTransformRegr$" = "percentile",
      "^InferenceSurvival(Univ|Multi)KKRankRegrIVWC$" = "percentile",
      "^InferenceIncidMultiKKClogitIVWC$" = "percentile",
      "^InferenceIncidMultiKKClogitPlusGLMMIVWC$" = "percentile",
      "^InferenceIncidMultiKKClogitCombinedLikelihood$" = "percentile",
      "^InferenceIncidMultiKKClogitPlusGLMMCombinedLikelihood$" = "percentile",
      "^InferenceOrdinalUnivKKCondPropOddsRegr$" = "percentile",
      "^InferenceOrdinalMultiAdjCatLogitRegr$" = "percentile",
      "^InferenceSurvival(Univ|Multi)KKStratCoxCombinedLikelihood$" = "percentile",
      "^InferenceCountMultiPoissonRegr$" = "percentile",
      "^InferenceCountMultiQuasiPoissonRegr$" = "percentile",
      "^InferencePropMultiKKQuantileRegrCombinedLikelihood$" = "percentile",
      "^InferenceSurvivalMultiDepCensTransformRegr$" = "percentile",
      "^InferencePropMultiKKQuantileRegrIVWC$" = "percentile",
      "^InferenceOrdinalMultiCumulProbitRegr$" = "percentile",
      "^InferenceOrdinalMultiPartialProportionalOddsRegr$" = "percentile",
      "^InferencePropMultiZeroOneInflatedBetaRegr$" = "percentile",
      "^InferencePropMultiFractionalLogit$" = "percentile",
      "^InferenceCountMultiHurdleNegBinRegr$" = "percentile",
      "^InferenceContinMultiRobustRegr$" = "percentile"
    ),
    design_class_overrides = list(
      FixedDesignBlockedCluster = c(
        "^InferenceContinMultiRobustRegr$" = "percentile",
        "^InferenceContinMultLin$" = "percentile",
        "^InferenceContinMultOLS$" = "percentile"
      )
    )
  )
}

edi_env$bootstrap_dispatch_policy_config = get_bootstrap_dispatch_policy()

#' Update the parallel dispatch policy
#'
#' EDI uses an empirical, blocklist-first dispatch policy to decide when an
#' inference routine should be forced serial even if multiple cores are
#' available. This function lets the user update that default policy without
#' editing package internals.
#'
#' @details
#' The policy can be updated in two ways:
#' \itemize{
#'   \item Pass a named list to merge with the built-in default policy
#'   configuration. Supported top-level keys are \code{bootstrap} and
#'   \code{rand_ci}, and each key may contain
#'   \code{serial_inference_class_patterns} and
#'   \code{serial_response_types}.
#'   \item Pass a custom function with signature
#'   \code{function(inference_class, response_type, operation)} that returns a
#'   list with at least \code{force_serial} and \code{reason}.
#' }
#'
#' Use \code{reset = TRUE} to restore the built-in default policy. Do not expect
#' a universal "more cores is faster" rule.
#'
#' @param policy Either \code{NULL}, a named list of policy overrides, or a
#'   custom function. If \code{NULL} and \code{reset = FALSE}, the current
#'   effective policy configuration is returned invisibly.
#' @param reset If \code{TRUE}, restore the built-in default policy and remove
#'   any custom function override.
#'
#' @return Invisible \code{NULL} or the current policy configuration.
#'
#' @export
set_parallel_dispatch_policy = function(policy = NULL, reset = FALSE) {
  checkmate::assertFlag(reset)

  if (isTRUE(reset)) {
    edi_env$parallel_dispatch_policy_config = get_parallel_dispatch_policy()
    edi_env$parallel_dispatch_policy_override = NULL
    return(invisible(edi_env$parallel_dispatch_policy_config))
  }

  if (is.null(policy)) {
    return(invisible(edi_env$parallel_dispatch_policy_config))
  }

  if (is.function(policy)) {
    edi_env$parallel_dispatch_policy_override = policy
    return(invisible(NULL))
  }

  checkmate::assertList(policy, names = "named")
  current_config = edi_env$parallel_dispatch_policy_config
  for (nm in names(policy)) {
    if (!nm %in% names(current_config)) {
      stop("Unknown policy section: ", nm, call. = FALSE)
    }
    if (!is.list(policy[[nm]])) {
      stop("Policy section '", nm, "' must be a list.", call. = FALSE)
    }
    current_config[[nm]] = utils::modifyList(current_config[[nm]], policy[[nm]])
  }
  edi_env$parallel_dispatch_policy_config = current_config
  edi_env$parallel_dispatch_policy_override = NULL
  invisible(NULL)
}

# Internal helper to get the global fork cluster
get_global_fork_cluster = function() {
  edi_env$global_fork_cluster
}

# Internal helper to get the global mirai core count
get_global_mirai_cores = function() {
  edi_env$global_mirai_num_cores
}

# Internal helper to get the current core count budget
get_num_cores = function() {
  if (!is.null(edi_env$num_cores_override)) return(edi_env$num_cores_override)
  cl = get_global_fork_cluster()
  if (!is.null(cl)) return(length(cl))
  mirai_cores = get_global_mirai_cores()
  if (!is.null(mirai_cores)) return(mirai_cores)
  1L
}

# Internal helper for empirical parallel dispatch policy.
# This is intentionally conservative and only forces serial execution for
# inference families that have shown repeat multicore regressions in the package
# benchmark suite. Everything else remains eligible for the usual warmup-based
# parallel heuristics. Do not expect a universal "more cores is faster" rule.
edi_parallel_dispatch_policy = function(inference_class, response_type, operation) {
  if (is.function(edi_env$parallel_dispatch_policy_override)) {
    return(edi_env$parallel_dispatch_policy_override(inference_class, response_type, operation))
  }

  inference_class = as.character(inference_class[[1]])
  response_type = as.character(response_type[[1]])
  operation = as.character(operation[[1]])

  config = edi_env$parallel_dispatch_policy_config
  matches_any = function(value, patterns) {
    if (is.null(patterns) || length(patterns) == 0L) return(FALSE)
    any(vapply(patterns, function(pattern) grepl(pattern, value, perl = TRUE), logical(1)))
  }

  reason = NULL
  if (identical(operation, "bootstrap")) {
    op_cfg = config$bootstrap
    if (matches_any(inference_class, op_cfg$serial_inference_class_patterns) ||
        matches_any(response_type, op_cfg$serial_response_types)) {
      reason = "bootstrap is forced serial by benchmark policy"
    }
  } else if (identical(operation, "rand_ci")) {
    op_cfg = config$rand_ci
    if (matches_any(inference_class, op_cfg$serial_inference_class_patterns) ||
        matches_any(response_type, op_cfg$serial_response_types)) {
      reason = "randomization confidence intervals are forced serial by benchmark policy"
    }
  }

  list(
    force_serial = !is.null(reason),
    reason = reason,
    inference_class = inference_class,
    response_type = response_type,
    operation = operation
  )
}

edi_bootstrap_dispatch_policy = function(inference_class, object = NULL) {
  config = edi_env$bootstrap_dispatch_policy_config
  inference_class = as.character(inference_class[[1]])
  overrides = config$inference_class_overrides
  design_overrides = config$design_class_overrides
  default_type = config$default_type
  if (is.null(default_type)) default_type = "bca"
  if (!is.null(object) && !is.null(design_overrides) && length(design_overrides) > 0L) {
    des_obj = tryCatch(object$.__enclos_env__$private$des_obj, error = function(e) NULL)
    if (!is.null(des_obj)) {
      for (design_class in names(design_overrides)) {
        if (is(des_obj, design_class)) {
          design_map = design_overrides[[design_class]]
          for (pattern in names(design_map)) {
            if (is.na(pattern) || pattern == "") next
            if (grepl(pattern, inference_class, perl = TRUE)) {
              return(tolower(design_map[[pattern]]))
            }
          }
        }
      }
    }
  }
  if (!is.null(overrides) && length(overrides) > 0L) {
    for (pattern in names(overrides)) {
      if (is.na(pattern) || pattern == "") next
      if (grepl(pattern, inference_class, perl = TRUE)) {
        return(tolower(overrides[[pattern]]))
      }
    }
  }
  tolower(default_type)
}

# Internal helper
set_package_threads = function(num_cores) {
  # Ensure it's an integer
  num_cores = as.integer(num_cores)

  # Check if we are already at the target to avoid slow Sys.setenv and other calls
  if (!is.null(options(".edi_last_set_threads")[[1]]) && options(".edi_last_set_threads")[[1]] == num_cores) {
    return(invisible(NULL))
  }

  # R packages with global thread setters
	  if (check_package_installed("data.table")) {
	    data.table::setDTthreads(num_cores)
	  }
	  if (check_package_installed("fixest")) {
	    suppressWarnings(try(fixest::setFixest_nthreads(num_cores), silent = TRUE))
	  }

  # Environment variables for OpenMP and BLAS/LAPACK
  # This helps prevent thread explosion in child processes
  # that call multi-threaded native libraries.
  # Sys.setenv is relatively slow, so we only do it if needed.
  Sys.setenv(OMP_NUM_THREADS =    num_cores)
  Sys.setenv(MKL_NUM_THREADS =    num_cores)
  Sys.setenv(OPENBLAS_NUM_THREADS =   num_cores)
  Sys.setenv(VECLIB_MAXIMUM_THREADS = num_cores)
  Sys.setenv(NUMEXPR_NUM_THREADS =    num_cores)
  options(".edi_last_set_threads" =   num_cores)
  invisible(NULL)
}
