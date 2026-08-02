#' Mixin for Off-Optimum Likelihood Evaluation
#'
#' A Pattern-1 mixin (plain list with \code{$public} and \code{$private} slots)
#' providing evaluation of the negative log-likelihood and the information
#' matrix at an arbitrary parameter vector \code{theta}, not only at a
#' completed fit's optimum. The core likelihood-test machinery in
#' \code{InferenceAsympLik} only ever evaluates \code{neg_loglik}/
#' \code{observed_information}/\code{fisher_information} at a \code{fit}
#' object returned by the native optimizer, because that is all Wald/score/
#' gradient/LR testing needs. Off-optimum evaluation is a distinct, strictly
#' additional capability needed by anything that must walk the likelihood
#' surface away from the MLE -- e.g. a Firth/Jeffreys penalty
#' \code{log|I(theta)|}, profile likelihoods, or derivative-free optimization
#' of a penalized objective.
#'
#' This is opt-in, not part of the standard \code{InferenceAsympLik}
#' composition: most concrete classes' \code{get_likelihood_test_spec()}
#' closures are hard-wired to extract \code{theta} from a \code{fit} object
#' (e.g. \code{c(fit$b, log(fit$theta_hat))}) rather than accepting an
#' arbitrary \code{theta} directly, even when the underlying native evaluator
#' they call (e.g. \code{get_negbin_regression_hessian_cpp(X, y, theta)}) is
#' already a genuine function of \code{theta}. Splicing in this mixin does not
#' by itself grant off-optimum support -- a class must additionally add two
#' fields to the list returned by \code{get_likelihood_test_spec()}:
#'
#' \itemize{
#'   \item \code{neg_loglik_at(theta)}: numeric scalar, the negative
#'     log-likelihood evaluated at \code{theta}.
#'   \item \code{information_at(theta, source = c("observed", "fisher"))}:
#'     p x p matrix, the requested information matrix evaluated at
#'     \code{theta}.
#' }
#'
#' Classes that omit these two spec fields simply do not support off-optimum
#' evaluation; \code{supports_off_optimum_likelihood_eval()} reports FALSE and
#' the evaluator methods below error clearly rather than silently evaluating
#' at the wrong point (e.g. falling back to the fitted MLE).
#'
#' Consumers must provide private methods \code{get_likelihood_test_spec()}
#' and \code{get_default_information_source()} (both already required by
#' \code{InferenceMixinInformationMatrix}, which every \code{InferenceAsympLik}
#' subclass already carries).
#'
#' Splice into a class with
#' \code{private = c(InferenceMixinOffOptimumLikelihoodEval$private, list(...))}.
#'
#' @noRd
InferenceMixinOffOptimumLikelihoodEval = list(
	public = list(),
	private = list(
		is_a_off_optimum_likelihood_eval = function() TRUE,
		supports_off_optimum_likelihood_eval = function(spec = NULL){
			if (is.null(spec)) {
				spec = private$get_likelihood_test_spec()
			}
			!is.null(spec) && is.function(spec$neg_loglik_at) && is.function(spec$information_at)
		},

		evaluate_neg_loglik_at = function(theta, spec = NULL){
			if (is.null(spec)) {
				spec = private$get_likelihood_test_spec()
			}
			if (is.null(spec) || !is.function(spec$neg_loglik_at)) {
				stop(class(self)[1], " does not support off-optimum negative log-likelihood evaluation.", call. = FALSE)
			}
			spec$neg_loglik_at(theta)
		},

		evaluate_information_at = function(theta, spec = NULL, source = "auto"){
			if (is.null(spec)) {
				spec = private$get_likelihood_test_spec()
			}
			if (is.null(spec) || !is.function(spec$information_at)) {
				stop(class(self)[1], " does not support off-optimum information-matrix evaluation.", call. = FALSE)
			}
			if (identical(source, "auto")) {
				source = private$get_default_information_source()
				if (identical(source, "legacy")) {
					source = "observed"
				}
			}
			spec$information_at(theta, source = source)
		},

		# Convenience combining the two evaluators into a single scalar
		# Firth/Jeffreys-penalized negative log-likelihood,
		# neg_loglik(theta) + 0.5 * log|I(theta)|. Returns Inf (rather than
		# erroring) when the information matrix is singular or otherwise
		# produces a non-finite log-determinant, so this can be handed
		# directly to an optimizer -- including a derivative-free one --
		# without a separate feasibility check at every candidate theta.
		evaluate_penalized_neg_loglik_at = function(theta, spec = NULL, source = "auto"){
			if (is.null(spec)) {
				spec = private$get_likelihood_test_spec()
			}
			neg_loglik = private$evaluate_neg_loglik_at(theta, spec = spec)
			information = private$evaluate_information_at(theta, spec = spec, source = source)
			log_det = tryCatch(
				as.numeric(determinant(as.matrix(information), logarithm = TRUE)$modulus),
				error = function(e) NA_real_
			)
			if (!is.finite(log_det)) {
				return(Inf)
			}
			neg_loglik + 0.5 * log_det
		}
	)
)
