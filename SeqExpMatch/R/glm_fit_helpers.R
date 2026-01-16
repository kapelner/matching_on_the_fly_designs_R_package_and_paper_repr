#' Fast Ordinary Least Squares Regression
#'
#' This function performs a fast Ordinary Least Squares (OLS) regression using Eigen.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric vector of the response variable.
#'
#' @return A list containing the following components:
#' \item{b}{A numeric vector of the estimated coefficients.}
#' \item{XtX}{The matrix X'X.}
#'
#' @name fast_ols_cpp
#' @rdname fast_ols_cpp
#' @export
NULL


fast_logistic_regression = function(Xmm, y){
  mod = fastLogisticRegressionWrap::fast_logistic_regression(
      Xmm = Xmm,
      ybin = as.numeric(y)
  )
  list(b = as.vector(mod$coefficients))
}

binomial_link_cache = binomial()

fast_logistic_regression_with_var = function(Xmm, y){
  # Compute inference on the second coefficient (treatment effect)
  mod = fastLogisticRegressionWrap::fast_logistic_regression(
      Xmm = Xmm,
      ybin = as.numeric(y),
      do_inference_on_var = 2
  )
  list(
    b = as.vector(mod$coefficients),
    ssq_b_2 = mod$se[2]^2  # Square the standard error to get variance
  )
}


fast_beta_regression <- function(Xmm, y,
                             start_phi = 10,
                             bounds_logphi = c(log(1e-3), log(1e4)),
                             control = stats::glm.control(epsilon=1e-8, maxit=100)) {
   
  weights <- rep(1, nrow(Xmm))

  # Use C++ logistic regression to find beta (quasi-likelihood, independent of phi)
  # This replaces the repetitive glm.fit calls inside the optimization loop
  mod_log = fast_logistic_regression_cpp(as.matrix(Xmm), as.numeric(y), maxit = control$maxit, tol = control$epsilon)
  b = mod_log$b
  mu = 1 / (1 + exp(-drop(Xmm %*% b)))
  # guard against numerical drift to 0/1
  mu <- pmin(pmax(mu, 1e-12), 1 - 1e-12)

  # objective: negative log-likelihood profiled over beta (beta is fixed now)
  obj <- function(logphi) {
    phi = exp(logphi)
    # return NEGATIVE log-likelihood
    -beta_loglik_cpp(y, mu, phi, wt = weights)
  }

  opt <- stats::nlminb(start = log(start_phi), objective = obj,
                lower = bounds_logphi[1], upper = bounds_logphi[2])

  phi_hat <- as.numeric(exp(opt$par))
  
  # Construct weights for the object matching beta_family/glm.fit behavior
  # weights = (1+phi) * mu * (1-mu)
  w_final = (1 + phi_hat) * mod_log$w

  out <- list(
    b = b,
    phi = phi_hat,
    converged_inner = TRUE,
    w = w_final
  )
  class(out) <- "beta_glm_profile"
  out
}

fast_beta_regression_with_var <- function(Xmm, y,
                             start_phi = 10,
                             bounds_logphi = c(log(1e-3), log(1e4)),
                             control = stats::glm.control(epsilon=1e-8, maxit=100)) {
    fit_hat = fast_beta_regression(Xmm, y, start_phi = start_phi, bounds_logphi = bounds_logphi, control = control) 
    XtWX <- crossprod(sqrt(fit_hat$w) * Xmm)
    fit_hat$ssq_b_2 = eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(XtWX, 2)
    fit_hat
}

fast_coxph_regression = function(Xmm, y, dead){
	mod = glmnet::glmnet(Xmm, survival::Surv(y, dead), family = "cox", lambda = 0)
	list(b = stats::coef(mod))
}

fast_negbin_regression <- function(Xmm, y) {
  list(
    b = fast_neg_bin_with_censoring_cpp(
      X = Xmm,
      y = as.integer(y),
      dead = rep(1L, length(y))
    )$b
  )  
}

fast_negbin_regression_with_var <- function(Xmm, y) {
  res <- fast_neg_bin_with_censoring_with_sd_cpp(
    X = Xmm,
    y = as.integer(y),
    dead = rep(1L, length(y))
  )
  
  # Extract beta_vcov
  p_beta <- length(res$b)
  beta_vcov_matrix <- as.matrix(solve(res$hess_fisher_info_matrix)[1:p_beta, 1:p_beta, drop = FALSE])
  
  list(
    b = res$b,
    ssq_b_2 = beta_vcov_matrix[2, 2]
  )
}




































# beta_family <- function(link = "logit", phi = 10) {
#   linkobj <- stats::make.link(link)
  
#   variance <- function(mu) {
#     mu * (1 - mu) / (1 + phi)
#   }
  
#   dev.resids <- function(y, mu, wt) {
#     # negative twice log-likelihood contribution
# #    2 * wt * (lbeta(mu * phi, (1 - mu) * phi) -
# #              (mu * phi - 1) * log(y) -
# #              ((1 - mu) * phi - 1) * log(1 - y))
# 	beta_dev_resids_cpp(y, mu, phi, wt)
#   }
  
#   aic <- function(y, n, mu, wt, dev) {
#     # -2*logLik + 2*edf
# #    -2 * sum(wt * (
# #      lgamma(phi) - lgamma(mu * phi) - lgamma((1 - mu) * phi) +
# #      (mu * phi - 1) * log(y) + ((1 - mu) * phi - 1) * log(1 - y)
# #    )) + 2 * (length(mu) + 1)
#     beta_aic_cpp(y, mu, phi, wt)
#   }
  
#   mu.eta <- linkobj$mu.eta
  
#   structure(
#     list(
#       family = "Beta",
#       link = linkobj$name,
#       linkfun = linkobj$linkfun,
#       linkinv = linkobj$linkinv,
#       variance = variance,
#       dev.resids = dev.resids,
#       aic = aic,
#       mu.eta = mu.eta,
#       initialize = expression({
#         if (any(y <= 0 | y >= 1))
#           stop("y values must be in (0,1) for beta regression")
#         mustart <- (y + 0.5) / 2  # crude initialization
#       })
#     ),
#     class = "family"
#   )
# }

#beta_loglik <- function(y, mu, phi, wt = 1) {
#  sum(wt * (
#    lgamma(phi) - lgamma(mu * phi) - lgamma((1 - mu) * phi) +
#      (mu * phi - 1) * log(y) + ((1 - mu) * phi - 1) * log1p(-y)
#  ))
#}

# fast_beta_regression_mle_r <- function(Xmm, y, start_phi = 10) {
#   # Get starting values for beta from a quick logistic regression
#   start_beta <- fast_logistic_regression(Xmm, y)$b
  
#   # Call the full MLE in C++ without computing standard errors
#   mod_cpp = fast_beta_regression_mle(y, Xmm, start_beta = start_beta, start_phi = start_phi, compute_std_errs = FALSE)
  
#   out <- list(
#     b = mod_cpp$coefficients,
#     phi = mod_cpp$phi
#   )
#   out
# }

# fast_beta_regression_mle_r_with_var <- function(Xmm, y, start_phi = 10) {
#   # Get starting values for beta from a quick logistic regression
#   start_beta <- fast_logistic_regression(Xmm, y)$b
  
#   # Call the full MLE in C++ and compute standard errors
#   mod_cpp = fast_beta_regression_mle(y, Xmm, start_beta = start_beta, start_phi = start_phi, compute_std_errs = TRUE)
  
#   out <- list(
#     b = mod_cpp$coefficients,
#     phi = mod_cpp$phi,
#     # ssq_b_2 is the squared standard error of the treatment effect, which is the second coefficient
#     # The std_errs from C++ includes standard errors for all coefficients and log(phi)
#     # We assume the second element of std_errs corresponds to the treatment effect std error
#     ssq_b_2 = mod_cpp$std_errs[2]^2 
#   )
#   class(out) <- "beta_glm_mle"
#   out
# }


#fast_glm_with_var = function(Xmm, y, glm_function){
#	mod = glm_function(Xmm, y)
#	XtWX = eigen_Xt_times_diag_w_times_X_cpp(Xmm, mod$w)	
#	mod$ssq_b_2 = eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(XtWX, 2)
#	mod	
#}
#
#fast_glm_nb <- function(X, y, maxit = 50, tol = 1e-8, trace = FALSE) {
#  stopifnot(is.matrix(X), is.numeric(y), length(y) == nrow(X))
#
#  n <- length(y)
#  p <- ncol(X)
#
#  beta <- rep(0, p)
#  avg_y = mean(y)
#  theta <- avg_y^2  / (var(y) - avg_y) #method of moments estimate to start
#  for (i in 1 : maxit) {
#    mu <- exp(X %*% beta)
#    W <- mu / (1 + mu / theta)
#    z <- X %*% beta + (y - mu) / mu
#    fit <- lm.wfit(X, z, w = as.vector(W))
#    beta_new <- fit$coefficients
#
#    if (any(is.na(beta_new))) stop("NA in coefficients; possibly singular matrix")
#	
#	opt <- optim(par = theta, fn = neg_loglik_nb_cpp, beta = beta_new, X = X, y = y, method = "L-BFGS-B", lower = 1e-8)
#    theta_new <- opt$par
#
#    if (trace) cat(sprintf("Iter %d: logLik=%.4f  theta=%.4f\n", i, loglik_nb(beta_new, theta_new), theta_new))
#
#    if (max(abs(beta_new - beta)) < tol && abs(theta_new - theta) < tol)
#      break
#
#    beta <- beta_new
#    theta <- theta_new
#  }
#
#  eta <- as.vector(X %*% beta)
#  mu <- exp(eta)
#
#  # Standard errors via observed information
#  W <- mu / (1 + mu / theta)
##  cov_beta <- tryCatch(solve(crossprod(X, X * W)), error = function(e) matrix(NA, p, p))
##  se <- sqrt(diag(cov_beta))
#
#  list(
#    b = beta,
#    ssq_b_2 = eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(crossprod(X, X * W), 2)
#  )
#}

	
#	loglik <- function(n, th, mu, y, w) sum(w * (lgamma(th + 
#        y) - lgamma(th) - lgamma(y + 1) + th * log(th) + y * 
#        log(mu + (y == 0)) - (th + y) * log(th + mu)))
#    link <- log
#    fam0 <- if (missing(init.theta)) 
#        do.call("poisson", list(link = link))
#    else do.call("negative.binomial", list(theta = init.theta, 
#        link = link))
#    mf <- Call <- match.call()
#    m <- match(c("formula", "data", "subset", "weights", "na.action", 
#        "etastart", "mustart", "offset"), names(mf), 0)
#    mf <- mf[c(1, m)]
#    mf$drop.unused.levels <- TRUE
#    mf[[1L]] <- quote(stats::model.frame)
#    mf <- eval.parent(mf)
#    Terms <- attr(mf, "terms")
#    if (method == "model.frame") 
#        return(mf)
#    Y <- model.response(mf, "numeric")
#    X <- if (!is.empty.model(Terms)) 
#        model.matrix(Terms, mf, contrasts)
#    else matrix( NROW(Y), 0)
#    w <- model.weights(mf)
#    if (!length(w)) 
#        w <- rep(1, nrow(mf))
#    else if (any(w < 0)) 
#        stop("negative weights not allowed")
#    offset <- model.offset(mf)
#    mustart <- model.extract(mf, "mustart")
#    etastart <- model.extract(mf, "etastart")
#    n <- length(Y)
#    if (!missing(method)) {
#        if (!exists(method, mode = "function")) 
#            stop(gettextf("unimplemented method: %s", sQuote(method)), 
#                domain = NA)
#        glm.fitter <- get(method)
#    }
#    else {
#        method <- "glm.fit"
#        glm.fitter <- stats::glm.fit
#    }
#    if (control$trace > 1) 
#        message("Initial fit:")
#    fit <- glm.fitter(x = X, y = Y, weights = w, start = start, 
#        etastart = etastart, mustart = mustart, offset = offset, 
#        family = fam0, control = list(maxit = control$maxit, 
#            epsilon = control$epsilon, trace = control$trace > 
#                1), intercept = attr(Terms, "intercept") > 0)
#    class(fit) <- c("glm", "lm")
#    mu <- fit$fitted.values
#    th <- as.vector(theta.ml(Y, mu, sum(w), w, limit = control$maxit, 
#        trace = control$trace > 2))
#    if (control$trace > 1) 
#        message(gettextf("Initial value for 'theta': %f", signif(th)), 
#            domain = NA)
#    fam <- do.call("negative.binomial", list(theta = th, link = link))
#    iter <- 0
#    d1 <- sqrt(2 * max(1, fit$df.residual))
#    d2 <- del <- 1
#    g <- fam$linkfun
#    Lm <- loglik(n, th, mu, Y, w)
#    Lm0 <- Lm + 2 * d1
#    while ((iter <- iter + 1) <= control$maxit && (abs(Lm0 - 
#        Lm)/d1 + abs(del)/d2) > control$epsilon) {
#        eta <- g(mu)
#        fit <- glm.fitter(x = X, y = Y, weights = w, etastart = eta, 
#            offset = offset, family = fam, control = list(maxit = control$maxit, 
#                epsilon = control$epsilon, trace = control$trace > 
#                  1), intercept = attr(Terms, "intercept") > 
#                0)
#        t0 <- th
#        th <- theta.ml(Y, mu, sum(w), w, limit = control$maxit, 
#            trace = control$trace > 2)
#        fam <- do.call("negative.binomial", list(theta = th, 
#            link = link))
#        mu <- fit$fitted.values
#        del <- t0 - th
#        Lm0 <- Lm
#        Lm <- loglik(n, th, mu, Y, w)
#        if (control$trace) {
#            Ls <- loglik(n, th, Y, Y, w)
#            Dev <- 2 * (Ls - Lm)
#            message(sprintf("Theta(%d) = %f, 2(Ls - Lm) = %f", 
#                iter, signif(th), signif(Dev)), domain = NA)
#        }
#    }
#    if (!is.null(attr(th, "warn"))) 
#        fit$th.warn <- attr(th, "warn")
#    if (iter > control$maxit) {
#        warning("alternation limit reached")
#        fit$th.warn <- gettext("alternation limit reached")
#    }
#    if (length(offset) && attr(Terms, "intercept")) {
#        null.deviance <- if (length(Terms)) 
#            glm.fitter(X[, "(Intercept)", drop = FALSE], Y, w, 
#                offset = offset, family = fam, control = list(maxit = control$maxit, 
#                  epsilon = control$epsilon, trace = control$trace > 
#                    1), intercept = TRUE)$deviance
#        else fit$deviance
#        fit$null.deviance <- null.deviance
#    }
#    class(fit) <- c("negbin", "glm", "lm")
#    fit$terms <- Terms
#    fit$formula <- as.vector(attr(Terms, "formula"))
#    Call$init.theta <- signif(as.vector(th), 10)
#    Call$link <- link
#    fit$call <- Call
#    if (model) 
#        fit$model <- mf
#    fit$na.action <- attr(mf, "na.action")
#    if (x) 
#        fit$x <- X
#    if (!y) 
#        fit$y <- NULL
#    fit$theta <- as.vector(th)
#    fit$SE.theta <- attr(th, "SE")
#    fit$twologlik <- as.vector(2 * Lm)
#    fit$aic <- -fit$twologlik + 2 * fit$rank + 2
#    fit$contrasts <- attr(X, "contrasts")
#    fit$xlevels <- .getXlevels(Terms, mf)
#    fit$method <- method
#    fit$control <- control
#    fit$offset <- offset
#    fit



