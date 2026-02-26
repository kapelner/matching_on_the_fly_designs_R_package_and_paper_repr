#' Fast Ordinary Least Squares (OLS) Regression (C++ Backend)
#'
#' This function provides a highly optimized implementation of Ordinary Least Squares (OLS)
#' regression using the Eigen C++ library for numerical computations. It is designed for speed
#' when only the regression coefficients are needed, without variance-covariance matrices
#' or other statistical inference components.
#'
#' @param X A numeric matrix of predictor variables. It is assumed that an intercept column
#'   (e.g., a column of ones) is already included in \code{X} if desired.
#' @param y A numeric vector of the response variable.
#'
#' @return A list containing the following components:
#' \describe{
#' \item{b}{A numeric vector of the estimated regression coefficients.}
#' }
#'
#' @examples
#' \dontrun{
#' # Generate some sample data
#' set.seed(123)
#' n <- 100
#' p <- 3
#' X_mat <- cbind(1, matrix(rnorm(n * (p - 1)), ncol = p - 1)) # X with intercept
#' y_vec <- X_mat %*% c(1, 0.5, -0.2) + rnorm(n, sd = 0.5)
#'
#' # Fit OLS using fast_ols_cpp
#' fit_ols <- fast_ols_cpp(X = X_mat, y = y_vec)
#' print(fit_ols$b)
#'
#' # Compare with standard R lm()
#' lm_fit <- lm(y_vec ~ X_mat - 1) # -1 because X_mat already has intercept
#' print(coef(lm_fit))
#' }
#' @name fast_ols_cpp
#' @rdname fast_ols_cpp
#' @export
NULL

#' Fast Ordinary Least Squares (OLS) Regression with Variance (C++ Backend)
#'
#' This function provides a highly optimized implementation of Ordinary Least Squares (OLS)
#' regression using the Eigen C++ library, including the calculation of the
#' variance-covariance matrix of the coefficients.
#'
#' @param X A numeric matrix of predictor variables. It is assumed that an intercept column
#'   (e.g., a column of ones) is already included in \code{X} if desired.
#' @param y A numeric vector of the response variable.
#' @param j This function will compute the variance of the jth coefficient estimator. Default is 2.
#'
#' @return A list containing the following components:
#' \describe{
#' \item{b}{A numeric vector of the estimated regression coefficients.}
#' \item{ssq_b_j}{The variance of the jth coefficient estimator.}
#' }
#'
#' @examples
#' \dontrun{
#' # Generate some sample data
#' set.seed(123)
#' n <- 100
#' p <- 3
#' X_mat <- cbind(1, matrix(rnorm(n * (p - 1)), ncol = p - 1)) # X with intercept
#' y_vec <- X_mat %*% c(1, 0.5, -0.2) + rnorm(n, sd = 0.5)
#'
#' # Fit OLS with variance using fast_ols_with_var_cpp
#' fit_ols_var <- fast_ols_with_var_cpp(X = X_mat, y = y_vec)
#' print(fit_ols_var$b)
#' print(fit_ols_var$ssq_b_j)
#' }
#'
#' @name fast_ols_with_var_cpp
#' @rdname fast_ols_with_var_cpp
#' @export
NULL

#' Fast Logistic Regression (R Wrapper)
#'
#' This function provides a fast implementation of logistic regression by wrapping
#' a C++ backend function. It is designed to efficiently estimate regression coefficients
#' for binary outcomes without computing variance components.
#'
#' @param Xmm A numeric matrix of predictor variables. It is assumed that an intercept column
#'   (e.g., a column of ones) is already included in \code{Xmm} if desired.
#' @param y A numeric vector of the response variable, expected to be binary (0 or 1).
#'
#' @return A list containing the following component:
#' \describe{
#' \item{b}{A numeric vector of the estimated logistic regression coefficients.}
#' }
#'
#' @examples
#' \dontrun{
#' # Generate some sample data
#' set.seed(123)
#' n <- 100
#' p <- 2
#' X_log <- cbind(1, matrix(rnorm(n * (p - 1)), ncol = p - 1)) # X with intercept
#' beta_log <- c(0.5, 1)
#' lin_pred <- X_log %*% beta_log
#' prob <- exp(lin_pred) / (1 + exp(lin_pred))
#' y_log <- rbinom(n, 1, prob)
#'
#' # Fit logistic regression using fast_logistic_regression
#' fit_logistic <- fast_logistic_regression(Xmm = X_log, y = y_log)
#' print(fit_logistic$b)
#'
#' # Compare with standard R glm()
#' glm_fit <- glm(y_log ~ X_log - 1, family = binomial) # -1 because X_log already has intercept
#' print(coef(glm_fit))
#' }
#'
#' @export
fast_logistic_regression = function(Xmm, y){
  na_b = function() list(b = rep(NA_real_, ncol(Xmm)))
  tryCatch({
    mod = suppressWarnings(fastLogisticRegressionWrap::fast_logistic_regression(
        Xmm = Xmm,
        ybin = as.numeric(y)
    ))
    b = as.vector(mod$coefficients)
    if (!all(is.finite(b))) return(na_b())
    list(b = b)
  }, error = function(e) {
    mod_canonical = suppressWarnings(
      stats::glm.fit(x = Xmm, y = as.numeric(y), family = stats::binomial())
    )
    # Non-convergence (e.g. complete separation) produces unreliable large-magnitude
    # coefficients; signal failure with NA so bootstrap callers can discard the sample.
    if (!isTRUE(mod_canonical$converged)) return(na_b())
    b = as.vector(mod_canonical$coefficients)
    if (!all(is.finite(b))) return(na_b())
    list(b = b)
  })
}

#' Fast Logistic Regression with Variance Calculation (R Wrapper)
#'
#' This function provides a fast implementation of logistic regression, wrapping a C++ backend.
#' It estimates regression coefficients for binary outcomes and specifically computes the
#' variance of the second coefficient (assumed to be the treatment effect).
#'
#' @param Xmm A numeric matrix of predictor variables. It is assumed that an intercept column
#'   (e.g., a column of ones) is already included in \code{Xmm} if desired.
#' @param y A numeric vector of the response variable, expected to be binary (0 or 1).
#' @return A list containing the following components:
#' \describe{
#' \item{b}{A numeric vector of the obtained logistic regression coefficients.}
#' \item{ssq_b_2}{The squared standard error (variance) of the second estimated coefficient,
#'   which typically corresponds to the treatment effect.}
#' }
#'
#' @examples
#' \dontrun{
#' # Generate some sample data
#' set.seed(123)
#' n <- 100
#' p <- 2
#' X_log_var <- cbind(1, rnorm(n)) # X with intercept and one predictor
#' beta_log_var <- c(0.5, 1) # First is intercept, second is treatment effect
#' lin_pred_var <- X_log_var %*% beta_log_var
#' prob_var <- exp(lin_pred_var) / (1 + exp(lin_pred_var))
#' y_log_var <- rbinom(n, 1, prob_var)
#'
#' # Fit logistic regression with variance using fast_logistic_regression_with_var
#' fit_logistic_var <- fast_logistic_regression_with_var(Xmm = X_log_var, y = y_log_var)
#' print(fit_logistic_var$b)
#' print(fit_logistic_var$ssq_b_2)
#'
#' # Compare with standard R glm()
#' glm_fit_var <- glm(y_log_var ~ X_log_var - 1, family = binomial)
#' # Extract squared standard error of the second coefficient (treatment effect)
#' glm_ssq_b_2 <- summary(glm_fit_var)$coefficients[2, 2]^2
#' print(coef(glm_fit_var))
#' print(glm_ssq_b_2)
#' }
#'
#' @export
fast_logistic_regression_with_var = function(Xmm, y){
  # Compute inference on the second coefficient (treatment effect)
  tryCatch({
    mod = suppressWarnings(fastLogisticRegressionWrap::fast_logistic_regression(
        Xmm = Xmm,
        ybin = as.numeric(y),
        do_inference_on_var = 2
    ))
    
    # Check for NaNs in standard errors
    if (any(is.na(mod$se)) || any(!is.finite(mod$se))) {
      stop("NaNs or non-finite values detected in standard errors from fastLogisticRegressionWrap")
    }
    
    list(
      b = as.vector(mod$coefficients),
      ssq_b_2 = mod$se[2]^2  # Square the standard error to get variance
    )
  }, error = function(e) {
    # Fallback to standard R glm if fast version fails
    #warning("fastLogisticRegressionWrap failed, falling back to stats::glm. Error: ", e$message)
    # Using glm.fit for speed
    mod_canonical = stats::glm.fit(x = Xmm, y = as.numeric(y), family = stats::binomial())
    
    # Extract variance of second coefficient (treatment effect)
    # The variance-covariance matrix is (X'WX)^-1
    R <- qr.R(mod_canonical$qr)
    Rinv <- backsolve(R, diag(ncol(R)))
    vcov <- Rinv %*% t(Rinv)
    
    list(
      b = as.vector(mod_canonical$coefficients),
      ssq_b_2 = vcov[2, 2]
    )
  })
}

#' Weibull Regression using survival package internals (fast, coefficients only)
#'
#' This function performs Weibull regression using survival::survreg.fit for speed.
#'
#' @param y Survival times.
#' @param dead Event indicator (1 = event, 0 = censored).
#' @param X A numeric matrix of predictor variables. It is assumed that an intercept column
#'   (e.g., a column of ones) is already included in \code{X} if desired.
#' @return A list containing the following components:
#' \describe{
#' \item{coefficients}{A numeric vector of the estimated Weibull regression coefficients, including the intercept.}
#' \item{log_sigma}{The logarithm of the scale parameter from the Weibull distribution.}
#' \item{vcov}{The variance-covariance matrix of the estimated coefficients.}
#' }
#' @examples
#' \dontrun{
#' # Generate some sample survival data
#' set.seed(123)
#' n <- 100
#' X_weibull <- matrix(rnorm(n), ncol = 1) # One predictor without intercept
#' true_beta <- 0.5
#' true_scale <- 1
#' linear_predictor <- 1 + X_weibull %*% true_beta # Intercept + predictor
#' time_true <- exp(linear_predictor + rweibull(n, shape = 1/true_scale, scale = 1))
#' censor_time <- rexp(n, rate = 0.1)
#' y_weibull <- pmin(time_true, censor_time)
#' dead_weibull <- as.numeric(time_true <= censor_time)
#'
#' # Fit Weibull regression using fast_weibull_regression
#' fit_weibull <- fast_weibull_regression(y = y_weibull, dead = dead_weibull, X = X_weibull)
#' print(fit_weibull$coefficients)
#' print(fit_weibull$log_sigma)
#'
#' # Compare with standard R survreg
#' library(survival)
#' surv_obj <- Surv(y_weibull, dead_weibull)
#' df_weibull <- data.frame(time = y_weibull, status = dead_weibull, X = X_weibull)
#' survreg_fit <- survreg(surv_obj ~ X, data = df_weibull, dist = "weibull")
#' print(coef(survreg_fit))
#' print(log(survreg_fit$scale))
#' }
#' @export
fast_weibull_regression = function(y, dead, X){
  Xmm = as.matrix(X)
  # Check if Xmm has an intercept column (a column of all ones)
  # Assuming intercept is the first column if present
  if (NCOL(Xmm) > 0 && all(Xmm[, 1] == 1)) {
    # If an intercept is present, remove it as survreg adds one automatically
    Xmm_no_intercept = Xmm[, -1, drop = FALSE]
  } else {
    Xmm_no_intercept = Xmm
  }

  # Drop linearly dependent columns before passing to survreg
  Xmm_no_intercept = drop_linearly_dependent_cols(Xmm_no_intercept)$M

  # Use survreg (not survreg.fit) to handle parameter defaults properly
  if (NCOL(Xmm_no_intercept) > 0) {
    # Preserve or create column names
    original_colnames = colnames(Xmm_no_intercept)
    if (is.null(original_colnames)) {
      original_colnames = paste0("X", 1:NCOL(Xmm_no_intercept))
    }

    # Create a data frame for survreg
    df = as.data.frame(Xmm_no_intercept)
    colnames(df) = original_colnames
    df$y = y
    df$dead = dead
    # Wrap column names in backticks to handle special characters
    backticked_colnames = paste0("`", original_colnames, "`")
    formula_str = paste("survival::Surv(y, dead) ~", paste(backticked_colnames, collapse = " + "))
    mod <- survival::survreg(as.formula(formula_str), data = df, dist = "weibull")

    # Extract coefficients and preserve names
    coefficients = as.vector(mod$coefficients)
    names(coefficients) = c("(Intercept)", original_colnames)
  } else {
    # Intercept-only model
    mod <- survival::survreg(survival::Surv(y, dead) ~ 1, dist = "weibull")
    coefficients = as.vector(mod$coefficients)
    names(coefficients) = "(Intercept)"
  }

  vcov = mod$var
  std_errs = if (is.matrix(vcov)) sqrt(diag(vcov)) else rep(NA_real_, length(coefficients) + 1)
  log_sigma = log(mod$scale)
  neg_log_lik = if (!is.null(mod$loglik) && length(mod$loglik) >= 2) -mod$loglik[2] else NA_real_

  # Throw (rather than silently return NaN) so callers like the bootstrap tryCatch can handle failure
  if (any(!is.finite(coefficients))) {
    stop("fast_weibull_regression: survreg returned non-finite coefficients (convergence failure)")
  }
  if (is.matrix(vcov) && any(!is.finite(diag(vcov)))) {
    stop("fast_weibull_regression: survreg returned non-finite variance-covariance (convergence failure)")
  }

  list(
    coefficients = coefficients,
    log_sigma = log_sigma,
    std_errs = std_errs,
    vcov = vcov,
    neg_log_lik = neg_log_lik,
    b = coefficients,
    ssq_b_2 = if (is.matrix(vcov) && nrow(vcov) >= 2) vcov[2, 2] else NA_real_
  )
}

# Internal helper for beta regression safety.
sanitize_beta_response = function(y){
  y = as.numeric(y)
  if (any(!is.finite(y))){
    stop("y must be finite for beta regression")
  }
  eps = .Machine$double.eps
  if (any(y <= 0 | y >= 1)){
    n = length(y)
    if (n > 1){
      y = (y * (n - 1) + 0.5) / n
    }
    y = pmin(pmax(y, eps), 1 - eps)
  }
  y
}

#' Fast Beta Regression (R Wrapper)
#'
#' This function provides a fast implementation of beta regression by wrapping
#' a C++ backend function. It is designed to efficiently estimate regression coefficients
#' for response variables that are continuous and restricted to the (0, 1) interval.
#'
#' @param Xmm A numeric matrix of predictor variables. It is assumed that an intercept column
#'   (e.g., a column of ones) is already included in \code{Xmm} if desired.
#' @param y A numeric vector of the response variable, with values strictly between 0 and 1.
#' @param start_phi A numeric value, the starting value for the precision parameter phi.
#'   Defaults to 10.
#'
#' @return A list containing the following component:
#' \item{b}{A numeric vector of the estimated beta regression coefficients.}
#'
#' @export
fast_beta_regression = function(Xmm, y, start_phi = 10){
  y = sanitize_beta_response(y)
  tryCatch({
    list(b = fast_beta_regression_cpp(Xmm, y, start_phi = start_phi)$coefficients)
  }, error = function(e) {
    warning("fast_beta_regression_cpp failed, falling back to betareg. Error: ", e$message)
    # create a data frame for betareg, removing the intercept from Xmm
    data_df <- as.data.frame(cbind(y, Xmm[, -1, drop = FALSE]))
    # rename columns for formula
    colnames(data_df) <- c("y", paste0("x", 1:(ncol(Xmm)-1)))
    # fit model with control to suppress precision parameter warning
    tryCatch({
      suppressWarnings({
        fit <- betareg::betareg(y ~ ., data = data_df,
                                control = betareg::betareg.control(start = list(phi = start_phi)))
      })
      list(b = coef(fit))
    }, error = function(e2) {
      warning("betareg fallback failed, using OLS on logit(y). Error: ", e2$message)
      list(b = fast_ols_cpp(Xmm, logit(y))$b)
    })
  })
}

#' Fast Beta Regression with Variance Calculation (R Wrapper)
#'
#' This function provides a fast implementation of beta regression, wrapping a C++ backend.
#' It estimates regression coefficients for response variables that are continuous and
#' restricted to the (0, 1) interval, and specifically computes the variance
#' of the second coefficient (assumed to be the treatment effect).
#'
#' @param Xmm A numeric matrix of predictor variables. It is assumed that an intercept column
#'   (e.g., a column of ones) is already included in \code{Xmm} if desired.
#' @param y A numeric vector of the response variable, with values strictly between 0 and 1.
#' @param start_phi A numeric value, the starting value for the precision parameter phi.
#'   Defaults to 10.
#'
#' @return A list containing the following components:
#' \item{b}{A numeric vector of the estimated beta regression coefficients.}
#' \item{ssq_b_2}{The squared standard error (variance) of the second estimated coefficient,
#'   which typically corresponds to the treatment effect.}
#'
#' @importFrom stats vcov
#' @export
fast_beta_regression_with_var = function(Xmm, y, start_phi = 10){
  y = sanitize_beta_response(y)
  tryCatch({
    mod = fast_beta_regression_with_var_cpp(Xmm, y, start_phi = start_phi)
    list(b = mod$coefficients, ssq_b_2 = mod$vcov[2, 2])
  }, error = function(e) {
    warning("fast_beta_regression_with_var_cpp failed, falling back to betareg. Error: ", e$message)
    # create a data frame for betareg, removing the intercept from Xmm
    data_df <- as.data.frame(cbind(y, Xmm[, -1, drop = FALSE]))
    # rename columns for formula
    colnames(data_df) <- c("y", paste0("x", 1:(ncol(Xmm)-1)))
    # fit model with control to suppress precision parameter warning
    tryCatch({
      suppressWarnings({
        fit <- betareg::betareg(y ~ ., data = data_df,
                                control = betareg::betareg.control(start = list(phi = start_phi)))
      })
      # Get the variance of the second coefficient
      vcov_matrix <- vcov(fit)
      list(b = coef(fit), ssq_b_2 = vcov_matrix[2, 2])
    }, error = function(e2) {
      warning("betareg fallback failed, using OLS on logit(y). Error: ", e2$message)
      mod = fast_ols_with_var_cpp(Xmm, logit(y), j = 2L)
      list(b = mod$b, ssq_b_2 = mod$ssq_b_j)
    })
  })
}

#' Fast Cox Proportional Hazards Regression (R Wrapper)
#'
#' This function provides a fast implementation of Cox Proportional Hazards regression
#' by wrapping the \code{glmnet} package's Cox model with \code{lambda = 0}.
#' It is designed to efficiently estimate regression coefficients for time-to-event data
#' without penalization.
#'
#' @param Xmm A numeric matrix of predictor variables. It is assumed that an intercept term
#'   is handled implicitly by the Cox model and should not be included in \code{Xmm}.
#' @param y A numeric vector representing the observed time (event time or censoring time).
#' @param dead A numeric vector (0 or 1) indicating event status (1 for event, 0 for censored).
#'
#' @return A list containing the following component:
#' \item{b}{A numeric vector of the estimated Cox regression coefficients.}
#'
#' @export
fast_coxph_regression = function(Xmm, y, dead){
	mod = glmnet::glmnet(Xmm, survival::Surv(y, dead), family = "cox", lambda = 0)
	list(b = stats::coef(mod))
}

#' Fast Negative Binomial Regression (R Wrapper)
#'
#' This function provides a fast implementation of negative binomial regression by wrapping
#' a C++ backend function. It is designed to efficiently estimate regression coefficients
#' for count data, assuming all observations are "not censored" (i.e., `dead = 1`).
#'
#' @param Xmm A numeric matrix of predictor variables. It is assumed that an intercept column
#'   (e.g., a column of ones) is already included in \code{Xmm} if desired.
#' @param y A numeric vector of the response variable, representing count data.
#'
#' @return A list containing the following component:
#' \item{b}{A numeric vector of the estimated negative binomial regression coefficients.}
#'
#' @importFrom stats glm.fit
#' @importFrom MASS negative.binomial
#' @export
fast_negbin_regression <- function(Xmm, y) {
  # if (nrow(Xmm) > 300 & ncol(Xmm) > 5){ #R's canonical version is faster as IRLS is a different algorithm
  #     
  # } else {
    # list(
    #   b = fast_neg_bin_cpp(
    #     X = Xmm,
    #     y = as.integer(y)
    #   )$b
    # )
  # }
  list(b = coef(stats::glm.fit(
    x = Xmm,
    y = y,
    family = MASS::negative.binomial(theta = 1)  # initial theta
  )))
}

#' Fast Negative Binomial Regression with Variance Calculation (R Wrapper)
#'
#' This function provides a fast implementation of negative binomial regression, wrapping
#' a C++ backend. It is designed to efficiently estimate regression coefficients for count
#' data and computes the variance of the second coefficient (assumed to be the treatment effect),
#' assuming all observations are "not censored" (i.e., `dead = 1`).
#'
#' @param Xmm A numeric matrix of predictor variables. It is assumed that an intercept column
#'   (e.g., a column of ones) is already included in \code{Xmm} if desired.
#' @param y A numeric vector of the response variable, representing count data.
#'
#' @return A list containing the following components:
#' \item{b}{A numeric vector of the estimated negative binomial regression coefficients.}
#' \item{ssq_b_2}{The squared standard error (variance) of the second estimated coefficient,
#'   which typically corresponds to the treatment effect.}
#'
#' @importFrom stats glm.fit
#' @importFrom MASS negative.binomial
#' @export
fast_negbin_regression_with_var <- function(Xmm, y) {
  # Use glm.fit for efficiency
  mod <- stats::glm.fit(
    x = Xmm,
    y = y,
    family = MASS::negative.binomial(theta = 1)  # initial theta
  )

  # Compute variance-covariance matrix from QR decomposition
  # The variance-covariance matrix is (X'WX)^-1 * dispersion
  # where W is the weight matrix from IRLS

  # Get the R matrix from QR decomposition
  R <- qr.R(mod$qr)

  # Compute (X'WX)^-1 = R^-1 %*% t(R^-1)
  Rinv <- backsolve(R, diag(ncol(R)))
  vcov_unscaled <- Rinv %*% t(Rinv)

  # Compute dispersion parameter
  # For negative binomial, dispersion = sum(weights * residuals^2) / df.residual
  dispersion <- sum(mod$weights * mod$residuals^2) / mod$df.residual

  # Scale the variance-covariance matrix by dispersion
  vcov <- vcov_unscaled * dispersion

  # Extract variance of second coefficient (treatment effect)
  ssq_b_2 <- vcov[2, 2]

  list(
    b = as.vector(mod$coefficients),
    ssq_b_2 = ssq_b_2
  )
}


































# fast_beta_regression <- function(Xmm, y,
#                              start_phi = 10,
#                              bounds_logphi = c(log(1e-3), log(1e4)),
#                              control = stats::glm.control(epsilon=1e-8, maxit=100)) {
   
#   weights <- rep(1, nrow(Xmm))

#   # Use C++ logistic regression to find beta (quasi-likelihood, independent of phi)
#   # This replaces the repetitive glm.fit calls inside the optimization loop
#   mod_log = fast_logistic_regression_cpp(as.matrix(Xmm), as.numeric(y), maxit = control$maxit, tol = control$epsilon)
#   b = mod_log$b
#   mu = 1 / (1 + exp(-drop(Xmm %*% b)))
#   # guard against numerical drift to 0/1
#   mu <- pmin(pmax(mu, 1e-12), 1 - 1e-12)

#   # objective: negative log-likelihood profiled over beta (beta is fixed now)
#   obj <- function(logphi) {
#     phi = exp(logphi)
#     # return NEGATIVE log-likelihood
#     -beta_loglik_cpp(y, mu, phi, wt = weights)
#   }

#   opt <- stats::nlminb(start = log(start_phi), objective = obj,
#                 lower = bounds_logphi[1], upper = bounds_logphi[2])

#   phi_hat <- as.numeric(exp(opt$par))
  
#   # Construct weights for the object matching beta_family/glm.fit behavior
#   # weights = (1+phi) * mu * (1-mu)
#   w_final = (1 + phi_hat) * mod_log$w

#   out <- list(
#     b = b,
#     phi = phi_hat,
#     converged_inner = TRUE,
#     w = w_final
#   )
#   class(out) <- "beta_glm_profile"
#   out
# }

# fast_beta_regression_with_var <- function(Xmm, y,
#                              start_phi = 10,
#                              bounds_logphi = c(log(1e-3), log(1e4)),
#                              control = stats::glm.control(epsilon=1e-8, maxit=100)) {
#     fit_hat = fast_beta_regression(Xmm, y, start_phi = start_phi, bounds_logphi = bounds_logphi, control = control) 
#     XtWX <- crossprod(sqrt(fit_hat$w) * Xmm)
#     fit_hat$ssq_b_2 = eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(XtWX, 2)
#     fit_hat
# }

# binomial_link_cache = binomial()

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
