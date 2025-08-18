
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
#    else matrix(, NROW(Y), 0)
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




fast_logistic_regression = function(Xmm, y){
	mod = glmnet::glmnet(
	    x = Xmm,
	    y = y,
	    family = "binomial",
	    alpha = 0,
	    lambda = 0,
	    standardize = FALSE,
	    intercept = FALSE
	)	
	list(b = mod$beta)
}

binomial_link_cache = binomial()

fast_logistic_regression_with_var = function(Xmm, y){
	mod = glm.fit(Xmm, y, family = binomial_link_cache)
	XtWX = eigen_Xt_times_diag_w_times_X_cpp(Xmm, mod$weights)
	list(b = mod$coefficients, ssq_b_2 = eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(XtWX, 2))
}