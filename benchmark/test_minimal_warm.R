
library(EDI)
library(microbenchmark)
library(data.table)

set.seed(42)
n = 200; p = 3
X = matrix(rnorm(n * p), n, p); X[, 1] = 1
beta = rnorm(p) * 0.5
eta = X %*% beta
y = rbinom(n, 1, 1/(1+exp(-eta)))

cat("--- Logistic IRLS Weight Start Test ---\n")
p_avg = mean(y)
w_init = rep(p_avg * (1 - p_avg), n)

bm = microbenchmark(
    Cold = fast_logistic_regression_with_var_cpp(X, y, optimization_alg = "irls"),
    Warm_Weights = {
        fast_logistic_regression_with_var_cpp(X, y, optimization_alg = "irls", warm_start_weights = w_init)
    },
    times = 5
)
print(bm)

cat("\n--- Poisson IRLS OLS Start Test ---\n")
y_p = rpois(n, exp(eta))
eta_init = log(y_p + 0.1)
b_init = fast_ols_cpp(X, eta_init)$coefficients
w_p_init = exp(X %*% b_init)

bm_p = microbenchmark(
    Cold = fast_poisson_regression_with_var_cpp(X, y_p, optimization_alg = "irls"),
    Warm_Full = {
        fast_poisson_regression_with_var_cpp(X, y_p, optimization_alg = "irls", warm_start_beta = b_init, warm_start_weights = w_p_init)
    },
    times = 5
)
print(bm_p)
