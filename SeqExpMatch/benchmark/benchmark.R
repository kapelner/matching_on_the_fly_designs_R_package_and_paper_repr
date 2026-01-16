library(SeqExpMatch)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(survival))

# Set seed for reproducibility
set.seed(123)

n = 50
p = 5

cat("Debugging Bootstrap NAs for CRD continuous SimpleMeanDiff\n")

# Define Designs
design_constructors = list(
  "CRD" = function(type) SeqDesignCRD$new(response_type = type)
)

# Inference classes mapped by response type
inference_map = list(
  "continuous" = list(
    "SimpleMeanDiff" = SeqDesignInferenceAllSimpleMeanDiff
  )
)

generate_data = function(n, p, type) {
  X = matrix(rnorm(n * p), nrow = n, ncol = p)
  X_dt = as.data.table(X)
  y = rnorm(n)
  dead = NULL
  list(X = X_dt, y = y, dead = dead)
}

des_name = "CRD"
resp_type = "continuous"
inf_name = "SimpleMeanDiff"

# Generate data
data = generate_data(n, p, resp_type)

# Instantiate design
des_completed = design_constructors[[des_name]](resp_type)
for (i in 1:n) {
  des_completed$add_subject_to_experiment_and_assign(data$X[i])
}
des_completed$add_all_subject_responses(data$y)

# Create an inference object
inf_obj = inference_map[[resp_type]][[inf_name]]$new(des_completed, num_cores = 1, verbose = TRUE)

# Directly get bootstrap samples and print summary
cat("\n--- Bootstrap Samples --- \n")
beta_hat_T_bs_debug = inf_obj$approximate_bootstrap_distribution_beta_hat_T(B = 101)
print(summary(beta_hat_T_bs_debug))


quit("no")

cat("\n\n--- Weibull Regression Benchmark ---\n")

# 1. Add Weibull data generation function
generate_weibull_data <- function(n, p, beta_true, lambda_true, nu_true, censoring_rate = 0.2) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X_dt <- as.data.table(X)

  # Generate Weibull survival times
  eta <- X %*% beta_true
  u <- runif(n)
  time_true <- (-log(u) / (lambda_true * exp(eta)))^(1/nu_true)

  # Generate censoring times
  censor_time <- rexp(n, rate = censoring_rate)

  # Observed time and event indicator
  time <- pmin(time_true, censor_time)
  dead <- as.integer(time_true <= censor_time)

  list(X = X_dt, y = time, dead = dead)
}

# Parameters for Weibull data
p_weibull = 3
beta_true_weibull = rep(0.5, p_weibull)
lambda_true_weibull = 0.1
nu_true_weibull = 1.5 # shape parameter

weibull_data = generate_weibull_data(n, p_weibull, beta_true_weibull, lambda_true_weibull, nu_true_weibull)

# Prepare data for C++ function
y_weibull <- weibull_data$y
dead_weibull <- weibull_data$dead
X_weibull <- as.matrix(weibull_data$X) # RcppEigen expects matrix

# Time the new fast_weibull_regression_with_var_cpp
cat("Timing fast_weibull_regression_with_var_cpp...\n")
time_weibull_cpp <- system.time({
  result_weibull_cpp <- fast_weibull_regression_with_var_cpp(
    y = y_weibull,
    dead = dead_weibull,
    X = X_weibull
  )
})
print(time_weibull_cpp)
print(result_weibull_cpp$coefficients)
print(result_weibull_cpp$log_sigma)
print(result_weibull_cpp$std_errs)
print(result_weibull_cpp$neg_log_lik)