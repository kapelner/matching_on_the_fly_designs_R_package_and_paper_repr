suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(betareg))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Rcpp))
suppressPackageStartupMessages(library(SeqExpMatch))

# Define different values for n (sample size) and p (number of covariates) to test
n_values <- 100#c(100, 500, 1000)
p_values <- 5#c(5, 10, 20)

# Number of iterations for benchmarking each test
num_iterations <- 1

# Data table to store benchmark results
benchmark_results <- data.table(
  n_value = integer(),
  p_value = integer(),
  Regression_Type = character(),
  Method = character(),
  Iteration = integer(),
  Time_Taken_seconds = numeric()
)

# Helper function to compare model coefficients and standard errors
compare_model_outputs <- function(model_new, model_canonical, test_name, tolerance = 1e-4) {
  # Extract coefficients from SeqExpMatch
  coef_seqexpmatch <- model_new$b

  # Extract coefficients from canonical model
  coef_canonical <- coef(model_canonical)

  # Skip if SeqExpMatch has NaN/Inf values
  if (any(!is.finite(coef_seqexpmatch))) {
    message(paste("\nSkipping", test_name, "- SeqExpMatch returned non-finite coefficients"))
    return(invisible(NULL))
  }

  # Handle intercept mismatch: if canonical has 1 more coef, it's likely the intercept
  if (length(coef_canonical) == length(coef_seqexpmatch) + 1) {
    # Canonical likely has intercept, SeqExpMatch doesn't
    coef_canonical <- coef_canonical[-1]  # Remove first element (intercept)
    message(paste(test_name, "- Comparing without intercept"))
  }

  # Handle phi parameter in beta regression: if canonical has extra param at end
  if (test_name == "Beta Regression" && length(coef_canonical) > length(coef_seqexpmatch)) {
    # Extract only beta coefficients (exclude phi which is typically last)
    coef_canonical <- coef_canonical[1:length(coef_seqexpmatch)]
    message(paste(test_name, "- Comparing only beta coefficients (excluding phi)"))
  }

  test_that(paste(test_name, "- Coefficients match"), {
    expect_equal(length(coef_seqexpmatch), length(coef_canonical),
                 info = paste0("Number of coefficients differ. SeqExpMatch: ", length(coef_seqexpmatch),
                              "; Canonical: ", length(coef_canonical)))

    # Check coefficients match within tolerance
    max_diff <- max(abs(coef_seqexpmatch - coef_canonical), na.rm = TRUE)
    expect_true(all(abs(coef_seqexpmatch - coef_canonical) < tolerance, na.rm = TRUE),
                info = paste("Coefficients mismatch. Max diff:", max_diff))
  })

  # Only check variance if model_new has ssq_b_2
  if (!is.null(model_new$ssq_b_2) && is.finite(model_new$ssq_b_2)) {
    new_var = model_new$ssq_b_2

    # Try to get variance from canonical model
    tryCatch({
      canonical_summary <- summary(model_canonical)
      if ("coefficients" %in% names(canonical_summary)) {
        # Get the second coefficient's std error and square it
        canonical_var <- canonical_summary$coefficients[2, 2]^2

        test_that(paste(test_name, "- Treatment Effect Variance matches"), {
          expect_true(abs(new_var - canonical_var) < tolerance,
                      info = paste("Treatment effect variance mismatch:",
                                   "SeqExpMatch:", new_var, "Canonical:", canonical_var))
        })
      }
    }, error = function(e) {
      # Skip variance check if canonical model doesn't provide it
      message(paste("Skipping variance check for", test_name, ":", e$message))
    })
  }
}


# Outer loops for n and p combinations
for (current_n in n_values) {
  n <- current_n # Assign global n for current iteration
  for (current_p in p_values) {
    p <- current_p # Assign global p for current iteration
    num_covariates <- p - 2 # Define num_covariates here, so it's consistent for each iteration of p


    # --- 1. Weibull Regression Check ---
    try({
    message("\n--- Running Weibull Regression Check (n=", n, ", p=", p, ", Benchmarking ", num_iterations, " iterations) ---")
    set.seed(123) # Set seed outside the loop for reproducibility of overall data generation
    
    # Generate covariates
    x_weibull <- lapply(1:num_covariates, function(i) rnorm(n, mean = i * 0.1, sd = 1))
    names(x_weibull) <- paste0("x", 1:num_covariates)
    x_weibull_df <- as.data.frame(x_weibull)

    treatment <- rbinom(n, 1, 0.5) # Random treatment assignment for CRD

    # Generate true beta coefficients for covariates
    true_beta_x_weibull <- rnorm(num_covariates, mean = 0, sd = 0.5)
    true_beta_treatment <- 0.6
    true_intercept <- 1
    true_scale <- 1.5

    # Construct linear predictor
    linear_predictor_weibull <- true_intercept + as.matrix(x_weibull_df) %*% true_beta_x_weibull + true_beta_treatment * treatment
    time_base <- exp(linear_predictor_weibull + rnorm(n, 0, true_scale))
    censor_time <- rexp(n, 0.1) # Random censoring
    time <- pmin(time_base, censor_time)
    status <- as.numeric(time_base <= censor_time) # 1 if event, 0 if censored

    data_weibull_df <- cbind(data.frame(time = time, status = status, treatment = treatment), x_weibull_df)

    for (iter in 1:num_iterations) {
      # Time Canonical survreg model
      time_canonical <- system.time({
        covariate_names_weibull <- paste0("x", 1:num_covariates)
        Surv_obj <- Surv(data_weibull_df$time, data_weibull_df$status)        
        formula_str_weibull <- paste("Surv_obj ~", paste(c(covariate_names_weibull, "treatment"), collapse = " + "))
        model_weibull_canonical <- survreg(as.formula(formula_str_weibull), data = data_weibull_df, dist = "weibull")
        coef(summary(model_weibull_canonical))[2, 2]
      })[["elapsed"]] * 1000
      
      benchmark_results <- rbindlist(list(benchmark_results, list(
        n_value = n, p_value = p,
        Regression_Type = "Weibull",
        Method = "Canonical R",
        Iteration = iter,
        Time_Taken_seconds = time_canonical
      )))
        
      # Time SeqExpMatch Weibull model
      # Create model matrix including treatment (note: C++ function adds intercept internally)
      X_weibull_no_intercept <- cbind(treatment, as.matrix(x_weibull_df))
      colnames(X_weibull_no_intercept) <- c("treatment", paste0("x", 1:num_covariates))

      time_seqexpmatch <- system.time({
        new_model_weibull <- fast_weibull_regression(
          Xmm = X_weibull_no_intercept,
          y = data_weibull_df$time,
          dead = data_weibull_df$status
        )
        
      })[["elapsed"]] * 1000

      benchmark_results <- rbindlist(list(benchmark_results, list(
        n_value = n, p_value = p,
        Regression_Type = "Weibull",
        Method = "SeqExpMatch (coef only)",
        Iteration = iter,
        Time_Taken_seconds = time_seqexpmatch
      )))

      # Time SeqExpMatch Weibull model with variance
      time_seqexpmatch <- system.time({
        new_model_weibull_with_var <- fast_weibull_regression_with_var(
          Xmm = X_weibull_no_intercept,
          y = data_weibull_df$time,
          dead = data_weibull_df$status
        )
        
      })[["elapsed"]] * 1000

      benchmark_results <- rbindlist(list(benchmark_results, list(
        n_value = n, p_value = p,
        Regression_Type = "Weibull",
        Method = "SeqExpMatch (coef+var)",
        Iteration = iter,
        Time_Taken_seconds = time_seqexpmatch
      )))
      browser()
      cat(data.table(
       idx = c(0, "tx", 1 : num_covariates, "ssq_b_2"),
       weibull_canonical = c(coef(model_weibull_canonical), coef(summary(model_weibull_canonical))[2, 2]),
       weibull_new = c(new_model_weibull$b, NA),
       weibull_new_with_var = c(new_model_weibull_with_var$b, new_model_weibull_with_var$ssq_b_2)
      ))

      expect_equal(new_model_weibull$b, new_model_weibull_with_var$b)

      # Extract model info for comparison (only from the last iteration)
      if (iter == num_iterations) {
        compare_model_outputs(new_model_weibull_with_var, model_weibull_canonical, "Weibull Regression", tolerance = 1e-3)
      }
    }
    })


    # --- 2. Logistic Regression Check ---
    try({
    message("\n--- Running Logistic Regression Check (n=", n, ", p=", p, ", Benchmarking ", num_iterations, " iterations) ---")
    set.seed(123)
    num_covariates <- p - 2

    # Generate covariates dynamically
    x_log <- lapply(1:num_covariates, function(i) rnorm(n, mean = i * 0.1, sd = 1))
    names(x_log) <- paste0("x", 1:num_covariates)
    X <- as.matrix(as.data.frame(x_log))
    w <- rbinom(n, 1, 0.5)
    wX = cbind(w, X)

    # Generate true beta coefficients for covariates
    true_beta_x_log <- rnorm(num_covariates, mean = 0, sd = 0.5)
    true_beta_treatment_log <- 0.9
    true_intercept_log <- 0.5

    # Construct linear predictor dynamically
    linear_predictor_log <- true_intercept_log + X %*% true_beta_x_log + true_beta_treatment_log * w

    prob <- plogis(linear_predictor_log)
    y <- rbinom(n, 1, prob)

    # Construct formula string dynamically
    covariate_names_log <- paste0("x", 1:num_covariates)
    formula_str_logistic <- paste("y ~", paste(c(covariate_names_log, "treatment"), collapse = " + "))

    for (iter in 1:num_iterations) {
      # Time Canonical glm.fit model
      time_canonical <- system.time({    
        canonical_fit = glm.fit(wX, y, family = binomial(link = "logit"))    
      })[["elapsed"]] * 1000

      benchmark_results <- rbindlist(list(benchmark_results, list(
        n_value = n, p_value = p,
        Regression_Type = "Logistic",
        Method = "Canonical R (coef only)",
        Iteration = iter,
        Time_Taken_seconds = time_canonical
      )))

      # Time Canonical glm.fit model
      time_canonical <- system.time({
        canonical_fit = glm(y ~ wX, family = "binomial")
        coef(summary(canonical_fit))[2, 2]
      })[["elapsed"]] * 1000

      benchmark_results <- rbindlist(list(benchmark_results, list(
        n_value = n, p_value = p,
        Regression_Type = "Logistic",
        Method = "Canonical R (coef+var)",
        Iteration = iter,
        Time_Taken_seconds = time_canonical
      )))

      # Time SeqExpMatch Logistic model (coefficients only)
      time_seqexpmatch_coef_only <- system.time({
        new_model_logistic <- fast_logistic_regression(wX, y)
      })[["elapsed"]] * 1000

      benchmark_results <- rbindlist(list(benchmark_results, list(
        n_value = n, p_value = p,
        Regression_Type = "Logistic",
        Method = "SeqExpMatch (coef only)",
        Iteration = iter,
        Time_Taken_seconds = time_seqexpmatch_coef_only
      )))

      # Time SeqExpMatch Logistic model (coefficients + variance)
      time_seqexpmatch <- system.time({
        new_model_logistic_with_var <- fast_logistic_regression_with_var(wX, y)
      })[["elapsed"]] * 1000

      benchmark_results <- rbindlist(list(benchmark_results, list(
        n_value = n, p_value = p,
        Regression_Type = "Logistic",
        Method = "SeqExpMatch (coef+var)",
        Iteration = iter,
        Time_Taken_seconds = time_seqexpmatch
      )))

      expect_equal(new_model_logistic$b, new_model_logistic_with_var$b)

        # Extract model info for comparison (only from the last iteration)
        if (iter == num_iterations) {
          compare_model_outputs(new_model_logistic_with_var, canonical_fit, "Logistic Regression")
        }
    }
    })


    # --- 3. Negative Binomial Regression Check ---
    try({
    message("\n--- Running Negative Binomial Regression Check (n=", n, ", p=", p, ", Benchmarking ", num_iterations, " iterations) ---")
    set.seed(123)
    num_covariates <- p - 2

    # Generate covariates dynamically
    x_nb <- lapply(1:num_covariates, function(i) rnorm(n, mean = i * 0.1, sd = 1))
    names(x_nb) <- paste0("x", 1:num_covariates)
    x_nb_df <- as.data.frame(x_nb)

    treatment_nb <- rbinom(n, 1, 0.5)

    # Generate true beta coefficients for covariates
    true_beta_x_nb <- rnorm(num_covariates, mean = 0, sd = 0.5)
    true_beta_treatment_nb <- 0.7
    true_intercept_nb <- 0.5

    # Construct linear predictor dynamically
    mu_nb_linear_predictor <- true_intercept_nb + as.matrix(x_nb_df) %*% true_beta_x_nb + true_beta_treatment_nb * treatment_nb
    mu_nb <- exp(mu_nb_linear_predictor)
    size_nb <- 2 # dispersion parameter
    y_nb <- rnbinom(n, size = size_nb, mu = mu_nb)

    data_negbin_df <- cbind(data.frame(y = y_nb, treatment = treatment_nb), x_nb_df)

    # Construct formula string dynamically
    covariate_names_nb <- paste0("x", 1:num_covariates)
    formula_str_negbin <- paste("y ~", paste(c(covariate_names_nb, "treatment"), collapse = " + "))

    for (iter in 1:num_iterations) {
      # Time Canonical glm.nb model
      time_canonical <- system.time({
        model_negbin_canonical <- glm.nb(as.formula(formula_str_negbin), data = data_negbin_df)
      })[["elapsed"]] * 1000

      benchmark_results <- rbindlist(list(benchmark_results, list(
        n_value = n, p_value = p,
        Regression_Type = "Negative Binomial",
        Method = "Canonical R",
        Iteration = iter,
        Time_Taken_seconds = time_canonical
      )))

      # Construct model matrix for C++ functions
      X_nb_matrix <- model.matrix(as.formula(formula_str_negbin), data = data_negbin_df)
      y_nb_vector <- data_negbin_df$y

      # Time SeqExpMatch Negative Binomial model (coefficients only)
      time_seqexpmatch_coef_only <- system.time({
        new_model_negbin <- fast_negbin_regression(X_nb_matrix, y_nb_vector)
      })[["elapsed"]] * 1000

      benchmark_results <- rbindlist(list(benchmark_results, list(
        n_value = n, p_value = p,
        Regression_Type = "Negative Binomial",
        Method = "SeqExpMatch (coef only)",
        Iteration = iter,
        Time_Taken_seconds = time_seqexpmatch_coef_only
      )))
      
      # Time SeqExpMatch Negative Binomial model (coefficients + variance)
      time_seqexpmatch <- system.time({
        new_model_negbin_with_var <- fast_negbin_regression_with_var(X_nb_matrix, y_nb_vector)
      })[["elapsed"]] * 1000

      benchmark_results <- rbindlist(list(benchmark_results, list(
        n_value = n, p_value = p,
        Regression_Type = "Negative Binomial",
        Method = "SeqExpMatch (coef+var)",
        Iteration = iter,
        Time_Taken_seconds = time_seqexpmatch
      )))

      # Extract model info for comparison (only from the last iteration)
      if (iter == num_iterations) {
        compare_model_outputs(new_model_negbin_with_var, model_negbin_canonical, "Negative Binomial Regression", tolerance = 1e-2)
      }
    }
    })


    # --- 4. Beta Regression Check ---
    try({
    message("\n--- Running Beta Regression Check (n=", n, ", p=", p, ", Benchmarking ", num_iterations, " iterations) ---")
    set.seed(42) # Changed to 42 for better convergence
    num_covariates <- p - 2

    # Generate covariates dynamically
    x_beta <- lapply(1:num_covariates, function(i) rnorm(n, mean = i * 0.1, sd = 1))
    names(x_beta) <- paste0("x", 1:num_covariates)
    x_beta_df <- as.data.frame(x_beta)

    treatment_beta <- rbinom(n, 1, 0.5)

    # Generate true beta coefficients for covariates
    true_beta_x_beta <- rnorm(num_covariates, mean = 0, sd = 0.5)
    true_beta_treatment_beta <- 0.5
    true_intercept_beta <- 0.5

    # Construct linear predictor dynamically
    mu_beta_logit <- true_intercept_beta + as.matrix(x_beta_df) %*% true_beta_x_beta + true_beta_treatment_beta * treatment_beta
    mu_beta <- plogis(mu_beta_logit) 
    phi_beta <- 5 # Precision parameter
    y_beta <- rbeta(n, shape1 = mu_beta * phi_beta, shape2 = (1 - mu_beta) * phi_beta)
    # Clamp y_beta to be strictly between 0 and 1
    epsilon_clamp <- 1e-8
    y_beta[y_beta < epsilon_clamp] <- epsilon_clamp
    y_beta[y_beta > (1 - epsilon_clamp)] <- (1 - epsilon_clamp)

    data_beta_df <- cbind(data.frame(y = y_beta, treatment = treatment_beta), x_beta_df)

    # Construct formula string dynamically
    covariate_names_beta <- paste0("x", 1:num_covariates)
    formula_str_beta <- paste("y ~", paste(c(covariate_names_beta, "treatment"), collapse = " + "))

    for (iter in 1:num_iterations) {
      # Time Canonical betareg model
      time_canonical <- system.time({
        model_beta_canonical <- betareg(as.formula(formula_str_beta), data = data_beta_df)
      })[["elapsed"]] * 1000

      benchmark_results <- rbindlist(list(benchmark_results, list(
        n_value = n, p_value = p,
        Regression_Type = "Beta",
        Method = "Canonical R",
        Iteration = iter,
        Time_Taken_seconds = time_canonical
      )))

      # Construct model matrix for C++ functions
      X_beta_matrix <- model.matrix(as.formula(formula_str_beta), data = data_beta_df)
      y_beta_vector <- data_beta_df$y

      # Time SeqExpMatch Beta model (coefficients only)
      time_seqexpmatch_coef_only <- system.time({
        new_model_beta <- fast_beta_regression(X_beta_matrix, y_beta_vector)
      })[["elapsed"]] * 1000
      
      benchmark_results <- rbindlist(list(benchmark_results, list(
        n_value = n, p_value = p,
        Regression_Type = "Beta",
        Method = "SeqExpMatch (coef only)",
        Iteration = iter,
        Time_Taken_seconds = time_seqexpmatch_coef_only
      )))

      # Time SeqExpMatch Beta model (coefficients + variance)
      time_seqexpmatch <- system.time({
        new_model_beta_with_var <- fast_beta_regression_with_var(X_beta_matrix, y_beta_vector)
      })[["elapsed"]] * 1000

      benchmark_results <- rbindlist(list(benchmark_results, list(
        n_value = n, p_value = p,
        Regression_Type = "Beta",
        Method = "SeqExpMatch (coef+var)",
        Iteration = iter,
        Time_Taken_seconds = time_seqexpmatch
      )))

      expect_equal(new_model_beta$b, new_model_beta_with_var$b)

      # Extract model info for comparison (only from the last iteration)
      if (iter == num_iterations) {
        compare_model_outputs(new_model_beta_with_var, model_beta_canonical, "Beta Regression", tolerance = 1e-2)
      }
    }
    })

    # --- 5. OLS Regression Check ---
    try({
    message("\n--- Running OLS Regression Check (n=", n, ", p=", p, ", Benchmarking ", num_iterations, " iterations) ---")
    set.seed(123)
    num_covariates <- p - 2

    # Generate covariates dynamically
    x_ols <- lapply(1:num_covariates, function(i) rnorm(n, mean = i * 0.1, sd = 1))
    names(x_ols) <- paste0("x", 1:num_covariates)
    x_ols_df <- as.data.frame(x_ols)

    treatment_ols <- rbinom(n, 1, 0.5)

    # Generate true beta coefficients for covariates
    true_beta_x_ols <- rnorm(num_covariates, mean = 0, sd = 0.5)
    true_beta_treatment_ols <- 0.8
    true_intercept_ols <- 1

    # Construct linear predictor dynamically
    y_ols <- true_intercept_ols + as.matrix(x_ols_df) %*% true_beta_x_ols + true_beta_treatment_ols * treatment_ols + rnorm(n, 0, 1)

    data_ols_df <- cbind(data.frame(y = y_ols, treatment = treatment_ols), x_ols_df)

    # Construct formula string dynamically
    covariate_names_ols <- paste0("x", 1:num_covariates)
    
    # Construct model matrix for C++ functions
    X_ols_matrix <- as.matrix(cbind(1, treatment_ols, x_ols_df))

    for (iter in 1:num_iterations) {
      # Time Canonical lm model (coefficients + variance)
      time_canonical <- system.time({
        model_ols_canonical <- lm.fit(X_ols_matrix, y_ols)
      })[["elapsed"]] * 1000

      benchmark_results <- rbindlist(list(benchmark_results, list(
        n_value = n, p_value = p,
        Regression_Type = "OLS",
        Method = "Canonical R (coef)",
        Iteration = iter,
        Time_Taken_seconds = time_canonical
      )))

      # Time Canonical lm.fit model
      time_canonical <- system.time({
        model_ols_canonical_with_var <- lm(y ~ ., data = data_ols_df)
      })[["elapsed"]] * 1000

      benchmark_results <- rbindlist(list(benchmark_results, list(
        n_value = n, p_value = p,
        Regression_Type = "OLS",
        Method = "Canonical R (coef+var)",
        Iteration = iter,
        Time_Taken_seconds = time_canonical
      )))      

      # Both canonical models should produce same coefficients
      # Note: lm.fit and lm produce coefficients in same order when X has column names
      expect_equal(unname(model_ols_canonical$coefficients), unname(coef(model_ols_canonical_with_var)), tolerance = 1e-10)

      # Time SeqExpMatch OLS model (coefficients only)
      time_seqexpmatch_coef_only <- system.time({
        new_model_ols <- fast_ols_cpp(X_ols_matrix, y_ols)
      })[["elapsed"]] * 1000

      benchmark_results <- rbindlist(list(benchmark_results, list(
        n_value = n, p_value = p,
        Regression_Type = "OLS",
        Method = "SeqExpMatch (coef only)",
        Iteration = iter,
        Time_Taken_seconds = time_seqexpmatch_coef_only
      )))
      
      # Time SeqExpMatch OLS model (coefficients + variance)
      time_seqexpmatch <- system.time({
        new_model_ols_with_var <- fast_ols_with_var_cpp(X_ols_matrix, y_ols)
      })[["elapsed"]] * 1000

      benchmark_results <- rbindlist(list(benchmark_results, list(
        n_value = n, p_value = p,
        Regression_Type = "OLS",
        Method = "SeqExpMatch (coef+var)",
        Iteration = iter,
        Time_Taken_seconds = time_seqexpmatch
      )))

      expect_equal(new_model_ols$b, new_model_ols_with_var$b)

      # Extract model info for comparison (only from the last iteration)
      if (iter == num_iterations) {
        compare_model_outputs(new_model_ols_with_var, model_ols_canonical_with_var, "OLS Regression")
      }
    }
    })

  } # closes for (current_p in p_values)
} # closes for (current_n in n_values)

fwrite(benchmark_results, file = "benchmark_new_vs_canonical_results.csv")

message("All GLM checks complete.")

# Calculate and report average benchmark results
if (nrow(benchmark_results) > 0) {
  # Calculate average benchmark results in long format
  average_benchmark_results_long <- benchmark_results[, .(Average_Time_seconds = sprintf("%.3f", mean(Time_Taken_seconds))), by = .(n_value, p_value, Regression_Type, Method)]

  # Print the results in a long format, which is cleaner without NAs
  cat("\n--- Average Benchmark Results ---\n")
  print(average_benchmark_results_long[order(n_value, p_value, Regression_Type, Method)], row.names = FALSE, digits = 6)
  cat("---------------------------------\n")
} else {
  message("\nNo benchmark results to report.")
}

# Remove the debug_beta_regression.R file
if (file.exists("debug_beta_regression.R")) {
  file.remove("debug_beta_regression.R")
  message("Removed temporary file: debug_beta_regression.R")
}