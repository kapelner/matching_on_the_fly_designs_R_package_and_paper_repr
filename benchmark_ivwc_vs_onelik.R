library(EDI)
library(dplyr)
library(tidyr)

# Parse command line arguments
args = commandArgs(trailingOnly = TRUE)
N_REP = 20
if (length(args) > 0) {
  N_REP = as.integer(args[1])
  if (is.na(N_REP)) {
    stop("Nrep must be an integer")
  }
}
cat("Using N_REP =", N_REP, "\n")

# Helper function to run a single comparison
run_comparison = function(response_type, ivwc_class, onelik_class, label, n = 100, Nrep = N_REP, betaT = 0.5) {
  cat("\n========================================================================\n")
  cat("Benchmarking:", label, "(", response_type, ")\n")
  cat("========================================================================\n")
  
  sim = SimulationFramework$new(
    response_type = response_type,
    design_classes = list(DesignSeqOneByOneKK21),
    inference_classes = list(ivwc_class, onelik_class),
    n = n,
    p = 5,
    data_type = "linear",
    Nrep = Nrep,
    betaT = betaT,
    verbose = FALSE
  )
  
  sim$run()
  summary_res = sim$summarize()
  
  # Add label and response_type for identification
  summary_res$path_label = label
  summary_res$response_type = response_type
  return(summary_res)
}

# Define the paths to benchmark
# Note: Some classes have different names for univariate vs multivariate or continuous vs proportion.
# We use representative ones here.

all_results = list()

# 1. Quantile Regression (Continuous)
all_results[[1]] = run_comparison("continuous", 
                                  InferenceContinKKQuantileRegrIVWC, 
                                  InferenceContinKKQuantileRegrOneLik, 
                                  "Quantile Regression")

# 2. Robust Regression (Continuous)
all_results[[2]] = run_comparison("continuous", 
                                  InferenceContinKKRobustRegrIVWC, 
                                  InferenceContinKKRobustRegrOneLik, 
                                  "Robust Regression")

# 3. Conditional Poisson (Count)
all_results[[3]] = run_comparison("count", 
                                  InferenceCountKKCPoissonIVWC, 
                                  InferenceCountKKCPoissonOneLik, 
                                  "Conditional Poisson")

# 4. Hurdle Poisson (Count)
all_results[[4]] = run_comparison("count", 
                                  InferenceCountKKHurdlePoissonIVWC, 
                                  InferenceCountKKHurdlePoissonOneLik, 
                                  "Hurdle Poisson")

# 5. Conditional Logistic (Incidence)
all_results[[5]] = run_comparison("incidence", 
                                  InferenceIncidKKClogitIVWC, 
                                  InferenceIncidKKClogitOneLik, 
                                  "Conditional Logistic")

# 6. Conditional Logistic Plus GLMM (Incidence)
all_results[[6]] = run_comparison("incidence", 
                                  InferenceIncidKKClogitPlusGLMMIVWC, 
                                  InferenceIncidKKClogitPlusGLMMOneLik, 
                                  "Clogit + GLMM")

# 7. Clayton Copula Weibull AFT (Survival)
all_results[[7]] = run_comparison("survival", 
                                  InferenceSurvivalKKClaytonCopulaIVWC, 
                                  InferenceSurvivalKKClaytonCopulaOneLik, 
                                  "Clayton Copula")

# 8. LWA Cox Partial Likelihood (Survival)
all_results[[8]] = run_comparison("survival", 
                                  InferenceSurvivalKKLWACoxIVWC, 
                                  InferenceSurvivalKKLWACoxOneLik, 
                                  "LWA Cox")

# 9. Stratified Cox Partial Likelihood (Survival)
all_results[[9]] = run_comparison("survival", 
                                  InferenceSurvivalKKStratCoxIVWC, 
                                  InferenceSurvivalKKStratCoxOneLik, 
                                  "Stratified Cox")

# Combine all results
final_res = do.call(rbind, all_results)

# Format and display the results focused on Power
# Note: simulate results uses 'inference' column name for the class name
benchmark_table = final_res %>%
  mutate(type = ifelse(grepl("IVWC", inference), "IVWC", "OneLik")) %>%
  select(path_label, response_type, type, power, MSE, n_est, n_pow) %>%
  arrange(path_label, type)

cat("\n\nFINAL BENCHMARK RESULTS (Power and MSE):\n")
print(as.data.frame(benchmark_table))

# Save to file
write.csv(benchmark_table, "benchmark_ivwc_vs_onelik_results.csv", row.names = FALSE)
cat("\nResults saved to benchmark_ivwc_vs_onelik_results.csv\n")
