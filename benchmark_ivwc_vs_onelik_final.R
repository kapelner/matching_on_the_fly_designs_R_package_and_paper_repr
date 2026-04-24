# Ensure local EDI is loaded
devtools::load_all("EDI")
library(dplyr)

# Parse command line arguments
args = commandArgs(trailingOnly = TRUE)
N_REP = 10
if (length(args) > 0) {
  N_REP = as.integer(args[1])
  if (is.na(N_REP)) {
    stop("Nrep must be an integer")
  }
}
cat("Using N_REP =", N_REP, "\n")

SAMPLE_SIZE = 100
EFFECT_SIZE = 0.5 

run_path_benchmark = function(label, resp_type, ivwc_cls, onelik_cls) {
  cat("\n", rep("-", 80), sep = "")
  cat("\nPATH: ", label, " [", resp_type, "]\n", sep = "")
  cat(rep("-", 80), "\n", sep = "")
  
  sim = SimulationFramework$new(
    response_type     = resp_type,
    design_classes    = list(DesignSeqOneByOneKK21),
    inference_classes = list(ivwc_cls, onelik_cls),
    n = SAMPLE_SIZE,
    p = 5,
    Nrep = N_REP,
    betaT = EFFECT_SIZE,
    verbose = TRUE,
    inf_types = c("asymp_pval")
  )
  
  # Run and capture results
  res = tryCatch({
    sim$run()
    # SimulationFramework$summarize() returns a data.table with:
    # design, inference, method, MSE, n_est, coverage, n_cov, ci_length, power, n_pow, design_params, inference_params
    sim$summarize()
  }, error = function(e) {
    message("Error in simulation for ", label, ": ", e$message)
    return(NULL)
  })
  
  if (!is.null(res)) {
    res$Path = label
    res$response_type = resp_type # Manually add since it's missing from summarize()
    res$Variant = ifelse(grepl("IVWC", res$inference), "IVWC", "OneLik")
    return(res %>% select(Path, response_type, Variant, power, MSE, n_est, n_pow))
    }
    return(NULL)
    }

# Define all comparison pairs
paths = list(
  list("Quantile (Contin)", "continuous", InferenceContinKKQuantileRegrIVWC,  InferenceContinKKQuantileRegrOneLik),
  list("Quantile (Prop)",   "proportion", InferencePropKKQuantileRegrIVWC,    InferencePropKKQuantileRegrOneLik),
  list("Robust Regr",       "continuous", InferenceContinKKRobustRegrIVWC,    InferenceContinKKRobustRegrOneLik),
  list("Cond Poisson",      "count",      InferenceCountKKCPoissonIVWC,       InferenceCountKKCPoissonOneLik),
  list("Hurdle Poisson",    "count",      InferenceCountKKHurdlePoissonIVWC,  InferenceCountKKHurdlePoissonOneLik),
  list("Clogit",            "incidence",  InferenceIncidKKClogitIVWC,         InferenceIncidKKClogitOneLik),
  list("Clogit + GLMM",     "incidence",  InferenceIncidKKClogitPlusGLMMIVWC, InferenceIncidKKClogitPlusGLMMOneLik),
  list("Clayton Copula",    "survival",   InferenceSurvivalKKClaytonCopulaIVWC, InferenceSurvivalKKClaytonCopulaOneLik),
  list("LWA Cox",           "survival",   InferenceSurvivalKKLWACoxIVWC,      InferenceSurvivalKKLWACoxOneLik),
  list("Stratified Cox",    "survival",   InferenceSurvivalKKStratCoxIVWC,    InferenceSurvivalKKStratCoxOneLik)
)

results_list = list()
for (p in paths) {
  results_list[[length(results_list) + 1]] = run_path_benchmark(p[[1]], p[[2]], p[[3]], p[[4]])
}

# Combine and Report
final_results = do.call(rbind, results_list)

cat("\n\n" , rep("=", 40), sep = "")
cat("\nSUMMARY OF POWER COMPARISON (IVWC vs OneLik)\n")
cat(rep("=", 40), "\n", sep = "")
print(as.data.frame(final_results), row.names = FALSE)

# Save to CSV
write.csv(final_results, "ivwc_vs_onelik_power_benchmark.csv", row.names = FALSE)
cat("\nResults saved to 'ivwc_vs_onelik_power_benchmark.csv'\n")
