suppressPackageStartupMessages(library(EDI))
suppressPackageStartupMessages(library(data.table))

args = commandArgs(trailingOnly = TRUE)
N_REP = 100
if (length(args) > 0) {
  N_REP = as.integer(args[1])
  if (is.na(N_REP)) stop("Nrep must be an integer")
}
cat("Using N_REP =", N_REP, "\n")

set_num_cores(1L)
on.exit(unset_num_cores(), add = TRUE)

inference_path_labels = c(  
  InferenceContinKKOLSIVWC = "OLS",
  InferenceContinKKOLSOneLik = "OLS",
  InferenceContinKKQuantileRegrIVWC = "Quantile (Contin)",
  InferenceContinKKQuantileRegrOneLik = "Quantile (Contin)",
  InferenceContinKKRobustRegrIVWC = "Robust Regr",
  InferenceContinKKRobustRegrOneLik = "Robust Regr",
  InferenceCountKKCPoissonIVWC = "Cond Poisson",
  InferenceCountKKCPoissonOneLik = "Cond Poisson",
  InferenceCountKKHurdlePoissonIVWC = "Hurdle Poisson",
  InferenceCountKKHurdlePoissonOneLik = "Hurdle Poisson",
  InferenceIncidKKClogitIVWC = "Clogit",
  InferenceIncidKKClogitOneLik = "Clogit",
  InferenceIncidKKClogitPlusGLMMIVWC = "Clogit + GLMM",
  InferenceIncidKKClogitPlusGLMMOneLik = "Clogit + GLMM",
  InferencePropKKQuantileRegrIVWC = "Quantile (Prop)",
  InferencePropKKQuantileRegrOneLik = "Quantile (Prop)",
  InferenceSurvivalKKClaytonCopulaIVWC = "Clayton Copula",
  InferenceSurvivalKKClaytonCopulaOneLik = "Clayton Copula",
  InferenceSurvivalKKLWACoxIVWC = "LWA Cox",
  InferenceSurvivalKKLWACoxOneLik = "LWA Cox",
  InferenceSurvivalKKStratCoxIVWC = "Stratified Cox",
  InferenceSurvivalKKStratCoxOneLik = "Stratified Cox"
)

sim = SimulationFramework$new(
  response_type = c("continuous", "count", "incidence", "proportion", "survival"),
  design_classes_and_params = list(
    DesignSeqOneByOneKK21
  ),
  inference_classes_and_params = list(
    InferenceContinKKOLSIVWC,
    InferenceContinKKOLSOneLik,
    InferenceContinKKQuantileRegrIVWC,
    InferenceContinKKQuantileRegrOneLik,
    InferenceContinKKRobustRegrIVWC,
    InferenceContinKKRobustRegrOneLik,
    InferenceCountKKCPoissonIVWC,
    InferenceCountKKCPoissonOneLik,
    InferenceCountKKHurdlePoissonIVWC,
    InferenceCountKKHurdlePoissonOneLik,
    InferenceIncidKKClogitIVWC,
    InferenceIncidKKClogitOneLik,
    InferenceIncidKKClogitPlusGLMMIVWC,
    InferenceIncidKKClogitPlusGLMMOneLik,
    InferencePropKKQuantileRegrIVWC,
    InferencePropKKQuantileRegrOneLik,
    # InferenceSurvivalKKClaytonCopulaIVWC,
    # InferenceSurvivalKKClaytonCopulaOneLik,
    InferenceSurvivalKKLWACoxIVWC,
    InferenceSurvivalKKLWACoxOneLik,
    InferenceSurvivalKKStratCoxIVWC,
    InferenceSurvivalKKStratCoxOneLik
  ),
  n = 100,
  p = 5,
  cond_exp_func_model = "linear",
  Nrep = N_REP,
  betaT = c(0, 0.25),
  num_cores = 5L,
  inference_types_and_params = list(
    asymp_pval = list()
  ),
  verbose = TRUE,
  turn_off_asserts_for_speed = FALSE,
  results_filename = "benchmark_ivwc_vs_onelik.csv.bz2",
  continue_from_last_result_row = TRUE
)

sim$run()
summary_res = sim$summarize()
benchmark_table = as.data.table(summary_res)[
  inference_type == "asymp_pval" & (n_pow > 0L | n_est > 0L)
][
  , `:=`(
    Path = unname(inference_path_labels[inference]),
    type = fifelse(grepl("IVWC", inference), "IVWC", "OneLik")
  )
][
  order(response_type, Path, type),
  .(Path, response_type, inference, type, power, size, MSE, n_est, n_pow, n_size)
]

cat("\n\nFINAL BENCHMARK RESULTS (Power and MSE):\n")
print(as.data.frame(benchmark_table), row.names = FALSE)

write.csv(benchmark_table, "benchmark_ivwc_vs_onelik_results.csv", row.names = FALSE)
cat("\nResults saved to benchmark_ivwc_vs_onelik_results.csv\n")
