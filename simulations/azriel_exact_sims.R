# simulations/azriel_exact_sims.R
#
# Compares exact / block-based incidence inference methods:
#   InferenceIncidAzriel           (asymp: block-stratified SE)
#   InferenceIncidExtendedRobins   (asymp: extended Robins block SE)
#   InferenceIncidenceExactBinomial (exact: matched-pair binomial)
#
# across four fixed designs (prob_T = 0.5):
#   FixedDesigniBCRD, FixedDesignBinaryMatch,
#   FixedDesignOptimalBlocks, FixedDesignBlocking
#
# for response_type = "incidence", data_type in {linear, nonlinear},
# n in {100, 200}, p in {1, 2, 5, 10}.
#
# Note: nonlinear (Friedman) requires p >= 5; lower-p runs are linear only.
# Note: Azriel and Robins require equal block sizes; they are silently skipped
#       for designs that do not satisfy that constraint.
# Note: ExactBinomial requires a matched-pair design; it is silently skipped
#       for non-matched designs.
#
# Usage:
#   devtools::load_all("EDI")
#   source("simulations/azriel_exact_sims.R")

library(EDI)

# ── Tunable parameters ────────────────────────────────────────────────────────
Nrep       = 500L   # Monte Carlo replications per cell
betaTs     = c(1, 0) # 0 → size / type-I error;  1 → power / coverage
alpha      = 0.05
ns         = c(128L, 256L)
ps         = c(1L, 2L, 5L, 10L)
data_types = c("linear", "nonlinear")
out_file   = "simulations/azriel_exact_sims_results.csv"

# ── Fixed: design classes ───────────────────────────────────────
design_classes = list(
  FixedDesigniBCRD,
  FixedDesignBinaryMatch,
  FixedDesignOptimalBlocks,
  FixedDesignOptimalBlocks,
  FixedDesignOptimalBlocks,
  FixedDesignOptimalBlocks,
  FixedDesignBlocking,
  FixedDesignBlocking,
  FixedDesignBlocking,
  FixedDesignBlocking
)
design_params = list(
  list(),
  list(),
  #FixedDesignOptimalBlocks
  list(B = 4),
  list(B = 8),
  list(B = 16),
  list(B = 32),
  #FixedDesignBlocking
  list(B_preferred = 4),
  list(B_preferred = 8),
  list(B_preferred = 16),
  list(B_preferred = 32)
)

# ── Fixed: inference classes ───────────────────────────────────────
inference_classes = list(
  InferenceAllSimpleMeanDiff,
  InferenceIncidAzriel,
  InferenceIncidExtendedRobins,
  InferenceIncidenceExactBinomial
)

# Only asymptotic (Azriel / Robins) and exact (ExactBinomial) types needed.
# Warnings about InferenceIncidAzriel / InferenceIncidExtendedRobins not
# inheriting InferenceExact are expected and suppressed below.
inf_types = c("asymp_ci", "asymp_pval", "exact_ci", "exact_pval")

# ── Run ───────────────────────────────────────────────────────────────────────
all_results = list()

for (betaT in betaTs) {
  for (n in ns) {
    for (p in ps) {
      for (data_type in data_types) {

        # Friedman nonlinear function needs at least 5 covariates
        if (data_type == "nonlinear" && p < 5L) next

        message(sprintf(
          "\n── betaT=%g  n=%d  p=%d  data_type=%s ───────────────────",
          betaT, n, p, data_type))

     
        sim = suppressWarnings(
          SimulationFramework$new(
            response_type     = "incidence",
            design_classes    = design_classes,
            inference_classes = inference_classes,
            n                 = n,
            p                 = p,
            data_type         = data_type,
            Nrep              = Nrep,
            betaT             = betaT,
            alpha             = alpha,
            inf_types         = inf_types,
            design_params     = design_params
          )
        )

        # Suppress the expected per-class exact-inheritance warnings during run
        suppressWarnings(sim$run())


        # print out the design
        ds = sim$get_designs()
        for (i in seq_along(ds)){
          message(sprintf(
            "\n── design #%d: %s, # blocks: %d \n",
            i, class(ds[[i]])[1], data.table::uniqueN(ds[[i]]$get_block_ids())))
          print(ds[[i]]$summarize_blocks())
        }        

        sm = sim$summarize()
        if (!is.null(sm) && nrow(sm) > 0L) {
          sm[, `:=`(betaT = betaT, n = n, p = p, data_type = data_type)]
          all_results[[length(all_results) + 1L]] = sm
        }
      }
    }
  }
}

# ── Collect and save ──────────────────────────────────────────────────────────
results_dt = rbindlist(all_results)

# Re-order columns for readability
setcolorder(results_dt, c("betaT", "n", "p", "data_type", "design", "inference", "method"))
setorder(results_dt, betaT, n, p, data_type, design, inference, method)

message("\n── Results ──────────────────────────────────────────────────────────")
print(results_dt, nrows = 200)

fwrite(results_dt, out_file)
message(sprintf("\nSaved to %s", out_file))
