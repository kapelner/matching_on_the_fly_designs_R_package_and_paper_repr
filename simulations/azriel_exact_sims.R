library(EDI)
suppressPackageStartupMessages(library(data.table))

# ── Tunable parameters ────────────────────────────────────────────────────────
Nrep       = 1000L   # Monte Carlo replications per cell
betaTs     = c(1, 0) # 0 → size / type-I error;  1 → power / coverage
alpha      = 0.05
ns         = c(64L, 128L, 256L)
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
  InferenceIncidenceWald,
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

          # ── Collect and save ──────────────────────────────────────────────────────────
          results_dt = rbindlist(all_results)

          # Re-order columns for readability
          setcolorder(results_dt, c("betaT", "n", "p", "data_type", "design", "inference", "method"))
          setorder(results_dt, betaT, n, p, data_type, design, inference, method)   
          
          fwrite(results_dt, out_file)
          message(sprintf("\nSaved to %s", out_file))

          print(as.data.frame(results_dt[, .(avg_power = mean(power, na.rm = TRUE), avg_coverage = mean(coverage, na.rm = TRUE)) , by = c("design", "inference", "n", "p")]))
        }
      }
    }
  }
}



