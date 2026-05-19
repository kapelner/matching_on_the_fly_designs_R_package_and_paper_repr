suppressPackageStartupMessages(library(EDI))
suppressPackageStartupMessages(library(data.table))

Nrep = 10000L   # Monte Carlo replications per cell

sim = SimulationFramework$new(
        Nrep                          = Nrep,
        num_cores                     = 12L,
        results_filename              = sprintf("simulations/cmh_exact_sims_plus_greedy_results_Nrep_%d.csv.bz2", Nrep),
        continue_from_last_result_row = TRUE,
        response_type                 = "incidence",
        n                             = c(64L, 128L, 256L),
        p                             = c(1L, 2L, 5L, 10L),
        random_X_draws                = FALSE,
        seed                          = 1984,
        betaT                         = c(1, 0), # 0 → size / type-I error;  1 → power / coverage
        alpha                         = 0.05,
        cond_exp_func_model           = c("linear"), #, "nonlinear" #nonlinear not as interesting for now
        design_classes_and_params     = list(
                                          DesignFixediBCRD,
                                          DesignFixedBinaryMatch,
                                          DesignFixedOptimalBlocks =   list(B = 4),
                                          DesignFixedOptimalBlocks =   list(B = 8),
                                          DesignFixedOptimalBlocks =   list(B = 16),
                                          DesignFixedOptimalBlocks =   list(B = 32),
                                          DesignFixedBlocking =        list(B_target = 4,  exact_num_blocks = TRUE),
                                          DesignFixedBlocking =        list(B_target = 8,  exact_num_blocks = TRUE),
                                          DesignFixedBlocking =        list(B_target = 16, exact_num_blocks = TRUE),
                                          DesignFixedBlocking =        list(B_target = 32, exact_num_blocks = TRUE),
                                          DesignFixedGreedy =          list(objective = "abs_sum_diff"),
                                          DesignFixedRerandomization = list(prop_acceptable = 0.01) 
                                        ),
        inference_classes_and_params  = list(
                                          InferenceIncidWald,
                                          InferenceIncidCMH,
                                          InferenceIncidExtendedRobins 
                                          #InferenceIncidExactBinomial
                                        ),
        inference_types_and_params    = list(
                                          asymp_ci   = list(),
                                          asymp_pval = list()
                                          # exact_ci   = list(),
                                          # exact_pval = list()
                                        ),
        keep_all_intermediate_data    = FALSE
      )  

suppressWarnings(sim$run())
# sm = sim$summarize()

# if (!is.null(sm) && nrow(sm) > 0L) {
#   print(as.data.frame(sm[, .(
#     avg_power = mean(power, na.rm = TRUE),
#     avg_coverage = mean(coverage, na.rm = TRUE)
#   ), by = c("design", "inference", "n", "p", "cond_exp_func_model")]))
# }
