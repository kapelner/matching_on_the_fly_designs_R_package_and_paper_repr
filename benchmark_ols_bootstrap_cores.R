rm(list = ls())
set.seed(1)
library(SeqExpMatch)
library(data.table)
library(dplyr)
library(mlbench)
library(AppliedPredictiveModeling)
library(PTE)

max_n_dataset = 150
source("package_tests/_dataset_load.R")

dataset_name = "diamonds"
response_type = "continuous"
design_type = "CRD"
B_samples = 10000
beta_T = 1
SD_NOISE = 0.1

apply_treatment_effect_and_noise = function(y_t, w_t, response_type){
    eps = rnorm(1, 0, SD_NOISE)
    bt = ifelse(w_t == 1, beta_T, 0)
    return(y_t + bt + eps)
}

D = datasets_and_response_models[[dataset_name]]
n = nrow(D$X)
dead = rep(1, n)
y = D$y_original[[response_type]]

cat("Initializing design ", design_type, " for dataset ", dataset_name, " (n = ", n, ")...
", sep="")

seq_des_obj = SeqDesignCRD$new(response_type = response_type, n = n)

for (t in 1 : n){
    w_t = seq_des_obj$add_subject_to_experiment_and_assign(D$X[t, ])
    y_t = apply_treatment_effect_and_noise(y[t], w_t, response_type)
    seq_des_obj$add_subject_response(t, y_t, dead[t])
}

timings = data.table(num_cores = integer(), duration_sec = numeric(), pval = numeric())

for (cores in 1:6) {
    cat(sprintf("
Benchmarking C++ OpenMP OLS bootstrap with %d core(s) and B = %d...
", cores, B_samples))
    
    seq_des_inf = SeqDesignInferenceContinMultOLS$new(seq_des_obj, num_cores = cores)
    
    start_time = proc.time()[["elapsed"]]
    
    pval = seq_des_inf$compute_bootstrap_two_sided_pval(B = B_samples, na.rm = TRUE)
    
    end_time = proc.time()[["elapsed"]]
    duration = end_time - start_time
    
    cat(sprintf("  Completed in %.3f seconds.
", duration))
    cat(sprintf("  Bootstrap P-value: %.5f
", pval))
    
    timings = rbind(timings, data.table(
        num_cores = cores, 
        duration_sec = duration, 
        pval = pval
    ))
}

cat("

--- C++ OLS Bootstrap Benchmark Results (ContinMultOLS) ---
")
print(timings)
