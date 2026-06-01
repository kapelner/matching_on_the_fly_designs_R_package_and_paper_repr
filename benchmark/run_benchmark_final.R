
library(EDI)
library(data.table)

N_VAL = 300
P_VAL = 5
B_VAL = 100 
R_VAL = 100 
J_VAL = 30  
PB_VAL = 20

RESULTS_CSV = "warm_starts_final_results.csv"

all_inf = grep("^Inference[A-Z]", ls("package:EDI"), value = TRUE)
inf_names = all_inf[!grepl("Abstract|Suite|Mixin|Compound|Likelihood$|Asymp$|Jackknife$|Bootstrap$|MLEorKMSummaryTable", all_inf)]

# Worker script content
worker_script = '
library(EDI)
library(data.table)
args = commandArgs(trailingOnly = TRUE)
cls_name = args[1]
N_VAL = as.numeric(args[2])
P_VAL = as.numeric(args[3])
B_VAL = as.numeric(args[4])
R_VAL = as.numeric(args[5])
J_VAL = as.numeric(args[6])
PB_VAL = as.numeric(args[7])

get_rt = function(cn) {
    if (grepl("Ordinal|Ridit|Jonck", cn)) return("ordinal")
    if (grepl("Incid", cn)) return("incidence")
    if (grepl("Count", cn)) return("count")
    if (grepl("Prop|Beta", cn)) return("proportion")
    if (grepl("Survival|KM|Cox|Weibull|LogRank", cn)) return("survival")
    return("continuous")
}

get_design_object = function(cn, rt, n) {
    if (grepl("KK", cn) || grepl("Exact|PairedSignTest", cn)) {
        return(DesignFixedBinaryMatch$new(response_type = rt, n = n))
    }
    if (grepl("ExtendedRobins", cn)) {
        return(DesignFixedBlocking$new(response_type = rt, n = n, m = rep(1:(n/2), each=2), equal_block_sizes = FALSE))
    }
    return(DesignFixedBernoulli$new(response_type = rt, n = n))
}

calc_s = function(tc, tw) {
    if (is.na(tc) || is.na(tw) || !is.finite(tc) || !is.finite(tw)) return("N/S")
    if (tc < 0.005) return("< 2%") 
    s = (tc - tw) / tc * 100
    if (s < 2 && s > -2) return("< 2%")
    if (s < 0) return(sprintf("%.1f%%", s))
    sprintf("+%.1f%%", s)
}

rt = get_rt(cls_name)

# Quick trial for dynamic N/P scaling if baseline cold start is too fast (sub-millisecond)
d_temp = get_design_object(cls_name, rt, N_VAL)
if (rt == "ordinal") d_temp$.__enclos_env__$private$ordinal_levels = as.character(1:5)
set.seed(42)
X_temp = as.data.frame(matrix(rnorm(N_VAL * P_VAL), N_VAL, P_VAL))
colnames(X_temp) = paste0("V", 2:(P_VAL+1))
d_temp$add_all_subjects_to_experiment(X_temp)
d_temp$assign_w_to_all_subjects()
y_temp = if(rt == "incidence") { rbinom(N_VAL, 1, 0.5) } else if(rt == "count") { rpois(N_VAL, 1) } else if(rt == "proportion") { runif(N_VAL) } else if(rt == "survival") { rexp(N_VAL) } else { sample(1:5, N_VAL, replace=TRUE) }
d_temp$add_all_subject_responses(y_temp)
inf_temp = get(cls_name)$new(d_temp)
inf_temp$.__enclos_env__$private$fit_warm_start_enabled = FALSE
t_single = tryCatch(system.time(inf_temp$compute_estimate())["elapsed"], error = function(e) NA)

if (!is.na(t_single)) {
    if (t_single < 0.00005) { # extremely fast closed-form models (Wilcoxon, mean difference, etc.)
        N_VAL = 20000
        P_VAL = 5
        B_VAL = 1000
        R_VAL = 1000
        J_VAL = 300
        PB_VAL = 100
    } else if (t_single < 0.0003) { # very fast models (OLS, simple GLM)
        N_VAL = 8000
        P_VAL = 35
        B_VAL = 500
        R_VAL = 500
        J_VAL = 150
        PB_VAL = 60
    } else if (t_single < 0.001) {
        N_VAL = 3500
        P_VAL = 20
        B_VAL = 250
        R_VAL = 250
        J_VAL = 80
        PB_VAL = 40
    } else if (t_single < 0.003) {
        N_VAL = 1500
        P_VAL = 12
        B_VAL = 150
        R_VAL = 150
        J_VAL = 50
        PB_VAL = 25
    }
}

if (grepl("KK", cls_name) || grepl("Exact|PairedSignTest", cls_name)) {
    N_VAL = min(N_VAL, 300)
}

if (grepl("Wilcox|Jonckheere", cls_name)) {
    N_VAL = min(N_VAL, 1000)
}

d = get_design_object(cls_name, rt, N_VAL)
if (rt == "ordinal") d$.__enclos_env__$private$ordinal_levels = as.character(1:5)

set.seed(42)
X = as.data.frame(matrix(rnorm(N_VAL * P_VAL), N_VAL, P_VAL))
colnames(X) = paste0("V", 2:(P_VAL+1))
d$add_all_subjects_to_experiment(X)
d$assign_w_to_all_subjects()

y = if(rt == "incidence") { rbinom(N_VAL, 1, 0.5) } else if(rt == "count") { rpois(N_VAL, 1) } else if(rt == "proportion") { runif(N_VAL) } else if(rt == "survival") { rexp(N_VAL) } else { sample(1:5, N_VAL, replace=TRUE) }
d$add_all_subject_responses(y)

# Warm Object Setup
inf_w = get(cls_name)$new(d)
inf_w$.__enclos_env__$private$fit_warm_start_enabled = TRUE
mle = tryCatch(inf_w$compute_estimate(), error = function(e) NULL)
if (!is.null(mle)) {
    pw = inf_w$.__enclos_env__$private
    if (is.list(mle)) {
        pw$set_fit_warm_start(
            mle$params %||% mle$b %||% mle$coefficients, 
            type = if (!is.null(mle$params)) "params" else "beta",
            fisher = mle$fisher_information %||% mle$observed_information %||% mle$XtWX,
            weights = mle$w %||% mle$mu,
            force_pd = TRUE
        )
    } else {
        pw$set_fit_warm_start(as.numeric(mle), "beta")
    }
    # Intentionally do NOT clear structural caches (cached_hardened_X_cov,
    # cached_design_matrix, etc.) populated by compute_estimate() above.
    # In production, compute_estimate() is always called before any resampling
    # method, so these caches are legitimately available and are part of the
    # measured advantage of operating in the warm state.
}

# Cold Object Setup
inf_c = get(cls_name)$new(d)
inf_c$.__enclos_env__$private$fit_warm_start_enabled = FALSE

# Mock context for Bootstrap/JK
mock_ctx = list(n_units = N_VAL, row_to_unit = 1:N_VAL)
inf_c$.__enclos_env__$private$current_bayesian_bootstrap_context = mock_ctx
inf_w$.__enclos_env__$private$current_bayesian_bootstrap_context = mock_ctx

# RAND
res_r = "N/S"
if (!(rt == "incidence" && is.null(inf_w$.__enclos_env__$private$custom_randomization_statistic_function))) {
    # Warm-up (r=1): equalize CPU instruction/data cache and R allocator state for both objects.
    # The r=1 cache key differs from r=R_VAL, so no result is reused in the timed call.
    tryCatch(inf_c$compute_rand_two_sided_pval(r=1L, show_progress=FALSE), error=function(e)NULL)
    tryCatch(inf_w$compute_rand_two_sided_pval(r=1L, show_progress=FALSE), error=function(e)NULL)
    t_rc = tryCatch(system.time(inf_c$compute_rand_two_sided_pval(r=R_VAL, show_progress=FALSE))["elapsed"], error=function(e) NA)
    t_rw = tryCatch(system.time(inf_w$compute_rand_two_sided_pval(r=R_VAL, show_progress=FALSE))["elapsed"], error=function(e) NA)
    res_r = calc_s(t_rc, t_rw)
}

# BOOT
# Warm-up (B=1): equalize CPU and allocator state.  B=1 cache key != B=B_VAL.
tryCatch(inf_c$compute_bootstrap_confidence_interval(B=1L, show_progress=FALSE), error=function(e)NULL)
tryCatch(inf_w$compute_bootstrap_confidence_interval(B=1L, show_progress=FALSE), error=function(e)NULL)
t_bc = tryCatch(system.time(inf_c$compute_bootstrap_confidence_interval(B=B_VAL, show_progress=FALSE))["elapsed"], error=function(e) NA)
t_bw = tryCatch(system.time(inf_w$compute_bootstrap_confidence_interval(B=B_VAL, show_progress=FALSE))["elapsed"], error=function(e) NA)
res_b = calc_s(t_bc, t_bw)

# JK
# Warm-up (1 iteration): populates cached_design_matrix and cached_hardened_X_cov on both
# objects equally, and seeds the R allocator free list with same-sized objects for both.
tryCatch({ww=rep(1,N_VAL);ww[1]=0;inf_c$compute_estimate_with_bootstrap_weights(ww,TRUE)}, error=function(e)NULL)
tryCatch({ww=rep(1,N_VAL);ww[1]=0;inf_w$compute_estimate_with_bootstrap_weights(ww,TRUE)}, error=function(e)NULL)
t_jc = tryCatch(system.time(for(k in 1:J_VAL){w=rep(1,N_VAL);w[k]=0;inf_c$compute_estimate_with_bootstrap_weights(w,TRUE)})["elapsed"], error=function(e) NA)
t_jw = tryCatch(system.time(for(k in 1:J_VAL){w=rep(1,N_VAL);w[k]=0;inf_w$compute_estimate_with_bootstrap_weights(w,TRUE)})["elapsed"], error=function(e) NA)
res_j = calc_s(t_jc, t_jw)

# PB
is_pb_sup = tryCatch(inf_w$.__enclos_env__$private$supports_lik_ratio_param_bootstrap(), error = function(e) FALSE)
res_p = "N/S"
if (isTRUE(is_pb_sup)) {
    # Warm-up (B=1): B=1 cache key != B=PB_VAL.
    tryCatch(inf_c$compute_lik_ratio_bootstrap_two_sided_pval(B=1L, show_progress=FALSE), error=function(e)NULL)
    tryCatch(inf_w$compute_lik_ratio_bootstrap_two_sided_pval(B=1L, show_progress=FALSE), error=function(e)NULL)
    t_pc = tryCatch(system.time(inf_c$compute_lik_ratio_bootstrap_two_sided_pval(B=PB_VAL, show_progress=FALSE))["elapsed"], error=function(e) NA)
    t_pw = tryCatch(system.time(inf_w$compute_lik_ratio_bootstrap_two_sided_pval(B=PB_VAL, show_progress=FALSE))["elapsed"], error=function(e) NA)
    res_p = calc_s(t_pc, t_pw)
}

cat(sprintf("%s,%s,%s,%s,%s\n", cls_name, res_r, res_b, res_j, res_p))
'

writeLines(worker_script, "benchmark/worker_v2.R")

cat("Starting definitive robust benchmark...\n")
if (!file.exists(RESULTS_CSV)) {
    fwrite(data.table(Path=character(), Rand=character(), Boot=character(), JK=character(), PB=character()), RESULTS_CSV)
}

done_paths = fread(RESULTS_CSV)$Path

for (cn in inf_names) {
    if (cn %in% done_paths) next
    cat(sprintf("Processing %s... ", cn))
    
    res = system2("Rscript", c("benchmark/worker_v2.R", cn, N_VAL, P_VAL, B_VAL, R_VAL, J_VAL, PB_VAL), stdout = TRUE, stderr = TRUE)
    
    # The last line should be the CSV row
    row_str = res[length(res)]
    if (grepl(",", row_str)) {
        write(row_str, RESULTS_CSV, append = TRUE)
        cat("Done.\n")
    } else {
        cat("Failed or Crashed.\n")
        # Log failure reason to stderr
        cat(paste(res, collapse="\n"), file=stderr())
        write(sprintf("%s,N/S,N/S,N/S,N/S", cn), RESULTS_CSV, append = TRUE)
    }
}
cat("Done.\n")
unlink("benchmark/worker_v2.R")
