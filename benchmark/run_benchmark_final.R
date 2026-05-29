
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
pkgload::load_all("EDI")
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
    if (grepl("Incid", cn)) return("incidence")
    if (grepl("Count", cn)) return("count")
    if (grepl("Prop|Beta", cn)) return("proportion")
    if (grepl("Survival|KM|Cox|Weibull|LogRank", cn)) return("survival")
    if (grepl("Ordinal|Ridit|Jonck", cn)) return("ordinal")
    return("continuous")
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
d = if(grepl("KK", cls_name)) DesignFixedBinaryMatch$new(rt, n=N_VAL) else DesignFixedBernoulli$new(rt, n=N_VAL)
if (rt == "ordinal") d$.__enclos_env__$private$ordinal_levels = as.character(1:5)

set.seed(42)
X = as.data.frame(matrix(rnorm(N_VAL * P_VAL), N_VAL, P_VAL))
colnames(X) = paste0("V", 2:(P_VAL+1))
d$add_all_subjects_to_experiment(X)
d$assign_w_to_all_subjects()

y = if(rt == "incidence") rbinom(N_VAL, 1, 0.5) 
    else if(rt == "count") rpois(N_VAL, 1) 
    else if(rt == "proportion") runif(N_VAL) 
    else if(rt == "survival") rexp(N_VAL) 
    else sample(1:5, N_VAL, replace=TRUE)
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
    t_rc = tryCatch(system.time(inf_c$compute_rand_two_sided_pval(r=R_VAL, show_progress=FALSE))["elapsed"], error=function(e) NA)
    t_rw = tryCatch(system.time(inf_w$compute_rand_two_sided_pval(r=R_VAL, show_progress=FALSE))["elapsed"], error=function(e) NA)
    res_r = calc_s(t_rc, t_rw)
}

# BOOT
t_bc = tryCatch(system.time(inf_c$compute_bootstrap_confidence_interval(B=B_VAL, show_progress=FALSE))["elapsed"], error=function(e) NA)
t_bw = tryCatch(system.time(inf_w$compute_bootstrap_confidence_interval(B=B_VAL, show_progress=FALSE))["elapsed"], error=function(e) NA)
res_b = calc_s(t_bc, t_bw)

# JK
t_jc = tryCatch(system.time(for(k in 1:J_VAL){w=rep(1,N_VAL);w[k]=0;inf_c$compute_estimate_with_bootstrap_weights(w,TRUE)})["elapsed"], error=function(e) NA)
t_jw = tryCatch(system.time(for(k in 1:J_VAL){w=rep(1,N_VAL);w[k]=0;inf_w$compute_estimate_with_bootstrap_weights(w,TRUE)})["elapsed"], error=function(e) NA)
res_j = calc_s(t_jc, t_jw)

# PB
is_pb_sup = tryCatch(inf_w$.__enclos_env__$private$supports_lik_ratio_param_bootstrap(), error = function(e) FALSE)
res_p = "N/S"
if (isTRUE(is_pb_sup)) {
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
