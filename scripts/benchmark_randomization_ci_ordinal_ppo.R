#!/usr/bin/env Rscript

lib_dir = Sys.getenv("SEQEXP_LIB", unset = "")
if (identical(lib_dir, "")) {
    stop("SEQEXP_LIB must point to the installed package library.")
}
label = Sys.getenv("SEQEXP_LABEL", unset = "portable")
num_cores = as.integer(Sys.getenv("SEQEXP_NUM_CORES", unset = "3"))
r = as.integer(Sys.getenv("SEQEXP_R", unset = "201"))
reps = as.integer(Sys.getenv("SEQEXP_REPS", unset = "3"))
force_mirai = identical(tolower(Sys.getenv("SEQEXP_FORCE_MIRAI", unset = "false")), "true")
inference_class = Sys.getenv("SEQEXP_INFERENCE_CLASS", unset = "InferenceOrdinalMultiPartialProportionalOddsRegr")
nonparallel = strsplit(Sys.getenv("SEQEXP_NONPARALLEL", unset = ""), ",", fixed = TRUE)[[1]]
nonparallel = nonparallel[nonparallel != ""]

.libPaths(c(lib_dir, .libPaths()))
library(EDI)

beta_T = 0.2
SD_NOISE = 0.1
max_n = 100L

X_design = na.omit(as.data.frame(airquality))
X_design = X_design[seq_len(max_n), c("Ozone", "Solar.R", "Wind", "Temp")]
y_base = cut(
    X_design$Temp,
    breaks = quantile(X_design$Temp, probs = seq(0, 1, length.out = 6), na.rm = TRUE),
    include.lowest = TRUE,
    labels = FALSE
)
y_base[is.na(y_base)] = 1L

stable_seed = function(...) {
    as.integer(sum(utf8ToInt(paste(..., collapse = "|"))) %% 2147483647L)
}

make_des_obj = function(seed_key) {
    set.seed(stable_seed("ordinal_ppo_build_benchmark", seed_key))
    des_obj = DesignSeqOneByOneKK14$new(response_type = "ordinal", n = nrow(X_design))
    for (t in seq_len(nrow(X_design))) {
        w_t = des_obj$add_one_subject_to_experiment_and_assign(X_design[t, , drop = FALSE])
        y_t = pmax(1L, as.integer(y_base[t] + (if (w_t == 1) beta_T else 0) + rnorm(1, 0, SD_NOISE)))
        des_obj$add_one_subject_response(t, y_t, dead = 1)
    }
    des_obj
}

set_num_cores(num_cores, force_mirai = force_mirai)
on.exit(unset_num_cores(), add = TRUE)

timings = numeric(reps)
cat(sprintf(
    "Benchmark: %s | class=%s | num_cores=%d | force_mirai=%s | r=%d | reps=%d\n",
    label, inference_class, num_cores, force_mirai, r, reps
))

for (i in seq_len(reps)) {
    des_obj = make_des_obj(paste(label, i, num_cores, r, sep = "|"))
    inf_obj = switch(
        inference_class,
        InferenceOrdinalMultiPartialProportionalOddsRegr = InferenceOrdinalMultiPartialProportionalOddsRegr$new(
            des_obj,
            nonparallel = nonparallel,
            verbose = FALSE
        ),
        InferenceOrdinalMultiOrderedProbitRegr = InferenceOrdinalMultiOrderedProbitRegr$new(
            des_obj,
            verbose = FALSE
        ),
        InferenceOrdinalMultiCumulProbitRegr = InferenceOrdinalMultiCumulProbitRegr$new(
            des_obj,
            verbose = FALSE
        ),
        stop("Unsupported inference class: ", inference_class)
    )
    gc(FALSE)
    t0 = proc.time()[["elapsed"]]
    ci = inf_obj$compute_confidence_interval_rand(
        alpha = 0.05,
        r = r,
        pval_epsilon = 0.05,
        show_progress = FALSE
    )
    timings[i] = round(proc.time()[["elapsed"]] - t0, 3)
    cat(sprintf(
        "rep %d: %.3fs CI=[%.4f, %.4f]\n",
        i, timings[i], ci[1], ci[2]
    ))
}

cat(sprintf("median: %.3fs\n", median(timings)))
