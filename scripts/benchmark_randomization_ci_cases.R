#!/usr/bin/env Rscript

lib_dir = Sys.getenv("EDI_LIB", unset = "")
if (identical(lib_dir, "")) {
    stop("EDI_LIB must point to the installed package library.")
}
label = Sys.getenv("EDI_LABEL", unset = "portable")
case_name = Sys.getenv("EDI_CASE", unset = "ordinal_ppo")
num_cores = as.integer(Sys.getenv("EDI_NUM_CORES", unset = "1"))
r = as.integer(Sys.getenv("EDI_R", unset = "201"))
reps = as.integer(Sys.getenv("EDI_REPS", unset = "3"))
force_mirai = identical(tolower(Sys.getenv("EDI_FORCE_MIRAI", unset = "false")), "true")

.libPaths(c(lib_dir, .libPaths()))
library(EDI)

stable_seed = function(...) {
    as.integer(sum(utf8ToInt(paste(..., collapse = "|"))) %% 2147483647L)
}

make_case = function(case_name, seed_key) {
    if (identical(case_name, "ordinal_ppo")) {
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

        set.seed(stable_seed("ordinal_ppo_build_benchmark", seed_key))
        des_obj = DesignSeqOneByOneKK14$new(response_type = "ordinal", n = nrow(X_design))
        for (t in seq_len(nrow(X_design))) {
            w_t = des_obj$add_one_subject_to_experiment_and_assign(X_design[t, , drop = FALSE])
            y_t = pmax(1L, as.integer(y_base[t] + (if (w_t == 1) beta_T else 0) + rnorm(1, 0, SD_NOISE)))
            des_obj$add_one_subject_response(t, y_t, dead = 1)
        }
        return(list(
            des_obj = des_obj,
            inference = function() InferenceOrdinalMultiPartialProportionalOddsRegr$new(
                des_obj,
                nonparallel = character(0),
                verbose = FALSE
            )
        ))
    }

    if (identical(case_name, "kk_compound_continuous")) {
        beta_T = 0.75
        SD_NOISE = 0.5
        max_n = 100L
        X_design = na.omit(as.data.frame(airquality))
        X_design = X_design[seq_len(max_n), c("Ozone", "Solar.R", "Wind", "Temp")]
        set.seed(stable_seed("kk_compound_continuous_build_benchmark", seed_key))
        des_obj = DesignSeqOneByOneKK14$new(response_type = "continuous", n = nrow(X_design))
        for (t in seq_len(nrow(X_design))) {
            w_t = des_obj$add_one_subject_to_experiment_and_assign(X_design[t, , drop = FALSE])
            y_t = X_design$Temp[t] + if (w_t == 1) beta_T else 0 + rnorm(1, 0, SD_NOISE)
            des_obj$add_one_subject_response(t, y_t, dead = 1)
        }
        return(list(
            des_obj = des_obj,
            inference = function() InferenceAllKKCompoundMeanDiff$new(des_obj, verbose = FALSE)
        ))
    }

    if (identical(case_name, "proportion_fractional_logit")) {
        x_dat = data.frame(
            x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
            x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
        )
        y_dat = c(0.10, 0.25, 0.20, 0.40, 0.35, 0.55, 0.60, 0.75)
        set.seed(stable_seed("proportion_fractional_logit_build_benchmark", seed_key))
        des_obj = DesignSeqOneByOneBernoulli$new(n = nrow(x_dat), response_type = "proportion", verbose = FALSE)
        for (i in seq_len(nrow(x_dat))) {
            des_obj$add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
        }
        des_obj$add_all_subject_responses(y_dat)
        return(list(
            des_obj = des_obj,
            inference = function() InferencePropMultiFractionalLogit$new(des_obj, verbose = FALSE)
        ))
    }

    if (identical(case_name, "count_poisson")) {
        x_dat = data.frame(
            x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
            x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
        )
        y_dat = c(0L, 1L, 1L, 2L, 2L, 3L, 3L, 4L)
        set.seed(stable_seed("count_poisson_build_benchmark", seed_key))
        des_obj = DesignSeqOneByOneBernoulli$new(n = nrow(x_dat), response_type = "count", verbose = FALSE)
        for (i in seq_len(nrow(x_dat))) {
            des_obj$add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
        }
        des_obj$add_all_subject_responses(y_dat)
        return(list(
            des_obj = des_obj,
            inference = function() InferenceCountUnivPoissonRegr$new(des_obj, verbose = FALSE)
        ))
    }

    stop("Unsupported case: ", case_name)
}

set_num_cores(num_cores, force_mirai = force_mirai)
on.exit(unset_num_cores(), add = TRUE)

timings = numeric(reps)
cat(sprintf(
    "Benchmark: %s | case=%s | num_cores=%d | force_mirai=%s | r=%d | reps=%d\n",
    label, case_name, num_cores, force_mirai, r, reps
))

for (i in seq_len(reps)) {
    bench = make_case(case_name, paste(label, i, num_cores, r, sep = "|"))
    des_obj = bench$des_obj
    inf_obj = bench$inference()
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
