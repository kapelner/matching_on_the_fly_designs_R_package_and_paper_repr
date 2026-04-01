
if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The 'pkgload' package is required to benchmark the current EDI source tree.")
}
pkgload::load_all("SeqExpMatch", quiet = TRUE)
library(EDI)

# Load remaining packages
packages = c("microbenchmark", "data.table", "dplyr", "PTE", "mlbench", "AppliedPredictiveModeling", "qgam", "mirai")
for (p in packages) {
    if (!require(p, character.only = TRUE)) {
        install.packages(p, repos = "https://cloud.r-project.org/", lib = "local_R_lib")
        library(p, character.only = TRUE)
    }
}

# Load datasets and setup as in simple_tests.R
max_n_dataset = 100 
source("package_tests/_dataset_load.R")

CORE_COUNTS = c(5, 2, 1)
FORCE_MIRAI_VALUES = c(TRUE, FALSE)

# Pre-create global clusters to avoid repeated 300ms startup penalty
# We benchmark one core count at a time, so we create/stop for each iteration of num_cores
# but since the current script structure is nested differently, let's just 
# handle it inside the loop or create once for the max cores.
# Actually, the user wants the objects to automatically pick it up.
# So we will wrap the num_cores loop logic.

# NOTE: Since Inference objects are recreated for every task, we want a persistent 
# cluster for the duration of each (response_type, inf_class, num_cores) triplet.
# But even better, let's just create it once for the largest num_cores and reuse it.
# However, the benchmark needs to test 5, 2, and 1 cores specifically.
# So we will manage the global cluster inside the num_cores loop.

MAX_N_DATASET = 100

custom_rand_stat = function(){
    yTs = private$des_obj_priv_int$y[private$des_obj_priv_int$w == 1]
    yCs = private$des_obj_priv_int$y[private$des_obj_priv_int$w == 0]
    if (length(yTs) < 2 || length(yCs) < 2) return(0)
    (mean(yTs) - mean(yCs)) / sqrt(var(yTs) / length(yTs) + var(yCs) / length(yCs))
}

r = 201
pval_epsilon = 0.05
beta_T = 0.2
SD_NOISE = 0.1

prepare_des_obj = function(response_type, dataset_name = "airquality", design_class = DesignSeqOneByOneKK14) {
    D = datasets_and_response_models[[dataset_name]]
    X_design = as.data.frame(D$X)
    n = min(nrow(X_design), MAX_N_DATASET)
    y = D$y_original[[response_type]]
    
    des_obj = design_class$new(response_type = response_type, n = n)
    
    for (t in 1:n) {
        w_t = des_obj$add_subject_to_experiment_and_assign(X_design[t, , drop = FALSE])
        if (response_type == "continuous") {
            y_t = y[t] + (if(w_t == 1) beta_T else 0) + rnorm(1, 0, SD_NOISE)
        } else if (response_type == "incidence") {
            y_t = as.integer(plogis(qlogis(pmax(0.01, pmin(0.99, y[t]))) + (if(w_t == 1) beta_T else 0)) > runif(1))
        } else if (response_type == "proportion") {
            y_t = plogis(qlogis(pmax(0.01, pmin(0.99, y[t]))) + (if(w_t == 1) beta_T else 0))
        } else if (response_type == "count") {
            y_t = as.integer(rpois(1, exp(log(pmax(0.1, y[t])) + (if(w_t == 1) beta_T else 0))))
        } else if (response_type == "survival") {
            y_t = pmax(0.1, y[t] * exp(if(w_t == 1) beta_T else 0))
        } else if (response_type == "ordinal") {
            y_t = pmax(1, as.integer(y[t] + (if(w_t == 1) beta_T else 0) + rnorm(1, 0, SD_NOISE)))
        }
        des_obj$add_subject_response(t, y_t, dead = 1)
    }
    return(des_obj)
}

stable_seed = function(...) {
    as.integer(sum(utf8ToInt(paste(..., collapse = "|"))) %% 2147483647L)
}

prepare_cold_des_obj = function(response_type, design_class, seed_key) {
    set.seed(stable_seed("benchmark_inference_all", response_type, seed_key))
    prepare_des_obj(response_type, design_class = design_class)
}

dispatch_policy = function(inference_class, response_type, operation) {
    get("edi_parallel_dispatch_policy", envir = asNamespace("EDI"))(
        inference_class = inference_class,
        response_type = response_type,
        operation = operation
    )
}

benchmark_dispatch_skip = function(inference_class, response_type, operation, num_cores) {
    num_cores > 1L && isTRUE(dispatch_policy(inference_class, response_type, operation)$force_serial)
}

skip_benchmark_step = function(label, inference_class, response_type, operation, num_cores) {
    cat(sprintf("%s (SKIP: forced serial dispatch policy)\n", label))
    flush.console()
    NA_real_
}

namespace_content = readLines("SeqExpMatch/NAMESPACE")
inf_classes = grep("^export\\(Inference", namespace_content, value = TRUE)
inf_classes = gsub("export\\(|\\)", "", inf_classes)

is_all = grepl("^InferenceAll", inf_classes)
is_ordinal = grepl("Ordinal", inf_classes)
categorized_classes = list(
    continuous = inf_classes[grepl("Contin|Bai", inf_classes) | is_all],
    incidence = inf_classes[grepl("Incid", inf_classes) | is_all],
    proportion = inf_classes[(grepl("Prop", inf_classes) & !is_ordinal) | is_all],
    count = inf_classes[grepl("Count", inf_classes) | is_all],
    survival = inf_classes[grepl("Survival", inf_classes) | is_all],
    ordinal = inf_classes[is_ordinal | is_all]
)

results_file = "benchmark_inference_results.csv"
if (file.exists(results_file)) file.remove(results_file)

MAX_BM_SECS = 200

bm_safe = function(label, expr, env = parent.frame(), num_cores_to_restore = NULL, force_mirai_to_restore = FALSE) {
    cat(sprintf("%s\n", label))
    flush.console()
    t_start = proc.time()[["elapsed"]]
    main_pid = Sys.getpid()

    # Launch a watchdog that kills all children of this process after MAX_BM_SECS.
    # This handles the case where mclapply children ignore setTimeLimit.
    watchdog_script = sprintf(
        "sleep %d && pkill -P %d -KILL 2>/dev/null; true",
        MAX_BM_SECS, main_pid
    )
    # Redirect watchdog stdin/stdout/stderr away from the inherited pipe;
    # otherwise system(intern=TRUE) would block until the watchdog finishes.
    watchdog_pid_str = system(
        sprintf("bash -c '%s' </dev/null >/dev/null 2>&1 & echo $!", watchdog_script),
        intern = TRUE
    )
    watchdog_pid = suppressWarnings(as.integer(watchdog_pid_str))

    kill_watchdog = function() {
        if (!is.na(watchdog_pid) && watchdog_pid > 0)
            system(sprintf("kill %d 2>/dev/null", watchdog_pid), ignore.stdout = TRUE)
    }

    restore_parallel_state = function() {
        if (is.null(num_cores_to_restore) || is.na(num_cores_to_restore) || num_cores_to_restore <= 1L) {
            return(invisible(NULL))
        }
        try(unset_num_cores(), silent = TRUE)
        try(set_num_cores(as.integer(num_cores_to_restore), force_mirai = force_mirai_to_restore), silent = TRUE)
        invisible(NULL)
    }

    res = tryCatch({
        setTimeLimit(elapsed = MAX_BM_SECS, transient = TRUE)
        on.exit({ setTimeLimit(elapsed = Inf, transient = FALSE); kill_watchdog() }, add = TRUE)
        eval(expr, envir = env)
        setTimeLimit(elapsed = Inf, transient = FALSE)
        kill_watchdog()
        t_end = proc.time()[["elapsed"]]
        duration = round(t_end - t_start, 3)
        # C++ code bypasses setTimeLimit; check manually
        if (duration >= MAX_BM_SECS) {
            restore_parallel_state()
            cat(sprintf(" (TIMEOUT >%.1fs) ", duration))
            return(Inf)
        }
        cat(sprintf(" (%.3fs) ", duration))
        duration
    }, error = function(e) {
        kill_watchdog()
        t_elapsed = round(proc.time()[["elapsed"]] - t_start, 3)
        restore_parallel_state()
        if (t_elapsed >= MAX_BM_SECS - 1 ||
            grepl("reached elapsed time limit", e$message, fixed = TRUE)) {
            cat(sprintf(" (TIMEOUT >%.1fs) ", t_elapsed))
            Inf
        } else {
            cat(sprintf(" (ERROR: %s) ", e$message))
            NA
        }
    })
    res
}

for (response_type in names(categorized_classes)) {
    cat(sprintf("\nResponse type: %s\n", response_type))
    flush.console()
    if (!(response_type %in% names(datasets_and_response_models$airquality$y_original))) next
    
    des_obj_kk14 = prepare_des_obj(response_type, design_class = DesignSeqOneByOneKK14)
    des_obj_kk21 = prepare_des_obj(response_type, design_class = DesignSeqOneByOneKK21)
    
    for (inf_class_name in unique(categorized_classes[[response_type]])) {
        if (inf_class_name == "InferenceAllKKCompoundMeanDiff") next
        
        cat(sprintf("  %s:\n", inf_class_name))
        flush.console()
        
        inf_class = tryCatch(get(inf_class_name), error = function(e) NULL)
        if (is.null(inf_class)) {
            cat("    FAILED (not found)\n")
            next
        }
        
        # Decide which design object to use
        current_design_class = if (grepl("TKK21", inf_class_name)) DesignSeqOneByOneKK21 else DesignSeqOneByOneKK14
        current_design_label = if (identical(current_design_class, DesignSeqOneByOneKK21)) "KK21" else "KK14"
        current_des_obj = if (identical(current_design_class, DesignSeqOneByOneKK21)) des_obj_kk21 else des_obj_kk14
        
        # Skip slow randomization CI for heavy models
        is_heavy_model = grepl("GLMM|GEE|QuantileRegr", inf_class_name)
        if (response_type != "continuous" && grepl("Wilcox", inf_class_name)) {
            is_heavy_model = TRUE
        }

        # Pre-screen with num_cores=1: skip if all operations time out (class too slow to benchmark)
        cat("    [pre-screen 1 core]: ")
        flush.console()
        inf_obj_screen = tryCatch(inf_class$new(current_des_obj, verbose = FALSE), error = function(e) NULL)
        screen_times = list(
            boot = if (!is.null(inf_obj_screen))
                bm_safe("boot", quote(inf_obj_screen$compute_bootstrap_two_sided_pval(B = r, na.rm = TRUE)))
                else NA,
            rand = if (!is.null(inf_obj_screen) && !is_heavy_model)
                bm_safe("rand", quote(inf_obj_screen$compute_two_sided_pval_for_treatment_effect_rand(r = r, show_progress = FALSE)))
                else NA
        )
        n_inf = sum(sapply(screen_times, function(x) !is.na(x) && is.finite(x)))
        if (n_inf == 0) {
            cat("  SKIP (all ops >10s at 1 core)\n")
            flush.console()
            next
        }
        cat("  OK\n")
        flush.console()

        for (force_mirai in FORCE_MIRAI_VALUES) {
            cat(sprintf("    force_mirai = %s\n", force_mirai))
            flush.console()
            for (num_cores in CORE_COUNTS) {
                cat(sprintf("      num_cores = %d: ", num_cores))
                flush.console()

                # Set the global number of cores
                set_num_cores(num_cores, force_mirai = force_mirai)

                # 1. boot (New cold design object and new inference object)
                boot_des_obj = prepare_cold_des_obj(
                    response_type,
                    current_design_class,
                    paste(response_type, current_design_label, inf_class_name, force_mirai, num_cores, "boot_cold", sep = "|")
                )
                inf_obj = tryCatch(inf_class$new(boot_des_obj, verbose = FALSE), error = function(e) NULL)
                boot_time = if (!is.null(inf_obj) && benchmark_dispatch_skip(inf_class_name, response_type, "bootstrap", num_cores)) {
                    skip_benchmark_step("boot", inf_class_name, response_type, "bootstrap", num_cores)
                } else if (!is.null(inf_obj)) {
                    bm_safe("boot", quote(inf_obj$compute_bootstrap_two_sided_pval(B = r, na.rm = TRUE)),
                            num_cores_to_restore = num_cores, force_mirai_to_restore = force_mirai)
                } else NA

                # 2. ci (New cold design object and new inference object)
                ci_des_obj = prepare_cold_des_obj(
                    response_type,
                    current_design_class,
                    paste(response_type, current_design_label, inf_class_name, force_mirai, num_cores, "ci_cold", sep = "|")
                )
                inf_obj = tryCatch(inf_class$new(ci_des_obj, verbose = FALSE), error = function(e) NULL)
                ci_time = if (!is.null(inf_obj) && !is_heavy_model &&
                               benchmark_dispatch_skip(inf_class_name, response_type, "rand_ci", num_cores)) {
                    skip_benchmark_step("ci", inf_class_name, response_type, "rand_ci", num_cores)
                } else if (!is.null(inf_obj) && !is_heavy_model) {
                    bm_safe("ci", quote(inf_obj$compute_confidence_interval_rand(r = r, pval_epsilon = pval_epsilon, show_progress = FALSE)),
                            num_cores_to_restore = num_cores, force_mirai_to_restore = force_mirai)
                } else NA

                # 3. custom_ci (New cold design object and new inference object)
                custom_ci_des_obj = prepare_cold_des_obj(
                    response_type,
                    current_design_class,
                    paste(response_type, current_design_label, inf_class_name, force_mirai, num_cores, "custom_ci_cold", sep = "|")
                )
                inf_obj = tryCatch(inf_class$new(custom_ci_des_obj, verbose = FALSE), error = function(e) NULL)
                custom_ci_time = if (!is.null(inf_obj) && !is_heavy_model &&
                                      benchmark_dispatch_skip(inf_class_name, response_type, "rand_ci", num_cores)) {
                    skip_benchmark_step("custom_ci", inf_class_name, response_type, "rand_ci", num_cores)
                } else if (!is.null(inf_obj) && !is_heavy_model) {
                    tryCatch(inf_obj$set_custom_randomization_statistic_function(custom_rand_stat), error = function(e) NULL)
                    bm_safe("custom_ci", quote(inf_obj$compute_confidence_interval_rand(r = r, pval_epsilon = pval_epsilon, show_progress = FALSE)),
                            num_cores_to_restore = num_cores, force_mirai_to_restore = force_mirai)
                } else NA

                # 4. rand (New cold design object and new inference object)
                rand_des_obj = prepare_cold_des_obj(
                    response_type,
                    current_design_class,
                    paste(response_type, current_design_label, inf_class_name, force_mirai, num_cores, "rand_cold", sep = "|")
                )
                inf_obj = tryCatch(inf_class$new(rand_des_obj, verbose = FALSE), error = function(e) NULL)
                rand_time = if (!is.null(inf_obj) && !is_heavy_model &&
                                 benchmark_dispatch_skip(inf_class_name, response_type, "rand_pval", num_cores)) {
                    skip_benchmark_step("rand", inf_class_name, response_type, "rand_pval", num_cores)
                } else if (!is.null(inf_obj) && !is_heavy_model) {
                    bm_safe("rand", quote(inf_obj$compute_two_sided_pval_for_treatment_effect_rand(r = r, show_progress = FALSE)),
                            num_cores_to_restore = num_cores, force_mirai_to_restore = force_mirai)
                } else NA

                # 5. ci_after_rand (REUSE warm design and inference object from rand)
                ci_after_rand_time = if (!is.null(inf_obj) && !is_heavy_model &&
                                          benchmark_dispatch_skip(inf_class_name, response_type, "rand_ci", num_cores)) {
                    skip_benchmark_step("ci_after_rand", inf_class_name, response_type, "rand_ci", num_cores)
                } else if (!is.null(inf_obj) && !is_heavy_model) {
                    bm_safe("ci_after_rand", quote(inf_obj$compute_confidence_interval_rand(r = r, pval_epsilon = pval_epsilon, show_progress = FALSE)),
                            num_cores_to_restore = num_cores, force_mirai_to_restore = force_mirai)
                } else NA

                # 6. asymp (Asymptotic p-value/CI)
                # Benchmark each on a fresh design object as well as a fresh inference object
                # so neither inference caches nor design-side all-subject-data caches are warm.
                asymp_seed_key = paste(
                    response_type, current_design_label, inf_class_name,
                    force_mirai, num_cores, "asymp_cold",
                    sep = "|"
                )
                asymp_des_obj = prepare_cold_des_obj(response_type, current_design_class, asymp_seed_key)
                inf_obj = tryCatch(inf_class$new(asymp_des_obj, verbose = FALSE), error = function(e) NULL)
                asymp_pval_time = if (!is.null(inf_obj)) {
                    bm_safe("asymp_pval", quote(inf_obj$compute_asymp_two_sided_pval_for_treatment_effect()),
                            num_cores_to_restore = num_cores, force_mirai_to_restore = force_mirai)
                } else NA
                asymp_des_obj = prepare_cold_des_obj(response_type, current_design_class, asymp_seed_key)
                inf_obj = tryCatch(inf_class$new(asymp_des_obj, verbose = FALSE), error = function(e) NULL)
                asymp_ci_time = if (!is.null(inf_obj)) {
                    bm_safe("asymp_ci", quote(inf_obj$compute_asymp_confidence_interval()),
                            num_cores_to_restore = num_cores, force_mirai_to_restore = force_mirai)
                } else NA

                # Unset cores
                unset_num_cores()

                row = data.table(
                    force_mirai = force_mirai,
                    num_cores = num_cores,
                    response_type = response_type,
                    inference_class = inf_class_name,
                    bootstrap_pval_time = boot_time,
                    rand_pval_time = rand_time,
                    rand_ci_time = ci_time,
                    rand_custom_ci_time = custom_ci_time,
                    rand_ci_after_rand_time = ci_after_rand_time,
                    asymp_pval_time = asymp_pval_time,
                    asymp_ci_time = asymp_ci_time
                )

                fwrite(row, results_file, append = file.exists(results_file))
                cat("done\n")
                flush.console()
            }
        }
    }
}

if (file.exists(results_file)) {
    benchmark_results = fread(results_file)
    summary_table = benchmark_results[, .(
        avg_bootstrap = mean(bootstrap_pval_time, na.rm = TRUE),
        avg_rand_pval = mean(rand_pval_time, na.rm = TRUE),
        avg_rand_ci = mean(rand_ci_time, na.rm = TRUE),
        avg_rand_custom_ci = mean(rand_custom_ci_time, na.rm = TRUE),
        avg_ci_after_rand = mean(rand_ci_after_rand_time, na.rm = TRUE),
        avg_asymp_pval = mean(asymp_pval_time, na.rm = TRUE),
        avg_asymp_ci = mean(asymp_ci_time, na.rm = TRUE)
    ), by = .(force_mirai, num_cores)]
    print(summary_table)
}
