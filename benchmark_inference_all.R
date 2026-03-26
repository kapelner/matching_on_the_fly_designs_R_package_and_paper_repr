
# Set library path
.libPaths(c("local_R_lib", .libPaths()))

# Load necessary packages
packages = c("SeqExpMatch", "microbenchmark", "data.table", "dplyr", "PTE", "mlbench", "AppliedPredictiveModeling", "qgam")
for (p in packages) {
    if (!require(p, character.only = TRUE)) {
        install.packages(p, repos = "https://cloud.r-project.org/", lib = "local_R_lib")
        library(p, character.only = TRUE)
    }
}

# Load datasets and setup as in simple_tests.R
max_n_dataset = 100 
source("package_tests/_dataset_load.R")
MAX_N_DATASET = 100

custom_rand_stat = function(){
    yTs = private$seq_des_obj_priv_int$y[private$seq_des_obj_priv_int$w == 1]
    yCs = private$seq_des_obj_priv_int$y[private$seq_des_obj_priv_int$w == 0]
    if (length(yTs) < 2 || length(yCs) < 2) return(0)
    (mean(yTs) - mean(yCs)) / sqrt(var(yTs) / length(yTs) + var(yCs) / length(yCs))
}

CORE_COUNTS = c(5, 2, 1) # User requested order
r = 201 
pval_epsilon = 0.05
beta_T = 0.2
SD_NOISE = 0.1

prepare_seq_des_obj = function(response_type, dataset_name = "airquality", design_class = SeqDesignKK14) {
    D = datasets_and_response_models[[dataset_name]]
    X_design = as.data.frame(D$X)
    n = min(nrow(X_design), MAX_N_DATASET)
    y = D$y_original[[response_type]]
    
    seq_des_obj = design_class$new(response_type = response_type, n = n)
    
    for (t in 1:n) {
        w_t = seq_des_obj$add_subject_to_experiment_and_assign(X_design[t, , drop = FALSE])
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
        seq_des_obj$add_subject_response(t, y_t, dead = 1)
    }
    return(seq_des_obj)
}

namespace_content = readLines("SeqExpMatch/NAMESPACE")
inf_classes = grep("^export\\(DesignInference", namespace_content, value = TRUE)
inf_classes = gsub("export\\(|\\)", "", inf_classes)

categorized_classes = list(
    continuous = inf_classes[grepl("Contin|All|Bai", inf_classes)],
    incidence = inf_classes[grepl("Incid|All", inf_classes)],
    proportion = inf_classes[grepl("Prop|All", inf_classes)],
    count = inf_classes[grepl("Count|All", inf_classes)],
    survival = inf_classes[grepl("Survival|All", inf_classes)],
    ordinal = inf_classes[grepl("Ordinal|All", inf_classes)]
)

results_file = "benchmark_inference_results.csv"
if (file.exists(results_file)) file.remove(results_file)

bm_safe = function(label, expr, env = parent.frame()) {
    cat(sprintf("%s\n", label))
    flush.console()
    t_start = proc.time()[["elapsed"]]
    res = tryCatch({
        eval(expr, envir = env)
        
        t_end = proc.time()[["elapsed"]]
        duration = round(t_end - t_start, 3)
        cat(sprintf(" (%.3fs) ", duration))
        duration
    }, error = function(e) {
        cat(sprintf(" (ERROR: %s) ", e$message))
        NA
    })
    res
}

for (response_type in names(categorized_classes)) {
    cat(sprintf("\nResponse type: %s\n", response_type))
    flush.console()
    if (!(response_type %in% names(datasets_and_response_models$airquality$y_original))) next
    
    seq_des_obj_kk14 = prepare_seq_des_obj(response_type, design_class = SeqDesignKK14)
    seq_des_obj_kk21 = prepare_seq_des_obj(response_type, design_class = SeqDesignKK21)
    
    for (inf_class_name in unique(categorized_classes[[response_type]])) {
        if (inf_class_name == "DesignInferenceAllKKCompoundMeanDiff") next
        
        cat(sprintf("  %s:\n", inf_class_name))
        flush.console()
        
        inf_class = tryCatch(get(inf_class_name), error = function(e) NULL)
        if (is.null(inf_class)) {
            cat("    FAILED (not found)\n")
            next
        }
        
        # Decide which design object to use
        current_seq_des_obj = if (grepl("TKK21", inf_class_name)) seq_des_obj_kk21 else seq_des_obj_kk14
        
        # Skip slow randomization CI for heavy models
        is_heavy_model = grepl("GLMM|GEE|QuantileRegr", inf_class_name)
        if (response_type != "continuous" && grepl("Wilcox", inf_class_name)) {
            is_heavy_model = TRUE
        }
        
        for (num_cores in CORE_COUNTS) {
            cat(sprintf("    num_cores = %d: ", num_cores))
            flush.console()
            
            # 1. boot (New object)
            inf_obj = tryCatch(inf_class$new(current_seq_des_obj, num_cores = num_cores, verbose = FALSE), error = function(e) NULL)
            boot_time = if (!is.null(inf_obj)) {
                bm_safe("boot", quote(inf_obj$compute_bootstrap_two_sided_pval(B = r, na.rm = TRUE)))
            } else NA
            
            # 2. ci (New object; save w before in case draw_ws_according_to_design corrupts it)
            w_original = current_seq_des_obj$.__enclos_env__$private$w
            inf_obj = tryCatch(inf_class$new(current_seq_des_obj, num_cores = num_cores, verbose = FALSE), error = function(e) NULL)
            ci_time = if (!is.null(inf_obj) && !is_heavy_model) {
                bm_safe("ci", quote(inf_obj$compute_confidence_interval_rand(r = r, pval_epsilon = pval_epsilon, show_progress = FALSE)))
            } else NA

            # 3. custom_ci (New object, reset design cache for cold measurement)
            current_seq_des_obj$.__enclos_env__$private$X = NULL
            current_seq_des_obj$.__enclos_env__$private$w = w_original
            inf_obj = tryCatch(inf_class$new(current_seq_des_obj, num_cores = num_cores, verbose = FALSE), error = function(e) NULL)
            custom_ci_time = if (!is.null(inf_obj) && !is_heavy_model) {
                tryCatch(inf_obj$set_custom_randomization_statistic_function(custom_rand_stat), error = function(e) NULL)
                bm_safe("custom_ci", quote(inf_obj$compute_confidence_interval_rand(r = r, pval_epsilon = pval_epsilon, show_progress = FALSE)))
            } else NA
            
            # 4. rand (New object)
            inf_obj = tryCatch(inf_class$new(current_seq_des_obj, num_cores = num_cores, verbose = FALSE), error = function(e) NULL)
            rand_time = if (!is.null(inf_obj) && !is_heavy_model) {
                bm_safe("rand", quote(inf_obj$compute_two_sided_pval_for_treatment_effect_rand(r = r, show_progress = FALSE)))
            } else NA
            
            # 5. ci_after_rand (REUSE object from rand)
            ci_after_rand_time = if (!is.null(inf_obj) && !is_heavy_model) {
                bm_safe("ci_after_rand", quote(inf_obj$compute_confidence_interval_rand(r = r, pval_epsilon = pval_epsilon, show_progress = FALSE)))
            } else NA
            
            # 6. asymp (Asymptotic p-value/CI)
            inf_obj = tryCatch(inf_class$new(current_seq_des_obj, num_cores = num_cores, verbose = FALSE), error = function(e) NULL)
            asymp_pval_time = if (!is.null(inf_obj)) {
                bm_safe("asymp_pval", quote(inf_obj$compute_asymp_two_sided_pval_for_treatment_effect()))
            } else NA
            asymp_ci_time = if (!is.null(inf_obj)) {
                bm_safe("asymp_ci", quote(inf_obj$compute_asymp_confidence_interval()))
            } else NA

            row = data.table(
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
    ), by = .(num_cores)]
    print(summary_table)
}
