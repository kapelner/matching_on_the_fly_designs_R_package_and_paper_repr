pacman::p_load(data.table, R.utils, ggplot2)

Nrep = 3000
tmpdir = tempfile()
dir.create(tmpdir)
results_file = sprintf("cmh_exact_sims_results_Nrep_%d.tar.bz2", Nrep)
utils::untar(results_file, exdir = tmpdir)

csv_file = file.path(tmpdir, sprintf("cmh_exact_sims_results_Nrep_%d.csv", Nrep))
raw_results_dt = data.table::fread(csv_file)
raw_results_dt[, reject := pval < 0.05]
raw_results_dt[, covers := ci_lo <= true_estimand & true_estimand <= ci_hi]
raw_results_dt[, ci_length := ci_hi - ci_lo]
results_dt = raw_results_dt[,
  .(
    pow_avg = mean(reject, na.rm = TRUE), 
    n_pow = sum(!is.na(reject)),
    cov_avg = mean(covers, na.rm = TRUE), 
    n_cov = sum(!is.na(covers)),
    len_avg = mean(ci_length, na.rm = TRUE), 
    len_sd = sd(ci_length, na.rm = TRUE), 
    n_len = sum(!is.na(ci_length))
  ),
  by = c("n", "p", "betaT", "design", "inference", "cond_exp_func_model")                         
]
#table(results_dt$design, results_dt$inference)
#table(results_dt$n, results_dt$p, results_dt$betaT)

# Wald CIs for power and coverage
z = qnorm(0.975)
results_dt[, `:=`(
  pow_a = pow_avg - z * sqrt(pow_avg * (1 - pow_avg) / n_pow),
  pow_b = pow_avg + z * sqrt(pow_avg * (1 - pow_avg) / n_pow),
  cov_a = cov_avg - z * sqrt(cov_avg * (1 - cov_avg) / n_cov),
  cov_b = cov_avg + z * sqrt(cov_avg * (1 - cov_avg) / n_cov),
  len_a = len_avg - z * len_sd / sqrt(n_len),
  len_b = len_avg + z * len_sd / sqrt(n_len)
)]

# Shorten labels
results_dt[, design_short    := gsub("FixedDesign", "", design)]
results_dt[, inference_short := inference]
results_dt[, inference_short := gsub("InferenceIncidence", "", inference_short)]
results_dt[, inference_short := gsub("InferenceIncid",     "", inference_short)]

# Row 1: iBCRD, OptimalBlocks B=4/8/16/32; Row 2: BinaryMatch, Blocking B_target=4/8/16/32
design_levels = c(
  "iBCRD",                    
  "OptimalBlocks (B=4)",     
  "OptimalBlocks (B=8)",   
  "OptimalBlocks (B=16)",  
  "OptimalBlocks (B=32)",   
  "BinaryMatch",  
  "Blocking (B_target=4)",  
  "Blocking (B_target=8)",  
  "Blocking (B_target=16)",    
  "Blocking (B_target=32)"
)
results_dt[, design_short := factor(design_short, levels = design_levels)]

make_plot = function(dat, metric, ylab, hline, title_str, filename_stem, plot = TRUE, save_PDF = FALSE, save_PNG = FALSE, free_y = FALSE, hline_linetype = "dashed") {
  dat = copy(dat)
  dat[, y := get(metric)]

  has_ci_cols = paste0(substr(metric, 1, 3), "_a") %in% names(dat)
  if (has_ci_cols) {
    dat[, y_lo := get(paste0(substr(metric, 1, 3), "_a"))]
    dat[, y_hi := get(paste0(substr(metric, 1, 3), "_b"))]
  }

  dodge = position_dodge(width = 0.4)
  plot_obj = ggplot(dat, aes(x = factor(n), y = y,
                      color = inference_short, group = inference_short)) +
    geom_point(size = 1.5, position = dodge) +
    facet_wrap(~ design_short, nrow = 2, ncol = 5, scales = if (free_y) "free_y" else "fixed") +
    scale_color_brewer(palette = "Set1") +
    labs(x = "n", y = ylab, color = "Inference", title = title_str) +
    theme_bw(base_size = 9) +
    theme(
      strip.text       = element_text(size = 10),
      legend.position  = "bottom",
      legend.text      = element_text(size = 12),
      legend.title     = element_text(size = 12),
      axis.title       = element_text(size = 13),
      axis.text        = element_text(size = 11),
      panel.grid.minor = element_blank()
    )
  if (!is.null(hline)) {
    plot_obj = plot_obj + geom_hline(yintercept = hline, linetype = hline_linetype, color = "gray60", linewidth = 0.4)
  }
  if (has_ci_cols) {
    plot_obj = plot_obj + geom_errorbar(aes(ymin = y_lo, ymax = y_hi), width = 0.2, linewidth = 0.4, position = dodge)
  }
  if (plot){
      plot(plot_obj)
  }
  if(save_PDF){
    ggsave(sprintf("%s.pdf", filename_stem), plot_obj, width = 12, height = 6)
  }
  if (save_PNG){
    ggsave(sprintf("%s.png", filename_stem), plot_obj, width = 12, height = 6, dpi = 150)
  }
  invisible(plot_obj)
}

Nrep = max(results_dt$n_pow, na.rm = TRUE)

for (p_ in unique(results_dt$p)) {
  for (dt_val in unique(results_dt$cond_exp_func_model)) {
    sub = results_dt[p == p_ & cond_exp_func_model == dt_val]
    stem_pow = sprintf("plot_power_p%s_%s_Nrep_%d",     p_, dt_val, Nrep)
    stem_cov = sprintf("plot_coverage_p%s_%s_Nrep_%d",  p_, dt_val, Nrep)
    stem_len = sprintf("plot_ci_length_p%s_%s_Nrep_%d", p_, dt_val, Nrep)
    stem_size = sprintf("plot_size_p%s_%s_Nrep_%d",     p_, dt_val, Nrep)
    if (nrow(sub[betaT == 1]) == 0L) next
    make_plot(sub[betaT == 1], "pow_avg",     "Power",    0.05,
              sprintf("Power (betaT=1) | p=%s, log_odds_model=%s",    p_, dt_val), stem_pow,  save_PDF = TRUE, plot = FALSE)
    make_plot(sub[betaT == 1], "cov_avg",  "Coverage", 0.95,
              sprintf("Coverage (betaT=1) | p=%s, log_odds_model=%s", p_, dt_val), stem_cov,  save_PDF = TRUE, plot = FALSE)
    make_plot(sub[betaT == 1], "len_avg", "CI Length", NULL,
              sprintf("CI Length (betaT=1) | p=%s, log_odds_model=%s", p_, dt_val), stem_len, save_PDF = TRUE, plot = FALSE, free_y = TRUE)
    if (nrow(sub[betaT == 0]) == 0L) next
    make_plot(sub[betaT == 0], "pow_avg",     "Size",     0.05,
              sprintf("Size (betaT=0) | p=%s, log_odds_model=%s",     p_, dt_val), stem_size, save_PDF = TRUE, plot = FALSE, hline_linetype = "dotted")
  }
}

