library(data.table)
library(ggplot2)

results_dt = fread("simulations/azriel_exact_sims_results.csv")

# Wald CIs for power and coverage
z = qnorm(0.975)
results_dt[, `:=`(
  pow_a = power - z * sqrt(power * (1 - power) / n_pow),
  pow_b = power + z * sqrt(power * (1 - power) / n_pow),
  cov_a = coverage - z * sqrt(coverage * (1 - coverage) / n_cov),
  cov_b = coverage + z * sqrt(coverage * (1 - coverage) / n_cov)
)]

# Shorten labels
results_dt[, design_short    := gsub("FixedDesign", "", design)]
results_dt[, inference_short := inference]
results_dt[, inference_short := gsub("InferenceIncidence", "", inference_short)]
results_dt[, inference_short := gsub("InferenceIncid",     "", inference_short)]

# Row 1: iBCRD, OptimalBlocks B=4/8/16/32; Row 2: BinaryMatch, Blocking B_preferred=4/8/16/32
design_levels = c(
  "iBCRD",                    
  "OptimalBlocks (B=4)",     
  "OptimalBlocks (B=8)",   
  "OptimalBlocks (B=16)",  
  "OptimalBlocks (B=32)",   
  "BinaryMatch",  
  "Blocking (B_preferred=4)",  
  "Blocking (B_preferred=8)",  
  "Blocking (B_preferred=16)",    
  "Blocking (B_preferred=32)"
)
results_dt[, design_short := factor(design_short, levels = design_levels)]

make_plot = function(dat, metric, ylab, hline, title_str, filename_stem, plot = TRUE, save_PDF = FALSE, save_PNG = FALSE) {
  dat = copy(dat)
  dat[, y    := get(metric)]
  dat[, y_lo := get(paste0(substr(metric, 1, 3), "_a"))]
  dat[, y_hi := get(paste0(substr(metric, 1, 3), "_b"))]

  dodge = position_dodge(width = 0.4)
  p = ggplot(dat, aes(x = factor(n), y = y,
                      color = inference_short, group = inference_short)) +
    geom_hline(yintercept = hline, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    geom_point(size = 1.5, position = dodge) +
    geom_errorbar(aes(ymin = y_lo, ymax = y_hi), width = 0.2, linewidth = 0.4, position = dodge) +
    facet_wrap(~ design_short, nrow = 2, ncol = 5) +
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
  if (plot){
      plot(p)
  } 
  if(save_PDF){
    ggsave(sprintf("simulations/%s.pdf", filename_stem), p, width = 12, height = 6)
  }
  if (save_PNG){
    ggsave(sprintf("simulations/%s.png", filename_stem), p, width = 12, height = 6, dpi = 150)
  }
  invisible(p)
}

Nrep = max(results_dt$n_pow, na.rm = TRUE)

for (p_val in unique(results_dt$p)) {
  for (dt_val in unique(results_dt$data_type)) {
    sub = results_dt[p == p_val & data_type == dt_val]
    if (nrow(sub) == 0L) next
    stem_pow = sprintf("plot_power_p%s_%s_Nrep_%d",    p_val, dt_val, Nrep)
    stem_cov = sprintf("plot_coverage_p%s_%s_Nrep_%d", p_val, dt_val, Nrep)
    make_plot(sub[betaT == 1], "power",    "Power",    0.05,
              sprintf("Power (betaT=1) | p=%s, log_odds_model=%s",    p_val, dt_val), stem_pow, save_PDF = TRUE, plot = FALSE)
    make_plot(sub[betaT == 0], "coverage", "Coverage", 0.95,
              sprintf("Coverage (betaT=1) | p=%s, log_odds_model=%s", p_val, dt_val), stem_cov, save_PDF = TRUE, plot = FALSE)
  }
}
