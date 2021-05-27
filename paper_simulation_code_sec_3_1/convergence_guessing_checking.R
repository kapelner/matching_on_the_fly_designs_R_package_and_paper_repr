
# lagged_indic_T = c(0, indic_T[1 : (n-1)])
# lagged_tx_prop = cumsum(lagged_indic_T) / (0 : (n - 1))
# decision = c("guess", ifelse(lagged_tx_prop[2 : n] == (1 : (n - 1)) / 2, "guess", ifelse(lagged_tx_prop[2 : n] < (1 : (n - 1)) / 2, "T", "C")))
# 
# 
# data.frame(
#   # cumsum_indicT = cumsum(indic_T), 
#   cumsum_indicT_lag = cumsum(lagged_indic_T), 
#   indic_T = indic_T, 
#   lagged_tx_prop= lagged_tx_prop, 
#   # decision = decision, 
#   indic_T_guesses = indic_T_guesses
# )

# 
# indic_T[lagged_tx_prop == 0.5]
# indic_T_guesses[lagged_tx_prop == 0.5]
# 
# sum(indic_T == indic_T_guesses)
guessing_res = rbind(
		master_results$conv_guessing_strategy_pct_correct_avg, 
		master_results$conv_guessing_strategy_pct_correct_sd / sqrt(Nsim_per_block),
		master_results$conv_guessing_strategy_pct_correct_after_n_0_avg, 
		master_results$conv_guessing_strategy_pct_correct_after_n_0_sd / sqrt(Nsim_per_block)
)
colnames(guessing_res) = master_results$algorithm
rownames(guessing_res) = c("tot_avg", "tot_se_avg", "matching_avg", "matching_se_avg")
guessing_res
