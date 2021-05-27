NUM_PLOTS_PER_RESP = 1
SIZE_OF_PLOT_PTS = 0.4
plot_response_model = function(nsim, betaT, final_reservoir_size){
	for (j in 1 : p){
		pdf(file = paste("response_", randomization_type, "_", response_model, "_n_", n, "_nsim_", nsim, "_alpha_", prob_match_cutoff_lambda, "_x_", j, ".pdf", sep = ""))
		
		plot(Xy[, j], Xy$y, type = "n", 
				main = paste("x vs y b:", round(betaT, 2), 
						"\nn =", n, "alpha =", prob_match_cutoff_lambda, "response:", response_model,
						"\nfinal reservoir size:", final_reservoir_size), ylab = "y", xlab = paste("x_", j, " (trt = red, con = blue)", sep = "")
		)
		#if we're pair matching...
		if (length(grep("pair", randomization_type)) > 0){
			#reservoir gets dots
			reservoir = Xy[Xy$match_indic == 0, ]
			reservoirT = reservoir[reservoir$indic_T == 1, ]
			reservoirC = reservoir[reservoir$indic_T == 0, ]
			points(reservoirT[, j], reservoirT$y, col = "red", pch = "+", cex = SIZE_OF_PLOT_PTS)
			points(reservoirC[, j], reservoirC$y, col = "blue", pch = "+", cex = SIZE_OF_PLOT_PTS)
			for (match_index in 1 : max(Xy$match_indic)){
				matched_pair = Xy[Xy$match_indic == match_index, ]
				matched_pairT = matched_pair[matched_pair$indic_T == 1, ]
				matched_pairC = matched_pair[matched_pair$indic_T == 0, ]
				letter = c(letters, LETTERS, "<", ">", "[", "]", "{", "}", "|", "(", ")", "=", "$", "%", "^", "&", "*", "\\", "/", ";")[match_index %% 70 + 1]
				points(matched_pairT[, j], matched_pairT$y, col = "red", pch = letter, cex = SIZE_OF_PLOT_PTS)
				points(matched_pairC[, j], matched_pairC$y, col = "blue", pch = letter, cex = SIZE_OF_PLOT_PTS)			
			}
		} else { 
			#ttest group - just reservoir
			XyT = Xy[Xy$indic_T == 1, ]
			XyC = Xy[Xy$indic_T == 0, ]		
			points(XyT[, j], XyT$y, col = "red", pch = "+", cex = SIZE_OF_PLOT_PTS)
			points(XyC[, j], XyC$y, col = "blue", pch = "+", cex = SIZE_OF_PLOT_PTS)	
		}
		
		
		dev.off()
	}
}

plot_balance_and_betas = function(randomization_type, n){
	pdf(file = paste("balances_", response_model, "_", randomization_type, "_n_", n, ".pdf", sep = ""))
	hist(balances, br = 50)
	dev.off()
	
	pdf(file = paste("beta_hat_Ts_", response_model, "_", randomization_type, "_n_", n, ".pdf", sep = ""))
	hist(beta_hat_T, br = 50)
	dev.off()
}