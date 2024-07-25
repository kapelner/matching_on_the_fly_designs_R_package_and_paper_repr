sim_results = read.csv("sim_results.csv")

#simulation constants
TEST_TYPES = c("Classic", "Linear", "Exact")
SCENARIOS = c("NL", "LI", "ZE")
ALLOCATION_METHODS = c("C", "E", "S", "M", "SM")
num_alloc_methods = length(ALLOCATION_METHODS)
Ns_in_experiment = c(50, 100, 200)
NUM_SIMs = 2000
lambda = 0.10

#convenience objects
powers = sim_results[sim_results$feature == "power", ]
balances = sim_results[sim_results$feature == "avg_max_std_diff_bal", ]
std_err_beta_Ts = sim_results[sim_results$feature == "std_err_beta_T", ]
sizes = sim_results[sim_results$feature == "power" & sim_results$beta_T == 0, ]
biases = sim_results[sim_results$feature == "avg_abs_bias" & sim_results$beta_T == 1 & sim_results$z_test == 1, ]


#### 
#### 
#### 
#### Regenerate Figure 1
#### 
#### 
#### 

draw_power_fig = function(lambda,
							pch_symbols = c(15, 17, 19),
							pch_symbols_z = c(0, 2, 1),
							cex_size = 2,
							make_pdfs = TRUE,
							margins_inner = c(2.5,2.5,0.1,0.7),
							margins = c(2.5,4.5,0.1,0.7),
							dot_size = 1.75,
							line_thickness = 4,
							allocation_methods = c("C", "E", "S", "M", "SM"), #ALLOCATION_METHODS = c("C", "E", "S", "M", "SPM", "SM")
							beta_T = 1){
	
	num_alloc_methods = length(allocation_methods)
	
	for (scenario in SCENARIOS){
		for (test_type in TEST_TYPES){
			if (test_type == "Classic"){
				par(mar = margins)
			} else {
				par(mar = margins_inner)
			}
			
			if (make_pdfs){
				pdf(paste("scenario_", scenario, "_test_type_", test_type, "_lambda_", gsub("\\.", "_", lambda), "_beta_T_", beta_T, ".PDF", sep = ""))
				par(mar = margins)
			} else {
				windows()
			}
			plot(1 : num_alloc_methods, rep(0, num_alloc_methods), ylim = c(0, 1), xlab = "", ylab = ifelse(test_type == "Classic", "Power", ""), cex.lab = cex_size, type = "n", xaxt = "n", cex.axis = cex_size)
			axis(1, at = 1 : num_alloc_methods, labels = allocation_methods, cex.axis = cex_size)
			
			for (allocation_method in allocation_methods){
				x_pos = which(allocation_method == allocation_methods)
				#get power
				# plot powers by n
				for (n_i in 1 : length(Ns_in_experiment)){
					power_runs = powers[powers$model_name == scenario & 
									powers$test_type == test_type & 
									powers$allocation_method == allocation_method &
									powers$beta_T == beta_T &
									powers$n == Ns_in_experiment[n_i], ]
					
					if (allocation_method %in% c("SM", "SPM", "CPM")){
						power_run = power_runs[power_runs$lambda == lambda & power_runs$z_test == 1, "val"]
					} else {
						power_run = power_runs["val"] #there's only one of them!
					}
					power_run = as.numeric(power_run)
					#let's also plot the t-test power point
					if (allocation_method == "SM" & test_type != "Exact"){					
						power_run_z_0 = power_runs[power_runs$z_test == 0 & power_runs$lambda == lambda, "val"]
						moe = 1.96 * sqrt(power_run_z_0 * (1 - power_run_z_0) / NUM_SIMs)
						segments(x_pos, power_run_z_0 - moe, y1 = power_run_z_0 + moe, lwd = line_thickness, col = "gray")
						points(x_pos, power_run_z_0, ylim = c(0, 1), cex = dot_size, pch = pch_symbols_z[n_i], col = "gray")
					} 
					#now do the actual point
					if (length(power_run) == 1){
						moe = 1.96 * sqrt(power_run * (1 - power_run) / NUM_SIMs)
						segments(x_pos, power_run - moe, y1 = power_run + moe, lwd = line_thickness)
						points(x_pos, power_run, ylim = c(0, 1), cex = dot_size, pch = pch_symbols[n_i])
					}
					
				}
			}
			if (make_pdfs){dev.off()}
		}
	}
	graphics.off()
}

#this will create a series of 9 PDF's which can be subfigured to create Figure 1
draw_power_fig(lambda = 0.10)





#### 
#### 
#### 
#### Regenerate Table 2: Balance and Efficiency
#### 
#### 
#### 

#first do balances results
balance_results = list()

for (n in Ns_in_experiment){
	balance_results[[as.character(n)]] = list()
	for (scenario in SCENARIOS){
		for (test_type in TEST_TYPES){
		for (allocation_method in ALLOCATION_METHODS){
				balance = balances[balances$allocation_method == allocation_method & balances$n == n, ]
				if (allocation_method %in% c("SM", "SPM", "CPM")){ #further limit by lambda here
					balance = balance[balance$lambda == lambda, ]
				}
				balance_results[[as.character(n)]][[allocation_method]] = mean(as.numeric(balance$val))
			}
		}
	}
}

#copy these results below into LaTeX to regenerate the table
balance_results


#now do sample relative efficiencies

s_sq_sm_over_s_sq_competitor = list()

for (n in Ns_in_experiment){
	s_sq_sm_over_s_sq_competitor[[as.character(n)]] = list()
	for (scenario in SCENARIOS){
		s_sq_sm_over_s_sq_competitor[[as.character(n)]][[scenario]] = list()
		for (test_type in TEST_TYPES){
			if (test_type == "Exact"){
				next
			}
			s_sq_sm_over_s_sq_competitor[[as.character(n)]][[scenario]][[test_type]] = list()

			for (allocation_method in ALLOCATION_METHODS){
				std_err_beta_T_sm = std_err_beta_Ts[std_err_beta_Ts$model_name == scenario & 
					std_err_beta_Ts$allocation_method == "SM" & 
					std_err_beta_Ts$n == n & 
					std_err_beta_Ts$test_type == test_type & 
					std_err_beta_Ts$z_test == 1 & 
					std_err_beta_Ts$beta_T == 1 & 
					std_err_beta_Ts$lambda == lambda, ]
				if (allocation_method == "SM"){ #it's relative to this
					next
				}				
				std_err_beta_T_alloc = std_err_beta_Ts[std_err_beta_Ts$model_name == scenario & 
					std_err_beta_Ts$allocation_method == allocation_method & 
					std_err_beta_Ts$test_type == test_type & 
					std_err_beta_Ts$beta_T == 1 &
					std_err_beta_Ts$n == n, ]

				rel_efficiency = std_err_beta_T_alloc$val^2 / std_err_beta_T_sm$val^2
				pval = pf(ifelse(rel_efficiency < 1, 1 / rel_efficiency, rel_efficiency), NUM_SIMs - 1, NUM_SIMs - 1, lower.tail = F)
				
				s_sq_sm_over_s_sq_competitor[[as.character(n)]][[scenario]][[test_type]][[allocation_method]] = c(round(rel_efficiency, 3), round(pval, 3))
			}
		}
	}
}

#copy these results below into LaTeX to regenerate the table
s_sq_sm_over_s_sq_competitor[["50"]][["NL"]]
s_sq_sm_over_s_sq_competitor[["50"]][["LI"]]
s_sq_sm_over_s_sq_competitor[["50"]][["ZE"]]

s_sq_sm_over_s_sq_competitor[["100"]][["NL"]]
s_sq_sm_over_s_sq_competitor[["100"]][["LI"]]
s_sq_sm_over_s_sq_competitor[["100"]][["ZE"]]

s_sq_sm_over_s_sq_competitor[["200"]][["NL"]]
s_sq_sm_over_s_sq_competitor[["200"]][["LI"]]
s_sq_sm_over_s_sq_competitor[["200"]][["ZE"]]


#### 
#### 
#### 
#### Regenerate Table 3: Sizes
#### 
#### 
#### 

size_results = list()

for (n in Ns_in_experiment){
	size_results[[as.character(n)]] = list()
	for (scenario in SCENARIOS){
		size_results[[as.character(n)]][[scenario]] = list()
		size_scenario = sizes[sizes$model_name == scenario, ]		
		for (test_type in TEST_TYPES){
			size_results[[as.character(n)]][[scenario]][[test_type]] = list()
			size_test_type = size_scenario[size_scenario$test_type == test_type, ]
			for (allocation_method in ALLOCATION_METHODS){
				size_alloc = size_test_type[size_test_type$allocation_method == allocation_method & size_test_type$n == n, ]
				if (allocation_method == "SM"){
					sizes_lambda = size_alloc[size_alloc$lambda == lambda, ]
					#okay so now we need to get the Z
					if (test_type == "Exact"){
						pval = prop.test(sizes_lambda$val * NUM_SIMs, NUM_SIMs, 0.05)$p.value
						size_results[[as.character(n)]][[scenario]][[test_type]][[allocation_method]] = c(round(sizes_lambda$val, 3), ifelse(pval * 162 < 0.05, "**", ifelse(pval < 0.05, "*", "NONE")))
					} else {
						
						size_z = sizes_lambda[sizes_lambda$z_test == 1, ]
						pval = prop.test(size_z$val * NUM_SIMs, NUM_SIMs, 0.05)$p.value
						size_results[[as.character(n)]][[scenario]][[test_type]][[paste(allocation_method, "z", sep = "")]] = c(round(size_z$val, 3), ifelse(pval * 162 < 0.05, "**", ifelse(pval < 0.05, "*", "NONE")))

						size_t = sizes_lambda[sizes_lambda$z_test == 0, ]
						pval = prop.test(size_t$val * NUM_SIMs, NUM_SIMs, 0.05)$p.value
						size_results[[as.character(n)]][[scenario]][[test_type]][[paste(allocation_method, "t", sep = "")]] = c(round(size_t$val, 3), ifelse(pval * 162 < 0.05, "**", ifelse(pval < 0.05, "*", "NONE")))

					}
					
				} else {
					pval = prop.test(size_alloc$val * NUM_SIMs, NUM_SIMs, 0.05)$p.value
					size_results[[as.character(n)]][[scenario]][[test_type]][[allocation_method]] = c(round(size_alloc$val, 3), ifelse(pval * 162 < 0.05, "**", ifelse(pval < 0.05, "*", "NONE")))
				}			
			}
		}
	}
}

#copy these results below into LaTeX to regenerate the table
size_results[["50"]][["NL"]]
size_results[["50"]][["LI"]]
size_results[["50"]][["ZE"]]

size_results[["100"]][["NL"]]
size_results[["100"]][["LI"]]
size_results[["100"]][["ZE"]]

size_results[["200"]][["NL"]]
size_results[["200"]][["LI"]]
size_results[["200"]][["ZE"]]


#### 
#### 
#### 
#### Regenerate Supplementary Materials: Table 1 - Biases
#### 
#### 
#### 

ALLOCATION_METHODS = c("C", "E", "S", "M", "SM", "SPM")

bias_results = list()

for (n in Ns_in_experiment){
	bias_results[[as.character(n)]] = list()
	for (scenario in SCENARIOS){
		bias_results[[as.character(n)]][[scenario]] = list()
		for (test_type in TEST_TYPES){
			if (test_type == "Exact"){
				next
			}
			bias_results[[as.character(n)]][[scenario]][[test_type]] = list()
			for (allocation_method in ALLOCATION_METHODS){
				bias_alloc = biases[biases$n == n & biases$model_name == scenario & biases$test_type == test_type & biases$allocation_method == allocation_method, ]
				if (allocation_method == "SM" || allocation_method == "SPM"){ #we need to just pull out the lambda value that's appropriate
					bias_alloc = bias_alloc[bias_alloc$lambda == lambda, ]					
				}
				bias_results[[as.character(n)]][[scenario]][[test_type]][[allocation_method]] = round(bias_alloc$val, 3)
			}
		}
	}
}

#copy these results below into LaTeX to regenerate the table
bias_results[["50"]][["NL"]]
bias_results[["50"]][["LI"]]
bias_results[["50"]][["ZE"]]

bias_results[["100"]][["NL"]]
bias_results[["100"]][["LI"]]
bias_results[["100"]][["ZE"]]

bias_results[["200"]][["NL"]]
bias_results[["200"]][["LI"]]
bias_results[["200"]][["ZE"]]

#### 
#### 
#### 
#### Regenerate Supplementary Materials: Figures 1-8 - Power results by lambda
#### 
#### 
#### 

LAMBDAS = c(0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.35, 0.5)
for (l in LAMBDAS){
	draw_power_fig(lambda = l, allocation_methods = c("C", "E", "S", "M", "SPM", "SM"))
}
