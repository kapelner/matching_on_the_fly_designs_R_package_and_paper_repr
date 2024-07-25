
efron_biased_coin_design = function(n, prob_trt, bias = 2 / 3){ #his personal favorite
	#initialize the indicator treatment vector
	indic_T = array(NA, n)
	
	for (i in 1 : n){
		if (sum(indic_T == 1, na.rm = TRUE) == sum(indic_T == 0, na.rm = TRUE)){
			indic_T[i] = rbinom(1, 1, prob_trt)
		} else if (sum(indic_T == 1, na.rm = TRUE) < sum(indic_T == 0, na.rm = TRUE)){
			indic_T[i] = rbinom(1, 1, bias)
		} else if (sum(indic_T == 1, na.rm = TRUE) > sum(indic_T == 0, na.rm = TRUE)){
			indic_T[i] = rbinom(1, 1, 1 - bias)
		}
	}	
	indic_T
}

indic_T = efron_biased_coin_design(n, prob_trt)

#create response vector
sys.source(paste("create_response_", response_model, ".R", sep = ""), envir = environment())

#create design matrix
Xy = as.data.frame(cbind(x_s, indic_T, y))
colnames(Xy) = c(paste0("x", 1 : p), "indic_T", "y")

#pull out yT, yC
yTs = Xy[Xy$indic_T == 1, "y"]
yCs = Xy[Xy$indic_T == 0, "y"]