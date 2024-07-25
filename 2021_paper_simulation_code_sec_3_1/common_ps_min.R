
#initialize the indicator treatment vector
indic_T = array(NA, n)
#3 categories for each variable...
cutoff_a = qnorm(1 / 3)
cutoff_b = qnorm(2 / 3)
blocks_x1 = array(NA, n)
blocks_x2 = array(NA, n)

sum_variance_scores_by_trt = function(trt){
	indic_T[i_match] = trt	
	score_tot = 0	
	
	for (block in unique(blocks_x1)){
		trts_in_level = indic_T[which(blocks_x1 == block)]
		num_T = length(trts_in_level[trts_in_level == 1])
		num_C = length(trts_in_level[trts_in_level == 0])
		#add the range (d_ik) to the total function (G)
		score_tot = score_tot + var(c(num_T, num_C)) 
	}
	for (block in unique(blocks_x2)){
		trts_in_level = indic_T[which(blocks_x2 == block)]
		num_T = length(trts_in_level[trts_in_level == 1])
		num_C = length(trts_in_level[trts_in_level == 0])
		#add the range (d_ik) to the total function (G)
		score_tot = score_tot + var(c(num_T, num_C)) 
	}
	score_tot
}


#now we're going to go through and do the matching
for (i_match in 1 : n){
	x_i = x_s[i_match, ]
	#first go through the decision tree to assign a block
	if (x_i[1] < cutoff_a){
		blocks_x1[i_match] = 1	
	} else if (x_i[1] < cutoff_b){
		blocks_x1[i_match] = 2		
	} else {
		blocks_x1[i_match] = 3
	}
	if (x_i[2] < cutoff_a){
		blocks_x2[i_match] = 1		
	} else if (x_i[2] < cutoff_b){
		blocks_x2[i_match] = 2
	} else {
		blocks_x2[i_match] = 3
	}
#	cat("i", i_match, "\n")
	#now assign T or C
	
	#we first need to go through all blocks and get the distance functions for both a hypothetical T and C
#	table(blocks_x1, indic_T)
#	table(blocks_x2, indic_T)
	
	sum_variances_control = sum_variance_scores_by_trt(0)
	sum_variances_trt = sum_variance_scores_by_trt(1)
	
	#now assign to the treatment that has minimum sum of total distances (if they're equal, do it randomly)
	if (sum_variances_control < sum_variances_trt){
		indic_T[i_match] = 0
	} else if (sum_variances_trt < sum_variances_control){
		indic_T[i_match] = 1
	} else {
		indic_T[i_match] = rbinom(1, 1, prob_trt)
	}
}

#create response vector
sys.source(paste("create_response_", response_model, ".R", sep = ""), envir = environment())

#create design matrix
Xy = data.frame(cbind(x_s, indic_T, y))
colnames(Xy) = c(paste0("x", 1 : p), "indic_T", "y")

#pull out yT, yC
yTs = Xy[Xy$indic_T == 1, "y"]
yCs = Xy[Xy$indic_T == 0, "y"]