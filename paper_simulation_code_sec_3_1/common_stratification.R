
#initialize the indicator treatment vector
indic_T = array(NA, n)
#3 categories for each variable... WARNING: not robust to changes in distribution of X!!!!
cutoff_a = qnorm(1 / 3)
cutoff_b = qnorm(2 / 3)
blocks = array(NA, n)


#now we're going to go through and do the matching
for (i_match in 1 : n){
	x_i = x_s[i_match, ]
	#first go through the decision tree to assign a block
	if (x_i[1] < cutoff_a){
		if (x_i[2] < cutoff_a){
			blocks[i_match] = "1-1"			
		} else if (x_i[2] < cutoff_b){
			blocks[i_match] = "1-2"
		} else {
			blocks[i_match] = "1-3"
		}		
	} else if (x_i[1] < cutoff_b){
		if (x_i[2] < cutoff_a){
			blocks[i_match] = "2-1"
		} else if (x_i[2] < cutoff_b){
			blocks[i_match] = "2-2"
		} else {
			blocks[i_match] = "2-3"
		}		
	} else {
		if (x_i[2] < cutoff_a){
			blocks[i_match] = "3-1"
		} else if (x_i[2] < cutoff_b){
			blocks[i_match] = "3-2"
		} else {
			blocks[i_match] = "3-3"
		}		
	}
#	cat("i", i_match, "\n")
	#now assign T or C
	block = blocks[i_match]
	previous_assigments = indic_T[1 : i_match][blocks == block]
	previous_assigments = previous_assigments[!is.na(previous_assigments)]
#	cat(previous_assigments, "\n")
	if (length(previous_assigments) == 0){ 
		#if there's no assignments, do it randomly
		indic_T[i_match] = rbinom(1, 1, prob_trt)
	} else {
		#otherwise, flip the previous assignment --- very generous
		indic_T[i_match] = 1 - previous_assigments[length(previous_assigments)]
	}
}

#create response vector
sys.source(paste("create_response_", response_model, ".R", sep = ""), envir = environment())

#create design matrix
Xy = data.frame(x_s, indic_T, blocks, y)
colnames(Xy) = c(paste0("x", 1 : p), "indic_T", "blocks", "y")
Xy$blocks = as.factor(Xy$blocks)

#pull out yT, yC
yTs = Xy[Xy$indic_T == 1, "y"]
yCs = Xy[Xy$indic_T == 0, "y"]