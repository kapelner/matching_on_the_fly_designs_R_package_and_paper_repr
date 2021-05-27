Xy$z = z #we want to keep information about what the true conditional means of the functions are

xTs = Xy[Xy$indic_T == 1, 1 : p]
xCs = Xy[Xy$indic_T == 0, 1 : p]
nC = nrow(xCs)
nT = nrow(xTs)
rownames(Xy)[Xy$indic_T == 0] = paste("c", 1 : nC, sep = "")
rownames(Xy)[Xy$indic_T == 1] = paste("t", 1 : nT, sep = "")

#given an experiment where all subjects have been allocated already, we want to match after the fact

#we want to tabulate all possible differences
mahalanobis_distances_sqd = matrix(Inf, nrow = nC, ncol = nT)
rownames(mahalanobis_distances_sqd) = paste("c", 1 : nC, sep = "")
colnames(mahalanobis_distances_sqd) = paste("t", 1 : nT, sep = "")

S_xs_inv_n = solve(var(rbind(xCs, xTs)))


F_crit =  qf(prob_match_cutoff_lambda, p, n - p)
T_cutoff_sq = p * (n - 1) / (n - p) * F_crit

for (i_0 in 1 : nC){
	for (j_0 in 1 : nT){
		x_diff = t(as.matrix(xCs[i_0, ] - xTs[j_0, ]))
		mahalanobis_distance_sqd = 1 / 2 * t(x_diff) %*% S_xs_inv_n %*% x_diff
		
		#now ensure some cannot be matched because they're too far apart in X space
		mahalanobis_distances_sqd[i_0, j_0] = ifelse(mahalanobis_distance_sqd > T_cutoff_sq, Inf, mahalanobis_distance_sqd)
	}
}

#round(mahalanobis_distances_sqd, 2)


#make (or reset) the match vector
Xy$match_indic = 0

#load library
tryCatch(library(optmatch, quietly = TRUE), error = function(e){install.packages("optmatch")}, finally = library(optmatch, quietly = TRUE))

#####optmatch sometimes does not produce matches
match_obj = tryCatch(
	pairmatch(mahalanobis_distances_sqd, data = Xy, remove.unmatchables = TRUE), 
	#if there's an error, transpose it
	error = function(e){
		return(tryCatch(
				pairmatch(t(mahalanobis_distances_sqd), data = Xy, remove.unmatchables = TRUE), 
				#if there's an error, transpose it
				error = function(e){NULL})
		)
	}
)

#sometimes optmatch just blows up, so just catch that here and skip to the next step
if (!is.null(match_obj)){
	control_matches = match_obj[paste("c", 1 : nC, sep = "")]
	treatment_matches = match_obj[paste("t", 1 : nT, sep = "")]
	
	
	#summary(match_obj)
	#matched(match_obj)
	#unmatched(match_obj)
	#matchfailed(match_obj)
	#sort(control_matches)
	#sort(treatment_matches)
	
	
	
	match_num = 0
	for (i_0 in 1 : nC){
		control_index = names(control_matches)[i_0]
		match = as.character(control_matches[i_0])
		if (is.na(match)){
			next
		}
		match_num = match_num + 1
		
		treatment_index = names(treatment_matches)[which(treatment_matches == match)]
	#	cat("c", control_index, "t", treatment_index, "match_num", match_num, "\n")
		#now mark it in Xy
		Xy[control_index, "match_indic"] = match_num
		Xy[treatment_index, "match_indic"] = match_num
	}
}
Xy[order(Xy$match_indic), ]




####GREEDY MATCHING ALWAYS WORKS BUT CAN BE INFINITELY BAD
#match_num = 0
#while (TRUE){
#	#if there's either no controls or treatments left, jet or if there's no more viable matches, jet
#	if (0 %in% dim(mahalanobis_distances_sqd) || min(mahalanobis_distances_sqd) == Inf){
#		break
#	}
#	match_indices = which(mahalanobis_distances_sqd == min(mahalanobis_distances_sqd), arr.ind = TRUE)
#		
#	matched_control = rownames(mahalanobis_distances_sqd)[match_indices[1]]
#	matched_treatment = colnames(mahalanobis_distances_sqd)[match_indices[2]]
#	
#	#now make a new match
#	match_num = match_num + 1
#	Xy[matched_control, "match_indic"] = match_num
#	Xy[matched_treatment, "match_indic"] = match_num
##	cat("made match", match_num, "c =", matched_control, "t =", matched_treatment, "\n")
#	
#	#remove match from dist matrix
#	mahalanobis_distances_sqd = mahalanobis_distances_sqd[-c(match_indices[1]), , drop = FALSE]
#	mahalanobis_distances_sqd = mahalanobis_distances_sqd[, -c(match_indices[2]), drop = FALSE]
##	print(mahalanobis_distances_sqd)
#}
#Xy[order(Xy$match_indic), ]


#now before we get to post processing, we should kill the block information. We only stratified to get better matches,
#we no longer need this information
Xy$blocks = NULL

#now we can do the post-processing the same as if it were sequentially matched
sys.source("common_matching_post_processing.R", envir = environment())
