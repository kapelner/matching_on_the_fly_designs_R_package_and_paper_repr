
#create design matrix
Xy = as.data.frame(cbind(x_s, indic_T, z, y, match_indic))
colnames(Xy) = c(paste0("x", 1 : p), "indic_T", "z", "y", "match_indic")

#pull out yT, yC, xT, xC
yTs = Xy[Xy$indic_T == 1, "y"]
yCs = Xy[Xy$indic_T == 0, "y"]
xTs = Xy[Xy$indic_T == 1, 1 : p]
xCs = Xy[Xy$indic_T == 0, 1 : p]

#and record other useful info
m = max(Xy$match_indic)

#now we want to just calculate the diffs inside matches and ignore the reservoir
ydiffs = array(NA, m)
zdiffs = array(NA, m)

if (m > 0){ #R runs the for loop if it's over the collection 1 : 0 => how stupid is that?
	for (match_id in 1 : m){
		yT = Xy[Xy$indic_T == 1 & Xy$match_indic == match_id, "y"]
		yC = Xy[Xy$indic_T == 0 & Xy$match_indic == match_id, "y"]
		ydiffs[match_id] = yT - yC
		zT = Xy[Xy$indic_T == 1 & Xy$match_indic == match_id, "z"]
		zC = Xy[Xy$indic_T == 0 & Xy$match_indic == match_id, "z"]
		zdiffs[match_id] = zT - zC	
	}
}

#the estimate is the one sample t test from Stat 101
if (length(ydiffs) <= 1){
	d_bar = 0
	beta_hat_T = d_bar
	pval = 1
	Rsq = 0
} else {
	d_bar = mean(ydiffs)
	ttest = t.test(ydiffs)
	beta_hat_T = d_bar
	pval = ttest$p.value
	Rsq = 1 - sum((ydiffs - beta_hat_T)^2) / sum((ydiffs - mean(ydiffs))^2)
}


#get reservoir data
Xyleft = Xy[Xy$match_indic == 0, ] #pull out the reservoir data
final_reservoir_size = nrow(Xyleft) / n

YleftT = Xyleft[Xyleft$indic_T == 1, ]$y #get the reservoir responses from the treatment
YleftC = Xyleft[Xyleft$indic_T == 0, ]$y #get the reservoir responses from the control
r_bar = mean(YleftT) - mean(YleftC) #compute the classic estimator from the reservoir: ybar_T - ybar_C

#get reservoir sample sizes
nR = nrow(Xyleft) #how many observations are there in the reservoir?
nRT = length(YleftT) #how many treatment observations are there in the reservoir?
nRC = length(YleftC) #how many control observations are there in the reservoir?

ssqR = (var(YleftT) * (nRT - 1) + var(YleftC) * (nRC - 1)) / (nR - 2) * (1 / nRT + 1 / nRC)
ssqD_bar = var(ydiffs) / m

w_star = ssqR / (ssqR + ssqD_bar)

Xymatched = Xy[Xy$match_indic > 0, ]
Xymatched = Xymatched[order(Xymatched$match_indic), ]

Xyleft = Xyleft[order(Xyleft$indic_T), ]