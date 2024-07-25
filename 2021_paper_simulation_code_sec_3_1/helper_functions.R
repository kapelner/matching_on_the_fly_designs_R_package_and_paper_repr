stepwise_weights = function(X_stepwise, y_stepwise){	
  p = ncol(X_stepwise)
  weights = array(NA, p)
  
  j_droppeds = c()
  #initialize cols to be all cols
  cols = setdiff(1 : p, j_droppeds)
	
  #iteratively add all variables and extract weights for matching later
  repeat {
  	#now get all the Rsqs and record the highest one i.e. the "best" covariate
  	rsqs = find_rsqs(X_stepwise, y_stepwise, cols)
  	j_max = which.max(rsqs)
  	weights[j_max] = rsqs[j_max]
  	j_droppeds = c(j_droppeds, j_max)
  	
  	#register the remaining cols
  	cols = setdiff(1 : p, j_droppeds)
  	#if there's none left, we jet
  	if (length(cols) == 0){
  		break
  	}
  	
  	#we now need to adjust the other covariates for the "best" covariate
  	x_best = X_stepwise[, j_max]
	X_st_w = cbind(1, x_best)
	X_st_w_tr = t(X_st_w)
	XtXinvXt_w_tr = X_st_w %*% solve(X_st_w_tr %*% X_st_w) %*% X_st_w_tr
  	for (j in cols){
#  		#predict x_j using x_best
#  		mod = lm(X_stepwise[, j] ~ x_best)
#  		#then set x_j equal to the residuals
#  		X_stepwise[, j] = summary(mod)$residuals
			
		#predict x_j using x_best
		yhat_st_w_j = XtXinvXt_w_tr %*% X_stepwise[, j]
		#then set x_j equal to the residuals
		X_stepwise[, j] = X_stepwise[, j] - yhat_st_w_j
  	}
  }
  
  #return normalized weights
  weights / sum(weights)
}

find_rsqs = function(X_stepwise, y_stepwise, cols){
	p = ncol(X_stepwise)
	Rsqs = array(NA, p)
	for (j in cols){
		Rsqs[j] = cor(y_stepwise, X_stepwise[, j])^2 #summary(lm(y_stepwise ~ X_stepwise[, j]))$r.squared
	}
	Rsqs
}

is_defined = function(sym) {
  sym = deparse(substitute(sym))
  env = parent.frame()
  exists(sym, env)
}

is_defined_or_NA = function(sym) {
  sym = deparse(substitute(sym))
  env = parent.frame()
  env = environment()
  
  ifelse(exists(sym, env), eval(parse(text = sym)), NA)
}

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}