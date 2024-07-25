

atkinson_assignment = function(X){
  X_with_T_and_int = cbind(rbinom(n, 1, prob_trt), 1, X)
  prob_Ts = array(NA, n)
  check_probTs = array(NA, n)
  for (t in (p + 2 + 1) : n){
    X_t_min_1 = X_with_T_and_int[1 : (t - 1), ]	
    
    M = (t - 1) * solve(t(X_t_min_1) %*% X_t_min_1)
    
    A = M[1, 2 : (p + 2)] %*% c(1, X[t, ]) 
    s_over_A_plus_one_sq = (M[1, 1] / A + 1)^2
    prob_T = s_over_A_plus_one_sq / (s_over_A_plus_one_sq + 1)
    prob_Ts[t] = prob_T
    
    n1 = sum(X_with_T_and_int[1 : (t-1), 1])
    n2 = (t - 1) - n1
    check_probTs[t] = n2^2 / (n1^2 + n2^2)
    X_with_T_and_int[t, 1] = rbinom(1, 1, prob_T)
  }	
  # prob_Ts
  # check_probTs
  # X_with_T_and_int[, 1] #indicT
  
  list(indicT = X_with_T_and_int[, 1], prob_Ts = prob_Ts, check_probTs = check_probTs)
}


n = 100
p = 10
prob_trt = 0.5
X = matrix(runif(n * p), nrow = n, ncol = p)
atkinson_assignment(X)
