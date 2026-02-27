#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_kk_compound_distr_parallel_cpp(
    const Eigen::VectorXd& y,
    const Eigen::MatrixXi& w_mat,
    const Eigen::MatrixXi& match_indic_mat,
    int num_cores) {
    
    int nsim = w_mat.cols();
    int n = y.size();
    NumericVector results(nsim);
    
#ifdef _OPENMP
    if (num_cores > 1) {
        omp_set_num_threads(num_cores);
    }
#endif

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < nsim; ++b) {
        Eigen::VectorXi w = w_mat.col(b);
        Eigen::VectorXi match_indic = match_indic_mat.col(b);
        
        int m = match_indic.maxCoeff();
        
        double d_bar = NA_REAL;
        double ssqD_bar = NA_REAL;
        
        if (m > 0) {
            Eigen::VectorXd diffs(m);
            for (int match_id = 1; match_id <= m; match_id++) {
                double yT = NA_REAL, yC = NA_REAL;
                bool foundT = false, foundC = false;
                for (int i = 0; i < n; i++) {
                    if (match_indic[i] == match_id) {
                        if (w[i] == 1) {
                            yT = y[i];
                            foundT = true;
                        } else {
                            yC = y[i];
                            foundC = true;
                        }
                        if (foundT && foundC) break;
                    }
                }
                diffs[match_id - 1] = yT - yC;
            }
            
            d_bar = diffs.mean();
            
            if (m > 1) {
                ssqD_bar = (diffs.array() - d_bar).square().sum() / (m - 1) / m;
            }
        }
        
        int nRT = 0, nRC = 0;
        double sum_T = 0, sum_C = 0;
        
        for (int i = 0; i < n; i++) {
            if (match_indic[i] == 0) {
                if (w[i] == 1) {
                    nRT++;
                    sum_T += y[i];
                } else {
                    nRC++;
                    sum_C += y[i];
                }
            }
        }
        
        double r_bar = NA_REAL;
        double ssqR = NA_REAL;
        
        if (nRT > 0 && nRC > 0) {
            double mean_T = sum_T / nRT;
            double mean_C = sum_C / nRC;
            r_bar = mean_T - mean_C;
            
            if (nRT > 1 && nRC > 1 && (nRT + nRC) > 2) {
                double sq_diff_T = 0, sq_diff_C = 0;
                for (int i = 0; i < n; i++) {
                    if (match_indic[i] == 0) {
                        if (w[i] == 1) {
                            sq_diff_T += (y[i] - mean_T) * (y[i] - mean_T);
                        } else {
                            sq_diff_C += (y[i] - mean_C) * (y[i] - mean_C);
                        }
                    }
                }
                
                double var_T = sq_diff_T / (nRT - 1);
                double var_C = sq_diff_C / (nRC - 1);
                int nR = nRT + nRC;
                
                ssqR = (var_T * (nRT - 1) + var_C * (nRC - 1)) / (nR - 2) * (1.0 / nRT + 1.0 / nRC);
            }
        }
        
        double beta_hat_T = NA_REAL;
        
        if (nRT <= 1 || nRC <= 1) {
            beta_hat_T = d_bar;
        } else if (m == 0) {
            beta_hat_T = r_bar;
        } else {
            if (!std::isfinite(ssqD_bar) || ssqD_bar <= 0) {
                beta_hat_T = r_bar; 
            } else if (!std::isfinite(ssqR) || ssqR <= 0) {
                beta_hat_T = d_bar; 
            } else {
                double w_star = ssqR / (ssqR + ssqD_bar);
                beta_hat_T = w_star * d_bar + (1.0 - w_star) * r_bar;
            }
        }
        
        results[b] = beta_hat_T;
    }
    
    return results;
}

// [[Rcpp::export]]
NumericVector compute_kk_compound_bootstrap_parallel_cpp(
    const Eigen::MatrixXd& y_mat,
    const Eigen::MatrixXi& w_mat,
    const Eigen::MatrixXi& match_indic_mat,
    int num_cores) {
    
    int nsim = w_mat.cols();
    int n = w_mat.rows();
    NumericVector results(nsim);
    
#ifdef _OPENMP
    if (num_cores > 1) {
        omp_set_num_threads(num_cores);
    }
#endif

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < nsim; ++b) {
        Eigen::VectorXd y = y_mat.col(b);
        Eigen::VectorXi w = w_mat.col(b);
        Eigen::VectorXi match_indic = match_indic_mat.col(b);
        
        int m = match_indic.maxCoeff();
        
        double d_bar = NA_REAL;
        double ssqD_bar = NA_REAL;
        
        if (m > 0) {
            Eigen::VectorXd diffs(m);
            for (int match_id = 1; match_id <= m; match_id++) {
                double yT = NA_REAL, yC = NA_REAL;
                bool foundT = false, foundC = false;
                for (int i = 0; i < n; i++) {
                    if (match_indic[i] == match_id) {
                        if (w[i] == 1) {
                            yT = y[i];
                            foundT = true;
                        } else {
                            yC = y[i];
                            foundC = true;
                        }
                        if (foundT && foundC) break;
                    }
                }
                diffs[match_id - 1] = yT - yC;
            }
            d_bar = diffs.mean();
            if (m > 1) {
                ssqD_bar = (diffs.array() - d_bar).square().sum() / (m - 1) / m;
            }
        }
        
        int nRT = 0, nRC = 0;
        double sum_T = 0, sum_C = 0;
        
        for (int i = 0; i < n; i++) {
            if (match_indic[i] == 0) {
                if (w[i] == 1) {
                    nRT++;
                    sum_T += y[i];
                } else {
                    nRC++;
                    sum_C += y[i];
                }
            }
        }
        
        double r_bar = NA_REAL;
        double ssqR = NA_REAL;
        
        if (nRT > 0 && nRC > 0) {
            double mean_T = sum_T / nRT;
            double mean_C = sum_C / nRC;
            r_bar = mean_T - mean_C;
            
            if (nRT > 1 && nRC > 1 && (nRT + nRC) > 2) {
                double sq_diff_T = 0, sq_diff_C = 0;
                for (int i = 0; i < n; i++) {
                    if (match_indic[i] == 0) {
                        if (w[i] == 1) {
                            sq_diff_T += (y[i] - mean_T) * (y[i] - mean_T);
                        } else {
                            sq_diff_C += (y[i] - mean_C) * (y[i] - mean_C);
                        }
                    }
                }
                double var_T = sq_diff_T / (nRT - 1);
                double var_C = sq_diff_C / (nRC - 1);
                int nR = nRT + nRC;
                ssqR = (var_T * (nRT - 1) + var_C * (nRC - 1)) / (nR - 2) * (1.0 / nRT + 1.0 / nRC);
            }
        }
        
        double beta_hat_T = NA_REAL;
        if (nRT <= 1 || nRC <= 1) {
            beta_hat_T = d_bar;
        } else if (m == 0) {
            beta_hat_T = r_bar;
        } else {
            if (!std::isfinite(ssqD_bar) || ssqD_bar <= 0) {
                beta_hat_T = r_bar;
            } else if (!std::isfinite(ssqR) || ssqR <= 0) {
                beta_hat_T = d_bar;
            } else {
                double w_star = ssqR / (ssqR + ssqD_bar);
                beta_hat_T = w_star * d_bar + (1.0 - w_star) * r_bar;
            }
        }
        results[b] = beta_hat_T;
    }
    
    return results;
}
