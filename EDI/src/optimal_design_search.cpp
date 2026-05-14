#include <RcppEigen.h>
#include <algorithm>
#include <vector>
#include <random>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix d_optimal_search_cpp(const Eigen::MatrixXd& P, int nsim, int n_T) {
    const int n = P.rows();
    const double* p_ptr = P.data();
    const Eigen::VectorXd p_diag = P.diagonal();
    IntegerMatrix w_mat(n, nsim);
    
    std::vector<int> indices(n);
    for(int i=0; i<n; ++i) indices[i] = i;
    
    std::random_device rd;
    std::mt19937 g(rd());

    for (int s = 0; s < nsim; ++s) {
        std::shuffle(indices.begin(), indices.end(), g);
        Eigen::VectorXd w(n);
        w.setZero();
        std::vector<int> t_idxs;
        std::vector<int> c_idxs;
        for (int i = 0; i < n; ++i) {
            if (i < n_T) {
                w(indices[i]) = 1.0;
                t_idxs.push_back(indices[i]);
            } else {
                c_idxs.push_back(indices[i]);
            }
        }

        Eigen::VectorXd Pw = P * w;
        bool improved = true;
        while (improved) {
            improved = false;
            double best_delta = -1e-10;
            int best_i = -1;
            int best_j = -1;
            int best_t_pos = -1;
            int best_c_pos = -1;

            const double* pw_ptr = Pw.data();
            for (int ti = 0; ti < (int)t_idxs.size(); ++ti) {
                const int i = t_idxs[ti];
                for (int cj = 0; cj < (int)c_idxs.size(); ++cj) {
                    const int j = c_idxs[cj];
                    const double delta = -2.0 * pw_ptr[i] + 2.0 * pw_ptr[j] +
                        p_diag[i] + p_diag[j] - 2.0 * p_ptr[static_cast<size_t>(j) * n + i];
                    
                    if (delta < best_delta) {
                        best_delta = delta;
                        best_i = i;
                        best_j = j;
                        best_t_pos = ti;
                        best_c_pos = cj;
                        improved = true;
                    }
                }
            }

            if (improved) {
                Pw -= P.col(best_i);
                Pw += P.col(best_j);
                
                // Swap indices in t_idxs and c_idxs
                t_idxs[best_t_pos] = best_j;
                c_idxs[best_c_pos] = best_i;
                w(best_i) = 0.0;
                w(best_j) = 1.0;
            }
        }

        for (int i = 0; i < n; ++i) {
            w_mat(i, s) = (int)w(i);
        }
    }

    return w_mat;
}

// [[Rcpp::export]]
IntegerMatrix a_optimal_search_cpp(const Eigen::MatrixXd& P, const Eigen::MatrixXd& H, int nsim, int n_T) {
    const int n = P.rows();
    const double* p_ptr = P.data();
    const double* h_ptr = H.data();
    const Eigen::VectorXd p_diag = P.diagonal();
    const Eigen::VectorXd h_diag = H.diagonal();
    IntegerMatrix w_mat(n, nsim);
    
    std::vector<int> indices(n);
    for(int i=0; i<n; ++i) indices[i] = i;
    
    std::random_device rd;
    std::mt19937 g(rd());

    for (int s = 0; s < nsim; ++s) {
        std::shuffle(indices.begin(), indices.end(), g);
        Eigen::VectorXd w(n);
        w.setZero();
        std::vector<int> t_idxs;
        std::vector<int> c_idxs;
        for (int i = 0; i < n; ++i) {
            if (i < n_T) {
                w(indices[i]) = 1.0;
                t_idxs.push_back(indices[i]);
            } else {
                c_idxs.push_back(indices[i]);
            }
        }

        Eigen::VectorXd Pw = P * w;
        Eigen::VectorXd Hw = H * w;
        double wPw = w.dot(Pw);
        double wHw = w.dot(Hw);
        
        // Objective: (wHw + 1) / (n_T - wPw)
        double obj_curr = (wHw + 1.0) / (n_T - wPw);

        bool improved = true;
        while (improved) {
            improved = false;
            double best_obj = obj_curr - 1e-12;
            int best_i = -1;
            int best_j = -1;
            int best_t_pos = -1;
            int best_c_pos = -1;
            double best_wPw = wPw;
            double best_wHw = wHw;
            const double* pw_ptr = Pw.data();
            const double* hw_ptr = Hw.data();

            for (int ti = 0; ti < (int)t_idxs.size(); ++ti) {
                const int i = t_idxs[ti];
                for (int cj = 0; cj < (int)c_idxs.size(); ++cj) {
                    const int j = c_idxs[cj];
                    
                    const double delta_wPw = -2.0 * pw_ptr[i] + 2.0 * pw_ptr[j] +
                        p_diag[i] + p_diag[j] - 2.0 * p_ptr[static_cast<size_t>(j) * n + i];
                    const double delta_wHw = -2.0 * hw_ptr[i] + 2.0 * hw_ptr[j] +
                        h_diag[i] + h_diag[j] - 2.0 * h_ptr[static_cast<size_t>(j) * n + i];
                    
                    double next_wPw = wPw + delta_wPw;
                    double next_wHw = wHw + delta_wHw;
                    
                    double denom = n_T - next_wPw;
                    if (denom <= 1e-10) continue; // Should not happen for valid designs
                    
                    double next_obj = (next_wHw + 1.0) / denom;
                    
                    if (next_obj < best_obj) {
                        best_obj = next_obj;
                        best_i = i;
                        best_j = j;
                        best_t_pos = ti;
                        best_c_pos = cj;
                        best_wPw = next_wPw;
                        best_wHw = next_wHw;
                        improved = true;
                    }
                }
            }

            if (improved) {
                obj_curr = best_obj;
                Pw -= P.col(best_i);
                Pw += P.col(best_j);
                Hw -= H.col(best_i);
                Hw += H.col(best_j);
                wPw = best_wPw;
                wHw = best_wHw;
                
                t_idxs[best_t_pos] = best_j;
                c_idxs[best_c_pos] = best_i;
                w(best_i) = 0.0;
                w(best_j) = 1.0;
            }
        }

        for (int i = 0; i < n; ++i) {
            w_mat(i, s) = (int)w(i);
        }
    }

    return w_mat;
}
