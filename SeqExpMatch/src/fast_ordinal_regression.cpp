#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

class OrdinalRegression {
private:
    const MatrixXd m_X;
    const VectorXd m_y;
    const std::vector<double> m_levels;
    const int m_n;
    const int m_p;
    const int m_K;

public:
    OrdinalRegression(const MatrixXd& X, const VectorXd& y) :
        m_X(X), m_y(y), m_n(X.rows()), m_p(X.cols()),
        m_levels(init_levels(y)), m_K(m_levels.size()) {}

    static std::vector<double> init_levels(const VectorXd& y) {
        std::vector<double> levels;
        for (int i = 0; i < y.size(); ++i) {
            if (std::find(levels.begin(), levels.end(), y[i]) == levels.end()) {
                levels.push_back(y[i]);
            }
        }
        std::sort(levels.begin(), levels.end());
        return levels;
    }

    double neg_log_likelihood(const VectorXd& params) const {
        int n_alpha = m_K - 1;
        VectorXd alpha = params.head(n_alpha);
        VectorXd beta = params.tail(m_p);
        
        // Ensure alpha is sorted (constraint)
        for (int k = 1; k < n_alpha; ++k) {
            if (alpha[k] <= alpha[k-1]) return 1e10; 
        }

        VectorXd eta = m_X * beta;
        double nll = 0.0;
        
        auto logit_inv = [](double x) { 
            if (x > 20) return 1.0;
            if (x < -20) return 0.0;
            return 1.0 / (1.0 + std::exp(-x)); 
        };

        for (int i = 0; i < m_n; ++i) {
            int yi_idx = -1;
            for (int k = 0; k < m_K; ++k) {
                if (m_y[i] == m_levels[k]) {
                    yi_idx = k;
                    break;
                }
            }

            double p_upper = (yi_idx == m_K - 1) ? 1.0 : logit_inv(alpha[yi_idx] - eta[i]);
            double p_lower = (yi_idx == 0) ? 0.0 : logit_inv(alpha[yi_idx - 1] - eta[i]);
            
            double prob = std::max(1e-12, p_upper - p_lower);
            nll -= std::log(prob);
        }
        return nll;
    }

    MatrixXd hessian(const VectorXd& params) const {
        int n_params = params.size();
        MatrixXd H(n_params, n_params);
        double h = 1e-4;

        for (int i = 0; i < n_params; ++i) {
            for (int j = i; j < n_params; ++j) {
                VectorXd p_pp = params; p_pp[i] += h; p_pp[j] += h;
                VectorXd p_pm = params; p_pm[i] += h; p_pm[j] -= h;
                VectorXd p_mp = params; p_mp[i] -= h; p_mp[j] += h;
                VectorXd p_mm = params; p_mm[i] -= h; p_mm[j] -= h;

                double f_pp = neg_log_likelihood(p_pp);
                double f_pm = neg_log_likelihood(p_pm);
                double f_mp = neg_log_likelihood(p_mp);
                double f_mm = neg_log_likelihood(p_mm);

                H(i, j) = (f_pp - f_pm - f_mp + f_mm) / (4.0 * h * h);
                H(j, i) = H(i, j);
            }
        }
        return H;
    }
};

// Simple solver using Newton-Raphson as we have a small number of parameters (usually)
// [[Rcpp::export]]
List fast_ordinal_regression_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int maxit = 100, double tol = 1e-6) {
    OrdinalRegression model(X, y);
    int p = X.cols();
    int K = OrdinalRegression::init_levels(y).size();
    int n_alpha = K - 1;
    int n_params = n_alpha + p;

    VectorXd params(n_params);
    // Initialize alpha
    for (int k = 0; k < n_alpha; ++k) {
        params[k] = -1.0 + 2.0 * (k + 1) / K;
    }
    // Initialize beta
    params.tail(p).setZero();

    double h = 1e-4;
    for (int iter = 0; iter < maxit; ++iter) {
        VectorXd grad(n_params);
        for (int i = 0; i < n_params; ++i) {
            VectorXd p_plus = params; p_plus[i] += h;
            VectorXd p_minus = params; p_minus[i] -= h;
            grad[i] = (model.neg_log_likelihood(p_plus) - model.neg_log_likelihood(p_minus)) / (2.0 * h);
        }

        if (grad.norm() < tol) break;

        MatrixXd H = model.hessian(params);
        FullPivLU<MatrixXd> lu(H);
        if (!lu.isInvertible()) break;

        VectorXd step = lu.solve(grad);
        
        // Simple line search
        double alpha_step = 1.0;
        double current_nll = model.neg_log_likelihood(params);
        while (alpha_step > 1e-4) {
            VectorXd next_params = params - alpha_step * step;
            if (model.neg_log_likelihood(next_params) < current_nll) {
                params = next_params;
                break;
            }
            alpha_step *= 0.5;
        }
        if (alpha_step <= 1e-4) break;
    }

    return List::create(
        Named("b") = params.tail(p),
        Named("alpha") = params.head(n_alpha),
        Named("n_params") = n_params,
        Named("params") = params
    );
}

// [[Rcpp::export]]
List fast_ordinal_regression_with_var_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) {
    List res = fast_ordinal_regression_cpp(X, y);
    VectorXd params = res["params"];
    OrdinalRegression model(X, y);
    MatrixXd H = model.hessian(params);
    
    FullPivLU<MatrixXd> lu(H);
    if (!lu.isInvertible()) {
        return List::create(Named("b") = res["b"], Named("ssq_b_2") = NA_REAL);
    }
    
    MatrixXd cov_mat = lu.inverse();
    int p = X.cols();
    int n_alpha = params.size() - p;
    
    // We want the variance of the first covariate after alphas (which is often the treatment)
    // In our params, it's at index n_alpha
    double ssq_b_2 = (p >= 1) ? cov_mat(n_alpha, n_alpha) : NA_REAL;

    return List::create(
        Named("b") = res["b"],
        Named("ssq_b_2") = ssq_b_2
    );
}

// [[Rcpp::export]]
List expand_continuation_ratio_data_cpp(const Eigen::VectorXi& y, const Eigen::VectorXi& w, const Eigen::VectorXi& strata, int K) {
    int n = y.size();
    int n_alpha = K - 1;
    int num_strata = strata.maxCoeff();
    
    std::vector<int> y_stack;
    std::vector<int> w_stack;
    std::vector<int> strata_stack;
    
    for (int i = 0; i < n; ++i) {
        int yi = y[i];
        int wi = w[i];
        int si = strata[i];
        
        for (int j = 1; j <= std::min(yi, n_alpha); ++j) {
            y_stack.push_back((yi == j) ? 1 : 0);
            w_stack.push_back(wi);
            strata_stack.push_back(si + (j - 1) * num_strata);
        }
    }
    
    return List::create(
        Named("y") = wrap(y_stack),
        Named("w") = wrap(w_stack),
        Named("strata") = wrap(strata_stack)
    );
}

// [[Rcpp::export]]
List expand_adjacent_category_data_cpp(const Eigen::VectorXi& y, const Eigen::VectorXi& w, const Eigen::VectorXi& strata, int K) {
    int n = y.size();
    int n_alpha = K - 1;
    int num_strata = strata.maxCoeff();
    
    std::vector<int> y_stack;
    std::vector<int> w_stack;
    std::vector<int> strata_stack;
    
    for (int i = 0; i < n; ++i) {
        int yi = y[i];
        int wi = w[i];
        int si = strata[i];
        
        for (int j = 1; j <= n_alpha; ++j) {
            if (yi == j || yi == j + 1) {
                y_stack.push_back((yi == j + 1) ? 1 : 0);
                w_stack.push_back(wi);
                strata_stack.push_back(si + (j - 1) * num_strata);
            }
        }
    }
    
    return List::create(
        Named("y") = wrap(y_stack),
        Named("w") = wrap(w_stack),
        Named("strata") = wrap(strata_stack)
    );
}
