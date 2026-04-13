#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

class OrdinalCLLRegression {
private:
    const MatrixXd m_X;
    const VectorXd m_y;
    const std::vector<double> m_levels;
    const int m_n;
    const int m_p;
    const int m_K;

public:
    OrdinalCLLRegression(const MatrixXd& X, const VectorXd& y) :
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
        
        auto cloglog_inv = [](double x) { 
            if (x > 5) return 1.0;
            if (x < -20) return 0.0;
            return 1.0 - std::exp(-std::exp(x)); 
        };

        for (int i = 0; i < m_n; ++i) {
            int yi_idx = -1;
            for (int k = 0; k < m_K; ++k) {
                if (m_y[i] == m_levels[k]) {
                    yi_idx = k;
                    break;
                }
            }

            // cloglog(P(Y <= k)) = alpha_k + eta_i
            // P(Y <= k) = 1 - exp(-exp(alpha_k + eta_i))
            double p_upper = (yi_idx == m_K - 1) ? 1.0 : cloglog_inv(alpha[yi_idx] + eta[i]);
            double p_lower = (yi_idx == 0) ? 0.0 : cloglog_inv(alpha[yi_idx - 1] + eta[i]);
            
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

// [[Rcpp::export]]
List fast_ordinal_cloglog_regression_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int maxit = 100, double tol = 1e-6) {
    OrdinalCLLRegression model(X, y);
    int p = X.cols();
    int K = OrdinalCLLRegression::init_levels(y).size();
    if (K < 2) return List::create();
    int n_alpha = K - 1;
    int n_params = n_alpha + p;

    VectorXd params(n_params);
    // Initialize alpha
    for (int k = 0; k < n_alpha; ++k) {
        params[k] = -1.0 + 2.0 * (k + 1) / K;
    }
    // Initialize beta
    params.tail(p).setZero();

    bool converged = false;
    double h = 1e-4;
    for (int iter = 0; iter < maxit; ++iter) {
        VectorXd grad(n_params);
        for (int i = 0; i < n_params; ++i) {
            VectorXd p_plus = params; p_plus[i] += h;
            VectorXd p_minus = params; p_minus[i] -= h;
            grad[i] = (model.neg_log_likelihood(p_plus) - model.neg_log_likelihood(p_minus)) / (2.0 * h);
        }

        if (grad.norm() < tol) {
            converged = true;
            break;
        }

        MatrixXd H = model.hessian(params);
        FullPivLU<MatrixXd> lu(H);
        if (!lu.isInvertible()) break;

        VectorXd step = lu.solve(grad);
        
        // Simple line search
        double alpha_step = 1.0;
        double current_nll = model.neg_log_likelihood(params);
        bool accepted = false;
        while (alpha_step > 1e-4) {
            VectorXd next_params = params - alpha_step * step;
            if (model.neg_log_likelihood(next_params) < current_nll) {
                params = next_params;
                accepted = true;
                break;
            }
            alpha_step *= 0.5;
        }
        if (!accepted) {
            break;
        }
        if ((alpha_step * step).norm() < tol) {
            converged = true;
            break;
        }
    }

    return List::create(
        Named("b") = params.tail(p),
        Named("alpha") = params.head(n_alpha),
        Named("n_params") = n_params,
        Named("params") = params,
        Named("converged") = converged
    );
}

// [[Rcpp::export]]
List fast_ordinal_cloglog_regression_with_var_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) {
    List res = fast_ordinal_cloglog_regression_cpp(X, y);
    if (res.size() == 0) return List::create(Named("b") = NumericVector::create(NA_REAL), Named("ssq_b_2") = NA_REAL);
    
    VectorXd params = res["params"];
    bool converged = res["converged"];
    OrdinalCLLRegression model(X, y);
    MatrixXd H = model.hessian(params);
    
    FullPivLU<MatrixXd> lu(H);
    if (!lu.isInvertible()) {
        return List::create(Named("b") = res["b"], Named("ssq_b_2") = NA_REAL, Named("converged") = converged);
    }
    
    MatrixXd cov_mat = lu.inverse();
    int p = X.cols();
    int n_alpha = params.size() - p;
    
    double ssq_b_2 = (p >= 1) ? cov_mat(n_alpha, n_alpha) : NA_REAL;

    return List::create(
        Named("b") = res["b"],
        Named("ssq_b_2") = ssq_b_2,
        Named("converged") = converged
    );
}
