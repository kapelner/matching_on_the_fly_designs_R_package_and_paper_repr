#include <RcppEigen.h>
#include <optimization/LBFGS.h>
#include <Rmath.h>
#include "_helper_functions.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace LBFGSpp;

namespace {

// Helper for log-sum-exp trick to compute log(exp(a) + exp(b) - exp(c))
// Used for Clayton Copula logA: log(exp(theta*H1) + exp(theta*H2) - 1)
// where 1 = exp(0).
inline double log_sum_exp_clayton(double a, double b) {
    double m = a;
    if (b > m) m = b;
    if (0.0 > m) m = 0.0;
    double inner = std::exp(a - m) + std::exp(b - m) - std::exp(-m);
    if (inner <= 0) return m + std::log(std::numeric_limits<double>::min());
    return m + std::log(inner);
}

// -----------------------------------------------------------------------------
// Clayton Copula Weibull AFT
// -----------------------------------------------------------------------------

class ClaytonWeibullLikelihood {
private:
    const Eigen::VectorXd& m_y;
    const Eigen::VectorXd& m_dead;
    const Eigen::MatrixXd& m_X;
    const Eigen::MatrixXi& m_pair_idx;
    const Eigen::VectorXi& m_singleton_rows;
    const bool m_has_pairs;
    const bool m_has_singletons;
    const int m_n;
    const int m_p;
    const Eigen::VectorXd m_log_y;

public:
    ClaytonWeibullLikelihood(const Eigen::VectorXd& y, 
                             const Eigen::VectorXd& dead, 
                             const Eigen::MatrixXd& X,
                             const Eigen::MatrixXi& pair_idx,
                             const Eigen::VectorXi& singleton_rows) :
        m_y(y), m_dead(dead), m_X(X), m_pair_idx(pair_idx), 
        m_singleton_rows(singleton_rows),
        m_has_pairs(pair_idx.rows() > 0),
        m_has_singletons(singleton_rows.size() > 0),
        m_n(y.size()), m_p(X.cols()),
        m_log_y(y.array().log().matrix()) {}

    double operator()(const Eigen::VectorXd& params, Eigen::VectorXd& grad) {
        double log_sigma = params[m_p];
        double log_theta = params[m_p + 1];
        double sigma = std::exp(log_sigma);
        double theta = std::exp(log_theta);
        Eigen::VectorXd beta = params.head(m_p);

        Eigen::VectorXd eta = m_X * beta;
        Eigen::VectorXd log_H = (m_log_y - eta) / sigma;
        Eigen::VectorXd H(m_n);
        Eigen::VectorXd log_f(m_n);
        
        for (int i = 0; i < m_n; ++i) {
            double lH = log_H[i];
            if (lH > 700.0) lH = 700.0;
            H[i] = std::exp(lH);
            log_f[i] = lH - log_sigma - m_log_y[i] - H[i];
        }

        double loglik = 0.0;
        grad.setZero();
        
        Eigen::VectorXd d_loglik_d_eta = Eigen::VectorXd::Zero(m_n);
        double d_loglik_d_log_sigma = 0.0;
        double d_loglik_d_log_theta = 0.0;

        if (m_has_pairs) {
            for (int k = 0; k < m_pair_idx.rows(); ++k) {
                int i1 = m_pair_idx(k, 0);
                int i2 = m_pair_idx(k, 1);
                double h1 = H[i1];
                double h2 = H[i2];
                double d1 = m_dead[i1];
                double d2 = m_dead[i2];
                
                double logA = log_sum_exp_clayton(theta * h1, theta * h2);
                double A = std::exp(logA);
                
                // Common terms for derivatives
                double dA_d_h1 = theta * std::exp(theta * h1);
                double dA_d_h2 = theta * std::exp(theta * h2);
                double dA_d_theta = h1 * std::exp(theta * h1) + h2 * std::exp(theta * h2);

                if (d1 < 0.5 && d2 < 0.5) { // mask00
                    loglik -= (1.0 / theta) * logA;
                    // d/d_theta (-1/theta * logA) = 1/theta^2 * logA - 1/theta * 1/A * dA_d_theta
                    d_loglik_d_log_theta += (logA / (theta * theta) - dA_d_theta / (theta * A)) * theta;
                    
                    double d_ll_d_A = -1.0 / (theta * A);
                    d_loglik_d_eta[i1] += d_ll_d_A * dA_d_h1 * (-H[i1] / sigma);
                    d_loglik_d_eta[i2] += d_ll_d_A * dA_d_h2 * (-H[i2] / sigma);
                    d_loglik_d_log_sigma += (d_ll_d_A * dA_d_h1 * (-H[i1] * log_H[i1]) + 
                                             d_ll_d_A * dA_d_h2 * (-H[i2] * log_H[i2]));
                } else if (d1 > 0.5 && d2 < 0.5) { // mask10
                    loglik += log_f[i1] + (-1.0/theta - 1.0) * logA + (theta + 1.0) * h1;
                    
                    d_loglik_d_log_theta += (logA / (theta * theta) + (-1.0/theta - 1.0) * dA_d_theta / A + h1) * theta;
                    
                    double d_ll_d_h1 = (theta + 1.0) + (-1.0/theta - 1.0) * dA_d_h1 / A;
                    double d_ll_d_h2 = (-1.0/theta - 1.0) * dA_d_h2 / A;
                    
                    // From log_f[i1]: d/d_eta = (H[i1] - 1)/sigma, d/d_log_sigma = -1 + (H[i1] - 1)*log_H[i1]
                    d_loglik_d_eta[i1] += (H[i1] - 1.0) / sigma + d_ll_d_h1 * (-H[i1] / sigma);
                    d_loglik_d_eta[i2] += d_ll_d_h2 * (-H[i2] / sigma);
                    d_loglik_d_log_sigma += -1.0 + (H[i1] - 1.0) * log_H[i1] + 
                                            d_ll_d_h1 * (-H[i1] * log_H[i1]) + 
                                            d_ll_d_h2 * (-H[i2] * log_H[i2]);
                } else if (d1 < 0.5 && d2 > 0.5) { // mask01
                    loglik += log_f[i2] + (-1.0/theta - 1.0) * logA + (theta + 1.0) * h2;
                    
                    d_loglik_d_log_theta += (logA / (theta * theta) + (-1.0/theta - 1.0) * dA_d_theta / A + h2) * theta;
                    
                    double d_ll_d_h1 = (-1.0/theta - 1.0) * dA_d_h1 / A;
                    double d_ll_d_h2 = (theta + 1.0) + (-1.0/theta - 1.0) * dA_d_h2 / A;
                    
                    d_loglik_d_eta[i1] += d_ll_d_h1 * (-H[i1] / sigma);
                    d_loglik_d_eta[i2] += (H[i2] - 1.0) / sigma + d_ll_d_h2 * (-H[i2] / sigma);
                    d_loglik_d_log_sigma += -1.0 + (H[i2] - 1.0) * log_H[i2] + 
                                            d_ll_d_h1 * (-H[i1] * log_H[i1]) + 
                                            d_ll_d_h2 * (-H[i2] * log_H[i2]);
                } else { // mask11
                    loglik += std::log(theta + 1.0) + log_f[i1] + log_f[i2] + (-1.0/theta - 2.0) * logA + (theta + 1.0) * (h1 + h2);
                    
                    d_loglik_d_log_theta += (1.0/(theta + 1.0) + (1.0/(theta*theta)) * logA + (-1.0/theta - 2.0) * dA_d_theta / A + h1 + h2) * theta;
                    
                    double d_ll_d_h1 = (theta + 1.0) + (-1.0/theta - 2.0) * dA_d_h1 / A;
                    double d_ll_d_h2 = (theta + 1.0) + (-1.0/theta - 2.0) * dA_d_h2 / A;
                    
                    d_loglik_d_eta[i1] += (H[i1] - 1.0) / sigma + d_ll_d_h1 * (-H[i1] / sigma);
                    d_loglik_d_eta[i2] += (H[i2] - 1.0) / sigma + d_ll_d_h2 * (-H[i2] / sigma);
                    d_loglik_d_log_sigma += -2.0 + (H[i1] - 1.0) * log_H[i1] + (H[i2] - 1.0) * log_H[i2] + 
                                            d_ll_d_h1 * (-H[i1] * log_H[i1]) + 
                                            d_ll_d_h2 * (-H[i2] * log_H[i2]);
                }
            }
        }

        if (m_has_singletons) {
            for (int k = 0; k < m_singleton_rows.size(); ++k) {
                int i = m_singleton_rows[k];
                if (m_dead[i] > 0.5) {
                    loglik += log_f[i];
                    d_loglik_d_eta[i] += (H[i] - 1.0) / sigma;
                    d_loglik_d_log_sigma += -1.0 + (H[i] - 1.0) * log_H[i];
                } else {
                    loglik -= H[i];
                    d_loglik_d_eta[i] += H[i] / sigma;
                    d_loglik_d_log_sigma += H[i] * log_H[i];
                }
            }
        }

        grad.head(m_p) = - m_X.transpose() * d_loglik_d_eta;
        grad[m_p] = - d_loglik_d_log_sigma;
        grad[m_p + 1] = - d_loglik_d_log_theta;

        return -loglik;
    }
    
    Eigen::MatrixXd hessian(const Eigen::VectorXd& params) {
        int total_p = params.size();
        Eigen::MatrixXd H(total_p, total_p);
        H.setZero();
        double h = 1e-6;
        Eigen::VectorXd grad_at_params(total_p);
        operator()(params, grad_at_params);

        for (int i = 0; i < total_p; ++i) {
            Eigen::VectorXd p_plus = params;
            p_plus[i] += h;
            Eigen::VectorXd g_plus(total_p);
            operator()(p_plus, g_plus);
            H.col(i) = (g_plus - grad_at_params) / h;
        }
        H = (H + H.transpose()) / 2.0;
        return H;
    }
};

// -----------------------------------------------------------------------------
// Dependent Censoring Transformation Regression
// -----------------------------------------------------------------------------

class DepCensTransformLikelihood {
private:
    const Eigen::VectorXd& m_y;
    const Eigen::VectorXd& m_dead;
    const Eigen::MatrixXd& m_X;
    const int m_n;
    const int m_p;
    const Eigen::VectorXd m_log_y;

public:
    DepCensTransformLikelihood(const Eigen::VectorXd& y, 
                               const Eigen::VectorXd& dead, 
                               const Eigen::MatrixXd& X) :
        m_y(y), m_dead(dead), m_X(X), m_n(y.size()), m_p(X.cols()),
        m_log_y(y.array().log().matrix()) {}

    double operator()(const Eigen::VectorXd& params, Eigen::VectorXd& grad) {
        // params: [beta_event (p), beta_cens (p), log_sigma_event (1), log_sigma_cens (1), atanh_rho (1)]
        Eigen::VectorXd beta_event = params.head(m_p);
        Eigen::VectorXd beta_cens = params.segment(m_p, m_p);
        double log_sigma_event = params[2 * m_p];
        double log_sigma_cens = params[2 * m_p + 1];
        double atanh_rho = params[2 * m_p + 2];

        double sigma_event = std::exp(log_sigma_event);
        double sigma_cens = std::exp(log_sigma_cens);
        double rho = std::tanh(atanh_rho);
        double one_minus_rho_sq = std::max(1.0 - rho * rho, 1e-12);
        double sd_cond = std::sqrt(one_minus_rho_sq);

        Eigen::VectorXd mu_event = m_X * beta_event;
        Eigen::VectorXd mu_cens = m_X * beta_cens;
        Eigen::VectorXd z_event = (m_log_y - mu_event) / sigma_event;
        Eigen::VectorXd z_cens = (m_log_y - mu_cens) / sigma_cens;

        double loglik = 0.0;
        grad.setZero();
        
        Eigen::VectorXd d_ll_d_mu_event = Eigen::VectorXd::Zero(m_n);
        Eigen::VectorXd d_ll_d_mu_cens = Eigen::VectorXd::Zero(m_n);
        double d_ll_d_log_sigma_event = 0.0;
        double d_ll_d_log_sigma_cens = 0.0;
        double d_ll_d_atanh_rho = 0.0;

        for (int i = 0; i < m_n; ++i) {
            double ze = z_event[i];
            double zc = z_cens[i];
            double d = m_dead[i];

            double log_f_event = R::dnorm(ze, 0.0, 1.0, 1) - log_sigma_event - m_log_y[i];
            double log_f_cens = R::dnorm(zc, 0.0, 1.0, 1) - log_sigma_cens - m_log_y[i];

            double w_event = (rho * ze - zc) / sd_cond;
            double w_cens = (rho * zc - ze) / sd_cond;
            
            double log_S_cond_cens = R::pnorm(w_event, 0.0, 1.0, 1, 1);
            double log_S_cond_event = R::pnorm(w_cens, 0.0, 1.0, 1, 1);

            double li = d * (log_f_event + log_S_cond_cens) + (1.0 - d) * (log_f_cens + log_S_cond_event);
            loglik += li;

            // Derivatives
            double mill_event = std::exp(R::dnorm(w_event, 0.0, 1.0, 1) - log_S_cond_cens);
            double mill_cens = std::exp(R::dnorm(w_cens, 0.0, 1.0, 1) - log_S_cond_event);

            // d/d_ze w_event = rho / sd_cond
            // d/d_zc w_event = -1 / sd_cond
            // d/d_rho w_event = (ze * sd_cond - (rho*ze - zc) * (-rho/sd_cond)) / (1-rho^2)
            //                = (ze * (1-rho^2) + rho*(rho*ze - zc)) / (1-rho^2)^(3/2)
            //                = (ze - rho^2*ze + rho^2*ze - rho*zc) / (1-rho^2)^(3/2)
            //                = (ze - rho*zc) / (1-rho^2)^(3/2)
            
            double d_w_event_d_rho = (ze - rho * zc) / std::pow(one_minus_rho_sq, 1.5);
            double d_w_cens_d_rho = (zc - rho * ze) / std::pow(one_minus_rho_sq, 1.5);

            if (d > 0.5) {
                // d * log_f_event
                d_ll_d_mu_event[i] += ze / sigma_event;
                d_ll_d_log_sigma_event += (ze * ze - 1.0);
                
                // d * log_S_cond_cens
                d_ll_d_mu_event[i] += mill_event * (rho / sd_cond) * (-1.0 / sigma_event);
                d_ll_d_mu_cens[i] += mill_event * (-1.0 / sd_cond) * (-1.0 / sigma_cens);
                d_ll_d_log_sigma_event += mill_event * (rho / sd_cond) * (-ze);
                d_ll_d_log_sigma_cens += mill_event * (-1.0 / sd_cond) * (-zc);
                d_ll_d_atanh_rho += mill_event * d_w_event_d_rho * one_minus_rho_sq; // d_rho/d_atanh_rho = 1-rho^2
            } else {
                // (1-d) * log_f_cens
                d_ll_d_mu_cens[i] += zc / sigma_cens;
                d_ll_d_log_sigma_cens += (zc * zc - 1.0);

                // (1-d) * log_S_cond_event
                d_ll_d_mu_cens[i] += mill_cens * (rho / sd_cond) * (-1.0 / sigma_cens);
                d_ll_d_mu_event[i] += mill_cens * (-1.0 / sd_cond) * (-1.0 / sigma_event);
                d_ll_d_log_sigma_cens += mill_cens * (rho / sd_cond) * (-zc);
                d_ll_d_log_sigma_event += mill_cens * (-1.0 / sd_cond) * (-ze);
                d_ll_d_atanh_rho += mill_cens * d_w_cens_d_rho * one_minus_rho_sq;
            }
        }

        grad.head(m_p) = - m_X.transpose() * d_ll_d_mu_event;
        grad.segment(m_p, m_p) = - m_X.transpose() * d_ll_d_mu_cens;
        grad[2 * m_p] = - d_ll_d_log_sigma_event;
        grad[2 * m_p + 1] = - d_ll_d_log_sigma_cens;
        grad[2 * m_p + 2] = - d_ll_d_atanh_rho;

        return -loglik;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& params) {
        int total_p = params.size();
        Eigen::MatrixXd H(total_p, total_p);
        H.setZero();
        double h = 1e-6;
        Eigen::VectorXd grad_at_params(total_p);
        operator()(params, grad_at_params);

        for (int i = 0; i < total_p; ++i) {
            Eigen::VectorXd p_plus = params;
            p_plus[i] += h;
            Eigen::VectorXd g_plus(total_p);
            operator()(p_plus, g_plus);
            H.col(i) = (g_plus - grad_at_params) / h;
        }
        H = (H + H.transpose()) / 2.0;
        return H;
    }
};

} // namespace

// -----------------------------------------------------------------------------
// R-exposed functions
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List fast_clayton_weibull_aft_optim_cpp(
    const Eigen::VectorXd& y,
    const Eigen::VectorXd& dead,
    const Eigen::MatrixXd& X,
    const Eigen::MatrixXi& pair_idx,
    const Eigen::VectorXi& singleton_rows,
    const Eigen::VectorXd& start_params,
    int maxit = 2000,
    double reltol = 1e-9
) {
    ClaytonWeibullLikelihood fun(y, dead, X, pair_idx, singleton_rows);
    
    LBFGSParam<double> lbfgs_params;
    lbfgs_params.epsilon = reltol;
    lbfgs_params.max_iterations = maxit;
    
    LBFGSSolver<double> solver(lbfgs_params);
    Eigen::VectorXd params = start_params;
    double neg_ll;
    int niter = 0;
    try {
        niter = solver.minimize(fun, params, neg_ll);
    } catch (...) {
        return List::create(Named("converged") = false);
    }

    return List::create(
        Named("par") = params,
        Named("value") = neg_ll,
        Named("niter") = niter,
        Named("converged") = (niter < maxit),
        Named("hessian") = fun.hessian(params)
    );
}

// [[Rcpp::export]]
List fast_dep_cens_transform_optim_cpp(
    const Eigen::VectorXd& y,
    const Eigen::VectorXd& dead,
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& start_params,
    int maxit = 2000,
    double reltol = 1e-9
) {
    DepCensTransformLikelihood fun(y, dead, X);
    
    LBFGSParam<double> lbfgs_params;
    lbfgs_params.epsilon = reltol;
    lbfgs_params.max_iterations = maxit;
    
    LBFGSSolver<double> solver(lbfgs_params);
    Eigen::VectorXd params = start_params;
    double neg_ll;
    int niter = 0;
    try {
        niter = solver.minimize(fun, params, neg_ll);
    } catch (...) {
        return List::create(Named("converged") = false);
    }

    return List::create(
        Named("par") = params,
        Named("value") = neg_ll,
        Named("niter") = niter,
        Named("converged") = (niter < maxit),
        Named("hessian") = fun.hessian(params)
    );
}
