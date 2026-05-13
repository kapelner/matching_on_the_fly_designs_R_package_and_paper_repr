#include "_helper_functions.h"
#include <RcppEigen.h>
#include <Rmath.h>

using namespace Rcpp;

namespace {

// Helper for log1mexp(x) = log(1 - exp(x)) for x < 0
double log1mexp(double x) {
    if (x >= 0) return -std::numeric_limits<double>::infinity();
    if (x > -0.6931471805599453) return std::log(-std::expm1(x));
    return std::log1p(-std::exp(x));
}

// Helper for log1pexp(x) = log(1 + exp(x))
double log1pexp(double x) {
    if (x > 0) return x + std::log1p(std::exp(-x));
    return std::log1p(std::exp(x));
}

class ZeroAugmentedPoisson {
private:
    const Eigen::VectorXd m_y;
    const Eigen::MatrixXd m_X;
    const Eigen::MatrixXd m_Xzi;
    const int m_n;
    const int m_p_cond;
    const int m_p_zi;
    const bool m_is_hurdle;

public:
    ZeroAugmentedPoisson(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, const Eigen::MatrixXd& Xzi, bool is_hurdle) :
        m_y(y), m_X(X), m_Xzi(Xzi), m_n(X.rows()), 
        m_p_cond(X.cols()), m_p_zi(Xzi.cols()), m_is_hurdle(is_hurdle) {}

    double operator()(const Eigen::VectorXd& params, Eigen::VectorXd& grad) {
        Eigen::VectorXd beta_cond = params.head(m_p_cond);
        Eigen::VectorXd beta_zi = params.tail(m_p_zi);

        Eigen::VectorXd eta_cond = m_X * beta_cond;
        Eigen::VectorXd eta_zi = m_Xzi * beta_zi;

        double neg_ll = 0;
        grad.setZero();
        Eigen::VectorXd grad_cond = grad.head(m_p_cond);
        Eigen::VectorXd grad_zi = grad.tail(m_p_zi);

        for (int i = 0; i < m_n; ++i) {
            double lambda = std::exp(std::min(eta_cond[i], 700.0));
            double pi = 1.0 / (1.0 + std::exp(-std::max(std::min(eta_zi[i], 700.0), -700.0)));
            // pi is P(structural zero) for ZIP, P(zero) for Hurdle

            if (m_is_hurdle) {
                // Hurdle Poisson
                if (m_y[i] == 0) {
                    neg_ll -= std::log(std::max(pi, 1e-15));
                    // d/d_eta_zi log(pi) = 1 - pi
                    grad_zi -= m_Xzi.row(i).transpose() * (1.0 - pi);
                } else {
                    double log_1_minus_pi = -log1pexp(eta_zi[i]);
                    double log_tp = m_y[i] * std::min(eta_cond[i], 700.0) - lambda - log1mexp(-lambda);
                    neg_ll -= (log_1_minus_pi + log_tp);
                    
                    // d/d_eta_zi log(1-pi) = -pi
                    grad_zi += m_Xzi.row(i).transpose() * pi;
                    
                    // d/d_eta_cond log_tp = y - lambda / (1 - exp(-lambda))
                    double exp_ml = std::exp(-lambda);
                    double d_tp = m_y[i] - lambda / (1.0 - exp_ml);
                    grad_cond -= m_X.row(i).transpose() * d_tp;
                }
            } else {
                // Zero-Inflated Poisson
                if (m_y[i] == 0) {
                    double exp_ml = std::exp(-lambda);
                    double prob_zero = pi + (1.0 - pi) * exp_ml;
                    neg_ll -= std::log(std::max(prob_zero, 1e-15));
                    
                    // d/d_eta_zi log(prob_zero) = [pi(1-pi)(1-exp(-lambda))] / prob_zero
                    double d_zi = (pi * (1.0 - pi) * (1.0 - exp_ml)) / prob_zero;
                    grad_zi -= m_Xzi.row(i).transpose() * d_zi;
                    
                    // d/d_eta_cond log(prob_zero) = [(1-pi) * (-lambda * exp(-lambda))] / prob_zero
                    double d_cond = ((1.0 - pi) * (-lambda * exp_ml)) / prob_zero;
                    grad_cond -= m_X.row(i).transpose() * d_cond;
                } else {
                    double log_1_minus_pi = -log1pexp(eta_zi[i]);
                    double log_poisson = m_y[i] * std::min(eta_cond[i], 700.0) - lambda; // ignore factorial
                    neg_ll -= (log_1_minus_pi + log_poisson);
                    
                    // d/d_eta_zi log(1-pi) = -pi
                    grad_zi += m_Xzi.row(i).transpose() * pi;
                    
                    // d/d_eta_cond log_poisson = y - lambda
                    grad_cond -= m_X.row(i).transpose() * (m_y[i] - lambda);
                }
            }
        }
        
        grad.head(m_p_cond) = grad_cond;
        grad.tail(m_p_zi) = grad_zi;

        return neg_ll;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& params) {
        int total_p = m_p_cond + m_p_zi;
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total_p, total_p);
        Eigen::VectorXd beta_cond = params.head(m_p_cond);
        Eigen::VectorXd beta_zi = params.tail(m_p_zi);
        Eigen::VectorXd eta_cond = m_X * beta_cond;
        Eigen::VectorXd eta_zi = m_Xzi * beta_zi;

        for (int i = 0; i < m_n; ++i) {
            double eta_c = std::min(eta_cond[i], 700.0);
            double lambda = std::exp(eta_c);
            double eta_z = std::max(std::min(eta_zi[i], 700.0), -700.0);
            double pi = 1.0 / (1.0 + std::exp(-eta_z));
            double pi_prime = pi * (1.0 - pi);
            Eigen::VectorXd xc = m_X.row(i).transpose();
            Eigen::VectorXd xz = m_Xzi.row(i).transpose();

            if (m_is_hurdle) {
                double h_zi = pi_prime;
                H.bottomRightCorner(m_p_zi, m_p_zi).noalias() += h_zi * (xz * xz.transpose());
                if (m_y[i] > 0) {
                    double exp_ml = std::exp(-lambda);
                    double denom = std::max(1.0 - exp_ml, 1e-15);
                    double h_cond = lambda * (denom - lambda * exp_ml) / (denom * denom);
                    H.topLeftCorner(m_p_cond, m_p_cond).noalias() += h_cond * (xc * xc.transpose());
                }
            } else {
                if (m_y[i] == 0) {
                    double r = std::exp(-lambda);
                    double q = std::max(pi + (1.0 - pi) * r, 1e-15);
                    double pi_second = pi_prime * (1.0 - 2.0 * pi);
                    double q_z = pi_prime * (1.0 - r);
                    double q_zz = pi_second * (1.0 - r);
                    double q_c = -(1.0 - pi) * lambda * r;
                    double q_cc = (1.0 - pi) * r * (lambda * lambda - lambda);
                    double q_zc = pi_prime * lambda * r;

                    double h_cc = q_c * q_c / (q * q) - q_cc / q;
                    double h_zz = q_z * q_z / (q * q) - q_zz / q;
                    double h_cz = q_c * q_z / (q * q) - q_zc / q;

                    H.topLeftCorner(m_p_cond, m_p_cond).noalias() += h_cc * (xc * xc.transpose());
                    H.bottomRightCorner(m_p_zi, m_p_zi).noalias() += h_zz * (xz * xz.transpose());
                    H.topRightCorner(m_p_cond, m_p_zi).noalias() += h_cz * (xc * xz.transpose());
                } else {
                    H.topLeftCorner(m_p_cond, m_p_cond).noalias() += lambda * (xc * xc.transpose());
                    H.bottomRightCorner(m_p_zi, m_p_zi).noalias() += pi_prime * (xz * xz.transpose());
                }
            }
        }
        H.bottomLeftCorner(m_p_zi, m_p_cond) = H.topRightCorner(m_p_cond, m_p_zi).transpose();
        return H;
    }
};

} // namespace

// [[Rcpp::export]]
Eigen::VectorXd get_zero_augmented_poisson_score_cpp(const Eigen::MatrixXd& X,
													 const Eigen::VectorXd& y,
													 const Eigen::MatrixXd& Xzi,
													 const Eigen::VectorXd& params,
													 bool is_hurdle) {
	ZeroAugmentedPoisson fun(y, X, Xzi, is_hurdle);
	Eigen::VectorXd grad(params.size());
	fun(params, grad);
	return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_zero_augmented_poisson_hessian_cpp(const Eigen::MatrixXd& X,
													   const Eigen::VectorXd& y,
													   const Eigen::MatrixXd& Xzi,
													   const Eigen::VectorXd& params,
													   bool is_hurdle) {
	ZeroAugmentedPoisson fun(y, X, Xzi, is_hurdle);
	return -fun.hessian(params);
}

// [[Rcpp::export]]
List fast_zero_augmented_poisson_cpp(const Eigen::MatrixXd& X,
                                     const Eigen::VectorXd& y,
                                     const Eigen::MatrixXd& Xzi,
                                     bool is_hurdle,
                                     Nullable<NumericVector> start_params = R_NilValue,
                                     bool estimate_only = false,
                                     int maxit = 1000,
                                     double tol = 1e-6,
                                     Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                     std::string optimization_alg = "newton_raphson",
                                     Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    int p_cond = X.cols();
    int p_zi = Xzi.cols();
    int total_p = p_cond + p_zi;

    Eigen::VectorXd params(total_p);
    if (start_params.isNotNull()) {
        params = as<Eigen::VectorXd>(NumericVector(start_params));
    } else {
        params.setZero();
        // Naive start for intercept
        double mean_y = y.mean();
        if (mean_y > 0) params[0] = std::log(mean_y);

        double prop_zero = 0;
        for(int i=0; i<y.size(); ++i) if(y[i] == 0) prop_zero++;
        prop_zero /= y.size();
            if (prop_zero > 0 && prop_zero < 1) params[p_cond] = std::log(prop_zero / (1.0 - prop_zero));
    }
    FixedParamSpec fixed_spec = make_fixed_param_spec(total_p, fixed_idx, fixed_values);

    ZeroAugmentedPoisson fun(y, X, Xzi, is_hurdle);
    Eigen::MatrixXd H_start;
    const Eigen::MatrixXd* h_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        H_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        h_ptr = &H_start;
    }

    LikelihoodFitResult fit;
    try {
        fit = optimize_fixed_likelihood(fun, params, fixed_spec, maxit, tol, optimization_alg, "newton_raphson", 0, h_ptr);
    } catch (...) {
        return List::create(Named("converged") = false);
    }
    params = fit.params;

    if (estimate_only) {
        return List::create(
            Named("coefficients") = List::create(
                Named("cond") = params.head(p_cond),
                Named("zi") = params.tail(p_zi)
            ),
            Named("converged") = fit.converged,
            Named("neg_ll") = fit.value
        );
    }

    Eigen::MatrixXd H = fun.hessian(params);
    Eigen::MatrixXd H_free = subset_matrix(H, fixed_spec.free_idx, fixed_spec.free_idx);
    Eigen::MatrixXd cov_free = H_free.inverse();
    Eigen::MatrixXd vcov = expand_free_covariance(total_p, fixed_spec, cov_free, true);

    return List::create(
        Named("coefficients") = List::create(
            Named("cond") = params.head(p_cond),
            Named("zi") = params.tail(p_zi)
        ),
        Named("vcov") = vcov,
        Named("converged") = fit.converged,
        Named("neg_ll") = fit.value
    );
}
