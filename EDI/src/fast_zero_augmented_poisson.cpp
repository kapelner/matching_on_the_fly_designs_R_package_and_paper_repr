#include "_helper_functions.h"
#include <RcppEigen.h>
#include <optimization/LBFGS.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace LBFGSpp;

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
    const Eigen::MatrixXd m_Xcond;
    const Eigen::MatrixXd m_Xzi;
    const int m_n;
    const int m_p_cond;
    const int m_p_zi;
    const bool m_is_hurdle;

public:
    ZeroAugmentedPoisson(const Eigen::VectorXd& y, const Eigen::MatrixXd& Xcond, const Eigen::MatrixXd& Xzi, bool is_hurdle) :
        m_y(y), m_Xcond(Xcond), m_Xzi(Xzi), m_n(Xcond.rows()), 
        m_p_cond(Xcond.cols()), m_p_zi(Xzi.cols()), m_is_hurdle(is_hurdle) {}

    double operator()(const Eigen::VectorXd& params, Eigen::VectorXd& grad) {
        Eigen::VectorXd beta_cond = params.head(m_p_cond);
        Eigen::VectorXd beta_zi = params.tail(m_p_zi);

        Eigen::VectorXd eta_cond = m_Xcond * beta_cond;
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
                    grad_cond -= m_Xcond.row(i).transpose() * d_tp;
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
                    grad_cond -= m_Xcond.row(i).transpose() * d_cond;
                } else {
                    double log_1_minus_pi = -log1pexp(eta_zi[i]);
                    double log_poisson = m_y[i] * std::min(eta_cond[i], 700.0) - lambda; // ignore factorial
                    neg_ll -= (log_1_minus_pi + log_poisson);
                    
                    // d/d_eta_zi log(1-pi) = -pi
                    grad_zi += m_Xzi.row(i).transpose() * pi;
                    
                    // d/d_eta_cond log_poisson = y - lambda
                    grad_cond -= m_Xcond.row(i).transpose() * (m_y[i] - lambda);
                }
            }
        }
        
        grad.head(m_p_cond) = grad_cond;
        grad.tail(m_p_zi) = grad_zi;

        return neg_ll;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& params) {
        int total_p = m_p_cond + m_p_zi;
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

// [[Rcpp::export]]
List fast_zero_augmented_poisson_cpp(const Eigen::VectorXd& y, 
                                     const Eigen::MatrixXd& Xcond, 
                                     const Eigen::MatrixXd& Xzi, 
                                     bool is_hurdle,
                                     Nullable<NumericVector> start_params = R_NilValue,
                                     bool estimate_only = false,
                                     int maxit = 1000, 
                                     double tol = 1e-6) {
    int p_cond = Xcond.cols();
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

    ZeroAugmentedPoisson fun(y, Xcond, Xzi, is_hurdle);
    LBFGSParam<double> lbfgs_params;
    lbfgs_params.epsilon = tol;
    lbfgs_params.max_iterations = maxit;

    LBFGSSolver<double> solver(lbfgs_params);
    double neg_ll;
    int niter = 0;
    try {
        niter = solver.minimize(fun, params, neg_ll);
    } catch (...) {
        return List::create(Named("converged") = false);
    }

    if (estimate_only) {
        return List::create(
            Named("coefficients") = List::create(
                Named("cond") = params.head(p_cond),
                Named("zi") = params.tail(p_zi)
            ),
            Named("converged") = (niter < maxit),
            Named("neg_ll") = neg_ll
        );
    }

    Eigen::MatrixXd H = fun.hessian(params);
    Eigen::MatrixXd vcov = H.inverse();

    return List::create(
        Named("coefficients") = List::create(
            Named("cond") = params.head(p_cond),
            Named("zi") = params.tail(p_zi)
        ),
        Named("vcov") = vcov,
        Named("converged") = (niter < maxit),
        Named("neg_ll") = neg_ll
    );
}
