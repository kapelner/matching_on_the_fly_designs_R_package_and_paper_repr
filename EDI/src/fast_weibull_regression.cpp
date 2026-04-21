#include "_helper_functions.h"
#include <RcppEigen.h>
#include <Rmath.h>

using namespace Rcpp;

namespace {

class WeibullAFTLikelihood {
private:
    const Eigen::VectorXd& m_y;
    const Eigen::VectorXd& m_dead;
    const Eigen::MatrixXd& m_X;
    const int m_n;
    const int m_p;
    const Eigen::VectorXd m_log_y;

public:
    WeibullAFTLikelihood(const Eigen::VectorXd& y, 
                         const Eigen::VectorXd& dead, 
                         const Eigen::MatrixXd& X) :
        m_y(y), m_dead(dead), m_X(X), m_n(y.size()), m_p(X.cols()),
        m_log_y(y.array().log().matrix()) {}

    double operator()(const Eigen::VectorXd& params, Eigen::VectorXd& grad) {
        // params: [beta (p), log_sigma (1)]
        Eigen::VectorXd beta = params.head(m_p);
        double log_sigma = params[m_p];
        double sigma = std::exp(log_sigma);

        Eigen::VectorXd eta = m_X * beta;
        Eigen::VectorXd w = (m_log_y - eta) / sigma;
        Eigen::VectorXd exp_w(m_n);
        
        double loglik = 0.0;
        grad.setZero();
        Eigen::VectorXd d_ll_d_eta = Eigen::VectorXd::Zero(m_n);
        double d_ll_d_log_sigma = 0.0;

        for (int i = 0; i < m_n; ++i) {
            double wi = w[i];
            if (wi > 700.0) wi = 700.0;
            double ewi = std::exp(wi);
            exp_w[i] = ewi;

            if (m_dead[i] > 0.5) {
                // Event: log(f(y)) = w - log(sigma) - log(y) - exp(w)
                loglik += wi - log_sigma - m_log_y[i] - ewi;
                d_ll_d_eta[i] += (ewi - 1.0) / sigma;
                d_ll_d_log_sigma += (ewi - 1.0) * wi - 1.0;
            } else {
                // Censored: log(S(y)) = -exp(w)
                loglik -= ewi;
                d_ll_d_eta[i] += ewi / sigma;
                d_ll_d_log_sigma += ewi * wi;
            }
        }

        grad.head(m_p) = - m_X.transpose() * d_ll_d_eta;
        grad[m_p] = - d_ll_d_log_sigma;

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

// [[Rcpp::export]]
Eigen::VectorXd get_weibull_regression_score_cpp(const Eigen::VectorXd& y,
                                                 const Eigen::VectorXd& dead,
                                                 const Eigen::MatrixXd& X,
                                                 const Eigen::VectorXd& params) {
    WeibullAFTLikelihood fun(y, dead, X);
    Eigen::VectorXd grad(params.size());
    fun(params, grad);
    return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_weibull_regression_hessian_cpp(const Eigen::VectorXd& y,
                                                   const Eigen::VectorXd& dead,
                                                   const Eigen::MatrixXd& X,
                                                   const Eigen::VectorXd& params) {
    WeibullAFTLikelihood fun(y, dead, X);
    return -fun.hessian(params);
}

// [[Rcpp::export]]
List fast_weibull_regression_cpp(const Eigen::VectorXd& y, 
                                 const Eigen::VectorXd& dead, 
                                 const Eigen::MatrixXd& X, 
                                 Nullable<NumericVector> start_params = R_NilValue,
                                 bool estimate_only = false,
                                 int maxit = 1000, 
                                 double tol = 1e-6,
                                 Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                 Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                 std::string optimization_alg = "lbfgs") {
    int p = X.cols();
    Eigen::VectorXd params(p + 1);
    
    if (start_params.isNotNull()) {
        params = as<Eigen::VectorXd>(NumericVector(start_params));
    } else {
        // Initial OLS on log(y) for beta, and rough estimate for log_sigma
        Eigen::VectorXd log_y = y.array().log().matrix();
        Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(X);
        params.head(p) = cod.solve(log_y);
        Eigen::VectorXd resid = log_y - X * params.head(p);
        double std_resid = std::sqrt(resid.squaredNorm() / (y.size() - p));
        params[p] = std::log(std_resid * 0.7797); // 0.7797 is ~sqrt(6)/pi
    }
    FixedParamSpec fixed_spec = make_fixed_param_spec(p + 1, fixed_idx, fixed_values);

    WeibullAFTLikelihood fun(y, dead, X);
    LikelihoodFitResult fit;
    try {
        fit = optimize_fixed_likelihood(fun, params, fixed_spec, maxit, tol, optimization_alg, "lbfgs");
    } catch (...) {
        return List::create(Named("converged") = false);
    }
    params = fit.params;

    if (estimate_only) {
        return List::create(
            Named("coefficients") = params.head(p),
            Named("log_sigma") = params[p],
            Named("converged") = fit.converged,
            Named("neg_ll") = fit.value
        );
    }

    Eigen::MatrixXd H = fun.hessian(params);
    Eigen::MatrixXd H_free = subset_matrix(H, fixed_spec.free_idx, fixed_spec.free_idx);
    Eigen::MatrixXd cov_free = H_free.inverse();
    Eigen::MatrixXd vcov = expand_free_covariance(p + 1, fixed_spec, cov_free, true);

    return List::create(
        Named("coefficients") = params.head(p),
        Named("log_sigma") = params[p],
        Named("vcov") = vcov,
        Named("converged") = fit.converged,
        Named("neg_ll") = fit.value
    );
}
