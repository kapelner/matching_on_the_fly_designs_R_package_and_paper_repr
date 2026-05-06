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
        const Eigen::ArrayXd w = ((m_log_y - eta) / sigma).array().min(700.0);
        const Eigen::ArrayXd exp_w = w.exp();
        const Eigen::ArrayXd dead = m_dead.array();
        const double loglik = (dead * (w - log_sigma - m_log_y.array()) - exp_w).sum();
        const Eigen::VectorXd d_ll_d_eta = ((exp_w - dead) / sigma).matrix();
        const double d_ll_d_log_sigma = (exp_w * w - dead * (w + 1.0)).sum();

        grad.setZero();

        grad.head(m_p) = - m_X.transpose() * d_ll_d_eta;
        grad[m_p] = - d_ll_d_log_sigma;

        return -loglik;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& params) {
        int total_p = params.size();
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total_p, total_p);
        Eigen::VectorXd beta = params.head(m_p);
        double sigma = std::exp(params[m_p]);
        Eigen::VectorXd eta = m_X * beta;
        const Eigen::ArrayXd w = ((m_log_y - eta) / sigma).array().min(700.0);
        const Eigen::ArrayXd exp_w = w.exp();
        const Eigen::VectorXd beta_weights = (exp_w / (sigma * sigma)).matrix();
        const Eigen::VectorXd cross_weights =
            ((exp_w * (w + 1.0) - m_dead.array()) / sigma).matrix();

        H.topLeftCorner(m_p, m_p).noalias() = weighted_crossprod(m_X, beta_weights);
        H.topRightCorner(m_p, 1).noalias() = m_X.transpose() * cross_weights;
        H(m_p, m_p) = (exp_w * (w.square() + w) - m_dead.array() * w).sum();
        H.bottomLeftCorner(1, m_p) = H.topRightCorner(m_p, 1).transpose();
        return H;
    }
};

} // namespace

//' @title Compute Weibull Regression Score
//' @description Calculates the score vector (gradient of the log-likelihood) for a Weibull AFT regression model.
//' @param y A numeric vector of survival times.
//' @param dead A numeric vector of event indicators.
//' @param X A numeric matrix of predictors.
//' @param params A numeric vector of parameters [beta, log_sigma].
//' @return A numeric vector representing the score.
//' @export
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

//' @title Compute Weibull Regression Hessian
//' @description Calculates the Hessian matrix (second derivatives of the log-likelihood) for a Weibull AFT regression model.
//' @param y A numeric vector of survival times.
//' @param dead A numeric vector of event indicators.
//' @param X A numeric matrix of predictors.
//' @param params A numeric vector of parameters [beta, log_sigma].
//' @return A numeric matrix representing the Hessian.
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd get_weibull_regression_hessian_cpp(const Eigen::VectorXd& y,
                                                   const Eigen::VectorXd& dead,
                                                   const Eigen::MatrixXd& X,
                                                   const Eigen::VectorXd& params) {
    WeibullAFTLikelihood fun(y, dead, X);
    return -fun.hessian(params);
}

//' @title Fast Weibull Regression (C++)
//' @description High-performance Weibull Accelerated Failure Time (AFT) regression fitting.
//' @param y A numeric vector of survival times.
//' @param dead A numeric vector of event indicators.
//' @param X A numeric matrix of predictors.
//' @param start_params Optional starting values for [beta, log_sigma].
//' @param estimate_only If TRUE, only return coefficients and likelihood.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @return A list containing coefficients, log_sigma, vcov, and convergence status.
//' @export
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
                                 std::string optimization_alg = "newton_raphson") {
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
        fit = optimize_fixed_likelihood(fun, params, fixed_spec, maxit, tol, optimization_alg, "newton_raphson");
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
