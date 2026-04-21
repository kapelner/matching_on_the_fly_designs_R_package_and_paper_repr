#include "_helper_functions.h"
#include <RcppEigen.h>

using namespace Rcpp;

namespace {

inline double log1pexp_logistic(double x) {
    if (x > 0.0) return x + std::log1p(std::exp(-x));
    return std::log1p(std::exp(x));
}

class LogisticNegLogLik {
private:
    const Eigen::MatrixXd& m_X;
    const Eigen::VectorXd& m_y;
    const Eigen::VectorXd& m_weights;
    const bool m_use_weights;
    const int m_n;

public:
    LogisticNegLogLik(const Eigen::MatrixXd& X,
                      const Eigen::VectorXd& y,
                      const Eigen::VectorXd& weights) :
        m_X(X), m_y(y), m_weights(weights), m_use_weights(weights.size() == X.rows()), m_n(X.rows()) {}

    double operator()(const Eigen::VectorXd& beta, Eigen::VectorXd& grad) {
        Eigen::VectorXd eta = m_X * beta;
        Eigen::VectorXd mu = (1.0 / (1.0 + (-eta).array().exp())).matrix();
        Eigen::VectorXd wt = m_use_weights ? m_weights : Eigen::VectorXd::Ones(m_n);

        double neg_ll = 0.0;
        for (int i = 0; i < m_n; ++i) {
            neg_ll += wt[i] * (log1pexp_logistic(eta[i]) - m_y[i] * eta[i]);
        }
        grad = m_X.transpose() * wt.cwiseProduct(mu - m_y);
        return neg_ll;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& beta) {
        Eigen::VectorXd eta = m_X * beta;
        Eigen::VectorXd mu = (1.0 / (1.0 + (-eta).array().exp())).matrix();
        Eigen::VectorXd w = mu.array() * (1.0 - mu.array());
        if (m_use_weights) w = w.cwiseProduct(m_weights);
        return m_X.transpose() * w.asDiagonal() * m_X;
    }
};

} // namespace

// Internal pure C++ logic
ModelResult fast_logistic_regression_internal(const Eigen::MatrixXd& X, 
                                              const Eigen::VectorXd& y, 
                                              const Eigen::VectorXd& weights = Eigen::VectorXd(),
                                              int maxit = 100, 
                                              double tol = 1e-8,
                                              Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                              Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                              std::string optimization_alg = "irls") {
    int n = X.rows();
    int p = X.cols();
    bool use_weights = (weights.size() == n);
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);

    std::string alg = normalize_optimizer_algorithm(optimization_alg, "irls", true);
    if (alg != "irls") {
        Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
        beta = apply_fixed_values(beta, fixed_spec);
        LogisticNegLogLik fun(X, y, weights);
        LikelihoodFitResult fit = (alg == "lbfgs") ?
            optimize_fixed_likelihood_lbfgs(fun, beta, fixed_spec, maxit, tol) :
            optimize_fixed_likelihood_newton(fun, beta, fixed_spec, maxit, tol);

        ModelResult res;
        res.b = fit.params;
        Eigen::VectorXd eta = X * res.b;
        res.mu = (1.0 / (1.0 + (-eta).array().exp())).matrix();
        Eigen::MatrixXd info_full = fun.hessian(res.b);
        res.XtWX = info_full;
        res.converged = fit.converged;
        return res;
    }

    const int p_free = fixed_spec.free_idx.size();
    Eigen::MatrixXd X_free(n, p_free);
    for (int j = 0; j < p_free; ++j) X_free.col(j) = X.col(fixed_spec.free_idx[j]);
    Eigen::VectorXd beta_free = Eigen::VectorXd::Zero(p_free);
    Eigen::VectorXd eta_fixed = Eigen::VectorXd::Zero(n);
    for (int j = 0; j < fixed_spec.fixed_idx.size(); ++j) {
        eta_fixed.noalias() += X.col(fixed_spec.fixed_idx[j]) * fixed_spec.fixed_values[j];
    }
    
    ModelResult res;
    res.b = Eigen::VectorXd::Zero(p);
    res.mu = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd w(n);
    Eigen::VectorXd z(n);

    for (int iter = 0; iter < maxit; iter++) {
        Eigen::VectorXd eta = eta_fixed + X_free * beta_free;
        for(int i=0; i<n; i++) {
            double val = 1.0 / (1.0 + std::exp(-eta[i]));
            res.mu[i] = val;
            double base_w = std::max(val * (1.0 - val), 1e-10);
            w[i] = use_weights ? weights[i] * base_w : base_w;
        }

        z = eta + (y - res.mu).cwiseQuotient((res.mu.array() * (1.0 - res.mu.array())).matrix());
        // Note: The adjusted dependent variable z in IRLS for binomial is eta + (y - mu) / [mu(1-mu)]
        // The weights in W are weights_i * mu_i * (1-mu_i)
        
        Eigen::VectorXd z_adj = z - eta_fixed;
        Eigen::MatrixXd XtW = X_free.transpose() * w.asDiagonal();
        Eigen::MatrixXd XtWX_free = XtW * X_free;
        Eigen::VectorXd XtWz = XtW * z_adj;

        Eigen::VectorXd beta_new = XtWX_free.ldlt().solve(XtWz);

        if ((beta_free - beta_new).norm() < tol) {
            beta_free = beta_new;
            res.converged = true;
            break;
        }
        beta_free = beta_new;
    }

    for (int j = 0; j < p_free; ++j) res.b[fixed_spec.free_idx[j]] = beta_free[j];
    for (int j = 0; j < fixed_spec.fixed_idx.size(); ++j) res.b[fixed_spec.fixed_idx[j]] = fixed_spec.fixed_values[j];

    // Final update of mu and XtWX
    Eigen::VectorXd eta = X * res.b;
    for(int i=0; i<n; i++) {
        double val = 1.0 / (1.0 + std::exp(-eta[i]));
        res.mu[i] = val;
        double base_w = std::max(val * (1.0 - val), 1e-10);
        w[i] = use_weights ? weights[i] * base_w : base_w;
    }
    Eigen::MatrixXd info_free = X_free.transpose() * w.asDiagonal() * X_free;
    res.XtWX = expand_free_covariance(p, fixed_spec, info_free, false);

    return res;
}

// [[Rcpp::export]]
Eigen::VectorXd get_logistic_regression_score_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& beta) {
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd mu = (1.0 / (1.0 + (-eta).array().exp())).matrix();
    return X.transpose() * (y - mu);
}

// [[Rcpp::export]]
Eigen::MatrixXd get_logistic_regression_hessian_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) {
    int n = X.rows();
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd w(n);
    for(int i=0; i<n; i++) {
        double mu_i = 1.0 / (1.0 + std::exp(-eta[i]));
        w[i] = mu_i * (1.0 - mu_i);
    }
    return - (X.transpose() * w.asDiagonal() * X);
}

// [[Rcpp::export]]
Eigen::VectorXd get_logistic_regression_weighted_score_cpp(const Eigen::MatrixXd& X,
                                                           const Eigen::VectorXd& y,
                                                           const Eigen::VectorXd& weights,
                                                           const Eigen::VectorXd& beta) {
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd mu = (1.0 / (1.0 + (-eta).array().exp())).matrix();
    return X.transpose() * weights.cwiseProduct(y - mu);
}

// [[Rcpp::export]]
Eigen::MatrixXd get_logistic_regression_weighted_hessian_cpp(const Eigen::MatrixXd& X,
                                                             const Eigen::VectorXd& weights,
                                                             const Eigen::VectorXd& beta) {
    int n = X.rows();
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd w(n);
    for(int i=0; i<n; i++) {
        double mu_i = 1.0 / (1.0 + std::exp(-eta[i]));
        w[i] = weights[i] * mu_i * (1.0 - mu_i);
    }
    return - (X.transpose() * w.asDiagonal() * X);
}

// [[Rcpp::export]]
List fast_logistic_regression_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int maxit = 100, double tol = 1e-8,
                                  Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                  Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                  std::string optimization_alg = "irls") {
    ModelResult res = fast_logistic_regression_internal(X, y, Eigen::VectorXd(), maxit, tol, fixed_idx, fixed_values, optimization_alg);
    return List::create(
        Named("b") = res.b,
        Named("w") = res.mu.array() * (1.0 - res.mu.array()) // return weights as expected by old code
    );
}

// [[Rcpp::export]]
List fast_logistic_regression_weighted_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& weights, int maxit = 100, double tol = 1e-8,
                                           Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                           Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                           std::string optimization_alg = "irls") {
    ModelResult res = fast_logistic_regression_internal(X, y, weights, maxit, tol, fixed_idx, fixed_values, optimization_alg);
    return List::create(
        Named("b") = res.b,
        Named("mu") = res.mu,
        Named("XtWX") = res.XtWX,
        Named("converged") = res.converged
    );
}

// [[Rcpp::export]]
List fast_logistic_regression_with_var_cpp(const Eigen::MatrixXd& Xmm, const Eigen::VectorXd& y, int j = 2,
                                           Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                           Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                           std::string optimization_alg = "irls") {
    ModelResult res = fast_logistic_regression_internal(Xmm, y, Eigen::VectorXd(), 100, 1e-8, fixed_idx, fixed_values, optimization_alg);
    FixedParamSpec fixed_spec = make_fixed_param_spec(Xmm.cols(), fixed_idx, fixed_values);
    Eigen::MatrixXd info_free = subset_matrix(res.XtWX, fixed_spec.free_idx, fixed_spec.free_idx);
    Eigen::MatrixXd cov_free = info_free.inverse();
    Eigen::MatrixXd vcov = expand_free_covariance(Xmm.cols(), fixed_spec, cov_free, true);
    res.ssq_b_j = (j > 0 && j <= Xmm.cols()) ? vcov(j - 1, j - 1) : NA_REAL;
    res.ssq_b_2 = (Xmm.cols() >= 2) ? vcov(1, 1) : NA_REAL;

    return List::create(
        Named("b") = res.b,
        Named("ssq_b_j") = res.ssq_b_j,
        Named("ssq_b_2") = res.ssq_b_2,
        Named("vcov") = vcov
    );
}
