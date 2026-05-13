#include "_helper_functions.h"
#include <RcppEigen.h>
#include <Rmath.h>

using namespace Rcpp;

namespace {

double smart_negbin_theta_start_from_beta(const Eigen::MatrixXd& X,
                                          const Eigen::VectorXi& y,
                                          const Eigen::VectorXd& beta,
                                          double legacy_theta_start) {
    const int n = X.rows();
    const int p = X.cols();
    const int df = std::max(1, n - p);
    if (beta.size() != p || !beta.allFinite()) return legacy_theta_start;

    Eigen::VectorXd eta = (X * beta).array().min(700.0).matrix();
    Eigen::VectorXd mu = eta.array().exp().max(1e-8).matrix();
    double alpha_sum = 0.0;

    for (int i = 0; i < n; ++i) {
        const double yi = static_cast<double>(y[i]);
        const double mui = mu[i];
        alpha_sum += ((yi - mui) * (yi - mui) - mui) / (mui * mui);
    }

    const double alpha_hat = alpha_sum / static_cast<double>(df);
    if (!std::isfinite(alpha_hat) || alpha_hat <= 0.0) return legacy_theta_start;

    const double theta_hat = 1.0 / alpha_hat;
    if (!std::isfinite(theta_hat) || theta_hat <= 0.0) return legacy_theta_start;
    return std::max(0.1, theta_hat);
}

class NBLogLik {
private:
    const Eigen::MatrixXd m_X;
    const Eigen::VectorXi m_y;
    const int m_n;
    const int m_p;

public:
    NBLogLik(const Eigen::MatrixXd& X, const Eigen::VectorXi& y) :
        m_X(X), m_y(y), m_n(X.rows()), m_p(X.cols()) {}

    double operator()(const Eigen::VectorXd& params, Eigen::VectorXd& grad) {
        Eigen::VectorXd beta = params.head(m_p);
        double log_theta = params[m_p];
        double theta = std::exp(log_theta);

        Eigen::VectorXd eta = m_X * beta;
        Eigen::VectorXd mu = eta.array().exp();

        double neg_ll = 0.0;
        Eigen::VectorXd score_beta = Eigen::VectorXd::Zero(m_p);
        double score_log_theta = 0.0;

        for (int i = 0; i < m_n; ++i) {
            double mu_i = mu[i];
            double yi = m_y[i];
            neg_ll -= R::dnbinom_mu(yi, theta, mu_i, true);
            double coef = yi - mu_i * (yi + theta) / (theta + mu_i);
            score_beta += coef * m_X.row(i).transpose();
            double dlogf_dtheta = R::digamma(yi + theta) - R::digamma(theta) + std::log(theta) - std::log(theta + mu_i) + 1.0 - (yi + theta) / (theta + mu_i);
            score_log_theta += theta * dlogf_dtheta;
        }
        grad.head(m_p) = -score_beta;
        grad[m_p] = -score_log_theta;
        return neg_ll;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& params) {
        int total_p = m_p + 1;
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total_p, total_p);
        Eigen::VectorXd beta = params.head(m_p);
        double theta = std::exp(params[m_p]);
        Eigen::VectorXd eta = m_X * beta;

        for (int i = 0; i < m_n; ++i) {
            double mu_i = std::exp(eta[i]);
            double yi = m_y[i];
            double denom = theta + mu_i;
            Eigen::VectorXd x = m_X.row(i).transpose();

            double w_beta = mu_i * theta * (yi + theta) / (denom * denom);
            H.topLeftCorner(m_p, m_p).noalias() += w_beta * (x * x.transpose());

            double d_score_beta_d_log_theta = theta * mu_i * (yi - mu_i) / (denom * denom);
            H.topRightCorner(m_p, 1).noalias() -= d_score_beta_d_log_theta * x;

            double A = R::digamma(yi + theta) - R::digamma(theta) +
                std::log(theta) - std::log(denom) + 1.0 - (yi + theta) / denom;
            double dA_dtheta = R::trigamma(yi + theta) - R::trigamma(theta) +
                1.0 / theta - 1.0 / denom + (yi - mu_i) / (denom * denom);
            H(m_p, m_p) -= theta * A + theta * theta * dA_dtheta;
        }
        H.bottomLeftCorner(1, m_p) = H.topRightCorner(m_p, 1).transpose();
        return H;
    }
};

ModelResult fast_neg_bin_internal(const Eigen::MatrixXd& X,
                                  const Eigen::VectorXi& y,
                                  Rcpp::Nullable<Rcpp::NumericVector> start_params = R_NilValue,
                                  bool smart_start = true,
                                  int maxit = 1000,
                                  double eps_g = 1e-5,
                                  Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                  Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                  std::string optimization_alg = "newton_raphson",
                                  Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    int p = X.cols();
    ModelResult res;
    Eigen::VectorXd params = Eigen::VectorXd::Zero(p + 1);
    const Eigen::VectorXd y_double = y.cast<double>();
    const double mean_y = y_double.mean();
    const double var_y = (y_double.array() - mean_y).square().sum() /
        static_cast<double>(std::max(1, static_cast<int>(y.size()) - 1));
    const double theta_start = (var_y > mean_y && mean_y > 0.0) ?
        std::max(0.1, mean_y * mean_y / (var_y - mean_y)) : 10.0;
    Eigen::VectorXd legacy_params = Eigen::VectorXd::Zero(p + 1);
    if (p > 0 && X.col(0).array().isApprox(Eigen::ArrayXd::Ones(X.rows()), 1e-12)) {
        legacy_params[0] = std::log(std::max(mean_y, 1e-8));
    }
    legacy_params[p] = std::log(theta_start);
    FixedParamSpec fixed_spec = make_fixed_param_spec(p + 1, fixed_idx, fixed_values);
    if (start_params.isNotNull()) {
        params = as<Eigen::VectorXd>(NumericVector(start_params));
        if (params.size() != p + 1) stop("start_params must have length equal to the number of model parameters");
    } else if (smart_start) {
        Eigen::VectorXd beta_smart = ols_start_beta_on_log1p(X, y_double);
        Eigen::VectorXd beta_start = vector_is_usable_start(beta_smart, p) ? beta_smart : legacy_params.head(p);
        params.head(p) = beta_start;
        params[p] = std::log(smart_negbin_theta_start_from_beta(X, y, beta_start, theta_start));
    } else {
        params = legacy_params;
    }
    params = apply_fixed_values(params, fixed_spec);

    NBLogLik fun(X, y);
    
    Eigen::MatrixXd H_start;
    const Eigen::MatrixXd* h_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        H_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        h_ptr = &H_start;
    }
    
    LikelihoodFitResult fit = optimize_fixed_likelihood(fun, params, fixed_spec, maxit, eps_g, optimization_alg, "newton_raphson", 0, h_ptr);
    params = fit.params;

    res.b = params.head(p);
    res.dispersion = std::exp(params[p]); // theta
    res.XtWX = fun.hessian(params); // Hessian
    res.iterations = fit.niter;
    res.converged = fit.converged;
    res.sigma2_hat = -fit.value; // using sigma2_hat to store logLik temporarily
    return res;
}

} // namespace

//' @title Compute Negative Binomial Regression Score
//' @description Calculates the score vector (gradient of the log-likelihood) for a negative binomial regression model.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (non-negative integers).
//' @param params A numeric vector of parameters [beta, log_theta].
//' @return A numeric vector representing the score.
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::VectorXd get_negbin_regression_score_cpp(const Eigen::MatrixXd& X,
                                                const Eigen::VectorXi& y,
                                                const Eigen::VectorXd& params) {
    NBLogLik fun(X, y);
    Eigen::VectorXd grad(params.size());
    fun(params, grad);
    return -grad;
}

//' @title Compute Negative Binomial Regression Hessian
//' @description Calculates the Hessian matrix (second derivatives of the log-likelihood) for a negative binomial regression model.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param params A numeric vector of parameters [beta, log_theta].
//' @return A numeric matrix representing the Hessian.
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd get_negbin_regression_hessian_cpp(const Eigen::MatrixXd& X,
                                                  const Eigen::VectorXi& y,
                                                  const Eigen::VectorXd& params) {
    NBLogLik fun(X, y);
    return -fun.hessian(params);
}

//' @title Fast Negative Binomial Regression with Variance (C++)
//' @description Negative binomial regression fitting with full variance-covariance matrix.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (non-negative integers).
//' @param start_params Optional starting values for coefficients and dispersion.
//' @param smart_start Logical. If TRUE, use an initial OLS-based guess.
//' @param maxit Maximum number of iterations.
//' @param eps_f Convergence tolerance for function value.
//' @param eps_g Convergence tolerance for gradient.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients, theta, vcov, and convergence status.
//' @export
//' @keywords internal
//' @examples
//' X = matrix(rnorm(100), 10, 10)
//' y = rpois(10, 2)
//' fast_neg_bin_with_var_cpp(X, y)
// [[Rcpp::export]]
List fast_neg_bin_with_var_cpp(Eigen::MatrixXd X,
                                Eigen::VectorXi y,
                                Nullable<NumericVector> start_params = R_NilValue,
                                bool smart_start = false,
                                int maxit = 1000,
                                double eps_f = 1e-8,
                                double eps_g = 1e-5,
                                Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                std::string optimization_alg = "newton_raphson",
                                Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    ModelResult res = fast_neg_bin_internal(X, y, start_params, smart_start, maxit, eps_g, fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info);
    FixedParamSpec fixed_spec = make_fixed_param_spec(X.cols() + 1, fixed_idx, fixed_values);
    Eigen::MatrixXd H_free = subset_matrix(res.XtWX, fixed_spec.free_idx, fixed_spec.free_idx);
    Eigen::MatrixXd cov_free = H_free.inverse();
    Eigen::MatrixXd vcov = expand_free_covariance(X.cols() + 1, fixed_spec, cov_free, true);
    return List::create(
        Named("b") = res.b,
        Named("theta_hat") = res.dispersion,
        Named("logLik") = res.sigma2_hat,
        Named("converged") = res.converged,
        Named("iterations") = res.iterations,
        Named("hess_fisher_info_matrix") = res.XtWX,
        Named("vcov") = vcov
    );
}

//' @title Fast Negative Binomial Regression (C++)
//' @description High-performance negative binomial regression fitting.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param start_params Optional starting values for coefficients and dispersion.
//' @param smart_start Logical. If TRUE, use an initial OLS-based guess.
//' @param maxit Maximum number of iterations.
//' @param eps_f Convergence tolerance for function value.
//' @param eps_g Convergence tolerance for gradient.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients, theta, and convergence status.
//' @export
//' @keywords internal
//' @examples
//' X = matrix(rnorm(100), 10, 10)
//' y = rpois(10, 2)
//' fast_neg_bin_cpp(X, y)
// [[Rcpp::export]]
List fast_neg_bin_cpp(Eigen::MatrixXd X,
                        Eigen::VectorXi y,
                        Nullable<NumericVector> start_params = R_NilValue,
                        bool smart_start = false,
                        int maxit = 1000,
                        double eps_f = 1e-8,
                        double eps_g = 1e-5,
                        Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                        Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                        std::string optimization_alg = "newton_raphson",
                        Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    ModelResult res = fast_neg_bin_internal(X, y, start_params, smart_start, maxit, eps_g, fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info);
    return List::create(
        Named("b") = res.b,
        Named("theta_hat") = res.dispersion,
        Named("logLik") = res.sigma2_hat,
        Named("converged") = res.converged,
        Named("iterations") = res.iterations
    );
}
