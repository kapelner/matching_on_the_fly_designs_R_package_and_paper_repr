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

        double* H_data = H.data();
        for (int i = 0; i < m_n; ++i) {
            double mu_i = std::exp(eta[i]);
            double yi = m_y[i];
            double denom = theta + mu_i;
            const double* xi = m_X.data() + i;  // xi[j * m_n] == X(i,j)

            double w_beta = mu_i * theta * (yi + theta) / (denom * denom);
            for (int c = 0; c < m_p; ++c) {
                const double wxi_c = w_beta * xi[c * m_n];
                for (int r = 0; r <= c; ++r)
                    H_data[r + c * total_p] += wxi_c * xi[r * m_n];
            }

            double d_score_beta_d_log_theta = theta * mu_i * (yi - mu_i) / (denom * denom);
            for (int r = 0; r < m_p; ++r)
                H_data[r + m_p * total_p] -= d_score_beta_d_log_theta * xi[r * m_n];

            double A = R::digamma(yi + theta) - R::digamma(theta) +
                std::log(theta) - std::log(denom) + 1.0 - (yi + theta) / denom;
            double dA_dtheta = R::trigamma(yi + theta) - R::trigamma(theta) +
                1.0 / theta - 1.0 / denom + (yi - mu_i) / (denom * denom);
            H(m_p, m_p) -= theta * A + theta * theta * dA_dtheta;
        }
        for (int c = 0; c < total_p; ++c)
            for (int r = 0; r < c; ++r)
                H_data[c + r * total_p] = H_data[r + c * total_p];
        return H;
    }

    Eigen::MatrixXd expected_hessian(const Eigen::VectorXd& params) {
        int total_p = m_p + 1;
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total_p, total_p);
        Eigen::VectorXd beta = params.head(m_p);
        double theta = std::exp(params[m_p]);
        Eigen::VectorXd eta = m_X * beta;

        double* H_data = H.data();
        for (int i = 0; i < m_n; ++i) {
            double mu_i = std::exp(eta[i]);
            double denom = theta + mu_i;
            const double* xi = m_X.data() + i;

            double w_beta = mu_i * theta / denom;
            for (int c = 0; c < m_p; ++c) {
                const double wxi_c = w_beta * xi[c * m_n];
                for (int r = 0; r <= c; ++r)
                    H_data[r + c * total_p] += wxi_c * xi[r * m_n];
            }

            const double e_trigamma = expected_trigamma_y_plus_theta(mu_i, theta);
            H(m_p, m_p) += -theta * theta * (
                e_trigamma - R::trigamma(theta) + 1.0 / theta - 1.0 / denom
            );
        }

        for (int c = 0; c < total_p; ++c)
            for (int r = 0; r < c; ++r)
                H_data[c + r * total_p] = H_data[r + c * total_p];
        return H;
    }

private:
    static double expected_trigamma_y_plus_theta(double mu, double theta) {
        const double prob_success = theta / (theta + mu);
        double pk = std::exp(theta * std::log(prob_success));
        double sum = pk * R::trigamma(theta);
        double cdf = pk;
        const double ratio_base = mu / (theta + mu);
        const double mean = mu;
        const double sd = std::sqrt(mu + mu * mu / theta);
        const int min_iter = static_cast<int>(std::ceil(mean + 10.0 * sd));
        const int max_iter = 100000;

        for (int k = 0; k < max_iter; ++k) {
            pk *= (static_cast<double>(k) + theta) / static_cast<double>(k + 1) * ratio_base;
            const int y = k + 1;
            sum += pk * R::trigamma(static_cast<double>(y) + theta);
            cdf += pk;
            if (y > min_iter && pk < 1e-14 && 1.0 - cdf < 1e-12) break;
        }
        return sum;
    }
};

ModelResult fast_neg_bin_internal(const Eigen::MatrixXd& X,
                                  const Eigen::VectorXi& y,
                                  Rcpp::Nullable<Rcpp::NumericVector> warm_start_params = R_NilValue,
                                  bool smart_cold_start = true,
                                  int maxit = 1000,
                                  double eps_g = 1e-6,
                                  Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                  Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                  std::string optimization_alg = "lbfgs",
                                  Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                  bool estimate_only = false) {
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
    if (warm_start_params.isNotNull()) {
        params = as<Eigen::VectorXd>(NumericVector(warm_start_params));
        if (params.size() != p + 1) stop("warm_start_params must have length equal to the number of model parameters");
    } else if (smart_cold_start) {
        Eigen::VectorXd beta_smart = ols_smart_cold_start_beta_on_log1p(X, y_double);
        Eigen::VectorXd beta_start = vector_is_usable_start(beta_smart, p) ? beta_smart : legacy_params.head(p);
        params.head(p) = beta_start;
        params[p] = std::log(smart_negbin_theta_start_from_beta(X, y, beta_start, theta_start));
    } else {
        params = legacy_params;
    }
    params = apply_fixed_values(params, fixed_spec);

    NBLogLik fun(X, y);
    
    Eigen::MatrixXd H_start_val;
    const Eigen::MatrixXd* h_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        H_start_val = as<Eigen::MatrixXd>(warm_start_fisher_info);
        h_ptr = &H_start_val;
    } else if (smart_cold_start) {
        H_start_val = edi_opt::negbin_smart_hessian(X, params, y);
        h_ptr = &H_start_val;
    }
    
    LikelihoodFitResult fit = optimize_fixed_likelihood(fun, params, fixed_spec, maxit, eps_g, optimization_alg, "lbfgs", 0, h_ptr);
    params = fit.params;

    res.b = params.head(p);
    res.dispersion = std::exp(params[p]); // theta
    res.XtWX = estimate_only ? Eigen::MatrixXd::Zero(p+1, p+1) : fun.expected_hessian(params);
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
//' @param warm_start_params Optional starting values for coefficients and dispersion. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
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
                                Nullable<NumericVector> warm_start_params = R_NilValue,
                                bool smart_cold_start = false,
                                int maxit = 1000,
                                double eps_f = 1e-8,
                                double eps_g = 1e-6,
                                Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                std::string optimization_alg = "lbfgs",
                                Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                bool estimate_only = false) {
    ModelResult res = fast_neg_bin_internal(X, y, warm_start_params, smart_cold_start, maxit, eps_g, fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info, estimate_only);
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
//' @param warm_start_params Optional starting values for coefficients and dispersion. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
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
                        Nullable<NumericVector> warm_start_params = R_NilValue,
                        bool smart_cold_start = false,
                        int maxit = 1000,
                        double eps_f = 1e-8,
                        double eps_g = 1e-6,
                        Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                        Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                        std::string optimization_alg = "lbfgs",
                        Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                        bool estimate_only = false) {
    ModelResult res = fast_neg_bin_internal(X, y, warm_start_params, smart_cold_start, maxit, eps_g, fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info, estimate_only);
    return List::create(
        Named("b") = res.b,
        Named("theta_hat") = res.dispersion,
        Named("logLik") = res.sigma2_hat,
        Named("converged") = res.converged,
        Named("iterations") = res.iterations,
        Named("fisher_information") = res.XtWX
    );
}
