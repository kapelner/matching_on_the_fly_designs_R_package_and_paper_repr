#include "_helper_functions.h"
#include "result_map_rcpp.h"
#include <RcppEigen.h>
#include <Rmath.h>
#include <unordered_map>

using namespace Rcpp;

namespace {

double smart_negbin_theta_start_from_beta(const Eigen::Ref<const Eigen::MatrixXd>& X,
                                          const Eigen::Ref<const Eigen::VectorXi>& y,
                                          const Eigen::Ref<const Eigen::VectorXd>& beta,
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
    const Eigen::Ref<const Eigen::MatrixXd> m_X;
    const Eigen::Ref<const Eigen::VectorXi> m_y;
    const Eigen::VectorXd m_weights;
    const int m_n;
    const int m_p;

    // Per-observation slot into distinct-y tables (built once at construction).
    std::vector<int>    m_y_slot;
    std::vector<double> m_distinct_y;
    std::vector<double> m_lgamma_y1;        // lgamma(y+1) = log(y!) per distinct y
    std::vector<double> m_lgamma_yptheta;   // preallocated; filled per operator() call
    std::vector<double> m_digamma_yptheta;  // preallocated; filled per operator() and hessian() call
    std::vector<double> m_trigamma_yptheta; // preallocated; filled per hessian() call
    Eigen::VectorXd     m_coef_vec;         // per-obs gradient scalar; GEMV after obs loop

public:
    NBLogLik(const Eigen::Ref<const Eigen::MatrixXd>& X, const Eigen::Ref<const Eigen::VectorXi>& y) :
        NBLogLik(X, y, Eigen::VectorXd::Ones(X.rows())) {}

    NBLogLik(const Eigen::Ref<const Eigen::MatrixXd>& X, const Eigen::Ref<const Eigen::VectorXi>& y,
             const Eigen::Ref<const Eigen::VectorXd>& weights) :
        m_X(X), m_y(y), m_weights(weights), m_n(X.rows()), m_p(X.cols()), m_y_slot(X.rows(), -1),
        m_coef_vec(X.rows())
    {
        std::unordered_map<int, int> seen;
        for (int i = 0; i < m_n; ++i) {
            const int yi = m_y[i];
            auto it = seen.find(yi);
            if (it == seen.end()) {
                const int slot = static_cast<int>(m_distinct_y.size());
                seen[yi] = slot;
                m_distinct_y.push_back(static_cast<double>(yi));
                m_lgamma_y1.push_back(std::lgamma(static_cast<double>(yi) + 1.0));
                m_y_slot[i] = slot;
            } else {
                m_y_slot[i] = it->second;
            }
        }
        const int nd = static_cast<int>(m_distinct_y.size());
        m_lgamma_yptheta.resize(nd);
        m_digamma_yptheta.resize(nd);
        m_trigamma_yptheta.resize(nd);
    }

    double operator()(const Eigen::VectorXd& params, Eigen::VectorXd& grad) {
        const Eigen::VectorXd beta = params.head(m_p);
        const double log_theta = params[m_p];
        const double theta = std::exp(log_theta);

        const Eigen::VectorXd eta = m_X * beta;

        // Hoist theta-only constants out of the observation loop.
        const double digamma_theta = fast_digamma(theta);
        const double lgamma_theta  = std::lgamma(theta);
        const double log_theta_val = std::log(theta);

        // Fill preallocated per-distinct-y tables for lgamma(y+theta) and digamma(y+theta).
        const int nd = static_cast<int>(m_distinct_y.size());
        for (int k = 0; k < nd; ++k) {
            const double ypt = m_distinct_y[k] + theta;
            m_lgamma_yptheta[k]  = std::lgamma(ypt);
            m_digamma_yptheta[k] = fast_digamma(ypt);
        }

        double neg_ll = 0.0;
        double score_log_theta = 0.0;

        for (int i = 0; i < m_n; ++i) {
            const double eta_i  = eta[i];
            const double mu_i   = std::exp(eta_i);
            const double yi     = m_distinct_y[m_y_slot[i]];
            const double denom  = theta + mu_i;
            const double log_denom = std::log(denom);
            const int slot = m_y_slot[i];
            const double w_i = m_weights[i];

            // log dnbinom_mu via explicit formula — avoids R::dnbinom_mu overhead.
            neg_ll -= w_i * (m_lgamma_yptheta[slot] - lgamma_theta - m_lgamma_y1[slot]
                    + theta * (log_theta_val - log_denom)
                    + yi * (eta_i - log_denom));

            m_coef_vec[i] = w_i * (yi - mu_i * (yi + theta) / denom);

            const double dlogf = m_digamma_yptheta[slot] - digamma_theta
                               + log_theta_val - log_denom + 1.0 - (yi + theta) / denom;
            score_log_theta += w_i * theta * dlogf;
        }
        grad.head(m_p).noalias() = -(m_X.transpose() * m_coef_vec);
        grad[m_p] = -score_log_theta;
        return neg_ll;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& params) {
        int total_p = m_p + 1;
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total_p, total_p);
        Eigen::VectorXd beta = params.head(m_p);
        double theta = std::exp(params[m_p]);
        Eigen::VectorXd eta = m_X * beta;

        const double digamma_theta = fast_digamma(theta);
        const double trigamma_theta = fast_trigamma(theta);
        const double log_theta = std::log(theta);
        const int nd = static_cast<int>(m_distinct_y.size());
        Eigen::ArrayXd ypt_arr(nd);
        for (int k = 0; k < nd; ++k) {
            const double ypt = m_distinct_y[k] + theta;
            ypt_arr[k] = ypt;
            m_digamma_yptheta[k] = fast_digamma(ypt);
        }
        Eigen::Map<Eigen::ArrayXd>(m_trigamma_yptheta.data(), nd) = fast_trigamma_vec(ypt_arr);

        double* H_data = H.data();
        for (int i = 0; i < m_n; ++i) {
            double mu_i = std::exp(eta[i]);
            const int slot = m_y_slot[i];
            double yi = m_distinct_y[slot];
            double denom = theta + mu_i;
            const double* xi = m_X.data() + i;  // xi[j * m_n] == X(i,j)

            const double w_i = m_weights[i];
            double w_beta = w_i * mu_i * theta * (yi + theta) / (denom * denom);
            for (int c = 0; c < m_p; ++c) {
                const double wxi_c = w_beta * xi[c * m_n];
                for (int r = 0; r <= c; ++r)
                    H_data[r + c * total_p] += wxi_c * xi[r * m_n];
            }

            double d_score_beta_d_log_theta = w_i * theta * mu_i * (yi - mu_i) / (denom * denom);
            for (int r = 0; r < m_p; ++r)
                H_data[r + m_p * total_p] -= d_score_beta_d_log_theta * xi[r * m_n];

            double A = m_digamma_yptheta[slot] - digamma_theta +
                log_theta - std::log(denom) + 1.0 - (yi + theta) / denom;
            double dA_dtheta = m_trigamma_yptheta[slot] - trigamma_theta +
                1.0 / theta - 1.0 / denom + (yi - mu_i) / (denom * denom);
            H(m_p, m_p) -= w_i * (theta * A + theta * theta * dA_dtheta);
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
        const double trigamma_theta = fast_trigamma(theta);

        double* H_data = H.data();
        for (int i = 0; i < m_n; ++i) {
            double mu_i = std::exp(eta[i]);
            double denom = theta + mu_i;
            const double* xi = m_X.data() + i;

            const double w_i = m_weights[i];
            double w_beta = w_i * mu_i * theta / denom;
            for (int c = 0; c < m_p; ++c) {
                const double wxi_c = w_beta * xi[c * m_n];
                for (int r = 0; r <= c; ++r)
                    H_data[r + c * total_p] += wxi_c * xi[r * m_n];
            }

            const double e_trigamma = expected_trigamma_y_plus_theta(mu_i, theta, trigamma_theta);
            H(m_p, m_p) += -w_i * theta * theta * (
                e_trigamma - trigamma_theta + 1.0 / theta - 1.0 / denom
            );
        }

        for (int c = 0; c < total_p; ++c)
            for (int r = 0; r < c; ++r)
                H_data[c + r * total_p] = H_data[r + c * total_p];
        return H;
    }

private:
    static double expected_trigamma_y_plus_theta(double mu,
                                                  double theta,
                                                  double trigamma_theta) {
        const double prob_success = theta / (theta + mu);
        double pk = std::exp(theta * std::log(prob_success));
        double trigamma_yptheta = trigamma_theta;
        double sum = pk * trigamma_yptheta;
        double cdf = pk;
        const double ratio_base = mu / (theta + mu);
        const double mean = mu;
        const double sd = std::sqrt(mu + mu * mu / theta);
        const int min_iter = static_cast<int>(std::ceil(mean + 10.0 * sd));
        const int max_iter = 100000;

        for (int k = 0; k < max_iter; ++k) {
            pk *= (static_cast<double>(k) + theta) / static_cast<double>(k + 1) * ratio_base;
            const int y = k + 1;
            const double x = static_cast<double>(k) + theta;
            trigamma_yptheta -= 1.0 / (x * x);
            sum += pk * trigamma_yptheta;
            cdf += pk;
            if (y > min_iter && pk < 1e-14 && 1.0 - cdf < 1e-12) break;
        }
        return sum;
    }
};

ModelResult fast_neg_bin_internal(const Eigen::Ref<const Eigen::MatrixXd>& X,
                                  const Eigen::Ref<const Eigen::VectorXi>& y,
                                  std::optional<Eigen::VectorXd> warm_start_params = std::nullopt,
                                  bool smart_cold_start = true,
                                  int maxit = 1000,
                                  double eps_g = 1e-6,
                                  std::optional<Eigen::VectorXi> fixed_idx = std::nullopt,
                                  std::optional<Eigen::VectorXd> fixed_values = std::nullopt,
                                  std::string optimization_alg = "lbfgs",
                                  std::optional<Eigen::MatrixXd> warm_start_fisher_info = std::nullopt,
                                  bool estimate_only = false,
                                  const Eigen::VectorXd* weights = nullptr) {
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
    if (warm_start_params.has_value()) {
        params = *warm_start_params;
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

    Eigen::VectorXd weights_work = weights == nullptr ? Eigen::VectorXd::Ones(X.rows()) : *weights;
    NBLogLik fun(X, y, weights_work);

    Eigen::MatrixXd H_start_val;
    const Eigen::MatrixXd* h_ptr = nullptr;
    if (warm_start_fisher_info.has_value()) {
        H_start_val = *warm_start_fisher_info;
        h_ptr = &H_start_val;
    } else if (smart_cold_start) {
        H_start_val = edi_opt::negbin_smart_hessian(X, params, y);
        h_ptr = &H_start_val;
    }

    LikelihoodFitResult fit = optimize_fixed_likelihood(fun, params, fixed_spec, maxit, eps_g, optimization_alg, "lbfgs", 0, h_ptr);
    params = fit.params;

    res.b = params.head(p);
    res.dispersion = std::exp(params[p]); // theta
    res.XtWX = estimate_only ? Eigen::MatrixXd::Zero(p+1, p+1) : fun.hessian(params);
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
Eigen::VectorXd get_negbin_regression_score_cpp(SEXP X_sexp,
                                                SEXP y_sexp,
                                                SEXP params_sexp) {
    NumericMatrix X_r(X_sexp);
    IntegerVector y_r(y_sexp);
    NumericVector params_r(params_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXi> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());

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
Eigen::MatrixXd get_negbin_regression_hessian_cpp(SEXP X_sexp,
                                                  SEXP y_sexp,
                                                  SEXP params_sexp) {
    NumericMatrix X_r(X_sexp);
    IntegerVector y_r(y_sexp);
    NumericVector params_r(params_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXi> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());

    NBLogLik fun(X, y);
    return -fun.hessian(params);
}

// [[Rcpp::export]]
Eigen::MatrixXd get_negbin_regression_expected_hessian_cpp(SEXP X_sexp,
                                                           SEXP y_sexp,
                                                           SEXP params_sexp) {
    NumericMatrix X_r(X_sexp);
    IntegerVector y_r(y_sexp);
    NumericVector params_r(params_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXi> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());

    NBLogLik fun(X, y);
    return fun.expected_hessian(params);
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
List fast_neg_bin_with_var_cpp(SEXP X_sexp,
                                SEXP y_sexp,
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
    NumericMatrix X_r(X_sexp);
    IntegerVector y_r(y_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXi> y(y_r.begin(), y_r.size());

    ModelResult res = fast_neg_bin_internal(
        X, y,
        nullable_to_optional<Eigen::VectorXd>(warm_start_params),
        smart_cold_start, maxit, eps_g,
        nullable_to_optional<Eigen::VectorXi>(fixed_idx),
        nullable_to_optional<Eigen::VectorXd>(fixed_values),
        optimization_alg,
        nullable_to_optional<Eigen::MatrixXd>(warm_start_fisher_info),
        estimate_only);
    FixedParamSpec fixed_spec = make_fixed_param_spec(
        X.cols() + 1,
        nullable_to_optional<Eigen::VectorXi>(fixed_idx),
        nullable_to_optional<Eigen::VectorXd>(fixed_values));
    Eigen::MatrixXd H_free = subset_matrix(res.XtWX, fixed_spec.free_idx, fixed_spec.free_idx);
    Eigen::MatrixXd cov_free = H_free.inverse();
    Eigen::MatrixXd vcov = expand_free_covariance(X.cols() + 1, fixed_spec, cov_free, true);
    return edi::to_rcpp_list(edi::ResultMap()
        .set("b", res.b)
        .set("theta_hat", res.dispersion)
        .set("logLik", res.sigma2_hat)
        .set("converged", res.converged)
        .set("iterations", res.iterations)
        .set("hess_fisher_info_matrix", res.XtWX)
        .set("vcov", vcov));
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
List fast_neg_bin_cpp(SEXP X_sexp,
                        SEXP y_sexp,
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
    NumericMatrix X_r(X_sexp);
    IntegerVector y_r(y_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXi> y(y_r.begin(), y_r.size());

    ModelResult res = fast_neg_bin_internal(
        X, y,
        nullable_to_optional<Eigen::VectorXd>(warm_start_params),
        smart_cold_start, maxit, eps_g,
        nullable_to_optional<Eigen::VectorXi>(fixed_idx),
        nullable_to_optional<Eigen::VectorXd>(fixed_values),
        optimization_alg,
        nullable_to_optional<Eigen::MatrixXd>(warm_start_fisher_info),
        estimate_only);
    return edi::to_rcpp_list(edi::ResultMap()
        .set("b", res.b)
        .set("theta_hat", res.dispersion)
        .set("logLik", res.sigma2_hat)
        .set("converged", res.converged)
        .set("iterations", res.iterations)
        .set("fisher_information", res.XtWX));
}

//' @title Fast Weighted Negative Binomial Regression (C++)
//' @description High-performance negative binomial regression fitting with nonnegative row weights.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param weights A nonnegative numeric vector of row weights.
//' @param warm_start_params Optional starting values for coefficients and dispersion. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param maxit Maximum number of iterations.
//' @param eps_f Convergence tolerance for function value.
//' @param eps_g Convergence tolerance for gradient.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @param estimate_only If TRUE, skip Fisher information calculation.
//' @return A list containing coefficients, theta, and convergence status.
//' @export
//' @keywords internal
//' @examples
//' X = matrix(rnorm(100), 10, 10)
//' y = rpois(10, 2)
//' fast_neg_bin_weighted_cpp(X, y, weights = rep(1, 10))
// [[Rcpp::export]]
List fast_neg_bin_weighted_cpp(SEXP X_sexp,
                        SEXP y_sexp,
                        SEXP weights_sexp,
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
    NumericMatrix X_r(X_sexp);
    IntegerVector y_r(y_sexp);
    NumericVector weights_r(weights_sexp);
    if (weights_r.size() != X_r.nrow()) {
        stop("weights length must equal nrow(X)");
    }
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXi> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::VectorXd> weights(weights_r.begin(), weights_r.size());
    if ((weights.array() < 0.0).any() || !weights.allFinite() || weights.sum() <= 0.0) {
        stop("weights must be finite, nonnegative, and have positive sum");
    }
    Eigen::VectorXd weights_vec = weights;

    ModelResult res = fast_neg_bin_internal(
        X, y,
        nullable_to_optional<Eigen::VectorXd>(warm_start_params),
        smart_cold_start, maxit, eps_g,
        nullable_to_optional<Eigen::VectorXi>(fixed_idx),
        nullable_to_optional<Eigen::VectorXd>(fixed_values),
        optimization_alg,
        nullable_to_optional<Eigen::MatrixXd>(warm_start_fisher_info),
        estimate_only, &weights_vec);
    return edi::to_rcpp_list(edi::ResultMap()
        .set("b", res.b)
        .set("theta_hat", res.dispersion)
        .set("logLik", res.sigma2_hat)
        .set("converged", res.converged)
        .set("iterations", res.iterations)
        .set("fisher_information", res.XtWX));
}
