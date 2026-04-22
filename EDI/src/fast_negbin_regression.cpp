#include "_helper_functions.h"
#include <RcppEigen.h>
#include <Rmath.h>

using namespace Rcpp;

namespace {

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
                                  int maxit = 1000,
                                  double eps_g = 1e-5,
                                  Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                  Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                  std::string optimization_alg = "newton_raphson") {
    int p = X.cols();
    ModelResult res;
    Eigen::VectorXd params = Eigen::VectorXd::Zero(p + 1);
    if (optimization_alg == "newton_raphson" || optimization_alg == "newton" || optimization_alg == "nr") {
        double mean_y = y.cast<double>().mean();
        if (p > 0 && X.col(0).array().isApprox(Eigen::ArrayXd::Ones(X.rows()), 1e-12)) {
            params[0] = std::log(std::max(mean_y, 1e-8));
        }
        double var_y = (y.cast<double>().array() - mean_y).square().sum() /
            static_cast<double>(std::max(1, static_cast<int>(y.size()) - 1));
        double theta_start = (var_y > mean_y && mean_y > 0.0) ?
            std::max(0.1, mean_y * mean_y / (var_y - mean_y)) : 10.0;
        params[p] = std::log(theta_start);
    }
    FixedParamSpec fixed_spec = make_fixed_param_spec(p + 1, fixed_idx, fixed_values);

    NBLogLik fun(X, y);
    LikelihoodFitResult fit = optimize_fixed_likelihood(fun, params, fixed_spec, maxit, eps_g, optimization_alg, "newton_raphson");
    params = fit.params;

    res.b = params.head(p);
    res.dispersion = std::exp(params[p]); // theta
    res.XtWX = fun.hessian(params); // Hessian
    res.converged = fit.converged;
    res.sigma2_hat = -fit.value; // using sigma2_hat to store logLik temporarily
    return res;
}

} // namespace

// [[Rcpp::export]]
Eigen::VectorXd get_negbin_regression_score_cpp(const Eigen::MatrixXd& X,
                                                const Eigen::VectorXi& y,
                                                const Eigen::VectorXd& params) {
    NBLogLik fun(X, y);
    Eigen::VectorXd grad(params.size());
    fun(params, grad);
    return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_negbin_regression_hessian_cpp(const Eigen::MatrixXd& X,
                                                  const Eigen::VectorXi& y,
                                                  const Eigen::VectorXd& params) {
    NBLogLik fun(X, y);
    return -fun.hessian(params);
}

// [[Rcpp::export]]
List fast_neg_bin_with_var_cpp(Eigen::MatrixXd X,
                                Eigen::VectorXi y,
                                int maxit = 1000,
                                double eps_f = 1e-8,
                                double eps_g = 1e-5,
                                Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                std::string optimization_alg = "newton_raphson") {
    ModelResult res = fast_neg_bin_internal(X, y, maxit, eps_g, fixed_idx, fixed_values, optimization_alg);
    FixedParamSpec fixed_spec = make_fixed_param_spec(X.cols() + 1, fixed_idx, fixed_values);
    Eigen::MatrixXd H_free = subset_matrix(res.XtWX, fixed_spec.free_idx, fixed_spec.free_idx);
    Eigen::MatrixXd cov_free = H_free.inverse();
    Eigen::MatrixXd vcov = expand_free_covariance(X.cols() + 1, fixed_spec, cov_free, true);
    return List::create(
        Named("b") = res.b,
        Named("theta_hat") = res.dispersion,
        Named("logLik") = res.sigma2_hat,
        Named("converged") = res.converged,
        Named("hess_fisher_info_matrix") = res.XtWX,
        Named("vcov") = vcov
    );
}

// [[Rcpp::export]]
List fast_neg_bin_cpp(Eigen::MatrixXd X,
                        Eigen::VectorXi y,
                        int maxit = 1000,
                        double eps_f = 1e-8,
                        double eps_g = 1e-5,
                        Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                        Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                        std::string optimization_alg = "newton_raphson") {
    ModelResult res = fast_neg_bin_internal(X, y, maxit, eps_g, fixed_idx, fixed_values, optimization_alg);
    return List::create(
        Named("b") = res.b,
        Named("theta_hat") = res.dispersion,
        Named("logLik") = res.sigma2_hat,
        Named("converged") = res.converged
    );
}
