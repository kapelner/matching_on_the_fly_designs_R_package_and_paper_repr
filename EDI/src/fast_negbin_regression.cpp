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

ModelResult fast_neg_bin_internal(const Eigen::MatrixXd& X,
                                  const Eigen::VectorXi& y,
                                  int maxit = 1000,
                                  double eps_g = 1e-5,
                                  Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                  Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                  std::string optimization_alg = "lbfgs") {
    int p = X.cols();
    ModelResult res;
    Eigen::VectorXd params = Eigen::VectorXd::Zero(p + 1);
    params[p] = 0.0; 
    FixedParamSpec fixed_spec = make_fixed_param_spec(p + 1, fixed_idx, fixed_values);

    NBLogLik fun(X, y);
    LikelihoodFitResult fit = optimize_fixed_likelihood(fun, params, fixed_spec, maxit, eps_g, optimization_alg, "lbfgs");
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
                                std::string optimization_alg = "lbfgs") {
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
                        std::string optimization_alg = "lbfgs") {
    ModelResult res = fast_neg_bin_internal(X, y, maxit, eps_g, fixed_idx, fixed_values, optimization_alg);
    return List::create(
        Named("b") = res.b,
        Named("theta_hat") = res.dispersion,
        Named("logLik") = res.sigma2_hat,
        Named("converged") = res.converged
    );
}
