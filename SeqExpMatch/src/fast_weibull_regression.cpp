#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

class WeibullRegression {
private:
    const Eigen::VectorXd m_y;
    const Eigen::VectorXi m_dead;
    const Eigen::MatrixXd m_X;
    const int m_n;
    const int m_p;

public:
    WeibullRegression(const Eigen::VectorXd& y, const Eigen::VectorXi& dead, const Eigen::MatrixXd& X) : 
        m_y(y), m_dead(dead), m_X(X), m_n(X.rows()), m_p(X.cols()) {}

    void gradient(const Eigen::VectorXd& params, Eigen::VectorXd& grad) {
        Eigen::VectorXd beta = params.head(m_p);
        double log_sigma = params(m_p);
        double sigma = std::exp(log_sigma);

        Eigen::VectorXd eta = m_X * beta;
        Eigen::VectorXd z = (m_y.array().log().matrix() - eta) / sigma;
        Eigen::VectorXd exp_z = z.array().exp().matrix();
        
        Eigen::VectorXd d_neg_ll_d_z = exp_z - m_dead.cast<double>();

        grad.head(m_p) = (-1.0 / sigma) * m_X.transpose() * d_neg_ll_d_z;
        grad(m_p) = - (z.array() * d_neg_ll_d_z.array()).sum() + m_dead.sum();
    }

    double operator()(const Eigen::VectorXd& params, Eigen::VectorXd& grad) {
        // params contains beta (p) and log(sigma) (1)
        Eigen::VectorXd beta = params.head(m_p);
        double log_sigma = params[m_p];
        double sigma = std::exp(log_sigma);

        Eigen::VectorXd eta = m_X * beta;
        Eigen::VectorXd z = (m_y.array().log().matrix() - eta) / sigma;
        Eigen::VectorXd exp_z = z.array().exp().matrix();

        // Negative Log-likelihood: sum(exp(z) - dead*z + dead*log_sigma)
        double neg_ll = (exp_z.array() - m_dead.cast<double>().array() * z.array()).sum();
        neg_ll += m_dead.sum() * log_sigma;

        // Gradient
        gradient(params, grad);
        
        return neg_ll;
    }
    
    void hessian(const Eigen::VectorXd& params, Eigen::MatrixXd& H) {
        Eigen::VectorXd beta = params.head(m_p);
        double log_sigma = params(m_p);
        double sigma = std::exp(log_sigma);

        Eigen::VectorXd eta = m_X * beta;
        Eigen::VectorXd z = (m_y.array().log().matrix() - eta) / sigma;
        Eigen::VectorXd exp_z = z.array().exp().matrix();

        H.setZero();
        
        Eigen::MatrixXd W = (exp_z.array() / (sigma * sigma)).matrix().asDiagonal();
        H.topLeftCorner(m_p, m_p) = m_X.transpose() * W * m_X;

        Eigen::VectorXd h_bs = (exp_z.array() * z.array() + exp_z.array() - m_dead.cast<double>().array()).matrix();
        H.topRightCorner(m_p, 1) = (1.0 / sigma) * m_X.transpose() * h_bs;
        H.bottomLeftCorner(1, m_p) = H.topRightCorner(m_p, 1).transpose();

        H(m_p, m_p) = (exp_z.array() * z.array().square() + exp_z.array() * z.array() - m_dead.cast<double>().array() * z.array()).sum();
    }
};

//' Fast Weibull Regression using Rcpp and Newton-Raphson
//'
//' @param y Numeric vector of survival times.
//' @param dead Integer vector of event indicators (1=event, 0=censored).
//' @param X Model matrix.
//' @return A list with coefficients, log(sigma), standard errors, vcov, etc.
//' @export
// [[Rcpp::export]]
List fast_weibull_regression(NumericVector y, IntegerVector dead, Eigen::MatrixXd X) {
    // We need to add an intercept term to X
    Eigen::MatrixXd X_intercept(X.rows(), X.cols() + 1);
    X_intercept.col(0).setOnes();
    X_intercept.rightCols(X.cols()) = X;

    Eigen::VectorXd y_eigen = Rcpp::as<Eigen::VectorXd>(y);
    Eigen::VectorXi dead_eigen = Rcpp::as<Eigen::VectorXi>(dead);

    int p = X_intercept.cols();
    
    // Initial parameters
    Eigen::VectorXd params = Eigen::VectorXd::Zero(p + 1);
    params[p] = log(1.0); // log(sigma) = log(1.0) => sigma = 1

    WeibullRegression fun(y_eigen, dead_eigen, X_intercept);

    // Newton-Raphson optimization
    int max_iter = 100;
    double tol = 1e-6;
    int niter = 0;
    for (int i = 0; i < max_iter; ++i) {
        Eigen::VectorXd grad(p + 1);
        Eigen::MatrixXd H_optim(p + 1, p + 1);

        fun(params, grad);
        fun.hessian(params, H_optim);
        
        Eigen::VectorXd delta = (H_optim + 1e-6 * Eigen::MatrixXd::Identity(p + 1, p + 1)).colPivHouseholderQr().solve(grad);
        params -= 0.5 * delta; // Apply damping factor
        niter++;

        if (delta.norm() < tol) {
            break;
        }
    }
    
    Eigen::VectorXd beta = params.head(p);
    double log_sigma = params[p];

    // Final Hessian and standard errors
    Eigen::MatrixXd H_final(p + 1, p + 1);
    fun.hessian(params, H_final);
    Eigen::MatrixXd cov_mat = (H_final + 1e-6 * Eigen::MatrixXd::Identity(p + 1, p + 1)).inverse();    
    Eigen::VectorXd std_errs = cov_mat.diagonal().array().sqrt();
    
    // Calculate final neg_ll
    Eigen::VectorXd grad_final(p + 1);
    double neg_ll = fun(params, grad_final);

    return List::create(
        Named("coefficients") = beta,
        Named("log_sigma") = log_sigma,
        Named("std_errs") = std_errs,
        Named("vcov") = cov_mat,
        Named("n_iter") = niter,
        Named("neg_log_lik") = neg_ll
    );
}

//' Fast Weibull Regression with variance using Rcpp and Newton-Raphson
//'
//' @param y Numeric vector of survival times.
//' @param dead Integer vector of event indicators (1=event, 0=censored).
//' @param X Model matrix.
//' @return A list with coefficients, log(sigma), standard errors, vcov, etc.
//' @export
// [[Rcpp::export]]
List fast_weibull_regression_with_var_cpp(NumericVector y, IntegerVector dead, Eigen::MatrixXd X) {
    return fast_weibull_regression(y, dead, X);
}