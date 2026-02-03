#include <Rcpp.h>
#include <RcppEigen.h>
#include <optimization/LBFGS.h>
#include <Rmath.h> // For R::digamma and R::trigamma

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace LBFGSpp;

// Functor for digamma
struct DigammaFunctor {
    double operator()(double x) const {
        return R::digamma(x);
    }
};

class BetaRegression {
private:
    const Eigen::VectorXd m_y;
    const Eigen::MatrixXd m_X;
    const int m_n;
    const int m_p;
    const Eigen::VectorXd m_log_y;
    const Eigen::VectorXd m_log1_y;

public:
    BetaRegression(const Eigen::VectorXd& y, const Eigen::MatrixXd& X) : 
        m_y(y), m_X(X), m_n(X.rows()), m_p(X.cols()),
        m_log_y(y.array().log().matrix()),
        m_log1_y((1.0 - y.array()).log().matrix()) {}

    double operator()(const Eigen::VectorXd& params, Eigen::VectorXd& grad) {
        Eigen::VectorXd beta = params.head(m_p);
        double log_phi = params[m_p];
        double phi = std::exp(log_phi);

        Eigen::VectorXd eta = m_X * beta;
        Eigen::VectorXd mu = (1.0 / (1.0 + (-eta).array().exp())).matrix();
        
        // Clamp mu to avoid numerical instability
        double epsilon = 1e-8;
        for(int i=0; i<m_n; ++i) {
            if (mu[i] < epsilon) mu[i] = epsilon;
            if (mu[i] > 1.0 - epsilon) mu[i] = 1.0 - epsilon;
        }

        // Negative log-likelihood
        Eigen::VectorXd mu_phi = mu.array() * phi;
        Eigen::VectorXd one_minus_mu_phi = (1.0 - mu.array()) * phi;

        double neg_ll = - (
            (m_n * R::lgammafn(phi)) -
            mu_phi.unaryExpr([](double x){ return R::lgammafn(x); }).sum() -
            one_minus_mu_phi.unaryExpr([](double x){ return R::lgammafn(x); }).sum() +
            ((mu_phi.array() - 1.0) * m_log_y.array()).sum() +
            ((one_minus_mu_phi.array() - 1.0) * m_log1_y.array()).sum()
        );
        
        // Gradient
        grad.resize(m_p + 1);

        Eigen::VectorXd d_mu_d_eta = mu.array() * (1.0 - mu.array());
        
        // d(-ll)/dmu = phi * [digamma(mu*phi) - digamma((1-mu)*phi) - log(y) + log(1-y)]
        Eigen::VectorXd d_neg_ll_d_mu = (
            mu_phi.unaryExpr(DigammaFunctor()).array() - 
            one_minus_mu_phi.unaryExpr(DigammaFunctor()).array() - 
            m_log_y.array() + m_log1_y.array()
        ) * phi;

        // d(-ll)/dbeta = X^T * (d(-ll)/dmu * dmu/deta)
        grad.head(m_p) = m_X.transpose() * (d_neg_ll_d_mu.array() * d_mu_d_eta.array()).matrix();
        
        // d(-ll)/dlog_phi = [ -n*digamma(phi) + sum(mu*digamma(mu*phi) + (1-mu)*digamma((1-mu)*phi) - mu*log(y) - (1-mu)*log(1-y)) ] * phi
        double d_neg_ll_d_phi = (
            -m_n * R::digamma(phi) + 
            (mu.array() * mu_phi.unaryExpr(DigammaFunctor()).array()).sum() +
            ((1.0 - mu.array()) * one_minus_mu_phi.unaryExpr(DigammaFunctor()).array()).sum() -
            (mu.array() * m_log_y.array()).sum() -
            ((1.0 - mu.array()) * m_log1_y.array()).sum()
        );
        grad[m_p] = d_neg_ll_d_phi * phi;

        return neg_ll;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& params) {
        Eigen::MatrixXd H(m_p + 1, m_p + 1);
        H.setZero();

        double h = 1e-6; // Step size for finite difference
        
        Eigen::VectorXd grad_at_params(m_p + 1);
        operator()(params, grad_at_params);

        for (int i = 0; i < m_p + 1; ++i) {
            Eigen::VectorXd params_plus_h_i = params;
            params_plus_h_i[i] += h;
            Eigen::VectorXd grad_plus_h_i(m_p + 1);
            operator()(params_plus_h_i, grad_plus_h_i);

            H.col(i) = (grad_plus_h_i - grad_at_params) / h;
        }
        H = (H + H.transpose()) / 2.0;
        
        return H;
    }
};

//' Fast Beta Regression using Rcpp and L-BFGS for joint MLE
//'
//' @param X Model matrix.
//' @param y Numeric vector of proportions (0 < y < 1).
//' @param start_beta Optional starting values for beta coefficients.
//' @param start_phi Optional starting value for phi.
//' @param compute_std_errs Logical, whether to compute and return standard errors (currently ignored in this function, use _with_var_cpp for SEs).
//' @return A list with coefficients, phi, and other optimization details.
//' @export
// [[Rcpp::export]]
List fast_beta_regression_cpp(Eigen::MatrixXd X, 
                              NumericVector y,                                
                              Nullable<NumericVector> start_beta = R_NilValue,
                              double start_phi = 10.0,
                              bool compute_std_errs = false) {
                                  
    Eigen::VectorXd y_eigen = Rcpp::as<Eigen::VectorXd>(y);
    int p = X.cols();
    
    // Initial parameters
    Eigen::VectorXd params(p + 1);
    if (start_beta.isNotNull()) {
        params.head(p) = Rcpp::as<Eigen::VectorXd>(start_beta);
    } else {
        params.head(p).setZero();
    }
    params[p] = std::log(start_phi);

    BetaRegression fun(y_eigen, X);
    LBFGSParam<double> lbfgs_params;
    lbfgs_params.epsilon = 1e-6;
    lbfgs_params.max_iterations = 1000;

    LBFGSSolver<double> solver(lbfgs_params);
    double neg_ll;
    int niter = solver.minimize(fun, params, neg_ll);

    Eigen::VectorXd beta = params.head(p);
    double phi = std::exp(params[p]);

    return List::create(
        Named("coefficients") = beta,
        Named("phi") = phi,
        Named("n_iter") = niter,
        Named("neg_log_lik") = neg_ll,
        Named("params") = params
    );
}

//' Fast Beta Regression with variance using Rcpp and L-BFGS
//'
//' @param X Model matrix.
//' @param y Numeric vector of proportions (0 < y < 1).
//' @param start_beta Optional starting values for beta coefficients.
//' @param start_phi Optional starting value for phi.
//' @param compute_std_errs Logical, whether to compute and return standard errors.
//' @return A list with coefficients, phi, std_errs, vcov, etc.
//' @export
// [[Rcpp::export]]
List fast_beta_regression_with_var_cpp(Eigen::MatrixXd X,
                                     NumericVector y,                                       
                                     Nullable<NumericVector> start_beta = R_NilValue,
                                     double start_phi = 10.0,
                                     bool compute_std_errs = true) {
    
    List fit = fast_beta_regression_cpp(X, y, start_beta, start_phi);
    Eigen::VectorXd params = fit["params"];
    
    Eigen::VectorXd y_eigen = Rcpp::as<Eigen::VectorXd>(y);
    BetaRegression fun(y_eigen, X);
    Eigen::MatrixXd H = fun.hessian(params);
    Eigen::MatrixXd cov_mat = H.inverse();
    
    fit["vcov"] = cov_mat;
    fit["std_errs"] = cov_mat.diagonal().array().sqrt();
    fit.erase(fit.findName("params"));

    return fit;
}
