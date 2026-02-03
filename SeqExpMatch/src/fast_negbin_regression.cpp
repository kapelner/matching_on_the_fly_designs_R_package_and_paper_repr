#include "_helper_functions.h"

using namespace Eigen;
using namespace Rcpp;

// Negative binomial log-likelihood with gradient (NO CENSORING)
class NBLogLik {
private:
    const MatrixXd& X;
    const VectorXi& y;
    int p;

public:
    NBLogLik(const MatrixXd& X_, const VectorXi& y_)
        : X(X_), y(y_), p(X_.cols()) {}

    double f_grad(const VectorXd& pars, VectorXd& grad) {
        //Rcpp::Rcout << "pars: " << pars.transpose() << std::endl;
        const VectorXd beta = pars.head(p);
        const double log_theta = pars[p];
        const double theta = std::exp(log_theta);

        VectorXd eta = X * beta;
        VectorXd mu = eta.array().exp();

        double nll = 0.0;
        VectorXd score_beta = VectorXd::Zero(p);
        double score_log_theta = 0.0;

        for (int i = 0; i < X.rows(); ++i) {
            double mu_i = mu[i];
            double r = theta;
            double yi = y[i];

            // Log-likelihood contribution
            nll -= R::dnbinom_mu(yi, r, mu_i, true);

            // Gradient w.r.t. beta
            double coef = yi - mu_i * (yi + r) / (r + mu_i);
            score_beta += coef * X.row(i).transpose();

            // Gradient w.r.t. log(theta)
            double dlogf_dr =
                R::digamma(yi + r) - R::digamma(r) +
                std::log(r) - std::log(r + mu_i) +
                1.0 - (yi + r) / (r + mu_i);
            score_log_theta += r * dlogf_dr;
        }

        grad.head(p) = -score_beta;
        grad[p] = -score_log_theta;

        return nll;
    }
    
    Eigen::MatrixXd hessian(const Eigen::VectorXd& params) {
        Eigen::MatrixXd H(p + 1, p + 1);
        H.setZero();

        double h = 1e-5; // Step size for finite difference
        
        Eigen::VectorXd grad_at_params(p + 1);
        f_grad(params, grad_at_params);

        for (int i = 0; i < p + 1; ++i) {
            Eigen::VectorXd params_plus_h_i = params;
            params_plus_h_i[i] += h;
            Eigen::VectorXd grad_plus_h_i(p + 1);
            f_grad(params_plus_h_i, grad_plus_h_i);

            H.col(i) = (grad_plus_h_i - grad_at_params) / h;
        }
        H = (H + H.transpose()) / 2.0;
        
        return H;
    }
};

//' Fast Negative Binomial Regression with standard deviations using Rcpp and L-BFGS
//'
//' @param X Design matrix.
//' @param y Response vector.
//' @param maxit Maximum iterations.
//' @param eps_f Convergence tolerance for function value.
//' @param eps_g Convergence tolerance for gradient.
//' @return A list with coefficients, theta, log-likelihood, and Hessian.
//' @export
// [[Rcpp::export]]
List fast_neg_bin_with_var_cpp(const Eigen::MatrixXd& X,
                               const Eigen::VectorXi& y,
                               int maxit = 100,
                               double eps_f = 1e-8,
                               double eps_g = 1e-5) {
    NBLogLik fn(X, y);

    int p = X.cols();
    VectorXd beta_and_theta = VectorXd::Zero(p + 1); //initialize theta so we start with the geometric distribution

    double learning_rate = 1e-5;
    double f_prev = std::numeric_limits<double>::infinity();
    int status = 1;

    VectorXd grad(p + 1);

    for (int i = 0; i < maxit; ++i) {
        double f_curr = fn.f_grad(beta_and_theta, grad);
        if (std::isinf(f_curr) || std::isnan(f_curr)) {
            Rcpp::warning("Objective function is NaN or Inf. Stopping.");
            status = 1;
            break;
        }

        beta_and_theta -= learning_rate * grad;

        if (i > 0 && std::abs(f_prev - f_curr) < eps_f) {
            status = 0;
            break;
        }
        f_prev = f_curr;
        if (i == maxit - 1) {
			Rcpp::warning("Max iterations reached without convergence.");
            status = 1;
		}
    }
    VectorXd grad_final(p + 1);
    double fmin = fn.f_grad(beta_and_theta, grad_final);

    VectorXd b = beta_and_theta.head(p);
    double theta_hat = std::exp(beta_and_theta[p]);

    // Calculate Hessian numerically for debugging
    MatrixXd hess_fisher_info_matrix = fn.hessian(beta_and_theta);

    return List::create(
        Named("b") = b,
        Named("theta_hat") = theta_hat,
        Named("logLik") = -fmin,
        Named("converged") = (status == 0),
        Named("hess_fisher_info_matrix") = hess_fisher_info_matrix
    );
}

// [[Rcpp::export]]
List fast_neg_bin_cpp(const Eigen::MatrixXd& X,
                      const Eigen::VectorXi& y,
                      int maxit = 100,
                      double eps_f = 1e-8,
                      double eps_g = 1e-5) {

    List res = fast_neg_bin_with_var_cpp(X, y, maxit, eps_f, eps_g);

    // Return only the needed elements (exclude hess_fisher_info_matrix)
    return List::create(
        Named("b") = res["b"],
        Named("theta_hat") = res["theta_hat"],
        Named("logLik") = res["logLik"],
        Named("converged") = res["converged"]
    );
}