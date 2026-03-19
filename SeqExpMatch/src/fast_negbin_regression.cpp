#include <Rcpp.h>
#include <RcppEigen.h>
#include <optimization/LBFGS.h>
#include <Rmath.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace LBFGSpp;

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
            
            // Negative log-likelihood
            neg_ll -= R::dnbinom_mu(yi, theta, mu_i, true);

            // Gradient components
            double coef = yi - mu_i * (yi + theta) / (theta + mu_i);
            score_beta += coef * m_X.row(i).transpose();

            double dlogf_dtheta =
                R::digamma(yi + theta) - R::digamma(theta) +
                std::log(theta) - std::log(theta + mu_i) +
                1.0 - (yi + theta) / (theta + mu_i);
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

// [[Rcpp::export]]
List fast_neg_bin_with_var_cpp(Eigen::MatrixXd X,
                                Eigen::VectorXi y,
                                int maxit = 1000,
                                double eps_f = 1e-8,
                                double eps_g = 1e-5) {
    int p = X.cols();
    Eigen::VectorXd params = Eigen::VectorXd::Zero(p + 1);
    params[p] = 0.0; // log(theta) = 0 -> theta = 1

    NBLogLik fun(X, y);
    LBFGSParam<double> lbfgs_params;
    lbfgs_params.epsilon = eps_g;
    lbfgs_params.max_iterations = maxit;

    LBFGSSolver<double> solver(lbfgs_params);
    double neg_ll;
    int niter = solver.minimize(fun, params, neg_ll);

    Eigen::VectorXd beta = params.head(p);
    double theta = std::exp(params[p]);
    Eigen::MatrixXd H = fun.hessian(params);

    return List::create(
        Named("b") = beta,
        Named("theta_hat") = theta,
        Named("logLik") = -neg_ll,
        Named("converged") = (niter < maxit),
        Named("hess_fisher_info_matrix") = H,
        Named("niter") = niter
    );
}

// [[Rcpp::export]]
List fast_neg_bin_cpp(Eigen::MatrixXd X,
                        Eigen::VectorXi y,
                        int maxit = 1000,
                        double eps_f = 1e-8,
                        double eps_g = 1e-5) {
    List res = fast_neg_bin_with_var_cpp(X, y, maxit, eps_f, eps_g);
    return List::create(
        Named("b") = res["b"],
        Named("theta_hat") = res["theta_hat"],
        Named("logLik") = res["logLik"],
        Named("converged") = res["converged"]
    );
}
