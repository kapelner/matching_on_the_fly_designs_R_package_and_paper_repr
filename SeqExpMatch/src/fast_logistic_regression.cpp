#include "_helper_functions.h"
#include <RcppEigen.h>

using namespace Rcpp;

// Internal pure C++ logic
ModelResult fast_logistic_regression_internal(const Eigen::MatrixXd& X, 
                                              const Eigen::VectorXd& y, 
                                              const Eigen::VectorXd& weights = Eigen::VectorXd(),
                                              int maxit = 100, 
                                              double tol = 1e-8) {
    int n = X.rows();
    int p = X.cols();
    bool use_weights = (weights.size() == n);
    
    ModelResult res;
    res.b = Eigen::VectorXd::Zero(p);
    res.mu = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd w(n);
    Eigen::VectorXd z(n);

    for (int iter = 0; iter < maxit; iter++) {
        Eigen::VectorXd eta = X * res.b;
        for(int i=0; i<n; i++) {
            double val = 1.0 / (1.0 + std::exp(-eta[i]));
            res.mu[i] = val;
            double base_w = std::max(val * (1.0 - val), 1e-10);
            w[i] = use_weights ? weights[i] * base_w : base_w;
        }

        z = eta + (y - res.mu).cwiseQuotient(use_weights ? (res.mu.array() * (1.0 - res.mu.array())).matrix() : w);
        // Note: The adjusted dependent variable z in IRLS for binomial is eta + (y - mu) / [mu(1-mu)]
        // The weights in W are weights_i * mu_i * (1-mu_i)
        
        Eigen::MatrixXd XtW = X.transpose() * w.asDiagonal();
        res.XtWX = XtW * X;
        Eigen::VectorXd XtWz = XtW * z;

        Eigen::VectorXd beta_new = res.XtWX.ldlt().solve(XtWz);

        if ((res.b - beta_new).norm() < tol) {
            res.b = beta_new;
            res.converged = true;
            break;
        }
        res.b = beta_new;
    }

    // Final update of mu and XtWX
    Eigen::VectorXd eta = X * res.b;
    for(int i=0; i<n; i++) {
        double val = 1.0 / (1.0 + std::exp(-eta[i]));
        res.mu[i] = val;
        double base_w = std::max(val * (1.0 - val), 1e-10);
        w[i] = use_weights ? weights[i] * base_w : base_w;
    }
    res.XtWX = X.transpose() * w.asDiagonal() * X;

    return res;
}

// [[Rcpp::export]]
List fast_logistic_regression_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int maxit = 100, double tol = 1e-8) {
    ModelResult res = fast_logistic_regression_internal(X, y, Eigen::VectorXd(), maxit, tol);
    return List::create(
        Named("b") = res.b,
        Named("w") = res.mu.array() * (1.0 - res.mu.array()) // return weights as expected by old code
    );
}

// [[Rcpp::export]]
List fast_logistic_regression_weighted_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& weights, int maxit = 100, double tol = 1e-8) {
    ModelResult res = fast_logistic_regression_internal(X, y, weights, maxit, tol);
    return List::create(
        Named("b") = res.b,
        Named("mu") = res.mu,
        Named("XtWX") = res.XtWX,
        Named("converged") = res.converged
    );
}

// [[Rcpp::export]]
List fast_logistic_regression_with_var_cpp(const Eigen::MatrixXd& Xmm, const Eigen::VectorXd& y, int j = 2) {
    ModelResult res = fast_logistic_regression_internal(Xmm, y);
    res.ssq_b_j = compute_diagonal_inverse_entry(res.XtWX, j);
    res.ssq_b_2 = (Xmm.cols() >= 2) ? compute_diagonal_inverse_entry(res.XtWX, 2) : NA_REAL;

    return List::create(
        Named("b") = res.b,
        Named("ssq_b_j") = res.ssq_b_j,
        Named("ssq_b_2") = res.ssq_b_2
    );
}
