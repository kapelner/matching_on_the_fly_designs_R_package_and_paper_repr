#include "_helper_functions.h"
#include <RcppEigen.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// --- Weight Functions ---

// Huber weight function
double huber_w(double r, double c) {
    double abs_r = std::abs(r);
    if (abs_r <= c) return 1.0;
    return c / abs_r;
}

// Tukey's Bisquare (Biweight) weight function
double bisquare_w(double r, double c) {
    double abs_r = std::abs(r);
    if (abs_r <= c) {
        double tmp = 1.0 - (r / c) * (r / c);
        return tmp * tmp;
    }
    return 0.0;
}

// --- Internal IRLS Logic ---

struct RobustModelResult {
    Eigen::VectorXd b;
    Eigen::VectorXd w;
    Eigen::MatrixXd XtWX;
    double scale;
    int iterations;
    bool converged;
    double ssq_b_j;

    RobustModelResult() : scale(NA_REAL), iterations(0), converged(false), ssq_b_j(NA_REAL) {}
};

RobustModelResult fast_robust_regression_internal(
    const Eigen::MatrixXd& X, 
    const Eigen::VectorXd& y, 
    Nullable<NumericVector> start_beta = R_NilValue,
    std::string method = "MM",
    double c = 1.345, // Huber constant
    double c_bisquare = 4.685, // Bisquare constant
    int maxit = 50,
    double tol = 1e-7,
    double scale_est = -1.0 // If negative, compute MAD
) {
    int n = X.rows();
    int p = X.cols();
    RobustModelResult res;

    // 1. Initial estimate
    if (start_beta.isNotNull()) {
        res.b = as<Eigen::VectorXd>(NumericVector(start_beta));
    } else {
        Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(X);
        res.b = cod.solve(y);
    }

    Eigen::VectorXd r = y - X * res.b;
    
    // 2. Scale estimation (MAD of residuals)
    if (scale_est < 0) {
        std::vector<double> abs_r(n);
        for (int i = 0; i < n; ++i) abs_r[i] = std::abs(r[i]);
        std::sort(abs_r.begin(), abs_r.end());
        double median_abs_r = (n % 2 == 0) ? (abs_r[n / 2 - 1] + abs_r[n / 2]) / 2.0 : abs_r[n / 2];
        res.scale = median_abs_r / 0.6745;
    } else {
        res.scale = scale_est;
    }

    if (res.scale < 1e-10) res.scale = 1e-10;

    // 3. IRLS loop
    res.w = Eigen::VectorXd::Ones(n);
    Eigen::VectorXd b_old = res.b;
    
    for (int iter = 1; iter <= maxit; ++iter) {
        res.iterations = iter;
        
        // Update weights
        for (int i = 0; i < n; ++i) {
            double u = r[i] / res.scale;
            if (method == "M") {
                res.w[i] = huber_w(u, c);
            } else {
                res.w[i] = bisquare_w(u, c_bisquare);
            }
        }

        // Solve Weighted Least Squares
        Eigen::MatrixXd XtWX = X.transpose() * res.w.asDiagonal() * X;
        Eigen::VectorXd XtWy = X.transpose() * res.w.asDiagonal() * y;
        
        Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
        if (ldlt.info() != Eigen::Success) break; // Numerical failure
        
        res.b = ldlt.solve(XtWy);
        r = y - X * res.b;

        // Check convergence
        if ((res.b - b_old).norm() / (res.b.norm() + 1e-10) < tol) {
            res.converged = true;
            res.XtWX = XtWX;
            break;
        }
        b_old = res.b;
    }

    return res;
}

// [[Rcpp::export]]
List fast_robust_regression_cpp(
    const Eigen::MatrixXd& X, 
    const Eigen::VectorXd& y, 
    Nullable<NumericVector> start_beta = R_NilValue,
    std::string method = "MM",
    int j = 2,
    double c = 1.345,
    int maxit = 50,
    double tol = 1e-7
) {
    RobustModelResult res = fast_robust_regression_internal(X, y, start_beta, method, c, 4.685, maxit, tol);
    
    if (!res.converged && res.XtWX.rows() == 0) {
        res.XtWX = X.transpose() * res.w.asDiagonal() * X;
    }

    int n = X.rows();
    int p = X.cols();
    Eigen::VectorXd r = y - X * res.b;
    
    double ssq_j = NA_REAL;
    if (res.converged || res.iterations == maxit) {
        Eigen::VectorXd psi_r(n);
        double sum_psi_prime = 0;
        
        for (int i = 0; i < n; ++i) {
            double u = r[i] / res.scale;
            double abs_u = std::abs(u);
            if (method == "M") {
                if (abs_u <= c) {
                    psi_r[i] = r[i];
                    sum_psi_prime += 1.0;
                } else {
                    psi_r[i] = c * res.scale * (u > 0 ? 1.0 : -1.0);
                    sum_psi_prime += 0.0;
                }
            } else {
                double c_b = 4.685;
                if (abs_u <= c_b) {
                    double tmp = 1.0 - (u/c_b)*(u/c_b);
                    psi_r[i] = r[i] * tmp * tmp;
                    sum_psi_prime += tmp * (1.0 - 5.0 * (u/c_b)*(u/c_b));
                } else {
                    psi_r[i] = 0;
                    sum_psi_prime += 0;
                }
            }
        }
        
        double m = sum_psi_prime / n;
        Eigen::MatrixXd XtX = X.transpose() * X;
        Eigen::MatrixXd invXtX = XtX.inverse();
        
        double sum_psi_sq = psi_r.squaredNorm();
        double factor = (n / (double(n - p))) * sum_psi_sq / (n * m * m);
        
        if (j > 0 && j <= p) {
            ssq_j = factor * invXtX(j-1, j-1);
        }
    }

    return List::create(
        Named("coefficients") = res.b,
        Named("scale") = res.scale,
        Named("converged") = res.converged,
        Named("iterations") = res.iterations,
        Named("ssq_b_j") = ssq_j
    );
}
