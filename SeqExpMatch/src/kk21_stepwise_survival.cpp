// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace Eigen;

// Forward declaration from _helper_functions.cpp
List fast_ols_with_var_cpp(const Eigen::MatrixXd& Xmm,
                           const Eigen::VectorXd& y,
                           int j_for_var = 2);

// Helper: Weibull AFT MLE using Newton-Raphson for scale parameter
// Model: log(T) = intercept + X*beta + scale*epsilon, epsilon ~ standard extreme value
// For censored data, we use profile likelihood approach
static bool weibull_aft_fit(
    const VectorXd& log_y,      // log survival times
    const VectorXd& delta,      // censoring indicator (1=event, 0=censored)
    const MatrixXd& X,          // design matrix (without intercept)
    VectorXd& beta,             // output: regression coefficients
    double& scale,              // output: scale parameter
    double& se_beta1,           // output: standard error of first covariate
    int maxit = 50,
    double tol = 1e-6
) {
    int n = log_y.size();
    int p = X.cols();

    if (n < p + 2) {
        return false;
    }

    // Count events
    double n_events = delta.sum();
    if (n_events < 2) {
        return false;
    }

    // Initialize with OLS on log(y)
    MatrixXd Xfull(n, p + 1);
    Xfull.col(0).setOnes();
    Xfull.rightCols(p) = X;

    // OLS solution for initial beta
    MatrixXd XtX = Xfull.transpose() * Xfull;
    ColPivHouseholderQR<MatrixXd> qr(XtX);
    if (qr.rank() < p + 1) {
        return false;
    }
    VectorXd Xty = Xfull.transpose() * log_y;
    VectorXd beta_full = qr.solve(Xty);

    // Initial residuals
    VectorXd resid = log_y - Xfull * beta_full;

    // Initial scale estimate (robust)
    double scale_init = 0.0;
    for (int i = 0; i < n; ++i) {
        if (delta(i) > 0.5) {  // uncensored
            scale_init += resid(i) * resid(i);
        }
    }
    scale_init = std::sqrt(scale_init / std::max(1.0, n_events - p - 1));
    if (scale_init < 0.01) scale_init = 0.01;
    if (scale_init > 10.0) scale_init = 10.0;
    scale = scale_init;

    // Newton-Raphson iteration for scale (profile likelihood)
    for (int iter = 0; iter < maxit; ++iter) {
        // Compute standardized residuals
        VectorXd z = resid / scale;

        // Compute weights for IRLS
        VectorXd w(n);
        VectorXd adj(n);
        for (int i = 0; i < n; ++i) {
            double exp_z = std::exp(z(i));
            if (delta(i) > 0.5) {  // uncensored
                w(i) = exp_z;
                adj(i) = z(i) - 1.0 + exp_z;
            } else {  // censored
                w(i) = exp_z;
                adj(i) = exp_z;
            }
        }

        // Update beta using weighted least squares
        VectorXd sqrt_w = w.array().sqrt();
        MatrixXd Xw = Xfull.array().colwise() * sqrt_w.array();
        VectorXd yw = (log_y.array() - scale * adj.array() / w.array()) * sqrt_w.array();

        MatrixXd XwXw = Xw.transpose() * Xw;
        ColPivHouseholderQR<MatrixXd> qr_w(XwXw);
        if (qr_w.rank() < p + 1) {
            return false;
        }
        VectorXd beta_new = qr_w.solve(Xw.transpose() * yw);

        // Update residuals
        VectorXd resid_new = log_y - Xfull * beta_new;

        // Update scale using MLE
        double sum_d = 0.0;
        double sum_wz = 0.0;
        for (int i = 0; i < n; ++i) {
            double z_new = resid_new(i) / scale;
            double exp_z_new = std::exp(z_new);
            if (delta(i) > 0.5) {
                sum_d += 1.0;
                sum_wz += z_new * exp_z_new - z_new + exp_z_new;
            } else {
                sum_wz += exp_z_new;
            }
        }

        // Scale update (simplified)
        double scale_new = scale;
        if (sum_d > 0) {
            // Newton step for scale
            double score = 0.0;
            double info = 0.0;
            for (int i = 0; i < n; ++i) {
                double z_i = resid_new(i) / scale;
                double exp_z_i = std::exp(z_i);
                if (delta(i) > 0.5) {
                    score += -1.0/scale + z_i/scale - z_i*exp_z_i/scale;
                    info += 1.0/(scale*scale);
                } else {
                    score += -z_i*exp_z_i/scale;
                }
            }
            if (std::fabs(info) > 1e-10) {
                scale_new = scale - score / info;
                scale_new = std::max(0.01, std::min(10.0, scale_new));
            }
        }

        // Check convergence
        double diff_beta = (beta_new - beta_full).norm();
        double diff_scale = std::fabs(scale_new - scale);

        beta_full = beta_new;
        scale = scale_new;
        resid = resid_new;

        if (diff_beta < tol && diff_scale < tol) {
            break;
        }
    }

    // Extract beta (without intercept)
    beta = beta_full.tail(p);

    // Compute standard error of first covariate coefficient
    VectorXd z = resid / scale;
    VectorXd w(n);
    for (int i = 0; i < n; ++i) {
        double exp_z = std::exp(z(i));
        w(i) = exp_z / (scale * scale);
    }

    MatrixXd Xw = Xfull.array().colwise() * w.array().sqrt();
    MatrixXd info_mat = Xw.transpose() * Xw;

    FullPivLU<MatrixXd> lu(info_mat);
    if (!lu.isInvertible()) {
        return false;
    }
    MatrixXd cov_mat = lu.inverse();

    if (cov_mat(1, 1) <= 0) {
        return false;
    }
    se_beta1 = std::sqrt(cov_mat(1, 1));

    return true;
}

// [[Rcpp::export]]
NumericVector kk21_stepwise_survival_weights_cpp(
    const NumericMatrix& X,      // Covariate matrix (n x p)
    const NumericVector& y,      // Survival times
    const NumericVector& delta,  // Censoring indicator (1=event, 0=censored)
    const NumericVector& w       // Treatment indicator
) {
    int n = X.nrow();
    int p = X.ncol();
    NumericVector weights(p, NA_REAL);
    double eps = std::numeric_limits<double>::epsilon();

    if (n == 0 || p == 0) {
        return weights;
    }

    // Prepare data
    Eigen::Map<Eigen::MatrixXd> X_map(as<Eigen::Map<Eigen::MatrixXd>>(X));
    Eigen::VectorXd y_vec = as<Eigen::VectorXd>(y);
    Eigen::VectorXd delta_vec = as<Eigen::VectorXd>(delta);
    Eigen::VectorXd w_vec = as<Eigen::VectorXd>(w);

    // Compute log(y), handling zeros
    Eigen::VectorXd log_y(n);
    for (int i = 0; i < n; ++i) {
        log_y(i) = std::log(std::max(y_vec(i), 1e-10));
    }

    // Stepwise selection
    std::vector<int> selected;
    selected.reserve(p);
    std::vector<bool> used(p, false);

    for (int step = 0; step < p; ++step) {
        double best_stat = -1.0;
        int best_j = -1;
        int k = static_cast<int>(selected.size());

        for (int j = 0; j < p; ++j) {
            if (used[j]) {
                continue;
            }

            // Build design matrix: [x_j, selected covariates, w]
            Eigen::MatrixXd Xmm(n, k + 2);
            Xmm.col(0) = X_map.col(j);
            for (int idx = 0; idx < k; ++idx) {
                Xmm.col(1 + idx) = X_map.col(selected[idx]);
            }
            Xmm.col(k + 1) = w_vec;

            double stat = 0.0;

            // Try Weibull AFT first
            Eigen::VectorXd beta;
            double scale, se_beta1;
            bool weibull_ok = weibull_aft_fit(log_y, delta_vec, Xmm, beta, scale, se_beta1, 30, 1e-5);

            if (weibull_ok && se_beta1 > 0 && std::isfinite(se_beta1)) {
                stat = std::fabs(beta(0) / se_beta1);
            }

            // Fall back to OLS on log(y) if Weibull failed
            if (!weibull_ok || !std::isfinite(stat) || stat <= 0) {
                // OLS fallback: use log(y) ~ x_j + selected + w
                Eigen::MatrixXd Xols(n, k + 3);
                Xols.col(0).setOnes();
                Xols.col(1) = X_map.col(j);
                for (int idx = 0; idx < k; ++idx) {
                    Xols.col(2 + idx) = X_map.col(selected[idx]);
                }
                Xols.col(k + 2) = w_vec;

                List mod = fast_ols_with_var_cpp(Xols, log_y, 2);
                Eigen::VectorXd b = mod["b"];
                double ssq = as<double>(mod["ssq_b_j"]);
                stat = std::fabs(b(1) / std::sqrt(ssq));
            }

            if (!std::isfinite(stat)) {
                stat = 0.0;
            }
            if (stat > best_stat) {
                best_stat = stat;
                best_j = j;
            }
        }

        if (best_j < 0) {
            break;
        }

        weights[best_j] = best_stat;
        used[best_j] = true;
        selected.push_back(best_j);
    }

    return weights;
}
