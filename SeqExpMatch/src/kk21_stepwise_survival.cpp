// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;
using namespace Eigen;

// Helper: Fit Weibull AFT on the reduced design matrix XS (intercept in col 0)
// using alternating IRLS (beta) and Newton-Raphson (scale). beta and scale are warm-started.
// On success, fills W_diag = exp(z_hat), resid = delta - W_diag, and XtWX_ldlt.
// Returns false if rank-deficient or too few events.
static bool weibull_reduced_fit_for_score_test(
    const Eigen::MatrixXd& XS,
    const Eigen::VectorXd& log_y,
    const Eigen::VectorXd& delta,
    Eigen::VectorXd& beta,
    double& scale,
    Eigen::VectorXd& W_diag,
    Eigen::VectorXd& resid,
    Eigen::LDLT<Eigen::MatrixXd>& XtWX_ldlt,
    int maxit = 50,
    double tol = 1e-6
) {
    int n = XS.rows();
    int m = XS.cols();
    if (n <= m + 1) return false;
    if (delta.sum() < 2) return false;

    scale = std::max(0.01, std::min(10.0, scale));
    W_diag.resize(n);

    for (int iter = 0; iter < maxit; ++iter) {
        // Standardized log-residuals z_i = (log(t_i) - XS*beta) / scale
        Eigen::VectorXd z = (log_y - XS * beta) / scale;
        Eigen::VectorXd w(n), adj(n);
        for (int i = 0; i < n; ++i) {
            double zi    = std::max(-20.0, std::min(20.0, z(i)));
            double exp_z = std::exp(zi);
            w(i)   = exp_z;
            adj(i) = (delta(i) > 0.5) ? (zi - 1.0 + exp_z) : exp_z;
        }

        // IRLS: weighted least squares update for beta
        Eigen::VectorXd sqrt_w = w.array().sqrt();
        Eigen::MatrixXd Xw = XS.array().colwise() * sqrt_w.array();
        Eigen::VectorXd yw(n);
        for (int i = 0; i < n; ++i)
            yw(i) = (log_y(i) - scale * adj(i) / w(i)) * sqrt_w(i);

        Eigen::LDLT<Eigen::MatrixXd> ldlt_inner((Xw.transpose() * Xw).eval());
        if (ldlt_inner.info() != Eigen::Success) return false;
        Eigen::VectorXd beta_new = ldlt_inner.solve(Xw.transpose() * yw);

        // Newton-Raphson step for scale
        Eigen::VectorXd z_new = (log_y - XS * beta_new) / scale;
        double score_s = 0.0, info_s = 0.0;
        for (int i = 0; i < n; ++i) {
            double zi     = std::max(-20.0, std::min(20.0, z_new(i)));
            double exp_zi = std::exp(zi);
            if (delta(i) > 0.5) {
                score_s += -1.0/scale + zi/scale - zi*exp_zi/scale;
                info_s  +=  1.0/(scale*scale);
            } else {
                score_s += -zi*exp_zi/scale;
            }
        }
        double scale_new = scale;
        if (info_s > 1e-10) {
            scale_new = scale - score_s / info_s;
            scale_new = std::max(0.01, std::min(10.0, scale_new));
        }

        double diff = (beta_new - beta).norm() + std::fabs(scale_new - scale);
        beta  = beta_new;
        scale = scale_new;
        if (diff < tol) break;
    }

    // Final W and score-test residuals
    Eigen::VectorXd z_fin = (log_y - XS * beta) / scale;
    for (int i = 0; i < n; ++i) {
        double zi = std::max(-20.0, std::min(20.0, z_fin(i)));
        W_diag(i) = std::exp(zi);
    }
    resid = delta - W_diag;

    Eigen::MatrixXd XtWX = XS.transpose() * W_diag.asDiagonal() * XS;
    XtWX_ldlt = XtWX.ldlt();
    return XtWX_ldlt.info() == Eigen::Success;
}

// Optimized stepwise Weibull AFT using the score test (LM test) from the reduced model.
// At each step k, instead of fitting (p-k) full Weibull AFT regressions of size n x (k+3), we:
//   1. Fit ONE Weibull AFT on the reduced model XS = [1, X_sel, w] (with warm start)
//   2. For each candidate x_j compute the Weibull score test statistic:
//        stat_j = |x_j' r| / sqrt(I_jj.S)
//      where r_i = delta_i - W_hat_i,  W_hat_i = exp(z_hat_i),  z_hat_i = (log(t_i) - XS*beta_hat) / sigma_hat
//      and I_jj.S = x_j'Wx_j - (XS'Wx_j)'(XS'WXS)^{-1}(XS'Wx_j)
//      Note: sigma_hat cancels in the statistic.
// Falls back to FWL OLS on log(y) if the Weibull fit fails.
// Asymptotically equivalent to the Wald t-statistic from the full Weibull model.
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

    if (n == 0 || p == 0) return weights;

    Eigen::Map<Eigen::MatrixXd> X_map(as<Eigen::Map<Eigen::MatrixXd>>(X));
    Eigen::VectorXd y_vec   = as<Eigen::VectorXd>(y);
    Eigen::VectorXd delta_v = as<Eigen::VectorXd>(delta);
    Eigen::VectorXd w_vec   = as<Eigen::VectorXd>(w);

    Eigen::VectorXd log_y(n);
    for (int i = 0; i < n; ++i)
        log_y(i) = std::log(std::max(y_vec(i), 1e-10));

    // ---- Weibull score-test path ----
    // Base model XS = [1, w]
    int m = 2;
    Eigen::MatrixXd XS(n, m);
    XS.col(0).setOnes();
    XS.col(1) = w_vec;

    // Initialize beta from OLS on log(y) ~ XS; scale from event residuals
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(m);
    double scale = 1.0;
    {
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr((XS.transpose() * XS).eval());
        if (qr.rank() >= m) {
            beta = qr.solve(XS.transpose() * log_y);
            Eigen::VectorXd r0 = log_y - XS * beta;
            double sse = 0.0, nd = 0.0;
            for (int i = 0; i < n; ++i)
                if (delta_v(i) > 0.5) { sse += r0(i)*r0(i); nd++; }
            if (nd > m)
                scale = std::max(0.01, std::min(10.0, std::sqrt(sse / (nd - m))));
        }
    }

    Eigen::VectorXd W_diag, resid;
    Eigen::LDLT<Eigen::MatrixXd> XtWX_ldlt;
    bool fit_ok = weibull_reduced_fit_for_score_test(
        XS, log_y, delta_v, beta, scale, W_diag, resid, XtWX_ldlt);

    // ---- FWL OLS on log(y) fallback state (lazy init) ----
    const double eps = std::numeric_limits<double>::epsilon();
    Eigen::MatrixXd Q(n, p + 2);
    int qm = 0;
    Eigen::VectorXd e_y_fwl;
    double c_fwl = 0.0;
    bool fwl_initialized = false;

    std::vector<int> selected_so_far;
    selected_so_far.reserve(p);
    std::vector<bool> used(p, false);

    for (int step = 0; step < p; ++step) {

        if (fit_ok) {
            // ---- Weibull score test ----
            Eigen::MatrixXd WXS = W_diag.asDiagonal() * XS;  // n x m

            double best_stat = -1.0;
            int best_j = -1;

            for (int j = 0; j < p; ++j) {
                if (used[j]) continue;
                Eigen::VectorXd xj = X_map.col(j);

                // Score: U_j = x_j' r  (sigma_hat cancels)
                double Uj = xj.dot(resid);

                // Adjusted information: I_jj.S = x_j'Wx_j - (XS'Wx_j)'(XS'WXS)^{-1}(XS'Wx_j)
                Eigen::VectorXd v_j  = WXS.transpose() * xj;
                Eigen::VectorXd h_j  = XtWX_ldlt.solve(v_j);
                double xjWxj = (W_diag.array() * xj.array().square()).sum();
                double Ijj   = xjWxj - v_j.dot(h_j);
                if (Ijj <= 0.0) continue;

                double stat = std::fabs(Uj) / std::sqrt(Ijj);
                if (!R_finite(stat)) stat = 0.0;
                if (stat > best_stat) { best_stat = stat; best_j = j; }
            }

            if (best_j < 0) break;
            weights[best_j] = best_stat;
            used[best_j]    = true;
            selected_so_far.push_back(best_j);

            // Append x_{best_j} to XS, extend beta with 0, keep scale, refit
            Eigen::MatrixXd XS_new(n, m + 1);
            XS_new.leftCols(m) = XS;
            XS_new.col(m)      = X_map.col(best_j);
            XS = XS_new;
            m++;
            Eigen::VectorXd beta_new(m);
            beta_new.head(m - 1) = beta;
            beta_new[m - 1]      = 0.0;
            beta = beta_new;

            fit_ok = weibull_reduced_fit_for_score_test(
                XS, log_y, delta_v, beta, scale, W_diag, resid, XtWX_ldlt);

        } else {
            // ---- FWL OLS on log(y) fallback ----

            // Lazy init: build orthonormal basis from [1, w, already-selected covariates]
            if (!fwl_initialized) {
                Q.col(0).setConstant(1.0 / std::sqrt(static_cast<double>(n)));
                qm = 1;

                // Add w
                {
                    Eigen::VectorXd u = w_vec - Q.leftCols(qm) * (Q.leftCols(qm).transpose() * w_vec);
                    double nm = u.norm();
                    if (nm > 1e-8) { Q.col(qm) = u / nm; qm++; }
                }
                // Add covariates already selected by Weibull path (if any)
                for (int sj : selected_so_far) {
                    Eigen::VectorXd v = X_map.col(sj);
                    Eigen::VectorXd u = v - Q.leftCols(qm) * (Q.leftCols(qm).transpose() * v);
                    double nm = u.norm();
                    if (nm > 1e-8) { Q.col(qm) = u / nm; qm++; }
                }
                e_y_fwl = log_y - Q.leftCols(qm) * (Q.leftCols(qm).transpose() * log_y);
                c_fwl = e_y_fwl.squaredNorm();
                fwl_initialized = true;
            }

            int df = n - qm - 1;
            if (df <= 0) break;

            double best_stat = -1.0;
            int best_j = -1;

            for (int j = 0; j < p; ++j) {
                if (used[j]) continue;
                Eigen::VectorXd xj = X_map.col(j);
                Eigen::VectorXd e_j = xj - Q.leftCols(qm) * (Q.leftCols(qm).transpose() * xj);
                double b = e_j.squaredNorm();
                if (b < eps) continue;

                double a    = e_j.dot(e_y_fwl);
                double sse  = c_fwl - a * a / b;
                if (sse < 0.0) sse = 0.0;
                double sigma2 = sse / static_cast<double>(df);
                if (sigma2 <= 0.0 || !std::isfinite(sigma2)) continue;

                double stat = std::fabs(a) / std::sqrt(b * sigma2);
                if (!R_finite(stat)) stat = 0.0;
                if (stat > best_stat) { best_stat = stat; best_j = j; }
            }

            if (best_j < 0) break;
            weights[best_j] = best_stat;
            used[best_j]    = true;
            selected_so_far.push_back(best_j);

            // Update FWL basis with x_{best_j}
            if (qm < p + 2) {
                Eigen::VectorXd xbest = X_map.col(best_j);
                Eigen::VectorXd u = xbest - Q.leftCols(qm) * (Q.leftCols(qm).transpose() * xbest);
                double nm = u.norm();
                if (nm > 1e-8) {
                    Q.col(qm) = u / nm;
                    double proj = Q.col(qm).dot(e_y_fwl);
                    e_y_fwl -= proj * Q.col(qm);
                    c_fwl   -= proj * proj;
                    if (c_fwl < 0.0) c_fwl = 0.0;
                    qm++;
                }
            }
        }
    }

    return weights;
}
