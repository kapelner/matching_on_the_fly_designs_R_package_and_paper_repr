#include "_helper_functions.h"
#include <RcppEigen.h>
#include <cmath>
#include <vector>
#include <string>

using namespace Rcpp;
using namespace Eigen;

ModelResult fast_logistic_regression_internal(const Eigen::MatrixXd& X_eigen,
                                              const Eigen::VectorXd& y_eigen,
                                              const Eigen::VectorXd& weights_eigen = Eigen::VectorXd(),
                                              int maxit = 100,
                                              double tol = 1e-8,
                                              Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                              Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                              std::string optimization_alg = "irls");
ModelResult fast_poisson_regression_internal(const Eigen::MatrixXd& X,
                                             const Eigen::VectorXd& y,
                                             const Eigen::VectorXd& weights = Eigen::VectorXd(),
                                             int maxit = 100,
                                             double tol = 1e-8,
                                             Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                             Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                             std::string optimization_alg = "lbfgs");

enum class GEEFamily { GAUSSIAN, BINOMIAL, POISSON };

struct GEEResult {
    VectorXd beta;
    double alpha;
    MatrixXd vcov;
    double quasi_loglik;
    bool converged;
    int niter;
};

inline double gee_variance(double mu, GEEFamily family) {
    if (family == GEEFamily::BINOMIAL) return std::max(mu * (1.0 - mu), 1e-10);
    if (family == GEEFamily::POISSON)  return std::max(mu, 1e-10);
    return 1.0;
}

inline double gee_link_inv(double eta, GEEFamily family) {
    if (family == GEEFamily::BINOMIAL) return 1.0 / (1.0 + std::exp(-std::min(std::max(eta, -20.0), 20.0)));
    if (family == GEEFamily::POISSON)  return std::exp(std::min(eta, 20.0));
    return eta;
}

inline double gee_dmu_deta(double mu, GEEFamily family) {
    if (family == GEEFamily::BINOMIAL) return mu * (1.0 - mu);
    if (family == GEEFamily::POISSON)  return mu;
    return 1.0;
}

inline double gee_quasi_loglik_contrib(double y, double mu, GEEFamily family) {
    if (family == GEEFamily::BINOMIAL) return y * std::log(std::max(mu, 1e-10)) + (1.0 - y) * std::log(std::max(1.0 - mu, 1e-10));
    if (family == GEEFamily::POISSON)  return y * std::log(std::max(mu, 1e-10)) - mu;
    return -0.5 * (y - mu) * (y - mu);
}

inline VectorXd gee_solve_system(const MatrixXd& A, const VectorXd& b) {
    MatrixXd A_sym = (A + A.transpose()) / 2.0;
    Eigen::LDLT<MatrixXd> ldlt(A_sym);
    if (ldlt.info() == Eigen::Success) {
        VectorXd x = ldlt.solve(b);
        if (x.allFinite()) return x;
    }
    return A.colPivHouseholderQr().solve(b);
}

inline MatrixXd gee_inverse_system(const MatrixXd& A) {
    return covariance_from_information((A + A.transpose()) / 2.0);
}

inline double gee_estimate_exchangeable_alpha(const VectorXd& resid, const VectorXd& mu,
                                              const std::vector<int>& grp_start, const std::vector<int>& grp_size,
                                              GEEFamily family) {
    double num = 0.0, den = 0.0;
    for (int gi = 0; gi < grp_start.size(); ++gi) {
        if (grp_size[gi] < 2) continue;
        int start = grp_start[gi];
        for (int a = 0; a < grp_size[gi]; ++a) {
            for (int b = a + 1; b < grp_size[gi]; ++b) {
                int i = start + a;
                int j = start + b;
                double ri = resid[i] / std::sqrt(gee_variance(mu[i], family));
                double rj = resid[j] / std::sqrt(gee_variance(mu[j], family));
                num += ri * rj;
                den += 1.0;
            }
        }
    }
    if (den <= 0.0) return 0.0;
    return std::min(std::max(num / den, -0.95), 0.95);
}

inline VectorXd gee_fit_independence_glm(const MatrixXd& X, const VectorXd& y, GEEFamily family, int maxit = 100, double tol = 1e-8) {
    ModelResult fit = (family == GEEFamily::BINOMIAL) ?
        fast_logistic_regression_internal(X, y, Eigen::VectorXd(), maxit, tol, R_NilValue, R_NilValue, "irls") :
        fast_poisson_regression_internal(X, y, Eigen::VectorXd(), maxit, tol, R_NilValue, R_NilValue, "irls");
    if (fit.b.size() == X.cols() && fit.b.allFinite()) {
        return fit.b;
    }
    VectorXd beta = VectorXd::Zero(X.cols());
    double y_mean = std::max(y.mean(), 1e-8);
    beta[0] = (family == GEEFamily::POISSON) ? std::log(y_mean) : std::log(std::min(std::max(y_mean, 1e-6), 1.0 - 1e-6) / (1.0 - std::min(std::max(y_mean, 1e-6), 1.0 - 1e-6)));
    return beta;
}

GEEResult gee_pairs_singletons_cpp_impl(const MatrixXd& X, const VectorXd& y, const std::vector<int>& grp_start, const std::vector<int>& grp_size, GEEFamily family, int maxit = 100, double tol = 1e-8) {
    const int n = X.rows(), p = X.cols(), G = grp_start.size();
    for (int gi = 0; gi < G; ++gi) {
        if (grp_size[gi] > 2) stop("gee_pairs_singletons_cpp: cluster %d has size %d (only singletons and pairs are supported)", gi, grp_size[gi]);
    }
    VectorXd beta = VectorXd::Zero(p);
    if (family == GEEFamily::GAUSSIAN) {
        MatrixXd XtX = X.transpose() * X;
        VectorXd Xty = X.transpose() * y;
        beta = gee_solve_system(XtX, Xty);
    } else {
        beta = gee_fit_independence_glm(X, y, family);
    }

    double alpha = 0.0; bool converged = false; int iter = 0;
    for (; iter < maxit; ++iter) {
        VectorXd mu(n), resid(n);
        for (int i = 0; i < n; ++i) { mu[i] = gee_link_inv(X.row(i).dot(beta), family); resid[i] = y[i] - mu[i]; }
        alpha = gee_estimate_exchangeable_alpha(resid, mu, grp_start, grp_size, family);

        MatrixXd Bread = MatrixXd::Zero(p, p); VectorXd Score = VectorXd::Zero(p);
        for (int gi = 0; gi < G; ++gi) {
            int sz = grp_size[gi], start = grp_start[gi];
            if (sz == 1) {
                double mui = mu[start], vi = gee_variance(mui, family), di = gee_dmu_deta(mui, family);
                double w = di * di / vi;
                for(int r=0; r<p; r++) {
                    for(int c=0; c<=r; c++) Bread(r, c) += w * X(start, r) * X(start, c);
                    Score[r] += (di / vi) * resid[start] * X(start, r);
                }
            } else {
                int i1 = start, i2 = start + 1;
                double m1 = mu[i1], m2 = mu[i2], v1 = gee_variance(m1, family), v2 = gee_variance(m2, family);
                double d1 = gee_dmu_deta(m1, family), d2 = gee_dmu_deta(m2, family);
                double s1 = std::sqrt(v1), s2 = std::sqrt(v2), detV = v1 * v2 * (1.0 - alpha * alpha);
                double V11 = v2/detV, V22 = v1/detV, V12 = -(alpha*s1*s2)/detV;
                double r1 = resid[i1], r2 = resid[i2];
                for(int r=0; r<p; r++) {
                    double d1xr = d1 * X(i1, r), d2xr = d2 * X(i2, r);
                    for(int c=0; c<=r; c++) {
                        double d1xc = d1 * X(i1, c), d2xc = d2 * X(i2, c);
                        Bread(r, c) += V11*d1xr*d1xc + V22*d2xr*d2xc + V12*(d1xr*d2xc + d2xr*d1xc);
                    }
                    Score[r] += d1xr*(V11*r1 + V12*r2) + d2xr*(V12*r1 + V22*r2);
                }
            }
        }
        for(int r=0; r<p; r++) for(int c=0; c<r; c++) Bread(c, r) = Bread(r, c);
        VectorXd delta = gee_solve_system(Bread, Score);
        if (!delta.allFinite()) break;
        beta += delta;
        if (delta.norm() < tol) { converged = true; break; }
    }

    MatrixXd Bread = MatrixXd::Zero(p, p), Meat = MatrixXd::Zero(p, p);
    double quasi_loglik = 0.0; VectorXd mu(n), resid(n);
    for (int i = 0; i < n; ++i) { mu[i] = gee_link_inv(X.row(i).dot(beta), family); resid[i] = y[i] - mu[i]; quasi_loglik += gee_quasi_loglik_contrib(y[i], mu[i], family); }
    alpha = gee_estimate_exchangeable_alpha(resid, mu, grp_start, grp_size, family);

    for (int gi = 0; gi < G; ++gi) {
        int sz = grp_size[gi], start = grp_start[gi];
        if (sz == 1) {
            double mui = mu[start], vi = gee_variance(mui, family), di = gee_dmu_deta(mui, family);
            double Vi_inv = 1.0 / vi, w = di * di * Vi_inv;
            for(int r=0; r<p; r++) {
                double sir = (di * Vi_inv * resid[start]) * X(start, r);
                for(int c=0; c<=r; c++) {
                    Bread(r, c) += w * X(start, r) * X(start, c);
                    Meat(r, c) += sir * (di * Vi_inv * resid[start]) * X(start, c);
                }
            }
        } else {
            int i1 = start, i2 = start + 1;
            double v1 = gee_variance(mu[i1], family), v2 = gee_variance(mu[i2], family), d1 = gee_dmu_deta(mu[i1], family), d2 = gee_dmu_deta(mu[i2], family);
            double s1 = std::sqrt(v1), s2 = std::sqrt(v2), detV = v1 * v2 * (1.0 - alpha * alpha);
            double V11 = v2/detV, V22 = v1/detV, V12 = -(alpha*s1*s2)/detV, r1 = resid[i1], r2 = resid[i2];
            for(int r=0; r<p; r++) {
                double d1xr = d1 * X(i1, r), d2xr = d2 * X(i2, r), sir = d1xr*(V11*r1 + V12*r2) + d2xr*(V12*r1 + V22*r2);
                for(int c=0; c<=r; c++) {
                    double d1xc = d1 * X(i1, c), d2xc = d2 * X(i2, c), sic = d1xc*(V11*r1 + V12*r2) + d2xc*(V12*r1 + V22*r2);
                    Bread(r, c) += V11*d1xr*d1xc + V22*d2xr*d2xc + V12*(d1xr*d2xc + d2xr*d1xc);
                    Meat(r, c) += sir * sic;
                }
            }
        }
    }
    for(int r=0; r<p; r++) for(int c=0; c<r; c++) { Bread(c, r) = Bread(r, c); Meat(c, r) = Meat(r, c); }
    MatrixXd BI = gee_inverse_system(Bread);
    return {beta, alpha, BI * Meat * BI, quasi_loglik, converged, iter};
}

// [[Rcpp::export]]
List gee_pairs_singletons_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXi& group_id, std::string family_str, int maxit = 100, double tol = 1e-8) {
    GEEFamily family = GEEFamily::GAUSSIAN;
    if (family_str == "binomial") family = GEEFamily::BINOMIAL;
    else if (family_str == "poisson") family = GEEFamily::POISSON;
    int n = y.size(); std::vector<int> ord(n); for(int i=0; i<n; i++) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](int a, int b){ return group_id[a] < group_id[b]; });
    MatrixXd X_s(n, X.cols()); VectorXd y_s(n);
    for(int i=0; i<n; i++) { X_s.row(i) = X.row(ord[i]); y_s[i] = y[ord[i]]; }
    std::vector<int> grp_start, grp_size; int prev = -1;
    for(int i=0; i<n; i++) {
        int g = group_id[ord[i]];
        if (g != prev) { grp_start.push_back(i); grp_size.push_back(1); prev = g; }
        else grp_size.back()++;
    }
    GEEResult res = gee_pairs_singletons_cpp_impl(X_s, y_s, grp_start, grp_size, family, maxit, tol);
    return List::create(Named("beta")=res.beta, Named("alpha")=res.alpha, Named("vcov")=res.vcov, Named("quasi_loglik")=res.quasi_loglik, Named("converged")=res.converged, Named("niter")=res.niter);
}
