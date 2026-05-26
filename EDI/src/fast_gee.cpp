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
                                              Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                              bool smart_cold_start = true,
                                              int maxit = 100, 
                                              double tol = 1e-8,
                                              Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                              Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                              std::string optimization_alg = "irls",
                                              Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                              Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                              bool estimate_only = false);
ModelResult fast_poisson_regression_internal(const Eigen::MatrixXd& X,
                                             const Eigen::VectorXd& y,
                                             const Eigen::VectorXd& weights = Eigen::VectorXd(),
                                             Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                             bool smart_cold_start = true,
                                             int maxit = 100,
                                             double tol = 1e-8,
                                             Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                             Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                             std::string optimization_alg = "lbfgs",
                                             Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                             Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                             bool estimate_only = false);

enum class GEEFamily { GAUSSIAN, BINOMIAL, POISSON };

struct GEEResult {
    VectorXd beta;
    double alpha;
    MatrixXd vcov;
    double quasi_loglik;
    VectorXd score;
    MatrixXd bread;
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

inline double gee_obs_weight(const VectorXd& weights, int i) {
    return (weights.size() > 0) ? std::max(weights[i], 0.0) : 1.0;
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
                                              GEEFamily family,
                                              const VectorXd& weights = VectorXd()) {
    double num = 0.0, den = 0.0;
    for (int gi = 0; gi < grp_start.size(); ++gi) {
        if (grp_size[gi] < 2) continue;
        int warm_start_params = grp_start[gi];
        for (int a = 0; a < grp_size[gi]; ++a) {
            for (int b = a + 1; b < grp_size[gi]; ++b) {
                int i = warm_start_params + a;
                int j = warm_start_params + b;
                double ri = resid[i] / std::sqrt(gee_variance(mu[i], family));
                double rj = resid[j] / std::sqrt(gee_variance(mu[j], family));
                double wij = std::sqrt(gee_obs_weight(weights, i) * gee_obs_weight(weights, j));
                num += wij * ri * rj;
                den += wij;
            }
        }
    }
    if (den <= 0.0) return 0.0;
    return std::min(std::max(num / den, -0.95), 0.95);
}

inline VectorXd gee_fit_independence_glm(const MatrixXd& X, const VectorXd& y, GEEFamily family, int maxit = 100, double tol = 1e-8) {
    ModelResult fit = (family == GEEFamily::BINOMIAL) ?
        fast_logistic_regression_internal(X, y, Eigen::VectorXd(), R_NilValue, true, maxit, tol, R_NilValue, R_NilValue, "irls", R_NilValue, R_NilValue) :
        fast_poisson_regression_internal(X, y, Eigen::VectorXd(), R_NilValue, true, maxit, tol, R_NilValue, R_NilValue, "irls", R_NilValue, R_NilValue);
    if (fit.b.size() == X.cols() && fit.b.allFinite()) {
        return fit.b;
    }
    VectorXd beta = VectorXd::Zero(X.cols());
    double y_mean = std::max(y.mean(), 1e-8);
    beta[0] = (family == GEEFamily::POISSON) ? std::log(y_mean) : std::log(std::min(std::max(y_mean, 1e-6), 1.0 - 1e-6) / (1.0 - std::min(std::max(y_mean, 1e-6), 1.0 - 1e-6)));
    return beta;
}

GEEResult gee_pairs_singletons_cpp_impl(const MatrixXd& X, const VectorXd& y, 
                                        const std::vector<int>& grp_start, 
                                        const std::vector<int>& grp_size, 
                                        GEEFamily family, 
                                        const VectorXd& weights = VectorXd(),
                                        Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                        Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                        int maxit = 100, double tol = 1e-8) {
    const int n = X.rows(), p = X.cols(), G = grp_start.size();
    for (int gi = 0; gi < G; ++gi) {
        if (grp_size[gi] > 2) stop("gee_pairs_singletons_cpp: cluster %d has size %d (only singletons and pairs are supported)", gi, grp_size[gi]);
    }
    VectorXd beta = VectorXd::Zero(p);
    bool use_warm_start = false;
    if (warm_start_beta.isNotNull()) {
        NumericVector sb(warm_start_beta);
        if (sb.size() == p) {
            beta = as<VectorXd>(sb);
            use_warm_start = true;
        }
    }
    
    if (!use_warm_start) {
        if (family == GEEFamily::GAUSSIAN) {
            MatrixXd XtX = X.transpose() * X;
            VectorXd Xty = X.transpose() * y;
            beta = gee_solve_system(XtX, Xty);
        } else {
            beta = gee_fit_independence_glm(X, y, family);
        }
    }

    double alpha = 0.0; bool converged = false; int iter = 0;
    for (; iter < maxit; ++iter) {
        VectorXd mu(n), resid(n);
        for (int i = 0; i < n; ++i) { mu[i] = gee_link_inv(X.row(i).dot(beta), family); resid[i] = y[i] - mu[i]; }
        alpha = gee_estimate_exchangeable_alpha(resid, mu, grp_start, grp_size, family, weights);

        MatrixXd Bread = MatrixXd::Zero(p, p); VectorXd Score = VectorXd::Zero(p);
        
        if (iter == 0 && warm_start_fisher_info.isNotNull()) {
            NumericMatrix wf(warm_start_fisher_info);
            if (wf.rows() == p && wf.cols() == p) {
                Bread = as<MatrixXd>(wf);
            }
        }
        
        bool bread_is_zero = Bread.isZero();

        for (int gi = 0; gi < G; ++gi) {
            int sz = grp_size[gi], warm_start_params = grp_start[gi];
            if (sz == 1) {
                double mui = mu[warm_start_params], vi = gee_variance(mui, family), di = gee_dmu_deta(mui, family);
                double vi_inv = 1.0 / vi;
                double wi = gee_obs_weight(weights, warm_start_params);
                double w = wi * di * di * vi_inv;
                if (bread_is_zero) {
                    for(int r=0; r<p; r++) {
                        for(int c=0; c<=r; c++) Bread(r, c) += w * X(warm_start_params, r) * X(warm_start_params, c);
                    }
                }
                for(int r=0; r<p; r++) Score[r] += wi * (di * vi_inv) * resid[warm_start_params] * X(warm_start_params, r);
            } else {
                int i1 = warm_start_params, i2 = warm_start_params + 1;
                double m1 = mu[i1], m2 = mu[i2], v1 = gee_variance(m1, family), v2 = gee_variance(m2, family);
                double d1 = gee_dmu_deta(m1, family), d2 = gee_dmu_deta(m2, family);
                double s1 = std::sqrt(v1), s2 = std::sqrt(v2), detV = v1 * v2 * (1.0 - alpha * alpha);
                double V11 = v2/detV, V22 = v1/detV, V12 = -(alpha*s1*s2)/detV;
                double r1 = resid[i1], r2 = resid[i2];
                double sw1 = std::sqrt(gee_obs_weight(weights, i1));
                double sw2 = std::sqrt(gee_obs_weight(weights, i2));
                for(int r=0; r<p; r++) {
                    double d1xr = sw1 * d1 * X(i1, r), d2xr = sw2 * d2 * X(i2, r);
                    if (bread_is_zero) {
                        for(int c=0; c<=r; c++) {
                            double d1xc = sw1 * d1 * X(i1, c), d2xc = sw2 * d2 * X(i2, c);
                            Bread(r, c) += V11*d1xr*d1xc + V22*d2xr*d2xc + V12*(d1xr*d2xc + d2xr*d1xc);
                        }
                    }
                    Score[r] += d1xr*(V11*(sw1 * r1) + V12*(sw2 * r2)) + d2xr*(V12*(sw1 * r1) + V22*(sw2 * r2));
                }
            }
        }
        if (bread_is_zero) {
            for(int r=0; r<p; r++) for(int c=0; c<r; c++) Bread(c, r) = Bread(r, c);
        }
        VectorXd delta = gee_solve_system(Bread, Score);
        if (!delta.allFinite()) break;
        beta += delta;
        if (delta.norm() < tol) { converged = true; break; }
        bread_is_zero = true; // Compute Bread in subsequent iterations
    }

    MatrixXd Bread = MatrixXd::Zero(p, p); 
    MatrixXd Meat = MatrixXd::Zero(p, p);
    VectorXd ScoreFinal = VectorXd::Zero(p);
    double quasi_loglik = 0.0; VectorXd mu(n), resid(n);
    for (int i = 0; i < n; ++i) { mu[i] = gee_link_inv(X.row(i).dot(beta), family); resid[i] = y[i] - mu[i]; quasi_loglik += gee_obs_weight(weights, i) * gee_quasi_loglik_contrib(y[i], mu[i], family); }
    alpha = gee_estimate_exchangeable_alpha(resid, mu, grp_start, grp_size, family, weights);

    for (int gi = 0; gi < G; ++gi) {
        int sz = grp_size[gi], gs = grp_start[gi];
        if (sz == 1) {
            double mui = mu[gs], vi = gee_variance(mui, family), di = gee_dmu_deta(mui, family);
            double wi = gee_obs_weight(weights, gs);
            double Vi_inv = 1.0 / vi, w = wi * di * di * Vi_inv;
            double si_val = wi * di * Vi_inv * resid[gs];
            for(int r=0; r<p; r++) {
                double sir = si_val * X(gs, r);
                ScoreFinal[r] += sir;
                for(int c=0; c<=r; c++) {
                    Bread(r, c) += w * X(gs, r) * X(gs, c);
                    Meat(r, c) += sir * si_val * X(gs, c);
                }
            }
        } else {
            int i1 = gs, i2 = gs + 1;
            double v1 = gee_variance(mu[i1], family), v2 = gee_variance(mu[i2], family), d1 = gee_dmu_deta(mu[i1], family), d2 = gee_dmu_deta(mu[i2], family);
            double s1 = std::sqrt(v1), s2 = std::sqrt(v2), detV = v1 * v2 * (1.0 - alpha * alpha);
            double V11 = v2/detV, V22 = v1/detV, V12 = -(alpha*s1*s2)/detV, r1 = resid[i1], r2 = resid[i2];
            double sw1 = std::sqrt(gee_obs_weight(weights, i1)), sw2 = std::sqrt(gee_obs_weight(weights, i2));
            for(int r=0; r<p; r++) {
                double d1xr = sw1 * d1 * X(i1, r), d2xr = sw2 * d2 * X(i2, r), sir = d1xr*(V11*(sw1 * r1) + V12*(sw2 * r2)) + d2xr*(V12*(sw1 * r1) + V22*(sw2 * r2));
                ScoreFinal[r] += sir;
                for(int c=0; c<=r; c++) {
                    double d1xc = sw1 * d1 * X(i1, c), d2xc = sw2 * d2 * X(i2, c), sic = d1xc*(V11*(sw1 * r1) + V12*(sw2 * r2)) + d2xc*(V12*(sw1 * r1) + V22*(sw2 * r2));
                    Bread(r, c) += V11*d1xr*d1xc + V22*d2xr*d2xc + V12*(d1xr*d2xc + d2xr*d1xc);
                    Meat(r, c) += sir * sic;
                }
            }
        }
    }
    for(int r=0; r<p; r++) for(int c=0; c<r; c++) { Bread(c, r) = Bread(r, c); Meat(c, r) = Meat(r, c); }
    MatrixXd BI = gee_inverse_system(Bread);
    return {beta, alpha, BI * Meat * BI, quasi_loglik, ScoreFinal, Bread, converged, std::min(maxit, iter + 1)};
}

// [[Rcpp::export]]
List gee_pairs_singletons_weighted_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXi& group_id, std::string family_str,
                                       const Eigen::VectorXd& weights,
                                       Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                       Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                       int maxit = 100, double tol = 1e-8) {
    if (weights.size() != y.size()) {
        stop("weights must have length equal to length(y)");
    }
    GEEFamily family = GEEFamily::GAUSSIAN;
    if (family_str == "binomial") family = GEEFamily::BINOMIAL;
    else if (family_str == "poisson") family = GEEFamily::POISSON;
    int n = y.size(); std::vector<int> ord(n); for(int i=0; i<n; i++) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](int a, int b){ return group_id[a] < group_id[b]; });
    MatrixXd X_s(n, X.cols()); VectorXd y_s(n), w_s(n);
    for(int i=0; i<n; i++) { X_s.row(i) = X.row(ord[i]); y_s[i] = y[ord[i]]; w_s[i] = weights[ord[i]]; }
    std::vector<int> grp_start, grp_size; int prev = -1;
    for(int i=0; i<n; i++) {
        int g = group_id[ord[i]];
        if (g != prev) { grp_start.push_back(i); grp_size.push_back(1); prev = g; }
        else grp_size.back()++;
    }
    GEEResult res = gee_pairs_singletons_cpp_impl(X_s, y_s, grp_start, grp_size, family, w_s, warm_start_beta, warm_start_fisher_info, maxit, tol);
    return List::create(
        Named("beta")=res.beta,
        Named("alpha")=res.alpha,
        Named("vcov")=res.vcov,
        Named("quasi_loglik")=res.quasi_loglik,
        Named("score")=res.score,
        Named("fisher_information")=res.bread,
        Named("converged")=res.converged,
        Named("niter")=res.niter
    );
}

// [[Rcpp::export]]
List gee_pairs_singletons_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXi& group_id, std::string family_str, 
                               Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                               Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                               int maxit = 100, double tol = 1e-8) {
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
    GEEResult res = gee_pairs_singletons_cpp_impl(X_s, y_s, grp_start, grp_size, family, VectorXd(), warm_start_beta, warm_start_fisher_info, maxit, tol);
    return List::create(
        Named("beta")=res.beta, 
        Named("alpha")=res.alpha, 
        Named("vcov")=res.vcov, 
        Named("quasi_loglik")=res.quasi_loglik, 
        Named("score")=res.score,
        Named("fisher_information")=res.bread,
        Named("converged")=res.converged, 
        Named("niter")=res.niter
    );
}
