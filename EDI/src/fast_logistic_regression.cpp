#include "_helper_functions.h"
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace {

inline double plogis_manual(double x) {
    if (x > 20.0) return 1.0;
    if (x < -20.0) return 0.0;
    return 1.0 / (1.0 + std::exp(-x));
}

// Manual LLT (Cholesky) decomposition and solver to avoid Eigen templates
// A is p*p symmetric, positive definite. b is p*1.
bool solve_llt_raw(double* x, const double* A, const double* b, int p) {
    std::vector<double> L(p * p, 0.0);
    for (int i = 0; i < p; i++) {
        for (int j = 0; j <= i; j++) {
            double s = 0;
            for (int k = 0; k < j; k++) s += L[i * p + k] * L[j * p + k];
            if (i == j) {
                double val = A[i * p + i] - s;
                if (val <= 1e-15) return false;
                L[i * p + j] = std::sqrt(val);
            } else {
                L[i * p + j] = (A[i * p + j] - s) / L[j * p + j];
            }
        }
    }
    std::vector<double> y_tmp(p);
    for (int i = 0; i < p; i++) {
        double s = 0;
        for (int k = 0; k < i; k++) s += L[i * p + k] * y_tmp[k];
        y_tmp[i] = (b[i] - s) / L[i * p + i];
    }
    for (int i = p - 1; i >= 0; i--) {
        double s = 0;
        for (int k = i + 1; k < p; k++) s += L[k * p + i] * x[k];
        x[i] = (y_tmp[i] - s) / L[i * p + i];
    }
    return true;
}

} // namespace

// Internal pure C++ logic
ModelResult fast_logistic_regression_internal(const Eigen::MatrixXd& X_eigen, 
                                              const Eigen::VectorXd& y_eigen, 
                                              const Eigen::VectorXd& weights_eigen = Eigen::VectorXd(),
                                              int maxit = 100, 
                                              double tol = 1e-8,
                                              Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                              Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                              std::string optimization_alg = "irls") {
    int n = X_eigen.rows();
    int p = X_eigen.cols();
    bool use_weights = (weights_eigen.size() == n);
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);

    int p_free = fixed_spec.free_idx.size();
    std::vector<double> beta_full(p, 0.0);
    for(size_t j=0; j<fixed_spec.fixed_idx.size(); j++) beta_full[fixed_spec.fixed_idx[j]] = fixed_spec.fixed_values[j];
    
    std::vector<double> beta_free(p_free, 0.0);
    for(int j=0; j<p_free; j++) beta_free[j] = beta_full[fixed_spec.free_idx[j]];

    std::vector<double> eta_fixed(n, 0.0);
    for(size_t k=0; k<fixed_spec.fixed_idx.size(); k++) {
        int idx = fixed_spec.fixed_idx[k];
        double val = fixed_spec.fixed_values[k];
        for(int i=0; i<n; i++) eta_fixed[i] += X_eigen(i, idx) * val;
    }

    std::vector<double> X_free_raw(n * p_free);
    for (int j = 0; j < p_free; ++j) {
        int col_idx = fixed_spec.free_idx[j];
        for(int i=0; i<n; i++) X_free_raw[i * p_free + j] = X_eigen(i, col_idx);
    }

    std::vector<double> mu(n);
    std::vector<double> w_diag(n);
    bool converged = false;

    int n_threads = 1;
    #ifdef _OPENMP
    n_threads = omp_get_max_threads();
    #endif

    std::vector<double> score_free(p_free);
    std::vector<double> final_XtWX(p_free * p_free);

    // Pre-allocate thread-local storage once
    std::vector<double> XtWX_threads(n_threads * p_free * p_free, 0.0);
    std::vector<double> score_threads(n_threads * p_free, 0.0);

    for (int iter = 0; iter < maxit; iter++) {
        #pragma omp parallel for
        for(int i=0; i<n; i++) {
            double eta_i = eta_fixed[i];
            const double* x_ptr = &X_free_raw[i * p_free];
            for(int j=0; j<p_free; j++) eta_i += x_ptr[j] * beta_free[j];
            mu[i] = plogis_manual(eta_i);
            double base_w = std::max(mu[i] * (1.0 - mu[i]), 1e-10);
            w_diag[i] = use_weights ? weights_eigen[i] * base_w : base_w;
        }

        std::fill(XtWX_threads.begin(), XtWX_threads.end(), 0.0);
        std::fill(score_threads.begin(), score_threads.end(), 0.0);

        #pragma omp parallel
        {
            int tid = 0;
            #ifdef _OPENMP
            tid = omp_get_thread_num();
            #endif
            double* t_xtwx = &XtWX_threads[tid * p_free * p_free];
            double* t_score = &score_threads[tid * p_free];

            #pragma omp for
            for(int i=0; i<n; i++) {
                double diff = y_eigen[i] - mu[i];
                if (use_weights) diff *= weights_eigen[i];
                double wi = w_diag[i];
                const double* x_ptr = &X_free_raw[i * p_free];
                for(int r=0; r<p_free; r++) {
                    double x_ir = x_ptr[r];
                    t_score[r] += x_ir * diff;
                    double xirwi = x_ir * wi;
                    for(int c=0; c<=r; c++) {
                        t_xtwx[r * p_free + c] += xirwi * x_ptr[c];
                    }
                }
            }
        }

        std::fill(final_XtWX.begin(), final_XtWX.end(), 0.0);
        std::fill(score_free.begin(), score_free.end(), 0.0);
        for (int t = 0; t < n_threads; ++t) {
            for(int r=0; r<p_free; r++) {
                score_free[r] += score_threads[t * p_free + r];
                for(int c=0; c<=r; c++) final_XtWX[r * p_free + c] += XtWX_threads[t * p_free * p_free + r * p_free + c];
            }
        }
        for (int r = 0; r < p_free; ++r) {
            for (int c = r + 1; c < p_free; ++c) final_XtWX[r * p_free + c] = final_XtWX[c * p_free + r];
        }

        std::vector<double> delta(p_free);
        if (!solve_llt_raw(delta.data(), final_XtWX.data(), score_free.data(), p_free)) break;
        
        double norm_delta_sq = 0.0;
        for(int j=0; j<p_free; j++) {
            beta_free[j] += delta[j];
            norm_delta_sq += delta[j] * delta[j];
        }
        if (std::sqrt(norm_delta_sq) < tol) { converged = true; break; }
    }

    ModelResult res;
    res.b.resize(p);
    for(int j=0; j<p; j++) res.b[j] = beta_full[j];
    for(int j=0; j<p_free; j++) res.b[fixed_spec.free_idx[j]] = beta_free[j];
    res.mu.resize(n);
    for(int i=0; i<n; i++) res.mu[i] = mu[i];
    
    Eigen::MatrixXd info_free_eigen = Eigen::MatrixXd::Zero(p_free, p_free);
    for(int r=0; r<p_free; r++) {
        for(int c=0; c<p_free; c++) info_free_eigen(r, c) = final_XtWX[r * p_free + c];
    }
    res.XtWX = expand_free_covariance(p, fixed_spec, info_free_eigen, false);
    res.converged = converged;
    return res;
}

// [[Rcpp::export]]
Eigen::VectorXd get_logistic_regression_score_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& beta) {
    int n = X.rows();
    int p = X.cols();
    Eigen::VectorXd res = Eigen::VectorXd::Zero(p);
    for(int i=0; i<n; i++) {
        double eta_i = 0.0;
        for(int j=0; j<p; j++) eta_i += X(i, j) * beta[j];
        double mu_i = plogis_manual(eta_i);
        double diff = y[i] - mu_i;
        for(int j=0; j<p; j++) res[j] += X(i, j) * diff;
    }
    return res;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_logistic_regression_hessian_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) {
    int n = X.rows();
    int p = X.cols();
    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(p, p);
    for(int i=0; i<n; i++) {
        double eta_i = 0.0;
        for(int j=0; j<p; j++) eta_i += X(i, j) * beta[j];
        double mu_i = plogis_manual(eta_i);
        double wi = mu_i * (1.0 - mu_i);
        for(int r=0; r<p; r++) {
            double xirwi = X(i, r) * wi;
            for(int c=0; c<=r; c++) res(r, c) += xirwi * X(i, c);
        }
    }
    for (int r = 0; r < p; ++r) {
        for (int c = r + 1; c < p; ++c) res(r, c) = res(c, r);
    }
    return -res;
}

// [[Rcpp::export]]
Eigen::VectorXd get_logistic_regression_weighted_score_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& weights, const Eigen::VectorXd& beta) {
    int n = X.rows();
    int p = X.cols();
    Eigen::VectorXd res = Eigen::VectorXd::Zero(p);
    for(int i=0; i<n; i++) {
        double eta_i = 0.0;
        for(int j=0; j<p; j++) eta_i += X(i, j) * beta[j];
        double mu_i = plogis_manual(eta_i);
        double diff = weights[i] * (y[i] - mu_i);
        for(int j=0; j<p; j++) res[j] += X(i, j) * diff;
    }
    return res;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_logistic_regression_weighted_hessian_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& weights, const Eigen::VectorXd& beta) {
    int n = X.rows();
    int p = X.cols();
    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(p, p);
    for(int i=0; i<n; i++) {
        double eta_i = 0.0;
        for(int j=0; j<p; j++) eta_i += X(i, j) * beta[j];
        double mu_i = plogis_manual(eta_i);
        double wi = weights[i] * mu_i * (1.0 - mu_i);
        for(int r=0; r<p; r++) {
            double xirwi = X(i, r) * wi;
            for(int c=0; c<=r; c++) res(r, c) += xirwi * X(i, c);
        }
    }
    for (int r = 0; r < p; ++r) {
        for (int c = r + 1; c < p; ++c) res(r, c) = res(c, r);
    }
    return -res;
}

// [[Rcpp::export]]
List fast_logistic_regression_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int maxit = 100, double tol = 1e-8,
                                  Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                  Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                  std::string optimization_alg = "irls") {
    ModelResult res = fast_logistic_regression_internal(X, y, Eigen::VectorXd(), maxit, tol, fixed_idx, fixed_values, optimization_alg);
    Eigen::VectorXd weights_vec(X.rows());
    for(int i=0; i<X.rows(); i++) weights_vec[i] = res.mu[i] * (1.0 - res.mu[i]);
    return List::create(
        Named("b") = res.b,
        Named("w") = weights_vec
    );
}

// [[Rcpp::export]]
List fast_logistic_regression_weighted_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& weights, int maxit = 100, double tol = 1e-8,
                                           Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                           Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                           std::string optimization_alg = "irls") {
    ModelResult res = fast_logistic_regression_internal(X, y, weights, maxit, tol, fixed_idx, fixed_values, optimization_alg);
    return List::create(
        Named("b") = res.b,
        Named("mu") = res.mu,
        Named("XtWX") = res.XtWX,
        Named("converged") = res.converged
    );
}

// [[Rcpp::export]]
List fast_logistic_regression_with_var_cpp(const Eigen::MatrixXd& Xmm, const Eigen::VectorXd& y, int j = 2,
                                           Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                           Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                           std::string optimization_alg = "irls") {
    ModelResult res = fast_logistic_regression_internal(Xmm, y, Eigen::VectorXd(), 100, 1e-8, fixed_idx, fixed_values, optimization_alg);
    FixedParamSpec fixed_spec = make_fixed_param_spec(Xmm.cols(), fixed_idx, fixed_values);
    
    int p_free = fixed_spec.free_idx.size();
    Eigen::MatrixXd info_free = subset_matrix(res.XtWX, fixed_spec.free_idx, fixed_spec.free_idx);
    std::vector<double> info_free_v(p_free * p_free);
    for(int r=0; r<p_free; r++) for(int c=0; c<p_free; c++) info_free_v[r * p_free + c] = info_free(r, c);

    Eigen::MatrixXd cov_free = Eigen::MatrixXd::Zero(p_free, p_free);
    for(int col=0; col<p_free; col++) {
        std::vector<double> b(p_free, 0.0);
        b[col] = 1.0;
        std::vector<double> x_sol(p_free);
        solve_llt_raw(x_sol.data(), info_free_v.data(), b.data(), p_free);
        for(int r=0; r<p_free; r++) cov_free(r, col) = x_sol[r];
    }
    
    Eigen::MatrixXd vcov = expand_free_covariance(Xmm.cols(), fixed_spec, cov_free, true);
    res.ssq_b_j = (j > 0 && j <= Xmm.cols()) ? vcov(j - 1, j - 1) : NA_REAL;
    res.ssq_b_2 = (Xmm.cols() >= 2) ? vcov(1, 1) : NA_REAL;
    Eigen::VectorXd score = get_logistic_regression_score_cpp(Xmm, y, res.b);
    Eigen::MatrixXd information = res.XtWX;
    Eigen::VectorXd eta = Xmm * res.b;
    double neg_loglik = 0.0;
    for (int i = 0; i < eta.size(); ++i) {
        double log_one_plus_exp_eta = (eta[i] > 0.0) ?
            eta[i] + std::log1p(std::exp(-eta[i])) :
            std::log1p(std::exp(eta[i]));
        neg_loglik += log_one_plus_exp_eta - y[i] * eta[i];
    }

    return List::create(
        Named("b") = res.b,
        Named("params") = res.b,
        Named("ssq_b_j") = res.ssq_b_j,
        Named("ssq_b_2") = res.ssq_b_2,
        Named("vcov") = vcov,
        Named("score") = score,
        Named("observed_information") = information,
        Named("fisher_information") = information,
        Named("information") = information,
        Named("information_type") = "fisher",
        Named("hessian") = -information,
        Named("neg_loglik") = neg_loglik,
        Named("neg_ll") = neg_loglik,
        Named("loglik") = R_finite(neg_loglik) ? -neg_loglik : NA_REAL,
        Named("converged") = res.converged
    );
}
