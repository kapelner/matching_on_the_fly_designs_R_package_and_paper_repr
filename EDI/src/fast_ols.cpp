#include "_helper_functions.h"
#include <RcppEigen.h>

using namespace Rcpp;

// Internal pure C++ logic
ModelResult fast_ols_internal(const Eigen::MatrixXd& X,
                              const Eigen::VectorXd& y,
                              Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                              Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                              bool estimate_only = false) {
    const int n = X.rows();
    const int p = X.cols();
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);
    const int p_free = fixed_spec.free_idx.size();
    
    Eigen::VectorXd y_adj = y;
    for (int j = 0; j < fixed_spec.fixed_idx.size(); ++j) {
        y_adj.noalias() -= X.col(fixed_spec.fixed_idx[j]) * fixed_spec.fixed_values[j];
    }

    ModelResult res;
    res.b = Eigen::VectorXd::Zero(p);
    for (int j = 0; j < fixed_spec.fixed_idx.size(); ++j) res.b[fixed_spec.fixed_idx[j]] = fixed_spec.fixed_values[j];

    if (p_free > 0) {
        // Build XtX and Xty directly to avoid X_free allocation
        Eigen::MatrixXd XtX_free = Eigen::MatrixXd::Zero(p_free, p_free);
        Eigen::VectorXd Xty_free = Eigen::VectorXd::Zero(p_free);
        
        for (int j = 0; j < p_free; ++j) {
            int col_j = fixed_spec.free_idx[j];
            Xty_free[j] = X.col(col_j).dot(y_adj);
            for (int i = 0; i <= j; ++i) {
                int col_i = fixed_spec.free_idx[i];
                double val = X.col(col_j).dot(X.col(col_i));
                XtX_free(i, j) = val;
                if (i != j) XtX_free(j, i) = val;
            }
        }

        Eigen::LDLT<Eigen::MatrixXd> ldlt(XtX_free);
        Eigen::VectorXd beta_free;
        if (ldlt.info() == Eigen::Success) {
            beta_free = ldlt.solve(Xty_free);
        } else {
            // Fall back to COD if XtX is not well-behaved
            Eigen::MatrixXd X_free(n, p_free);
            for (int j = 0; j < p_free; ++j) X_free.col(j) = X.col(fixed_spec.free_idx[j]);
            Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(X_free);
            beta_free = cod.solve(y_adj);
        }
        
        for (int j = 0; j < p_free; ++j) res.b[fixed_spec.free_idx[j]] = beta_free[j];
        if (!estimate_only) res.XtWX = expand_free_covariance(p, fixed_spec, XtX_free, false);
    } else {
        if (!estimate_only) res.XtWX = Eigen::MatrixXd::Zero(p, p);
    }

    if (!res.b.allFinite()) {
        res.b = Eigen::VectorXd::Constant(p, NA_REAL);
    }
    return res;
}

//' @title Fast Ordinary Least Squares (C++)
//' @description High-performance OLS fitting using Eigen's Complete Orthogonal Decomposition.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @return A list containing coefficients and the XtX matrix.
//' @export
//' @keywords internal
//' @examples
//' X = matrix(rnorm(100), 10, 10)
//' y = rnorm(10)
//' fast_ols_cpp(X, y)
// [[Rcpp::export]]
List fast_ols_cpp(const Eigen::MatrixXd& X,
                  const Eigen::VectorXd& y,
                  Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                  Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue) {
    ModelResult res = fast_ols_internal(X, y, fixed_idx, fixed_values, true);
    return List::create(
        Named("b") = res.b
    );
}

//' @title Fast OLS with Variance (C++)
//' @description OLS fitting with variance-covariance matrix and error variance estimation.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param j 1-based index of the parameter for which to return specific variance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @return A list containing coefficients, vcov, ssq_b_j, and sigma2_hat.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_ols_with_var_cpp(SEXP X_sexp, SEXP y_sexp, int j = 2,
                           Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                           Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue) {
    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());

    const int n = X.rows();
    const int p = X.cols();
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);
    const int p_free = fixed_spec.free_idx.size();

    Eigen::VectorXd y_adj = y;
    for (int k = 0; k < fixed_spec.fixed_idx.size(); ++k) {
        y_adj.noalias() -= X.col(fixed_spec.fixed_idx[k]) * fixed_spec.fixed_values[k];
    }

    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    for (int k = 0; k < fixed_spec.fixed_idx.size(); ++k) beta[fixed_spec.fixed_idx[k]] = fixed_spec.fixed_values[k];

    Eigen::MatrixXd XtX_free(p_free, p_free);
    Eigen::VectorXd Xty_free(p_free);
    Eigen::MatrixXd X_free(n, p_free);

    for (int k = 0; k < p_free; ++k) X_free.col(k) = X.col(fixed_spec.free_idx[k]);
    
    XtX_free.setZero().selfadjointView<Eigen::Lower>().rankUpdate(X_free.transpose());
    XtX_free.triangularView<Eigen::Upper>() = XtX_free.transpose();
    Xty_free.noalias() = X_free.transpose() * y_adj;

    Eigen::LDLT<Eigen::MatrixXd> ldlt(XtX_free);
    bool converged = (ldlt.info() == Eigen::Success);
    Eigen::VectorXd beta_free;
    
    if (converged) {
        beta_free = ldlt.solve(Xty_free);
    } else {
        Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(X_free);
        beta_free = cod.solve(y_adj);
    }

    for (int k = 0; k < p_free; ++k) beta[fixed_spec.free_idx[k]] = beta_free[k];

    Eigen::VectorXd e = y_adj - X_free * beta_free;
    double sse = e.squaredNorm();
    double sigma2_hat = sse / (n - p_free);

    auto free_idx_of = [&](int k) -> int {
        for (int jj = 0; jj < p_free; ++jj)
            if (fixed_spec.free_idx[jj] == k) return jj;
        return -1;
    };

    double ssq_b_j = NA_REAL;
    int f_j = (j > 0 && j <= p) ? free_idx_of(j - 1) : -1;
    if (f_j >= 0 && converged) {
        Eigen::VectorXd unit = Eigen::VectorXd::Unit(p_free, f_j);
        ssq_b_j = sigma2_hat * ldlt.solve(unit)(f_j);
    }

    double ssq_b_2 = NA_REAL;
    int f_2 = (p >= 2) ? free_idx_of(1) : -1;
    if (f_2 >= 0 && converged) {
        if (f_2 == f_j) {
            ssq_b_2 = ssq_b_j;
        } else {
            Eigen::VectorXd unit = Eigen::VectorXd::Unit(p_free, f_2);
            ssq_b_2 = sigma2_hat * ldlt.solve(unit)(f_2);
        }
    }

    return List::create(
        Named("b") = beta,
        Named("XtX") = expand_free_covariance(p, fixed_spec, XtX_free, false),
        Named("ssq_b_j") = ssq_b_j,
        Named("ssq_b_2") = ssq_b_2,
        Named("sigma2_hat") = sigma2_hat,
        Named("converged") = converged
    );
}
