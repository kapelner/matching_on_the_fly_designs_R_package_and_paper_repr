#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

List compute_zhang_match_data_cpp(const NumericMatrix& X,
                                  const NumericVector& y,
                                  const IntegerVector& w,
                                  const IntegerVector& m_vec);

// [[Rcpp::export]]
List match_diffs_cpp(SEXP X_sexp,
                        SEXP y_sexp,
                        SEXP w_sexp,
                                                SEXP m_vec_sexp,
                                                int m) {
        NumericMatrix X_mat(X_sexp);
        NumericVector y_vec(y_sexp);
        IntegerVector w_int_vec(w_sexp);
        IntegerVector m_int_vec(m_vec_sexp);
        Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
        Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
        Eigen::Map<const Eigen::VectorXi> w(w_int_vec.begin(), w_int_vec.size());
        Eigen::Map<const Eigen::VectorXi> m_vec(m_int_vec.begin(), m_int_vec.size());

        List match_data = compute_zhang_match_data_cpp(
                wrap(X),
                wrap(y),
                wrap(w),
                wrap(m_vec)
        );

        return List::create(
        _["yTs_matched"] =              match_data["yTs_matched"],
        _["yCs_matched"] =              match_data["yCs_matched"],
        _["X_matched_diffs"] =  match_data["X_matched_diffs"]
        );
}

