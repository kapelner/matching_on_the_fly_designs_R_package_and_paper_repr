#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
List match_diffs_cpp(const Eigen::VectorXi& w, 
                 		const Eigen::VectorXi& match_indic, 
                 		const Eigen::VectorXd& y,
                 		const Eigen::MatrixXd& X,
                 		int m) {
  
  Eigen::VectorXd yTs_matched(m);
  Eigen::VectorXd yCs_matched(m);
  Eigen::MatrixXd X_matched_diffs(m, X.cols());
  
  for (int match_id = 1; match_id <= m; match_id++) {
    // placeholders
    double yT = NA_REAL, yC = NA_REAL;
    Eigen::RowVectorXd xmT(X.cols()), xmC(X.cols());
    xmT.setZero();
    xmC.setZero();
    bool foundT = false, foundC = false;
    
    for (int i = 0; i < w.size(); i++) {
      if (match_indic[i] == match_id) {
        if (w[i] == 1) {
          yT = y[i];
          xmT = X.row(i);
          foundT = true;
        } else {
          yC = y[i];
          xmC = X.row(i);
          foundC = true;
        }
        if (foundT && foundC) break;
      }
    }
    yTs_matched[match_id - 1] = yT;
    yCs_matched[match_id - 1] = yC;
    X_matched_diffs.row(match_id - 1) = xmT - xmC;
  }
  
  return List::create(
    _["yTs_matched"] = 		yTs_matched,
    _["yCs_matched"] = 		yCs_matched,
    _["X_matched_diffs"] = 	X_matched_diffs
  );
}
