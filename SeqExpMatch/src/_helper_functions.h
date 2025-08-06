#include <RcppEigen.h>
using namespace Rcpp;

typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;

double eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(Eigen::MatrixXd M, int j);
Eigen::MatrixXd eigen_Xt_times_diag_w_times_X_cpp(Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::VectorXd> w);
List fast_ols_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y);