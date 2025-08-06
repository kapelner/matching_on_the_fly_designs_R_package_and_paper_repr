#include <RcppEigen.h>

//typedef Eigen::Map<Eigen::MatrixXd> MapMat;
// Eigen::Map<Eigen::VectorXd> MapVec;

double eigen_compute_single_entry_of_diagonal_matrix_cpp(Eigen::Map<Eigen::MatrixXd> M, int j);
Eigen::MatrixXd eigen_Xt_times_diag_w_times_X_cpp(Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::VectorXd> w);