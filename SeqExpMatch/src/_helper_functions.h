#include <RcppEigen.h>

typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;

double eigen_compute_single_entry_of_diagonal_matrix_cpp(MapMat M, int j);
Eigen::MatrixXd eigen_Xt_times_diag_w_times_X_cpp(MapMat X, MapVec w);