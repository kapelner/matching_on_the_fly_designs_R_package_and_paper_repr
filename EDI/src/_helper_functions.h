#ifndef EDI_HELPERS_H
#define EDI_HELPERS_H

#include <RcppEigen.h>
#include <vector>

// Pure C++ result structure to avoid R List contention
struct ModelResult {
    Eigen::VectorXd b;
    Eigen::VectorXd mu;
    Eigen::MatrixXd XtWX;
    double ssq_b_j;
    double ssq_b_2;
    double dispersion;
    double sigma2_hat;
    bool converged;

    ModelResult() : ssq_b_j(NA_REAL), ssq_b_2(NA_REAL), dispersion(NA_REAL), sigma2_hat(NA_REAL), converged(false) {}
};

// Pure C++ internal helpers
double compute_diagonal_inverse_entry(const Eigen::MatrixXd& M, int j);

// R-facing exports
double eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(Eigen::MatrixXd M, int j);
Eigen::MatrixXd eigen_Xt_times_diag_w_times_X_cpp(Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::VectorXd> w);

#endif
