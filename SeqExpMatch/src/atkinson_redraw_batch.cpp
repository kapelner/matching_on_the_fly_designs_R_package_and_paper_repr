// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace Eigen;

// Helper: Find which columns have variation in rows 0..(n-1)
static std::vector<int> which_cols_vary_subset(const MatrixXd& X, int n_rows) {
    std::vector<int> varying_cols;
    int p = X.cols();

    for (int j = 0; j < p; ++j) {
        double first_val = X(0, j);
        for (int i = 1; i < n_rows; ++i) {
            if (X(i, j) != first_val) {
                varying_cols.push_back(j);
                break;
            }
        }
    }
    return varying_cols;
}

// Helper: Extract submatrix with specific rows and columns
static MatrixXd extract_submatrix(const MatrixXd& X, int n_rows, const std::vector<int>& cols) {
    if (cols.empty() || n_rows == 0) {
        return MatrixXd(n_rows, 0);
    }
    MatrixXd result(n_rows, cols.size());
    for (int i = 0; i < n_rows; ++i) {
        for (size_t j = 0; j < cols.size(); ++j) {
            result(i, j) = X(i, cols[j]);
        }
    }
    return result;
}

// Helper: Extract row with specific columns
static VectorXd extract_row_cols(const MatrixXd& X, int row, const std::vector<int>& cols) {
    VectorXd result(cols.size());
    for (size_t j = 0; j < cols.size(); ++j) {
        result(j) = X(row, cols[j]);
    }
    return result;
}

// Helper: Find linearly independent columns via QR with column pivoting
static std::vector<int> find_independent_cols(const MatrixXd& M, double tol = 1e-12) {
    if (M.cols() == 0 || M.rows() == 0) {
        return std::vector<int>();
    }

    ColPivHouseholderQR<MatrixXd> qr(M);
    qr.setThreshold(tol);
    int rank = qr.rank();

    std::vector<int> indep_cols;
    const auto& pivot = qr.colsPermutation().indices();
    for (int i = 0; i < rank; ++i) {
        indep_cols.push_back(pivot(i));
    }
    std::sort(indep_cols.begin(), indep_cols.end());
    return indep_cols;
}

// Helper: Compute single Atkinson weight for subject t (0-indexed)
// Returns the probability used for assignment (not the actual 0/1 assignment)
static double compute_atkinson_weight(
    const VectorXd& w_prev,      // Previous weights (length t)
    const MatrixXd& X_prev,      // Previous subjects' processed covariates (t x rank)
    const VectorXd& xt,          // Current subject's processed covariates (length rank)
    int rank_prev,
    int t                        // 1-indexed subject number
) {
    int rows = w_prev.size();
    int cols = rank_prev + 2;

    if (rows == 0 || cols < 2) {
        return 0.5;  // Will use CRD
    }

    // Build design matrix [w, 1, X_prev]
    MatrixXd XprevWT(rows, cols);
    for (int i = 0; i < rows; ++i) {
        XprevWT(i, 0) = w_prev(i);
        XprevWT(i, 1) = 1.0;
        for (int j = 0; j < rank_prev; ++j) {
            XprevWT(i, j + 2) = X_prev(i, j);
        }
    }

    // Compute (X'X)^-1
    MatrixXd XwtXw = XprevWT.transpose() * XprevWT;
    FullPivLU<MatrixXd> lu(XwtXw);
    if (!lu.isInvertible()) {
        return 0.5;  // Will use CRD
    }

    MatrixXd M = static_cast<double>(t - 1) * lu.inverse();

    // Extract row segment for computation
    VectorXd row_segment(rank_prev + 1);
    for (int j = 0; j < rank_prev + 1; ++j) {
        row_segment(j) = M(0, j + 1);
    }

    // Build xt vector with intercept
    VectorXd xt_full(rank_prev + 1);
    xt_full(0) = 1.0;
    for (int j = 0; j < rank_prev; ++j) {
        xt_full(j + 1) = xt(j);
    }

    double A = row_segment.dot(xt_full);
    if (A == 0 || !std::isfinite(A)) {
        return 0.5;  // Will use CRD
    }

    double val = M(0, 0) / A + 1.0;
    double s_over_A_plus_one_sq = val * val;
    double prob = s_over_A_plus_one_sq / (s_over_A_plus_one_sq + 1.0);
    prob = std::max(0.0, std::min(1.0, prob));

    return prob;
}

// [[Rcpp::export]]
NumericVector atkinson_redraw_batch_cpp(
    const Eigen::MatrixXd& X,     // Full covariate matrix (n x p), already numeric
    int n,                         // Number of subjects
    int p_raw,                     // Number of raw covariates (for early-subject CRD threshold)
    double prob_T = 0.5           // Treatment probability for CRD
) {
    NumericVector w(n);

    // Threshold for using CRD (early subjects)
    int crd_threshold = p_raw + 2 + 1;

    for (int t = 1; t <= n; ++t) {
        // For early subjects, use CRD
        if (t <= crd_threshold) {
            w[t - 1] = R::rbinom(1, prob_T);
            continue;
        }

        // Find which columns vary in rows 0..(t-1)
        std::vector<int> var_cols = which_cols_vary_subset(X, t);

        if (var_cols.empty()) {
            w[t - 1] = R::rbinom(1, prob_T);
            continue;
        }

        // Extract submatrix for rows 0..(t-2) with varying columns
        MatrixXd X_var = extract_submatrix(X, t - 1, var_cols);

        // Find linearly independent columns
        std::vector<int> indep_cols = find_independent_cols(X_var);

        if (indep_cols.empty()) {
            w[t - 1] = R::rbinom(1, prob_T);
            continue;
        }

        // Extract final processed matrix and vector
        MatrixXd X_prev(t - 1, indep_cols.size());
        for (int i = 0; i < t - 1; ++i) {
            for (size_t j = 0; j < indep_cols.size(); ++j) {
                X_prev(i, j) = X_var(i, indep_cols[j]);
            }
        }

        // Get xt_prev (current subject's row with same column processing)
        VectorXd xt_var = extract_row_cols(X, t - 1, var_cols);
        VectorXd xt_prev(indep_cols.size());
        for (size_t j = 0; j < indep_cols.size(); ++j) {
            xt_prev(j) = xt_var(indep_cols[j]);
        }

        // Get previous weights
        VectorXd w_prev(t - 1);
        for (int i = 0; i < t - 1; ++i) {
            w_prev(i) = w[i];
        }

        int rank_prev = X_prev.cols();

        // Compute Atkinson probability
        double prob = compute_atkinson_weight(w_prev, X_prev, xt_prev, rank_prev, t);

        // Draw assignment
        w[t - 1] = R::rbinom(1, prob);
    }

    return w;
}
