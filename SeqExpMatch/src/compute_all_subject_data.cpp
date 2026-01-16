// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// Helper: Find which columns have variation (not all same value)
static std::vector<int> which_cols_vary(const MatrixXd& X) {
    std::vector<int> varying_cols;
    int n = X.rows();
    int p = X.cols();

    for (int j = 0; j < p; ++j) {
        double first_val = X(0, j);
        for (int i = 1; i < n; ++i) {
            if (X(i, j) != first_val) {
                varying_cols.push_back(j);
                break;
            }
        }
    }
    return varying_cols;
}

// Helper: Extract columns by index
static MatrixXd extract_cols(const MatrixXd& X, const std::vector<int>& cols) {
    if (cols.empty()) {
        return MatrixXd(X.rows(), 0);
    }
    MatrixXd result(X.rows(), cols.size());
    for (size_t j = 0; j < cols.size(); ++j) {
        result.col(j) = X.col(cols[j]);
    }
    return result;
}

// Helper: Extract row as vector
static VectorXd extract_row_cols(const MatrixXd& X, int row, const std::vector<int>& cols) {
    VectorXd result(cols.size());
    for (size_t j = 0; j < cols.size(); ++j) {
        result(j) = X(row, cols[j]);
    }
    return result;
}

// Helper: Scale columns (z-score standardization)
static MatrixXd scale_columns(const MatrixXd& X) {
    int n = X.rows();
    int p = X.cols();
    MatrixXd X_scaled(n, p);

    for (int j = 0; j < p; ++j) {
        double mean = X.col(j).mean();
        double var = (X.col(j).array() - mean).square().sum() / (n - 1);
        double sd = std::sqrt(var);

        if (sd > 0) {
            X_scaled.col(j) = (X.col(j).array() - mean) / sd;
        } else {
            // Constant column: set to 0 after scaling
            X_scaled.col(j).setZero();
        }
    }
    return X_scaled;
}

// Helper: Drop linearly dependent columns using QR with column pivoting
// Returns indices of linearly independent columns
static std::vector<int> find_independent_cols(const MatrixXd& M, double tol = 1e-12) {
    if (M.cols() == 0) {
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
    // Sort to maintain original column order
    std::sort(indep_cols.begin(), indep_cols.end());
    return indep_cols;
}

// Helper: Compute matrix rank
static int matrix_rank(const MatrixXd& M, double tol = 1e-12) {
    if (M.rows() == 0 || M.cols() == 0) {
        return 0;
    }
    ColPivHouseholderQR<MatrixXd> qr(M);
    qr.setThreshold(tol);
    return qr.rank();
}

// Helper struct to hold result of processing a subset
struct SubsetResult {
    MatrixXd Xint;
    VectorXd xt;
    int rank;
    std::vector<int> col_indices; // 1-based indices for R
};

// Process a subset of rows with optional scaling
static SubsetResult process_subset(
    const MatrixXd& X,           // Full matrix (only first t rows are valid)
    int t,                       // Current subject (1-indexed, so row t-1 is current)
    const std::vector<int>& row_indices, // 0-indexed row indices for this subset
    bool scaled,
    double tol = 1e-12
) {
    SubsetResult result;
    int p = X.cols();

    if (row_indices.empty()) {
        result.Xint = MatrixXd(0, p);
        result.rank = 0;
        result.xt = X.row(t - 1).transpose();
        // Return all column indices 1..p as technically none were dropped due to dependence (just no data)
        // But logic below suggests we should probably return empty or all?
        // Existing R code implies we might want all cols if we have no rows?
        // Let's stick to what valid columns are there.
        for(int j=0; j<p; ++j) result.col_indices.push_back(j + 1);
        return result;
    }

    int n_subset = row_indices.size();

    // Extract subset rows
    MatrixXd Xsub(n_subset, p);
    for (int i = 0; i < n_subset; ++i) {
        Xsub.row(i) = X.row(row_indices[i]);
    }

    // Find columns with variation in this subset
    std::vector<int> var_cols = which_cols_vary(Xsub);

    if (var_cols.empty()) {
        // No varying columns
        result.Xint = MatrixXd(n_subset, 0);
        result.rank = 0;
        result.xt = VectorXd(0);
        return result;
    }

    // Extract only varying columns
    MatrixXd Xvar = extract_cols(Xsub, var_cols);
    VectorXd xt_var = extract_row_cols(X, t - 1, var_cols);
    
    std::vector<int> indep_cols;

    if (scaled && n_subset > 1) {
        // Append current subject's row for scaling
        MatrixXd Xvar_with_t(n_subset + 1, var_cols.size());
        Xvar_with_t.topRows(n_subset) = Xvar;
        Xvar_with_t.row(n_subset) = xt_var.transpose();

        // Scale all together
        Xvar_with_t = scale_columns(Xvar_with_t);

        // Handle NaN from constant columns (set to 0)
        Xvar_with_t = Xvar_with_t.unaryExpr([](double v) {
            return std::isnan(v) ? 0.0 : v;
        });

        // Find independent columns after scaling
        indep_cols = find_independent_cols(Xvar_with_t, tol);

        if (indep_cols.empty()) {
            result.Xint = MatrixXd(n_subset, 0);
            result.rank = 0;
            result.xt = VectorXd(0);
            return result;
        }

        // Extract independent columns
        MatrixXd Xfinal = extract_cols(Xvar_with_t, indep_cols);
        result.rank = matrix_rank(Xfinal, tol);

        // Split back: rows 0..n_subset-1 are Xint, row n_subset is xt
        result.Xint = Xfinal.topRows(n_subset);
        result.xt = Xfinal.row(n_subset).transpose();
    } else {
        // Not scaled: find independent columns directly
        indep_cols = find_independent_cols(Xvar, tol);

        if (indep_cols.empty()) {
            result.Xint = MatrixXd(n_subset, 0);
            result.rank = 0;
            result.xt = VectorXd(0);
            return result;
        }

        result.Xint = extract_cols(Xvar, indep_cols);
        result.rank = matrix_rank(result.Xint, tol);

        // xt uses the same column indices
        VectorXd xt_indep(indep_cols.size());
        for (size_t j = 0; j < indep_cols.size(); ++j) {
            xt_indep(j) = xt_var(indep_cols[j]);
        }
        result.xt = xt_indep;
    }

    // Map back to original column indices (1-based)
    for (int idx : indep_cols) {
        result.col_indices.push_back(var_cols[idx] + 1);
    }

    return result;
}

// [[Rcpp::export]]
List compute_all_subject_data_cpp(
    const Eigen::MatrixXd& X,           // Full covariate matrix (n x p), but only first t rows are valid
    int t,                               // Current subject index (1-indexed)
    const IntegerVector& i_all_y_present_R, // R indices (1-indexed) where y is present
    double rank_tol = 1e-12
) {
    // Convert R indices (1-indexed) to C++ indices (0-indexed)
    std::vector<int> i_past;          // indices 0 to t-2 (subjects 1 to t-1)
    std::vector<int> i_all;           // indices 0 to t-1 (subjects 1 to t)
    std::vector<int> i_all_y_present; // indices where y is present

    // Build i_past and i_all
    for (int i = 0; i < t; ++i) {
        i_all.push_back(i);
        if (i < t - 1) {
            i_past.push_back(i);
        }
    }

    // Convert R indices to 0-indexed
    for (int k = 0; k < i_all_y_present_R.size(); ++k) {
        i_all_y_present.push_back(i_all_y_present_R[k] - 1);
    }

    // Handle t == 1 case specially (first subject)
    if (t == 1) {
        VectorXd xt = X.row(0).transpose();
        MatrixXd Xall = X.row(0);

        // Check if current subject has y present
        bool has_y = std::find(i_all_y_present.begin(), i_all_y_present.end(), 0) != i_all_y_present.end();

        MatrixXd X_all_with_y_scaled;
        int rank_all_with_y_scaled;
        IntegerVector all_col_indices(X.cols());
        for(int j=0; j<X.cols(); ++j) all_col_indices[j] = j + 1;

        if (has_y) {
            X_all_with_y_scaled = Xall;
            rank_all_with_y_scaled = xt.size();
        } else {
            X_all_with_y_scaled = MatrixXd(0, X.cols());
            rank_all_with_y_scaled = 0;
        }

        return List::create(
            Named("X_prev") = MatrixXd(0, X.cols()),
            Named("rank_prev") = 0,
            Named("xt_prev") = VectorXd(0),
            Named("cols_prev") = all_col_indices,
            Named("X_all") = Xall,
            Named("cols_all") = all_col_indices,
            Named("X_all_scaled") = Xall,
            Named("xt_all_scaled") = xt,
            Named("cols_all_scaled") = all_col_indices,
            Named("X_all_with_y_scaled") = X_all_with_y_scaled,
            Named("rank_all_with_y_scaled") = rank_all_with_y_scaled,
            Named("cols_all_with_y_scaled") = all_col_indices
        );
    }

    // Process the 4 needed variants
    // 1. past_info: i_past (unscaled) - for X_prev, rank_prev, xt_prev
    SubsetResult past_info = process_subset(X, t, i_past, false, rank_tol);

    // 2. all_info: i_all (unscaled) - only need X_all
    SubsetResult all_info = process_subset(X, t, i_all, false, rank_tol);

    // 3. all_info_scaled: i_all (scaled) - for X_all_scaled, xt_all_scaled
    SubsetResult all_info_scaled = process_subset(X, t, i_all, true, rank_tol);

    // 4. all_info_with_y_scaled: i_all_y_present (scaled) - for X_all_with_y_scaled, rank_all_with_y_scaled
    SubsetResult all_info_with_y_scaled = process_subset(X, t, i_all_y_present, true, rank_tol);

    return List::create(
        Named("X_prev") = past_info.Xint,
        Named("rank_prev") = past_info.rank,
        Named("xt_prev") = past_info.xt,
        Named("cols_prev") = past_info.col_indices,
        Named("X_all") = all_info.Xint,
        Named("cols_all") = all_info.col_indices,
        Named("X_all_scaled") = all_info_scaled.Xint,
        Named("xt_all_scaled") = all_info_scaled.xt,
        Named("cols_all_scaled") = all_info_scaled.col_indices,
        Named("X_all_with_y_scaled") = all_info_with_y_scaled.Xint,
        Named("rank_all_with_y_scaled") = all_info_with_y_scaled.rank,
        Named("cols_all_with_y_scaled") = all_info_with_y_scaled.col_indices
    );
}
