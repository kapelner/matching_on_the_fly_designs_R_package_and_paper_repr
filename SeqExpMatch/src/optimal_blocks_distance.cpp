#include <Rcpp.h>
using namespace Rcpp;

namespace {

enum DistanceCode {
  DIST_EUCLIDEAN_SQ = 1,
  DIST_SUM_ABS_DIFF = 2,
  DIST_MAHAL = 3
};

NumericMatrix invert_matrix_cpp(const NumericMatrix& A) {
  const int p = A.nrow();
  if (A.ncol() != p) {
    stop("Matrix inversion requires a square matrix.");
  }

  NumericMatrix aug(p, 2 * p);
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p; ++j) {
      aug(i, j) = A(i, j);
      aug(i, p + j) = (i == j) ? 1.0 : 0.0;
    }
  }

  for (int col = 0; col < p; ++col) {
    int pivot_row = col;
    double pivot_abs = std::abs(aug(col, col));
    for (int row = col + 1; row < p; ++row) {
      const double cand = std::abs(aug(row, col));
      if (cand > pivot_abs) {
        pivot_abs = cand;
        pivot_row = row;
      }
    }
    if (pivot_abs < 1e-12) {
      stop("Mahalanobis covariance matrix is numerically singular even after ridge regularization.");
    }
    if (pivot_row != col) {
      for (int j = 0; j < 2 * p; ++j) {
        const double tmp = aug(col, j);
        aug(col, j) = aug(pivot_row, j);
        aug(pivot_row, j) = tmp;
      }
    }

    const double pivot = aug(col, col);
    for (int j = 0; j < 2 * p; ++j) {
      aug(col, j) /= pivot;
    }

    for (int row = 0; row < p; ++row) {
      if (row == col) {
        continue;
      }
      const double factor = aug(row, col);
      if (factor == 0.0) {
        continue;
      }
      for (int j = 0; j < 2 * p; ++j) {
        aug(row, j) -= factor * aug(col, j);
      }
    }
  }

  NumericMatrix inv(p, p);
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p; ++j) {
      inv(i, j) = aug(i, p + j);
    }
  }
  return inv;
}

NumericMatrix compute_mahal_inverse_covariance(const NumericMatrix& X) {
  const int n = X.nrow();
  const int p = X.ncol();
  NumericVector means(p);
  NumericMatrix cov(p, p);

  for (int j = 0; j < p; ++j) {
    double total = 0.0;
    for (int i = 0; i < n; ++i) {
      total += X(i, j);
    }
    means[j] = total / static_cast<double>(n);
  }

  if (n > 1) {
    for (int a = 0; a < p; ++a) {
      for (int b = a; b < p; ++b) {
        double acc = 0.0;
        for (int i = 0; i < n; ++i) {
          const double xa = X(i, a) - means[a];
          const double xb = X(i, b) - means[b];
          acc += xa * xb;
        }
        const double cov_ab = acc / static_cast<double>(n - 1);
        cov(a, b) = cov_ab;
        cov(b, a) = cov_ab;
      }
    }
  } else {
    for (int j = 0; j < p; ++j) {
      cov(j, j) = 1.0;
    }
  }

  for (int a = 0; a < p; ++a) {
    for (int b = 0; b < p; ++b) {
      if (!R_finite(cov(a, b))) {
        cov(a, b) = 0.0;
      }
    }
    cov(a, a) += 1e-8;
  }

  return invert_matrix_cpp(cov);
}

NumericMatrix distance_matrix_custom_cpp_impl(const NumericMatrix& X, const Function& dist_fn) {
  const int n = X.nrow();
  const int p = X.ncol();
  NumericMatrix D(n, n);
  for (int i = 0; i < n; ++i) {
    D(i, i) = 0.0;
    for (int j = i + 1; j < n; ++j) {
      NumericVector xi(p);
      NumericVector xj(p);
      for (int k = 0; k < p; ++k) {
        xi[k] = X(i, k);
        xj[k] = X(j, k);
      }
      SEXP d_obj = dist_fn(xi, xj);
      if (!Rf_isNumeric(d_obj) || Rf_length(d_obj) != 1) {
        stop("Custom dist must return one finite nonnegative number per subject pair.");
      }
      const double d_ij = as<double>(d_obj);
      if (!R_finite(d_ij) || d_ij < 0.0) {
        stop("Custom dist must return one finite nonnegative number per subject pair.");
      }
      D(i, j) = d_ij;
      D(j, i) = d_ij;
    }
  }
  return D;
}

}  // namespace

// [[Rcpp::export]]
NumericMatrix distance_matrix_euclidean_sq_cpp(const NumericMatrix& X) {
  const int n = X.nrow();
  const int p = X.ncol();
  NumericMatrix D(n, n);
  for (int i = 0; i < n; ++i) {
    D(i, i) = 0.0;
    for (int j = i + 1; j < n; ++j) {
      double acc = 0.0;
      for (int k = 0; k < p; ++k) {
        const double diff = X(i, k) - X(j, k);
        acc += diff * diff;
      }
      D(i, j) = acc;
      D(j, i) = acc;
    }
  }
  return D;
}

// [[Rcpp::export]]
NumericMatrix distance_matrix_custom_cpp(const NumericMatrix& X, const Function& dist_fn) {
  return distance_matrix_custom_cpp_impl(X, dist_fn);
}

// [[Rcpp::export]]
NumericMatrix distance_matrix_sum_abs_diff_cpp(const NumericMatrix& X) {
  const int n = X.nrow();
  const int p = X.ncol();
  NumericMatrix D(n, n);
  for (int i = 0; i < n; ++i) {
    D(i, i) = 0.0;
    for (int j = i + 1; j < n; ++j) {
      double acc = 0.0;
      for (int k = 0; k < p; ++k) {
        acc += std::abs(X(i, k) - X(j, k));
      }
      D(i, j) = acc;
      D(j, i) = acc;
    }
  }
  return D;
}

// [[Rcpp::export]]
NumericMatrix distance_matrix_mahal_cpp(const NumericMatrix& X) {
  const int n = X.nrow();
  const int p = X.ncol();
  NumericMatrix D(n, n);
  NumericMatrix S_inv = compute_mahal_inverse_covariance(X);
  NumericVector diff(p);
  NumericVector tmp(p);
  for (int i = 0; i < n; ++i) {
    D(i, i) = 0.0;
    for (int j = i + 1; j < n; ++j) {
      for (int k = 0; k < p; ++k) {
        diff[k] = X(i, k) - X(j, k);
      }
      for (int a = 0; a < p; ++a) {
        double acc = 0.0;
        for (int b = 0; b < p; ++b) {
          acc += S_inv(a, b) * diff[b];
        }
        tmp[a] = acc;
      }
      double dist = 0.0;
      for (int k = 0; k < p; ++k) {
        dist += diff[k] * tmp[k];
      }
      D(i, j) = dist;
      D(j, i) = dist;
    }
  }
  return D;
}

// [[Rcpp::export]]
NumericMatrix optimal_blocks_distance_matrix_cpp(const NumericMatrix& X,
                                                 const int dist_code,
                                                 Nullable<Function> dist_fn = R_NilValue) {
  switch (dist_code) {
  case DIST_EUCLIDEAN_SQ:
    return distance_matrix_euclidean_sq_cpp(X);
  case DIST_SUM_ABS_DIFF:
    return distance_matrix_sum_abs_diff_cpp(X);
  case DIST_MAHAL:
    return distance_matrix_mahal_cpp(X);
  default:
    if (dist_fn.isNull()) {
      stop("Custom distance dispatch requires a distance function.");
    }
    return distance_matrix_custom_cpp_impl(X, as<Function>(dist_fn));
  }
}
