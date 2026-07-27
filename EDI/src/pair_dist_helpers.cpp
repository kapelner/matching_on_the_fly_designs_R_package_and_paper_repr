#include <RcppEigen.h>
#include <limits>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
NumericMatrix compute_pair_averages_cpp(const NumericMatrix& X,
										const IntegerVector& m_vec,
										const int m) {
	int n = X.nrow();
	int p = X.ncol();

	NumericMatrix pair_avg(m, p);
	IntegerVector counts(m);
	const double* X_ptr = X.begin();
	double* pair_ptr = pair_avg.begin();
	const int* m_ptr = m_vec.begin();

	for (int i = 0; i < n; ++i) {
		const int id = m_ptr[i];
		if (id > 0 && id <= m) {
			counts[id - 1] += 1;
		}
	}

	for (int j = 0; j < p; ++j) {
		const double* X_col = X_ptr + static_cast<size_t>(j) * n;
		double* avg_col = pair_ptr + static_cast<size_t>(j) * m;
		for (int i = 0; i < n; ++i) {
			const int id = m_ptr[i];
			if (id > 0 && id <= m) {
				avg_col[id - 1] += X_col[i];
			}
		}
	}

	for (int i = 0; i < m; ++i) {
	if (counts[i] > 0) {
		double denom = static_cast<double>(counts[i]);
		for (int j = 0; j < p; ++j) {
		pair_avg(i, j) /= denom;
		}
	} else {
		for (int j = 0; j < p; ++j) {
		pair_avg(i, j) = NA_REAL;
		}
	}
	}

	return pair_avg;
}

// [[Rcpp::export]]
NumericMatrix compute_pair_distance_matrix_cpp(const NumericMatrix& pair_avg,
                                                 const NumericVector& weights) {
	int m = pair_avg.nrow();
	int p = pair_avg.ncol();
	bool use_weights = weights.size() == p;

	NumericMatrix dist_mat(m, m);
	double inf = std::numeric_limits<double>::infinity();
	const double* weight_ptr = use_weights ? weights.begin() : nullptr;
	double* dist_ptr = dist_mat.begin();
	Eigen::Map<const Eigen::MatrixXd> pair_col(pair_avg.begin(), m, p);
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pair_row = pair_col;

	for (int i = 0; i < m; ++i) {
		dist_ptr[i + static_cast<size_t>(m) * i] = inf;
	}

	for (int i = 0; i < m; ++i) {
		const double* row_i = pair_row.data() + static_cast<size_t>(i) * p;
		for (int j = i + 1; j < m; ++j) {
			const double* row_j = pair_row.data() + static_cast<size_t>(j) * p;
			double sum = 0.0;
#pragma omp simd reduction(+:sum)
			for (int k = 0; k < p; ++k) {
				const double diff = row_i[k] - row_j[k];
				const double w = use_weights ? weight_ptr[k] : 1.0;
				sum += w * diff * diff;
			}
			dist_ptr[i + static_cast<size_t>(m) * j] = sum;
			dist_ptr[j + static_cast<size_t>(m) * i] = sum;
		}
	}

	return dist_mat;
}

// [[Rcpp::export]]
double compute_lambda_squ_cpp(const NumericVector& d_i,
								const IntegerMatrix& halves) {
	int n_halves = halves.nrow();
	if (n_halves == 0) {
	return 0.0;
	}

	double sum = 0.0;
	for (int i = 0; i < n_halves; ++i) {
	int id1 = halves(i, 0);
	int id2 = halves(i, 1);
	if (id1 > 0 && id2 > 0 && id1 <= d_i.size() && id2 <= d_i.size()) {
		sum += d_i[id1 - 1] * d_i[id2 - 1];
	}
	}

	return sum / static_cast<double>(n_halves);
}
