#include <RcppEigen.h>
#include <vector>
#include <array>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix bootstrap_m_indices_cpp(
	const IntegerVector& m_vec,
	const IntegerVector& i_reservoir,
	int n_reservoir,
	int m,
	int B
) {
	int row_length = n_reservoir + 2 * m;
	IntegerMatrix result(B, row_length);

	std::vector< std::array<int, 2> > match_pairs(m);
	std::vector<int> count(m, 0);
	for (int idx = 0; idx < m_vec.size(); ++idx) {
	int match_id = m_vec[idx];
	if (match_id > 0 && match_id <= m) {
		int pos = count[match_id - 1]++;
		if (pos < 2) {
		match_pairs[match_id - 1][pos] = idx + 1;
		}
	}
	}

	for (int row = 0; row < B; ++row) {
	if (n_reservoir > 0) {
		for (int j = 0; j < n_reservoir; ++j) {
		int idx = static_cast<int>(R::runif(0.0, 1.0) * n_reservoir);
		if (idx == n_reservoir) idx = n_reservoir - 1;
		result(row, j) = i_reservoir[idx];
		}
	}

	for (int k = 0; k < m; ++k) {
		int match_id = static_cast<int>(R::runif(0.0, 1.0) * m);
		if (match_id == m) match_id = m - 1;
		auto pair = match_pairs[match_id];
		result(row, n_reservoir + 2 * k) = pair[0];
		result(row, n_reservoir + 2 * k + 1) = pair[1];
	}
	}
	return result;
}

List compute_zhang_match_data_cpp(const IntegerVector& w,
								  const IntegerVector& m_vec,
								  const NumericVector& y,
								  const NumericMatrix& X);

// [[Rcpp::export]]
List match_stats_from_indices_cpp(
	const NumericVector& y,
	const NumericVector& w,
	const NumericMatrix& X,
	const IntegerVector& original_m_vec, // Changed name
	const IntegerVector& i_b,
	int m
) {
	int n_rows = i_b.size();
	int p = X.ncol();
	NumericVector y_sample(n_rows);
	NumericVector w_sample(n_rows);
	NumericMatrix X_sample(n_rows, p);
	IntegerVector m_vec_sample(n_rows); // This will hold the resampled m_vec

	for (int i = 0; i < n_rows; ++i) {
	int idx = i_b[i] - 1;
	y_sample[i] = y[idx];
	w_sample[i] = w[idx];
	m_vec_sample[i] = original_m_vec[idx]; // Corrected resampling
	for (int j = 0; j < p; ++j) {
		X_sample(i, j) = X(idx, j);
	}
	}

	List match_data = compute_zhang_match_data_cpp(
		as<IntegerVector>(w_sample),
		m_vec_sample,
		y_sample,
		X_sample
	);

	return match_data;
}
