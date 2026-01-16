#include <Rcpp.h>
#include <unordered_map>

using namespace Rcpp;

namespace {

struct CountInfo {
	int count;
	int first_index;
};

inline void update_info(CountInfo& info, int idx) {
	if (info.count == 0) {
		info.first_index = idx;
	}
	++info.count;
}

struct SexpHash {
	std::size_t operator()(SEXP x) const {
		return std::hash<void*>()(static_cast<void*>(x));
	}
};

struct SexpEq {
	bool operator()(SEXP a, SEXP b) const {
		return a == b;
	}
};

Rcpp::RObject sample_mode_integer(const IntegerVector& data, bool is_factor) {
	int n = data.size();
	if (n == 0) {
		return data;
	}

	std::unordered_map<int, CountInfo> counts;
	counts.reserve(static_cast<size_t>(n));
	CountInfo na_info{0, -1};

	for (int i = 0; i < n; ++i) {
		int val = data[i];
		if (val == NA_INTEGER) {
			update_info(na_info, i);
			continue;
		}
		auto it = counts.find(val);
		if (it == counts.end()) {
			counts.emplace(val, CountInfo{1, i});
		} else {
			++it->second.count;
		}
	}

	int best_count = -1;
	int best_first = n + 1;
	int best_val = NA_INTEGER;
	bool best_is_na = false;

	if (na_info.count > 0) {
		best_count = na_info.count;
		best_first = na_info.first_index;
		best_is_na = true;
	}

	for (const auto& kv : counts) {
		const CountInfo& info = kv.second;
		if (info.count > best_count || (info.count == best_count && info.first_index < best_first)) {
			best_count = info.count;
			best_first = info.first_index;
			best_val = kv.first;
			best_is_na = false;
		}
	}

	IntegerVector res(1);
	res[0] = best_is_na ? NA_INTEGER : best_val;
	if (is_factor) {
		res.attr("class") = data.attr("class");
		res.attr("levels") = data.attr("levels");
	}
	return res;
}

Rcpp::RObject sample_mode_logical(const LogicalVector& data) {
	int n = data.size();
	if (n == 0) {
		return data;
	}

	std::unordered_map<int, CountInfo> counts;
	counts.reserve(static_cast<size_t>(n));
	CountInfo na_info{0, -1};

	for (int i = 0; i < n; ++i) {
		int val = data[i];
		if (val == NA_LOGICAL) {
			update_info(na_info, i);
			continue;
		}
		auto it = counts.find(val);
		if (it == counts.end()) {
			counts.emplace(val, CountInfo{1, i});
		} else {
			++it->second.count;
		}
	}

	int best_count = -1;
	int best_first = n + 1;
	int best_val = NA_LOGICAL;
	bool best_is_na = false;

	if (na_info.count > 0) {
		best_count = na_info.count;
		best_first = na_info.first_index;
		best_is_na = true;
	}

	for (const auto& kv : counts) {
		const CountInfo& info = kv.second;
		if (info.count > best_count || (info.count == best_count && info.first_index < best_first)) {
			best_count = info.count;
			best_first = info.first_index;
			best_val = kv.first;
			best_is_na = false;
		}
	}

	LogicalVector res(1);
	res[0] = best_is_na ? NA_LOGICAL : best_val;
	return res;
}

Rcpp::RObject sample_mode_numeric(const NumericVector& data) {
	int n = data.size();
	if (n == 0) {
		return data;
	}

	std::unordered_map<double, CountInfo> counts;
	counts.reserve(static_cast<size_t>(n));
	CountInfo na_info{0, -1};
	CountInfo nan_info{0, -1};

	for (int i = 0; i < n; ++i) {
		double val = data[i];
		if (R_IsNA(val)) {
			update_info(na_info, i);
			continue;
		}
		if (R_IsNaN(val)) {
			update_info(nan_info, i);
			continue;
		}
		auto it = counts.find(val);
		if (it == counts.end()) {
			counts.emplace(val, CountInfo{1, i});
		} else {
			++it->second.count;
		}
	}

	int best_count = -1;
	int best_first = n + 1;
	int best_kind = 0; // 0 normal, 1 NA, 2 NaN
	double best_val = NA_REAL;

	if (na_info.count > 0) {
		best_count = na_info.count;
		best_first = na_info.first_index;
		best_kind = 1;
	}
	if (nan_info.count > 0 &&
		(nan_info.count > best_count || (nan_info.count == best_count && nan_info.first_index < best_first))) {
		best_count = nan_info.count;
		best_first = nan_info.first_index;
		best_kind = 2;
	}

	for (const auto& kv : counts) {
		const CountInfo& info = kv.second;
		if (info.count > best_count || (info.count == best_count && info.first_index < best_first)) {
			best_count = info.count;
			best_first = info.first_index;
			best_kind = 0;
			best_val = kv.first;
		}
	}

	NumericVector res(1);
	if (best_kind == 1) {
		res[0] = NA_REAL;
	} else if (best_kind == 2) {
		res[0] = R_NaN;
	} else {
		res[0] = best_val;
	}
	return res;
}

Rcpp::RObject sample_mode_character(const CharacterVector& data) {
	int n = data.size();
	if (n == 0) {
		return data;
	}

	std::unordered_map<SEXP, CountInfo, SexpHash, SexpEq> counts;
	counts.reserve(static_cast<size_t>(n));
	CountInfo na_info{0, -1};

	for (int i = 0; i < n; ++i) {
		SEXP val = STRING_ELT(data, i);
		if (val == NA_STRING) {
			update_info(na_info, i);
			continue;
		}
		auto it = counts.find(val);
		if (it == counts.end()) {
			counts.emplace(val, CountInfo{1, i});
		} else {
			++it->second.count;
		}
	}

	int best_count = -1;
	int best_first = n + 1;
	SEXP best_val = NA_STRING;
	bool best_is_na = false;

	if (na_info.count > 0) {
		best_count = na_info.count;
		best_first = na_info.first_index;
		best_is_na = true;
	}

	for (const auto& kv : counts) {
		const CountInfo& info = kv.second;
		if (info.count > best_count || (info.count == best_count && info.first_index < best_first)) {
			best_count = info.count;
			best_first = info.first_index;
			best_val = kv.first;
			best_is_na = false;
		}
	}

	CharacterVector res(1);
	res[0] = best_is_na ? NA_STRING : best_val;
	return res;
}

} // namespace

// [[Rcpp::export]]
SEXP sample_mode_cpp(SEXP data) {
	switch (TYPEOF(data)) {
	case INTSXP:
		return sample_mode_integer(IntegerVector(data), Rf_isFactor(data));
	case LGLSXP:
		return sample_mode_logical(LogicalVector(data));
	case REALSXP:
		return sample_mode_numeric(NumericVector(data));
	case STRSXP:
		return sample_mode_character(CharacterVector(data));
	default:
		stop("sample_mode_cpp: unsupported type");
	}

	return R_NilValue;
}
