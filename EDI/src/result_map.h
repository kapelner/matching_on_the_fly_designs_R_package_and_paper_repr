#ifndef EDI_RESULT_MAP_H
#define EDI_RESULT_MAP_H

// Rcpp-free core boundary: no <Rcpp.h> here, only Eigen types. Model-fitting
// core functions build a ResultMap; the thin exported wrapper converts it to
// an Rcpp::List via result_map_rcpp.h's to_rcpp_list() (Python's pybind11
// wrapper does the analogous thing via result_map_pybind.h's to_py_dict()).
// See package_metadata/sexp_removal_rcppeigen_conversion_spec.md and
// package_metadata/python_bindings_package_spec.md.
//
// EDI_CORE_ONLY: see fast_gamma_functions.h for the convention. This header
// only ever needed Eigen::VectorXd/MatrixXd (never SEXP), so the two
// branches differ only in which header supplies those types.
#ifdef EDI_CORE_ONLY
#include <Eigen/Dense>
#else
#include <RcppEigen.h>
#endif
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace edi {

// std::monostate represents R_NilValue / Python None (typed "not available").
using ResultValue = std::variant<
	std::monostate,
	bool,
	int,
	double,
	std::string,
	Eigen::VectorXd,
	Eigen::MatrixXd
>;

// Eigen expressions (e.g. par.head(p), a block, X.transpose() * y) are not
// ResultValue themselves -- they're lazy expression templates. Materialize
// them into a concrete VectorXd/MatrixXd at the point of insertion, which is
// exactly the right place to do it (this is the R/Python boundary).
template <typename Derived>
inline ResultValue to_result_value(const Eigen::MatrixBase<Derived>& v) {
	if constexpr (Derived::ColsAtCompileTime == 1 || Derived::RowsAtCompileTime == 1) {
		return ResultValue{Eigen::VectorXd(v)};
	} else {
		return ResultValue{Eigen::MatrixXd(v)};
	}
}

inline ResultValue to_result_value(ResultValue v) {
	return v;
}

class ResultMap {
public:
	// Reserves up front so a chain of .set() calls doesn't pay for vector
	// reallocations as it grows (observed field counts across the 33
	// model-fitting files range 6-39; 100 is deliberately generous headroom
	// for future growth, since a single reserve()'d allocation this size is
	// negligible either way -- 20 vs 100 entries made no measurable
	// difference in practice). This was NOT the dominant cost, though: most
	// of the ResultMap-vs-List::create construction gap (~8.7us for a
	// 13-field result) comes from a structural double-copy of Eigen-typed
	// fields (to_result_value() copies into the variant, then
	// to_rcpp_list()'s wrap() copies again into the SEXP) inherent to
	// ResultMap's deliberate Rcpp-free two-stage design -- reserve() only
	// recovers the reallocation share of that gap (~16% of it, measured).
	// See package_metadata/perf_experiments_final.md TODO-130 follow-up.
	ResultMap() { entries_.reserve(100); }

	template <typename T>
	ResultMap& set(std::string name, const T& value) {
		entries_.emplace_back(std::move(name), to_result_value(value));
		return *this;
	}

	const std::vector<std::pair<std::string, ResultValue>>& entries() const {
		return entries_;
	}

	// Lookup for the rare case where one core function's ResultMap is a
	// producer another core function needs to read back from internally
	// (e.g. a *_with_var variant reusing a plain fit's result) -- not needed
	// by the common case (build once, convert once at the R/Python
	// boundary), so this stays a linear scan over the (small) entry list
	// rather than adding an index the common case would pay to maintain.
	template <typename T>
	const T* get_if(const std::string& name) const {
		for (const auto& [n, v] : entries_) {
			if (n == name) return std::get_if<T>(&v);
		}
		return nullptr;
	}

	bool has(const std::string& name) const {
		for (const auto& [n, v] : entries_) {
			if (n == name) return true;
		}
		return false;
	}

private:
	std::vector<std::pair<std::string, ResultValue>> entries_;
};

} // namespace edi

#endif
