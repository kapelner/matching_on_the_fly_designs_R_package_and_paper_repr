#ifndef EDI_RESULT_MAP_RCPP_H
#define EDI_RESULT_MAP_RCPP_H

// The only place that includes both result_map.h and Rcpp.h -- never
// included by a core (Rcpp-free) .cpp/.h file, only by the thin
// // [[Rcpp::export]] wrapper functions.

#include "result_map.h"
#include <Rcpp.h>

namespace edi {

inline Rcpp::List to_rcpp_list(const ResultMap& m) {
	Rcpp::List out;
	Rcpp::CharacterVector names;
	for (const auto& [name, value] : m.entries()) {
		out.push_back(std::visit([](auto&& v) -> SEXP {
			using T = std::decay_t<decltype(v)>;
			if constexpr (std::is_same_v<T, std::monostate>) {
				return R_NilValue;
			} else {
				return Rcpp::wrap(v);
			}
		}, value));
		names.push_back(name);
	}
	out.names() = names;
	return out;
}

} // namespace edi

#endif
