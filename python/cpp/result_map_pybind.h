#ifndef EDI_RESULT_MAP_PYBIND_H
#define EDI_RESULT_MAP_PYBIND_H

// The pybind11-side twin of EDI/src/result_map_rcpp.h -- converts a
// portable edi::ResultMap (built by an EDI/src/*_internal function, no
// Rcpp/R dependency) into a Python dict. Mirrors that file's std::visit
// shape exactly so the two boundary layers stay in lockstep; see
// package_metadata/python_bindings_package_spec.md's "Result Conversion".

#include "result_map.h"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace edi {

inline py::dict to_py_dict(const ResultMap& m) {
	py::dict out;
	for (const auto& [name, value] : m.entries()) {
		out[py::str(name)] = std::visit([](auto&& v) -> py::object {
			using T = std::decay_t<decltype(v)>;
			if constexpr (std::is_same_v<T, std::monostate>) {
				return py::none();
			} else {
				return py::cast(v);
			}
		}, value);
	}
	return out;
}

} // namespace edi

#endif
