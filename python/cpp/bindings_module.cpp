// Single pybind11 module entry point. Each family gets its own
// bindings_<family>.cpp with a bind_<family>(py::module_&) function (mirrors
// EDI/src's own per-model-family .cpp layout); this file just wires them
// together into one compiled extension.

#include <pybind11/pybind11.h>

namespace py = pybind11;

void bind_fast_math(py::module_& m);
void bind_glmm(py::module_& m);
void bind_continuous(py::module_& m);
void bind_binary(py::module_& m);
void bind_count(py::module_& m);
void bind_proportion(py::module_& m);
void bind_ordinal(py::module_& m);
void bind_incidence(py::module_& m);
void bind_survival(py::module_& m);

PYBIND11_MODULE(_core, m) {
    m.doc() = "EDI kernels (pybind11) -- compiled directly from EDI/src, no copies";

    bind_fast_math(m);
    bind_glmm(m);
    bind_continuous(m);
    bind_binary(m);
    bind_count(m);
    bind_proportion(m);
    bind_ordinal(m);
    bind_incidence(m);
    bind_survival(m);
}
