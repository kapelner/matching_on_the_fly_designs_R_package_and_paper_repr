"""Parity test for edi_kernels.mn_ci against an R-generated fixture.

This does NOT re-verify the statistical correctness of the Miettinen-Nurminen
score CI itself -- that's the same mn_ci_internal object code the R
package's own test suite already exercises via EDI:::mn_ci_cpp (both call
the identical mn_ci_internal in EDI/src/miettinen_nurminen_speedups.cpp;
mn_ci_cpp is now a thin Rcpp::NumericVector wrapper around it). What this
test covers: the pybind11 binding layer itself (argument marshaling,
defaults, tuple conversion) and the EDI_CORE_ONLY compiled path.

The expected values below were computed once via:
    EDI:::mn_ci_cpp(x_t, n_t, x_c, n_c, p_t, p_c, alpha, pval_epsilon)
in R on x_t=30, n_t=100, x_c=18, n_c=95 (see
package_metadata/python_bindings_package_spec.md "Testing").
"""
import pytest

from edi_kernels import mn_ci

ATOL = 1e-9
RTOL = 1e-9

X_T, N_T, X_C, N_C = 30.0, 100.0, 18.0, 95.0
P_T, P_C = X_T / N_T, X_C / N_C

# R reference, EDI 1.0.0, computed 2026-08-05 via:
#   EDI:::mn_ci_cpp(30, 100, 18, 95, 0.3, 18/95, 0.05, 1e-7)
R_LOWER = -0.0109324923824448
R_UPPER = 0.22947596992664


def test_matches_r_fixture():
    lower, upper = mn_ci(X_T, N_T, X_C, N_C, P_T, P_C)
    assert lower == pytest.approx(R_LOWER, abs=ATOL, rel=RTOL)
    assert upper == pytest.approx(R_UPPER, abs=ATOL, rel=RTOL)


def test_default_alpha_and_pval_epsilon_match_r_defaults():
    # R's compute_asymp_confidence_interval() defaults to alpha=0.05,
    # pval_epsilon=1e-7 -- the binding's defaults must match so callers who
    # omit them get the same fixture above.
    lower, upper = mn_ci(X_T, N_T, X_C, N_C, P_T, P_C)
    lower_explicit, upper_explicit = mn_ci(X_T, N_T, X_C, N_C, P_T, P_C, alpha=0.05, pval_epsilon=1e-7)
    assert (lower, upper) == (lower_explicit, upper_explicit)


def test_lower_less_than_upper():
    lower, upper = mn_ci(X_T, N_T, X_C, N_C, P_T, P_C)
    assert lower < upper


def test_ci_contains_point_estimate():
    lower, upper = mn_ci(X_T, N_T, X_C, N_C, P_T, P_C)
    assert lower <= (P_T - P_C) <= upper


def test_zero_n_t_matches_r_fixture():
    # n_t=0 does NOT propagate to NaN here: mn_ci_internal's bisection
    # (find_bound) treats a non-finite p-value from mn_pvalue_cpp as p=0.0
    # and keeps searching rather than short-circuiting, so both bounds
    # collapse near p_c instead. Confirmed this is EDI's actual existing
    # behavior (not a binding bug) via:
    #   EDI:::mn_ci_cpp(0, 0, 18, 95, 0, 18/95, 0.05, 1e-7)
    # in R, which returns the same values.
    lower, upper = mn_ci(0.0, 0.0, X_C, N_C, 0.0, P_C)
    assert lower == pytest.approx(-0.189473684210527, abs=ATOL, rel=RTOL)
    assert upper == pytest.approx(-0.189473684210526, abs=ATOL, rel=RTOL)
