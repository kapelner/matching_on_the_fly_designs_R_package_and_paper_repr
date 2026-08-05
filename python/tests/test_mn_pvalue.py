"""Parity test for edi_kernels.mn_pvalue against an R-generated fixture.

This does NOT re-verify the statistical correctness of the Miettinen-Nurminen
score p-value itself -- that's the same mn_pvalue_cpp object code the R
package's own test suite already exercises (mn_pvalue_cpp is entirely
scalar-double math with zero Rcpp-type dependency, so the EDI_CORE_ONLY and
R builds compile the identical function body). What this test covers: the
pybind11 binding layer itself and the EDI_CORE_ONLY compiled path.

The expected value below was computed once via:
    EDI:::mn_pvalue_cpp(x_t, n_t, x_c, n_c, delta, p_t, p_c)
in R on x_t=30, n_t=100, x_c=18, n_c=95, delta=0.
"""
import math

import pytest

from edi_kernels import mn_pvalue

ATOL = 1e-9
RTOL = 1e-9

X_T, N_T, X_C, N_C = 30.0, 100.0, 18.0, 95.0
P_T, P_C = X_T / N_T, X_C / N_C

# R reference, EDI 1.0.0, computed 2026-08-05 via:
#   EDI:::mn_pvalue_cpp(30, 100, 18, 95, 0, 0.3, 18/95)
R_PVALUE = 0.0740542399588677


def test_matches_r_fixture():
    got = mn_pvalue(X_T, N_T, X_C, N_C, 0.0, P_T, P_C)
    assert got == pytest.approx(R_PVALUE, abs=ATOL, rel=RTOL)


def test_pvalue_in_unit_interval():
    got = mn_pvalue(X_T, N_T, X_C, N_C, 0.0, P_T, P_C)
    assert 0.0 <= got <= 1.0


def test_delta_equal_to_observed_diff_gives_pvalue_near_one():
    delta = P_T - P_C
    got = mn_pvalue(X_T, N_T, X_C, N_C, delta, P_T, P_C)
    assert got == pytest.approx(1.0, abs=1e-6)


def test_zero_n_returns_nan():
    got = mn_pvalue(0.0, 0.0, X_C, N_C, 0.0, 0.0, P_C)
    assert math.isnan(got)


def test_out_of_range_delta_returns_nan():
    assert math.isnan(mn_pvalue(X_T, N_T, X_C, N_C, -1.5, P_T, P_C))
    assert math.isnan(mn_pvalue(X_T, N_T, X_C, N_C, 1.5, P_T, P_C))
