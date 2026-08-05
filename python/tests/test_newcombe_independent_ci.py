"""Parity test for edi_kernels.newcombe_independent_ci against an
R-generated fixture.

This does NOT re-verify the statistical correctness of the Newcombe hybrid
score interval itself -- that's the same newcombe_independent_ci_internal
object code the R package's own test suite already exercises via
EDI:::newcombe_independent_ci_cpp (both call the identical
newcombe_independent_ci_internal in EDI/src/newcombe_speedups.cpp;
newcombe_independent_ci_cpp is now a thin Rcpp::NumericVector wrapper around
it). What this test covers: the pybind11 binding layer itself and the
EDI_CORE_ONLY compiled path.

The expected values below were computed once via:
    EDI:::newcombe_independent_ci_cpp(x1, n1, x2, n2, alpha)
in R on x1=30, n1=100, x2=18, n2=95, alpha=0.05.
"""
import math

import pytest

from edi_kernels import newcombe_independent_ci

ATOL = 1e-9
RTOL = 1e-9

X1, N1, X2, N2 = 30.0, 100.0, 18.0, 95.0

# R reference, EDI 1.0.0, computed 2026-08-05 via:
#   EDI:::newcombe_independent_ci_cpp(30, 100, 18, 95, 0.05)
R_LOWER = -0.0107855464549578
R_UPPER = 0.226971499864617


def test_matches_r_fixture():
    lower, upper = newcombe_independent_ci(X1, N1, X2, N2)
    assert lower == pytest.approx(R_LOWER, abs=ATOL, rel=RTOL)
    assert upper == pytest.approx(R_UPPER, abs=ATOL, rel=RTOL)


def test_default_alpha_matches_r_default():
    # R's compute_asymp_confidence_interval() defaults to alpha=0.05.
    default = newcombe_independent_ci(X1, N1, X2, N2)
    explicit = newcombe_independent_ci(X1, N1, X2, N2, alpha=0.05)
    assert default == explicit


def test_lower_less_than_upper():
    lower, upper = newcombe_independent_ci(X1, N1, X2, N2)
    assert lower < upper


def test_ci_contains_point_estimate():
    lower, upper = newcombe_independent_ci(X1, N1, X2, N2)
    est = X1 / N1 - X2 / N2
    assert lower <= est <= upper


def test_narrower_alpha_gives_wider_interval():
    lower95, upper95 = newcombe_independent_ci(X1, N1, X2, N2, alpha=0.05)
    lower99, upper99 = newcombe_independent_ci(X1, N1, X2, N2, alpha=0.01)
    assert (upper99 - lower99) > (upper95 - lower95)


def test_zero_n_returns_nan():
    lower, upper = newcombe_independent_ci(0.0, 0.0, X2, N2)
    assert math.isnan(lower)
    assert math.isnan(upper)
