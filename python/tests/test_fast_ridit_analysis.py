"""Parity test for edi_kernels.fast_ridit_analysis against an R-generated
fixture.

This does NOT re-verify the statistical correctness of Ridit analysis
itself -- that's the same fast_ridit_scores_from_ref_indices object code
the R package's own test suite already exercises via
EDI:::fast_ridit_analysis_cpp (both call the identical algorithm; only
RiditScoreData's internal storage changed from Rcpp::NumericVector to
std::vector<double> so a separate Python binding translation unit could
call it, and fast_ridit_analysis_result was added as a portable sibling of
fast_ridit_analysis_cpp avoiding its SEXP/Rcpp::List boundary). What this
test covers: the pybind11 binding layer itself and the EDI_CORE_ONLY
compiled path.

The expected values below were computed once via:
    EDI:::fast_ridit_analysis_cpp(w, y, "control")
in R on a synthetic n=45 dataset generated with set.seed(654).
"""
import math

import numpy as np
import pytest

from edi_kernels import fast_ridit_analysis

ATOL = 1e-9
RTOL = 1e-9

Y = np.array([
    2, 4, 3, 2, 1, 2, 3, 2, 3, 4, 4, 3, 1, 2, 1, 4, 4, 1, 3, 1, 3, 1, 2, 4,
    2, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 2, 2, 2, 3, 4, 3, 1, 4,
], dtype=np.int32)
W = np.array([1, 0] * 22 + [1], dtype=np.int32)

# R reference, EDI 1.0.0, computed 2026-08-05 via:
#   EDI:::fast_ridit_analysis_cpp(w, y, "control")
R_MEAN_RIDIT_T = 0.558300395256917
R_MEAN_RIDIT_C = 0.5
R_ESTIMATE = 0.0583003952569171
R_SE = 0.0487241695783548
R_LEVELS = [1, 2, 3, 4]


def test_matches_r_fixture():
    res = fast_ridit_analysis(W, Y, "control")
    assert res["mean_ridit_t"] == pytest.approx(R_MEAN_RIDIT_T, abs=ATOL, rel=RTOL)
    assert res["mean_ridit_c"] == pytest.approx(R_MEAN_RIDIT_C, abs=ATOL, rel=RTOL)
    assert res["estimate"] == pytest.approx(R_ESTIMATE, abs=ATOL, rel=RTOL)
    assert res["se"] == pytest.approx(R_SE, abs=ATOL, rel=RTOL)
    assert res["levels"] == R_LEVELS


def test_default_reference_is_control():
    default = fast_ridit_analysis(W, Y)
    explicit = fast_ridit_analysis(W, Y, "control")
    assert default["estimate"] == pytest.approx(explicit["estimate"], abs=ATOL, rel=RTOL)


def test_reference_group_mean_ridit_is_half_by_construction():
    # A group's mean ridit relative to its OWN empirical distribution is
    # always exactly 0.5 -- true for whichever arm is chosen as reference.
    res_control = fast_ridit_analysis(W, Y, "control")
    assert res_control["mean_ridit_c"] == pytest.approx(0.5, abs=ATOL)

    res_treatment = fast_ridit_analysis(W, Y, "treatment")
    assert res_treatment["mean_ridit_t"] == pytest.approx(0.5, abs=ATOL)


def test_estimate_is_treatment_mean_ridit_minus_half():
    res = fast_ridit_analysis(W, Y, "control")
    assert res["estimate"] == pytest.approx(res["mean_ridit_t"] - 0.5, abs=ATOL, rel=RTOL)


def test_empty_reference_group_returns_nan_fields():
    # Unlike R's fast_ridit_analysis_cpp (which returns an empty list when
    # the reference group is empty), this binding always returns a
    # fully-populated dict -- NaN scalars and empty list fields -- matching
    # the NaN-populated-dict convention this session's other bindings use
    # for degenerate inputs (e.g. mn_ci at n=0), rather than an empty dict.
    w_all_treat = np.ones_like(W)
    res = fast_ridit_analysis(w_all_treat, Y, "control")
    assert math.isnan(res["mean_ridit_t"])
    assert math.isnan(res["mean_ridit_c"])
    assert math.isnan(res["estimate"])
    assert math.isnan(res["se"])
    assert res["scores"] == []
    assert res["levels"] == []
    assert res["ref_p"] == []
