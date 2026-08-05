"""Parity test for edi_kernels.wilcox_hl_point_estimate against an
R-generated fixture.

This does NOT re-verify the statistical correctness of the Hodges-Lehmann
estimator itself -- that's the same hl_from_groups object code the R
package's own test suite already exercises via
EDI:::wilcox_hl_point_estimate_cpp (both call the identical hl_from_groups
in EDI/src/fast_wilcox_hl.cpp, re-exposed under external linkage as
wilcox_hl_point_estimate_result so this binding's translation unit can call
it). What this test covers: the pybind11 binding layer itself and the
EDI_CORE_ONLY compiled path.

The expected value below was computed once via:
    EDI:::wilcox_hl_point_estimate_cpp(w, y)
in R on a synthetic n=40 dataset generated with set.seed(789).
"""
import numpy as np
import pytest

from edi_kernels import wilcox_hl_point_estimate

ATOL = 1e-9
RTOL = 1e-9

Y = np.array([
    6.048193421526, 0.47846423155274, 4.96064056382295, 5.36627978723278,
    4.27729703654138, 4.03103201572166, 3.66737410663191, 4.65106144208571,
    2.97808066532158, 6.47939211319428, 4.19540738662459, 2.99443935059055,
    4.64456158477191, 4.02420087236944, 6.85581478069316, 3.451149550331,
    5.84574224923325, 3.78607322234924, 5.41874995346491, 3.44532561986474,
    3.59589976010749, 6.36693834391202, 3.28462485435676, 5.73552167273596,
    2.14061594104969, 3.97535035486544, 4.46377925364063, 4.60153575528796,
    6.71314851504582, 4.66642127575769, 4.24865686884112, 2.73004625387865,
    6.21862950782745, 6.13171063321983, 3.52790282464051, 2.88289620617719,
    2.36427145227111, 4.59122042740955, 3.79722212455877, 6.32120523178529,
])
W = np.array([1, 0] * 20, dtype=np.int32)

# R reference, EDI 1.0.0, computed 2026-08-05 via:
#   EDI:::wilcox_hl_point_estimate_cpp(w, y)
R_HL = 0.0650472158359303


def test_matches_r_fixture():
    got = wilcox_hl_point_estimate(Y, W)
    assert got == pytest.approx(R_HL, abs=ATOL, rel=RTOL)


def test_matches_manual_median_of_pairwise_diffs():
    # Independent check against the textbook (not the fast bisection) HL
    # definition: median of all pairwise treatment-minus-control differences.
    yt = Y[W == 1]
    yc = Y[W == 0]
    diffs = (yt[:, None] - yc[None, :]).ravel()
    assert wilcox_hl_point_estimate(Y, W) == pytest.approx(np.median(diffs), abs=ATOL, rel=RTOL)


def test_empty_group_returns_nan():
    import math
    w_all_treat = np.ones_like(W)
    assert math.isnan(wilcox_hl_point_estimate(Y, w_all_treat))


def test_non_finite_y_dropped():
    y_with_nan = Y.copy()
    y_with_nan[0] = np.nan
    got_with_nan = wilcox_hl_point_estimate(y_with_nan, W)
    got_dropped = wilcox_hl_point_estimate(Y[1:], W[1:])
    assert got_with_nan == pytest.approx(got_dropped, abs=ATOL, rel=RTOL)
