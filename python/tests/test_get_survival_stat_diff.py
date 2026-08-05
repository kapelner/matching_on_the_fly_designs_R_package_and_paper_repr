"""Parity test for edi_kernels.get_survival_stat_diff against an
R-generated fixture.

This does NOT re-verify the statistical correctness of the KM-median-
difference statistic itself -- that's the same algorithm (KM curve, R's
survival::quantile.survfit step-function median semantics) the R package's
own test suite already exercises via EDI:::get_survival_stat_diff (both
call the identical logic in EDI/src/fast_survival_stats.cpp, ported to
get_survival_stat_diff_result/get_survival_stat_for_group_result so this
binding's translation unit can call it directly without an SEXP/wrap()
round-trip). What this test covers: the pybind11 binding layer itself and
the EDI_CORE_ONLY compiled path.

The expected value below was computed once via:
    EDI:::get_survival_stat_diff(y, dead, w, "median")
in R on a synthetic n=50 dataset generated with set.seed(321).
"""
import math

import numpy as np
import pytest

from edi_kernels import get_survival_stat_diff

ATOL = 1e-9
RTOL = 1e-9

Y = np.array([
    0.330243668478817, 1.42688344410601, 2.51039010464194, 2.11573323599043,
    1.78822659557811, 1.22638287954032, 0.422884277999401, 2.28902306213024,
    1.06535389833152, 6.43074203994045, 0.368687500245869, 3.99870373664384,
    0.530952102504671, 2.61501431829276, 1.71032162747326, 0.561836042441428,
    0.543871022760868, 0.673052502288582, 1.00495299510658, 0.172406522503754,
    4.14066198078876, 0.841013251803815, 0.703319482505322, 0.608175668979316,
    8.2061071801455, 0.689319854602218, 0.163812545128167, 1.28198249172419,
    1.59524546555224, 2.54871021713012, 2.12777091740611, 3.04149094952231,
    0.244976702881798, 0.15137149207294, 1.63156884155216, 8.84666130641356,
    5.25673054506011, 0.647397911176085, 0.709341395646334, 2.87318276969676,
    1.48813046692851, 1.92634124859216, 0.706980625167489, 4.17351907093949,
    4.22120489527412, 0.64476263243705, 5.51034237046666, 0.948096717707813,
    4.56539801855512, 0.608004671898981,
])
DEAD = np.array([
    1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
    1, 1,
], dtype=np.int32)
W = np.array([1, 0] * 25, dtype=np.int32)

# R reference, EDI 1.0.0, computed 2026-08-05 via:
#   EDI:::get_survival_stat_diff(y, dead, w, "median")
R_KM_DIFF = 0.10673367511481


def test_matches_r_fixture():
    got = get_survival_stat_diff(Y, DEAD, W, "median")
    assert got == pytest.approx(R_KM_DIFF, abs=ATOL, rel=RTOL)


def test_default_requested_stat_is_median():
    default = get_survival_stat_diff(Y, DEAD, W)
    explicit = get_survival_stat_diff(Y, DEAD, W, "median")
    assert default == pytest.approx(explicit, abs=ATOL, rel=RTOL)


def test_zero_diff_when_groups_identical():
    # Same data in both arms (interleaved so each unique time appears in
    # both groups) must give a KM-median difference of exactly zero.
    n = len(Y) // 2
    y_dup = np.concatenate([Y[:n], Y[:n]])
    dead_dup = np.concatenate([DEAD[:n], DEAD[:n]])
    w_dup = np.array([1] * n + [0] * n, dtype=np.int32)
    got = get_survival_stat_diff(y_dup, dead_dup, w_dup, "median")
    assert got == pytest.approx(0.0, abs=1e-9)


def test_empty_group_returns_nan():
    w_all_treat = np.ones_like(W)
    assert math.isnan(get_survival_stat_diff(Y, DEAD, w_all_treat, "median"))
