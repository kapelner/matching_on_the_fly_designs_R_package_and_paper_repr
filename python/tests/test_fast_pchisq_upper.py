import numpy as np
import pytest
from scipy.stats import chi2

from edi_kernels import fast_pchisq_upper


@pytest.mark.parametrize("statistic,df", [
    (3.84, 1),
    (0.0, 1),
    (10.0, 3),
    (50.0, 20),
])
def test_matches_scipy_chi2_sf(statistic, df):
    got = fast_pchisq_upper(statistic, df)
    expected = chi2.sf(statistic, df)
    assert got == pytest.approx(expected, abs=1e-9, rel=1e-9)


def test_nonpositive_statistic_is_one():
    assert fast_pchisq_upper(0.0, 5) == 1.0
    assert fast_pchisq_upper(-1.0, 5) == 1.0


def test_invalid_df_is_nan():
    assert np.isnan(fast_pchisq_upper(1.0, 0.0))
    assert np.isnan(fast_pchisq_upper(1.0, -1.0))
