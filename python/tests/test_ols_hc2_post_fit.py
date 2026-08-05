"""Parity test for edi_kernels.ols_hc2_post_fit against an R-generated
fixture.

This does NOT re-verify the statistical correctness of the HC2 sandwich
standard error itself -- that's the same summarize_with_vcov object code
the R package's own test suite already exercises via
EDI:::ols_hc2_post_fit_cpp (both call the identical bread/hat/omega/meat
algorithm; only summarize_with_vcov's error path changed from Rcpp::stop()
to std::invalid_argument and its return type from Rcpp::List to a plain
struct, and ols_hc2_post_fit_result fuses ols_hc2_setup_cpp +
ols_hc2_post_fit_precomputed_cpp into one call avoiding
ols_hc2_post_fit_cpp's internal SEXP round-trip). What this test covers:
the pybind11 binding layer itself and the EDI_CORE_ONLY compiled path.

This kernel backs InferenceContinLin (Lin (2013) covariate-adjusted OLS):
X_fit there is [1, w, Xc, Xc*w] (intercept, treatment, centered covariates,
treatment-by-centered-covariate interactions), j_treat=2. The expected
values below were computed once via:
    coef_hat <- lm.fit(Xint, y)$coefficients
    EDI:::ols_hc2_post_fit_cpp(Xint, y, coef_hat, 2L)
in R on a synthetic n=60 dataset generated with set.seed(111).
"""
import numpy as np
import pytest

from edi_kernels import fast_ols, ols_hc2_post_fit

ATOL = 1e-9
RTOL = 1e-9

Y = np.array([
    2.0960997839798, 2.25908515607663, 0.614002522497192, 0.330797252285849,
    2.25839380868223, 0.525359000780545, 1.93967267870741, 1.19211574500912,
    1.28620613746427, 0.748141139993753, 3.05770353773665, 1.53483804804903,
    1.069925620492, 1.02505735131063, 2.06833894033676, -0.273124113677753,
    -0.725783324344603, 0.584658553951541, 1.23759849938716, 0.608901717527366,
    2.33183728305345, 1.7066242321227, 1.04187937695559, 1.31783414397417,
    2.4553509456443, 0.988430865724889, 2.12872861647617, 1.74806449095798,
    1.59033222456738, 1.01435759444172, 0.0918207217461808, -0.692582920995227,
    1.8014601059065, 0.514071232118755, -0.998673744778482, 1.56349264422455,
    0.914204812445817, 1.43585780043915, 0.628154323055885, 1.28781663087841,
    1.28051085951854, 0.564628840247903, 1.66503555962024, 1.73784819948674,
    2.19266934893766, 1.46632247421565, 2.16746857543284, 1.70535381550428,
    3.08657209032454, 1.03892841370038, 1.9172047047306, 1.38344279742332,
    1.50457424253425, 1.02292772096137, 1.1034800752667, 0.764274121615753,
    2.5142500727475, 1.89615349154569, 1.99281797347785, 0.518341829982467,
])
W = np.array([1, 0] * 30, dtype=float)
X1 = np.array([
    0.235220711613698, -0.330735871626691, -0.311623823978762, -2.30234565842938,
    -0.170876044613865, 0.140278225103231, -1.49742665556512, -1.01018841905138,
    -0.948475604952447, -0.493962217235108, -0.17367412797058, -0.406598780040652,
    1.84563626380708, 0.394054109974501, 0.797528501261517, -1.56666536018558,
    -0.0858510088226061, -0.359139481237783, -1.19360896656373, 0.364186737022221,
    0.361662451588275, 0.346964370113566, 0.189736526689275, -0.159576805538779,
    0.326549238431867, 0.59825420239349, -1.84153430006245, 2.71805559954175,
    0.191244390457294, -1.30129606527832, -3.11321730129512, -0.94135739527934,
    1.40025878166735, -1.62047002834908, -2.26599595710602, 1.16299358918063,
    -0.116155040696488, 0.334256009374765, -0.620858106369138, -1.30984490848087,
    -1.175726041705, -1.12121553281381, -1.36190448365479, 0.481124578468714,
    0.741971625532911, 0.0278246252153275, 0.331379710818421, 0.644114131991377,
    2.48566156143796, 1.95998170769623, 0.19166338489197, 1.55254427210515,
    0.914242286889689, 0.358625374801395, 0.175095636090836, -0.847267768953972,
    0.978231657400864, 1.80586825914887, 0.122914802352916, -0.1297720262759,
])
X2 = np.array([
    -0.216428659149337, 1.44647816735024, 0.409709801556001, 0.910916571988728,
    1.43035816657513, -0.381291955622273, 0.202307176881315, -0.806199194081451,
    0.294634183706975, 1.40488308312069, 1.02376684645418, 0.476126064485472,
    -0.670330329898876, 0.159234321208676, -0.382715383134029, 0.935762593244525,
    -0.631532267741634, -0.0983060774786364, 1.03198498311984, 0.387808429457905,
    -1.25612931333732, -0.786952730090235, 0.429811548667884, -0.37641621920597,
    -1.21622906636522, 1.02927851411841, 0.430396999024907, -1.24557402009338,
    -0.602728486487283, 0.660069387382248, 2.05074952611361, 0.490808178818426,
    -1.73147941878002, 0.71088365940818, 0.0138229114622208, -1.40104159736621,
    1.25912366544167, -0.127477518375642, -0.729386513651465, -1.21136136021004,
    0.599619744841827, -1.16032952722025, 0.439093434435602, 0.204853736631967,
    -0.699181340101498, -0.926625676673078, -1.01348237524073, 0.604987059771553,
    1.73440013278809, -0.349850530637584, 1.19918550882917, 1.02390099805067,
    -0.0462751493457263, 1.3865084257746, -2.19527271893927, -0.127556912959721,
    -1.31699684187755, -1.42693514077612, 0.879441630777448, 0.393991575061551,
])

X = np.column_stack([X1, X2])
XC = X - X.mean(axis=0)
X_FIT = np.column_stack([np.ones(len(Y)), W, XC, XC * W[:, None]])

R_BETA_HAT = 0.525943573983647
R_SE = 0.181738156858422
R_SSQ_HAT = 0.0330287576582965


def test_matches_r_fixture():
    fit = fast_ols(X_FIT, Y, estimate_only=True)
    coef_hat = np.asarray(fit["b"])
    post = ols_hc2_post_fit(X_FIT, Y, coef_hat, j_treat=2)
    assert post["beta_hat"] == pytest.approx(R_BETA_HAT, abs=ATOL, rel=RTOL)
    assert post["se"] == pytest.approx(R_SE, abs=ATOL, rel=RTOL)
    assert post["ssq_hat"] == pytest.approx(R_SSQ_HAT, abs=ATOL, rel=RTOL)


def test_default_j_treat_is_2():
    fit = fast_ols(X_FIT, Y, estimate_only=True)
    coef_hat = np.asarray(fit["b"])
    default = ols_hc2_post_fit(X_FIT, Y, coef_hat)
    explicit = ols_hc2_post_fit(X_FIT, Y, coef_hat, j_treat=2)
    assert default["beta_hat"] == pytest.approx(explicit["beta_hat"], abs=ATOL, rel=RTOL)


def test_beta_hat_matches_coef_hat_at_j_treat():
    fit = fast_ols(X_FIT, Y, estimate_only=True)
    coef_hat = np.asarray(fit["b"])
    post = ols_hc2_post_fit(X_FIT, Y, coef_hat, j_treat=2)
    assert post["beta_hat"] == pytest.approx(coef_hat[1], abs=ATOL, rel=RTOL)


def test_result_shapes():
    fit = fast_ols(X_FIT, Y, estimate_only=True)
    coef_hat = np.asarray(fit["b"])
    post = ols_hc2_post_fit(X_FIT, Y, coef_hat, j_treat=2)
    p = X_FIT.shape[1]
    assert post["vcov"].shape == (p, p)
    assert post["std_err"].shape == (p,)
    assert post["z_vals"].shape == (p,)


def test_se_is_sqrt_ssq_hat():
    fit = fast_ols(X_FIT, Y, estimate_only=True)
    coef_hat = np.asarray(fit["b"])
    post = ols_hc2_post_fit(X_FIT, Y, coef_hat, j_treat=2)
    assert post["se"] == pytest.approx(np.sqrt(post["ssq_hat"]), abs=ATOL, rel=RTOL)


def test_out_of_bounds_j_treat_raises():
    fit = fast_ols(X_FIT, Y, estimate_only=True)
    coef_hat = np.asarray(fit["b"])
    with pytest.raises(Exception):
        ols_hc2_post_fit(X_FIT, Y, coef_hat, j_treat=999)
