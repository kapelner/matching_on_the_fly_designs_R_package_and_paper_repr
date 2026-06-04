// Gaussian Linear Mixed Model for KK designs
// Model: y_ij = X_ij β + u_i + ε_ij
//   u_i ~ N(0, σ²_b)  (random intercept per matched pair / singleton)
//   ε_ij ~ N(0, σ²_e) (residual)
//
// Parameter vector: par = (β[0..p-1], log_σ_e, log_σ_b)
//   log_σ_e = log(σ_e),  log_σ_b = log(σ_b)
//
// Groups: matched pairs have size 2; reservoir subjects have size 1.

#include "_helper_functions.h"
#include <RcppEigen.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>

using namespace Rcpp;

static const double LOG2PI = std::log(2.0 * M_PI);

namespace {

// ── Pre-computed, design-fixed quantities per group ─────────────────────────
struct GroupInfo {
    int    size;       // m_i  (1 or 2)
    int    warm_start_params;      // index into sorted-obs arrays
};

// ── Core data object (passed by const-ref to all functions) ──────────────────
struct LMMData {
    // Observations sorted so that obs within each group are contiguous
    Eigen::VectorXd y_s;         // n sorted responses
    Eigen::MatrixXd X_s;         // n × p sorted design (includes intercept at col 0)
    std::vector<GroupInfo> grps; // one entry per group
    int n;
    int p;
    int G;  // number of groups

    LMMData(const Eigen::Ref<const VectorXd>& y,
            const Eigen::Ref<const MatrixXd>& X,
            const std::vector<int>& gid)   // 0-based group ids, length n
        : n(y.size()), p(X.cols())
    {
        // Sort observations by group id
        std::vector<int> ord(n);
        std::iota(ord.begin(), ord.end(), 0);
        std::stable_sort(ord.begin(), ord.end(),
                         [&](int a, int b){ return gid[a] < gid[b]; });

        y_s.resize(n);
        X_s.resize(n, p);
        for (int i = 0; i < n; ++i) {
            y_s[i]   = y[ord[i]];
            X_s.row(i) = X.row(ord[i]);
        }

        // Build group list
        int prev = -1, gi = 0;
        for (int i = 0; i < n; ++i) {
            int g = gid[ord[i]];
            if (g != prev) {
                GroupInfo info;
                info.warm_start_params = i;
                info.size  = 1;
                grps.push_back(info);
                prev = g;
                gi++;
            } else {
                grps.back().size++;
            }
        }
        G = (int)grps.size();
    }
};

// ── Neg log-likelihood with analytical gradient ──────────────────────────────
//
// neg_ll = (n/2)log(2π)
//        + Σ_g [ (m_g-1)/2 · log(v_e) + 1/2 · log(a_g) ]
//        + 1/(2 v_e) · Σ_g [ Q_g - (v_b/a_g)·S_g² ]
//
// where a_g = v_e + m_g · v_b,  S_g = Σ r,  Q_g = Σ r²,  r = y - Xβ
//
// Analytical gradients:
//   ∂/∂β       = -(1/v_e)·X^T·r + (v_b/v_e) Σ_g (S_g/a_g)·Σ_{j∈g} x_j
//   ∂/∂log_σ_e = Σ_g(m_g-1) + v_e·Σ_g(1/a_g) - (1/v_e)·Σ_g W_g
//   ∂/∂log_σ_b = Σ_g(m_g·v_b/a_g) - v_b·Σ_g(S_g²/a_g²)
//
// where W_g = Q_g - (v_b/a_g)·S_g²

double neg_ll_and_grad(const LMMData& dat,
                       const Eigen::Ref<const Eigen::VectorXd>& par,
                       Eigen::Ref<VectorXd> grad)
{
    const int p = dat.p, n = dat.n;
    const Eigen::VectorXd beta = par.head(p);
    const double lse = par[p];       // log σ_e
    const double lsb = par[p + 1];   // log σ_b

    const double v_e = std::exp(2.0 * lse);
    const double v_b = std::exp(2.0 * lsb);

    if (!std::isfinite(v_e) || !std::isfinite(v_b) || v_e < 1e-300)
        return 1e300;

    // Residuals
    Eigen::VectorXd r = dat.y_s - dat.X_s * beta;

    // Accumulators for the two variance-parameter gradients
    double d_lse = 0.0;   // will add contributions below
    double d_lsb = 0.0;

    // Accumulate: weighted residual for β gradient
    // ∂/∂β contribution from group g:
    //   -(1/v_e)·r_g  +  (v_b/v_e)·(S_g/a_g)·1_g
    // We accumulate into grad_beta using the sorted obs.

    Eigen::VectorXd grad_beta = Eigen::VectorXd::Zero(p);

    double neg_ll = (n * 0.5) * LOG2PI;

    for (int gi = 0; gi < dat.G; ++gi) {
        const GroupInfo& g = dat.grps[gi];
        const int m = g.size;
        const int s = g.warm_start_params;
        const double a = v_e + m * v_b;

        double S = 0.0, Q = 0.0;
        for (int j = s; j < s + m; ++j) {
            S += r[j];
            Q += r[j] * r[j];
        }

        const double W      = Q - (v_b / a) * S * S;
        const double inv_a  = 1.0 / a;

        // --- neg_ll contribution ---
        neg_ll += 0.5 * ((m - 1) * 2.0 * lse + std::log(a));  // det term
        neg_ll += W / (2.0 * v_e);                              // quadratic term

        // --- gradient for β ---
        const double coeff_S = (v_b / (v_e * a)) * S;
        for (int j = s; j < s + m; ++j) {
            const double w_j = -r[j] / v_e + coeff_S;
            // grad_beta -= w_j * x_j  (remember: neg_ll, so ∂neg_ll/∂β = -score)
            grad_beta -= w_j * dat.X_s.row(j).transpose();
        }

        // --- gradient for log σ_e ---
        // ∂neg_ll/∂lse = Σ_g [(m_g-1) + v_e/a_g + v_b·S_g²/a_g² - W_g/v_e]
        // Derivation: ∂/∂lse of [(m-1)·lse + ½log(a) + W/(2v_e)]
        //   det:  (m-1) + v_e/a   (chain ∂v_e/∂lse = 2v_e)
        //   quad: v_b·S²/a²  (from ∂W/∂lse through a) − W/v_e
        d_lse += (m - 1) + v_e * inv_a + v_b * (S * S) * (inv_a * inv_a) - W / v_e;

        // --- gradient for log σ_b ---
        // ∂neg_ll/∂lsb = Σ_g(m·v_b/a) - v_b·Σ_g(S²/a²)
        d_lsb += (m * v_b * inv_a) - v_b * (S * S) * (inv_a * inv_a);
    }

    // Fix the gradient sign for β (we subtracted above, which was the score; now it's already
    // the neg_ll gradient direction after the sign flip in the loop)
    grad.resize(p + 2);
    grad.head(p) = -grad_beta;   // grad_beta accumulated score; neg_ll grad = -score
    grad[p]     = d_lse;
    grad[p + 1] = d_lsb;

    return neg_ll;
}

Eigen::MatrixXd lmm_analytic_hessian(const LMMData& dat,
                                     const Eigen::Ref<const Eigen::VectorXd>& par)
{
    const int p = dat.p;
    const int k = (int)par.size();
    const Eigen::VectorXd beta = par.head(p);
    const double v_e = std::exp(2.0 * par[p]);
    const double v_b = std::exp(2.0 * par[p + 1]);

    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(k, k);
    if (!std::isfinite(v_e) || !std::isfinite(v_b) || v_e < 1e-300) {
        H.setConstant(NA_REAL);
        return H;
    }

    const Eigen::VectorXd r_all = dat.y_s - dat.X_s * beta;

    for (int gi = 0; gi < dat.G; ++gi) {
        const GroupInfo& g = dat.grps[gi];
        const int m = g.size;
        const int s = g.warm_start_params;
        const Eigen::MatrixXd Xg = dat.X_s.middleRows(s, m);
        const Eigen::VectorXd rg = r_all.segment(s, m);

        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(m, m);
        Eigen::MatrixXd J = Eigen::MatrixXd::Ones(m, m);
        Eigen::MatrixXd V = v_e * I + v_b * J;
        Eigen::LDLT<Eigen::MatrixXd> ldlt(V);
        if (ldlt.info() != Eigen::Success) {
            H.setConstant(NA_REAL);
            return H;
        }
        Eigen::MatrixXd P = ldlt.solve(I);

        Eigen::MatrixXd dV_e = 2.0 * v_e * I;
        Eigen::MatrixXd dV_b = 2.0 * v_b * J;
        Eigen::MatrixXd d2V_ee = 4.0 * v_e * I;
        Eigen::MatrixXd d2V_bb = 4.0 * v_b * J;

        Eigen::MatrixXd dP_e = -P * dV_e * P;
        Eigen::MatrixXd dP_b = -P * dV_b * P;

        H.topLeftCorner(p, p).noalias() += Xg.transpose() * P * Xg;
        H.topRightCorner(p, 1).noalias() += -Xg.transpose() * dP_e * rg;
        H.block(0, p + 1, p, 1).noalias() += -Xg.transpose() * dP_b * rg;

        Eigen::MatrixXd dV[2] = {dV_e, dV_b};
        Eigen::MatrixXd dP[2] = {dP_e, dP_b};
        Eigen::MatrixXd d2V[2][2] = {
            {d2V_ee, Eigen::MatrixXd::Zero(m, m)},
            {Eigen::MatrixXd::Zero(m, m), d2V_bb}
        };

        for (int a = 0; a < 2; ++a) {
            for (int b = a; b < 2; ++b) {
                Eigen::MatrixXd term_mat =
                    dP[a] * dV[b] * P +
                    P * d2V[a][b] * P +
                    P * dV[b] * dP[a];
                const double h_ab =
                    0.5 * (dP[a] * dV[b] + P * d2V[a][b]).trace()
                    - 0.5 * (rg.transpose() * term_mat * rg)(0, 0);
                H(p + a, p + b) += h_ab;
                if (a != b) H(p + b, p + a) += h_ab;
            }
        }
    }

    H.bottomLeftCorner(2, p) = H.topRightCorner(p, 2).transpose();
    return (H + H.transpose()) / 2.0;
}

// ── Wrapper satisfying LBFGSpp operator() signature ─────────────────────────
class GaussianLMMObjective {
public:
    const LMMData& dat;
    explicit GaussianLMMObjective(const LMMData& d) : dat(d) {}

    double operator()(const Eigen::VectorXd& par, Eigen::VectorXd& grad) {
        double nll = neg_ll_and_grad(dat, par, grad);
        if (!std::isfinite(nll) || nll >= 1e299) {
            grad.setZero();
            return 1e300;
        }
        return nll;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& par) const {
        return lmm_analytic_hessian(dat, par);
    }
};

// ── Hessian of neg_ll (for Fisher info / vcov) ─────────────────────────────
Eigen::MatrixXd lmm_fisher_hessian(const LMMData& dat,
                                   const Eigen::Ref<const Eigen::VectorXd>& par,
                                   double h_rel = 1e-4)
{
    (void)h_rel;
    return lmm_analytic_hessian(dat, par);
}

// ── Starting values: OLS β, residual-based σ_e, σ_b = σ_e/2 ─────────────────
Eigen::VectorXd make_start(const LMMData& dat)
{
    const int p = dat.p;
    // OLS: β = (X^T X)^{-1} X^T y
    Eigen::MatrixXd XtX = dat.X_s.transpose() * dat.X_s;
    Eigen::VectorXd Xty = dat.X_s.transpose() * dat.y_s;
    Eigen::VectorXd beta = XtX.ldlt().solve(Xty);

    // Residual SD
    Eigen::VectorXd res = dat.y_s - dat.X_s * beta;
    double sigma2_e = res.squaredNorm() / std::max(dat.n - p, 1);
    double sigma_e  = std::sqrt(std::max(sigma2_e, 1e-8));

    Eigen::VectorXd warm_start_params(p + 2);
    warm_start_params.head(p)   = beta;
    warm_start_params[p]     = std::log(sigma_e);        // log σ_e
    warm_start_params[p + 1] = std::log(sigma_e * 0.5);  // log σ_b ≈ σ_e/2
    return warm_start_params;
}

} // anonymous namespace

// ── R-exported: fit Gaussian LMM ─────────────────────────────────────────────
// [[Rcpp::export]]
List fast_gaussian_lmm_cpp(
    SEXP X_r,       // n × p, intercept in col 0, treatment in col 1
    SEXP y_r,
    SEXP group_id_r, // 1-based group IDs (length n)
    Rcpp::Nullable<Rcpp::NumericVector> warm_start_params = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
    bool  estimate_only = false,
    int   maxit  = 300,
    double eps_g = 1e-6,
    Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
    std::string optimization_alg = "lbfgs",
    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue
) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    IntegerVector group_id_int(group_id_r);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    Eigen::Map<const Eigen::VectorXi> group_id(group_id_int.begin(), group_id_int.size());

    const int n = y.size(), p = X.cols();

    // Convert group_id to 0-based sorted integer IDs
    std::vector<int> gid(n);
    {
        // Map R group ids (any positive integers) to 0-based consecutive ints
        const int* gid_ptr = group_id.data();
        std::vector<int> gid_r(gid_ptr, gid_ptr + n);
        std::vector<int> uniq = gid_r;
        std::sort(uniq.begin(), uniq.end());
        uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
        for (int i = 0; i < n; ++i)
            gid[i] = (int)(std::lower_bound(uniq.begin(), uniq.end(), gid_r[i]) - uniq.begin());
    }

    LMMData dat(y, X, gid);
    GaussianLMMObjective obj(dat);

    // Starting point
    Eigen::VectorXd par = make_start(dat);
    if (warm_start_params.isNotNull()) {
        NumericVector sp(warm_start_params);
        if (sp.size() == p + 2)
            for (int i = 0; i < p + 2; ++i) par[i] = sp[i];
    } else if (warm_start_beta.isNotNull()) {
        VectorXd sb = as<VectorXd>(warm_start_beta);
        if (sb.size() == p + 2) {
            par = sb;
        } else if (sb.size() == p) {
            par.head(p) = sb;
        }
    }
    FixedParamSpec fixed_spec = make_fixed_param_spec(p + 2, fixed_idx, fixed_values);

    Eigen::MatrixXd info_start;
    const Eigen::MatrixXd* info_start_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        info_start_ptr = &info_start;
    }

    double neg_ll = 1e300;
    int niter = maxit;
    bool converged = false;
    try {
        LikelihoodFitResult fit = optimize_fixed_likelihood(obj, par, fixed_spec, maxit, eps_g, optimization_alg, "lbfgs", 0, info_start_ptr);
        par = fit.params;
        neg_ll = fit.value;
        niter = fit.niter;
        converged = std::isfinite(neg_ll) && fit.converged;
    } catch (...) {
        converged = false;
    }

    // Return par names: β[0..p-1], log_sigma_e, log_sigma_b
    NumericVector b_r(p + 2);
    for (int i = 0; i < p + 2; ++i) b_r[i] = par[i];
    CharacterVector b_names(p + 2);
    for (int i = 0; i < p; ++i) b_names[i] = "b" + std::to_string(i);
    b_names[p]   = "log_sigma_e";
    b_names[p+1] = "log_sigma_b";
    b_r.names() = b_names;

    if (estimate_only) {
        return List::create(
            Named("b")         = b_r,
            Named("ssq_b_T")   = NA_REAL,
            Named("neg_loglik")= neg_ll,
            Named("converged") = converged,
            Named("niter")     = niter
        );
    }

    // Variance-covariance via Hessian of neg_ll
    Eigen::MatrixXd H = lmm_fisher_hessian(dat, par);
    Eigen::MatrixXd H_free = subset_matrix(H, fixed_spec.free_idx, fixed_spec.free_idx);
    Eigen::LDLT<Eigen::MatrixXd> ldlt(H_free);

    double ssq_b_T = NA_REAL;
    NumericMatrix vcov_r(p + 2, p + 2);
    std::fill(vcov_r.begin(), vcov_r.end(), NA_REAL);

    if (ldlt.info() == Eigen::Success) {
        Eigen::MatrixXd V_free = ldlt.solve(Eigen::MatrixXd::Identity(H_free.rows(), H_free.cols()));
        Eigen::MatrixXd V = expand_free_covariance(p + 2, fixed_spec, V_free, true);
        if (V.allFinite()) {
            ssq_b_T = V(1, 1);   // treatment is index 1 (after intercept at 0)
            for (int i = 0; i < p + 2; ++i)
                for (int j = 0; j < p + 2; ++j)
                    vcov_r(i, j) = V(i, j);
        }
    }

    return List::create(
        Named("b")         = b_r,
        Named("params")    = b_r,
        Named("ssq_b_T")   = ssq_b_T,
        Named("vcov")      = vcov_r,
        Named("neg_loglik")= neg_ll,
        Named("converged") = converged,
        Named("niter")     = niter,
        Named("fisher_information") = H
    );
}

// ── R-exported: score (gradient of log_lik) at arbitrary par ─────────────────
// [[Rcpp::export]]
NumericVector get_gaussian_lmm_score_cpp(
    SEXP X_r,
    SEXP y_r,
    SEXP group_id_r,
    SEXP par_sexp
) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    IntegerVector group_id_int(group_id_r);
    NumericVector par_r(par_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    Eigen::Map<const Eigen::VectorXi> group_id(group_id_int.begin(), group_id_int.size());
    Eigen::Map<const Eigen::VectorXd> par(par_r.begin(), par_r.size());

    const int n = y.size();
    std::vector<int> gid(n);
    {
        const int* gid_ptr = group_id.data();
        std::vector<int> gid_r(gid_ptr, gid_ptr + n);
        std::vector<int> uniq = gid_r;
        std::sort(uniq.begin(), uniq.end());
        uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
        for (int i = 0; i < n; ++i)
            gid[i] = (int)(std::lower_bound(uniq.begin(), uniq.end(), gid_r[i]) - uniq.begin());
    }
    LMMData dat(y, X, gid);
    const int k = X.cols() + 2;  // p betas + log_sigma_e + log_sigma_b
    if (par.size() != k) {
        return NumericVector(k, NA_REAL);
    }
    Eigen::VectorXd grad = Eigen::VectorXd::Zero(k);
    neg_ll_and_grad(dat, par, grad);
    // score = -grad of neg_ll
    Eigen::VectorXd score = -grad;
    return wrap(score);
}

// ── R-exported: observed Fisher information (Hessian of neg_ll) at par ───────
// [[Rcpp::export]]
NumericMatrix get_gaussian_lmm_fisher_cpp(
    SEXP X_r,
    SEXP y_r,
    SEXP group_id_r,
    SEXP par_sexp,
    double h_rel = 1e-4
) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    IntegerVector group_id_int(group_id_r);
    NumericVector par_r(par_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    Eigen::Map<const Eigen::VectorXi> group_id(group_id_int.begin(), group_id_int.size());
    Eigen::Map<const Eigen::VectorXd> par(par_r.begin(), par_r.size());

    const int n = y.size();
    std::vector<int> gid(n);
    {
        const int* gid_ptr = group_id.data();
        std::vector<int> gid_r(gid_ptr, gid_ptr + n);
        std::vector<int> uniq = gid_r;
        std::sort(uniq.begin(), uniq.end());
        uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
        for (int i = 0; i < n; ++i)
            gid[i] = (int)(std::lower_bound(uniq.begin(), uniq.end(), gid_r[i]) - uniq.begin());
    }
    LMMData dat(y, X, gid);
    Eigen::MatrixXd H = lmm_fisher_hessian(dat, par, h_rel);
    return wrap(H);
}
