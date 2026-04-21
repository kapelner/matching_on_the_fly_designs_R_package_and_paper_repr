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
    int    start;      // index into sorted-obs arrays
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

    LMMData(const Eigen::VectorXd& y,
            const Eigen::MatrixXd& X,
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
                info.start = i;
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
                       const Eigen::VectorXd& par,
                       Eigen::VectorXd& grad)
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
        const int s = g.start;
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

// ── Wrapper satisfying LBFGSpp operator() signature ─────────────────────────
class GaussianLMMObjective {
public:
    const LMMData& dat;
    explicit GaussianLMMObjective(const LMMData& d) : dat(d) {}

    double operator()(const Eigen::VectorXd& par, Eigen::VectorXd& grad) {
        double nll = neg_ll_and_grad(dat, par, grad);
        if (!std::isfinite(nll)) {
            grad.setZero();
            return 1e300;
        }
        return nll;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& par) const {
        const int k = par.size();
        Eigen::MatrixXd H(k, k);
        const double h_rel = 1e-4;

        for (int j = 0; j < k; ++j) {
            const double h = h_rel * std::max(1.0, std::abs(par[j]));
            Eigen::VectorXd pp = par, pm = par;
            pp[j] += h;
            pm[j] -= h;
            Eigen::VectorXd gp(k), gm(k);
            neg_ll_and_grad(dat, pp, gp);
            neg_ll_and_grad(dat, pm, gm);
            H.col(j) = (gp - gm) / (2.0 * h);
        }
        return (H + H.transpose()) / 2.0;
    }
};

// ── Numerical Hessian of neg_ll (for Fisher info / vcov) ────────────────────
Eigen::MatrixXd numerical_hessian(const LMMData& dat,
                                  const Eigen::VectorXd& par,
                                  double h_rel = 1e-4)
{
    const int k = par.size();
    Eigen::MatrixXd H(k, k);
    Eigen::VectorXd dummy(k);

    for (int j = 0; j < k; ++j) {
        const double h = h_rel * std::max(1.0, std::abs(par[j]));
        Eigen::VectorXd pp = par, pm = par;
        pp[j] += h;  pm[j] -= h;
        Eigen::VectorXd gp(k), gm(k);
        neg_ll_and_grad(dat, pp, gp);
        neg_ll_and_grad(dat, pm, gm);
        H.col(j) = (gp - gm) / (2.0 * h);
    }
    return (H + H.transpose()) / 2.0;
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

    Eigen::VectorXd start(p + 2);
    start.head(p)   = beta;
    start[p]     = std::log(sigma_e);        // log σ_e
    start[p + 1] = std::log(sigma_e * 0.5);  // log σ_b ≈ σ_e/2
    return start;
}

} // anonymous namespace

// ── R-exported: fit Gaussian LMM ─────────────────────────────────────────────
// [[Rcpp::export]]
List fast_gaussian_lmm_cpp(
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& X,       // n × p, intercept in col 0, treatment in col 1
    const Eigen::VectorXi& group_id, // 1-based group IDs (length n)
    Nullable<NumericVector> start_par = R_NilValue,
    bool  estimate_only = false,
    int   maxit  = 300,
    double eps_g = 1e-6,
    Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
    std::string optimization_alg = "lbfgs"
) {
    const int n = y.size(), p = X.cols();

    // Convert group_id to 0-based sorted integer IDs
    std::vector<int> gid(n);
    {
        // Map R group ids (any positive integers) to 0-based consecutive ints
        std::vector<int> gid_r(group_id.data(), group_id.data() + n);
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
    if (start_par.isNotNull()) {
        NumericVector sp(start_par);
        if (sp.size() == p + 2)
            for (int i = 0; i < p + 2; ++i) par[i] = sp[i];
    }
    FixedParamSpec fixed_spec = make_fixed_param_spec(p + 2, fixed_idx, fixed_values);

    double neg_ll = 1e300;
    int niter = maxit;
    bool converged = false;
    try {
        LikelihoodFitResult fit = optimize_fixed_likelihood(obj, par, fixed_spec, maxit, eps_g, optimization_alg, "lbfgs", 100);
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
    Eigen::MatrixXd H = numerical_hessian(dat, par);
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
        Named("ssq_b_T")   = ssq_b_T,
        Named("vcov")      = vcov_r,
        Named("neg_loglik")= neg_ll,
        Named("converged") = converged,
        Named("niter")     = niter
    );
}

// ── R-exported: score (gradient of log_lik) at arbitrary par ─────────────────
// [[Rcpp::export]]
NumericVector get_gaussian_lmm_score_cpp(
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& X,
    const Eigen::VectorXi& group_id,
    const Eigen::VectorXd& par
) {
    const int n = y.size();
    std::vector<int> gid(n);
    {
        std::vector<int> gid_r(group_id.data(), group_id.data() + n);
        std::vector<int> uniq = gid_r;
        std::sort(uniq.begin(), uniq.end());
        uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
        for (int i = 0; i < n; ++i)
            gid[i] = (int)(std::lower_bound(uniq.begin(), uniq.end(), gid_r[i]) - uniq.begin());
    }
    LMMData dat(y, X, gid);
    Eigen::VectorXd grad;
    neg_ll_and_grad(dat, par, grad);
    // score = -grad of neg_ll
    Eigen::VectorXd score = -grad;
    return wrap(score);
}

// ── R-exported: observed Fisher information (Hessian of neg_ll) at par ───────
// [[Rcpp::export]]
NumericMatrix get_gaussian_lmm_fisher_cpp(
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& X,
    const Eigen::VectorXi& group_id,
    const Eigen::VectorXd& par,
    double h_rel = 1e-4
) {
    const int n = y.size();
    std::vector<int> gid(n);
    {
        std::vector<int> gid_r(group_id.data(), group_id.data() + n);
        std::vector<int> uniq = gid_r;
        std::sort(uniq.begin(), uniq.end());
        uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
        for (int i = 0; i < n; ++i)
            gid[i] = (int)(std::lower_bound(uniq.begin(), uniq.end(), gid_r[i]) - uniq.begin());
    }
    LMMData dat(y, X, gid);
    Eigen::MatrixXd H = numerical_hessian(dat, par, h_rel);
    return wrap(H);
}
