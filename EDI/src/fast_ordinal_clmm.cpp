#include <RcppEigen.h>
#include "_helper_functions.h"
#include "_glmm_links.h"
#include "_glmm_engine.h"

using namespace Rcpp;
using namespace glmm;

namespace {

// Ordinal model with thresholds
template <typename Link>
class OrdinalLikelihood {
    int K;
public:
    OrdinalLikelihood(int K_) : K(K_) {}

    int n_model_params() const { return K - 1; }

    inline void get_alpha_bounds(const Eigen::VectorXd& par, int y_ir, double& alpha_lo, double& alpha_up) const {
        alpha_lo = 0.0;
        alpha_up = 0.0;
        const int upper_idx = (y_ir >= K) ? -1 : (y_ir - 1);
        const int lower_idx = (y_ir <= 1) ? -1 : (y_ir - 2);
        const int last_idx = std::max(lower_idx, upper_idx);
        if (last_idx < 0) return;
        double alpha_curr = par[0];
        if (lower_idx == 0) alpha_lo = alpha_curr;
        if (upper_idx == 0) alpha_up = alpha_curr;
        for (int k = 1; k <= last_idx; ++k) {
            alpha_curr += std::exp(par[k]);
            if (k == lower_idx) alpha_lo = alpha_curr;
            if (k == upper_idx) alpha_up = alpha_curr;
        }
    }

    double log_prob(double y_val, double eta, const Eigen::VectorXd& par) const {
        int y_ir = static_cast<int>(y_val);
        double alpha_lo, alpha_up;
        get_alpha_bounds(par, y_ir, alpha_lo, alpha_up);
        double F_up = (y_ir >= K) ? 1.0 : Link::cdf(alpha_up - eta);
        double F_lo = (y_ir <= 1) ? 0.0 : Link::cdf(alpha_lo - eta);
        return std::log(std::max(1e-15, F_up - F_lo));
    }

    double log_prob_derivs(double y_val, double eta, const Eigen::VectorXd& par, double& de, Eigen::VectorXd& dp) const {
        int y_ir = static_cast<int>(y_val);
        int na = K - 1;

        double alpha_lo, alpha_up;
        get_alpha_bounds(par, y_ir, alpha_lo, alpha_up);

        double F_up = (y_ir >= K) ? 1.0 : Link::cdf(alpha_up - eta);
        double F_lo = (y_ir <= 1) ? 0.0 : Link::cdf(alpha_lo - eta);
        double prob = std::max(1e-15, F_up - F_lo);
        double lp = std::log(prob);

        // pdf_from_cdf avoids re-calling cdf() inside pdf() — F_up/F_lo already computed.
        const double g_up = (y_ir >= K) ? 0.0 : Link::pdf_from_cdf(alpha_up - eta, F_up) / prob;
        const double g_lo = (y_ir <= 1) ? 0.0 : Link::pdf_from_cdf(alpha_lo - eta, F_lo) / prob;

        dp.setZero(na);
        de = -(g_up - g_lo);
        if (na > 0) dp[0] = g_up - g_lo;

        if (y_ir <= 1) return lp;

        const int lower_idx = y_ir - 2;
        const double base_suffix = (y_ir >= K) ? -g_lo : (g_up - g_lo);
        for (int j = 1; j <= lower_idx; ++j) {
            dp[j] = base_suffix * std::exp(par[j]);
        }
        if (y_ir < K) {
            const int upper_idx = y_ir - 1;
            dp[upper_idx] = g_up * std::exp(par[upper_idx]);
        }
        return lp;
    }
};

template <typename Link>
List run_fast_ordinal_clmm(const GLMMData& dat, int K, int j_T, bool estimate_only, int maxit, double eps_g, 
                           const Eigen::Ref<const Eigen::VectorXd>& start_full, const std::string& optimization_alg,
                           const FixedParamSpec& fixed_spec,
                           const Eigen::MatrixXd* info_start_ptr = nullptr) {
    OrdinalLikelihood<Link> model(K);
    GLMMObjective<OrdinalLikelihood<Link>> obj(dat, model);
    
    Eigen::VectorXd par = start_full;
    double neg_ll = NA_REAL;
    bool converged = false;
    
    try {
        LikelihoodFitResult fit = optimize_fixed_likelihood(obj, par, fixed_spec, maxit, eps_g, optimization_alg, "lbfgs", 0, info_start_ptr);
        par = fit.params;
        neg_ll = fit.value;
        converged = std::isfinite(neg_ll) && fit.converged;
    } catch (...) {
        converged = false;
    }

    int total = par.size();
    int na = K - 1;
    int j_T_full = na + j_T;
    double ssq_b_T = NA_REAL;

    if (!estimate_only && converged) {
        FixedParameterFunctor<GLMMObjective<OrdinalLikelihood<Link>>> fixed_obj(obj, fixed_spec, par);
        Eigen::VectorXd params_free = subset_vector(par, fixed_spec.free_idx);
        Eigen::MatrixXd H_free = fixed_obj.hessian(params_free);
        Eigen::LDLT<Eigen::MatrixXd> ldlt(H_free);
        if (ldlt.info() == Eigen::Success) {
            Eigen::MatrixXd inv_free = ldlt.solve(Eigen::MatrixXd::Identity(H_free.rows(), H_free.cols()));
            if (inv_free.allFinite()) {
                Eigen::MatrixXd inv = expand_free_covariance(total, fixed_spec, inv_free, true);
                if (j_T_full < total) ssq_b_T = inv(j_T_full, j_T_full);
            }
        }
    }

    return List::create(
        Named("b")          = par.segment(na, dat.p),
        Named("alpha")      = par.head(na),
        Named("log_sigma")  = par[total - 1],
        Named("ssq_b_T")    = ssq_b_T,
        Named("converged")  = converged,
        Named("neg_loglik") = neg_ll
    );
}

} // namespace

// [[Rcpp::export]]
List fast_ordinal_clmm_cpp(
    const Rcpp::NumericMatrix& X,
    const Rcpp::IntegerVector& y,
    const Rcpp::IntegerVector& group_id,
    int K,
    int j_T,
    std::string link = "logit",
    bool estimate_only = false,
    int n_gh = 20,
    double max_abs_log_sigma = 8.0,
    int maxit = 300,
    double eps_g = 1e-6,
    Rcpp::Nullable<Rcpp::NumericVector> warm_start_params = R_NilValue,
    std::string optimization_alg = "lbfgs",
    Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue
) {
    Eigen::Map<const Eigen::MatrixXd> map_X(X.begin(), X.rows(), X.cols());
    Eigen::Map<const Eigen::VectorXi> map_y(y.begin(), y.size());
    Eigen::Map<const Eigen::VectorXi> map_group_id(group_id.begin(), group_id.size());

    int n = map_X.rows();
    int p = map_X.cols();
    int na = K - 1;
    int total = na + p + 1;

    Eigen::VectorXd y_v(n);
    std::vector<int> gid_v(n);
    for (int i = 0; i < n; ++i) {
        y_v[i] = static_cast<double>(map_y[i]);
        gid_v[i] = map_group_id[i];
    }

    GLMMData dat(map_X, y_v, gid_v, n_gh, max_abs_log_sigma);
    FixedParamSpec fixed_spec = make_fixed_param_spec(total, fixed_idx, fixed_values);

    Eigen::VectorXd start_full(total);
    if (warm_start_params.isNotNull()) {
        Rcpp::NumericVector sv(warm_start_params);
        if (sv.size() == total) {
            for (int i = 0; i < total; ++i) start_full[i] = sv[i];
        } else {
            start_full.setZero();
            start_full[total - 1] = -3.0;
        }
    } else {
        start_full.setZero();
        start_full[total - 1] = -3.0;
    }
    start_full = apply_fixed_values(start_full, fixed_spec);

    Eigen::MatrixXd info_start;
    const Eigen::MatrixXd* info_start_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        info_start_ptr = &info_start;
    }

    if (link == "logit")   return run_fast_ordinal_clmm<LogitLink>(dat, K, j_T, estimate_only, maxit, eps_g, start_full, optimization_alg, fixed_spec, info_start_ptr);
    if (link == "probit")  return run_fast_ordinal_clmm<ProbitLink>(dat, K, j_T, estimate_only, maxit, eps_g, start_full, optimization_alg, fixed_spec, info_start_ptr);
    if (link == "cauchit") return run_fast_ordinal_clmm<CauchitLink>(dat, K, j_T, estimate_only, maxit, eps_g, start_full, optimization_alg, fixed_spec, info_start_ptr);
    if (link == "cloglog") return run_fast_ordinal_clmm<CloglogLink>(dat, K, j_T, estimate_only, maxit, eps_g, start_full, optimization_alg, fixed_spec, info_start_ptr);

    Rcpp::stop("Unknown link: " + link);
    return List::create();
}
