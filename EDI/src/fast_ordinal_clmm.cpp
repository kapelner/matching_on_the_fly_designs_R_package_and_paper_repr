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

    // Recover thresholds from log-difference parameterization
    inline std::vector<double> get_alpha(const Eigen::VectorXd& par) const {
        int na = K - 1;
        std::vector<double> alpha(na);
        alpha[0] = par[0];
        for (int k = 1; k < na; ++k) alpha[k] = alpha[k - 1] + std::exp(par[k]);
        return alpha;
    }

    double log_prob(double y_val, double eta, const Eigen::VectorXd& par) const {
        int y_ir = static_cast<int>(y_val);
        std::vector<double> alpha = get_alpha(par);
        double F_up = (y_ir >= K) ? 1.0 : Link::cdf(alpha[y_ir - 1] - eta);
        double F_lo = (y_ir <= 1) ? 0.0 : Link::cdf(alpha[y_ir - 2] - eta);
        return std::log(std::max(1e-15, F_up - F_lo));
    }

    double log_prob_derivs(double y_val, double eta, const Eigen::VectorXd& par, double& de, Eigen::VectorXd& dp) const {
        int y_ir = static_cast<int>(y_val);
        int na = K - 1;
        std::vector<double> alpha = get_alpha(par);
        
        double F_up = (y_ir >= K) ? 1.0 : Link::cdf(alpha[y_ir - 1] - eta);
        double F_lo = (y_ir <= 1) ? 0.0 : Link::cdf(alpha[y_ir - 2] - eta);
        double prob = std::max(1e-15, F_up - F_lo);
        double lp = std::log(prob);

        double f_up = (y_ir >= K) ? 0.0 : Link::pdf(alpha[y_ir - 1] - eta);
        double f_lo = (y_ir <= 1) ? 0.0 : Link::pdf(alpha[y_ir - 2] - eta);

        // Gradient w.r.t eta
        de = -(f_up - f_lo) / prob;

        // Gradient w.r.t alpha
        Eigen::VectorXd d_alpha = Eigen::VectorXd::Zero(na);
        if (y_ir < K) d_alpha[y_ir - 1] += f_up / prob;
        if (y_ir > 1) d_alpha[y_ir - 2] -= f_lo / prob;

        // Chain rule for log-diff parameterization
        dp.setZero(na);
        dp[0] = d_alpha.sum();
        for (int j = 1; j < na; ++j) {
            double s = 0.0;
            for (int k = j; k < na; ++k) s += d_alpha[k];
            dp[j] = s * std::exp(par[j]);
        }

        return lp;
    }
};

template <typename Link>
List run_fast_ordinal_clmm(const GLMMData& dat, int K, int j_T, bool estimate_only, int maxit, double eps_g, 
                           const Eigen::VectorXd& start_full, const std::string& optimization_alg,
                           const FixedParamSpec& fixed_spec) {
    OrdinalLikelihood<Link> model(K);
    GLMMObjective<OrdinalLikelihood<Link>> obj(dat, model);
    
    Eigen::VectorXd par = start_full;
    double neg_ll = NA_REAL;
    bool converged = false;
    
    try {
        LikelihoodFitResult fit = optimize_fixed_likelihood(obj, par, fixed_spec, maxit, eps_g, optimization_alg, "lbfgs");
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
    const Eigen::MatrixXd& X,
    const Eigen::VectorXi& y,
    const Eigen::VectorXi& group_id,
    int K,
    int j_T,
    std::string link = "logit",
    bool estimate_only = false,
    int n_gh = 20,
    double max_abs_log_sigma = 8.0,
    int maxit = 300,
    double eps_g = 1e-6,
    Rcpp::Nullable<Rcpp::NumericVector> start = R_NilValue,
    std::string optimization_alg = "lbfgs",
    Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue
) {
    int n = X.rows();
    int p = X.cols();
    int na = K - 1;
    int total = na + p + 1;

    Eigen::VectorXd y_v(n);
    std::vector<int> gid_v(n);
    for (int i = 0; i < n; ++i) {
        y_v[i] = static_cast<double>(y[i]);
        gid_v[i] = group_id[i];
    }

    GLMMData dat(X, y_v, gid_v, n_gh, max_abs_log_sigma);
    FixedParamSpec fixed_spec = make_fixed_param_spec(total, fixed_idx, fixed_values);

    Eigen::VectorXd start_full(total);
    if (start.isNotNull()) {
        Rcpp::NumericVector sv(start);
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

    if (link == "logit")   return run_fast_ordinal_clmm<LogitLink>(dat, K, j_T, estimate_only, maxit, eps_g, start_full, optimization_alg, fixed_spec);
    if (link == "probit")  return run_fast_ordinal_clmm<ProbitLink>(dat, K, j_T, estimate_only, maxit, eps_g, start_full, optimization_alg, fixed_spec);
    if (link == "cauchit") return run_fast_ordinal_clmm<CauchitLink>(dat, K, j_T, estimate_only, maxit, eps_g, start_full, optimization_alg, fixed_spec);
    if (link == "cloglog") return run_fast_ordinal_clmm<CloglogLink>(dat, K, j_T, estimate_only, maxit, eps_g, start_full, optimization_alg, fixed_spec);

    Rcpp::stop("Unknown link: " + link);
    return List::create();
}
