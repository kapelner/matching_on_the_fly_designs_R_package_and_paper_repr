// Zero-Inflated Negative Binomial (ZINB) regression.
//
// Model:
//   P(Y=0)   = pi + (1-pi) * (theta/(theta+mu))^theta
//   P(Y=y>0) = (1-pi) * NegBin(y; mu, theta)
//
// Parameter vector: [beta_cond(p_cond), beta_zi(p_zi), log_theta]
//   mu  = exp(eta_cond),  pi = sigmoid(eta_zi),  theta = exp(log_theta)
//
// Analytic gradient; analytic Hessian.

#include "_helper_functions.h"
#include <RcppEigen.h>
#include <Rmath.h>
#include <cmath>

using namespace Rcpp;

namespace {

// log(1 + exp(x)) — numerically stable
inline double lse_zinb(double x) {
    if (x > 0.0) return x + std::log1p(std::exp(-x));
    return std::log1p(std::exp(x));
}

// sigmoid(x)
inline double sigmoid_zinb(double x) {
    if (x >  35.0) return 1.0;
    if (x < -35.0) return 0.0;
    return 1.0 / (1.0 + std::exp(-x));
}

class ZeroInflatedNegBin {
    const Eigen::VectorXd m_y;
    const Eigen::MatrixXd m_Xc;
    const Eigen::MatrixXd m_Xz;
    const int m_n, m_pc, m_pz;

public:
    ZeroInflatedNegBin(const Eigen::VectorXd& y,
                       const Eigen::MatrixXd& Xc,
                       const Eigen::MatrixXd& Xz)
        : m_y(y), m_Xc(Xc), m_Xz(Xz),
          m_n(y.size()), m_pc(Xc.cols()), m_pz(Xz.cols()) {}

    double operator()(const Eigen::VectorXd& par, Eigen::VectorXd& grad) {
        const Eigen::VectorXd bc  = par.head(m_pc);
        const Eigen::VectorXd bz  = par.segment(m_pc, m_pz);
        const double log_theta    = par[m_pc + m_pz];
        const double theta        = std::exp(std::min(log_theta, 700.0));

        const Eigen::VectorXd eta_c = m_Xc * bc;
        const Eigen::VectorXd eta_z = m_Xz * bz;

        grad.setZero(m_pc + m_pz + 1);

        double nll = 0.0;
        double d_log_theta = 0.0;

        for (int i = 0; i < m_n; ++i) {
            const double mu  = std::exp(std::min(eta_c[i], 700.0));
            const double pi  = sigmoid_zinb(eta_z[i]);
            const double A   = theta + mu;
            const double p0  = std::pow(theta / A, theta);  // NB prob at y=0

            if (m_y[i] == 0.0) {
                // q = pi + (1-pi)*p0
                const double q = std::max(pi + (1.0 - pi) * p0, 1e-300);
                nll -= std::log(q);

                // d/d eta_c
                const double dc = (1.0 - pi) * theta * p0 * mu / (A * q);
                grad.head(m_pc).noalias() += dc * m_Xc.row(i).transpose();

                // d/d eta_z
                const double dz = -(1.0 - p0) * pi * (1.0 - pi) / q;
                grad.segment(m_pc, m_pz).noalias() += dz * m_Xz.row(i).transpose();

                // d/d log_theta = theta * d/d theta
                // d p0/d theta = p0 * (log(theta/A) + mu/A)
                const double dp0_dtheta = p0 * (std::log(theta / A) + mu / A);
                const double dq_dtheta  = (1.0 - pi) * dp0_dtheta;
                d_log_theta += -theta * dq_dtheta / q;

            } else {
                // NLL_i = -log(1-pi) - lgamma(y+theta) + lgamma(theta) + lgamma(y+1)
                //         - theta*log(theta/A) - y*log(mu/A)
                const double yi = m_y[i];
                nll += lse_zinb(eta_z[i])   // = -log(1-pi)
                     - R::lgammafn(yi + theta)
                     + R::lgammafn(theta)
                     + R::lgammafn(yi + 1.0)
                     - theta * std::log(theta / A)
                     - yi * std::log(mu / A);

                // d/d eta_c
                const double dc = mu * (yi + theta) / A - yi;
                grad.head(m_pc).noalias() += dc * m_Xc.row(i).transpose();

                // d/d eta_z
                grad.segment(m_pc, m_pz).noalias() += pi * m_Xz.row(i).transpose();

                // d/d log_theta = theta * (digamma(theta) - digamma(yi+theta) - log(theta/A) + (yi-mu)/A)
                const double dt = theta * (R::digamma(theta) - R::digamma(yi + theta)
                                           - std::log(theta / A) + (yi - mu) / A);
                d_log_theta += dt;
            }
        }
        grad[m_pc + m_pz] = d_log_theta;
        return nll;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& par) {
        const int total = m_pc + m_pz + 1;
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total, total);

        const Eigen::VectorXd bc  = par.head(m_pc);
        const Eigen::VectorXd bz  = par.segment(m_pc, m_pz);
        const double log_theta    = par[m_pc + m_pz];
        const double theta        = std::exp(std::min(log_theta, 700.0));

        const Eigen::VectorXd eta_c = m_Xc * bc;
        const Eigen::VectorXd eta_z = m_Xz * bz;

        for (int i = 0; i < m_n; ++i) {
            const double mu = std::exp(std::min(eta_c[i], 700.0));
            const double pi = sigmoid_zinb(eta_z[i]);
            const double A  = theta + mu;
            const double p0 = std::pow(theta / A, theta);

            const Eigen::VectorXd xc = m_Xc.row(i);
            const Eigen::VectorXd xz = m_Xz.row(i);

            if (m_y[i] == 0.0) {
                const double q = std::max(pi + (1.0 - pi) * p0, 1e-300);
                
                const double dq_dmu = (1.0 - pi) * p0 * (-theta / A);
                const double dq_detac = dq_dmu * mu;
                
                const double dpi_detaz = pi * (1.0 - pi);
                const double dq_detaz = (1.0 - p0) * dpi_detaz;
                
                const double dp0_dtheta = p0 * (std::log(theta / A) + mu / A);
                const double dq_dtheta = (1.0 - pi) * dp0_dtheta;
                const double dq_dlogtheta = dq_dtheta * theta;

                const double d2q_dmu2 = (1.0 - pi) * p0 * (theta / (A * A)) * (theta + 1.0);
                const double d2q_detac2 = d2q_dmu2 * mu * mu + dq_dmu * mu;
                
                const double d2pi_detaz2 = pi * (1.0 - pi) * (1.0 - 2.0 * pi);
                const double d2q_detaz2 = (1.0 - p0) * d2pi_detaz2;
                
                const double d2p0_dtheta2 = p0 * (std::pow(std::log(theta / A) + mu / A, 2) + 
                                               (1.0 / theta - 2.0 / A + theta / (A * A)));
                const double d2q_dtheta2 = (1.0 - pi) * d2p0_dtheta2;
                const double d2q_dlogtheta2 = d2q_dtheta2 * theta * theta + dq_dtheta * theta;

                const double d2q_detac_detaz = -mu * (theta / A) * p0 * dpi_detaz;
                const double d2q_detac_dlogtheta = (1.0 - pi) * theta * (
                    dp0_dtheta * (-theta / A) + p0 * (-mu / (A * A))
                );
                const double d2q_detaz_dlogtheta = -dpi_detaz * dp0_dtheta * theta;

                const double inv_q = 1.0 / q;
                const double inv_q2 = inv_q * inv_q;

                const double h_cc = (inv_q2 * dq_detac * dq_detac - inv_q * d2q_detac2);
                const double h_zz = (inv_q2 * dq_detaz * dq_detaz - inv_q * d2q_detaz2);
                const double h_tt = (inv_q2 * dq_dlogtheta * dq_dlogtheta - inv_q * d2q_dlogtheta2);
                const double h_cz = (inv_q2 * dq_detac * dq_detaz - inv_q * d2q_detac_detaz);
                const double h_ct = (inv_q2 * dq_detac * dq_dlogtheta - inv_q * d2q_detac_dlogtheta);
                const double h_zt = (inv_q2 * dq_detaz * dq_dlogtheta - inv_q * d2q_detaz_dlogtheta);

                for (int r = 0; r < m_pc; ++r) {
                    for (int c = 0; r >= c && c < m_pc; ++c) H(r, c) += h_cc * xc[r] * xc[c];
                    for (int c = 0; c < m_pz; ++c) H(r, m_pc + c) += h_cz * xc[r] * xz[c];
                    H(r, total - 1) += h_ct * xc[r];
                }
                for (int r = 0; r < m_pz; ++r) {
                    for (int c = 0; r >= c && c < m_pz; ++c) H(m_pc + r, m_pc + c) += h_zz * xz[r] * xz[c];
                    H(m_pc + r, total - 1) += h_zt * xz[r];
                }
                H(total - 1, total - 1) += h_tt;

            } else {
                const double yi = m_y[i];
                const double h_zz = pi * (1.0 - pi);
                const double h_cc = theta * mu * (yi + theta) / (A * A);
                const double dnb_dtheta = R::digamma(theta) - R::digamma(yi + theta) - std::log(theta / A) + (yi - mu) / A;
                const double d2nb_dtheta2 = R::trigamma(theta) - R::trigamma(yi + theta) - 1.0 / theta + 2.0 / A - (yi + theta) / (A * A);
                const double h_tt = (d2nb_dtheta2 * theta * theta + dnb_dtheta * theta);
                const double h_ct = (-mu * (yi - mu) / (A * A)) * theta;

                for (int r = 0; r < m_pc; ++r) {
                    for (int c = 0; r >= c && c < m_pc; ++c) H(r, c) += h_cc * xc[r] * xc[c];
                    H(r, total - 1) += h_ct * xc[r];
                }
                for (int r = 0; r < m_pz; ++r) {
                    for (int c = 0; r >= c && c < m_pz; ++c) H(m_pc + r, m_pc + c) += h_zz * xz[r] * xz[c];
                }
                H(total - 1, total - 1) += h_tt;
            }
        }
        for (int r = 0; r < total; ++r) {
            for (int c = 0; c < r; ++c) H(c, r) = H(r, c);
        }
        return H;
    }
};

} // namespace

// [[Rcpp::export]]
Eigen::VectorXd get_zinb_score_cpp(const Eigen::MatrixXd& X,
                                   const Eigen::VectorXd& y,
                                   const Eigen::MatrixXd& Xzi,
                                   const Eigen::VectorXd& params) {
    ZeroInflatedNegBin fun(y, X, Xzi);
    return likelihood_score(fun, params);
}

// [[Rcpp::export]]
Eigen::MatrixXd get_zinb_hessian_cpp(const Eigen::MatrixXd& X,
                                     const Eigen::VectorXd& y,
                                     const Eigen::MatrixXd& Xzi,
                                     const Eigen::VectorXd& params) {
    ZeroInflatedNegBin fun(y, X, Xzi);
    return -fun.hessian(params);
}

// [[Rcpp::export]]
double get_zinb_neg_loglik_cpp(const Eigen::MatrixXd& X,
                               const Eigen::VectorXd& y,
                               const Eigen::MatrixXd& Xzi,
                               const Eigen::VectorXd& params) {
    ZeroInflatedNegBin fun(y, X, Xzi);
    return likelihood_value(fun, params);
}

// [[Rcpp::export]]
List fast_zinb_cpp(
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& Xzi,
    Nullable<NumericVector> start_params = R_NilValue,
    bool estimate_only = false,
    int maxit = 1000,
    double tol = 1e-6,
    std::string optimization_alg = "newton_raphson",
    Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue
) {
    const int pc = X.cols();
    const int pz = Xzi.cols();
    const int total = pc + pz + 1;

    Eigen::VectorXd par(total);
    if (start_params.isNotNull()) {
        par = as<Eigen::VectorXd>(NumericVector(start_params));
    } else {
        par.setZero();
        // Init cond intercept from mean of positive observations; zi starts at 0;
        // log_theta = 0 (theta=1). Starting zi at logit(prop_zeros) is wrong
        // because it attributes all zeros to the structural component, biasing the
        // optimizer toward a poor local minimum.
        double sum_pos = 0.0; int cnt_pos = 0;
        for (int i = 0; i < y.size(); ++i)
            if (y[i] > 0.0) { sum_pos += y[i]; ++cnt_pos; }
        if (cnt_pos > 0) par[0] = std::log(sum_pos / cnt_pos);
        // par[pc..pc+pz-1] = 0 (zi), par[pc+pz] = 0 (log_theta) already set
    }

    ZeroInflatedNegBin fun(y, X, Xzi);
    FixedParamSpec fixed_spec = make_fixed_param_spec(total, fixed_idx, fixed_values);
    
    Eigen::MatrixXd H_start;
    const Eigen::MatrixXd* h_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        H_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        h_ptr = &H_start;
    }
    
    LikelihoodFitResult fit;
    try {
        fit = optimize_fixed_likelihood(fun, par, fixed_spec, maxit, tol, optimization_alg, "newton_raphson", 0, h_ptr);
    } catch (...) {
        return List::create(Named("converged") = false);
    }
    par = fit.params;
    Eigen::VectorXd score = get_zinb_score_cpp(X, y, Xzi, par);
    Eigen::MatrixXd information = fun.hessian(par);

    if (estimate_only) {
        return List::create(
            Named("params") = par,
            Named("coefficients") = List::create(
                Named("cond") = par.head(pc),
                Named("zi")   = par.segment(pc, pz)
            ),
            Named("converged") = fit.converged,
            Named("neg_ll")    = fit.value,
            Named("neg_loglik") = fit.value,
            Named("loglik") = R_finite(fit.value) ? -fit.value : NA_REAL,
            Named("fisher_information") = information
        );
    }

    Eigen::MatrixXd vcov = expand_free_covariance(total, fixed_spec, covariance_from_information(subset_matrix(information, fixed_spec.free_idx, fixed_spec.free_idx)), true);

    return List::create(
        Named("params") = par,
        Named("coefficients") = List::create(
            Named("cond") = par.head(pc),
            Named("zi")   = par.segment(pc, pz)
        ),
        Named("vcov")      = vcov,
        Named("score")     = score,
        Named("observed_information") = information,
        Named("information") = information,
        Named("information_type") = "observed",
        Named("hessian")   = -information,
        Named("converged") = fit.converged,
        Named("neg_ll")    = fit.value,
        Named("neg_loglik") = fit.value,
        Named("loglik") = R_finite(fit.value) ? -fit.value : NA_REAL,
        Named("fisher_information") = information
    );
}
