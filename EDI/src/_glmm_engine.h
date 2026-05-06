#ifndef EDI_GLMM_ENGINE_H
#define EDI_GLMM_ENGINE_H

#include <RcppEigen.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <limits>
#include "_helper_functions.h"

namespace glmm {

struct GHRule {
    Eigen::VectorXd nodes;
    Eigen::VectorXd log_norm_weights;
};

inline GHRule gauss_hermite_rule(int n) {
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n, n);
    for (int i = 0; i < n - 1; ++i) {
        const double v = std::sqrt((i + 1.0) / 2.0);
        J(i, i + 1) = v;
        J(i + 1, i) = v;
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(J);
    GHRule rule;
    rule.nodes = es.eigenvalues();
    rule.log_norm_weights = (std::sqrt(M_PI) * es.eigenvectors().row(0).array().square()).log()
                            - 0.5 * std::log(M_PI);
    return rule;
}

inline double log_sum_exp(const Eigen::VectorXd& x) {
    const double m = x.maxCoeff();
    if (!std::isfinite(m)) return m;
    return m + std::log((x.array() - m).exp().sum());
}

struct GLMMData {
    Eigen::MatrixXd X_s;
    Eigen::VectorXd y_s;
    std::vector<int> grp_start, grp_size;
    int n, p, G;
    GHRule gh;
    double max_abs_log_sigma;

    GLMMData(const Eigen::MatrixXd& X, 
             const Eigen::VectorXd& y, 
             const std::vector<int>& gid,
             int n_gh, double max_abs_log_sigma_)
        : n(X.rows()), p(X.cols()), gh(gauss_hermite_rule(n_gh)), max_abs_log_sigma(max_abs_log_sigma_) {
        
        std::vector<int> ord(n);
        std::iota(ord.begin(), ord.end(), 0);
        std::stable_sort(ord.begin(), ord.end(), [&](int a, int b){ return gid[a] < gid[b]; });

        X_s.resize(n, p);
        y_s.resize(n);
        for(int i = 0; i < n; ++i) {
            X_s.row(i) = X.row(ord[i]);
            y_s[i] = y[ord[i]];
        }

        int prev = -1;
        for(int i = 0; i < n; ++i) {
            if(gid[ord[i]] != prev) {
                grp_start.push_back(i);
                grp_size.push_back(1);
                prev = gid[ord[i]];
            } else {
                grp_size.back()++;
            }
        }
        G = static_cast<int>(grp_start.size());
    }
};

inline double sigma_penalty(double log_sigma, double center = 5.0, double scale = 10.0) {
    const double d = std::abs(log_sigma) - center;
    if (d <= 0.0) return 0.0;
    return scale * d * d;
}

inline double sigma_penalty_grad(double log_sigma, double center = 5.0, double scale = 10.0) {
    const double d = std::abs(log_sigma) - center;
    if (d <= 0.0) return 0.0;
    return 2.0 * scale * d * (log_sigma > 0 ? 1.0 : -1.0);
}

inline double sigma_penalty_hessian(double log_sigma, double center = 5.0, double scale = 10.0) {
    const double d = std::abs(log_sigma) - center;
    if (d <= 0.0) return 0.0;
    return 2.0 * scale;
}

// A generic GLMM Objective that uses Gauss-Hermite quadrature.
// Model template must provide:
// - int n_model_params()
// - void unpack(const VectorXd& par, ...)
// - double log_prob(double y, double eta, ...)
// - void log_prob_derivs(double y, double eta, double& d_eta, VectorXd& d_par, ...)
// - void log_prob_hessians(double y, double eta, double& d2_eta, MatrixXd& d2_par, MatrixXd& d2_par_eta, ...)
template <typename Model>
class GLMMObjective {
    const GLMMData& dat;
    Model model;

public:
    GLMMObjective(const GLMMData& d, const Model& m) : dat(d), model(m) {}

    double value(const Eigen::VectorXd& par) const {
        const int nm = model.n_model_params();
        const double log_sigma = par[nm + dat.p];
        if (!std::isfinite(log_sigma) || std::abs(log_sigma) > dat.max_abs_log_sigma) return 1e100;
        const double sigma = std::exp(log_sigma);
        const double pen = sigma_penalty(log_sigma);

        const Eigen::VectorXd beta = par.segment(nm, dat.p);
        const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * dat.gh.nodes;
        const int nn = static_cast<int>(b_vals.size());

        double total_ll = 0.0;
        for (int gi = 0; gi < dat.G; ++gi) {
            const int start = dat.grp_start[gi];
            const int sz = dat.grp_size[gi];
            const Eigen::VectorXd eta0 = dat.X_s.middleRows(start, sz) * beta;

            Eigen::VectorXd log_terms(nn);
            for (int k = 0; k < nn; ++k) {
                double ll = dat.gh.log_norm_weights[k];
                for (int r = 0; r < sz; ++r) {
                    ll += model.log_prob(dat.y_s[start + r], eta0[r] + b_vals[k], par);
                }
                log_terms[k] = ll;
            }
            const double ll_g = log_sum_exp(log_terms);
            if (!std::isfinite(ll_g)) return 1e100;
            total_ll += ll_g;
        }
        return -total_ll + pen;
    }

    double operator()(const Eigen::VectorXd& par, Eigen::VectorXd& grad) {
        const int nm = model.n_model_params();
        const double log_sigma = par[nm + dat.p];
        const double sigma = std::exp(log_sigma);
        const double pen = sigma_penalty(log_sigma);
        
        const Eigen::VectorXd beta = par.segment(nm, dat.p);
        const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * dat.gh.nodes;
        const int nn = static_cast<int>(b_vals.size());

        double total_nll = pen;
        grad.setZero();
        grad[nm + dat.p] += sigma_penalty_grad(log_sigma);

        for (int gi = 0; gi < dat.G; ++gi) {
            const int start = dat.grp_start[gi];
            const int sz = dat.grp_size[gi];
            const Eigen::MatrixXd Xg = dat.X_s.middleRows(start, sz);
            const Eigen::VectorXd eta0 = Xg * beta;

            Eigen::VectorXd log_terms(nn);
            std::vector<Eigen::VectorXd> dL_dp_sum_nodes(nn, Eigen::VectorXd::Zero(nm));
            std::vector<Eigen::VectorXd> dL_de_nodes(nn, Eigen::VectorXd::Zero(sz));

            for (int k = 0; k < nn; ++k) {
                double ll = dat.gh.log_norm_weights[k];
                for (int r = 0; r < sz; ++r) {
                    double de;
                    Eigen::VectorXd dp(nm);
                    double lp = model.log_prob_derivs(dat.y_s[start + r], eta0[r] + b_vals[k], par, de, dp);
                    ll += lp;
                    dL_dp_sum_nodes[k] += dp;
                    dL_de_nodes[k][r] = de;
                }
                log_terms[k] = ll;
            }

            const double ll_g = log_sum_exp(log_terms);
            total_nll -= ll_g;

            for (int k = 0; k < nn; ++k) {
                double pk = std::exp(log_terms[k] - ll_g);
                if (pk < 1e-15) continue;

                const Eigen::VectorXd dLi_db = Xg.transpose() * dL_de_nodes[k];
                grad.head(nm) -= pk * dL_dp_sum_nodes[k];
                grad.segment(nm, dat.p) -= pk * dLi_db;
                grad[nm + dat.p] -= pk * dL_de_nodes[k].sum() * std::sqrt(2.0) * dat.gh.nodes[k] * sigma;
            }
        }
        return total_nll;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& par) {
        return numerical_hessian(*this, par); // Fallback to numerical Hessian for simplicity in generic engine
    }
};

} // namespace glmm

#endif
