#ifndef EDI_ORDINAL_FIXED_LINK_HELPERS_H
#define EDI_ORDINAL_FIXED_LINK_HELPERS_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace edi_ordinal {

enum class Link {
    Logit,
    Probit,
    Cloglog,
    Cauchit
};

inline std::vector<double> init_levels(const Eigen::VectorXd& y) {
    std::vector<double> levels(y.data(), y.data() + y.size());
    std::sort(levels.begin(), levels.end());
    levels.erase(std::unique(levels.begin(), levels.end()), levels.end());
    return levels;
}

inline double cdf(Link link, double z) {
    switch (link) {
    case Link::Logit:
        if (z >= 0.0) {
            const double e = std::exp(-z);
            return 1.0 / (1.0 + e);
        } else {
            const double e = std::exp(z);
            return e / (1.0 + e);
        }
    case Link::Probit:
        return R::pnorm5(z, 0.0, 1.0, 1, 0);
    case Link::Cloglog:
        if (z > 5.0) return 1.0;
        if (z < -37.0) return 0.0;
        return 1.0 - std::exp(-std::exp(z));
    case Link::Cauchit:
        return 0.5 + std::atan(z) / M_PI;
    }
    return NA_REAL;
}

inline double pdf(Link link, double z) {
    switch (link) {
    case Link::Logit: {
        const double F = cdf(link, z);
        return F * (1.0 - F);
    }
    case Link::Probit:
        return R::dnorm4(z, 0.0, 1.0, 0);
    case Link::Cloglog:
        if (z > 5.0 || z < -37.0) return 0.0;
        return std::exp(z - std::exp(z));
    case Link::Cauchit:
        return 1.0 / (M_PI * (1.0 + z * z));
    }
    return NA_REAL;
}

inline double pdf_derivative(Link link, double z) {
    switch (link) {
    case Link::Logit: {
        const double F = cdf(link, z);
        const double f = F * (1.0 - F);
        return f * (1.0 - 2.0 * F);
    }
    case Link::Probit: {
        const double f = R::dnorm4(z, 0.0, 1.0, 0);
        return -z * f;
    }
    case Link::Cloglog: {
        if (z > 5.0 || z < -37.0) return 0.0;
        const double ez = std::exp(z);
        const double f = std::exp(z - ez);
        return f * (1.0 - ez);
    }
    case Link::Cauchit: {
        const double denom = 1.0 + z * z;
        return -2.0 * z / (M_PI * denom * denom);
    }
    }
    return NA_REAL;
}

inline int level_index(const std::vector<double>& levels, double y) {
    for (int k = 0; k < static_cast<int>(levels.size()); ++k) {
        if (y == levels[k]) return k;
    }
    return -1;
}

class FixedOrdinalRegression {
private:
    const Eigen::MatrixXd m_X;
    const Eigen::VectorXd m_y;
    const Eigen::VectorXd m_weights;
    const std::vector<double> m_levels;
    const int m_n;
    const int m_p;
    const int m_K;
    const Link m_link;
    const double m_eta_sign;
    const bool m_use_weights;
    mutable Eigen::VectorXd m_scratch_dq;
    mutable Eigen::MatrixXd m_scratch_d2q;
    mutable Eigen::VectorXd m_scratch_v;

    inline double obs_weight(int i) const {
        return m_use_weights ? std::max(m_weights[i], 0.0) : 1.0;
    }

    bool validate_params(const Eigen::VectorXd& params) const {
        const int n_alpha = m_K - 1;
        if (params.size() != n_alpha + m_p) return false;
        for (int k = 1; k < n_alpha; ++k) {
            if (params[k] <= params[k - 1]) return false;
        }
        return true;
    }

    void add_endpoint_derivatives(int alpha_idx,
                                  double z,
                                  double endpoint_sign,
                                  const Eigen::RowVectorXd& x,
                                  Eigen::VectorXd& dq,
                                  Eigen::MatrixXd& d2q) const {
        m_scratch_v.setZero();
        m_scratch_v[alpha_idx] = 1.0;
        m_scratch_v.tail(m_p) = m_eta_sign * x.transpose();

        const double f = pdf(m_link, z);
        const double fp = pdf_derivative(m_link, z);
        dq.noalias() += endpoint_sign * f * m_scratch_v;
        d2q.noalias() += endpoint_sign * fp * (m_scratch_v * m_scratch_v.transpose());
    }

    void add_endpoint_gradient(int alpha_idx,
                               double z,
                               double endpoint_sign,
                               const Eigen::RowVectorXd& x,
                               Eigen::VectorXd& dq) const {
        const int n_alpha = m_K - 1;
        const double f = pdf(m_link, z);
        dq[alpha_idx] += endpoint_sign * f;
        dq.tail(m_p).noalias() += endpoint_sign * f * m_eta_sign * x.transpose();
    }

public:
    FixedOrdinalRegression(const Eigen::MatrixXd& X,
                           const Eigen::VectorXd& y,
                           Link link,
                           double eta_sign,
                           const Eigen::VectorXd& weights = Eigen::VectorXd()) :
        m_X(X), m_y(y), m_weights(weights), m_levels(init_levels(y)), m_n(X.rows()), m_p(X.cols()),
        m_K(m_levels.size()), m_link(link), m_eta_sign(eta_sign), m_use_weights(weights.size() == X.rows()),
        m_scratch_dq((m_K - 1) + m_p),
        m_scratch_d2q((m_K - 1) + m_p, (m_K - 1) + m_p),
        m_scratch_v((m_K - 1) + m_p) {}

    static std::vector<double> init_levels_static(const Eigen::VectorXd& y) {
        return init_levels(y);
    }

    double neg_log_likelihood(const Eigen::VectorXd& params) const {
        if (!validate_params(params)) return 1e10;
        const int n_alpha = m_K - 1;
        const Eigen::VectorXd alpha = params.head(n_alpha);
        const Eigen::VectorXd beta = params.tail(m_p);
        const Eigen::VectorXd eta = m_X * beta;
        double nll = 0.0;

        for (int i = 0; i < m_n; ++i) {
            const int yi_idx = level_index(m_levels, m_y[i]);
            if (yi_idx < 0) return 1e10;
            const double p_upper = (yi_idx == m_K - 1) ? 1.0 : cdf(m_link, alpha[yi_idx] + m_eta_sign * eta[i]);
            const double p_lower = (yi_idx == 0) ? 0.0 : cdf(m_link, alpha[yi_idx - 1] + m_eta_sign * eta[i]);
            const double prob = std::max(1e-12, p_upper - p_lower);
            nll -= obs_weight(i) * std::log(prob);
        }
        return nll;
    }

    double operator()(const Eigen::VectorXd& params, Eigen::VectorXd& grad) const {
        const int n_params = params.size();
        grad = Eigen::VectorXd::Zero(n_params);
        if (!validate_params(params)) return 1e10;

        const int n_alpha = m_K - 1;
        const Eigen::VectorXd alpha = params.head(n_alpha);
        const Eigen::VectorXd beta = params.tail(m_p);
        const Eigen::VectorXd eta = m_X * beta;
        double nll = 0.0;

        for (int i = 0; i < m_n; ++i) {
            const int yi_idx = level_index(m_levels, m_y[i]);
            if (yi_idx < 0) return 1e10;

            const double p_upper = (yi_idx == m_K - 1) ? 1.0 : cdf(m_link, alpha[yi_idx] + m_eta_sign * eta[i]);
            const double p_lower = (yi_idx == 0) ? 0.0 : cdf(m_link, alpha[yi_idx - 1] + m_eta_sign * eta[i]);
            const double prob = std::max(1e-12, p_upper - p_lower);
            m_scratch_dq.setZero();
            const double wi = obs_weight(i);

            if (yi_idx < m_K - 1) {
                add_endpoint_gradient(yi_idx, alpha[yi_idx] + m_eta_sign * eta[i], 1.0, m_X.row(i), m_scratch_dq);
            }
            if (yi_idx > 0) {
                add_endpoint_gradient(yi_idx - 1, alpha[yi_idx - 1] + m_eta_sign * eta[i], -1.0, m_X.row(i), m_scratch_dq);
            }
            nll -= wi * std::log(prob);
            grad.noalias() -= wi * m_scratch_dq / prob;
        }
        return nll;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& params) const {
        const int n_params = params.size();
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n_params, n_params);
        if (!validate_params(params)) {
            H.diagonal().array() = 1e10;
            return H;
        }

        const int n_alpha = m_K - 1;
        const Eigen::VectorXd alpha = params.head(n_alpha);
        const Eigen::VectorXd beta = params.tail(m_p);
        const Eigen::VectorXd eta = m_X * beta;

        for (int i = 0; i < m_n; ++i) {
            const int yi_idx = level_index(m_levels, m_y[i]);
            if (yi_idx < 0) continue;

            const double p_upper = (yi_idx == m_K - 1) ? 1.0 : cdf(m_link, alpha[yi_idx] + m_eta_sign * eta[i]);
            const double p_lower = (yi_idx == 0) ? 0.0 : cdf(m_link, alpha[yi_idx - 1] + m_eta_sign * eta[i]);
            const double prob = std::max(1e-12, p_upper - p_lower);
            m_scratch_dq.setZero();
            m_scratch_d2q.setZero();
            const double wi = obs_weight(i);

            if (yi_idx < m_K - 1) {
                add_endpoint_derivatives(yi_idx, alpha[yi_idx] + m_eta_sign * eta[i], 1.0, m_X.row(i), m_scratch_dq, m_scratch_d2q);
            }
            if (yi_idx > 0) {
                add_endpoint_derivatives(yi_idx - 1, alpha[yi_idx - 1] + m_eta_sign * eta[i], -1.0, m_X.row(i), m_scratch_dq, m_scratch_d2q);
            }
            H.noalias() += wi * ((m_scratch_dq * m_scratch_dq.transpose()) / (prob * prob) - m_scratch_d2q / prob);
        }
        return 0.5 * (H + H.transpose());
    }
};

} // namespace edi_ordinal

#endif
