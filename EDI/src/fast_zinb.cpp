// Zero-Inflated Negative Binomial (ZINB) regression.
//
// Model:
//   P(Y=0)   = pi + (1-pi) * (theta/(theta+mu))^theta
//   P(Y=y>0) = (1-pi) * NegBin(y; mu, theta)
//
// Parameter vector: [beta_cond(p_cond), beta_zi(p_zi), log_theta]
//   mu  = exp(eta_cond),  pi = sigmoid(eta_zi),  theta = exp(log_theta)
//
// Analytic gradient; numerical Hessian via numerical_hessian().

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
        return numerical_hessian(*this, par);
    }
};

} // namespace

// [[Rcpp::export]]
List fast_zinb_cpp(
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& Xcond,
    const Eigen::MatrixXd& Xzi,
    Nullable<NumericVector> start_params = R_NilValue,
    bool estimate_only = false,
    int maxit = 1000,
    double tol = 1e-6,
    std::string optimization_alg = "newton_raphson"
) {
    const int pc = Xcond.cols();
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

    ZeroInflatedNegBin fun(y, Xcond, Xzi);
    LikelihoodFitResult fit;
    try {
        fit = optimize_likelihood(fun, par, maxit, tol, optimization_alg, "newton_raphson");
    } catch (...) {
        return List::create(Named("converged") = false);
    }
    par = fit.params;

    if (estimate_only) {
        return List::create(
            Named("coefficients") = List::create(
                Named("cond") = par.head(pc),
                Named("zi")   = par.segment(pc, pz)
            ),
            Named("converged") = fit.converged,
            Named("neg_ll")    = fit.value
        );
    }

    Eigen::MatrixXd H = fun.hessian(par);
    Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
    Eigen::MatrixXd vcov = Eigen::MatrixXd::Constant(total, total, NA_REAL);
    if (ldlt.info() == Eigen::Success) {
        Eigen::MatrixXd inv = ldlt.solve(Eigen::MatrixXd::Identity(total, total));
        if (inv.allFinite()) vcov = inv;
    }

    return List::create(
        Named("coefficients") = List::create(
            Named("cond") = par.head(pc),
            Named("zi")   = par.segment(pc, pz)
        ),
        Named("vcov")      = vcov,
        Named("converged") = fit.converged,
        Named("neg_ll")    = fit.value
    );
}
