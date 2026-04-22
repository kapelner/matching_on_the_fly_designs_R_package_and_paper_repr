#ifndef EDI_HELPERS_H
#define EDI_HELPERS_H

#include <RcppEigen.h>
#include <optimization/LBFGS.h>
#include <vector>
#include <set>
#include <limits>
#include <cmath>
#include <string>

// Pure C++ result structure to avoid R List contention
struct ModelResult {
    Eigen::VectorXd b;
    Eigen::VectorXd mu;
    Eigen::MatrixXd XtWX;
    double ssq_b_j;
    double ssq_b_2;
    double dispersion;
    double sigma2_hat;
    bool converged;

    ModelResult() : ssq_b_j(NA_REAL), ssq_b_2(NA_REAL), dispersion(NA_REAL), sigma2_hat(NA_REAL), converged(false) {}
};

// Pure C++ internal helpers
double compute_diagonal_inverse_entry(const Eigen::MatrixXd& M, int j);

// R-facing exports
double eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(Eigen::MatrixXd M, int j);
Eigen::MatrixXd eigen_Xt_times_diag_w_times_X_cpp(Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::VectorXd> w);

struct FixedParamSpec {
    Eigen::VectorXi fixed_idx;
    Eigen::VectorXi free_idx;
    Eigen::VectorXd fixed_values;
    bool has_fixed;

    FixedParamSpec() : has_fixed(false) {}
};

inline FixedParamSpec make_fixed_param_spec(
    int n_params,
    Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue
) {
    FixedParamSpec spec;

    if (fixed_idx.isNull()) {
        spec.free_idx.resize(n_params);
        for (int i = 0; i < n_params; ++i) spec.free_idx[i] = i;
        return spec;
    }

    Rcpp::IntegerVector fixed_idx_r(fixed_idx);
    if (fixed_idx_r.size() == 0) {
        spec.free_idx.resize(n_params);
        for (int i = 0; i < n_params; ++i) spec.free_idx[i] = i;
        return spec;
    }
    if (fixed_values.isNull()) {
        Rcpp::stop("fixed_values must be supplied when fixed_idx is non-empty");
    }
    Rcpp::NumericVector fixed_values_r(fixed_values);
    if (fixed_values_r.size() != fixed_idx_r.size()) {
        Rcpp::stop("fixed_idx and fixed_values must have the same length");
    }

    std::vector<int> fixed_zero_based(fixed_idx_r.size());
    std::set<int> seen;
    for (int i = 0; i < fixed_idx_r.size(); ++i) {
        int idx = fixed_idx_r[i] - 1;
        if (idx < 0 || idx >= n_params) {
            Rcpp::stop("fixed_idx entries must be one-based indices within the parameter vector");
        }
        if (!seen.insert(idx).second) {
            Rcpp::stop("fixed_idx cannot contain duplicate entries");
        }
        if (!R_finite(fixed_values_r[i])) {
            Rcpp::stop("fixed_values must be finite");
        }
        fixed_zero_based[i] = idx;
    }

    spec.has_fixed = true;
    spec.fixed_idx.resize(fixed_zero_based.size());
    spec.fixed_values.resize(fixed_zero_based.size());
    std::vector<int> is_fixed(n_params, 0);
    for (int i = 0; i < static_cast<int>(fixed_zero_based.size()); ++i) {
        spec.fixed_idx[i] = fixed_zero_based[i];
        spec.fixed_values[i] = fixed_values_r[i];
        is_fixed[fixed_zero_based[i]] = 1;
    }

    int n_free = n_params - fixed_zero_based.size();
    if (n_free <= 0) {
        Rcpp::stop("at least one parameter must remain free");
    }
    spec.free_idx.resize(n_free);
    int k = 0;
    for (int i = 0; i < n_params; ++i) {
        if (!is_fixed[i]) spec.free_idx[k++] = i;
    }

    return spec;
}

inline Eigen::VectorXd subset_vector(const Eigen::VectorXd& x, const Eigen::VectorXi& idx) {
    Eigen::VectorXd out(idx.size());
    for (int i = 0; i < idx.size(); ++i) out[i] = x[idx[i]];
    return out;
}

inline Eigen::MatrixXd subset_matrix(const Eigen::MatrixXd& M, const Eigen::VectorXi& row_idx, const Eigen::VectorXi& col_idx) {
    Eigen::MatrixXd out(row_idx.size(), col_idx.size());
    for (int i = 0; i < row_idx.size(); ++i) {
        for (int j = 0; j < col_idx.size(); ++j) {
            out(i, j) = M(row_idx[i], col_idx[j]);
        }
    }
    return out;
}

inline Eigen::VectorXd apply_fixed_values(Eigen::VectorXd params, const FixedParamSpec& spec) {
    for (int i = 0; i < spec.fixed_idx.size(); ++i) {
        params[spec.fixed_idx[i]] = spec.fixed_values[i];
    }
    return params;
}

inline Eigen::VectorXd expand_free_params(const Eigen::VectorXd& free_params,
                                          const Eigen::VectorXd& full_template,
                                          const FixedParamSpec& spec) {
    Eigen::VectorXd full = apply_fixed_values(full_template, spec);
    for (int i = 0; i < spec.free_idx.size(); ++i) {
        full[spec.free_idx[i]] = free_params[i];
    }
    return full;
}

inline Eigen::MatrixXd expand_free_covariance(int n_params,
                                              const FixedParamSpec& spec,
                                              const Eigen::MatrixXd& cov_free,
                                              bool fixed_as_na = true) {
    Eigen::MatrixXd cov(n_params, n_params);
    if (fixed_as_na) {
        cov.setConstant(NA_REAL);
    } else {
        cov.setZero();
    }
    for (int i = 0; i < spec.free_idx.size(); ++i) {
        for (int j = 0; j < spec.free_idx.size(); ++j) {
            cov(spec.free_idx[i], spec.free_idx[j]) = cov_free(i, j);
        }
    }
    return cov;
}

inline double plogis_safe(double x) {
    if (x >= 0.0) { const double z = std::exp(-x); return 1.0 / (1.0 + z); }
    const double z = std::exp(x); return z / (1.0 + z);
}

template <typename Functor>
inline Eigen::VectorXd numerical_gradient(const Functor& fun, const Eigen::VectorXd& par) {
    const int n = par.size();
    Eigen::VectorXd grad(n);
    for (int i = 0; i < n; ++i) {
        double h = 1e-5 * std::max(1.0, std::abs(par[i]));
        Eigen::VectorXd p_plus = par, p_minus = par;
        p_plus[i] += h; p_minus[i] -= h;
        grad[i] = (fun.value(p_plus) - fun.value(p_minus)) / (2.0 * h);
    }
    return grad;
}

template <typename Functor>
inline Eigen::MatrixXd numerical_hessian(Functor& fun, const Eigen::VectorXd& par) {
    const int n = par.size();
    Eigen::MatrixXd hess(n, n);
    for (int i = 0; i < n; ++i) {
        double h = 1e-4 * std::max(1.0, std::abs(par[i]));
        Eigen::VectorXd p_plus = par, p_minus = par;
        p_plus[i] += h; p_minus[i] -= h;
        
        Eigen::VectorXd g_plus(n), g_minus(n);
        fun(p_plus, g_plus);
        fun(p_minus, g_minus);
        
        hess.col(i) = (g_plus - g_minus) / (2.0 * h);
    }
    return (hess + hess.transpose()) / 2.0;
}

struct LikelihoodFitResult {
    Eigen::VectorXd params;
    double value;
    int niter;
    bool converged;

    LikelihoodFitResult() :
        value(std::numeric_limits<double>::quiet_NaN()),
        niter(0),
        converged(false) {}
};

inline std::string normalize_optimizer_algorithm(const std::string& optimization_alg,
                                                 const std::string& default_optimization_alg,
                                                 bool allow_irls) {
    std::string alg = optimization_alg.empty() ? default_optimization_alg : optimization_alg;
    if (alg == "nr" || alg == "newton" || alg == "newton-raphson") {
        alg = "newton_raphson";
    } else if (alg == "l-bfgs" || alg == "L-BFGS" || alg == "LBFGS") {
        alg = "lbfgs";
    }

    if (alg == "lbfgs" || alg == "newton_raphson" || (allow_irls && alg == "irls")) {
        return alg;
    }
    if (allow_irls) {
        Rcpp::stop("optimization_alg must be one of 'irls', 'lbfgs', or 'newton_raphson'");
    }
    Rcpp::stop("optimization_alg must be one of 'lbfgs' or 'newton_raphson'");
}

template <typename LikelihoodFunctor>
inline double likelihood_value(LikelihoodFunctor& fun,
                               const Eigen::VectorXd& params) {
    Eigen::VectorXd grad(params.size());
    return fun(params, grad);
}

template <typename FullFunctor>
class FixedParameterFunctor {
private:
    FullFunctor& m_fun;
    const FixedParamSpec& m_spec;
    const Eigen::VectorXd& m_full_template;

public:
    FixedParameterFunctor(FullFunctor& fun,
                          const FixedParamSpec& spec,
                          const Eigen::VectorXd& full_template) :
        m_fun(fun), m_spec(spec), m_full_template(full_template) {}

    double operator()(const Eigen::VectorXd& free_params, Eigen::VectorXd& grad_free) {
        Eigen::VectorXd full_params = expand_free_params(free_params, m_full_template, m_spec);
        Eigen::VectorXd grad_full(full_params.size());
        double val = m_fun(full_params, grad_full);
        grad_free = subset_vector(grad_full, m_spec.free_idx);
        return val;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& free_params) {
        Eigen::VectorXd full_params = expand_free_params(free_params, m_full_template, m_spec);
        Eigen::MatrixXd H_full = m_fun.hessian(full_params);
        return subset_matrix(H_full, m_spec.free_idx, m_spec.free_idx);
    }

    Eigen::VectorXd expand(const Eigen::VectorXd& free_params) const {
        return expand_free_params(free_params, m_full_template, m_spec);
    }
};

template <typename LikelihoodFunctor>
inline LikelihoodFitResult optimize_likelihood_lbfgs(LikelihoodFunctor& fun,
                                                     Eigen::VectorXd params,
                                                     int maxit,
                                                     double tol,
                                                     int max_linesearch = 0) {
    LBFGSpp::LBFGSParam<double> lbfgs_params;
    lbfgs_params.epsilon = tol;
    lbfgs_params.max_iterations = maxit;
    if (max_linesearch > 0) lbfgs_params.max_linesearch = max_linesearch;

    LBFGSpp::LBFGSSolver<double> solver(lbfgs_params);
    LikelihoodFitResult fit;
    fit.params = params;
    fit.niter = solver.minimize(fun, fit.params, fit.value);
    fit.converged = (fit.niter < maxit);
    return fit;
}

template <typename LikelihoodFunctor>
inline LikelihoodFitResult optimize_likelihood_newton(LikelihoodFunctor& fun,
                                                      Eigen::VectorXd params,
                                                      int maxit,
                                                      double tol) {
    LikelihoodFitResult fit;
    fit.params = params;

    for (int iter = 0; iter < maxit; ++iter) {
        Eigen::VectorXd grad(params.size());
        double current_value = fun(params, grad);
        if (!std::isfinite(current_value) || !grad.allFinite()) break;
        if (grad.norm() < tol) {
            fit.value = current_value;
            fit.niter = iter;
            fit.converged = true;
            fit.params = params;
            return fit;
        }

        Eigen::MatrixXd H = fun.hessian(params);
        if (!H.allFinite()) break;
        Eigen::FullPivLU<Eigen::MatrixXd> lu(H);
        if (!lu.isInvertible()) break;
        Eigen::VectorXd step = lu.solve(grad);
        if (!step.allFinite()) break;

        double step_scale = 1.0;
        bool accepted = false;
        while (step_scale > 1e-4) {
            Eigen::VectorXd candidate = params - step_scale * step;
            double candidate_value = likelihood_value(fun, candidate);
            if (std::isfinite(candidate_value) && candidate_value < current_value) {
                params = candidate;
                fit.value = candidate_value;
                accepted = true;
                break;
            }
            step_scale *= 0.5;
        }
        if (!accepted) break;

        fit.niter = iter + 1;
        if ((step_scale * step).norm() < tol) {
            fit.converged = true;
            fit.params = params;
            return fit;
        }
    }

    fit.params = params;
    fit.value = likelihood_value(fun, params);
    fit.converged = false;
    return fit;
}

template <typename LikelihoodFunctor>
inline LikelihoodFitResult optimize_likelihood_newton_then_lbfgs(LikelihoodFunctor& fun,
                                                                 Eigen::VectorXd params,
                                                                 int maxit,
                                                                 double tol,
                                                                 int max_linesearch = 0) {
    LikelihoodFitResult newton_fit = optimize_likelihood_newton(fun, params, maxit, tol);
    if (newton_fit.converged) {
        return newton_fit;
    }

    try {
        LikelihoodFitResult lbfgs_fit = optimize_likelihood_lbfgs(fun, newton_fit.params, maxit, tol, max_linesearch);
        if (lbfgs_fit.converged) {
            return lbfgs_fit;
        }
        if (std::isfinite(lbfgs_fit.value) &&
            (!std::isfinite(newton_fit.value) || lbfgs_fit.value < newton_fit.value)) {
            return lbfgs_fit;
        }
    } catch (...) {
        // Keep the best damped-Newton result when the fallback optimizer also fails.
    }

    return newton_fit;
}

template <typename FullFunctor>
inline LikelihoodFitResult optimize_fixed_likelihood_lbfgs(FullFunctor& fun,
                                                           Eigen::VectorXd params,
                                                           const FixedParamSpec& fixed_spec,
                                                           int maxit,
                                                           double tol,
                                                           int max_linesearch = 0) {
    params = apply_fixed_values(params, fixed_spec);
    Eigen::VectorXd params_free = subset_vector(params, fixed_spec.free_idx);
    FixedParameterFunctor<FullFunctor> fixed_fun(fun, fixed_spec, params);
    LikelihoodFitResult fit = optimize_likelihood_lbfgs(fixed_fun, params_free, maxit, tol, max_linesearch);
    fit.params = fixed_fun.expand(fit.params);
    return fit;
}

template <typename FullFunctor>
inline LikelihoodFitResult optimize_fixed_likelihood_newton(FullFunctor& fun,
                                                            Eigen::VectorXd params,
                                                            const FixedParamSpec& fixed_spec,
                                                            int maxit,
                                                            double tol) {
    params = apply_fixed_values(params, fixed_spec);
    Eigen::VectorXd params_free = subset_vector(params, fixed_spec.free_idx);
    FixedParameterFunctor<FullFunctor> fixed_fun(fun, fixed_spec, params);
    LikelihoodFitResult fit = optimize_likelihood_newton(fixed_fun, params_free, maxit, tol);
    fit.params = fixed_fun.expand(fit.params);
    return fit;
}

template <typename FullFunctor>
inline LikelihoodFitResult optimize_fixed_likelihood(FullFunctor& fun,
                                                     Eigen::VectorXd params,
                                                     const FixedParamSpec& fixed_spec,
                                                     int maxit,
                                                     double tol,
                                                     const std::string& optimization_alg,
                                                     const std::string& default_optimization_alg,
                                                     int max_linesearch = 0) {
    std::string alg = normalize_optimizer_algorithm(optimization_alg, default_optimization_alg, false);
    if (alg == "lbfgs") {
        return optimize_fixed_likelihood_lbfgs(fun, params, fixed_spec, maxit, tol, max_linesearch);
    }
    params = apply_fixed_values(params, fixed_spec);
    Eigen::VectorXd params_free = subset_vector(params, fixed_spec.free_idx);
    FixedParameterFunctor<FullFunctor> fixed_fun(fun, fixed_spec, params);
    LikelihoodFitResult fit = optimize_likelihood_newton_then_lbfgs(fixed_fun, params_free, maxit, tol, max_linesearch);
    fit.params = fixed_fun.expand(fit.params);
    return fit;
}

template <typename LikelihoodFunctor>
inline LikelihoodFitResult optimize_likelihood(LikelihoodFunctor& fun,
                                               Eigen::VectorXd params,
                                               int maxit,
                                               double tol,
                                               const std::string& optimization_alg,
                                               const std::string& default_optimization_alg,
                                               int max_linesearch = 0) {
    std::string alg = normalize_optimizer_algorithm(optimization_alg, default_optimization_alg, false);
    if (alg == "lbfgs") {
        return optimize_likelihood_lbfgs(fun, params, maxit, tol, max_linesearch);
    }
    return optimize_likelihood_newton_then_lbfgs(fun, params, maxit, tol, max_linesearch);
}

#endif
