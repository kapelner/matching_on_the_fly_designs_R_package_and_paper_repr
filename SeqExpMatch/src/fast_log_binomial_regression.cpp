// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>
#include <limits>

using namespace Rcpp;

namespace {

constexpr double kEps = 1e-8;
constexpr double kMinMu = kEps;
constexpr double kMaxMu = 1.0 - kEps;
constexpr double kMaxEtaLog = -kEps;

enum class BinomialConstrainedLink {
  kLog = 0,
  kIdentity = 1
};

inline double clamp_prob(double mu) {
  if (mu <= kMinMu) return kMinMu;
  if (mu >= kMaxMu) return kMaxMu;
  return mu;
}

inline double safe_mu_from_eta(double eta, BinomialConstrainedLink link_type) {
  if (link_type == BinomialConstrainedLink::kLog) {
    if (eta >= kMaxEtaLog) return std::exp(kMaxEtaLog);
    return std::exp(eta);
  }
  return clamp_prob(eta);
}

double loglik_constrained_binomial(const Eigen::MatrixXd& X,
                                   const Eigen::VectorXd& y,
                                   const Eigen::VectorXd& beta,
                                   BinomialConstrainedLink link_type) {
  const int n = X.rows();
  Eigen::VectorXd eta = X * beta;
  double ll = 0.0;
  for (int i = 0; i < n; ++i) {
    if (link_type == BinomialConstrainedLink::kLog) {
      if (eta[i] >= kMaxEtaLog) return R_NegInf;
    } else {
      if (eta[i] <= kMinMu || eta[i] >= kMaxMu) return R_NegInf;
    }
    const double mu = safe_mu_from_eta(eta[i], link_type);
    ll += y[i] * std::log(mu) + (1.0 - y[i]) * std::log1p(-mu);
  }
  return ll;
}

bool all_finite_vec(const Eigen::VectorXd& x) {
  for (int i = 0; i < x.size(); ++i) {
    if (!R_finite(x[i])) return false;
  }
  return true;
}

bool all_finite_mat(const Eigen::MatrixXd& X) {
  for (int j = 0; j < X.cols(); ++j) {
    for (int i = 0; i < X.rows(); ++i) {
      if (!R_finite(X(i, j))) return false;
    }
  }
  return true;
}

List fit_constrained_binomial_cpp_impl(const Eigen::MatrixXd& X,
                                       const Eigen::VectorXd& y,
                                       BinomialConstrainedLink link_type,
                                       int maxit,
                                       double tol) {
  const int n = X.rows();
  const int p = X.cols();
  if (y.size() != n) stop("dimension mismatch in constrained binomial regression");

  Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
  const double y_mean = std::min(std::max(y.mean(), kMinMu), kMaxMu);
  beta[0] = (link_type == BinomialConstrainedLink::kLog) ? std::log(y_mean) : y_mean;

  bool converged = false;
  Eigen::VectorXd mu = Eigen::VectorXd::Constant(n, y_mean);
  Eigen::VectorXd w = Eigen::VectorXd::Constant(
    n,
    (link_type == BinomialConstrainedLink::kLog) ?
      y_mean / (1.0 - y_mean) :
      1.0 / (y_mean * (1.0 - y_mean))
  );

  for (int iter = 0; iter < maxit; ++iter) {
    Eigen::VectorXd eta = X * beta;
    for (int i = 0; i < n; ++i) {
      if (link_type == BinomialConstrainedLink::kLog) {
        if (eta[i] >= kMaxEtaLog) eta[i] = kMaxEtaLog;
      } else {
        eta[i] = clamp_prob(eta[i]);
      }
      mu[i] = safe_mu_from_eta(eta[i], link_type);
      if (link_type == BinomialConstrainedLink::kLog) {
        w[i] = std::max(mu[i] / std::max(1.0 - mu[i], kEps), kEps);
      } else {
        w[i] = std::max(1.0 / std::max(mu[i] * (1.0 - mu[i]), kEps), kEps);
      }
    }

    Eigen::VectorXd z = (link_type == BinomialConstrainedLink::kLog) ?
      eta + (y - mu).cwiseQuotient(mu.array().max(kEps).matrix()) :
      y;

    Eigen::MatrixXd XtW = X.transpose() * w.asDiagonal();
    Eigen::MatrixXd XtWX = XtW * X;
    Eigen::VectorXd XtWz = XtW * z;

    Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
    if (ldlt.info() != Eigen::Success) {
      return List::create(_["b"] = beta, _["mu_hat"] = mu, _["working_weights"] = w, _["converged"] = false);
    }

    Eigen::VectorXd beta_target = ldlt.solve(XtWz);
    if (ldlt.info() != Eigen::Success || !all_finite_vec(beta_target)) {
      return List::create(_["b"] = beta, _["mu_hat"] = mu, _["working_weights"] = w, _["converged"] = false);
    }

    const double ll_curr = loglik_constrained_binomial(X, y, beta, link_type);
    double step = 1.0;
    Eigen::VectorXd beta_new = beta;
    bool accepted = false;
    while (step >= 1e-8) {
      beta_new = beta + step * (beta_target - beta);
      const double ll_new = loglik_constrained_binomial(X, y, beta_new, link_type);
      if (R_finite(ll_new) && ll_new >= ll_curr - 1e-10) {
        accepted = true;
        break;
      }
      step *= 0.5;
    }
    if (!accepted) break;
    if ((beta_new - beta).norm() < tol) {
      beta = beta_new;
      converged = true;
      break;
    }
    beta = beta_new;
  }

  Eigen::VectorXd eta = X * beta;
  for (int i = 0; i < n; ++i) {
    if (link_type == BinomialConstrainedLink::kLog) {
      if (eta[i] >= kMaxEtaLog) eta[i] = kMaxEtaLog;
    } else {
      eta[i] = clamp_prob(eta[i]);
    }
    mu[i] = safe_mu_from_eta(eta[i], link_type);
    if (link_type == BinomialConstrainedLink::kLog) {
      w[i] = std::max(mu[i] / std::max(1.0 - mu[i], kEps), kEps);
    } else {
      w[i] = std::max(1.0 / std::max(mu[i] * (1.0 - mu[i]), kEps), kEps);
    }
  }

  return List::create(
    _["b"] = beta,
    _["mu_hat"] = mu,
    _["working_weights"] = w,
    _["converged"] = converged && all_finite_vec(beta) && all_finite_vec(mu) && all_finite_vec(w)
  );
}

List fit_constrained_binomial_with_var_cpp_impl(const Eigen::MatrixXd& Xmm,
                                                const Eigen::VectorXd& y,
                                                BinomialConstrainedLink link_type,
                                                int j,
                                                int maxit,
                                                double tol) {
  List fit = fit_constrained_binomial_cpp_impl(Xmm, y, link_type, maxit, tol);
  const bool converged = as<bool>(fit["converged"]);
  Eigen::VectorXd beta = fit["b"];
  Eigen::VectorXd w = fit["working_weights"];

  if (!converged || !all_finite_vec(beta) || !all_finite_vec(w)) {
    return List::create(
      _["b"] = beta,
      _["vcov"] = NumericMatrix(0, 0),
      _["std_err"] = NumericVector(0),
      _["z_vals"] = NumericVector(0),
      _["ssq_b_j"] = NA_REAL,
      _["converged"] = false
    );
  }

  Eigen::MatrixXd XtWX = Xmm.transpose() * w.asDiagonal() * Xmm;
  Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
  if (ldlt.info() != Eigen::Success) {
    return List::create(
      _["b"] = beta,
      _["vcov"] = NumericMatrix(0, 0),
      _["std_err"] = NumericVector(0),
      _["z_vals"] = NumericVector(0),
      _["ssq_b_j"] = NA_REAL,
      _["converged"] = false
    );
  }

  Eigen::MatrixXd vcov = ldlt.solve(Eigen::MatrixXd::Identity(Xmm.cols(), Xmm.cols()));
  if (ldlt.info() != Eigen::Success || !all_finite_mat(vcov)) {
    return List::create(
      _["b"] = beta,
      _["vcov"] = NumericMatrix(0, 0),
      _["std_err"] = NumericVector(0),
      _["z_vals"] = NumericVector(0),
      _["ssq_b_j"] = NA_REAL,
      _["converged"] = false
    );
  }

  vcov = 0.5 * (vcov + vcov.transpose());
  Eigen::VectorXd std_err(Xmm.cols());
  Eigen::VectorXd z_vals(Xmm.cols());
  for (int k = 0; k < Xmm.cols(); ++k) {
    const double var_k = vcov(k, k);
    std_err[k] = (R_finite(var_k) && var_k >= 0.0) ? std::sqrt(var_k) : NA_REAL;
    z_vals[k] = (R_finite(std_err[k]) && std_err[k] > 0.0) ? beta[k] / std_err[k] : NA_REAL;
  }

  const int j0 = j - 1;
  const double ssq_b_j = (j0 >= 0 && j0 < vcov.rows()) ? vcov(j0, j0) : NA_REAL;

  return List::create(
    _["b"] = beta,
    _["vcov"] = vcov,
    _["std_err"] = std_err,
    _["z_vals"] = z_vals,
    _["ssq_b_j"] = ssq_b_j,
    _["converged"] = true
  );
}

}  // namespace

// [[Rcpp::export]]
List fast_log_binomial_regression_cpp(const Eigen::MatrixXd& X,
                                      const Eigen::VectorXd& y,
                                      int maxit = 100,
                                      double tol = 1e-8) {
  return fit_constrained_binomial_cpp_impl(X, y, BinomialConstrainedLink::kLog, maxit, tol);
}

// [[Rcpp::export]]
List fast_log_binomial_regression_with_var_cpp(const Eigen::MatrixXd& Xmm,
                                               const Eigen::VectorXd& y,
                                               int j = 2,
                                               int maxit = 100,
                                               double tol = 1e-8) {
  return fit_constrained_binomial_with_var_cpp_impl(Xmm, y, BinomialConstrainedLink::kLog, j, maxit, tol);
}

// [[Rcpp::export]]
List fast_identity_binomial_regression_cpp(const Eigen::MatrixXd& X,
                                           const Eigen::VectorXd& y,
                                           int maxit = 100,
                                           double tol = 1e-8) {
  return fit_constrained_binomial_cpp_impl(X, y, BinomialConstrainedLink::kIdentity, maxit, tol);
}

// [[Rcpp::export]]
List fast_identity_binomial_regression_with_var_cpp(const Eigen::MatrixXd& Xmm,
                                                    const Eigen::VectorXd& y,
                                                    int j = 2,
                                                    int maxit = 100,
                                                    double tol = 1e-8) {
  return fit_constrained_binomial_with_var_cpp_impl(Xmm, y, BinomialConstrainedLink::kIdentity, j, maxit, tol);
}
