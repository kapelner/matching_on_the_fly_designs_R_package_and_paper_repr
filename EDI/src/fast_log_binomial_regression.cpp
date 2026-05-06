// [[Rcpp::depends(RcppEigen)]]
#include "_helper_functions.h"
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
                                       double tol,
                                       Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                       Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue) {
  const int n = X.rows();
  const int p = X.cols();
  if (y.size() != n) stop("dimension mismatch in constrained binomial regression");
  FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);
  const int p_free = fixed_spec.free_idx.size();
  Eigen::MatrixXd X_free(n, p_free);
  for (int j = 0; j < p_free; ++j) X_free.col(j) = X.col(fixed_spec.free_idx[j]);
  Eigen::VectorXd eta_fixed = Eigen::VectorXd::Zero(n);
  for (int j = 0; j < fixed_spec.fixed_idx.size(); ++j) {
    eta_fixed.noalias() += X.col(fixed_spec.fixed_idx[j]) * fixed_spec.fixed_values[j];
  }

  Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
  const double y_mean = std::min(std::max(y.mean(), kMinMu), kMaxMu);
  beta[0] = (link_type == BinomialConstrainedLink::kLog) ? std::log(y_mean) : y_mean;
  beta = apply_fixed_values(beta, fixed_spec);
  Eigen::VectorXd beta_free = subset_vector(beta, fixed_spec.free_idx);

  bool converged = false;
  Eigen::VectorXd mu = Eigen::VectorXd::Constant(n, y_mean);
  Eigen::VectorXd w = Eigen::VectorXd::Constant(
    n,
    (link_type == BinomialConstrainedLink::kLog) ?
      y_mean / (1.0 - y_mean) :
      1.0 / (y_mean * (1.0 - y_mean))
  );

  for (int iter = 0; iter < maxit; ++iter) {
    Eigen::VectorXd eta = eta_fixed + X_free * beta_free;
    if (link_type == BinomialConstrainedLink::kLog) {
      eta = eta.array().min(kMaxEtaLog).matrix();
      mu = eta.array().exp().matrix();
      w = (mu.array() / (1.0 - mu.array()).max(kEps)).max(kEps).matrix();
    } else {
      eta = eta.array().max(kMinMu).min(kMaxMu).matrix();
      mu = eta;
      w = (1.0 / (mu.array() * (1.0 - mu.array())).max(kEps)).max(kEps).matrix();
    }

    Eigen::VectorXd z = (link_type == BinomialConstrainedLink::kLog) ?
      eta + (y - mu).cwiseQuotient(mu.array().max(kEps).matrix()) :
      y;
    Eigen::VectorXd z_adj = z - eta_fixed;

    Eigen::MatrixXd XtWX = weighted_crossprod(X_free, w);
    Eigen::VectorXd XtWz = weighted_crossprod_rhs(X_free, w, z_adj);

    Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
    if (ldlt.info() != Eigen::Success) {
      return List::create(_["b"] = beta, _["mu_hat"] = mu, _["working_weights"] = w, _["converged"] = false);
    }

    Eigen::VectorXd beta_free_target = ldlt.solve(XtWz);
    if (ldlt.info() != Eigen::Success || !all_finite_vec(beta_free_target)) {
      return List::create(_["b"] = beta, _["mu_hat"] = mu, _["working_weights"] = w, _["converged"] = false);
    }

    const double ll_curr = loglik_constrained_binomial(X, y, beta, link_type);
    double step = 1.0;
    Eigen::VectorXd beta_new = beta;
    Eigen::VectorXd beta_free_new = beta_free;
    bool accepted = false;
    while (step >= 1e-8) {
      beta_free_new = beta_free + step * (beta_free_target - beta_free);
      beta_new = expand_free_params(beta_free_new, beta, fixed_spec);
      const double ll_new = loglik_constrained_binomial(X, y, beta_new, link_type);
      if (R_finite(ll_new) && ll_new >= ll_curr - 1e-10) {
        accepted = true;
        break;
      }
      step *= 0.5;
    }
    if (!accepted) break;
    if ((beta_free_new - beta_free).norm() < tol) {
      beta = beta_new;
      beta_free = beta_free_new;
      converged = true;
      break;
    }
    beta = beta_new;
    beta_free = beta_free_new;
  }

  Eigen::VectorXd eta = X * beta;
  if (link_type == BinomialConstrainedLink::kLog) {
    eta = eta.array().min(kMaxEtaLog).matrix();
    mu = eta.array().exp().matrix();
    w = (mu.array() / (1.0 - mu.array()).max(kEps)).max(kEps).matrix();
  } else {
    eta = eta.array().max(kMinMu).min(kMaxMu).matrix();
    mu = eta;
    w = (1.0 / (mu.array() * (1.0 - mu.array())).max(kEps)).max(kEps).matrix();
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
                                                double tol,
                                                Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                                Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue) {
  List fit = fit_constrained_binomial_cpp_impl(Xmm, y, link_type, maxit, tol, fixed_idx, fixed_values);
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

  FixedParamSpec fixed_spec = make_fixed_param_spec(Xmm.cols(), fixed_idx, fixed_values);
  Eigen::MatrixXd X_free(Xmm.rows(), fixed_spec.free_idx.size());
  for (int col = 0; col < fixed_spec.free_idx.size(); ++col) X_free.col(col) = Xmm.col(fixed_spec.free_idx[col]);
  Eigen::MatrixXd XtWX = weighted_crossprod(X_free, w);
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

  Eigen::MatrixXd vcov_free = ldlt.solve(Eigen::MatrixXd::Identity(fixed_spec.free_idx.size(), fixed_spec.free_idx.size()));
  if (ldlt.info() != Eigen::Success || !all_finite_mat(vcov_free)) {
    return List::create(
      _["b"] = beta,
      _["vcov"] = NumericMatrix(0, 0),
      _["std_err"] = NumericVector(0),
      _["z_vals"] = NumericVector(0),
      _["ssq_b_j"] = NA_REAL,
      _["converged"] = false
    );
  }

  Eigen::MatrixXd vcov = expand_free_covariance(Xmm.cols(), fixed_spec, 0.5 * (vcov_free + vcov_free.transpose()), true);
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

Eigen::VectorXd constrained_binomial_score_cpp_impl(const Eigen::MatrixXd& X,
													const Eigen::VectorXd& y,
													const Eigen::VectorXd& beta,
													BinomialConstrainedLink link_type) {
	const int p = beta.size();
	Eigen::VectorXd score(p);
	const double h = 1e-6;
	for (int j = 0; j < p; ++j) {
		Eigen::VectorXd bp = beta;
		Eigen::VectorXd bm = beta;
		bp[j] += h;
		bm[j] -= h;
		score[j] = (loglik_constrained_binomial(X, y, bp, link_type) -
					loglik_constrained_binomial(X, y, bm, link_type)) / (2.0 * h);
	}
	return score;
}

Eigen::MatrixXd constrained_binomial_hessian_cpp_impl(const Eigen::MatrixXd& X,
													  const Eigen::VectorXd& y,
													  const Eigen::VectorXd& beta,
													  BinomialConstrainedLink link_type) {
	const int p = beta.size();
	Eigen::MatrixXd H(p, p);
	const double h = 1e-4;
	for (int i = 0; i < p; ++i) {
		for (int j = i; j < p; ++j) {
			Eigen::VectorXd bpp = beta; bpp[i] += h; bpp[j] += h;
			Eigen::VectorXd bpm = beta; bpm[i] += h; bpm[j] -= h;
			Eigen::VectorXd bmp = beta; bmp[i] -= h; bmp[j] += h;
			Eigen::VectorXd bmm = beta; bmm[i] -= h; bmm[j] -= h;
			H(i, j) = (loglik_constrained_binomial(X, y, bpp, link_type) -
					   loglik_constrained_binomial(X, y, bpm, link_type) -
					   loglik_constrained_binomial(X, y, bmp, link_type) +
					   loglik_constrained_binomial(X, y, bmm, link_type)) / (4.0 * h * h);
			H(j, i) = H(i, j);
		}
	}
	return H;
}

}  // namespace

//' @title Compute Log-Binomial Regression Score
//' @description Calculates the score vector (gradient of the log-likelihood) for a log-binomial regression model.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param beta A numeric vector of coefficients.
//' @return A numeric vector representing the score.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd get_log_binomial_regression_score_cpp(const Eigen::MatrixXd& X,
													  const Eigen::VectorXd& y,
													  const Eigen::VectorXd& beta) {
	return constrained_binomial_score_cpp_impl(X, y, beta, BinomialConstrainedLink::kLog);
}

//' @title Compute Log-Binomial Regression Hessian
//' @description Calculates the Hessian matrix (second derivatives of the log-likelihood) for a log-binomial regression model.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param beta A numeric vector of coefficients.
//' @return A numeric matrix representing the Hessian.
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd get_log_binomial_regression_hessian_cpp(const Eigen::MatrixXd& X,
														const Eigen::VectorXd& y,
														const Eigen::VectorXd& beta) {
	return constrained_binomial_hessian_cpp_impl(X, y, beta, BinomialConstrainedLink::kLog);
}

//' @title Compute Identity-Binomial Regression Score
//' @description Calculates the score vector for a binomial regression model with an identity link.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param beta A numeric vector of coefficients.
//' @return A numeric vector representing the score.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd get_identity_binomial_regression_score_cpp(const Eigen::MatrixXd& X,
														   const Eigen::VectorXd& y,
														   const Eigen::VectorXd& beta) {
	return constrained_binomial_score_cpp_impl(X, y, beta, BinomialConstrainedLink::kIdentity);
}

//' @title Compute Identity-Binomial Regression Hessian
//' @description Calculates the Hessian matrix for a binomial regression model with an identity link.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param beta A numeric vector of coefficients.
//' @return A numeric matrix representing the Hessian.
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd get_identity_binomial_regression_hessian_cpp(const Eigen::MatrixXd& X,
															 const Eigen::VectorXd& y,
															 const Eigen::VectorXd& beta) {
	return constrained_binomial_hessian_cpp_impl(X, y, beta, BinomialConstrainedLink::kIdentity);
}

//' @title Fast Log-Binomial Regression (C++)
//' @description High-performance log-binomial regression fitting using Fisher scoring.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @return A list containing coefficients and fitted values.
//' @export
// [[Rcpp::export]]
List fast_log_binomial_regression_cpp(const Eigen::MatrixXd& X,
                                      const Eigen::VectorXd& y,
                                      int maxit = 100,
                                      double tol = 1e-8,
                                      Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                      Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue) {
  return fit_constrained_binomial_cpp_impl(X, y, BinomialConstrainedLink::kLog, maxit, tol, fixed_idx, fixed_values);
}

//' @title Fast Log-Binomial Regression with Variance (C++)
//' @description Log-binomial regression with variance-covariance matrix and standard errors.
//' @param Xmm A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param j 1-based index of the parameter for which to return specific variance.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @return A list containing coefficients, vcov, and standard errors.
//' @export
// [[Rcpp::export]]
List fast_log_binomial_regression_with_var_cpp(const Eigen::MatrixXd& Xmm,
                                               const Eigen::VectorXd& y,
                                               int j = 2,
                                               int maxit = 100,
                                               double tol = 1e-8,
                                               Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                               Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue) {
  return fit_constrained_binomial_with_var_cpp_impl(Xmm, y, BinomialConstrainedLink::kLog, j, maxit, tol, fixed_idx, fixed_values);
}

//' @title Fast Identity-Binomial Regression (C++)
//' @description High-performance binomial regression with identity link using Fisher scoring.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @return A list containing coefficients and fitted values.
//' @export
// [[Rcpp::export]]
List fast_identity_binomial_regression_cpp(const Eigen::MatrixXd& X,
                                           const Eigen::VectorXd& y,
                                           int maxit = 100,
                                           double tol = 1e-8,
                                           Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                           Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue) {
  return fit_constrained_binomial_cpp_impl(X, y, BinomialConstrainedLink::kIdentity, maxit, tol, fixed_idx, fixed_values);
}

//' @title Fast Identity-Binomial Regression with Variance (C++)
//' @description Binomial regression with identity link, providing variance-covariance matrix and standard errors.
//' @param Xmm A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param j 1-based index of the parameter for which to return specific variance.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @return A list containing coefficients, vcov, and standard errors.
//' @export
// [[Rcpp::export]]
List fast_identity_binomial_regression_with_var_cpp(const Eigen::MatrixXd& Xmm,
                                                    const Eigen::VectorXd& y,
                                                    int j = 2,
                                                    int maxit = 100,
                                                    double tol = 1e-8,
                                                    Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                                    Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue) {
  return fit_constrained_binomial_with_var_cpp_impl(Xmm, y, BinomialConstrainedLink::kIdentity, j, maxit, tol, fixed_idx, fixed_values);
}
