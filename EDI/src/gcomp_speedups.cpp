// [[Rcpp::depends(RcppEigen)]]
#include "_helper_functions.h"
#include <unordered_map>

using namespace Rcpp;

namespace {

inline double plogis_stable(double x) {
  if (x >= 0.0) {
    const double z = std::exp(-x);
    return 1.0 / (1.0 + z);
  }
  const double z = std::exp(x);
  return z / (1.0 + z);
}

inline Eigen::ArrayXd plogis_array(const Eigen::ArrayXd& eta) {
  const Eigen::Array<bool, Eigen::Dynamic, 1> nonnegative = (eta >= 0.0);
  const Eigen::ArrayXd pos = 1.0 / (1.0 + (-eta).exp());
  const Eigen::ArrayXd neg_exp = eta.exp();
  const Eigen::ArrayXd neg = neg_exp / (1.0 + neg_exp);
  return nonnegative.select(pos, neg);
}

Eigen::MatrixXd cluster_meat(const Eigen::MatrixXd& X_fit,
                             const Eigen::VectorXd& resid,
                             const IntegerVector& cluster_id) {
  const int n = X_fit.rows();
  const int p = X_fit.cols();
  if (cluster_id.size() != n) {
    stop("dimension mismatch in cluster_meat");
  }

  Eigen::MatrixXd meat = Eigen::MatrixXd::Zero(p, p);
  std::unordered_map<int, Eigen::VectorXd> cluster_scores;
  cluster_scores.reserve(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    const int id = cluster_id[i];
    auto it = cluster_scores.find(id);
    if (it == cluster_scores.end()) {
      it = cluster_scores.emplace(id, Eigen::VectorXd::Zero(p)).first;
    }
    it->second.noalias() += X_fit.row(i).transpose() * resid[i];
  }
  for (const auto& entry : cluster_scores) {
    const auto& score_g = entry.second;
    meat.noalias() += score_g * score_g.transpose();
  }
  return meat;
}

}  // namespace

// [[Rcpp::export]]
List gcomp_logistic_post_fit_cpp(const Eigen::MatrixXd& X_fit,
                                 const Eigen::VectorXd& y,
                                 const Eigen::VectorXd& coef_hat,
                                 const Eigen::VectorXd& mu_hat,
                                 int j_treat) {
  const int n = X_fit.rows();
  const int p = X_fit.cols();
  const int j_treat0 = j_treat - 1;

  if (j_treat0 < 0 || j_treat0 >= p) {
    stop("treatment column index is out of bounds");
  }
  if (y.size() != n || coef_hat.size() != p || mu_hat.size() != n) {
    stop("dimension mismatch in gcomp_logistic_post_fit_cpp");
  }

  for (int i = 0; i < n; ++i) {
    const double mu_i = mu_hat[i];
    if (!R_finite(mu_i) || mu_i <= 0.0 || mu_i >= 1.0) {
      stop("non-finite or boundary fitted values");
    }
  }
  const Eigen::VectorXd W = (mu_hat.array() * (1.0 - mu_hat.array())).matrix();

  const Eigen::MatrixXd XtWX = weighted_crossprod(X_fit, W);
  Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
  if (ldlt.info() != Eigen::Success) {
    stop("failed to factorize X'WX");
  }

  const Eigen::MatrixXd bread = ldlt.solve(Eigen::MatrixXd::Identity(p, p));
  if (ldlt.info() != Eigen::Success) {
    stop("failed to invert X'WX");
  }

  const Eigen::VectorXd score_resid = y - mu_hat;
  const Eigen::VectorXd resid_sq = score_resid.array().square().matrix();
  Eigen::MatrixXd meat = weighted_crossprod(X_fit, resid_sq);
  Eigen::MatrixXd vcov_robust = bread * meat * bread;
  vcov_robust = 0.5 * (vcov_robust + vcov_robust.transpose());

  for (int j = 0; j < p; ++j) {
    for (int k = 0; k < p; ++k) {
      if (!R_finite(vcov_robust(j, k))) {
        stop("non-finite robust covariance");
      }
    }
  }

  const double ssq_treat = vcov_robust(j_treat0, j_treat0);
  if (!R_finite(ssq_treat) || ssq_treat <= 0.0) {
    stop("non-positive treatment variance");
  }

  const Eigen::VectorXd eta = X_fit * coef_hat;
  const Eigen::VectorXd eta_base = eta - coef_hat[j_treat0] * X_fit.col(j_treat0);

  const Eigen::ArrayXd risk1_arr = plogis_array(eta_base.array() + coef_hat[j_treat0]);
  const Eigen::ArrayXd risk0_arr = plogis_array(eta_base.array());
  const Eigen::VectorXd risk1_i = risk1_arr.matrix();
  const Eigen::VectorXd risk0_i = risk0_arr.matrix();

  const double risk1 = risk1_i.mean();
  const double risk0 = risk0_i.mean();

  const double inv_n = 1.0 / static_cast<double>(n);
  const Eigen::VectorXd wt1 = (risk1_arr * (1.0 - risk1_arr) * inv_n).matrix();
  const Eigen::VectorXd wt0 = (risk0_arr * (1.0 - risk0_arr) * inv_n).matrix();
  Eigen::VectorXd grad1 = X_fit.transpose() * wt1;
  Eigen::VectorXd grad0 = X_fit.transpose() * wt0;
  grad1[j_treat0] = wt1.sum();
  grad0[j_treat0] = 0.0;

  const double rd = risk1 - risk0;
  const Eigen::VectorXd grad_rd = grad1 - grad0;
  const Eigen::VectorXd bread_grad_rd = ldlt.solve(grad_rd);
  const double var_rd = (bread_grad_rd.transpose() * meat * bread_grad_rd)(0, 0);
  const double se_rd = (R_finite(var_rd) && var_rd >= 0.0) ? std::sqrt(var_rd) : NA_REAL;

  double log_rr = NA_REAL;
  double rr = NA_REAL;
  double se_log_rr = NA_REAL;
  if (risk1 > 0.0 && risk0 > 0.0) {
    log_rr = std::log(risk1) - std::log(risk0);
    rr = std::exp(log_rr);
    const Eigen::VectorXd grad_log_rr = grad1 / risk1 - grad0 / risk0;
    const Eigen::VectorXd bread_grad_log_rr = ldlt.solve(grad_log_rr);
    const double var_log_rr = (bread_grad_log_rr.transpose() * meat * bread_grad_log_rr)(0, 0);
    if (R_finite(var_log_rr) && var_log_rr >= 0.0) {
      se_log_rr = std::sqrt(var_log_rr);
    }
  }

  Eigen::VectorXd std_err(p);
  Eigen::VectorXd z_vals(p);
  for (int j = 0; j < p; ++j) {
    const double var_j = vcov_robust(j, j);
    std_err[j] = (R_finite(var_j) && var_j >= 0.0) ? std::sqrt(var_j) : NA_REAL;
    z_vals[j] = (R_finite(std_err[j]) && std_err[j] > 0.0) ? coef_hat[j] / std_err[j] : NA_REAL;
  }

  return List::create(
    _["vcov"] = vcov_robust,
    _["std_err"] = std_err,
    _["z_vals"] = z_vals,
    _["risk1"] = risk1,
    _["risk0"] = risk0,
    _["rd"] = rd,
    _["se_rd"] = se_rd,
    _["log_rr"] = log_rr,
    _["rr"] = rr,
    _["se_log_rr"] = se_log_rr
  );
}

// [[Rcpp::export]]
List gcomp_fractional_logit_post_fit_cpp(const Eigen::MatrixXd& X_fit,
                                         const Eigen::VectorXd& y,
                                         const Eigen::VectorXd& coef_hat,
                                         const Eigen::VectorXd& mu_hat,
                                         int j_treat) {
  List base = gcomp_logistic_post_fit_cpp(X_fit, y, coef_hat, mu_hat, j_treat);
  return List::create(
    _["vcov"] = base["vcov"],
    _["std_err"] = base["std_err"],
    _["z_vals"] = base["z_vals"],
    _["mean1"] = base["risk1"],
    _["mean0"] = base["risk0"],
    _["md"] = base["rd"],
    _["se_md"] = base["se_rd"]
  );
}

// [[Rcpp::export]]
List gcomp_logistic_cluster_post_fit_cpp(const Eigen::MatrixXd& X_fit,
                                         const Eigen::VectorXd& y,
                                         const Eigen::VectorXd& coef_hat,
                                         const Eigen::VectorXd& mu_hat,
                                         const IntegerVector& cluster_id,
                                         int j_treat) {
  const int n = X_fit.rows();
  const int p = X_fit.cols();
  const int j_treat0 = j_treat - 1;

  if (j_treat0 < 0 || j_treat0 >= p) {
    stop("treatment column index is out of bounds");
  }
  if (y.size() != n || coef_hat.size() != p || mu_hat.size() != n) {
    stop("dimension mismatch in gcomp_logistic_cluster_post_fit_cpp");
  }

  for (int i = 0; i < n; ++i) {
    const double mu_i = mu_hat[i];
    if (!R_finite(mu_i) || mu_i <= 0.0 || mu_i >= 1.0) {
      stop("non-finite or boundary fitted values");
    }
  }
  const Eigen::VectorXd W = (mu_hat.array() * (1.0 - mu_hat.array())).matrix();

  const Eigen::MatrixXd XtWX = weighted_crossprod(X_fit, W);
  Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
  if (ldlt.info() != Eigen::Success) {
    stop("failed to factorize X'WX");
  }

  const Eigen::MatrixXd bread = ldlt.solve(Eigen::MatrixXd::Identity(p, p));
  if (ldlt.info() != Eigen::Success) {
    stop("failed to invert X'WX");
  }

  const Eigen::VectorXd score_resid = y - mu_hat;
  Eigen::MatrixXd meat = cluster_meat(X_fit, score_resid, cluster_id);
  Eigen::MatrixXd vcov_robust = bread * meat * bread;
  vcov_robust = 0.5 * (vcov_robust + vcov_robust.transpose());

  for (int j = 0; j < p; ++j) {
    for (int k = 0; k < p; ++k) {
      if (!R_finite(vcov_robust(j, k))) {
        stop("non-finite robust covariance");
      }
    }
  }

  const double ssq_treat = vcov_robust(j_treat0, j_treat0);
  if (!R_finite(ssq_treat) || ssq_treat <= 0.0) {
    stop("non-positive treatment variance");
  }

  const Eigen::VectorXd eta = X_fit * coef_hat;
  const Eigen::VectorXd eta_base = eta - coef_hat[j_treat0] * X_fit.col(j_treat0);

  const Eigen::ArrayXd risk1_arr = plogis_array(eta_base.array() + coef_hat[j_treat0]);
  const Eigen::ArrayXd risk0_arr = plogis_array(eta_base.array());
  const Eigen::VectorXd risk1_i = risk1_arr.matrix();
  const Eigen::VectorXd risk0_i = risk0_arr.matrix();

  const double risk1 = risk1_i.mean();
  const double risk0 = risk0_i.mean();

  const double inv_n = 1.0 / static_cast<double>(n);
  const Eigen::VectorXd wt1 = (risk1_arr * (1.0 - risk1_arr) * inv_n).matrix();
  const Eigen::VectorXd wt0 = (risk0_arr * (1.0 - risk0_arr) * inv_n).matrix();
  Eigen::VectorXd grad1 = X_fit.transpose() * wt1;
  Eigen::VectorXd grad0 = X_fit.transpose() * wt0;
  grad1[j_treat0] = wt1.sum();
  grad0[j_treat0] = 0.0;

  const double rd = risk1 - risk0;
  const Eigen::VectorXd grad_rd = grad1 - grad0;
  const Eigen::VectorXd bread_grad_rd = ldlt.solve(grad_rd);
  const double var_rd = (bread_grad_rd.transpose() * meat * bread_grad_rd)(0, 0);
  const double se_rd = (R_finite(var_rd) && var_rd >= 0.0) ? std::sqrt(var_rd) : NA_REAL;

  double log_rr = NA_REAL;
  double rr = NA_REAL;
  double se_log_rr = NA_REAL;
  if (risk1 > 0.0 && risk0 > 0.0) {
    log_rr = std::log(risk1) - std::log(risk0);
    rr = std::exp(log_rr);
    const Eigen::VectorXd grad_log_rr = grad1 / risk1 - grad0 / risk0;
    const Eigen::VectorXd bread_grad_log_rr = ldlt.solve(grad_log_rr);
    const double var_log_rr = (bread_grad_log_rr.transpose() * meat * bread_grad_log_rr)(0, 0);
    if (R_finite(var_log_rr) && var_log_rr >= 0.0) {
      se_log_rr = std::sqrt(var_log_rr);
    }
  }

  Eigen::VectorXd std_err(p);
  Eigen::VectorXd z_vals(p);
  for (int j = 0; j < p; ++j) {
    const double var_j = vcov_robust(j, j);
    std_err[j] = (R_finite(var_j) && var_j >= 0.0) ? std::sqrt(var_j) : NA_REAL;
    z_vals[j] = (R_finite(std_err[j]) && std_err[j] > 0.0) ? coef_hat[j] / std_err[j] : NA_REAL;
  }

  return List::create(
    _["vcov"] = vcov_robust,
    _["std_err"] = std_err,
    _["z_vals"] = z_vals,
    _["risk1"] = risk1,
    _["risk0"] = risk0,
    _["rd"] = rd,
    _["se_rd"] = se_rd,
    _["log_rr"] = log_rr,
    _["rr"] = rr,
    _["se_log_rr"] = se_log_rr
  );
}

// [[Rcpp::export]]
List gcomp_ordinal_proportional_odds_post_fit_cpp(const Eigen::MatrixXd& X_fit,
                                                  const Eigen::VectorXd& coef_hat,
                                                  const Eigen::VectorXd& alpha_hat,
                                                  int j_treat) {
  const int n = X_fit.rows();
  const int K_minus_1 = alpha_hat.size();
  const int j_treat0 = j_treat - 1;

  Eigen::VectorXd eta_base = X_fit * coef_hat - coef_hat[j_treat0] * X_fit.col(j_treat0);
  Eigen::VectorXd eta1 = eta_base.array() + coef_hat[j_treat0];
  Eigen::VectorXd eta0 = eta_base;

  auto compute_mean = [&](const Eigen::VectorXd& eta_vec) {
    Eigen::ArrayXd mean = Eigen::ArrayXd::Ones(n);
    for (int k = 0; k < K_minus_1; ++k) {
      mean += 1.0 - plogis_array(Eigen::ArrayXd::Constant(n, alpha_hat[k]) - eta_vec.array());
    }
    return mean.mean();
  };

  double mean1 = compute_mean(eta1);
  double mean0 = compute_mean(eta0);

  return List::create(
    _["mean1"] = mean1,
    _["mean0"] = mean0,
    _["md"] = mean1 - mean0
  );
}
