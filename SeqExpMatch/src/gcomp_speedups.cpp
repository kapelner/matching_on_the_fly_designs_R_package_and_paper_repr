// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>
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

Eigen::MatrixXd cluster_meat(const Eigen::MatrixXd& X_fit,
                             const Eigen::VectorXd& resid,
                             const IntegerVector& cluster_id) {
  const int n = X_fit.rows();
  const int p = X_fit.cols();
  if (cluster_id.size() != n) {
    stop("dimension mismatch in cluster_meat");
  }

  std::unordered_map<int, int> cluster_lookup;
  std::vector<Eigen::VectorXd> cluster_scores;
  cluster_scores.reserve(static_cast<std::size_t>(n));

  for (int i = 0; i < n; ++i) {
    const int id = cluster_id[i];
    auto it = cluster_lookup.find(id);
    int pos;
    if (it == cluster_lookup.end()) {
      pos = static_cast<int>(cluster_scores.size());
      cluster_lookup.emplace(id, pos);
      cluster_scores.push_back(Eigen::VectorXd::Zero(p));
    } else {
      pos = it->second;
    }
    cluster_scores[static_cast<std::size_t>(pos)].noalias() += X_fit.row(i).transpose() * resid[i];
  }

  Eigen::MatrixXd meat = Eigen::MatrixXd::Zero(p, p);
  for (const auto& score_g : cluster_scores) {
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

  Eigen::VectorXd W(n);
  for (int i = 0; i < n; ++i) {
    const double mu_i = mu_hat[i];
    if (!R_finite(mu_i) || mu_i <= 0.0 || mu_i >= 1.0) {
      stop("non-finite or boundary fitted values");
    }
    W[i] = mu_i * (1.0 - mu_i);
  }

  const Eigen::MatrixXd XtWX = X_fit.transpose() * W.asDiagonal() * X_fit;
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
  Eigen::MatrixXd meat = X_fit.transpose() * resid_sq.asDiagonal() * X_fit;
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

  Eigen::VectorXd risk1_i(n);
  Eigen::VectorXd risk0_i(n);
  for (int i = 0; i < n; ++i) {
    risk1_i[i] = plogis_stable(eta_base[i] + coef_hat[j_treat0]);
    risk0_i[i] = plogis_stable(eta_base[i]);
  }

  const double risk1 = risk1_i.mean();
  const double risk0 = risk0_i.mean();

  Eigen::VectorXd grad1 = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd grad0 = Eigen::VectorXd::Zero(p);
  const double inv_n = 1.0 / static_cast<double>(n);
  for (int i = 0; i < n; ++i) {
    const double wt1 = risk1_i[i] * (1.0 - risk1_i[i]) * inv_n;
    const double wt0 = risk0_i[i] * (1.0 - risk0_i[i]) * inv_n;
    for (int j = 0; j < p; ++j) {
      const double x_ij = X_fit(i, j);
      grad1[j] += (j == j_treat0 ? 1.0 : x_ij) * wt1;
      grad0[j] += (j == j_treat0 ? 0.0 : x_ij) * wt0;
    }
  }

  const double rd = risk1 - risk0;
  const Eigen::VectorXd grad_rd = grad1 - grad0;
  const double var_rd = (grad_rd.transpose() * vcov_robust * grad_rd)(0, 0);
  const double se_rd = (R_finite(var_rd) && var_rd >= 0.0) ? std::sqrt(var_rd) : NA_REAL;

  double log_rr = NA_REAL;
  double rr = NA_REAL;
  double se_log_rr = NA_REAL;
  if (risk1 > 0.0 && risk0 > 0.0) {
    log_rr = std::log(risk1) - std::log(risk0);
    rr = std::exp(log_rr);
    const Eigen::VectorXd grad_log_rr = grad1 / risk1 - grad0 / risk0;
    const double var_log_rr = (grad_log_rr.transpose() * vcov_robust * grad_log_rr)(0, 0);
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

  Eigen::VectorXd W(n);
  for (int i = 0; i < n; ++i) {
    const double mu_i = mu_hat[i];
    if (!R_finite(mu_i) || mu_i <= 0.0 || mu_i >= 1.0) {
      stop("non-finite or boundary fitted values");
    }
    W[i] = mu_i * (1.0 - mu_i);
  }

  const Eigen::MatrixXd XtWX = X_fit.transpose() * W.asDiagonal() * X_fit;
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

  Eigen::VectorXd risk1_i(n);
  Eigen::VectorXd risk0_i(n);
  for (int i = 0; i < n; ++i) {
    risk1_i[i] = plogis_stable(eta_base[i] + coef_hat[j_treat0]);
    risk0_i[i] = plogis_stable(eta_base[i]);
  }

  const double risk1 = risk1_i.mean();
  const double risk0 = risk0_i.mean();

  Eigen::VectorXd grad1 = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd grad0 = Eigen::VectorXd::Zero(p);
  const double inv_n = 1.0 / static_cast<double>(n);
  for (int i = 0; i < n; ++i) {
    const double wt1 = risk1_i[i] * (1.0 - risk1_i[i]) * inv_n;
    const double wt0 = risk0_i[i] * (1.0 - risk0_i[i]) * inv_n;
    for (int j = 0; j < p; ++j) {
      const double x_ij = X_fit(i, j);
      grad1[j] += (j == j_treat0 ? 1.0 : x_ij) * wt1;
      grad0[j] += (j == j_treat0 ? 0.0 : x_ij) * wt0;
    }
  }

  const double rd = risk1 - risk0;
  const Eigen::VectorXd grad_rd = grad1 - grad0;
  const double var_rd = (grad_rd.transpose() * vcov_robust * grad_rd)(0, 0);
  const double se_rd = (R_finite(var_rd) && var_rd >= 0.0) ? std::sqrt(var_rd) : NA_REAL;

  double log_rr = NA_REAL;
  double rr = NA_REAL;
  double se_log_rr = NA_REAL;
  if (risk1 > 0.0 && risk0 > 0.0) {
    log_rr = std::log(risk1) - std::log(risk0);
    rr = std::exp(log_rr);
    const Eigen::VectorXd grad_log_rr = grad1 / risk1 - grad0 / risk0;
    const double var_log_rr = (grad_log_rr.transpose() * vcov_robust * grad_log_rr)(0, 0);
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
