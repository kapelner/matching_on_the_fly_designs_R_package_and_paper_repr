#include "_helper_functions.h"
#include <RcppEigen.h>
#include <cmath>
#include <limits>

using namespace Rcpp;

// Forward declarations for functions from RcppExports
Rcpp::List fast_ols_with_var_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int j);

// [[Rcpp::export]]
NumericVector kk21_continuous_weights_cpp(const NumericMatrix& X,
                                          const NumericVector& y) {
  int n = X.nrow();
  int p = X.ncol();
  NumericVector weights(p);

  if (n < 2 || p == 0) {
    std::fill(weights.begin(), weights.end(), std::numeric_limits<double>::epsilon());
    return weights;
  }

  double sumy = 0.0;
  double sumy2 = 0.0;
  for (int i = 0; i < n; ++i) {
    double yi = y[i];
    sumy += yi;
    sumy2 += yi * yi;
  }
  double ybar = sumy / static_cast<double>(n);
  double eps = std::numeric_limits<double>::epsilon();

  for (int j = 0; j < p; ++j) {
    double sumx = 0.0;
    double sumx2 = 0.0;
    double sumxy = 0.0;
    for (int i = 0; i < n; ++i) {
      double x = X(i, j);
      sumx += x;
      sumx2 += x * x;
      sumxy += x * y[i];
    }
    double xbar = sumx / static_cast<double>(n);
    double varx = sumx2 - static_cast<double>(n) * xbar * xbar;

    if (varx <= eps || n <= 2) {
      weights[j] = eps;
      continue;
    }

    double b1 = (sumxy - static_cast<double>(n) * xbar * ybar) / varx;
    double b0 = ybar - b1 * xbar;
    double sse = sumy2 + static_cast<double>(n) * b0 * b0 + b1 * b1 * sumx2 +
      2.0 * b0 * b1 * sumx - 2.0 * b0 * sumy - 2.0 * b1 * sumxy;
    double denom = static_cast<double>(n - 2);
    if (denom <= 0.0 || sse <= 0.0) {
      weights[j] = eps;
      continue;
    }
    double sigma2 = sse / denom;
    double se_b1 = std::sqrt(sigma2 / varx);
    if (!R_finite(se_b1) || se_b1 <= 0.0) {
      weights[j] = eps;
    } else {
      weights[j] = std::fabs(b1 / se_b1);
    }
  }

  return weights;
}

// [[Rcpp::export]]
NumericVector kk21_logistic_weights_cpp(const NumericMatrix& X,
                                        const NumericVector& y,
                                        const int maxit = 100,
                                        const double tol = 1e-8) {
  int n = X.nrow();
  int p = X.ncol();
  NumericVector weights(p);
  double eps = std::numeric_limits<double>::epsilon();

  if (n == 0 || p == 0) {
    std::fill(weights.begin(), weights.end(), eps);
    return weights;
  }

  for (int j = 0; j < p; ++j) {
    double b0 = 0.0;
    double b1 = 0.0;
    bool singular = false;

    for (int iter = 0; iter < maxit; ++iter) {
      double S0 = 0.0;
      double S1 = 0.0;
      double S2 = 0.0;
      double Sz0 = 0.0;
      double Sz1 = 0.0;

      for (int i = 0; i < n; ++i) {
        double x = X(i, j);
        double eta = b0 + b1 * x;
        double p_i = 1.0 / (1.0 + std::exp(-eta));
        double w = p_i * (1.0 - p_i);
        if (w < 1e-10) {
          w = 1e-10;
        }
        double z = eta + (y[i] - p_i) / w;

        S0 += w;
        S1 += w * x;
        S2 += w * x * x;
        Sz0 += w * z;
        Sz1 += w * x * z;
      }

      double det = S0 * S2 - S1 * S1;
      if (det <= eps) {
        singular = true;
        break;
      }

      double b0_new = (Sz0 * S2 - Sz1 * S1) / det;
      double b1_new = (S0 * Sz1 - S1 * Sz0) / det;

      if (std::fabs(b0_new - b0) < tol && std::fabs(b1_new - b1) < tol) {
        b0 = b0_new;
        b1 = b1_new;
        break;
      }
      b0 = b0_new;
      b1 = b1_new;
    }

    if (singular) {
      weights[j] = eps;
      continue;
    }

    double S0 = 0.0;
    double S1 = 0.0;
    double S2 = 0.0;
    for (int i = 0; i < n; ++i) {
      double x = X(i, j);
      double eta = b0 + b1 * x;
      double p_i = 1.0 / (1.0 + std::exp(-eta));
      double w = p_i * (1.0 - p_i);
      if (w < 1e-10) {
        w = 1e-10;
      }
      S0 += w;
      S1 += w * x;
      S2 += w * x * x;
    }

    double det = S0 * S2 - S1 * S1;
    if (det <= eps || S0 <= 0.0) {
      weights[j] = eps;
      continue;
    }
    double var_b1 = S0 / det;
    double se_b1 = std::sqrt(var_b1);
    if (!R_finite(se_b1) || se_b1 <= 0.0) {
      weights[j] = eps;
    } else {
      weights[j] = std::fabs(b1 / se_b1);
    }
  }

  return weights;
}

// [[Rcpp::export]]
NumericVector kk21_stepwise_continuous_weights_cpp(const NumericMatrix& X,
                                                   const NumericVector& y,
                                                   const NumericVector& w) {
  int n = X.nrow();
  int p = X.ncol();
  NumericVector weights(p, NA_REAL);

  if (n == 0 || p == 0) {
    return weights;
  }

  Eigen::Map<Eigen::MatrixXd> X_map(as<Eigen::Map<Eigen::MatrixXd>>(X));
  Eigen::VectorXd y_vec = as<Eigen::VectorXd>(y);
  Eigen::VectorXd w_vec = as<Eigen::VectorXd>(w);

  std::vector<int> selected;
  selected.reserve(p);
  std::vector<bool> used(p, false);

  for (int step = 0; step < p; ++step) {
    double best_stat = -1.0;
    int best_j = -1;
    int k = static_cast<int>(selected.size());

    for (int j = 0; j < p; ++j) {
      if (used[j]) {
        continue;
      }

      Eigen::MatrixXd Xmm(n, k + 3);
      Xmm.col(0).setOnes();
      Xmm.col(1) = X_map.col(j);
      for (int idx = 0; idx < k; ++idx) {
        Xmm.col(2 + idx) = X_map.col(selected[idx]);
      }
      Xmm.col(k + 2) = w_vec;

      List mod = fast_ols_with_var_cpp(Xmm, y_vec, 2);
      Eigen::VectorXd b = mod["b"];
      double ssq = as<double>(mod["ssq_b_j"]);
      double stat = std::fabs(b(1) / std::sqrt(ssq));
      if (!R_finite(stat)) {
        stat = 0.0;
      }
      if (stat > best_stat) {
        best_stat = stat;
        best_j = j;
      }
    }

    if (best_j < 0) {
      break;
    }

    weights[best_j] = best_stat;
    used[best_j] = true;
    selected.push_back(best_j);
  }

  return weights;
}

// [[Rcpp::export]]
NumericVector kk21_stepwise_logistic_weights_cpp(const NumericMatrix& X,
                                                 const NumericVector& y,
                                                 const NumericVector& w) {
  int n = X.nrow();
  int p = X.ncol();
  NumericVector weights(p, NA_REAL);

  if (n == 0 || p == 0) {
    return weights;
  }

  Eigen::Map<Eigen::MatrixXd> X_map(as<Eigen::Map<Eigen::MatrixXd>>(X));
  Eigen::VectorXd y_vec = as<Eigen::VectorXd>(y);
  Eigen::VectorXd w_vec = as<Eigen::VectorXd>(w);

  std::vector<int> selected;
  selected.reserve(p);
  std::vector<bool> used(p, false);

  for (int step = 0; step < p; ++step) {
    double best_stat = -1.0;
    int best_j = -1;
    int k = static_cast<int>(selected.size());

    for (int j = 0; j < p; ++j) {
      if (used[j]) {
        continue;
      }

      Eigen::MatrixXd Xmm(n, k + 3);
      Xmm.col(0).setOnes();
      Xmm.col(1) = X_map.col(j);
      for (int idx = 0; idx < k; ++idx) {
        Xmm.col(2 + idx) = X_map.col(selected[idx]);
      }
      Xmm.col(k + 2) = w_vec;

      List mod = fast_logistic_regression_with_var_cpp(Xmm, y_vec);
      Eigen::VectorXd b = mod["b"];
      double ssq = as<double>(mod["ssq_b_2"]);
      double stat = std::fabs(b(1) / std::sqrt(ssq));
      if (!R_finite(stat)) {
        stat = 0.0;
      }
      if (stat > best_stat) {
        best_stat = stat;
        best_j = j;
      }
    }

    if (best_j < 0) {
      break;
    }

    weights[best_j] = best_stat;
    used[best_j] = true;
    selected.push_back(best_j);
  }

  return weights;
}

// Helper: Univariate beta regression for a single covariate
// Returns t-statistic for the covariate, or -1 if fitting fails
static double univariate_beta_tstat(
    const Eigen::VectorXd& y,
    const Eigen::VectorXd& x,
    int maxit = 100,
    double tol = 1e-8
) {
    int n = y.size();
    if (n < 3) return -1.0;

    // Build design matrix [1, x]
    Eigen::MatrixXd X(n, 2);
    X.col(0).setOnes();
    X.col(1) = x;

    // Fit logistic regression via IRLS (for mean model with logit link)
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(2);
    Eigen::VectorXd prob(n);
    Eigen::VectorXd w(n);

    for (int iter = 0; iter < maxit; ++iter) {
        Eigen::VectorXd eta = X * beta;
        for (int i = 0; i < n; ++i) {
            double val = 1.0 / (1.0 + std::exp(-eta[i]));
            val = std::max(1e-10, std::min(1.0 - 1e-10, val));
            prob[i] = val;
            w[i] = std::max(val * (1.0 - val), 1e-10);
        }

        Eigen::VectorXd z = eta + (y - prob).cwiseQuotient(w);
        Eigen::MatrixXd XtW = X.transpose() * w.asDiagonal();
        Eigen::MatrixXd XtWX = XtW * X;
        Eigen::VectorXd XtWz = XtW * z;

        Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
        if (ldlt.info() != Eigen::Success) return -1.0;

        Eigen::VectorXd beta_new = ldlt.solve(XtWz);
        if ((beta - beta_new).norm() < tol) {
            beta = beta_new;
            break;
        }
        beta = beta_new;
    }

    // Compute final mu and weights
    Eigen::VectorXd eta = X * beta;
    for (int i = 0; i < n; ++i) {
        double val = 1.0 / (1.0 + std::exp(-eta[i]));
        val = std::max(1e-12, std::min(1.0 - 1e-12, val));
        prob[i] = val;
        w[i] = std::max(val * (1.0 - val), 1e-10);
    }

    // Profile likelihood optimization for phi
    // Grid search over log(phi) from log(0.01) to log(1000)
    double best_phi = 10.0;
    double best_ll = -std::numeric_limits<double>::infinity();

    for (double log_phi = -2.0; log_phi <= 7.0; log_phi += 0.5) {
        double phi = std::exp(log_phi);
        double ll = 0.0;
        for (int i = 0; i < n; ++i) {
            double mu_i = prob[i];
            double yi = y[i];
            double mu_phi = mu_i * phi;
            double one_minus_mu_phi = (1.0 - mu_i) * phi;
            mu_phi = std::max(1e-12, mu_phi);
            one_minus_mu_phi = std::max(1e-12, one_minus_mu_phi);

            ll += R::lgammafn(phi) - R::lgammafn(mu_phi) - R::lgammafn(one_minus_mu_phi)
                + (mu_phi - 1.0) * std::log(yi)
                + (one_minus_mu_phi - 1.0) * std::log1p(-yi);
        }
        if (ll > best_ll) {
            best_ll = ll;
            best_phi = phi;
        }
    }

    // Final weights for variance computation: w_final = (1+phi) * mu * (1-mu)
    Eigen::VectorXd w_final(n);
    for (int i = 0; i < n; ++i) {
        w_final[i] = (1.0 + best_phi) * w[i];
    }

    // Compute XtWX and get variance of beta[1]
    Eigen::MatrixXd XtWX = X.transpose() * w_final.asDiagonal() * X;
    Eigen::FullPivLU<Eigen::MatrixXd> lu(XtWX);
    if (!lu.isInvertible()) return -1.0;

    Eigen::MatrixXd cov_mat = lu.inverse();
    if (cov_mat(1, 1) <= 0) return -1.0;

    double se = std::sqrt(cov_mat(1, 1));
    if (!std::isfinite(se) || se <= 0) return -1.0;

    return std::fabs(beta(1) / se);
}

// Helper: Univariate negative binomial regression for a single covariate
// Returns t-statistic for the covariate, or -1 if fitting fails
static double univariate_negbin_tstat(
    const Eigen::VectorXd& y,
    const Eigen::VectorXd& x,
    int maxit = 50,
    double tol = 1e-6
) {
    int n = y.size();
    if (n < 3) return -1.0;

    // Build design matrix [1, x]
    Eigen::MatrixXd X(n, 2);
    X.col(0).setOnes();
    X.col(1) = x;

    // Initialize with Poisson regression coefficients (theta -> infinity)
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(2);

    // Initialize theta using method of moments
    double y_mean = y.mean();
    double y_var = (y.array() - y_mean).square().sum() / (n - 1);
    double theta = (y_var > y_mean) ? (y_mean * y_mean) / (y_var - y_mean) : 10.0;
    theta = std::max(0.01, std::min(1000.0, theta));

    Eigen::VectorXd mu(n);
    Eigen::VectorXd w(n);

    for (int outer = 0; outer < maxit; ++outer) {
        // E-step / IRLS for beta given theta
        for (int iter = 0; iter < 25; ++iter) {
            Eigen::VectorXd eta = X * beta;
            for (int i = 0; i < n; ++i) {
                mu[i] = std::exp(std::min(eta[i], 20.0)); // prevent overflow
                w[i] = mu[i] / (1.0 + mu[i] / theta);
                w[i] = std::max(w[i], 1e-10);
            }

            Eigen::VectorXd z = eta + (y - mu).cwiseQuotient(mu);
            Eigen::MatrixXd XtW = X.transpose() * w.asDiagonal();
            Eigen::MatrixXd XtWX = XtW * X;
            Eigen::VectorXd XtWz = XtW * z;

            Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
            if (ldlt.info() != Eigen::Success) return -1.0;

            Eigen::VectorXd beta_new = ldlt.solve(XtWz);
            if ((beta - beta_new).norm() < tol) {
                beta = beta_new;
                break;
            }
            beta = beta_new;
        }

        // Update mu with final beta
        Eigen::VectorXd eta = X * beta;
        for (int i = 0; i < n; ++i) {
            mu[i] = std::exp(std::min(eta[i], 20.0));
        }

        // M-step: update theta using one-step Newton-Raphson
        double score = 0.0;
        double info = 0.0;
        for (int i = 0; i < n; ++i) {
            double yi = y[i];
            double mui = mu[i];
            score += R::digamma(yi + theta) - R::digamma(theta)
                   + std::log(theta) - std::log(theta + mui) + 1.0
                   - (yi + theta) / (theta + mui);
            info += -R::trigamma(yi + theta) + R::trigamma(theta)
                  - 1.0 / theta + 2.0 / (theta + mui)
                  - (yi + theta) / ((theta + mui) * (theta + mui));
        }

        if (std::fabs(info) < 1e-10) break;
        double theta_new = theta - score / info;
        theta_new = std::max(0.01, std::min(1000.0, theta_new));

        if (std::fabs(theta_new - theta) < tol) {
            theta = theta_new;
            break;
        }
        theta = theta_new;
    }

    // Compute final weights and information matrix
    Eigen::VectorXd eta = X * beta;
    for (int i = 0; i < n; ++i) {
        mu[i] = std::exp(std::min(eta[i], 20.0));
        w[i] = mu[i] / (1.0 + mu[i] / theta);
        w[i] = std::max(w[i], 1e-10);
    }

    Eigen::MatrixXd XtWX = X.transpose() * w.asDiagonal() * X;
    Eigen::FullPivLU<Eigen::MatrixXd> lu(XtWX);
    if (!lu.isInvertible()) return -1.0;

    Eigen::MatrixXd cov_mat = lu.inverse();
    if (cov_mat(1, 1) <= 0) return -1.0;

    double se = std::sqrt(cov_mat(1, 1));
    if (!std::isfinite(se) || se <= 0) return -1.0;

    return std::fabs(beta(1) / se);
}

// [[Rcpp::export]]
NumericVector kk21_beta_weights_cpp(const NumericMatrix& X,
                                    const NumericVector& y) {
    int n = X.nrow();
    int p = X.ncol();
    NumericVector weights(p);
    double eps = std::numeric_limits<double>::epsilon();

    if (n < 3 || p == 0) {
        std::fill(weights.begin(), weights.end(), eps);
        return weights;
    }

    // Convert to Eigen
    Eigen::Map<Eigen::MatrixXd> X_map(as<Eigen::Map<Eigen::MatrixXd>>(X));
    Eigen::VectorXd y_vec = as<Eigen::VectorXd>(y);

    for (int j = 0; j < p; ++j) {
        double tstat = univariate_beta_tstat(y_vec, X_map.col(j));

        if (tstat > 0 && std::isfinite(tstat)) {
            weights[j] = tstat;
        } else {
            // Fall back to continuous weights on logit transform
            Eigen::VectorXd y_logit(n);
            for (int i = 0; i < n; ++i) {
                double yi = std::max(1e-10, std::min(1.0 - 1e-10, y_vec[i]));
                y_logit[i] = std::log(yi / (1.0 - yi));
            }

            // Simple OLS
            double sumy = y_logit.sum();
            double sumy2 = y_logit.squaredNorm();
            double ybar = sumy / n;
            double sumx = 0.0, sumx2 = 0.0, sumxy = 0.0;
            for (int i = 0; i < n; ++i) {
                double xi = X_map(i, j);
                sumx += xi;
                sumx2 += xi * xi;
                sumxy += xi * y_logit[i];
            }
            double xbar = sumx / n;
            double varx = sumx2 - n * xbar * xbar;

            if (varx <= eps || n <= 2) {
                weights[j] = eps;
                continue;
            }

            double b1 = (sumxy - n * xbar * ybar) / varx;
            double b0 = ybar - b1 * xbar;
            double sse = sumy2 + n * b0 * b0 + b1 * b1 * sumx2
                       + 2.0 * b0 * b1 * sumx - 2.0 * b0 * sumy - 2.0 * b1 * sumxy;

            if (n <= 2 || sse <= 0) {
                weights[j] = eps;
                continue;
            }

            double sigma2 = sse / (n - 2);
            double se_b1 = std::sqrt(sigma2 / varx);

            if (!std::isfinite(se_b1) || se_b1 <= 0) {
                weights[j] = eps;
            } else {
                weights[j] = std::fabs(b1 / se_b1);
            }
        }
    }

    return weights;
}

// [[Rcpp::export]]
NumericVector kk21_negbin_weights_cpp(const NumericMatrix& X,
                                      const NumericVector& y) {
    int n = X.nrow();
    int p = X.ncol();
    NumericVector weights(p);
    double eps = std::numeric_limits<double>::epsilon();

    if (n < 3 || p == 0) {
        std::fill(weights.begin(), weights.end(), eps);
        return weights;
    }

    // Convert to Eigen
    Eigen::Map<Eigen::MatrixXd> X_map(as<Eigen::Map<Eigen::MatrixXd>>(X));
    Eigen::VectorXd y_vec = as<Eigen::VectorXd>(y);

    // Pre-compute OLS quantities for fallback (using log(y+1))
    Eigen::VectorXd log_y1(n);
    for (int i = 0; i < n; ++i) {
        log_y1[i] = std::log1p(y_vec[i]);
    }
    double sumy = log_y1.sum();
    double sumy2 = log_y1.squaredNorm();
    double ybar = sumy / n;

    for (int j = 0; j < p; ++j) {
        double tstat = univariate_negbin_tstat(y_vec, X_map.col(j));

        if (tstat > 0 && std::isfinite(tstat)) {
            weights[j] = tstat;
            continue;
        }

        // Fall back to OLS on log(y+1)
        double sumx = 0.0, sumx2 = 0.0, sumxy = 0.0;
        for (int i = 0; i < n; ++i) {
            double x = X_map(i, j);
            sumx += x;
            sumx2 += x * x;
            sumxy += x * log_y1[i];
        }
        double xbar = sumx / n;
        double varx = sumx2 - n * xbar * xbar;

        if (varx <= eps || n <= 2) {
            weights[j] = eps;
            continue;
        }

        double b1 = (sumxy - n * xbar * ybar) / varx;
        double b0 = ybar - b1 * xbar;
        double sse = sumy2 + n * b0 * b0 + b1 * b1 * sumx2
                   + 2.0 * b0 * b1 * sumx - 2.0 * b0 * sumy - 2.0 * b1 * sumxy;

        if (n <= 2 || sse <= 0) {
            weights[j] = eps;
            continue;
        }

        double sigma2 = sse / (n - 2);
        double se_b1 = std::sqrt(sigma2 / varx);

        if (!std::isfinite(se_b1) || se_b1 <= 0) {
            weights[j] = eps;
        } else {
            weights[j] = std::fabs(b1 / se_b1);
        }
    }

    return weights;
}

// Helper: Univariate Weibull AFT regression for a single covariate
// Returns t-statistic for the covariate, or -1 if fitting fails
static double univariate_weibull_tstat(
    const Eigen::VectorXd& log_y,
    const Eigen::VectorXd& delta,
    const Eigen::VectorXd& x,
    int maxit = 30,
    double tol = 1e-5
) {
    int n = log_y.size();
    if (n < 4) return -1.0;

    double n_events = delta.sum();
    if (n_events < 2) return -1.0;

    // Build design matrix [1, x]
    Eigen::MatrixXd X(n, 2);
    X.col(0).setOnes();
    X.col(1) = x;

    // Initialize with OLS on log(y)
    Eigen::MatrixXd XtX = X.transpose() * X;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(XtX);
    if (qr.rank() < 2) return -1.0;

    Eigen::VectorXd beta = qr.solve(X.transpose() * log_y);
    Eigen::VectorXd resid = log_y - X * beta;

    // Initial scale estimate
    double scale = 0.0;
    for (int i = 0; i < n; ++i) {
        if (delta(i) > 0.5) {
            scale += resid(i) * resid(i);
        }
    }
    scale = std::sqrt(scale / std::max(1.0, n_events - 2));
    scale = std::max(0.01, std::min(10.0, scale));

    // Newton-Raphson for Weibull AFT
    for (int iter = 0; iter < maxit; ++iter) {
        Eigen::VectorXd z = resid / scale;
        Eigen::VectorXd w(n);
        Eigen::VectorXd adj(n);

        for (int i = 0; i < n; ++i) {
            double exp_z = std::exp(z(i));
            if (delta(i) > 0.5) {
                w(i) = exp_z;
                adj(i) = z(i) - 1.0 + exp_z;
            } else {
                w(i) = exp_z;
                adj(i) = exp_z;
            }
        }

        // Weighted least squares update
        Eigen::VectorXd sqrt_w = w.array().sqrt();
        Eigen::MatrixXd Xw = X.array().colwise() * sqrt_w.array();
        Eigen::VectorXd yw = (log_y.array() - scale * adj.array() / w.array()) * sqrt_w.array();

        Eigen::MatrixXd XwXw = Xw.transpose() * Xw;
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr_w(XwXw);
        if (qr_w.rank() < 2) return -1.0;

        Eigen::VectorXd beta_new = qr_w.solve(Xw.transpose() * yw);
        Eigen::VectorXd resid_new = log_y - X * beta_new;

        // Scale update
        double sum_d = 0.0;
        double scale_new = scale;
        for (int i = 0; i < n; ++i) {
            if (delta(i) > 0.5) sum_d += 1.0;
        }
        if (sum_d > 0) {
            double score = 0.0, info = 0.0;
            for (int i = 0; i < n; ++i) {
                double z_i = resid_new(i) / scale;
                double exp_z_i = std::exp(z_i);
                if (delta(i) > 0.5) {
                    score += -1.0/scale + z_i/scale - z_i*exp_z_i/scale;
                    info += 1.0/(scale*scale);
                } else {
                    score += -z_i*exp_z_i/scale;
                }
            }
            if (std::fabs(info) > 1e-10) {
                scale_new = scale - score / info;
                scale_new = std::max(0.01, std::min(10.0, scale_new));
            }
        }

        double diff = (beta_new - beta).norm() + std::fabs(scale_new - scale);
        beta = beta_new;
        scale = scale_new;
        resid = resid_new;

        if (diff < tol) break;
    }

    // Compute standard error of beta[1]
    Eigen::VectorXd z = resid / scale;
    Eigen::VectorXd w(n);
    for (int i = 0; i < n; ++i) {
        w(i) = std::exp(z(i)) / (scale * scale);
    }

    Eigen::MatrixXd Xw = X.array().colwise() * w.array().sqrt();
    Eigen::MatrixXd info_mat = Xw.transpose() * Xw;

    Eigen::FullPivLU<Eigen::MatrixXd> lu(info_mat);
    if (!lu.isInvertible()) return -1.0;

    Eigen::MatrixXd cov_mat = lu.inverse();
    if (cov_mat(1, 1) <= 0) return -1.0;

    double se = std::sqrt(cov_mat(1, 1));
    if (!std::isfinite(se) || se <= 0) return -1.0;

    return std::fabs(beta(1) / se);
}

// [[Rcpp::export]]
NumericVector kk21_survival_weights_cpp(const NumericMatrix& X,
                                        const NumericVector& y,
                                        const NumericVector& delta) {
    int n = X.nrow();
    int p = X.ncol();
    NumericVector weights(p);
    double eps = std::numeric_limits<double>::epsilon();

    if (n < 2 || p == 0) {
        std::fill(weights.begin(), weights.end(), eps);
        return weights;
    }

    // Convert to Eigen
    Eigen::Map<Eigen::MatrixXd> X_map(as<Eigen::Map<Eigen::MatrixXd>>(X));
    Eigen::VectorXd y_vec = as<Eigen::VectorXd>(y);
    Eigen::VectorXd delta_vec = as<Eigen::VectorXd>(delta);

    // Compute log(y)
    Eigen::VectorXd log_y(n);
    for (int i = 0; i < n; ++i) {
        log_y(i) = std::log(std::max(y_vec(i), 1e-10));
    }

    // Pre-compute OLS quantities for fallback
    double sumy = log_y.sum();
    double sumy2 = log_y.squaredNorm();
    double ybar = sumy / n;

    for (int j = 0; j < p; ++j) {
        // Try Weibull first
        double tstat = univariate_weibull_tstat(log_y, delta_vec, X_map.col(j));

        if (tstat > 0 && std::isfinite(tstat)) {
            weights[j] = tstat;
            continue;
        }

        // Fall back to OLS on log(y)
        double sumx = 0.0, sumx2 = 0.0, sumxy = 0.0;
        for (int i = 0; i < n; ++i) {
            double x = X_map(i, j);
            sumx += x;
            sumx2 += x * x;
            sumxy += x * log_y(i);
        }
        double xbar = sumx / n;
        double varx = sumx2 - n * xbar * xbar;

        if (varx <= eps || n <= 2) {
            weights[j] = eps;
            continue;
        }

        double b1 = (sumxy - n * xbar * ybar) / varx;
        double b0 = ybar - b1 * xbar;
        double sse = sumy2 + n * b0 * b0 + b1 * b1 * sumx2 +
                     2.0 * b0 * b1 * sumx - 2.0 * b0 * sumy - 2.0 * b1 * sumxy;

        if (n <= 2 || sse <= 0) {
            weights[j] = eps;
            continue;
        }

        double sigma2 = sse / (n - 2);
        double se_b1 = std::sqrt(sigma2 / varx);

        if (!std::isfinite(se_b1) || se_b1 <= 0) {
            weights[j] = eps;
        } else {
            weights[j] = std::fabs(b1 / se_b1);
        }
    }

    return weights;
}

// Helper: Multivariate beta regression returning t-statistic for coefficient at index coef_idx
// Design matrix X should include intercept as first column
static double multivariate_beta_tstat(
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& y,
    int coef_idx = 1,
    int maxit = 100,
    double tol = 1e-8
) {
    int n = X.rows();
    int p = X.cols();
    if (n < p + 1 || p < 2) return -1.0;

    // Fit logistic regression via IRLS (for mean model with logit link)
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    Eigen::VectorXd prob(n);
    Eigen::VectorXd w(n);

    for (int iter = 0; iter < maxit; ++iter) {
        Eigen::VectorXd eta = X * beta;
        for (int i = 0; i < n; ++i) {
            double val = 1.0 / (1.0 + std::exp(-eta[i]));
            val = std::max(1e-10, std::min(1.0 - 1e-10, val));
            prob[i] = val;
            w[i] = std::max(val * (1.0 - val), 1e-10);
        }

        Eigen::VectorXd z = eta + (y - prob).cwiseQuotient(w);
        Eigen::MatrixXd XtW = X.transpose() * w.asDiagonal();
        Eigen::MatrixXd XtWX = XtW * X;
        Eigen::VectorXd XtWz = XtW * z;

        Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
        if (ldlt.info() != Eigen::Success) return -1.0;

        Eigen::VectorXd beta_new = ldlt.solve(XtWz);
        if ((beta - beta_new).norm() < tol) {
            beta = beta_new;
            break;
        }
        beta = beta_new;
    }

    // Compute final mu and weights
    Eigen::VectorXd eta = X * beta;
    for (int i = 0; i < n; ++i) {
        double val = 1.0 / (1.0 + std::exp(-eta[i]));
        val = std::max(1e-12, std::min(1.0 - 1e-12, val));
        prob[i] = val;
        w[i] = std::max(val * (1.0 - val), 1e-10);
    }

    // Profile likelihood optimization for phi (grid search)
    double best_phi = 10.0;
    double best_ll = -std::numeric_limits<double>::infinity();

    for (double log_phi = -2.0; log_phi <= 7.0; log_phi += 0.5) {
        double phi = std::exp(log_phi);
        double ll = 0.0;
        for (int i = 0; i < n; ++i) {
            double mu_i = prob[i];
            double yi = y[i];
            double mu_phi = mu_i * phi;
            double one_minus_mu_phi = (1.0 - mu_i) * phi;
            mu_phi = std::max(1e-12, mu_phi);
            one_minus_mu_phi = std::max(1e-12, one_minus_mu_phi);

            ll += R::lgammafn(phi) - R::lgammafn(mu_phi) - R::lgammafn(one_minus_mu_phi)
                + (mu_phi - 1.0) * std::log(std::max(yi, 1e-12))
                + (one_minus_mu_phi - 1.0) * std::log1p(-std::min(yi, 1.0 - 1e-12));
        }
        if (ll > best_ll) {
            best_ll = ll;
            best_phi = phi;
        }
    }

    // Final weights for variance computation
    Eigen::VectorXd w_final(n);
    for (int i = 0; i < n; ++i) {
        w_final[i] = (1.0 + best_phi) * w[i];
    }

    // Compute XtWX and get variance of beta[coef_idx]
    Eigen::MatrixXd XtWX = X.transpose() * w_final.asDiagonal() * X;
    Eigen::FullPivLU<Eigen::MatrixXd> lu(XtWX);
    if (!lu.isInvertible()) return -1.0;

    Eigen::MatrixXd cov_mat = lu.inverse();
    if (coef_idx >= p || cov_mat(coef_idx, coef_idx) <= 0) return -1.0;

    double se = std::sqrt(cov_mat(coef_idx, coef_idx));
    if (!std::isfinite(se) || se <= 0) return -1.0;

    return std::fabs(beta(coef_idx) / se);
}

// Helper: Multivariate negative binomial regression returning t-statistic for coefficient at index coef_idx
// Design matrix X should include intercept as first column
static double multivariate_negbin_tstat(
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& y,
    int coef_idx = 1,
    int maxit = 50,
    double tol = 1e-6
) {
    int n = X.rows();
    int p = X.cols();
    if (n < p + 1 || p < 2) return -1.0;

    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);

    // Initialize theta using method of moments
    double y_mean = y.mean();
    double y_var = (y.array() - y_mean).square().sum() / (n - 1);
    double theta = (y_var > y_mean) ? (y_mean * y_mean) / (y_var - y_mean) : 10.0;
    theta = std::max(0.01, std::min(1000.0, theta));

    Eigen::VectorXd mu(n);
    Eigen::VectorXd w(n);

    for (int outer = 0; outer < maxit; ++outer) {
        // IRLS for beta given theta
        for (int iter = 0; iter < 25; ++iter) {
            Eigen::VectorXd eta = X * beta;
            for (int i = 0; i < n; ++i) {
                mu[i] = std::exp(std::min(eta[i], 20.0));
                w[i] = mu[i] / (1.0 + mu[i] / theta);
                w[i] = std::max(w[i], 1e-10);
            }

            Eigen::VectorXd z = eta + (y - mu).cwiseQuotient(mu);
            Eigen::MatrixXd XtW = X.transpose() * w.asDiagonal();
            Eigen::MatrixXd XtWX = XtW * X;
            Eigen::VectorXd XtWz = XtW * z;

            Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
            if (ldlt.info() != Eigen::Success) return -1.0;

            Eigen::VectorXd beta_new = ldlt.solve(XtWz);
            if ((beta - beta_new).norm() < tol) {
                beta = beta_new;
                break;
            }
            beta = beta_new;
        }

        // Update mu
        Eigen::VectorXd eta = X * beta;
        for (int i = 0; i < n; ++i) {
            mu[i] = std::exp(std::min(eta[i], 20.0));
        }

        // Update theta using Newton-Raphson
        double score = 0.0;
        double info = 0.0;
        for (int i = 0; i < n; ++i) {
            double yi = y[i];
            double mui = mu[i];
            score += R::digamma(yi + theta) - R::digamma(theta)
                   + std::log(theta) - std::log(theta + mui) + 1.0
                   - (yi + theta) / (theta + mui);
            info += -R::trigamma(yi + theta) + R::trigamma(theta)
                  - 1.0 / theta + 2.0 / (theta + mui)
                  - (yi + theta) / ((theta + mui) * (theta + mui));
        }

        if (std::fabs(info) < 1e-10) break;
        double theta_new = theta - score / info;
        theta_new = std::max(0.01, std::min(1000.0, theta_new));

        if (std::fabs(theta_new - theta) < tol) {
            theta = theta_new;
            break;
        }
        theta = theta_new;
    }

    // Compute final weights
    Eigen::VectorXd eta = X * beta;
    for (int i = 0; i < n; ++i) {
        mu[i] = std::exp(std::min(eta[i], 20.0));
        w[i] = mu[i] / (1.0 + mu[i] / theta);
        w[i] = std::max(w[i], 1e-10);
    }

    Eigen::MatrixXd XtWX = X.transpose() * w.asDiagonal() * X;
    Eigen::FullPivLU<Eigen::MatrixXd> lu(XtWX);
    if (!lu.isInvertible()) return -1.0;

    Eigen::MatrixXd cov_mat = lu.inverse();
    if (coef_idx >= p || cov_mat(coef_idx, coef_idx) <= 0) return -1.0;

    double se = std::sqrt(cov_mat(coef_idx, coef_idx));
    if (!std::isfinite(se) || se <= 0) return -1.0;

    return std::fabs(beta(coef_idx) / se);
}

// [[Rcpp::export]]
NumericVector kk21_stepwise_beta_weights_cpp(const NumericMatrix& X,
                                             const NumericVector& y,
                                             const NumericVector& w) {
    int n = X.nrow();
    int p = X.ncol();
    NumericVector weights(p, NA_REAL);

    if (n == 0 || p == 0) {
        return weights;
    }

    Eigen::Map<Eigen::MatrixXd> X_map(as<Eigen::Map<Eigen::MatrixXd>>(X));
    Eigen::VectorXd y_vec = as<Eigen::VectorXd>(y);
    Eigen::VectorXd w_vec = as<Eigen::VectorXd>(w);

    // Clamp y to (0,1)
    for (int i = 0; i < n; ++i) {
        y_vec[i] = std::max(1e-10, std::min(1.0 - 1e-10, y_vec[i]));
    }

    std::vector<int> selected;
    selected.reserve(p);
    std::vector<bool> used(p, false);

    for (int step = 0; step < p; ++step) {
        double best_stat = -1.0;
        int best_j = -1;
        int k = static_cast<int>(selected.size());

        for (int j = 0; j < p; ++j) {
            if (used[j]) {
                continue;
            }

            // Build design matrix: [1, x_j, selected_covs..., w]
            Eigen::MatrixXd Xmm(n, k + 3);
            Xmm.col(0).setOnes();
            Xmm.col(1) = X_map.col(j);
            for (int idx = 0; idx < k; ++idx) {
                Xmm.col(2 + idx) = X_map.col(selected[idx]);
            }
            Xmm.col(k + 2) = w_vec;

            double stat = multivariate_beta_tstat(Xmm, y_vec, 1);
            if (!R_finite(stat) || stat < 0) {
                stat = 0.0;
            }
            if (stat > best_stat) {
                best_stat = stat;
                best_j = j;
            }
        }

        if (best_j < 0) {
            break;
        }

        weights[best_j] = best_stat;
        used[best_j] = true;
        selected.push_back(best_j);
    }

    return weights;
}

// [[Rcpp::export]]
NumericVector kk21_stepwise_negbin_weights_cpp(const NumericMatrix& X,
                                               const NumericVector& y,
                                               const NumericVector& w) {
    int n = X.nrow();
    int p = X.ncol();
    NumericVector weights(p, NA_REAL);

    if (n == 0 || p == 0) {
        return weights;
    }

    Eigen::Map<Eigen::MatrixXd> X_map(as<Eigen::Map<Eigen::MatrixXd>>(X));
    Eigen::VectorXd y_vec = as<Eigen::VectorXd>(y);
    Eigen::VectorXd w_vec = as<Eigen::VectorXd>(w);

    std::vector<int> selected;
    selected.reserve(p);
    std::vector<bool> used(p, false);

    for (int step = 0; step < p; ++step) {
        double best_stat = -1.0;
        int best_j = -1;
        int k = static_cast<int>(selected.size());

        for (int j = 0; j < p; ++j) {
            if (used[j]) {
                continue;
            }

            // Build design matrix: [1, x_j, selected_covs..., w]
            Eigen::MatrixXd Xmm(n, k + 3);
            Xmm.col(0).setOnes();
            Xmm.col(1) = X_map.col(j);
            for (int idx = 0; idx < k; ++idx) {
                Xmm.col(2 + idx) = X_map.col(selected[idx]);
            }
            Xmm.col(k + 2) = w_vec;

            double stat = multivariate_negbin_tstat(Xmm, y_vec, 1);
            if (!R_finite(stat) || stat < 0) {
                stat = 0.0;
            }
            if (stat > best_stat) {
                best_stat = stat;
                best_j = j;
            }
        }

        if (best_j < 0) {
            break;
        }

        weights[best_j] = best_stat;
        used[best_j] = true;
        selected.push_back(best_j);
    }

    return weights;
}
