#include <RcppEigen.h>
#include <cmath>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

//' Sequential computation of both CI bounds (called from R-level parallelism)
//'
//' This function is kept for backwards compatibility but the outer parallelism
//' (running lower/upper bounds simultaneously) is now handled at the R level via
//' parallel::mclapply, which is safe. Calling R functions from OpenMP threads is
//' undefined behaviour in R and caused process crashes.
//'
//' @param pval_fn R function: pval_fn(nsim, delta, transform_responses, num_cores)
//' @param nsim_exact_test Number of randomization iterations
//' @param l_lower Initial lower bound for lower CI bound search
//' @param u_lower Initial upper bound for lower CI bound search (typically the estimate)
//' @param l_upper Initial lower bound for upper CI bound search (typically the estimate)
//' @param u_upper Initial upper bound for upper CI bound search
//' @param pval_th P-value threshold (typically alpha/2 for two-sided CI)
//' @param tol Tolerance for convergence (in p-value space)
//' @param transform_responses String: "none", "log", or "logit"
//' @param num_cores Passed through to pval_fn for inner parallelism
//'
//' @return Numeric vector of length 2: [lower_bound, upper_bound]
//'
//' @export
// [[Rcpp::export]]
NumericVector bisection_ci_parallel_cpp(
  Function pval_fn,
  int nsim_exact_test,
  double l_lower,
  double u_lower,
  double l_upper,
  double u_upper,
  double pval_th,
  double tol,
  std::string transform_responses,
  int num_cores = 1
) {
  NumericVector ci_bounds(2);

  // --- Lower CI bound ---
  {
    double l = l_lower;
    double u = u_lower;

    double pval_l = as<double>(pval_fn(nsim_exact_test, l, transform_responses, num_cores));
    double pval_u = as<double>(pval_fn(nsim_exact_test, u, transform_responses, num_cores));

    while (R_finite(pval_l) && R_finite(pval_u) && (pval_u - pval_l) > tol) {
      double m = (l + u) / 2.0;
      double pval_m = as<double>(pval_fn(nsim_exact_test, m, transform_responses, num_cores));

      if (!R_finite(pval_m)) { l = m; pval_l = 0.0; }
      else if (pval_m >= pval_th) { u = m; pval_u = pval_m; }
      else                        { l = m; pval_l = pval_m; }
    }
    ci_bounds[0] = l;
  }

  // --- Upper CI bound ---
  {
    double l = l_upper;
    double u = u_upper;

    double pval_l = as<double>(pval_fn(nsim_exact_test, l, transform_responses, num_cores));
    double pval_u = as<double>(pval_fn(nsim_exact_test, u, transform_responses, num_cores));

    while (R_finite(pval_l) && R_finite(pval_u) && (pval_u - pval_l) > tol) {
      double m = (l + u) / 2.0;
      double pval_m = as<double>(pval_fn(nsim_exact_test, m, transform_responses, num_cores));

      if (!R_finite(pval_m)) { u = m; pval_u = 0.0; }
      else if (pval_m >= pval_th) { l = m; pval_l = pval_m; }
      else                        { u = m; pval_u = pval_m; }
    }
    ci_bounds[1] = u;
  }

  return ci_bounds;
}


//' Single-threaded helper for computing one CI bound
//'
//' @keywords internal
// [[Rcpp::export]]
double bisection_ci_single_bound_cpp(
  Function pval_fn,
  int nsim_exact_test,
  double l,
  double u,
  double pval_th,
  double tol,
  std::string transform_responses,
  bool lower,
  int num_cores = 1
) {
  double pval_l = as<double>(pval_fn(nsim_exact_test, l, transform_responses, num_cores));
  double pval_u = as<double>(pval_fn(nsim_exact_test, u, transform_responses, num_cores));

  while (R_finite(pval_l) && R_finite(pval_u) && (pval_u - pval_l) > tol) {
    double m = (l + u) / 2.0;
    double pval_m = as<double>(pval_fn(nsim_exact_test, m, transform_responses, num_cores));

    if (!R_finite(pval_m)) {
      if (lower) { l = m; pval_l = 0.0; }
      else       { u = m; pval_u = 0.0; }
    } else if (pval_m >= pval_th && lower)  { u = m; pval_u = pval_m; }
    else if   (pval_m >= pval_th && !lower) { l = m; pval_l = pval_m; }
    else if   (lower)                       { l = m; pval_l = pval_m; }
    else                                    { u = m; pval_u = pval_m; }
  }

  return lower ? l : u;
}
