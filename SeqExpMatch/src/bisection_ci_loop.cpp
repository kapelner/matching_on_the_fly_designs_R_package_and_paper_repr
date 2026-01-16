#include <RcppEigen.h>
#include <cmath>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

//' Bisection loop for computing confidence interval bounds by inverting randomization tests
//'
//' This function implements the bisection algorithm to find CI bounds by inverting
//' the randomization test. It repeatedly calls the p-value computation function
//' until convergence.
//'
//' @param pval_fn R function that computes two-sided p-value given delta
//' @param nsim_exact_test Number of randomization iterations
//' @param l Initial lower bound
//' @param u Initial upper bound
//' @param pval_th P-value threshold (typically alpha/2 for two-sided CI)
//' @param tol Tolerance for convergence (in p-value space)
//' @param transform_responses String: "none", "log", or "logit"
//' @param lower Logical: TRUE for lower CI bound, FALSE for upper
//'
//' @return The CI bound value
//'
//' @export
// [[Rcpp::export]]
double bisection_ci_loop_cpp(
  Function pval_fn,
  int nsim_exact_test,
  double l,
  double u,
  double pval_th,
  double tol,
  std::string transform_responses,
  bool lower
) {
  // Compute initial p-values at bounds
  double pval_l = as<double>(pval_fn(nsim_exact_test, l, transform_responses));
  double pval_u = as<double>(pval_fn(nsim_exact_test, u, transform_responses));

  // Bisection loop
  while (true) {
    // Check convergence
    if ((pval_u - pval_l) <= tol) {
      return lower ? l : u;
    }

    // Compute midpoint
    double m = (l + u) / 2.0;
    double pval_m = as<double>(pval_fn(nsim_exact_test, m, transform_responses));

    // Update bounds based on bisection logic
    if (pval_m >= pval_th && lower) {
      u = m;
      pval_u = pval_m;
    } else if (pval_m >= pval_th && !lower) {
      l = m;
      pval_l = pval_m;
    } else if (lower) {
      l = m;
      pval_l = pval_m;
    } else { // !lower
      u = m;
      pval_u = pval_m;
    }
  }

  // Should never reach here due to while(true), but compiler needs return
  return lower ? l : u;
}
