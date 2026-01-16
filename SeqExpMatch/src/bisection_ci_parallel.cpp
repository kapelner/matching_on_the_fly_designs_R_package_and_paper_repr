#include <RcppEigen.h>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

//' Parallel computation of both CI bounds by inverting randomization tests
//'
//' This function computes both the lower and upper confidence interval bounds
//' in parallel using OpenMP. Each bound uses the bisection algorithm.
//'
//' @param pval_fn R function that computes two-sided p-value given delta
//' @param nsim_exact_test Number of randomization iterations
//' @param l_lower Initial lower bound for lower CI bound search
//' @param u_lower Initial upper bound for lower CI bound search (typically the estimate)
//' @param l_upper Initial lower bound for upper CI bound search (typically the estimate)
//' @param u_upper Initial upper bound for upper CI bound search
//' @param pval_th P-value threshold (typically alpha/2 for two-sided CI)
//' @param tol Tolerance for convergence (in p-value space)
//' @param transform_responses String: "none", "log", or "logit"
//' @param num_cores Number of cores to use (0 = auto-detect, 1 = no parallelization)
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
  int num_cores = 0
) {
  NumericVector ci_bounds(2);

  // Determine number of threads for outer parallelization (CI bounds)
  int outer_threads = 1;
#ifdef _OPENMP
  if (num_cores == 0) {
    outer_threads = std::min(2, omp_get_max_threads());
  } else if (num_cores > 1) {
    outer_threads = std::min(2, num_cores);
  }

  // Allocate remaining threads for inner parallelization (randomization)
  // Each CI bound search will use num_cores_inner threads for randomization
  int num_cores_inner = (num_cores > 2) ? (num_cores / 2) : 1;
#else
  int num_cores_inner = 1;
#endif

  // Compute both bounds in parallel
#pragma omp parallel num_threads(outer_threads) if(outer_threads > 1)
  {
#pragma omp sections
    {
      // Section 1: Compute lower CI bound
#pragma omp section
      {
        double l = l_lower;
        double u = u_lower;

        // Compute initial p-values at bounds
        SEXP pval_l_sexp = pval_fn(nsim_exact_test, l, transform_responses, num_cores_inner);
        SEXP pval_u_sexp = pval_fn(nsim_exact_test, u, transform_responses, num_cores_inner);
        double pval_l = as<double>(pval_l_sexp);
        double pval_u = as<double>(pval_u_sexp);

        // Bisection loop for lower bound
        while ((pval_u - pval_l) > tol) {
          double m = (l + u) / 2.0;
          SEXP pval_m_sexp = pval_fn(nsim_exact_test, m, transform_responses, num_cores_inner);
          double pval_m = as<double>(pval_m_sexp);

          if (pval_m >= pval_th) {
            u = m;
            pval_u = pval_m;
          } else {
            l = m;
            pval_l = pval_m;
          }
        }

        ci_bounds[0] = l;
      }

      // Section 2: Compute upper CI bound
#pragma omp section
      {
        double l = l_upper;
        double u = u_upper;

        // Compute initial p-values at bounds
        SEXP pval_l_sexp = pval_fn(nsim_exact_test, l, transform_responses, num_cores_inner);
        SEXP pval_u_sexp = pval_fn(nsim_exact_test, u, transform_responses, num_cores_inner);
        double pval_l = as<double>(pval_l_sexp);
        double pval_u = as<double>(pval_u_sexp);

        // Bisection loop for upper bound
        while ((pval_u - pval_l) > tol) {
          double m = (l + u) / 2.0;
          SEXP pval_m_sexp = pval_fn(nsim_exact_test, m, transform_responses, num_cores_inner);
          double pval_m = as<double>(pval_m_sexp);

          if (pval_m >= pval_th) {
            l = m;
            pval_l = pval_m;
          } else {
            u = m;
            pval_u = pval_m;
          }
        }

        ci_bounds[1] = u;
      }
    }
  }

  return ci_bounds;
}


//' Single-threaded helper for computing one CI bound
//'
//' Helper function that computes a single CI bound. Used internally.
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
  // Compute initial p-values at bounds
  SEXP pval_l_sexp = pval_fn(nsim_exact_test, l, transform_responses, num_cores);
  SEXP pval_u_sexp = pval_fn(nsim_exact_test, u, transform_responses, num_cores);
  double pval_l = as<double>(pval_l_sexp);
  double pval_u = as<double>(pval_u_sexp);

  // Bisection loop
  while ((pval_u - pval_l) > tol) {
    // Compute midpoint
    double m = (l + u) / 2.0;
    SEXP pval_m_sexp = pval_fn(nsim_exact_test, m, transform_responses, num_cores);
    double pval_m = as<double>(pval_m_sexp);

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

  return lower ? l : u;
}
