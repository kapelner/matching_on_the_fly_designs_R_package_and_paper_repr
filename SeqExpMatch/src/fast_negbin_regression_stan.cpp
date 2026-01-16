//// [[Rcpp::depends(RcppEigen)]]
//// [[Rcpp::depends(StanHeaders)]]
//#include "_helper_functions.h"
//#include <Rcpp.h>
//#include <RcppEigen.h>
//#include <stan/math.hpp>
//#include "optimization/LBFGS.h"  // Make sure this points to the correct path of LBFGSpp
//
//using namespace Rcpp;
//using namespace stan::math;
//using Eigen::VectorXd;
//using Eigen::MatrixXd;
//
//// Define functor for LBFGS++
//class NegBinCensoredFunctor {
//  private:
//	  const VectorXd& y;
//	  const VectorXd& dead;
//	  const MatrixXd& X;
//  public:
//	  NegBinCensoredFunctor(const VectorXd& y_, const VectorXd& dead_, const MatrixXd& X_) :
//	    y(y_), dead(dead_), X(X_) {}
//	
//	  double operator()(const VectorXd& x, VectorXd& grad) {
////	    using stan::math::var;
////	    using stan::math::gradient;
////	
////	    std::vector<var> theta(x.size());
////	    for (int i = 0; i < x.size(); ++i) theta[i] = x[i];
////	
////	    int K = X.cols();
////	    VectorXd b = x.head(K);
////	    var phi = x[K];
////	
////	    auto log_lik_var = [&]() {
////	      var ll = 0;
////	      for (int i = 0; i < y.size(); ++i) {
////	        var eta = dot_product(X.row(i), b);
////	        var mu = exp(eta);
////	        var alpha = phi;
////	
////	        if (dead[i] == 1) {
////	          ll += stan::math::neg_binomial_log(y[i], mu, alpha);
////	        } else {
////	          ll += log(1.0 - stan::math::neg_binomial_cdf(y[i], mu, alpha));
////	        }
////	      }
////	      return ll;
////	    };
////	
////	    var ll = log_lik_var();
////	    double val = ll.val();
////	
////	    stan::math::grad(ll.vi_);
////	    for (int i = 0; i < theta.size(); ++i) {
////	      grad[i] = theta[i].adj();
////	    }
////	
////	    return -val; // LBFGS minimizes, so return negative log-lik
//return 0;
//	  }
//};
//
//class Rosenbrock
//{
//private:
//    int n;
//public:
//    Rosenbrock(int n_) : n(n_) {}
//    double operator()(const VectorXd& x, VectorXd& grad)
//    {
//        double fx = 0.0;
//        for(int i = 0; i < n; i += 2)
//        {
//            double t1 = 1.0 - x[i];
//            double t2 = 10 * (x[i + 1] - x[i] * x[i]);
//            grad[i + 1] = 20 * t2;
//            grad[i]     = -2.0 * (x[i] * grad[i + 1] + t1);
//            fx += t1 * t1 + t2 * t2;
//        }
//        return fx;
//    }
//};
//
//// [[Rcpp::export]]
//List fast_stan_neg_bin_with_censoring_cpp(const MatrixXd& X,
//                         const VectorXd& y,
//                         const VectorXd& dead,
//                         const VectorXd& b_init,
//                         double tol = 1e-6,
//                         int maxit = 1000) {
////  NegBinCensoredFunctor negbin_censored_log_lik_function(y, dead, X);
//  int n = 10;
//  Rosenbrock ros(n);
//
//  LBFGSpp::LBFGSParam<double> param;
//  param.epsilon = tol;
//  param.max_iterations = maxit;
//  LBFGSpp::LBFGSSolver<double> solver(param);
//
//  double ll;
////  int niter = solver.minimize(negbin_censored_log_lik_function, b_init, ll);
//  VectorXd x = VectorXd::Zero(n);
//  int niter = solver.minimize(ros, x, ll);
//
//  return List::create(
//    _["b"] = b_init,
//    _["logLik"] = -ll,
//    _["niter"] = niter
//  );
//}
//
//// [[Rcpp::export]]
//List fast_stan_neg_bin_with_censoring_with_sd_cpp(const MatrixXd& X,
//                         const VectorXd& y,
//                         const VectorXd& dead,
//                         const VectorXd& b_init,
//                         double tol = 1e-6,
//                         int maxit = 1000) {
//
//  List mod = fast_stan_neg_bin_with_censoring_cpp(X, y, dead, b_init, tol, maxit);
//
//  return mod;
//}
