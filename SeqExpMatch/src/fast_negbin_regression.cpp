//#include "_helper_functions.h"
//#include <RcppNumerical.h>
//
//using namespace Numer;
//using namespace Eigen;
//using namespace Rcpp;
//
////requires tbb
////https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb-download.html
//
//// Negative binomial censored log-likelihood with gradient
//class NBCensoredLogLik : public MFuncGrad {
//	private:
//	  const MatrixXd& X;
//	  const VectorXi& y;
//	  const VectorXi& dead;
//	  int p;
//	
//	public:
//	  NBCensoredLogLik(const MatrixXd& X_,
//	                   const VectorXi& y_,
//	                   const VectorXi& dead_)
//	    : X(X_), y(y_), dead(dead_), p(X_.cols()) {}
//
//	  double f_grad(Constvec& pars, Refvec grad) {
//	    const VectorXd beta = pars.head(p);
//	    const double log_theta = pars[p];
//	    const double theta = std::exp(log_theta);
//	
//	    VectorXd eta = X * beta;
//	    VectorXd mu = eta.array().exp();
//	
//	    double nll = 0.0;
//	    VectorXd score_beta = VectorXd::Zero(p);
//	    double score_log_theta = 0.0;
//	
//	    for (int i = 0; i < X.rows(); ++i) {
//	      double mu_i = mu[i];
//	      double r = theta;
//	      double yi = y[i];
//	
//	      if (dead[i] == 1) {
//	        // Log-likelihood contribution
//	        nll -= R::dnbinom_mu(yi, r, mu_i, true);
//	
//	        // Gradient w.r.t. beta
//	        double coef = yi - mu_i * (yi + r) / (r + mu_i);
//	        score_beta += coef * X.row(i).transpose();
//	
//	        // Gradient w.r.t. log(theta)
//	        double dlogf_dr = 
//	          R::digamma(yi + r) - R::digamma(r) +
//	          std::log(r) - std::log(r + mu_i) +
//	          1.0 - (yi + r) / (r + mu_i);
//	        score_log_theta += r * dlogf_dr;
//	
//	      } else {
//	        // Right-censored log-likelihood
//	        double surv = 1.0 - R::pnbinom_mu(yi, r, mu_i, true, false);
//	        if (surv > 1e-12)
//	          nll -= std::log(surv);
//	        else
//	          nll += 1e6;  // Penalty for numerical underflow
//	        // Gradient is omitted for censored (this is for speed - it can be added in later for accuracy when there is many censored observations)
//	      }
//	    }
//	
//	    grad.head(p) = -score_beta;
//	    grad[p] = -score_log_theta;
//	
//	    return nll;
//	  }
//};
//
//
//// Negative Binomial log CCDF for right-censoring
//double log_ccdf_nb(int y, double mu, double phi) {
//  return std::log(1.0 - R::pnbinom_mu(y, phi, mu, 1, 0)); // lower.tail=TRUE, log=FALSE
//}
//
//// Analytic Hessian from uncensored observations
//MatrixXd analytic_hessian_uncensored(const MatrixXd &X,
//									const VectorXi &y,
//                                   	const VectorXi &dead,
//                                  	const VectorXd &b,
//                                  	const double theta_hat,
//                                  	double eps = 1e-5) {
//  int p_plus_one = X.cols();
//  MatrixXd H = MatrixXd::Zero(p_plus_one + 1, p_plus_one + 1);
//
//  for (int i = 0; i < y.size(); ++i) {
//    if (dead[i] == 0) continue;
//
//    double yi = (double)y[i];
//    double eta = X.row(i).dot(b);
//    double mu = std::exp(eta);
//    double denom = mu + theta_hat;
//
//    VectorXd x = X.row(i);
//    VectorXd dmu_dbeta = mu * x;
//    MatrixXd d2mu_dbetasq = mu * (x * x.transpose());
//
//    // ββ
//    H.topLeftCorner(p_plus_one, p_plus_one) +=
//      (yi / (mu * denom) - (yi + theta_hat) / (denom * denom)) * (dmu_dbeta * dmu_dbeta.transpose()) +
//      ((yi + theta_hat) / denom) * d2mu_dbetasq;
//
//    // βφ
//    H.topRightCorner(p_plus_one, 1) +=
//      (1.0 / denom - (yi + theta_hat) / (denom * denom)) * dmu_dbeta;
//
//    // φβ
//    H.bottomLeftCorner(1, p_plus_one) = H.topRightCorner(p_plus_one, 1).transpose();
//
//    // φφ
//    H(p_plus_one, p_plus_one) += R::trigamma(yi + theta_hat) - R::trigamma(theta_hat)
//               + std::pow(1.0 / denom - (yi + theta_hat) / (denom * denom), 2)
//               + (yi + theta_hat) / (denom * denom);
//  }
//
//  return H;
//}
//
//// Numeric Hessian from censored observations via finite differences
//MatrixXd numeric_hessian_censored(const MatrixXd &X,
//									const VectorXi &y,
//                                   	const VectorXi &dead,
//                                  	const VectorXd &b,
//                                  	const double theta_hat,
//                                  	double eps = 1e-5) {
//  int p_plus_two = b.size() + 1;
//  VectorXd theta(p_plus_two);
//  for (int i = 0; i < b.size(); ++i) {
//	theta[i] = b[i];
//  }
//  theta[p_plus_two] = theta_hat;
//  MatrixXd H = MatrixXd::Zero(p_plus_two, p_plus_two);
//
//  for (int i = 0; i < p_plus_two; ++i) {
//    for (int j = i; j < p_plus_two; ++j) {
//      VectorXd theta_pp = theta, theta_pm = theta, theta_mp = theta, theta_mm = theta;
//      theta_pp(i) += eps; theta_pp(j) += eps;
//      theta_pm(i) += eps; theta_pm(j) -= eps;
//      theta_mp(i) -= eps; theta_mp(j) += eps;
//      theta_mm(i) -= eps; theta_mm(j) -= eps;
//
//      double fpp = 0, fpm = 0, fmp = 0, fmm = 0;
//
//      for (int i = 0; i < y.size(); ++i) {
//        if (dead[i] == 1) continue;
//        VectorXd x = X.row(i);
//        int yi = y[i];
//
//        auto mu_phi = [&](const VectorXd &th) {
//          double mu = std::exp(x.dot(th.head(X.cols())));
//          double phi = th(X.cols());
//          return std::make_pair(mu, phi);
//        };
//
//        auto loglik = [&](const VectorXd &th) {
//          auto [mu, phi] = mu_phi(th);
//          return log_ccdf_nb(yi, mu, phi);
//        };
//
//        fpp += loglik(theta_pp);
//        fpm += loglik(theta_pm);
//        fmp += loglik(theta_mp);
//        fmm += loglik(theta_mm);
//      }
//
//      double hij = (fpp - fpm - fmp + fmm) / (4 * eps * eps);
//      H(i, j) = hij;
//      H(j, i) = hij;
//    }
//  }
//
//  return H;
//}
//
//// [[Rcpp::export]]
//List fast_neg_bin_with_censoring_cpp(const Eigen::MatrixXd& X,
//	                                 const Eigen::VectorXi& y,
//	                                 const Eigen::VectorXi& dead,
//	                                 int maxit = 100,
//	                                 double eps_f = 1e-8, 
//	                                 double eps_g = 1e-5) {
//  NBCensoredLogLik fn(X, y, dead);
//  
//  int p_plus_one = X.cols();
//  VectorXd beta_and_theta = VectorXd::Zero(p_plus_one + 1); //initialize theta so we start with the geometric distribution
//
//  VectorXd b_ols = fast_ols_cpp(X, y.cast<double>().array().log1p())["b"];
//  for (int j = 0; j < p_plus_one; j++){
//		beta_and_theta[j] = b_ols[j]; //start at OLS estimates to at least get the signs correct
//  }
//
//  double fmin = 0.0;
//  int status = 0;
//  status = optim_lbfgs(fn, beta_and_theta, fmin, maxit, eps_f, eps_g);
//
//  VectorXd b = beta_and_theta.head(p_plus_one);
//  double theta_hat = std::exp(beta_and_theta[p_plus_one + 1]);
//  
//  return List::create(
//	  Named("b") = 			b,
//	  Named("theta_hat") = 	theta_hat,
//	  Named("logLik") = 	-fmin,
//	  Named("converged") = 	(status == 0)
//  );
//}
//
//// [[Rcpp::export]]
//List fast_neg_bin_with_censoring_with_sd_cpp(const Eigen::MatrixXd& X,
//	                                 const Eigen::VectorXi& y,
//	                                 const Eigen::VectorXi& dead,
//	                                 int maxit = 100,
//	                                 double eps_f = 1e-8, 
//	                                 double eps_g = 1e-5) {
//	List mod = fast_neg_bin_with_censoring_cpp(X, y, dead, maxit, eps_f, eps_g);
//	VectorXd b = mod["b"];
//	double theta_hat = mod["theta_hat"];
//	
//	bool any_censoring = false;
//	bool any_dead = false;
//  	for (int i = 0; i < y.size(); i++){
//		if (dead[i] == 0){
//			any_censoring = true;
//		} else {
//			any_dead = true;
//		}
//	}
//	MatrixXd hess_fisher_info_matrix;	
//	if (any_dead & any_censoring){
//		hess_fisher_info_matrix = analytic_hessian_uncensored(X, y, dead, b, theta_hat) + 
//			numeric_hessian_censored(X, y, dead, b, theta_hat);
//	} else if (any_dead){
//		hess_fisher_info_matrix = analytic_hessian_uncensored(X, y, dead, b, theta_hat);
//	} else if (any_censoring){
//		hess_fisher_info_matrix = numeric_hessian_censored(X, y, dead, b, theta_hat);
//	}	
//	mod["ssq_b_2"] = eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(hess_fisher_info_matrix, 2); //R index
//	
//	return mod;									
//}