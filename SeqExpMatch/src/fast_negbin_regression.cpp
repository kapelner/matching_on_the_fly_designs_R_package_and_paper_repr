#include "_helper_functions.h"
#include <RcppNumerical.h>
using namespace Numer;
using namespace Rcpp;

// Negative binomial censored log-likelihood with gradient
class NBCensoredLogLik : public MFuncGrad {
	private:
	  const Eigen::MatrixXd& X;
	  const Eigen::VectorXi& y;
	  const Eigen::VectorXi& dead_;
	  int p;
	
	public:
	  NBCensoredLogLik(const Eigen::MatrixXd& X_,
	                   const Eigen::VectorXi& y_,
	                   const Eigen::VectorXi& dead_)
	    : X(X_), y(y_), dead_(dead_), p(X_.cols()) {}
	
	  // Negative log-likelihood
//	  double f(Constvec& pars) override {
//	    Eigen::VectorXd beta = pars.head(p);
//	    double log_theta = pars[p];
//	    double theta = std::exp(log_theta);
//	
//	    Eigen::VectorXd eta = X * beta;
//	    double ll = 0.0;
//	
//	    for (int i = 0; i < X.rows(); ++i) {
//	      double mu = std::exp(eta[i]);
//	      if (dead[i] == 1) {
//	        ll += R::dnbinom_mu(y[i], theta, mu, true);
//	      } else {
//	        double surv = 1.0 - R::pnbinom_mu(y[i], theta, mu, true, false);
//	        if (surv > 1e-12) {
//	          ll -= std::log(surv);
//	        } else {
//	          ll += 1e6;  // penalty for numerical instability
//	        }
//	      }
//	    }
//	    return -ll; // Negative log-likelihood for minimization
//	  }
	
	  // Gradient (automatic numerical gradient)
	  double f_grad(Constvec& pars, Refvec grad) {
//	    Eigen::VectorXd beta = pars.head(p);
//	    double log_theta = pars[p];
//	    double theta = std::exp(log_theta);
//	
//	    Eigen::VectorXd eta = X * beta;
//	    g.setZero(pars.size());
//	
//	    for (int i = 0; i < n; ++i) {
//	      double mu = std::exp(eta[i]);
//	      double r = theta;
//	
//	      if (dead[i] == 1) {
//	        double yi = y[i];
//	        double denom = r + mu;
//	        double dlogf_dmu = yi / mu - (yi + r) / denom;
//	        Eigen::VectorXd dmu_dbeta = mu * X.row(i).transpose();
//	        g.head(p) -= dlogf_dmu * dmu_dbeta;
//	
//	        double dlogf_dr = R::digamma(yi + r) - R::digamma(r) +
//	                          std::log(r) - std::log(denom) +
//	                          1 - (yi + r) / denom;
//	        g[p] -= r * dlogf_dr;  // chain rule for log(theta)
//	      } else {
//	        // Use numerical gradient for censored observation
//	        // since derivatives of CDF are intractable
//	        NumericalGradient<NBCensoredLogLik> numgrad(*this);
//	        Eigen::VectorXd g_i(pars.size());
//	        numgrad.grad(pars, g_i, i); // only point i perturbed
//	        g += g_i;
//	      }
		return 0.0;
	  }
};

// [[Rcpp::export]]
List fast_neg_bin_with_censoring_cpp(const Eigen::MatrixXd& X,
	                                 const Eigen::VectorXi& y,
	                                 const Eigen::VectorXi& dead,
	                                 int maxit = 100,
	                                 double eps_f = 1e-8, double eps_g = 1e-5) {
  int p = X.cols();
  Eigen::VectorXd beta = Eigen::VectorXd::Zero(p + 1);
  beta[p] = std::log(1.0); // Initial log(theta)

  NBCensoredLogLik fn(X, y, dead);

  double fmin = 0.0;
  int status = 0;
  //status = optim_lbfgs(fn, beta, fmin, maxit, eps_f, eps_g);

  Eigen::VectorXd beta_hat = beta.head(p);
  double theta_hat = std::exp(beta[p]);
  
  return List::create(
	  Named("coefficients") = 	beta_hat,
	  Named("theta") = 			theta_hat,
	  Named("logLik") = 		-fmin,
	  Named("converged") = 		(status == 0),
	  Named("message") = 		status == 0 ? "Converged" : "Not converged"
  );
}