#include "_helper_functions.h"
#include <RcppNumerical.h>
using namespace Rcpp;
using namespace Numer;

typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;

//see https://github.com/yixuan/RcppNumerical/blob/6ad26382a3414c248c9562c92985bb9e82fa1f04/src/fastLR.cpp
class LogisticReg: public MFuncGrad {
	private:
	   const MapMat X;
	   const MapVec Y;
	   const int n;
	   Eigen::VectorXd xbeta;  // contains X*beta
	   Eigen::VectorXd prob;   // contains log(1+exp(X*beta)) and 1/(1+exp(-X*beta))
	public:
	   LogisticReg(const MapMat x_, const MapVec y_) :
	       X(x_),
	       Y(y_),
	       n(X.rows()),
	       xbeta(n),
	       prob(n)
	   {}
	
	   double f_grad(Constvec& beta, Refvec grad){
	       // Negative log likelihood
	       //   sum(log(1 + exp(X * beta))) - y' * X * beta
	       xbeta.noalias() = X * beta;
	       const double yxbeta = Y.dot(xbeta);
	       // Calculate log(1 + exp(X * beta)), avoiding overflow
	       for(int i = 0; i < n; i++)
	           prob[i] = R::log1pexp(xbeta[i]);
	       const double f = prob.sum() - yxbeta;
	
	       // Gradient
	       //   X' * (p - y), p = exp(X * beta) / (1 + exp(X * beta))
	       //                   = exp(X * beta - log(1 + exp(X * beta)))
	       prob = (xbeta - prob).array().exp();
	       grad.noalias() = X.transpose() * (prob - Y);
	
	       return f;
	   }
	
	   Eigen::VectorXd current_xb() const { return xbeta; }
	   Eigen::VectorXd current_p()  const { return prob; }
};


// [[Rcpp::export]]
List fast_logistic_regression_cpp(NumericMatrix X, NumericVector y, NumericVector start, double eps_f = 1e-8, double eps_g = 1e-5, int maxit = 300) {
   const MapMat xx = Rcpp::as<MapMat>(X);
   const MapVec yy = Rcpp::as<MapVec>(y);
   // Negative log likelihood
   LogisticReg nll(xx, yy);
   // Initial guess
   NumericVector b = Rcpp::clone(start);
   MapVec beta(b.begin(), b.length());

   double fopt;
   int status = optim_lbfgs(nll, beta, fopt, maxit, eps_f, eps_g);

   return Rcpp::List::create(
       Named("b")      		= b,
       Named("converged")   = (status >= 0),
       Named("xx")         	= xx,
       Named("yy")         	= yy
   );
}

//////not finished yet!!!
// List fast_logistic_regression_with_sd_cpp(Rcpp::NumericMatrix X, Rcpp::NumericVector y, Rcpp::NumericVector start, double eps_f = 1e-8, double eps_g = 1e-5, int maxit = 300) {
//   List mod = fast_logistic_regression_cpp(X, y, start, eps_f, eps_g, maxit);
//   Eigen::VectorXd b = mod["b"];
// 
//   Eigen::VectorXd w(X.ncol());
//   // for (int j = 0; j < X.ncol(); j++){
//   //   Eigen::VectorXd X_dot_b_j = (X(, j) * b)(1, 1);
//   // 	w[j] = X_dot_b_j / ((1 + X_dot_b_j) * (1 + X_dot_b_j));
//   // }
//   //exp_X_dot_b = exp(X %*% mod$b)
//   //w = as.numeric(exp_Xmdot_b / (1 + exp_X_dot_b)^2)
// 
//   //const MapVec ww(w.data(), w.size());
//   Eigen::MatrixXd XmmtWmatXmm = eigen_Xt_times_diag_w_times_X(X, w); //t(Xmm) %*% diag(w) %*% Xmm
// 
//   mod["s_2"] = sqrt(eigen_compute_single_entry_of_diagonal_matrix(XmmtWmatXmm, 2));
// 
//   return mod;
// }

