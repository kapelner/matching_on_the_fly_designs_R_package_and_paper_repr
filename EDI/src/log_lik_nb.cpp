#include <RcppEigen.h>
#include <cmath>    // for std::log, std::exp, std::lgamma

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
double neg_loglik_nb_cpp(double theta,
							SEXP beta_sexp,
			                SEXP X_sexp,
			                SEXP y_sexp
	                 	) {
	using namespace Rcpp;
	NumericVector beta_r(beta_sexp);
	NumericMatrix X_r(X_sexp);
	IntegerVector y_r(y_sexp);
	Eigen::Map<const Eigen::VectorXd> beta(beta_r.begin(), beta_r.size());
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
	Eigen::Map<const Eigen::VectorXi> y(y_r.begin(), y_r.size());

	Eigen::VectorXd eta = X * beta;
	Eigen::VectorXd mu = eta.array().exp();

	double ll = 0.0;
	int n = y.size();

	for (int i = 0; i < n; i++) {
	double yi = static_cast<double>(y[i]);
	double mui = mu[i];

	ll += std::lgamma(yi + theta)
		- std::lgamma(theta)
		- std::lgamma(yi + 1.0)
		+ theta * std::log(theta)
		+ yi * std::log(mui)
		- (yi + theta) * std::log(mui + theta);
	}

	return -ll;
}
