#include "_helper_functions.h"
using namespace Rcpp;
using namespace Eigen;

// Fast optimizer for the combined conditional-Poisson + Poisson log-likelihood:
//   L_total = L_cond_Poisson(pairs) + L_Poisson(reservoir)
//
// Pair component (conditional Poisson = weighted binomial logistic):
//   eta_k = beta_T + X_diff_k' beta_xs
//   ll_pairs = sum(yT_k * eta_k - n_k * log(1 + exp(eta_k)))
//
// Reservoir component (marginal Poisson):
//   eta_i = beta_0 + w_i * beta_T + X_i' beta_xs
//   ll_res = sum(y_i * eta_i - exp(eta_i))
//
// Parameter layout: [beta_0 (0), beta_T (1), beta_xs (2..p+1)]  (size p+2)
// beta_T is at 0-based index 1; pass j=2 (1-based) to extract its variance.
//
// Uses Newton's method with the analytic Fisher-information matrix as Hessian.
// At each iteration:
//   grad = gradient of neg-log-lik (= -score)
//   H    = Fisher info = X_pairs_eff' W_pairs X_pairs_eff + X_res_eff' W_res X_res_eff
//   params -= H^{-1} * grad
//
// Converges in O(1) iterations near the optimum (quadratic convergence).
// At convergence H = observed information; vcov = H^{-1}; ssq_b_j = H^{-1}[1,1].
//

static List cpoisson_combined_score_info_cpp_impl(
	const Eigen::VectorXd& yT_v,
	const Eigen::VectorXd& n_k_v,
	const Eigen::MatrixXd& X_diff_v,
	const Eigen::VectorXd& y_r,
	const Eigen::VectorXd& w_r,
	const Eigen::MatrixXd& X_r,
	const Eigen::VectorXd& params
) {
	const int nd = (int)yT_v.size();
	const int nR = (int)y_r.size();
	const int p  = (int)X_diff_v.cols();
	const int np = p + 2;

	VectorXd score = VectorXd::Zero(np);
	MatrixXd info = MatrixXd::Zero(np, np);

	const double beta_0 = params[0];
	const double beta_T = params[1];
	const VectorXd beta_xs = params.tail(p);

	VectorXd eta_p = VectorXd::Constant(nd, beta_T);
	if (p > 0) eta_p.noalias() += X_diff_v * beta_xs;
	ArrayXd p_k_arr = 1.0 / (1.0 + (-eta_p.array()).exp());
	ArrayXd w_p_arr = n_k_v.array() * p_k_arr * (1.0 - p_k_arr);
	VectorXd score_p = (yT_v.array() - n_k_v.array() * p_k_arr).matrix();
	score[1] += score_p.sum();
	if (p > 0) score.tail(p).noalias() += X_diff_v.transpose() * score_p;

	VectorXd w_p_vec = w_p_arr.matrix();
	info(1, 1) += w_p_vec.sum();
	if (p > 0) {
		VectorXd Xdw = X_diff_v.transpose() * w_p_vec;
		info.block(1, 2, 1, p).noalias() += Xdw.transpose();
		info.block(2, 1, p, 1).noalias() += Xdw;
		info.block(2, 2, p, p).noalias() += X_diff_v.transpose() * w_p_vec.asDiagonal() * X_diff_v;
	}

	VectorXd eta_r = VectorXd::Constant(nR, beta_0) + beta_T * w_r;
	if (p > 0) eta_r.noalias() += X_r * beta_xs;
	VectorXd mu_r = eta_r.array().exp();
	VectorXd score_r = y_r - mu_r;
	score[0] += score_r.sum();
	score[1] += score_r.dot(w_r);
	if (p > 0) score.tail(p).noalias() += X_r.transpose() * score_r;

	const double mu_sum = mu_r.sum();
	const double muw_sum = mu_r.dot(w_r);
	info(0, 0) += mu_sum;
	info(0, 1) += muw_sum;
	info(1, 0) += muw_sum;
	info(1, 1) += (mu_r.array() * w_r.array().square()).sum();
	if (p > 0) {
		VectorXd Xrmu = X_r.transpose() * mu_r;
		VectorXd Xrwmu = X_r.transpose() * mu_r.cwiseProduct(w_r);
		info.block(0, 2, 1, p).noalias() += Xrmu.transpose();
		info.block(2, 0, p, 1).noalias() += Xrmu;
		info.block(1, 2, 1, p).noalias() += Xrwmu.transpose();
		info.block(2, 1, p, 1).noalias() += Xrwmu;
		info.block(2, 2, p, p).noalias() += X_r.transpose() * mu_r.asDiagonal() * X_r;
	}

	return List::create(Named("score") = score, Named("info") = info);
}

static double cpoisson_combined_neg_loglik_cpp_impl(
	const Eigen::VectorXd& yT_v,
	const Eigen::VectorXd& n_k_v,
	const Eigen::MatrixXd& X_diff_v,
	const Eigen::VectorXd& y_r,
	const Eigen::VectorXd& w_r,
	const Eigen::MatrixXd& X_r,
	const Eigen::VectorXd& params
) {
	const int nd = (int)yT_v.size();
	const int nR = (int)y_r.size();
	const int p  = (int)X_diff_v.cols();

	const double beta_0 = params[0];
	const double beta_T = params[1];
	const VectorXd beta_xs = params.tail(p);

	double loglik = 0.0;
	VectorXd eta_p = VectorXd::Constant(nd, beta_T);
	if (p > 0) eta_p.noalias() += X_diff_v * beta_xs;
	for (int i = 0; i < nd; ++i) {
		const double eta = eta_p[i];
		const double log_denom = eta > 0.0 ? eta + std::log1p(std::exp(-eta)) : std::log1p(std::exp(eta));
		loglik += yT_v[i] * eta - n_k_v[i] * log_denom;
	}

	VectorXd eta_r = VectorXd::Constant(nR, beta_0) + beta_T * w_r;
	if (p > 0) eta_r.noalias() += X_r * beta_xs;
	for (int i = 0; i < nR; ++i) {
		loglik += y_r[i] * eta_r[i] - std::exp(eta_r[i]) - R::lgammafn(y_r[i] + 1.0);
	}
	return -loglik;
}

// [[Rcpp::export]]
Eigen::VectorXd get_cpoisson_combined_score_cpp(
	const Eigen::VectorXd& yT_v,
	const Eigen::VectorXd& n_k_v,
	const Eigen::MatrixXd& X_diff_v,
	const Eigen::VectorXd& y_r,
	const Eigen::VectorXd& w_r,
	const Eigen::MatrixXd& X_r,
	const Eigen::VectorXd& params
) {
	List out = cpoisson_combined_score_info_cpp_impl(yT_v, n_k_v, X_diff_v, y_r, w_r, X_r, params);
	return Rcpp::as<Eigen::VectorXd>(out["score"]);
}

// [[Rcpp::export]]
Eigen::MatrixXd get_cpoisson_combined_hessian_cpp(
	const Eigen::VectorXd& yT_v,
	const Eigen::VectorXd& n_k_v,
	const Eigen::MatrixXd& X_diff_v,
	const Eigen::VectorXd& y_r,
	const Eigen::VectorXd& w_r,
	const Eigen::MatrixXd& X_r,
	const Eigen::VectorXd& params
) {
	List out = cpoisson_combined_score_info_cpp_impl(yT_v, n_k_v, X_diff_v, y_r, w_r, X_r, params);
	MatrixXd info = Rcpp::as<Eigen::MatrixXd>(out["info"]);
	return -info;
}

// [[Rcpp::export]]
List fast_cpoisson_combined_with_var_cpp(
	const Eigen::VectorXd& yT_v,       // treated count per valid pair (nd)
	const Eigen::VectorXd& n_k_v,      // total count per valid pair (nd)
	const Eigen::MatrixXd& X_diff_v,   // covariate diffs (nd x p; p=0 valid)
	const Eigen::VectorXd& y_r,        // reservoir outcomes (nR)
	const Eigen::VectorXd& w_r,        // reservoir treatment indicator (nR)
	const Eigen::MatrixXd& X_r,        // reservoir covariates (nR x p)
	int    maxit = 100,
	double tol   = 1e-8,
	Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue
) {
	const int nd = (int)yT_v.size();
	const int nR = (int)y_r.size();
	const int p  = (int)X_diff_v.cols();
	const int np = p + 2;              // beta_0, beta_T, beta_xs (p)

	// ---- Initialise params -----------------------------------------------
	VectorXd params = VectorXd::Zero(np);
	if (nR > 0) params[0] = std::log(std::max(1.0, y_r.mean()));
	FixedParamSpec fixed_spec = make_fixed_param_spec(np, fixed_idx, fixed_values);
	for (int k = 0; k < fixed_spec.fixed_idx.size(); ++k) {
		params[fixed_spec.fixed_idx[k]] = fixed_spec.fixed_values[k];
	}

	VectorXd grad(np);
	MatrixXd H(np, np);
	bool converged = false;

	for (int iter = 0; iter < maxit; ++iter) {
		const double  beta_0  = params[0];
		const double  beta_T  = params[1];
		const VectorXd beta_xs = params.tail(p);

		grad.setZero();
		H.setZero();

		// ---- Pair component (conditional Poisson) -------------------------
		// X_pairs_eff row k = [0, 1, X_diff_v[k,:]]  →  cols 1..np-1 only
		VectorXd eta_p = VectorXd::Constant(nd, beta_T);
		if (p > 0) eta_p.noalias() += X_diff_v * beta_xs;

		// Logistic probabilities and Fisher weights
		ArrayXd  p_k_arr  = 1.0 / (1.0 + (-eta_p.array()).exp());
		ArrayXd  w_p_arr  = n_k_v.array() * p_k_arr * (1.0 - p_k_arr);
		VectorXd resid_p  = (n_k_v.array() * p_k_arr - yT_v.array()).matrix();

		// Gradient contributions from pairs  (d(-L)/d beta_T and d beta_xs)
		grad[1] += resid_p.sum();
		if (p > 0) grad.tail(p).noalias() += X_diff_v.transpose() * resid_p;

		// Fisher block H[1..np-1, 1..np-1]:
		//   [1 | X_diff_v]' * diag(w_p) * [1 | X_diff_v]
		VectorXd w_p_vec = w_p_arr.matrix();
		H(1, 1) += w_p_vec.sum();
		if (p > 0) {
			VectorXd Xdw = X_diff_v.transpose() * w_p_vec;        // p-vector
			H.block(1, 2, 1, p).noalias() += Xdw.transpose();
			H.block(2, 1, p, 1).noalias() += Xdw;
			H.block(2, 2, p, p).noalias() += X_diff_v.transpose() * w_p_vec.asDiagonal() * X_diff_v;
		}

		// ---- Reservoir component (marginal Poisson) -----------------------
		// X_res_eff row i = [1, w_r[i], X_r[i,:]]
		VectorXd eta_r = VectorXd::Constant(nR, beta_0) + beta_T * w_r;
		if (p > 0) eta_r.noalias() += X_r * beta_xs;

		VectorXd mu_r   = eta_r.array().exp();
		VectorXd resid_r = mu_r - y_r;

		// Gradient contributions from reservoir
		grad[0] += resid_r.sum();
		grad[1] += resid_r.dot(w_r);
		if (p > 0) grad.tail(p).noalias() += X_r.transpose() * resid_r;

		// Fisher block H += [1|w_r|X_r]' * diag(mu_r) * [1|w_r|X_r]
		const double mu_sum  = mu_r.sum();
		const double muw_sum = mu_r.dot(w_r);
		H(0, 0) += mu_sum;
		H(0, 1) += muw_sum;
		H(1, 0) += muw_sum;
		H(1, 1) += (mu_r.array() * w_r.array().square()).sum();
		if (p > 0) {
			VectorXd Xrmu  = X_r.transpose() * mu_r;
			VectorXd Xrwmu = X_r.transpose() * mu_r.cwiseProduct(w_r);
			H.block(0, 2, 1, p).noalias() += Xrmu.transpose();
			H.block(2, 0, p, 1).noalias() += Xrmu;
			H.block(1, 2, 1, p).noalias() += Xrwmu.transpose();
			H.block(2, 1, p, 1).noalias() += Xrwmu;
			H.block(2, 2, p, p).noalias() += X_r.transpose() * mu_r.asDiagonal() * X_r;
		}

		// ---- Newton step over free parameters ----------------------------
		VectorXd delta_full = VectorXd::Zero(np);
		if (fixed_spec.free_idx.size() > 0) {
			MatrixXd H_free = subset_matrix(H, fixed_spec.free_idx, fixed_spec.free_idx);
			VectorXd grad_free(fixed_spec.free_idx.size());
			for (int k = 0; k < fixed_spec.free_idx.size(); ++k) grad_free[k] = grad[fixed_spec.free_idx[k]];
			VectorXd delta_free = H_free.ldlt().solve(grad_free);
			for (int k = 0; k < fixed_spec.free_idx.size(); ++k) {
				delta_full[fixed_spec.free_idx[k]] = delta_free[k];
				params[fixed_spec.free_idx[k]] -= delta_free[k];
			}
		}
		for (int k = 0; k < fixed_spec.fixed_idx.size(); ++k) {
			params[fixed_spec.fixed_idx[k]] = fixed_spec.fixed_values[k];
		}

		if (delta_full.norm() < tol) {
			converged = true;
			break;
		}
	}

	// ---- Extract Var(beta_T) from H^{-1}[1,1] (1-based index 2) ---------
	List final_si = cpoisson_combined_score_info_cpp_impl(yT_v, n_k_v, X_diff_v, y_r, w_r, X_r, params);
	VectorXd score = Rcpp::as<VectorXd>(final_si["score"]);
	MatrixXd info = Rcpp::as<MatrixXd>(final_si["info"]);
	MatrixXd info_free = subset_matrix(info, fixed_spec.free_idx, fixed_spec.free_idx);
	MatrixXd cov_free = covariance_from_information(info_free);
	MatrixXd vcov = expand_free_covariance(np, fixed_spec, cov_free, true);
	double ssq_b_j = (np >= 2) ? vcov(1, 1) : NA_REAL;
	double neg_loglik = cpoisson_combined_neg_loglik_cpp_impl(yT_v, n_k_v, X_diff_v, y_r, w_r, X_r, params);
	if (!converged || !std::isfinite(ssq_b_j)) {
		Rcout << "DEBUG C++: converged=" << converged << " ssq_b_j=" << ssq_b_j << std::endl;
	}

	return List::create(
		Named("b")         = params,
		Named("params")    = params,
		Named("ssq_b_j")   = ssq_b_j,
		Named("vcov")      = vcov,
		Named("score")     = score,
		Named("information") = info,
		Named("hessian")   = -info,
		Named("neg_loglik") = neg_loglik,
		Named("neg_ll")    = neg_loglik,
		Named("loglik")    = R_finite(neg_loglik) ? -neg_loglik : NA_REAL,
		Named("converged") = converged
	);
}
