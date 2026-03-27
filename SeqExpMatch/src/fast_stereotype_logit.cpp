#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

class StereotypeLogitRegression {
private:
    MatrixXd m_X;
    std::vector<int> m_y;
    int m_n;
    int m_p;
    int m_K;

public:
    StereotypeLogitRegression(const MatrixXd& X, const VectorXd& y) :
        m_X(X), m_n(X.rows()), m_p(X.cols()), m_K(0) {
        std::vector<double> levels = init_levels(y);
        m_K = static_cast<int>(levels.size());
        m_y.resize(m_n);
        for (int i = 0; i < m_n; ++i) {
            double yi = y[i];
            int idx = 0;
            while (idx < m_K && levels[idx] != yi) {
                ++idx;
            }
            m_y[i] = idx + 1;
        }
    }

    static std::vector<double> init_levels(const VectorXd& y) {
        std::vector<double> levels(y.data(), y.data() + y.size());
        std::sort(levels.begin(), levels.end());
        levels.erase(std::unique(levels.begin(), levels.end()), levels.end());
        return levels;
    }

    int num_categories() const { return m_K; }
    int num_alpha() const { return m_K - 1; }
    int num_gamma() const { return std::max(0, m_K - 2); }
    int num_params() const { return num_alpha() + m_p + num_gamma(); }

    VectorXd initialize_params() const {
        VectorXd params(num_params());
        params.setZero();

        std::vector<double> counts(m_K, 0.5);
        for (int i = 0; i < m_n; ++i) {
            counts[m_y[i] - 1] += 1.0;
        }
        for (int j = 1; j < m_K; ++j) {
            params[j - 1] = std::log(counts[j] / counts[0]);
        }
        return params;
    }

    void compute_scores(
        const VectorXd& gamma,
        std::vector<double>& score_vals,
        MatrixXd& dscore_dgamma
    ) const {
        score_vals.assign(m_K, 0.0);
        dscore_dgamma.setZero(m_K, num_gamma());

        if (num_gamma() == 0) {
            if (m_K >= 2) {
                score_vals[m_K - 1] = 1.0;
            }
            return;
        }

        VectorXd v = gamma.array().exp();
        double denom = 1.0 + v.sum();
        VectorXd cum_v(num_gamma());
        double running = 0.0;
        for (int r = 0; r < num_gamma(); ++r) {
            running += v[r];
            cum_v[r] = running;
        }

        for (int j = 2; j <= m_K - 1; ++j) {
            int interior_idx = j - 2;
            double c_j = cum_v[interior_idx];
            score_vals[j - 1] = c_j / denom;

            for (int r = 0; r < num_gamma(); ++r) {
                double dc = (r <= interior_idx) ? v[r] : 0.0;
                dscore_dgamma(j - 1, r) = (dc * denom - c_j * v[r]) / (denom * denom);
            }
        }
        score_vals[m_K - 1] = 1.0;
    }

    double loglik_grad(
        const VectorXd& params,
        VectorXd* grad = NULL
    ) const {
        const int n_alpha = num_alpha();
        const int n_gamma = num_gamma();

        VectorXd alpha = params.head(n_alpha);
        VectorXd beta = params.segment(n_alpha, m_p);
        VectorXd gamma = (n_gamma > 0) ? params.tail(n_gamma) : VectorXd(0);

        if (grad != NULL) {
            grad->setZero(num_params());
        }

        std::vector<double> score_vals;
        MatrixXd dscore_dgamma;
        compute_scores(gamma, score_vals, dscore_dgamma);

        std::vector<double> logits(m_K, 0.0);
        std::vector<double> probs(m_K, 0.0);
        double ll = 0.0;

        for (int i = 0; i < m_n; ++i) {
            double eta = (m_p > 0) ? m_X.row(i).dot(beta) : 0.0;
            logits[0] = 0.0;
            for (int j = 2; j <= m_K; ++j) {
                logits[j - 1] = alpha[j - 2] + score_vals[j - 1] * eta;
            }

            double max_logit = *std::max_element(logits.begin(), logits.end());
            double denom = 0.0;
            for (int j = 0; j < m_K; ++j) {
                probs[j] = std::exp(logits[j] - max_logit);
                denom += probs[j];
            }
            for (int j = 0; j < m_K; ++j) {
                probs[j] /= denom;
            }

            int yi = m_y[i] - 1;
            ll += logits[yi] - max_logit - std::log(denom);

            if (grad != NULL) {
                for (int j = 2; j <= m_K; ++j) {
                    (*grad)[j - 2] += ((yi == (j - 1)) ? 1.0 : 0.0) - probs[j - 1];
                }

                if (m_p > 0) {
                    double observed_score = score_vals[yi];
                    double expected_score = 0.0;
                    for (int j = 0; j < m_K; ++j) {
                        expected_score += probs[j] * score_vals[j];
                    }
                    grad->segment(n_alpha, m_p).noalias() += m_X.row(i).transpose() * (observed_score - expected_score);
                }

                if (n_gamma > 0) {
                    for (int r = 0; r < n_gamma; ++r) {
                        double observed_dscore = dscore_dgamma(yi, r);
                        double expected_dscore = 0.0;
                        for (int j = 0; j < m_K; ++j) {
                            expected_dscore += probs[j] * dscore_dgamma(j, r);
                        }
                        (*grad)[n_alpha + m_p + r] += eta * (observed_dscore - expected_dscore);
                    }
                }
            }
        }

        return ll;
    }
};

static MatrixXd numeric_hessian_from_gradient(
    const StereotypeLogitRegression& model,
    const VectorXd& params,
    double h = 1e-5
) {
    const int d = params.size();
    MatrixXd H(d, d);

    for (int j = 0; j < d; ++j) {
        VectorXd grad_plus(d), grad_minus(d);
        VectorXd p_plus = params;
        VectorXd p_minus = params;
        p_plus[j] += h;
        p_minus[j] -= h;
        model.loglik_grad(p_plus, &grad_plus);
        model.loglik_grad(p_minus, &grad_minus);
        H.col(j) = (grad_plus - grad_minus) / (2.0 * h);
    }
    return 0.5 * (H + H.transpose());
}

static MatrixXd pseudo_inverse_symmetric(const MatrixXd& A, double tol = 1e-8) {
    JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
    VectorXd sing = svd.singularValues();
    MatrixXd Dinv = MatrixXd::Zero(sing.size(), sing.size());
    for (int i = 0; i < sing.size(); ++i) {
        if (sing[i] > tol) {
            Dinv(i, i) = 1.0 / sing[i];
        }
    }
    return svd.matrixV() * Dinv * svd.matrixU().transpose();
}

static VectorXd set_beta_and_pack_nuisance(
    const VectorXd& nuisance,
    double beta_fixed,
    int n_alpha,
    int p,
    int n_gamma,
    int beta_index = 0
) {
    const int d = n_alpha + p + n_gamma;
    VectorXd params(d);
    int cursor = 0;

    if (n_alpha > 0) {
        params.head(n_alpha) = nuisance.segment(cursor, n_alpha);
        cursor += n_alpha;
    }

    for (int j = 0; j < p; ++j) {
        if (j == beta_index) {
            params[n_alpha + j] = beta_fixed;
        } else {
            params[n_alpha + j] = nuisance[cursor];
            ++cursor;
        }
    }

    if (n_gamma > 0) {
        params.tail(n_gamma) = nuisance.tail(n_gamma);
    }

    return params;
}

static VectorXd extract_nuisance_params(
    const VectorXd& params,
    int n_alpha,
    int p,
    int n_gamma,
    int beta_index = 0
) {
    VectorXd nuisance(params.size() - 1);
    int cursor = 0;

    if (n_alpha > 0) {
        nuisance.segment(cursor, n_alpha) = params.head(n_alpha);
        cursor += n_alpha;
    }

    for (int j = 0; j < p; ++j) {
        if (j == beta_index) {
            continue;
        }
        nuisance[cursor] = params[n_alpha + j];
        ++cursor;
    }

    if (n_gamma > 0) {
        nuisance.tail(n_gamma) = params.tail(n_gamma);
    }

    return nuisance;
}

static double nuisance_loglik_grad(
    const StereotypeLogitRegression& model,
    const VectorXd& nuisance,
    double beta_fixed,
    int n_alpha,
    int p,
    int n_gamma,
    VectorXd* grad = NULL,
    int beta_index = 0
) {
    VectorXd params = set_beta_and_pack_nuisance(nuisance, beta_fixed, n_alpha, p, n_gamma, beta_index);
    VectorXd full_grad(params.size());
    double ll = model.loglik_grad(params, grad == NULL ? NULL : &full_grad);
    if (grad != NULL) {
        *grad = extract_nuisance_params(full_grad, n_alpha, p, n_gamma, beta_index);
    }
    return ll;
}

static MatrixXd nuisance_numeric_hessian_from_gradient(
    const StereotypeLogitRegression& model,
    const VectorXd& nuisance,
    double beta_fixed,
    int n_alpha,
    int p,
    int n_gamma,
    double h = 1e-5,
    int beta_index = 0
) {
    const int d = nuisance.size();
    MatrixXd H(d, d);

    for (int j = 0; j < d; ++j) {
        VectorXd grad_plus(d), grad_minus(d);
        VectorXd p_plus = nuisance;
        VectorXd p_minus = nuisance;
        p_plus[j] += h;
        p_minus[j] -= h;
        nuisance_loglik_grad(model, p_plus, beta_fixed, n_alpha, p, n_gamma, &grad_plus, beta_index);
        nuisance_loglik_grad(model, p_minus, beta_fixed, n_alpha, p, n_gamma, &grad_minus, beta_index);
        H.col(j) = (grad_plus - grad_minus) / (2.0 * h);
    }
    return 0.5 * (H + H.transpose());
}

static VectorXd optimize_nuisance_given_beta(
    const StereotypeLogitRegression& model,
    const VectorXd& full_params_start,
    double beta_fixed,
    int n_alpha,
    int p,
    int n_gamma,
    int maxit = 50,
    double tol = 1e-8,
    int beta_index = 0
) {
    VectorXd nuisance = extract_nuisance_params(full_params_start, n_alpha, p, n_gamma, beta_index);
    const int d = nuisance.size();
    if (d == 0) {
        return nuisance;
    }

    VectorXd grad(d);
    for (int iter = 0; iter < maxit; ++iter) {
        double current_ll = nuisance_loglik_grad(model, nuisance, beta_fixed, n_alpha, p, n_gamma, &grad, beta_index);
        if (grad.norm() < tol) {
            break;
        }

        MatrixXd H = nuisance_numeric_hessian_from_gradient(
            model, nuisance, beta_fixed, n_alpha, p, n_gamma, 1e-5, beta_index
        );
        FullPivLU<MatrixXd> lu(H);
        if (!lu.isInvertible()) {
            break;
        }

        VectorXd step = lu.solve(grad);
        double scale = 1.0;
        bool accepted = false;
        while (scale > 1e-8) {
            VectorXd next_nuisance = nuisance - scale * step;
            double next_ll = nuisance_loglik_grad(
                model, next_nuisance, beta_fixed, n_alpha, p, n_gamma, NULL, beta_index
            );
            if (R_finite(next_ll) && next_ll > current_ll) {
                nuisance = next_nuisance;
                accepted = true;
                break;
            }
            scale *= 0.5;
        }

        if (!accepted || (scale * step).norm() < tol) {
            break;
        }
    }

    return nuisance;
}

static double profile_loglik_for_beta(
    const StereotypeLogitRegression& model,
    const VectorXd& full_params_start,
    double beta_fixed,
    int n_alpha,
    int p,
    int n_gamma,
    int beta_index = 0
) {
    VectorXd nuisance = optimize_nuisance_given_beta(
        model, full_params_start, beta_fixed, n_alpha, p, n_gamma, 50, 1e-8, beta_index
    );
    VectorXd params = set_beta_and_pack_nuisance(nuisance, beta_fixed, n_alpha, p, n_gamma, beta_index);
    return model.loglik_grad(params, NULL);
}

static VectorXd stereotype_newton_fit(
    const StereotypeLogitRegression& model,
    int maxit,
    double tol,
    bool* converged = NULL
) {
    VectorXd params = model.initialize_params();
    const int d = params.size();
    VectorXd grad(d);

    if (converged != NULL) {
        *converged = false;
    }

    for (int iter = 0; iter < maxit; ++iter) {
        double current_ll = model.loglik_grad(params, &grad);
        if (grad.norm() < tol) {
            if (converged != NULL) {
                *converged = true;
            }
            break;
        }

        MatrixXd H = numeric_hessian_from_gradient(model, params);
        FullPivLU<MatrixXd> lu(H);
        if (!lu.isInvertible()) {
            break;
        }

        VectorXd step = lu.solve(grad);
        double scale = 1.0;
        bool accepted = false;

        while (scale > 1e-8) {
            VectorXd next_params = params - scale * step;
            double next_ll = model.loglik_grad(next_params, NULL);
            if (R_finite(next_ll) && next_ll > current_ll) {
                params = next_params;
                accepted = true;
                break;
            }
            scale *= 0.5;
        }

        if (!accepted) {
            break;
        }

        if ((scale * step).norm() < tol) {
            if (converged != NULL) {
                *converged = true;
            }
            break;
        }
    }

    return params;
}

// [[Rcpp::export]]
List fast_stereotype_logit_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int maxit = 100, double tol = 1e-8) {
    StereotypeLogitRegression model(X, y);
    if (model.num_categories() < 2) {
        stop("Stereotype logistic regression requires at least two observed outcome categories.");
    }

    bool converged = false;
    VectorXd params = stereotype_newton_fit(model, maxit, tol, &converged);
    int n_alpha = model.num_alpha();
    int p = X.cols();

    return List::create(
        Named("b") = params.segment(n_alpha, p),
        Named("alpha") = params.head(n_alpha),
        Named("scores_raw") = params.tail(model.num_gamma()),
        Named("params") = params,
        Named("converged") = converged
    );
}

// [[Rcpp::export]]
List fast_stereotype_logit_with_var_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int maxit = 100, double tol = 1e-8) {
    StereotypeLogitRegression model(X, y);
    if (model.num_categories() < 2) {
        stop("Stereotype logistic regression requires at least two observed outcome categories.");
    }

    bool converged = false;
    VectorXd params = stereotype_newton_fit(model, maxit, tol, &converged);
    MatrixXd H = numeric_hessian_from_gradient(model, params);
    MatrixXd info = -H;
    FullPivLU<MatrixXd> lu(info);

    int n_alpha = model.num_alpha();
    int p = X.cols();

    MatrixXd cov;
    if (lu.isInvertible()) {
        cov = lu.inverse();
    } else {
        cov = pseudo_inverse_symmetric(info);
    }
    double ssq_b_1 = (p >= 1) ? cov(n_alpha, n_alpha) : NA_REAL;
    if ((!R_finite(ssq_b_1) || ssq_b_1 <= 0) && p >= 1) {
        double beta_hat = params[n_alpha];
        double h = std::max(1e-4, 1e-3 * (std::abs(beta_hat) + 1.0));
        double ll_0 = profile_loglik_for_beta(model, params, beta_hat, n_alpha, p, model.num_gamma(), 0);
        double ll_p = profile_loglik_for_beta(model, params, beta_hat + h, n_alpha, p, model.num_gamma(), 0);
        double ll_m = profile_loglik_for_beta(model, params, beta_hat - h, n_alpha, p, model.num_gamma(), 0);
        double info_beta = -(ll_p - 2.0 * ll_0 + ll_m) / (h * h);
        if (R_finite(info_beta) && info_beta > 0) {
            ssq_b_1 = 1.0 / info_beta;
        }
    }
    if (!R_finite(ssq_b_1) || ssq_b_1 <= 0) {
        ssq_b_1 = NA_REAL;
    }

    return List::create(
        Named("b") = params.segment(n_alpha, p),
        Named("ssq_b_1") = ssq_b_1,
        Named("converged") = converged
    );
}

// [[Rcpp::export]]
double fast_stereotype_profile_loglik_cpp(
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& y,
    double beta_fixed,
    int maxit = 100,
    double tol = 1e-8
) {
    StereotypeLogitRegression model(X, y);
    if (model.num_categories() < 2) {
        stop("Stereotype logistic regression requires at least two observed outcome categories.");
    }

    bool converged = false;
    VectorXd params = stereotype_newton_fit(model, maxit, tol, &converged);
    int n_alpha = model.num_alpha();
    int p = X.cols();
    int n_gamma = model.num_gamma();
    return profile_loglik_for_beta(model, params, beta_fixed, n_alpha, p, n_gamma, 0);
}

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]

//' Parallel Stereotype Logit Randomization Distribution
//'
//' @param y Numeric vector of response values (pre-null-shifted for treated).
//' @param X_covars Matrix of covariates (without intercept or treatment).
//' @param w_mat Integer matrix of permuted treatment assignments (n x nsim).
//' @param delta Null treatment effect (additive shift).
//' @param num_cores Number of OpenMP threads.
//' @return Numeric vector of length nsim with treatment coefficients.
// [[Rcpp::export]]
NumericVector compute_stereotype_logit_distr_parallel_cpp(
	const Eigen::VectorXd& y,
	const Eigen::MatrixXd& X_covars,
	const Rcpp::IntegerMatrix& w_mat,
	double delta,
	int num_cores
) {
	int nsim = w_mat.cols();
	int n = y.size();
	int p_covars = X_covars.cols();
	int p_full = p_covars + 1; // treatment + covars (no intercept — thresholds handle location)

	std::vector<double> results(nsim, NA_REAL);
	const int* w_ptr = w_mat.begin();

#ifdef _OPENMP
	omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(static)
	for (int b = 0; b < nsim; ++b) {
		const int* w_col = w_ptr + (size_t)b * n;

		Eigen::MatrixXd X_full(n, p_full);
		Eigen::VectorXd y_shifted(n);

		for (int i = 0; i < n; ++i) {
			X_full(i, 0) = (double)w_col[i];
			for (int k = 0; k < p_covars; ++k) {
				X_full(i, 1 + k) = X_covars(i, k);
			}
			y_shifted[i] = (w_col[i] == 1) ? y[i] + delta : y[i];
		}

		StereotypeLogitRegression model(X_full, y_shifted);
		if (model.num_categories() < 2) continue;

		bool converged = false;
		Eigen::VectorXd params = stereotype_newton_fit(model, 100, 1e-8, &converged);

		int n_alpha = model.num_alpha();
		// accept result even if not formally converged, matching R generate_mod behaviour
		if ((int)params.size() >= n_alpha + 1 && std::isfinite(params[n_alpha])) {
			results[b] = params[n_alpha];
		}
	}

	return wrap(results);
}
