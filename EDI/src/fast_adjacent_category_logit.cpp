#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

class AdjacentCategoryLogitRegression {
private:
    MatrixXd m_X;
    std::vector<int> m_y;
    int m_n;
    int m_p;
    int m_K;

public:
    AdjacentCategoryLogitRegression(const MatrixXd& X, const VectorXd& y) :
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

    int num_params() const {
        return (m_K - 1) + m_p;
    }

    int num_categories() const {
        return m_K;
    }

    double loglik_grad_hess(const VectorXd& params, VectorXd* grad = NULL, MatrixXd* hess = NULL) const {
        const int n_alpha = m_K - 1;
        const VectorXd alpha = params.head(n_alpha);
        const VectorXd beta = params.tail(m_p);

        if (grad != NULL) {
            grad->setZero(params.size());
        }
        if (hess != NULL) {
            hess->setZero(params.size(), params.size());
        }

        double ll = 0.0;
        std::vector<double> alpha_suffix(m_K, 0.0);
        for (int j = m_K - 2; j >= 0; --j) {
            alpha_suffix[j] = alpha_suffix[j + 1] + alpha[j];
        }

        std::vector<double> scores(m_K, 0.0);
        std::vector<double> probs(m_K, 0.0);
        std::vector<double> cdf(m_K - 1, 0.0);
        std::vector<double> prefix_first_moment(m_K - 1, 0.0);

        for (int i = 0; i < m_n; ++i) {
            double eta = (m_p > 0) ? m_X.row(i).dot(beta) : 0.0;
            double max_score = 0.0;

            for (int j = 0; j < m_K - 1; ++j) {
                scores[j] = alpha_suffix[j] - static_cast<double>(m_K - 1 - j) * eta;
            }
            scores[m_K - 1] = 0.0;
            max_score = *std::max_element(scores.begin(), scores.end());

            double denom = 0.0;
            for (int j = 0; j < m_K; ++j) {
                probs[j] = std::exp(scores[j] - max_score);
                denom += probs[j];
            }
            for (int j = 0; j < m_K; ++j) {
                probs[j] /= denom;
            }

            int y_i = m_y[i];
            ll += scores[y_i - 1] - max_score - std::log(denom);

            double ey = 0.0;
            double ey2 = 0.0;
            double running_cdf = 0.0;
            double running_first_moment = 0.0;
            for (int j = 0; j < m_K; ++j) {
                double y_val = static_cast<double>(j + 1);
                ey += y_val * probs[j];
                ey2 += y_val * y_val * probs[j];
                if (j < m_K - 1) {
                    running_cdf += probs[j];
                    running_first_moment += y_val * probs[j];
                    cdf[j] = running_cdf;
                    prefix_first_moment[j] = running_first_moment;
                }
            }
            double var_y = std::max(0.0, ey2 - ey * ey);

            if (grad != NULL) {
                for (int j = 0; j < m_K - 1; ++j) {
                    (*grad)[j] += ((y_i <= (j + 1)) ? 1.0 : 0.0) - cdf[j];
                }
                if (m_p > 0) {
                    grad->tail(m_p).noalias() += m_X.row(i).transpose() * (static_cast<double>(y_i) - ey);
                }
            }

            if (hess != NULL) {
                for (int j = 0; j < m_K - 1; ++j) {
                    for (int k = j; k < m_K - 1; ++k) {
                        double f_min = cdf[std::min(j, k)];
                        double val = -(f_min - cdf[j] * cdf[k]);
                        (*hess)(j, k) += val;
                        if (j != k) {
                            (*hess)(k, j) += val;
                        }
                    }
                }

                if (m_p > 0) {
                    for (int j = 0; j < m_K - 1; ++j) {
                        double cov_ind_y = prefix_first_moment[j] - ey * cdf[j];
                        VectorXd cross = m_X.row(i).transpose() * cov_ind_y;
                        hess->block(j, n_alpha, 1, m_p) += cross.transpose();
                        hess->block(n_alpha, j, m_p, 1) += cross;
                    }

                    hess->block(n_alpha, n_alpha, m_p, m_p).noalias() -= var_y * (m_X.row(i).transpose() * m_X.row(i));
                }
            }
        }

        return ll;
    }
};

static VectorXd adjacent_category_newton_fit(
    const AdjacentCategoryLogitRegression& model,
    int maxit,
    double tol,
    bool* converged = NULL
) {
    const int n_params = model.num_params();
    const int K = model.num_categories();
    const int n_alpha = K - 1;

    VectorXd params(n_params);
    params.setZero();

    VectorXd grad(n_params);
    MatrixXd hess(n_params, n_params);

    if (converged != NULL) {
        *converged = false;
    }

    for (int iter = 0; iter < maxit; ++iter) {
        double current_ll = model.loglik_grad_hess(params, &grad, &hess);
        if (grad.norm() < tol) {
            if (converged != NULL) {
                *converged = true;
            }
            break;
        }

        FullPivLU<MatrixXd> lu(hess);
        if (!lu.isInvertible()) {
            break;
        }

        VectorXd step = lu.solve(grad);
        double step_scale = 1.0;
        bool accepted = false;

        while (step_scale > 1e-8) {
            VectorXd next_params = params - step_scale * step;
            double next_ll = model.loglik_grad_hess(next_params, NULL, NULL);
            if (R_finite(next_ll) && next_ll > current_ll) {
                params = next_params;
                accepted = true;
                break;
            }
            step_scale *= 0.5;
        }

        if (!accepted) {
            break;
        }

        if ((step_scale * step).norm() < tol) {
            if (converged != NULL) {
                *converged = true;
            }
            break;
        }
    }

    return params;
}

// [[Rcpp::export]]
List fast_adjacent_category_logit_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int maxit = 100, double tol = 1e-8) {
    AdjacentCategoryLogitRegression model(X, y);
    if (model.num_categories() < 2) {
        stop("Adjacent-category logits require at least two observed outcome categories.");
    }

    bool converged = false;
    VectorXd params = adjacent_category_newton_fit(model, maxit, tol, &converged);
    int n_alpha = model.num_categories() - 1;

    return List::create(
        Named("b") = params.tail(X.cols()),
        Named("alpha") = params.head(n_alpha),
        Named("params") = params,
        Named("converged") = converged
    );
}

// [[Rcpp::export]]
List fast_adjacent_category_logit_with_var_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int maxit = 100, double tol = 1e-8) {
    AdjacentCategoryLogitRegression model(X, y);
    if (model.num_categories() < 2) {
        stop("Adjacent-category logits require at least two observed outcome categories.");
    }

    bool converged = false;
    VectorXd params = adjacent_category_newton_fit(model, maxit, tol, &converged);
    VectorXd grad(model.num_params());
    MatrixXd hess(model.num_params(), model.num_params());
    model.loglik_grad_hess(params, &grad, &hess);

    MatrixXd info = -hess;
    FullPivLU<MatrixXd> lu(info);

    if (!lu.isInvertible()) {
        return List::create(
            Named("b") = params.tail(X.cols()),
            Named("ssq_b_1") = NA_REAL,
            Named("converged") = converged
        );
    }

    MatrixXd cov = lu.inverse();
    int n_alpha = model.num_categories() - 1;
    double ssq_b_1 = (X.cols() >= 1) ? cov(n_alpha, n_alpha) : NA_REAL;

    return List::create(
        Named("b") = params.tail(X.cols()),
        Named("ssq_b_1") = ssq_b_1,
        Named("converged") = converged
    );
}

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]

//' Parallel Adjacent-Category Logit Randomization Distribution
//'
//' @param y Numeric vector of response values (pre-null-shifted for treated).
//' @param X_covars Matrix of covariates (without intercept or treatment).
//' @param w_mat Integer matrix of permuted treatment assignments (n x nsim).
//' @param delta Null treatment effect (additive shift).
//' @param num_cores Number of OpenMP threads.
//' @return Numeric vector of length nsim with treatment coefficients.
// [[Rcpp::export]]
NumericVector compute_adj_cat_logit_distr_parallel_cpp(
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

		AdjacentCategoryLogitRegression model(X_full, y_shifted);
		if (model.num_categories() < 2) continue;

		bool converged = false;
		Eigen::VectorXd params = adjacent_category_newton_fit(model, 100, 1e-8, &converged);

		int n_alpha = model.num_categories() - 1;
		// accept result even if not formally converged, matching R generate_mod behaviour
		if ((int)params.size() >= n_alpha + 1 && std::isfinite(params[n_alpha])) {
			results[b] = params[n_alpha];
		}
	}

	return wrap(results);
}
