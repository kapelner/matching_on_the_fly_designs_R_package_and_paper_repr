// Greedy experimental design search — native C++ (RcppEigen + OpenMP).
//
// Algorithm (per design):
//   1. Random balanced init: n/2 treated via Fisher-Yates.
//   2. Criterion vector  d = M * (2w − 1),  M is p×n.
//   3. For n_iter steps: pick random (treated i, control j), compute
//      Δ = 2(M[:,j] − M[:,i]), accept swap if f(d+Δ) < f(d).
//
// Objectives
//   abs_sum_diff : M = X_std'/n,  f(d) = ‖d‖₁  (fused expr, no temp vector)
//   mahal_dist   : M = L⁻¹X'/n   (L = Cholesky of Σ),
//                  f(d+Δ) = f + 2 d·Δ + ‖Δ‖²  (incremental — avoids forming d+Δ)
//
// Parallelism: outer loop over r designs is embarrassingly parallel (OMP).
//   Per-thread std::mt19937 RNGs are seeded from R's RNG before the parallel
//   region so set.seed() governs reproducibility at the R level.
//   Thread count follows whatever omp_set_num_threads() was last called with
//   (EDI sets this via set_num_cores() → set_omp_num_threads_cpp()).

#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppEigen)]]

#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <string>
#include <random>
#include <limits>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVectorXd;
using Eigen::LLT;
using Eigen::Success;

// [[Rcpp::export]]
Rcpp::IntegerMatrix greedy_design_search_cpp(
    const Eigen::Map<Eigen::MatrixXd> X_raw,
    const int                         r,
    const std::string&                objective,
    const int                         n_iter
) {
    const int n  = static_cast<int>(X_raw.rows());
    const int p  = static_cast<int>(X_raw.cols());
    const int nt = n / 2;

    Rcpp::IntegerMatrix result(n, r);

    // ── Seed per-thread RNGs from R's RNG (single-threaded, before parallel) ─
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    std::vector<uint32_t> seeds(static_cast<std::size_t>(nthreads));
    GetRNGstate();
    for (int t = 0; t < nthreads; t++)
        seeds[static_cast<std::size_t>(t)] =
            static_cast<uint32_t>(::unif_rand() * 4294967295.0);
    PutRNGstate();

    // ── No covariates: pure balanced randomisation ───────────────────────────
    if (p == 0) {
#pragma omp parallel
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            std::mt19937 rng(seeds[static_cast<std::size_t>(tid)]);
            std::vector<int> idx(n);
#pragma omp for schedule(static)
            for (int d = 0; d < r; d++) {
                std::iota(idx.begin(), idx.end(), 0);
                for (int i = n - 1; i > 0; i--)
                    std::swap(idx[static_cast<std::size_t>(i)],
                              idx[static_cast<std::size_t>(
                                  std::uniform_int_distribution<int>(0, i)(rng))]);
                for (int i = 0; i < nt; i++)
                    result(idx[static_cast<std::size_t>(i)], d) = 1;
            }
        }
        return result;
    }

    // ── Build M (p×n): d = M*(2w−1) encodes the balance criterion ───────────
    RowVectorXd col_means = X_raw.colwise().mean();
    MatrixXd    X         = X_raw.rowwise() - col_means;    // centred, n×p

    MatrixXd M(p, n);
    bool abs_mode = true;

    if (objective == "mahal_dist") {
        MatrixXd cov_mat = (X.transpose() * X) / std::max(1, n - 1);
        LLT<MatrixXd> llt(cov_mat);
        if (llt.info() == Success) {
            // d = L⁻¹ X'(2w−1)/n  →  f = ‖d‖²   (Σ = LL')
            M        = llt.matrixL().solve(X.transpose()) / static_cast<double>(n);
            abs_mode = false;
        }
        // singular covariance → fall through to abs_sum_diff
    }

    if (abs_mode) {
        VectorXd inv_sd(p);
        for (int j = 0; j < p; j++) {
            double var = X.col(j).squaredNorm() / std::max(1, n - 1);
            inv_sd[j]  = (var < 1e-24) ? 1.0 : 1.0 / std::sqrt(var);
        }
        M = (X * inv_sd.asDiagonal()).transpose() / static_cast<double>(n);
    }

    // ── Parallel greedy search over r independent designs ────────────────────
#pragma omp parallel
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        std::mt19937                     rng(seeds[static_cast<std::size_t>(tid)]);
        std::uniform_int_distribution<int> pick_nt(0, nt - 1);

        // Thread-local working buffers
        std::vector<int> w(n), treated(nt), control(nt), order(n);
        VectorXd dbl_w(n), d(p), delta(p);

#pragma omp for schedule(static)
        for (int des = 0; des < r; des++) {
            // ── Random balanced init (Fisher-Yates) ──────────────────────────
            std::iota(order.begin(), order.end(), 0);
            for (int i = n - 1; i > 0; i--)
                std::swap(order[static_cast<std::size_t>(i)],
                          order[static_cast<std::size_t>(
                              std::uniform_int_distribution<int>(0, i)(rng))]);
            std::fill(w.begin(), w.end(), 0);
            for (int i = 0; i < nt; i++) w[static_cast<std::size_t>(order[static_cast<std::size_t>(i)])] = 1;

            { int ti = 0, ci = 0;
              for (int i = 0; i < n; i++) {
                  if (w[static_cast<std::size_t>(i)]) treated[static_cast<std::size_t>(ti++)] = i;
                  else                                 control[static_cast<std::size_t>(ci++)] = i;
              } }

            for (int i = 0; i < n; i++) dbl_w[i] = w[static_cast<std::size_t>(i)] ? 1.0 : -1.0;
            d.noalias() = M * dbl_w;
            double f    = abs_mode ? d.lpNorm<1>() : d.squaredNorm();

            // ── Greedy swap iterations ───────────────────────────────────────
            for (int it = 0; it < n_iter; it++) {
                const int ti = pick_nt(rng);
                const int ci = pick_nt(rng);
                const int i  = treated[static_cast<std::size_t>(ti)];
                const int j  = control[static_cast<std::size_t>(ci)];

                // Δ = 2(M[:,j] − M[:,i])  — single fused SIMD pass
                delta.noalias() = 2.0 * (M.col(j) - M.col(i));

                double f_cand;
                if (abs_mode) {
                    // Fused: abs_mode avoids writing d+Δ; Eigen evaluates as
                    // a single vectorised pass over (d[k]+delta[k]).
                    f_cand = (d + delta).lpNorm<1>();
                } else {
                    // ‖d+Δ‖² = ‖d‖² + 2 d·Δ + ‖Δ‖²  — no temp vector needed
                    f_cand = f + 2.0 * d.dot(delta) + delta.squaredNorm();
                }

                if (f_cand < f) {
                    w[static_cast<std::size_t>(i)] = 0;
                    w[static_cast<std::size_t>(j)] = 1;
                    treated[static_cast<std::size_t>(ti)] = j;
                    control[static_cast<std::size_t>(ci)] = i;
                    d += delta;   // in-place — no copy of d
                    f  = f_cand;
                }
            }

            for (int i = 0; i < n; i++) result(i, des) = w[static_cast<std::size_t>(i)];
        }
    } // end omp parallel

    return result;
}
