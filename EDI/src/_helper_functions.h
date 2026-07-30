#ifndef EDI_HELPERS_H
#define EDI_HELPERS_H

// Thin Rcpp-glue wrapper around _helper_functions_core.h (portable, no Rcpp/R
// dependency). This file adds back just the handful of functions that
// inherently need Rcpp types (Rcpp::Nullable arguments, Rcpp::List return
// values) -- everything else lives in the core header. See
// package_metadata/python_bindings_package_spec.md's "Math-Utility Function
// Exports"/R-dependency-audit sections and perf_experiments_final.md for the
// rationale; this split lets a future Python/pybind11 build include
// _helper_functions_core.h directly with zero Rcpp/R headers on the include
// path, while every existing R-facing .cpp file keeps including this exact
// filename unchanged.

// R's CXXFLAGS includes -UNDEBUG which overrides PKG_CXXFLAGS's -DNDEBUG.
// Source-level #define takes effect after command-line flags, so re-defining
// NDEBUG HERE -- before RcppEigen.h/Eigen are first included anywhere in this
// translation unit -- disables C assert() for Eigen/LBFGSpp headers. Without
// this, Eigen's BDCSVD and LBFGSpp line search both fire abort() on numerical
// edge cases, crashing R. _helper_functions.h is the first include in nearly
// every .cpp file in this codebase, so this ordering (before <RcppEigen.h>
// below) is what actually matters; _helper_functions_core.h repeats the same
// guard for when it's included directly without going through this file.
#ifndef NDEBUG
#define NDEBUG
#endif

#include <RcppEigen.h>
#include "_helper_functions_core.h"

using Rcpp::as;
using Rcpp::List;
using Rcpp::Named;
using Rcpp::stop;

// Converts an Rcpp::Nullable<T> (the required R-facing type for an optional
// exported-function argument) into a std::optional<NativeT> so all the
// downstream logic in a function body can use std::optional idioms
// (has_value()/*opt) instead of Rcpp's isNotNull()/as<T>(). Keep the call at
// the top of a function, once per parameter — this is the only place
// Rcpp::Nullable should appear once a function has been migrated.
template <typename NativeT, typename RcppT>
inline std::optional<NativeT> nullable_to_optional(const Rcpp::Nullable<RcppT>& x) {
    if (x.isNull()) return std::nullopt;
    return Rcpp::as<NativeT>(x);
}

inline FixedParamSpec make_fixed_param_spec(
    int n_params,
    const Rcpp::Nullable<Rcpp::IntegerVector>& fixed_idx,
    const Rcpp::Nullable<Rcpp::NumericVector>& fixed_values
) {
    return make_fixed_param_spec(
        n_params,
        nullable_to_optional<Eigen::VectorXi>(fixed_idx),
        nullable_to_optional<Eigen::VectorXd>(fixed_values)
    );
}

inline FixedParamSpec make_fixed_param_spec(
    int n_params,
    SEXP fixed_idx,
    SEXP fixed_values
) {
    if (Rf_isNull(fixed_idx)) {
        return make_fixed_param_spec(n_params);
    }
    return make_fixed_param_spec(
        n_params,
        Rcpp::as<Eigen::VectorXi>(fixed_idx),
        Rf_isNull(fixed_values) ? std::optional<Eigen::VectorXd>() : std::optional<Eigen::VectorXd>(Rcpp::as<Eigen::VectorXd>(fixed_values))
    );
}

inline Rcpp::List make_uniform_likelihood_fit_result(const Eigen::VectorXd& params,
                                                     double neg_loglik,
                                                     bool converged,
                                                     const Eigen::VectorXd& score,
                                                     const Eigen::MatrixXd& observed_information,
                                                     bool estimate_only = false,
                                                     const Eigen::MatrixXd* fisher_information = nullptr,
                                                     const std::string& information_type = "observed") {
    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("params") = params,
        Rcpp::Named("neg_loglik") = neg_loglik,
        Rcpp::Named("neg_ll") = neg_loglik,
        Rcpp::Named("loglik") = R_finite(neg_loglik) ? -neg_loglik : NA_REAL,
        Rcpp::Named("converged") = converged
    );
    if (!estimate_only) {
        out["score"] = score;
        out["observed_information"] = observed_information;
        out["hessian"] = -observed_information;
        if (fisher_information != nullptr) {
            out["fisher_information"] = *fisher_information;
            out["information"] = *fisher_information;
            out["information_type"] = "fisher";
            out["vcov"] = covariance_from_information(*fisher_information);
        } else {
            out["information"] = observed_information;
            out["information_type"] = information_type;
            out["vcov"] = covariance_from_information(observed_information);
        }
    }
    return out;
}

inline Rcpp::List likelihood_ratio_test_from_negloglik(double unrestricted_neg_loglik,
                                                       double null_neg_loglik,
                                                       int df = 1) {
    double statistic = NA_REAL;
    double p_value = NA_REAL;
    if (R_finite(unrestricted_neg_loglik) && R_finite(null_neg_loglik) && df > 0) {
        statistic = std::max(0.0, 2.0 * (null_neg_loglik - unrestricted_neg_loglik));
        p_value = fast_pchisq_upper(statistic, static_cast<double>(df));
    }
    return Rcpp::List::create(
        Rcpp::Named("statistic") = statistic,
        Rcpp::Named("df") = df,
        Rcpp::Named("p_value") = p_value
    );
}

inline Rcpp::List score_test_from_score_information(const Eigen::VectorXd& score,
                                                    const Eigen::MatrixXd& information,
                                                    int tested_idx) {
    const int idx = tested_idx - 1;
    if (idx < 0 || idx >= score.size() || information.rows() != score.size() || information.cols() != score.size()) {
        Rcpp::stop("tested_idx must be a one-based index within the parameter vector");
    }

    std::vector<int> nuisance_idx_v;
    nuisance_idx_v.reserve(score.size() - 1);
    for (int i = 0; i < score.size(); ++i) {
        if (i != idx) nuisance_idx_v.push_back(i);
    }

    double info_eff = information(idx, idx);
    if (!information.allFinite()) {
        return Rcpp::List::create(
            Rcpp::Named("statistic") = NA_REAL,
            Rcpp::Named("df") = 1,
            Rcpp::Named("p_value") = NA_REAL,
            Rcpp::Named("score") = NA_REAL,
            Rcpp::Named("information_effective") = NA_REAL
        );
    }
    if (!nuisance_idx_v.empty()) {
        Eigen::VectorXi nuisance_idx(nuisance_idx_v.size());
        for (int i = 0; i < static_cast<int>(nuisance_idx_v.size()); ++i) nuisance_idx[i] = nuisance_idx_v[i];
        Eigen::MatrixXd I_nn = subset_matrix(information, nuisance_idx, nuisance_idx);
        Eigen::VectorXd I_nT = subset_matrix(information, nuisance_idx, Eigen::VectorXi::Constant(1, idx)).col(0);
        Eigen::VectorXd solved = covariance_from_information(I_nn) * I_nT;
        info_eff -= I_nT.dot(solved);
    }

    double statistic = NA_REAL;
    double p_value = NA_REAL;
    const double score_t = score[idx];
    if (R_finite(score_t) && R_finite(info_eff) && info_eff > 0.0) {
        statistic = score_t * score_t / info_eff;
        p_value = fast_pchisq_upper(statistic, 1.0);
    }

    return Rcpp::List::create(
        Rcpp::Named("statistic") = statistic,
        Rcpp::Named("df") = 1,
        Rcpp::Named("p_value") = p_value,
        Rcpp::Named("score") = score_t,
        Rcpp::Named("information_effective") = info_eff
    );
}

inline Rcpp::List gradient_test_from_restricted_score(const Eigen::VectorXd& score,
                                                      double unrestricted_estimate,
                                                      double null_value,
                                                      int tested_idx) {
    const int idx = tested_idx - 1;
    if (idx < 0 || idx >= score.size()) {
        Rcpp::stop("tested_idx must be a one-based index within the parameter vector");
    }

    double statistic = NA_REAL;
    double p_value = NA_REAL;
    const double score_t = score[idx];
    const double estimate_gap = unrestricted_estimate - null_value;

    if (R_finite(score_t) && R_finite(estimate_gap)) {
        statistic = score_t * estimate_gap;
        if (statistic < 0.0) {
            statistic = 0.0;
        }
        if (R_finite(statistic)) {
            p_value = fast_pchisq_upper(statistic, 1.0);
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("statistic") = statistic,
        Rcpp::Named("df") = 1,
        Rcpp::Named("p_value") = p_value,
        Rcpp::Named("score") = score_t,
        Rcpp::Named("estimate_gap") = estimate_gap
    );
}

#endif
