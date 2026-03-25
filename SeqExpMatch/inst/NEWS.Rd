\name{NEWS}
\title{News for Package "SeqExpMatch"}

\section{Changes in SeqExpMatch version 2.0.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Complete R6 class inheritance revolution for all inference methodologies.
      \item Implementation of the dual-strategy framework for designs with matched-pair strata:
        \itemize{
          \item \strong{IVWC}: Rigorous inverse-variance weighted combination of component estimators.
          \item \strong{Combined Likelihood}: Joint design-aware modeling of matched and unmatched data.
        }
      \item Support for Zhang (2026) Exact Inference for distribution-free p-values and confidence intervals (Incidence and Quantile Regression).
      \item Standardized all survival-based methodologies to Accelerated Failure Time (AFT) log-time ratios.
      \item Integrated modern statistical engines: \pkg{glmmTMB} (GLMM), \pkg{fixest} (Fast GLMs), \pkg{aftgee} (Rank AFT), \pkg{parfm} (Frailty), and \pkg{Rfit} (R-estimation).
      \item Support for Stratified Cox and Conditional Poisson regression.
      \item Added LWA-style marginal Cox and Clayton-copula paired-survival inference for KK designs under both IVWC and Combined Likelihood pooling, plus transformation-model survival inference with dependent censoring.
      \item Standalone log-rank inference powered by a new C++ score/variance routine plus automatic covariate-based strata construction.
      \item Added Lin-style adjusted estimators, robust regression, and quantile regression support for continuous outcomes across non-KK and KK designs.
      \item Added new incidence-response families: modified Poisson, marginal standardization / g-computation for risk differences and risk ratios, Miettinen-Nurminen and Newcombe risk-difference inference, log-binomial regression, and binomial identity-link regression.
      \item Added new count-response families for KK designs: hurdle-Poisson inference under both IVWC and Combined Likelihood pooling.
      \item Added new proportion-response families: fractional-logit and zero/one-inflated beta regression.
      \item Randomization confidence intervals now support ordinal GLM estimators (e.g., cumulative logit/probit) by inverting the exact permutation test for their shared treatment effect, following the Wang & Rosenberger inversion approach.
      \item Expanded ordinal support with cumulative probit, adjacent-category logit, stereotype logit/probit, ridit, paired sign, Jonckheere-Terpstra, partial proportional odds, conditional ordinal KK estimators, and KK ordinal mixed models with both logit and probit links.
      \item Zero-inflated/negative binomial and hurdle negative binomial regression classes plus Poisson/quasi-Poisson inference that now use high-speed Rcpp IRLS solvers with Huber-White variance.
      \item Automated missing data pipeline using \pkg{missRanger} and \pkg{missForest}.
      \item Full support for unequal treatment allocation probabilities (prob_T != 0.5).
    }
  }
  \subsection{PERFORMANCE AND OPTIMIZATION}{
    \itemize{
      \item Extreme speedups via C++/Eigen/OpenMP implementations of bootstrap and randomization distributions.
      \item Custom high-speed C++ IRLS solvers for Logistic, Beta, Negative Binomial, log-binomial, and identity-link binomial regression.
      \item Upgraded the Negative Binomial IRLS solver to use the L-BFGS algorithm for superior convergence stability.
      \item Enhanced \code{fast_beta_regression_with_var} to return the precision parameter \code{phi}.
      \item Added C++ helpers for ordinal g-computation, Jonckheere-Terpstra exact calculations, ridit analysis, and KK grouping / cluster setup for newer pooled models.
      \item Rcpp helper for computing informative survival strata IDs and fast log-rank statistics.
      \item Intelligent thread-control for parallelized inference.
      \item Implemented a specialized randomization-draw optimizer for Clayton Copula models that reuses initial MLE estimates and design configurations, achieving significant speedups in exact inference and confidence interval inversion.
      \item Optimized Weibull frailty models (Combined Likelihood and IVWC) by implementing a high-speed approximation for resampling draws, reducing bootstrap and randomization CI times by over 95\% for univariate settings.
      \item Accelerated GLMM and GEE randomization inference by bypassing expensive summary table computations during permutations.
      \item Implemented high-speed inference object duplication using \code{self$clone()} to bypass constructor overhead during large-scale resampling runs.
    }
  }
  \subsection{ROBUSTNESS, USABILITY, AND BUG FIXES}{
    \itemize{
      \item Standardized global indentation to tabs and modernized Roxygen documentation.
      \item Hidden internal abstract base classes for a streamlined public API.
      \item Robust error handling for small-sample degeneracy and non-convergence.
      \item Improved randomization-test caching so treatment-assignment draws are reused across repeated p-value and confidence-interval inversion calls, including custom randomization statistics.
      \item Cleaned up Rd generation issues and removed links to internal abstract classes from compiled help pages.
      \item Fixed a critical "attempt to apply non-function" error in \code{DesignInferenceIncidUnivKKNewcombeRiskDiff} by implementing the missing \code{pool_estimates_ivwc} method.
      \item Resolved \code{num_cores} assertion failures during \code{SeqDesign} initialization across all subclasses, ensuring correct handling of parallelization parameters.
      \item Corrected Rd documentation line-width issues by consolidating method examples and applying strategic line breaks in R6 documentation.
      \item Implemented missing C++ parallelization routines for Ridit analysis bootstrap and randomization distributions.
      \item Fixed broken test cases related to argument permutations and output structure in fast GLM solvers.
      \item Standardized the \code{DESCRIPTION} file with proper Author/Maintainer fields and fixed DCF formatting for CRAN compliance.
      \item Added automated \code{src} directory cleanup to the build process.
    }
  }
}

\section{Changes in SeqExpMatch version 1.0.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Support for general response types: incidence, proportion, count, and survival.
      \item Support for adding responses whenever they become recorded.
      \item Support for general data frames including factors, characters, and logicals.
      \item Integrated \pkg{missRanger} and \pkg{missForest} for missing data.
    }
  }
}

\section{Changes in SeqExpMatch version 0.1.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Initial Release.
      \item Support for continuous response types.
      \item Support for standard designs: Bernoulli, BCRD, Atkinson, KK14, and KK21.
    }
  }
}
