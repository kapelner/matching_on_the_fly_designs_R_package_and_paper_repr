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
      \item Automated missing data pipeline using \pkg{missRanger} and \pkg{missForest}.
      \item Full support for unequal treatment allocation probabilities (prob_T != 0.5).
    }
  }
  \subsection{PERFORMANCE AND OPTIMIZATION}{
    \itemize{
      \item Extreme speedups via C++/Eigen/OpenMP implementations of bootstrap and randomization distributions.
      \item Custom high-speed C++ IRLS solvers for Logistic, Beta, and Negative Binomial regression.
      \item Intelligent thread-control for parallelized inference.
    }
  }
  \subsection{ROBUSTNESS AND USABILITY}{
    \itemize{
      \item Standardized global indentation to tabs and modernized Roxygen documentation.
      \item Hidden internal abstract base classes for a streamlined public API.
      \item Robust error handling for small-sample degeneracy and non-convergence.
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
      \item Support for standard designs: CRD, BCRD, Atkinson, KK14, and KK21.
    }
  }
}
