% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inference_continuous_multivariate_ols.R
\name{SeqDesignInferenceContinMultOLS}
\alias{SeqDesignInferenceContinMultOLS}
\title{Simple Mean Difference Inference based on Maximum Likelihood}
\description{
The methods that support confidence intervals and testing for the mean difference
in all response types (except Weibull with censoring) sequential experimental design estimation and test object after the sequential design is completed.
}
\examples{

## ------------------------------------------------
## Method `SeqDesignInferenceContinMultOLS$compute_mle_confidence_interval`
## ------------------------------------------------

seq_des = SeqDesign$new(n = 6, p = 10, design = "CRD")
seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))

seq_des_inf = SeqDesignInferenceContMultOLS$new(seq_des, test_type = "MLE-or-KM-based")
seq_des_inf$compute_confidence_interval()
	
}
\section{Super classes}{
\code{SeqExpMatch::SeqDesignInference} -> \code{SeqExpMatch::SeqDesignInferenceMLEorKM} -> \code{SeqExpMatch::SeqDesignInferenceMLEorKMSummaryTable} -> \code{SeqDesignInferenceContinMultOLS}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-SeqDesignInferenceContinMultOLS-new}{\code{SeqDesignInferenceContinMultOLS$new()}}
\item \href{#method-SeqDesignInferenceContinMultOLS-compute_mle_confidence_interval}{\code{SeqDesignInferenceContinMultOLS$compute_mle_confidence_interval()}}
\item \href{#method-SeqDesignInferenceContinMultOLS-clone}{\code{SeqDesignInferenceContinMultOLS$clone()}}
}
}
\if{html}{\out{
<details><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesignInference" data-id="compute_beta_hat_T_randomization_distr_under_sharp_null"><a href='../../SeqExpMatch/html/SeqDesignInference.html#method-SeqDesignInference-compute_beta_hat_T_randomization_distr_under_sharp_null'><code>SeqExpMatch::SeqDesignInference$compute_beta_hat_T_randomization_distr_under_sharp_null()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesignInference" data-id="compute_confidence_interval_rand"><a href='../../SeqExpMatch/html/SeqDesignInference.html#method-SeqDesignInference-compute_confidence_interval_rand'><code>SeqExpMatch::SeqDesignInference$compute_confidence_interval_rand()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesignInference" data-id="compute_two_sided_pval_for_treatment_effect_rand"><a href='../../SeqExpMatch/html/SeqDesignInference.html#method-SeqDesignInference-compute_two_sided_pval_for_treatment_effect_rand'><code>SeqExpMatch::SeqDesignInference$compute_two_sided_pval_for_treatment_effect_rand()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesignInferenceMLEorKM" data-id="approximate_boostrap_distribution_beta_hat_T"><a href='../../SeqExpMatch/html/SeqDesignInferenceMLEorKM.html#method-SeqDesignInferenceMLEorKM-approximate_boostrap_distribution_beta_hat_T'><code>SeqExpMatch::SeqDesignInferenceMLEorKM$approximate_boostrap_distribution_beta_hat_T()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesignInferenceMLEorKM" data-id="compute_bootstrap_confidence_interval"><a href='../../SeqExpMatch/html/SeqDesignInferenceMLEorKM.html#method-SeqDesignInferenceMLEorKM-compute_bootstrap_confidence_interval'><code>SeqExpMatch::SeqDesignInferenceMLEorKM$compute_bootstrap_confidence_interval()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesignInferenceMLEorKMSummaryTable" data-id="compute_mle_two_sided_pval_for_treatment_effect"><a href='../../SeqExpMatch/html/SeqDesignInferenceMLEorKMSummaryTable.html#method-SeqDesignInferenceMLEorKMSummaryTable-compute_mle_two_sided_pval_for_treatment_effect'><code>SeqExpMatch::SeqDesignInferenceMLEorKMSummaryTable$compute_mle_two_sided_pval_for_treatment_effect()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesignInferenceMLEorKMSummaryTable" data-id="compute_treatment_estimate"><a href='../../SeqExpMatch/html/SeqDesignInferenceMLEorKMSummaryTable.html#method-SeqDesignInferenceMLEorKMSummaryTable-compute_treatment_estimate'><code>SeqExpMatch::SeqDesignInferenceMLEorKMSummaryTable$compute_treatment_estimate()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SeqDesignInferenceContinMultOLS-new"></a>}}
\if{latex}{\out{\hypertarget{method-SeqDesignInferenceContinMultOLS-new}{}}}
\subsection{Method \code{new()}}{
Initialize a sequential experimental design estimation and test object after the sequential design is completed.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SeqDesignInferenceContinMultOLS$new(seq_des_obj, num_cores = 1, verbose = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{seq_des_obj}}{A SeqDesign object whose entire n subjects are assigned and response y is recorded within.}

\item{\code{num_cores}}{The number of CPU cores to use to parallelize the sampling during randomization-based inference 
(which is very slow). The default is 1 for serial computation. This parameter is ignored
for \code{test_type = "MLE-or-KM-based"}.}

\item{\code{verbose}}{A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}

Compute confidence interval}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SeqDesignInferenceContinMultOLS-compute_mle_confidence_interval"></a>}}
\if{latex}{\out{\hypertarget{method-SeqDesignInferenceContinMultOLS-compute_mle_confidence_interval}{}}}
\subsection{Method \code{compute_mle_confidence_interval()}}{
Computes a 1-alpha level frequentist confidence interval differently for all response types, estimate types and test types.

Here we use the theory that MLE's computed for GLM's are asymptotically normal. 
Hence these confidence intervals are asymptotically valid and thus approximate for any sample size.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SeqDesignInferenceContinMultOLS$compute_mle_confidence_interval(alpha = 0.05)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{alpha}}{The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A (1 - alpha)-sized frequentist confidence interval for the treatment effect
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{seq_des = SeqDesign$new(n = 6, p = 10, design = "CRD")
seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))

seq_des_inf = SeqDesignInferenceContMultOLS$new(seq_des, test_type = "MLE-or-KM-based")
seq_des_inf$compute_confidence_interval()
	
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SeqDesignInferenceContinMultOLS-clone"></a>}}
\if{latex}{\out{\hypertarget{method-SeqDesignInferenceContinMultOLS-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SeqDesignInferenceContinMultOLS$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
