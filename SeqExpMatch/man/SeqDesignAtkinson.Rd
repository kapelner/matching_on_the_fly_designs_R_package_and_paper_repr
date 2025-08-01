% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/design_Atkinson.R
\name{SeqDesignAtkinson}
\alias{SeqDesignAtkinson}
\title{Atkinson's (1982) Covariate-Adjusted Biased Coin Sequential Design}
\description{
An R6 Class encapsulating the data and functionality for a sequential experimental design.
This class takes care of data intialization and sequential assignments. The class object
should be saved securely after each assignment e.g. on an encrypted cloud server.
}
\examples{

## ------------------------------------------------
## Method `SeqDesignAtkinson$new`
## ------------------------------------------------

seq_des = SeqDesign$new(response_type = "continuous")
 
}
\section{Super class}{
\code{SeqExpMatch::SeqDesign} -> \code{SeqDesignAtkinson}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-SeqDesignAtkinson-new}{\code{SeqDesignAtkinson$new()}}
\item \href{#method-SeqDesignAtkinson-clone}{\code{SeqDesignAtkinson$clone()}}
}
}
\if{html}{\out{
<details><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesign" data-id="add_all_subject_responses"><a href='../../SeqExpMatch/html/SeqDesign.html#method-SeqDesign-add_all_subject_responses'><code>SeqExpMatch::SeqDesign$add_all_subject_responses()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesign" data-id="add_subject_response"><a href='../../SeqExpMatch/html/SeqDesign.html#method-SeqDesign-add_subject_response'><code>SeqExpMatch::SeqDesign$add_subject_response()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesign" data-id="add_subject_to_experiment_and_assign"><a href='../../SeqExpMatch/html/SeqDesign.html#method-SeqDesign-add_subject_to_experiment_and_assign'><code>SeqExpMatch::SeqDesign$add_subject_to_experiment_and_assign()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesign" data-id="any_censoring"><a href='../../SeqExpMatch/html/SeqDesign.html#method-SeqDesign-any_censoring'><code>SeqExpMatch::SeqDesign$any_censoring()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesign" data-id="assert_even_allocation"><a href='../../SeqExpMatch/html/SeqDesign.html#method-SeqDesign-assert_even_allocation'><code>SeqExpMatch::SeqDesign$assert_even_allocation()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesign" data-id="assert_experiment_completed"><a href='../../SeqExpMatch/html/SeqDesign.html#method-SeqDesign-assert_experiment_completed'><code>SeqExpMatch::SeqDesign$assert_experiment_completed()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesign" data-id="assert_fixed_sample"><a href='../../SeqExpMatch/html/SeqDesign.html#method-SeqDesign-assert_fixed_sample'><code>SeqExpMatch::SeqDesign$assert_fixed_sample()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesign" data-id="check_experiment_completed"><a href='../../SeqExpMatch/html/SeqDesign.html#method-SeqDesign-check_experiment_completed'><code>SeqExpMatch::SeqDesign$check_experiment_completed()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesign" data-id="get_X_imp"><a href='../../SeqExpMatch/html/SeqDesign.html#method-SeqDesign-get_X_imp'><code>SeqExpMatch::SeqDesign$get_X_imp()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesign" data-id="get_X_raw"><a href='../../SeqExpMatch/html/SeqDesign.html#method-SeqDesign-get_X_raw'><code>SeqExpMatch::SeqDesign$get_X_raw()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesign" data-id="get_covariate_weights"><a href='../../SeqExpMatch/html/SeqDesign.html#method-SeqDesign-get_covariate_weights'><code>SeqExpMatch::SeqDesign$get_covariate_weights()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesign" data-id="get_dead"><a href='../../SeqExpMatch/html/SeqDesign.html#method-SeqDesign-get_dead'><code>SeqExpMatch::SeqDesign$get_dead()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesign" data-id="get_design"><a href='../../SeqExpMatch/html/SeqDesign.html#method-SeqDesign-get_design'><code>SeqExpMatch::SeqDesign$get_design()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesign" data-id="get_n"><a href='../../SeqExpMatch/html/SeqDesign.html#method-SeqDesign-get_n'><code>SeqExpMatch::SeqDesign$get_n()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesign" data-id="get_prob_T"><a href='../../SeqExpMatch/html/SeqDesign.html#method-SeqDesign-get_prob_T'><code>SeqExpMatch::SeqDesign$get_prob_T()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesign" data-id="get_response_type"><a href='../../SeqExpMatch/html/SeqDesign.html#method-SeqDesign-get_response_type'><code>SeqExpMatch::SeqDesign$get_response_type()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesign" data-id="get_t"><a href='../../SeqExpMatch/html/SeqDesign.html#method-SeqDesign-get_t'><code>SeqExpMatch::SeqDesign$get_t()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesign" data-id="get_w"><a href='../../SeqExpMatch/html/SeqDesign.html#method-SeqDesign-get_w'><code>SeqExpMatch::SeqDesign$get_w()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesign" data-id="get_y"><a href='../../SeqExpMatch/html/SeqDesign.html#method-SeqDesign-get_y'><code>SeqExpMatch::SeqDesign$get_y()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesign" data-id="is_fixed_sample_size"><a href='../../SeqExpMatch/html/SeqDesign.html#method-SeqDesign-is_fixed_sample_size'><code>SeqExpMatch::SeqDesign$is_fixed_sample_size()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeqExpMatch" data-topic="SeqDesign" data-id="print_current_subject_assignment"><a href='../../SeqExpMatch/html/SeqDesign.html#method-SeqDesign-print_current_subject_assignment'><code>SeqExpMatch::SeqDesign$print_current_subject_assignment()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SeqDesignAtkinson-new"></a>}}
\if{latex}{\out{\hypertarget{method-SeqDesignAtkinson-new}{}}}
\subsection{Method \code{new()}}{
Initialize a sequential experimental design
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SeqDesignAtkinson$new(
  response_type,
  prob_T = 0.5,
  include_is_missing_as_a_new_feature = TRUE,
  verbose = TRUE,
  n = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{response_type}}{The data type of response values which must be one of the following: 
"continuous", 
"incidence", 
"proportion", 
"count", 
"survival".
This package will enforce that all added responses via the \code{add_subject_response} method will be
of the appropriate type.}

\item{\code{prob_T}}{The probability of the treatment assignment. This defaults to \code{0.5}.}

\item{\code{include_is_missing_as_a_new_feature}}{If missing data is present in a variable, should we include another dummy variable for its
missingness in addition to imputing its value? If the feature is type factor, instead of creating
a new column, we allow missingness to be its own level. The default is \code{TRUE}.}

\item{\code{verbose}}{A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}.}

\item{\code{n}}{The sample size (if fixed). Default is \code{NULL} for not fixed.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new `SeqDesign` object of the specific type
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{seq_des = SeqDesign$new(response_type = "continuous")
 
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SeqDesignAtkinson-clone"></a>}}
\if{latex}{\out{\hypertarget{method-SeqDesignAtkinson-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SeqDesignAtkinson$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
