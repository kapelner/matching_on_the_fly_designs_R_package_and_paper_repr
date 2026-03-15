\name{SeqDesignInferenceOrdinalUnivCondContRatioRegr}
\alias{SeqDesignInferenceOrdinalUnivCondContRatioRegr}
\title{Univariate Conditional Continuation-Ratio Inference}
\usage{
SeqDesignInferenceOrdinalUnivCondContRatioRegr$new(
  seq_des_obj,
  num_cores = 1,
  verbose = FALSE
)
}
\arguments{
\item{seq_des_obj}{A SeqDesign object with an ordinal response.}
\item{num_cores}{Number of CPU cores.}
\item{verbose}{Whether to print progress messages.}
}
\description{
Fits a conditional continuation-ratio logit model for ordinal responses by
expanding the response into continuation trials and applying conditional
logistic regression across strata induced by matched sets or singleton
reservoir units.
}
