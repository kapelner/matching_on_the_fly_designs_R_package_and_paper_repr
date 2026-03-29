#' EDI
#'
#' Provides comprehensive support for many fixed and sequential experimental designs and many
#' infererential methods (parametric, nonparametric, exact) for response types continuous, 
#' incidence, count, proportion, survival  (with censoring) and ordinal. Supports automatic 
#' missing data imputation, parallelization, provides robustness fallbacks and is
#' optimized with C++.
#'
#' @name 		EDI
#' @title 		Experimental Design and Inference
#' @author 		Adam Kapelner \email{kapelner@@qc.cuny.edu}
#' @references         Adam Kapelner and Abba Krieger A Matching Procedure for Sequential
#'   Experiments that Iteratively Learns which Covariates Improve Power, Arxiv
#'   2010.05980
#' @keywords	design pair-matching
#' @import      checkmate
#' @import      data.table
#' @import      R6
#' @import      missForest
#' @import      missRanger
#' @import      Rcpp
#' @import      methods
#' @importFrom     stats coef cor glm.control make.link nlminb pnorm pt rnorm qnorm qf rbinom
#' var median sd model.matrix as.formula formula qt pchisq binom.test fisher.test
#' @importFrom	survival Surv survreg.control survfit survdiff
#' @importFrom	utils packageVersion
#' @importFrom	graphics hist
#' @useDynLib   EDI, .registration=TRUE
"_PACKAGE"
