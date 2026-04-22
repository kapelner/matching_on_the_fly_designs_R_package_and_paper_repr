
# Ensure we are in the project root if run from package_tests
if (basename(getwd()) == "package_tests") {
  setwd("..")
}

devtools::load_all("EDI")
ns = getNamespace("EDI")
all_objs = ls(ns)
inference_classes = all_objs[grep("^Inference", all_objs)]

# Filter for R6 classes
is_r6 = sapply(inference_classes, function(x) {
  obj = get(x, envir = ns)
  inherits(obj, "R6ClassGenerator")
})
inference_classes = inference_classes[is_r6]

# Base/Abstract classes to skip
base_classes = c(
  "Inference", "InferenceAsymp", "InferenceBoot", "InferenceRand",
  "InferenceRandCI", "InferenceExact", "InferenceKKPassThrough",
  "InferenceKKPassThroughCompound", "InferenceMLEorKMforGLMs",
  "InferenceSuite", "InferenceAbstractQuantileRandCI"
)
# Also skip classes with "Abstract" in name
concrete_classes = setdiff(inference_classes, base_classes)
concrete_classes = concrete_classes[!grepl("Abstract", concrete_classes)]

# Optimization metadata mapping
get_opt_metadata = function(res_type, base_name) {
  # Defaults
  via = "TODO"
  alg = "NA"
  
  if (base_name %in% c("BaiAdjustedT", "BaiAdjustedTKK14", "BaiAdjustedTKK21", "KKCompoundMeanDiff", "SimpleMeanDiff", "SimpleMeanDiffPooledVar", "Wald", "RiskDiff")) {
    via = "none/closed-form statistic"
  } else if (base_name %in% c("KKWilcoxIVWC", "SimpleWilcox", "Ridit", "GehanWilcox", "LogRank")) {
    via = "none/rank statistic"
  } else if (grepl("QuantileRegr", base_name)) {
    via = "quantreg"
  } else if (base_name == "KKGLMM") {
    via = "ourself/Rcpp GLMM"
    alg = "L-BFGS"
  } else if (base_name == "MultGLS") {
    via = "nlme"
  } else if (grepl("Exact", base_name)) {
    via = "none/exact enumeration"
  } else if (base_name %in% c("Azriel", "ExtendedRobins")) {
    via = "none/design-specific statistic"
  } else if (base_name == "KMDiff" || base_name == "RestrictedMeanDiff") {
    via = "none/KM statistic"
  } else if (base_name == "PairedSignTest") {
    via = "none/sign statistic"
  } else if (base_name == "JonckheereTerpstraTest") {
    via = "none/Jonckheere-Terpstra statistic"
  } else if (base_name == "Lin") {
    via = "ourself/closed-form linear estimator"
    alg = "closed-form QR solve"
  } else if (base_name == "OLS") {
    via = "ourself/Rcpp OLS"
    alg = "closed-form QR solve"
  } else if (grepl("RobustRegr", base_name)) {
    via = "ourself/Rcpp robust"
    alg = "IRLS"
  } else if (base_name == "LogRegr") {
    via = "ourself/Rcpp logistic"
    alg = "IRLS"
  } else if (base_name == "ModifiedPoisson" || base_name == "KKModifiedPoisson") {
    via = "ourself/Rcpp Poisson working likelihood"
    alg = "IRLS"
  } else if (base_name == "Poisson" || base_name == "RobustPoisson") {
    via = "ourself/Rcpp Poisson"
    alg = "IRLS"
  } else if (base_name == "QuasiPoisson") {
    via = "ourself/Rcpp quasi-Poisson"
    alg = "IRLS"
  } else if (base_name == "NegBin" || base_name == "ZeroInflatedNegBin") {
    via = "ourself/Rcpp negative binomial"
    alg = "L-BFGS"
  } else if (base_name == "HurdleNegBin") {
    via = "ourself/Rcpp hurdle negative binomial"
    alg = "IRLS plus L-BFGS"
  } else if (base_name %in% c("HurdlePoisson", "ZeroInflatedPoisson", "ZeroInflatedNegBin") || grepl("KKHurdlePoisson", base_name)) {
    via = "glmmTMB"
  } else if (base_name == "KKCPoissonIVWC") {
    via = "ourself/Rcpp conditional-Poisson plus negative binomial"
    alg = "weighted logistic IRLS plus L-BFGS"
  } else if (base_name == "KKCPoissonOneLik") {
    via = "ourself/Rcpp conditional-Poisson combined likelihood"
    alg = "Newton-Raphson"
  } else if (grepl("KKClogitPlusGLMM", base_name)) {
    via = "ourself/Rcpp conditional-logistic plus GLMM"
    alg = "L-BFGS"
  } else if (grepl("KKClogit", base_name)) {
    via = "ourself/Rcpp conditional-logistic"
    alg = "L-BFGS"
  } else if (grepl("KKGEE", base_name)) {
    via = "geepack"
  } else if (grepl("Newcombe", base_name) || grepl("Miettinen", base_name)) {
    via = "none/score or interval equations"
  } else if (base_name %in% c("LogBinomial", "BinomialIdentityRiskDiff")) {
    via = "ourself/Rcpp identity-binomial"
    alg = "constrained Fisher scoring with step-halving"
  } else if (base_name == "BetaRegr") {
    via = "ourself/Rcpp beta regression"
    alg = "L-BFGS"
  } else if (base_name == "FractionalLogit") {
    via = "ourself/Rcpp quasi-binomial-style objective"
    alg = "IRLS"
  } else if (base_name == "ZeroOneInflatedBetaRegr") {
    via = "ourself/Rcpp zero-one-inflated beta"
    alg = "L-BFGS"
  } else if (base_name == "CoxPHRegr" || grepl("LWACox", base_name) || grepl("StratCox", base_name)) {
    via = "ourself/Rcpp Cox partial likelihood"
    alg = "Newton-Raphson"
  } else if (base_name == "DepCensTransformRegr") {
    via = "ourself/Rcpp dependent-censoring transformation"
    alg = "L-BFGS"
  } else if (grepl("ClaytonCopula", base_name)) {
    via = "ourself/Rcpp Clayton copula Weibull AFT"
    alg = "L-BFGS"
  } else if (base_name == "WeibullRegr") {
    via = "ourself/Rcpp Weibull AFT"
    alg = "L-BFGS"
  } else if (grepl("WeibullFrailty", base_name)) {
    via = "survival"
  } else if (grepl("RankRegr", base_name)) {
    via = "none/rank-regression path"
  } else if (base_name %in% c("CauchitRegr", "CloglogRegr", "OrderedProbitRegr", "PropOddsRegr")) {
    via = "ourself/Rcpp cumulative link"
    alg = "Newton-Raphson with finite-difference derivatives and line search"
  } else if (grepl("PartialProportionalOdds", base_name)) {
    via = "VGAM"
  } else if (grepl("CLMM", base_name)) {
    via = "ordinal"
  } else if (grepl("WeibullFrailty", base_name)) {
    via = "survival"
  } else if (grepl("GComp", base_name)) {
    via = "ourself/nuisance likelihood"
  } else if (base_name == "MLEorKMSummaryTable") {
    via = "none/summary table"
  }
  
  irls_avail = "NA"
  if (grepl("ourself", via)) {
    if (base_name %in% c("LogRegr", "Poisson", "RobustPoisson", "QuasiPoisson", "ModifiedPoisson", "KKModifiedPoisson", "FractionalLogit", "RobustRegr", "OLS", "Lin") || grepl("RobustRegr", base_name)) {
       irls_avail = "Y"
    } else if (base_name == "HurdleNegBin") {
       irls_avail = "Y" # Hurdle models use IRLS for at least one component
    } else if (base_name %in% c("NegBin", "BetaRegr", "ZeroOneInflatedBetaRegr", "CoxPHRegr", "WeibullRegr", "DepCensTransformRegr") || grepl("ClaytonCopula", base_name) || grepl("LWACox", base_name) || grepl("StratCox", base_name) || grepl("Clogit", base_name)) {
       irls_avail = "N"
    } else if (base_name %in% c("CauchitRegr", "CloglogRegr", "OrderedProbitRegr", "PropOddsRegr")) {
       irls_avail = "N" # Usually NR for cumulative link
    } else if (base_name %in% c("LogBinomial", "BinomialIdentityRiskDiff")) {
       irls_avail = "Y" # Technically Fisher scoring is IRLS
    } else {
       irls_avail = "N"
    }
  }
  
  list(optimization_via = via, optim_alg = alg, irls_avail = irls_avail)
}

# Mapping function
get_info = function(class_name) {
  obj_gen = get(class_name, envir = ns)
  
  # Check methods including inheritance
  has_method = function(gen, method_name) {
    curr = gen
    while (!is.null(curr)) {
      if (method_name %in% names(curr$public_methods) || method_name %in% names(curr$private_methods)) return(TRUE)
      
      # Walk up inherit
      if (is.null(curr$inherit)) break
      inherit_obj = curr$inherit
      if (is.symbol(inherit_obj)) {
        name = as.character(inherit_obj)
        if (exists(name, envir = ns)) {
          curr = get(name, envir = ns)
        } else {
          break
        }
      } else if (is.environment(inherit_obj) && exists("classname", envir = inherit_obj)) {
        curr = inherit_obj
      } else {
        break
      }
    }
    FALSE
  }
  
  has_likelihood = has_method(obj_gen, "compute_likelihood") || 
                   has_method(obj_gen, "shared_combined_likelihood") ||
                   has_method(obj_gen, "fit_glmm") ||
                   has_method(obj_gen, "generate_mod")
  has_loglik     = has_method(obj_gen, "compute_log_likelihood") || 
                   has_method(obj_gen, "shared_combined_likelihood") ||
                   has_method(obj_gen, "fit_glmm") ||
                   has_method(obj_gen, "generate_mod")
  has_score      = has_method(obj_gen, "compute_score") || 
                   has_method(obj_gen, "shared_combined_likelihood") ||
                   has_method(obj_gen, "fit_glmm") ||
                   has_method(obj_gen, "generate_mod")
  
  # Check inheritance for Wald
  curr = obj_gen
  uses_wald = FALSE
  while (!is.null(curr$inherit)) {
    inherit_obj = curr$inherit
    if (is.symbol(inherit_obj)) {
       name = as.character(inherit_obj)
       if (exists(name, envir = ns)) {
         inherit_obj = get(name, envir = ns)
       } else {
         break
       }
    }
    if (is.environment(inherit_obj) && exists("classname", envir = inherit_obj)) {
       if (inherit_obj$classname == "InferenceAsymp") {
         uses_wald = TRUE
         break
       }
    } else {
       break
    }
    curr = inherit_obj
  }
  
  # Response type and normalized name
  res_type = "unknown"
  base_name = class_name
  
  if (grepl("^InferenceAll", class_name)) {
    res_type = "all"
    base_name = sub("^InferenceAll", "", class_name)
  } else if (grepl("^InferenceContin", class_name)) {
    res_type = "continuous"
    base_name = sub("^InferenceContin", "", class_name)
  } else if (grepl("^InferenceCount", class_name)) {
    res_type = "count"
    base_name = sub("^InferenceCount", "", class_name)
  } else if (grepl("^InferenceIncidence", class_name)) {
     res_type = "incidence"
     base_name = sub("^InferenceIncidence", "", class_name)
  } else if (grepl("^InferenceIncid", class_name)) {
    res_type = "incidence"
    base_name = sub("^InferenceIncid", "", class_name)
  } else if (grepl("^InferenceOrdinal", class_name)) {
    res_type = "ordinal"
    base_name = sub("^InferenceOrdinal", "", class_name)
  } else if (grepl("^InferenceProp", class_name)) {
    res_type = "proportion"
    base_name = sub("^InferenceProp", "", class_name)
  } else if (grepl("^InferenceSurvival", class_name)) {
    res_type = "survival"
    base_name = sub("^InferenceSurvival", "", class_name)
  } else if (grepl("^InferenceBaiAdjustedT", class_name)) {
    res_type = "all"
    base_name = sub("^Inference", "", class_name)
  } else if (grepl("^InferenceMLEorKMSummaryTable", class_name)) {
    res_type = "all"
    base_name = sub("^Inference", "", class_name)
  }
  
  opt = get_opt_metadata(res_type, base_name)
  
  list(
    response_type = res_type,
    inference_path = base_name,
    has_likelihood = if (has_likelihood) "Y" else "N",
    has_loglik = if (has_loglik) "Y" else "N",
    has_score = if (has_score) "Y" else "N",
    uses_wald_ci = if (uses_wald) "Y" else "N",
    uses_wald_test = if (uses_wald) "Y" else "N",
    optimization_via = opt$optimization_via,
    optim_alg = opt$optim_alg,
    irls_avail = opt$irls_avail
  )
}

results = lapply(concrete_classes, get_info)
df_new = do.call(rbind, lapply(results, as.data.frame))

# Load existing to preserve manually curated columns if possible
existing = read.csv("package_metadata/inference_paths_info.csv")

# Merge strategy: keep manual columns if they are not TODO and not NA
final_df = merge(df_new, existing[, c("response_type", "inference_path", "has_multi")], 
                 by = c("response_type", "inference_path"), all.x = TRUE)

# Fill NA in has_multi (if missing)
final_df$has_multi[is.na(final_df$has_multi)] = "N"

# Reorder columns
final_df = final_df[, c("response_type", "inference_path", "has_multi", "has_likelihood", "has_loglik", "has_score", "uses_wald_ci", "uses_wald_test", "optimization_via", "optim_alg", "irls_avail")]

# Sort
final_df = final_df[order(final_df$response_type, final_df$inference_path), ]

write.csv(final_df, "package_metadata/inference_paths_info.csv", row.names = FALSE, quote = FALSE)
cat("Refreshed package_metadata/inference_paths_info.csv with IRLS availability.\n")
