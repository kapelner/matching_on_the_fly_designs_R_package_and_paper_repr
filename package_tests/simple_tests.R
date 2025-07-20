#/c/Program\ Files/R/R-devel/bin/R.exe CMD INSTALL -l ~/AppData/Local/R/win-library/4.5/ SeqExpMatch/
rm(list = ls())
pacman::p_load(doParallel, PTE, datasets, qgam, mlbench, AppliedPredictiveModeling, dplyr, ggplot2, gridExtra, profvis, data.table, profvis)
library(SeqExpMatch)
options(error=recover)
# options(warn=2)
set.seed(1)

max_n_dataset = 500
source("_dataset_load.R")


D = datasets_and_response_models$boston

#try to create a CRD design
n = nrow(D$X)
response_type = "survival"
y = D$y_original[[response_type]]
seq_des_obj = SeqDesign$new(n = n, design = "KK14", response_type = response_type)

# profvis({
  for (t in 1 : n){
    w_t = seq_des_obj$add_subject_to_experiment_and_assign(D$X[t, ])
    seq_des_obj$add_subject_response(t, y[t])
  }
# })

seq_des_inf = SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj)
seq_des_inf$compute_treatment_estimate()
seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
seq_des_inf$compute_mle_confidence_interval(0.05)
seq_des_inf$compute_bootstrap_confidence_interval()
# profvis({
seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand()
# })


seq_des_inf = SeqDesignInferenceKKCompoundMeanDiff$new(seq_des_obj)
seq_des_inf$compute_treatment_estimate()
seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
seq_des_inf$compute_mle_confidence_interval(0.05)
seq_des_inf$compute_bootstrap_confidence_interval()
# profvis({
seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand()
# })

seq_des_inf = SeqDesignInferencePropUniBetaRegr$new(seq_des_obj)

seq_des_inf = SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj)
seq_des_inf$compute_treatment_estimate()
seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
seq_des_inf$compute_mle_confidence_interval()
# profvis({
seq_des_inf$compute_bootstrap_confidence_interval()
# })
# profvis({
seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand()
# })


seq_des_inf = SeqDesignInferenceSurvivalRestrictedMeanDiff$new(seq_des_obj)
seq_des_inf$compute_treatment_estimate()
seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
seq_des_inf$compute_mle_confidence_interval()
# profvis({
seq_des_inf$compute_bootstrap_confidence_interval()
# })
# profvis({
seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand()
# })


seq_des_inf = SeqDesignInferenceSurvivalKMDiff$new(seq_des_obj)
seq_des_inf$compute_treatment_estimate()
seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
# seq_des_inf$compute_mle_confidence_interval()
# profvis({
seq_des_inf$compute_bootstrap_confidence_interval()
# })
# profvis({
seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand()
# })

seq_des_inf = SeqDesignInferenceSurvivalUniWeibullRegr$new(seq_des_obj)
seq_des_inf$compute_treatment_estimate()
seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
seq_des_inf$compute_mle_confidence_interval()
# profvis({
seq_des_inf$compute_bootstrap_confidence_interval()
# })
# profvis({
seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand()
# })

seq_des_inf = SeqDesignInferenceSurvivalMultiWeibullRegr$new(seq_des_obj)
seq_des_inf$compute_treatment_estimate()
seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
seq_des_inf$compute_mle_confidence_interval()
# profvis({
seq_des_inf$compute_bootstrap_confidence_interval()
# })
# profvis({
seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand()
# })


seq_des_inf = SeqDesignInferenceSurvivalUniCoxPHRegr$new(seq_des_obj)
seq_des_inf$compute_treatment_estimate()
seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
seq_des_inf$compute_mle_confidence_interval()
# profvis({
seq_des_inf$compute_bootstrap_confidence_interval()
# })
# profvis({
seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(na.rm = TRUE)
# })

seq_des_inf = SeqDesignInferenceSurvivalMultiCoxPHRegr$new(seq_des_obj)
seq_des_inf$compute_treatment_estimate()
seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
seq_des_inf$compute_mle_confidence_interval()
# profvis({
seq_des_inf$compute_bootstrap_confidence_interval()
# })
# profvis({
seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(na.rm = TRUE)
# })

seq_des_inf = SeqDesignInferencePropUniBetaRegr$new(seq_des_obj)
seq_des_inf$compute_treatment_estimate()
seq_des_inf$compute_mle_confidence_interval()
seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
# profvis({
seq_des_inf$compute_bootstrap_confidence_interval(na.rm = TRUE)
# })
# profvis({
seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(na.rm = TRUE)
# })

seq_des_inf = SeqDesignInferencePropMultiBetaRegr$new(seq_des_obj)
seq_des_inf$compute_treatment_estimate()
seq_des_inf$compute_mle_confidence_interval()
seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
# profvis({
seq_des_inf$compute_bootstrap_confidence_interval(na.rm = TRUE)
# })
# profvis({
seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(na.rm = TRUE)
# })


# seq_des_inf = SeqDesignInferenceIncUnivLogRegr$new(seq_des_obj)
# seq_des_inf$compute_treatment_estimate()
# seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
# seq_des_inf$compute_mle_confidence_interval(0.05)
# profvis({
# seq_des_inf$compute_bootstrap_confidence_interval()
# })
# profvis({
# seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand()
# })

# seq_des_inf = SeqDesignInferenceIncMultiLogRegr$new(seq_des_obj)
# seq_des_inf$compute_treatment_estimate()
# seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
# seq_des_inf$compute_mle_confidence_interval(0.05)
# profvis({
#   seq_des_inf$compute_bootstrap_confidence_interval()
# })
# profvis({
#   seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand()
# })


# seq_des_inf = SeqDesignInferenceContMultOLSKK$new(seq_des_obj)
# seq_des_inf$compute_treatment_estimate()
# seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
# seq_des_inf$compute_mle_confidence_interval()
# profvis({
#   seq_des_inf$compute_bootstrap_confidence_interval()
# })
# profvis({
#   seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand()
# })
