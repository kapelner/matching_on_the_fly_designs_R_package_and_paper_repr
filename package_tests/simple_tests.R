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
dead = rep(1, n)
response_type = "continuous"
y = D$y_original[[response_type]]

#add some censoring if it's survival
if (response_type == "survival"){
  prob_censoring = 0.3
  for (i in 1 : n){
    if (runif(1) < prob_censoring){
      y[i] = runif(1, 0, y[i])
      dead[i] = 0
    }
  }  
}


#all current designs
seq_des_obj = SeqDesignCRD$new(response_type = response_type, n = n)
seq_des_obj = SeqDesigniBCRD$new(response_type = response_type, n = n)
seq_des_obj = SeqDesignEfron$new(response_type = response_type, n = n)
seq_des_obj = SeqDesignKK14$new(response_type = response_type, n = n)
seq_des_obj = SeqDesignKK21$new(response_type = response_type, n = n)
seq_des_obj = SeqDesignKK21stepwise$new(response_type = response_type, n = n)

# profvis({
  for (t in 1 : n){
    seq_des_obj$add_subject_to_experiment_and_assign(D$X[t, ])
    seq_des_obj$add_subject_response(t, y[t], dead[t])
  }
# })

#all response types
seq_des_inf = SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj)
seq_des_inf = SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj) #only KK

#continuous only
seq_des_inf = SeqDesignInferenceContinMultOLSKK$new(seq_des_obj)

#incidence only
seq_des_inf = SeqDesignInferenceIncidUnivLogRegr$new(seq_des_obj)
seq_des_inf = SeqDesignInferenceIncidMultiLogRegr$new(seq_des_obj)

#proportion only
seq_des_inf = SeqDesignInferencePropUniBetaRegr$new(seq_des_obj)
seq_des_inf = SeqDesignInferencePropMultiBetaRegr$new(seq_des_obj)

#count only
seq_des_inf = SeqDesignInferenceCountUnivNegBinRegr$new(seq_des_obj)
seq_des_inf = SeqDesignInferenceCountMultiNegBinRegr$new(seq_des_obj)

#survival only
seq_des_inf = SeqDesignInferenceSurvivalRestrictedMeanDiff$new(seq_des_obj)
seq_des_inf = SeqDesignInferenceSurvivalKMDiff$new(seq_des_obj)
seq_des_inf = SeqDesignInferenceSurvivalUniWeibullRegr$new(seq_des_obj)
seq_des_inf = SeqDesignInferenceSurvivalMultiWeibullRegr$new(seq_des_obj)
seq_des_inf = SeqDesignInferenceSurvivalUniCoxPHRegr$new(seq_des_obj)
seq_des_inf = SeqDesignInferenceSurvivalMultiCoxPHRegr$new(seq_des_obj)


seq_des_inf$compute_treatment_estimate()
seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
seq_des_inf$compute_mle_confidence_interval(0.05)
# head(seq_des_inf$approximate_bootstrap_distribution_beta_hat_T())
seq_des_inf$compute_bootstrap_confidence_interval()
seq_des_inf$compute_bootstrap_two_sided_pval()
# head(seq_des_inf$compute_beta_hat_T_randomization_distr_under_sharp_null())
seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand()
seq_des_inf$compute_confidence_interval_rand()
seq_des_inf$set_custom_randomization_statistic_function(function(){
  yTs = private$seq_des_obj_priv_int$y[private$seq_des_obj_priv_int$w == 1]
  yCs = private$seq_des_obj_priv_int$y[private$seq_des_obj_priv_int$w == 0]
  (mean(yTs) - mean(yCs)) / sqrt(var(yTs) / length(yTs) + var(yCs) / length(yCs))
})
seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand()
seq_des_inf$compute_confidence_interval_rand()
seq_des_inf$set_custom_randomization_statistic_function(NULL)

profvis({
seq_des_inf$compute_bootstrap_confidence_interval()
})
profvis({
seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand()
})
profvis({
  seq_des_inf$compute_confidence_interval_rand()
})
