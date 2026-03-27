library(SeqExpMatch)
library(data.table)

# Minimal data for incidence
set.seed(1)
n = 20
X = data.frame(x = rnorm(n))
y = rbinom(n, 1, 0.5)

seq_des = DesignSeqOneByOneBernoulli$new(response_type = "incidence", n = n)
for (i in 1:n){
    seq_des$add_subject_to_experiment_and_assign(X[i, , drop=FALSE])
}
seq_des$add_all_subject_responses(y)

inf = InferenceIncidExactZhang$new(seq_des)

cat("Testing compute_exact_two_sided_pval_for_treatment_effect()...
")
res = tryCatch({
    inf$compute_exact_two_sided_pval_for_treatment_effect()
}, error = function(e) {
    cat("Error caught:", e$message, "
")
    # Print traceback-like info if possible
    return(NULL)
})

if (!is.null(res)) {
    cat("Result:", res, "
")
}

cat("
Testing compute_exact_confidence_interval()...
")
res_ci = tryCatch({
    inf$compute_exact_confidence_interval(pval_epsilon = 0.1) # Fast tolerance
}, error = function(e) {
    cat("CI Error caught:", e$message, "
")
    return(NULL)
})

if (!is.null(res_ci)) {
    cat("CI Result:", res_ci, "
")
}
