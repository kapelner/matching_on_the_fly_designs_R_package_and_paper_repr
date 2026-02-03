library(SeqExpMatch)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(survival))

# Set seed for reproducibility
set.seed(123)

n = 50
p = 5

cat("Debugging Bootstrap NAs for CRD continuous SimpleMeanDiff\n")

# Define Designs
design_constructors = list(
  "CRD" = function(type) SeqDesignCRD$new(response_type = type)
)

# Inference classes mapped by response type
inference_map = list(
  "continuous" = list(
    "SimpleMeanDiff" = SeqDesignInferenceAllSimpleMeanDiff
  )
)

generate_data = function(n, p, type) {
  X = matrix(rnorm(n * p), nrow = n, ncol = p)
  X_dt = as.data.table(X)
  y = rnorm(n)
  dead = NULL
  list(X = X_dt, y = y, dead = dead)
}

des_name = "CRD"
resp_type = "continuous"
inf_name = "SimpleMeanDiff"

# Generate data
data = generate_data(n, p, resp_type)

# Instantiate design
des_completed = design_constructors[[des_name]](resp_type)
for (i in 1:n) {
  des_completed$add_subject_to_experiment_and_assign(data$X[i])
}
des_completed$add_all_subject_responses(data$y)

# Create an inference object
inf_obj = inference_map[[resp_type]][[inf_name]]$new(des_completed, num_cores = 1, verbose = TRUE)

# Directly get bootstrap samples and print summary
cat("\n--- Bootstrap Samples --- \n")
beta_hat_T_bs_debug = inf_obj$approximate_bootstrap_distribution_beta_hat_T(B = 101)
print(summary(beta_hat_T_bs_debug))


quit("no")