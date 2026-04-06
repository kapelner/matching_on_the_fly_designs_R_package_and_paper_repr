MockKKGEEFallback <- R6::R6Class(
	"MockKKGEEFallback",
	inherit = InferenceAbstractKKGEE,
	private = list(
		gee_response_type = function() "incidence",
		gee_family = function() binomial(link = "logit"),
		gee_predictors_df = function() data.frame(w = private$w),
		fit_gee_on_data = function(dat, id_sorted, std_err = TRUE) {
			cluster_sizes = table(id_sorted)
			if (any(cluster_sizes == 1L)) return(NULL)
			stats::lm(as.numeric(y) ~ w, data = dat)
		}
	)
)

test_that("KK GEE falls back to matched-only fit when reservoir singletons break the full fit", {
	des <- DesignSeqOneByOneKK14$new(n = 8, response_type = "incidence", verbose = FALSE)
	X <- data.frame(x1 = seq_len(8), x2 = seq_len(8) %% 2)
	for (i in seq_len(nrow(X))) {
		des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
	}
	des$add_all_subject_responses(c(0, 1, 0, 1, 0, 1, 0, 1))
	des$.__enclos_env__$private$m <- c(1L, 1L, 2L, 2L, 0L, 0L, 0L, 0L)

	inf <- MockKKGEEFallback$new(des, verbose = FALSE)
	priv <- inf$.__enclos_env__$private

	expect_true(priv$gee_has_reservoir())
	expect_null(priv$fit_gee(std_err = FALSE, include_reservoir = TRUE))
	mod_fb = priv$fit_gee_with_fallback(std_err = FALSE, estimate_only = TRUE)
	expect_true(!is.null(mod_fb))
	expect_true(is.finite(priv$extract_gee_treatment_estimate(mod_fb)))
})
