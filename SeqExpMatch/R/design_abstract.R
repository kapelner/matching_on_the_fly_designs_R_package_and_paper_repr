# A Design
#
# @description
# An abstract R6 Class encapsulating the data and functionality for an experimental design.
# This class takes care of data storage and response handling.
#
# @keywords internal
Design = R6::R6Class("Design",
	public = list(
		#
		# @description
		# Initialize an experimental design
		#
		# @param response_type 	The data type of response values which must be one of the following:
		# 							"continuous" (the default),
		# 							"incidence",
		# 							"proportion",
		# 							"count",
		# 							"survival",
		# 							"ordinal".
		# 							This package will enforce that all added responses via the \code{add_subject_response} method will be
		# 							of the appropriate type.
		# @param prob_T	The probability of the treatment assignment. This defaults to \code{0.5}.
		# @param include_is_missing_as_a_new_feature	If missing data is present in a variable, should we include another dummy variable for its
		# 												missingness in addition to imputing its value? If the feature is type factor, instead of creating
		# 												a new column, we allow missingness to be its own level. The default is \code{TRUE}.
		# @param n			The sample size (if fixed). Default is \code{NULL} for not fixed.
		# @param verbose	A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}.
		#
		# @return 			A new `Design` object
		#
		initialize = function(
				response_type = "continuous",
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				num_cores = 1,
				verbose = FALSE
			) {
			assertChoice(response_type, c("continuous", "incidence", "proportion", "count", "survival", "ordinal"))
			assertNumeric(prob_T, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
			assertFlag(include_is_missing_as_a_new_feature)
			assertFlag(verbose)
			assertCount(n, null.ok = TRUE)
			assertCount(num_cores, positive = TRUE)

			if (is.null(n)){
				private$fixed_sample = FALSE
			} else {
				n = as.integer(n)
				private$n = n
				private$fixed_sample = TRUE
			}

			private$prob_T = prob_T
			private$response_type = response_type
			private$include_is_missing_as_a_new_feature = include_is_missing_as_a_new_feature
			private$num_cores = num_cores
			set_package_threads(num_cores)
			private$verbose = verbose

			if (private$fixed_sample){
				private$y = 	rep(NA_real_, n)
				private$w = 	rep(NA_real_, n)
				private$dead =  rep(NA_real_, n)
			}

			if (private$verbose){
				cat(paste0("Initialized a ",
				class(self)[1],
				" experiment with response type ",
				response_type,
				" and ",
				ifelse(private$fixed_sample, "fixed sample", "not fixed sample"),
				 ".\n"))
			}
		},

		# @description
		# Add subject-specific measurements for the next subject entrant
		#
		# @param x_new 			A row of the data frame corresponding to the new subject to be added (must be type data.table).
		# @param allow_new_cols	Should we allow new/different features than previously seen in previous subjects in the
		# 							new subject's covariates? Default is \code{TRUE}.
		#
		add_subject = function(x_new, allow_new_cols = TRUE){
			assertClass(x_new, "data.frame")
			x_new = as.data.table(x_new)
			if (nrow(x_new) != 1){
				stop("You can only add one subject at a time.")
			}
			j_with_NAs = is.na(unlist(x_new))
			if (any(j_with_NAs) & private$t == 0){
				x_new = x_new[which(!j_with_NAs)]
				if (!allow_new_cols){
					warning("There is missing data in the first subject's covariate value(s). Setting the flag allow_new_cols = FALSE will disallow additional subjects")
				}
			}

			xnew_data_types = get_column_types_cpp(x_new)
			if ("ordered" %in% xnew_data_types){
				stop("Ordered factor data type is not supported; please convert to either an unordered factor or numeric.")
			}
			if ("Date" %in% xnew_data_types){
				stop("Date data type is not supported; please convert to numeric.")
			}

			if (private$t > 0){
				Xraw_data_types = get_column_types_cpp(private$Xraw)
				colnames_Xraw = names(private$Xraw)
				colnames_xnew = names(x_new)
				if (setequal(colnames_Xraw, colnames_xnew)){
					idx_data_types_that_changed = which(xnew_data_types != Xraw_data_types)
					if (length(idx_data_types_that_changed) > 0){
						for (e in idx_data_types_that_changed){
							warning("You entered data type ", xnew_data_types[e], " for attribute named ", colnames_Xraw[e], " that was previously entered with data type ", Xraw_data_types[e])
						}
					}
				} else {
					if (allow_new_cols){ #make NA's in appropriate places
						new_Xraw_cols = setdiff(colnames_xnew, colnames_Xraw)
						if (length(new_Xraw_cols) > 0){
							new_Xraw_col_types = xnew_data_types[new_Xraw_cols]
							for (j in 1 : length(new_Xraw_cols)){
								private$Xraw[, (new_Xraw_cols[j]) := switch(new_Xraw_col_types[j],
									character = NA_character_,
									factor =    NA_character_, #I think this is correct
									numeric =   NA_real_,
									logical =   NA_real_, #just let it be zero or one
									integer =   NA_real_  #I don't want to take the risk on a decimal popping up somewhere
								)]
							}
						}

						new_xnew_cols = setdiff(colnames_Xraw, colnames_xnew)
						if (length(new_xnew_cols) > 0){
							new_xnew_cols_types = Xraw_data_types[new_xnew_cols]
							for (j in 1 : length(new_xnew_cols)){
								x_new[, (new_xnew_cols[j]) := switch(new_xnew_cols_types[j],
									character = NA_character_,
									factor =    NA_character_, #I think this is correct
									numeric =   NA_real_,
									logical =   NA_real_, #just let it be zero or one
									integer =   NA_real_  #I don't want to take the risk on a decimal popping up somewhere
								)]
							}
						}


					} else {
						stop(paste(
							"The new subject vector has columns:\n  ",
							paste(colnames_xnew, collapse = ", "),
							"\nwhich are not the same as the current dataset's columns:\n  ",
							paste(colnames_Xraw, collapse = ", "),
							"\nIf you want to allow new columns on-the-fly, run this function again with the option\n  'allow_new_cols = TRUE'"
						))
					}
				}
			}
			#add new subject's measurements to the raw data frame (there should be the same exact columns even if there are new ones introduced)
			private$Xraw = rbindlist(list(private$Xraw, x_new))
			private$p_raw_t = ncol(private$Xraw)

			#iterate t
			private$t = private$t + 1L #t must be an integer for data.table's fast "set" function below to work

			#we only bother with imputation and model matrices if we have enough data otherwise it's a huge mess
			#thus, designs cannot utilize imputations nor model matrices until this condition is met
			#luckily, those are the designs implemented herein so we have complete control (if you are extending this package, you'll have to deal with this issue here)
			if (private$t > (ncol(private$Xraw) + 2) & private$uses_covariates){ #we only need to impute if we need the X's to make the allocation decisions
				private$covariate_impute_if_necessary_and_then_create_model_matrix()
			}
		},

		# @description
		# For CARA designs, add subject response for the a subject
		#
		# @param t 	 The subject index for which to attach a response (beginning with 1, ending with n). You cannot add responses
		# 				 for subjects that have not yet been added to the experiment via the \code{add_subject_to_experiment_and_assign} method.
		# @param y 	 The response value which must be appropriate for the response_type.
		# @param dead	 If the response is censored, enter 0 for this value. This is only necessary to specify for response type
		# 				 "survival" otherwise do not specify this argument (as it will default to 1).
		#
		add_subject_response = function(t, y, dead = 1) {
			assertNumeric(t, len = 1) #make sure it's length one here
			assertNumeric(y, len = 1) #make sure it's length one here
			assertNumeric(dead, len = 1) #make sure it's length one here
			assertChoice(dead, c(0, 1))
			assertCount(t, positive = TRUE)
			if (t > private$t){
				stop(paste("You cannot add response for subject", t, "when the most recent subjects' record added is", private$t))
			}

			if (length(private$y) >= t & !is.na(private$y[t])){
				warning(paste("Overwriting previous response for t =", t, "y[t] =", private$y[t]))
			}

			#deal with the myriad checks on the response value based on response_type
			if (private$response_type == "continuous"){
				assertNumeric(y, any.missing = FALSE)
			} else if (private$response_type == "incidence"){
				assertChoice(y, c(0, 1))
			} else if (private$response_type == "proportion"){
				assertNumeric(y, any.missing = FALSE, lower = 0, upper = 1)
			} else if (private$response_type == "count"){
				assertCount(y, na.ok = FALSE)
			} else if (private$response_type == "survival"){
				assertNumeric(y, any.missing = FALSE, lower = 0)
				if (y == 0){
					warning("0 survival responses not allowed --- recording .Machine$double.eps as its value instead")
					y = .Machine$double.eps
				}
			} else if (private$response_type == "ordinal"){
				if (is.factor(y)){
					assertFactor(y, ordered = TRUE, any.missing = FALSE)
					y = as.integer(y)
				} else {
					assertCount(y, positive = TRUE, na.ok = FALSE)
				}
			}
			if (dead == 0 & private$response_type != "survival"){
				stop("censored observations are only available for survival response types")
			}
			#finally, record the response value and the time at which it was recorded
			if (private$fixed_sample | t <= length(private$y)){
				private$y[t] = y
				private$dead[t] = dead
			} else if (t == length(private$y) + 1){
				private$y = c(private$y, y)
				private$dead = c(private$dead, dead)
			} else {
				stop("You cannot add a response for a subject that has not yet arrived when the sample size is not fixed in advance.")
			}
			private$y_i_t_i[[t]] = private$t
		},

		# @description
		# For non-CARA designs, add all subject responses (usually done at the end of the study)
		#
		# @param ys 		The responses as a numeric vector of length n
		# @param deads	    The binary vector of length n where 1 indicates the the subject
		# 					is dead (survival value is uncensored) and 0 indicates the subject is
		# 					alive (survival value is censored). This is only necessary for response type
		# 				 	"survival" otherwise do not specify and the value
		# 					will default to 1.
		#
		add_all_subject_responses = function(ys, deads = NULL) {
			if (is.null(deads)){
				deads = rep(1, private$t)
			}
			if (private$response_type == "ordinal" && is.factor(ys)){
				assertFactor(ys, len = private$t, ordered = TRUE, any.missing = FALSE)
			} else {
				assertNumeric(ys, len = private$t)
			}
			assertNumeric(deads, len = private$t)

			for (t in 1 : private$t){
				self$add_subject_response(t, ys[t], deads[t]) #piggy back on the checks therein
			}
		},

		# @description
		# Check if this design was initialized with a fixed sample size n
		#
		is_fixed_sample_size = function(){
			private$fixed_sample
		},

		# @description
		# Asserts if the experiment is completed (all n assignments are assigned
		# in the w vector and all n responses in the y vector are recorded), i.e. throws
		# descriptive error if the experiment is incomplete.
		#
		assert_experiment_completed = function(){
			#cat("Design assert_experiment_completed\n")
			if (private$fixed_sample & private$all_assignments_not_yet_allocated()){
				stop("This experiment is incomplete as all n assignments aren't administered yet.")
			}
			if (private$all_responses_not_yet_recorded()){
				stop("This experiment is incomplete as all responses aren't recorded yet.")
			}
		},

		# @description
		# Checks if the experiment is completed (all n assignments are assigned
		# in the w vector and all n responses in the y vector are recorded).
		#
		# @return	\code{TRUE} if experiment is complete, \code{FALSE} otherwise.
		#
		check_experiment_completed = function(){
			if (private$fixed_sample & private$all_assignments_not_yet_allocated(self$t)){
				FALSE
			} else if (private$all_responses_not_yet_recorded()){
				FALSE
			} else {
				TRUE
			}
		},

		# @description
		# Checks if the experiment has a 50-50 allocation to treatment and control
		#
		# @return	\code{TRUE} if 50-50, \code{FALSE} otherwise.
		#
		assert_even_allocation = function(){
			if (private$prob_T != 0.5){
		 		stop("This type of design currently only works with even treatment allocation, i.e. you must set prob_T = 0.5 upon initialization")
		 	}
		},

		#
		# @description
		# Checks if the experiment has a 50-50 allocation to treatment and control
		#
		# @return	\code{TRUE} if 50-50, \code{FALSE} otherwise.
		#
		assert_fixed_sample = function(){
			if (!private$fixed_sample){
				stop("This type of design currently only works with fixed sample, i.e., you must specify n upon initialization")
			}
		},

		# @description
		# Checks if the experiment has any censored responses
		#
		# @return	\code{TRUE} if there are any censored responses, \code{FALSE} otherwise.
		#
		any_censoring = function(){
			sum(private$dead) < length(private$dead)
		},

		# @description Get t
		#
		# @return 			The current number of subjects in this experiment (begins at zero).
		#
		get_t = function(){
			private$t
		},

		# @description Get raw X information
		#
		# @return 			A data frame (data.table object) of subject data with number of rows n (the number of subjects) and number of
		# 					columns p (the number of characteristics measured for each subject). This data frame is filled in
		# 					sequentially by the experimenter and thus will have data present for rows 1...t (i.e. the number of subjects in the
		# 					experiment currently) but otherwise will be missing.
		get_X_raw = function(){
			private$Xraw
		},

		# @description Get imputed X information
		#
		# @return 		Same as \code{Xraw} except with imputations for missing values (if necessary) and deletions of linearly dependent columns
		get_X_imp = function(){
			private$Ximp
		},

		# @description Get X matrix
		#
		# @return 			A numeric matrix of subject data with number of rows n (the number of subjects) and number of
		# 					columns p (the number of characteristics measured for each subject).
		get_X = function(){
			private$X
		},

		# @description Get y
		#
		# @return 			A numeric vector of subject responses with number of entries n (the number of subjects).
		get_y = function(){
			private$y
		},

		# @description Get w
		#
		# @return 			A binary vector of subject assignments with number of entries n (the number of subjects).
		# 					This vector is filled in sequentially by this package (similar to X) and will have assignments present for
		# 					entries 1...t (i.e. the number of subjects in the experiment currently) but otherwise will be missing.
		get_w = function(){
			private$w
		},

		# @description Get n, the sample size
		#
		# If n is fixed, it returns n, if n is not fixed, it returns the current number of subjects, t
		#
		# @return 			The number of subjects
		get_n = function(){
			ifelse(private$fixed_sample, private$n, private$t)
		},

		# @description Get dead
		#
		# @return 			A binary vector of whether the subject is dead with number of entries n (the number of subjects).
		get_dead = function(){
			private$dead
		},

		# @description Get probability of treatment
		#
		# @return 			The experimenter-specified probability a subject becomes treated to the treatment arm.
		get_prob_T = function(){
			private$prob_T
		},

		# @description Get response type
		#
		# @return 			The experimenter-specified response type
		get_response_type = function(){
			private$response_type
		},

		# @description Check whether this design uses blocking / matched-group structure
		#
		# @return 			\code{TRUE} if this design stores grouping information in \code{m}, \code{FALSE} otherwise.
		is_blocking_design = function(){
			!is.null(private$m)
		},

		# @description Check whether this blocking design has equal-sized blocks
		#
		# @return 			\code{TRUE} if all observed blocks have the same size, \code{FALSE} otherwise.
		is_complete_blocking_design = function(){
			if (!self$is_blocking_design()){
				return(FALSE)
			}
			block_sizes = as.integer(table(private$m))
			length(block_sizes) > 0 && length(unique(block_sizes)) == 1
		},

		# @description
		# Duplicate this design object
		#
		# @param verbose 	A flag indicating whether messages should be displayed to the user. Default is \code{FALSE}
		# @return 			A new `Design` object with the same data
		duplicate = function(verbose = FALSE){
			self$assert_experiment_completed() #can't duplicate without the experiment being done
			# Use the built-in R6 clone method (shallow by default) to bypass $new() logic.
			d = self$clone()
			d$.__enclos_env__$private$verbose = verbose
			d
		}
	),


	private = list(
		t = 0L,
		n = NULL,
		Xraw = data.table(),
		p_raw_t = NULL,
		Ximp = data.table(),
		X = NULL,
		w = numeric(),
		y = numeric(),
		dead = numeric(),
		m = NULL,
		prob_T = NULL,
		response_type = NULL,
		fixed_sample = NULL,
		num_cores = NULL,
		include_is_missing_as_a_new_feature = NULL,
		verbose = NULL,
		y_i_t_i = list(),	 #at what point during the experiment are the subjects recorded?
		uses_covariates = FALSE, #does this design use the covariates to make assignments? The default is FALSE

		redraw_w_according_to_design = function(){
			stop("Must be implemented by subclass.")
		},

		all_assignments_allocated = function(){
			if (private$fixed_sample){
				private$t >= private$n
			} else {
				FALSE
			}
		},

		all_assignments_not_yet_allocated = function(){
			!private$all_assignments_allocated()
		},

		all_responses_not_yet_recorded = function(){
			sum(!is.na(private$y)) != length(private$w)
		},

		covariate_impute_if_necessary_and_then_create_model_matrix = function(){
			#make a copy... sometimes the raw will be the same as the imputed if there are no imputations
			private$Ximp = copy(private$Xraw)

			column_has_missingness = columns_have_missingness_cpp(private$Xraw)
			if (any(column_has_missingness)){
				#deal with include_is_missing_as_a_new_feature here
				if (private$include_is_missing_as_a_new_feature){
					missing_cols_idx = which(column_has_missingness)
					if (length(missing_cols_idx) > 0){
						# Use C++ function to create missingness indicators efficiently
						missingness_indicators = create_missingness_indicators_cpp(private$Ximp, missing_cols_idx)

						# Add the new columns to Ximp
						for (col_name in names(missingness_indicators)) {
							private$Ximp[[col_name]] = missingness_indicators[[col_name]]
						}
					}
				}

				#we need to convert characters into factor for the imputation to work
				col_types = get_column_types_cpp(private$Ximp)
				idx_cols_to_convert_to_factor = which(col_types == "character")
				private$Ximp[, (idx_cols_to_convert_to_factor) := lapply(.SD, as.factor), .SDcols = idx_cols_to_convert_to_factor]

				#now do the imputation here by using missRanger (fast but fragile) and if that fails, use missForest (slow but more robust)
				private$Ximp = tryCatch({
										if (any(!is.na(private$y))){
											suppressWarnings(missRanger(cbind(private$Ximp, private$y[1 : nrow(private$Ximp)]), verbose = FALSE, num.threads = private$num_cores)[, 1 : ncol(private$Ximp)])
										} else {
											suppressWarnings(missRanger(private$Ximp, verbose = FALSE, num.threads = private$num_cores))
										}
									}, error = function(e){
										if (any(!is.na(private$y))){
											suppressWarnings(missForest(cbind(private$Ximp, private$y[1 : nrow(private$Ximp)]), num.threads = private$num_cores)$ximp[, 1 : ncol(private$Ximp)])
										} else {
											suppressWarnings(missForest(private$Ximp, num.threads = private$num_cores)$ximp)
										}
									}
								)
			}

			#now let's drop any columns that don't have any variation
			num_unique_values_per_column = count_unique_values_cpp(private$Ximp)
			private$Ximp = private$Ximp[, .SD, .SDcols = which(num_unique_values_per_column > 1)]

			#for nonblank data frames...
			if (ncol(private$Ximp) > 0){
				#now we need to update the numeric model matrix which may have expanded due to new factors, new missingness cols, etc
				private$X = model.matrix(~ ., data = private$Ximp)[, -1, drop = FALSE]
				# Ensure it is a numeric matrix (not character)
				if (is.character(private$X)){
					stop("model.matrix returned a character matrix - this should not happen.")
				}
				private$X = drop_linearly_dependent_cols(private$X)$M

				if (nrow(private$X) != nrow(private$Xraw) | nrow(private$X) != nrow(private$Ximp) | nrow(private$Ximp) != nrow(private$Xraw)){
					stop("improper sizing for the internal X representation")
				}
			} else {
				#blank covariate matrix
				private$X = matrix(NA, nrow = nrow(private$Xraw), ncol = 0)
			}
		},

		compute_all_subject_data = function(){
			i_present_y = which(!is.na(private$y))
			i_all = 1 : private$t
			i_all_y_present = intersect(i_all, i_present_y)

			# Call consolidated C++ function for all matrix computations
			cpp_result = compute_all_subject_data_cpp(
				as.matrix(private$X[1:private$t, , drop = FALSE]),
				private$t,
				as.integer(i_all_y_present)
			)

			# Restore column names
			X_names = colnames(private$X)

			if (length(cpp_result$cols_prev) > 0) {
				nms = X_names[cpp_result$cols_prev]
				colnames(cpp_result$X_prev) = nms
				names(cpp_result$xt_prev) = nms
			}

			if (length(cpp_result$cols_all) > 0) {
				colnames(cpp_result$X_all) = X_names[cpp_result$cols_all]
			}

			if (length(cpp_result$cols_all_scaled) > 0) {
				nms = X_names[cpp_result$cols_all_scaled]
				colnames(cpp_result$X_all_scaled) = nms
				names(cpp_result$xt_all_scaled) = nms
			}

			if (length(cpp_result$cols_all_with_y_scaled) > 0) {
				colnames(cpp_result$X_all_with_y_scaled) = X_names[cpp_result$cols_all_with_y_scaled]
			}

			# Add the simple array slices that don't need C++ optimization
			cpp_result$w_all_with_y_scaled = private$w[i_all_y_present]
			cpp_result$y_all = private$y[i_all_y_present]
			cpp_result$dead_all = private$dead[i_all_y_present]

			cpp_result
		},

		assign_wt_Bernoulli = function(){
			rbinom(1, 1, private$prob_T)
		},

		has_private_method = function(method_name) {
			exists(method_name, envir = self$.__enclos_env__$private, inherits = FALSE)
		}
	)
)
