#' A Sequential Design
#' 
#' @description
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#' This class takes care of data intialization and sequential assignments. The class object
#' should be saved securely after each assignment e.g. on an encrypted cloud server.
#' 
SeqDesign = R6::R6Class("SeqDesign",
	public = list(
		#' 				
		#' @description
		#' Initialize a sequential experimental design
		#' 
		#' @param response_type 	The data type of response values which must be one of the following: 
		#' 							"continuous" (the default),  
		#' 							"incidence", 
		#' 							"proportion", 
		#' 							"count", 
		#' 							"survival".
		#' 							This package will enforce that all added responses via the \code{add_subject_response} method will be
		#' 							of the appropriate type.
		#' @param prob_T	The probability of the treatment assignment. This defaults to \code{0.5}.
		#' @param include_is_missing_as_a_new_feature	If missing data is present in a variable, should we include another dummy variable for its
		#' 												missingness in addition to imputing its value? If the feature is type factor, instead of creating
		#' 												a new column, we allow missingness to be its own level. The default is \code{TRUE}.
		#' @param n			The sample size (if fixed). Default is \code{NULL} for not fixed.
		#' @param verbose	A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}.
		#' @param thin		For internal use only. Do not specify. You can thank R6's single constructor-only for this coding noise.
		#'
		#' @return 			A new `SeqDesign` object of the specific type
		#'  
		initialize = function(
				response_type = "continuous", 
				prob_T = 0.5, 
				include_is_missing_as_a_new_feature = TRUE, 
				n = NULL,
				verbose = FALSE,
				thin = FALSE
			) {
			if (!thin){
				assertChoice(response_type, c("continuous", "incidence", "proportion", "count", "survival"))
				assertNumeric(prob_T, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
				assertFlag(include_is_missing_as_a_new_feature)
				assertFlag(verbose)
				assertCount(n, null.ok = TRUE)
				
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
				private$verbose = verbose
				
				if (private$fixed_sample){
					private$y = 	rep(NA_real_, n)
					private$w = 	rep(NA_real_, n)
					private$dead =  rep(NA_real_, n)
				}
	
				if (private$verbose){
					cat(paste0("Intialized a ", 
					class(self)[1], 
					" experiment with response type ", 
					response_type, 
					" and ", 
					ifelse(private$fixed_sample, "fixed sample", "not fixed sample"),
					 ".\n"))
				}
			}				
		},
		
		#' @description
		#' Add subject-specific measurements for the next subject entrant and return this new subject's treatment assignment
		#' 
		#' @param x_new 			A row of the data frame corresponding to the new subject to be added (must be type data.table).
		#' @param allow_new_cols	Should we allow new/different features than previously seen in previous subjects in the 
		#' 							new subject's covariates? Default is \code{TRUE}.
		#' 
		#' @examples
		#' seq_des = SeqDesignCRD$new(n = 100, p = 10, design = "CRD", response_type = "continuous")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' 
		add_subject_to_experiment_and_assign = function(x_new, allow_new_cols = TRUE) {					
			assertDataFrame(x_new, nrows = 1)
			assertClass(x_new, "data.table")
			if (private$fixed_sample & private$all_assignments_allocated()){
				stop(paste("You cannot add any new subjects as all n =", private$n, "subjects have already been added."))
			}
			j_with_NAs = is.na(unlist(x_new))
			if (any(j_with_NAs) & private$t == 0){
				x_new = x_new[which(!j_with_NAs)]
				if (!allow_new_cols){
					warning("There is missing data in the first subject's covariate value(s). Setting the flag allow_new_cols = FALSE will disallow additional subjects")
				}					
			}
			
			xnew_data_types = sapply(x_new, function(xj) class(xj)[1])
			if ("ordered" %in% xnew_data_types | "Date" %in% xnew_data_types){
				stop("Ordered factor and Date data type is not supported; please convert to either an unordered factor or numeric.")
			}
			
			if (private$t > 0){
				Xraw_data_types = sapply(private$Xraw, function(xj) class(xj)[1]) #sapply(private$Xraw, class)
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
							"The new subject vectorhas columns:\n  ", 
							paste(colnames_xnew, collapse = ", "), 
							"\nwhich are not the same as the current dataset's columns:\n  ", 
							paste(colnames_Xrew, collapse = ", "),
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
			
			#now make the assignment
			if (private$fixed_sample){
				private$w[private$t] = private$assign_wt()
			} else {
				private$w = c(private$w, private$assign_wt())
			}
			private$w[private$t]
		},
		
		#' @description
		#' Prints the current assignment to screen. Should be called after \code{add_subject_to_experiment_and_assign}.
		#' 
		#' @examples
		#' seq_des = SeqDesign$new(n = 100, p = 10, design = "CRD", response_type = "continuous")
		#' 
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$print_current_subject_assignment()
		#' 
		print_current_subject_assignment = function(){
			cat("Subject number", private$t, "is assigned to", ifelse(private$w[private$t] == 1, "TREATMENT", "CONTROL"), "via design", class(self)[1], "\n")
		},
		
		#' @description
		#' For CARA designs, add subject response for the a subject
		#' 
		#' @param t 	 The subject index for which to attach a response (beginning with 1, ending with n). You cannot add responses
		#' 				 for subjects that have not yet been added to the experiment via the \code{add_subject_to_experiment_and_assign} method.
		#' @param y 	 The response value which must be appropriate for the response_type. 
		#' @param dead	 If the response is censored, enter 0 for this value. This is only necessary to specify for response type
		#' 				 "survival" otherwise do not specify this argument (as it will default to 1).
		#' @examples
		#' seq_des = SeqDesign$new(n = 100, p = 10, design = "KK21", response_type = "continuous")
		#' 
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' 
		#' seq_des$add_subject_response(4.71, 1)
		#' #works
		#' seq_des$add_subject_response(4.71, 2)
		#' #fails
		#' 
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
				assertNumeric(y, any.missing = FALSE, lower = 0, upper = 1) #the betareg package can handle 0's and 1's exactly (see their documentation) this seems to be why the statmod / numDeriv packages is also required 
			} else if (private$response_type == "count"){
				assertCount(y, na.ok = FALSE)
			} else if (private$response_type == "survival"){
				assertNumeric(y, any.missing = FALSE, lower = 0)
				if (y == 0){
					warning("0 survival responses not allowed --- recording .Machine$double.eps as its value instead")
					y = .Machine$double.eps
				}						
			}
			if (dead == 0 & private$response_type %in% c("continuous", "incidence", "proportion")){
				stop("censored observations are only available for count or survival response types")
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
		
		#' @description
		#' For non-CARA designs, add all subject responses (usually done at the end of the study)
		#' 
		#' @param ys 		The responses as a numeric vector of length n
		#' @param deads	    The binary vector of length n where 1 indicates the the subject 
		#' 					is dead (survival value is uncensored) and 0 indicates the subject is
		#' 					alive (survival value is censored). This is only necessary for response type
		#' 				 	"survival" otherwise do not specify and the value 
		#' 					will default to 1.
		#' @examples
		#' seq_des = SeqDesign$new(n = 6, p = 10, design = "CRD", response_type = "continuous")
		#' 
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' 
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 				
		add_all_subject_responses = function(ys, deads = NULL) {
			if (is.null(deads)){
				deads = rep(1, private$t)
			}
			assertNumeric(ys, len = private$t)
			assertNumeric(deads, len = private$t)
			
			for (t in 1 : private$t){
				self$add_subject_response(t, ys[t], deads[t])
			}
		},

		#' @description
		#' Check if this design was initialized with a fixed sample size n
		#' 
		#' @examples
		#' seq_des = SeqDesignCRD$new(n = 6, response_type = "continuous")
		#' seq_des$is_fixed_sample_size() #returns TRUE
		#' seq_des = SeqDesignCRD$new(response_type = "continuous")
		#' seq_des$is_fixed_sample_size() #returns FALSE
		#' 
		is_fixed_sample_size = function(){
			private$fixed_sample
		},
		
		#' @description
		#' Asserts if the experiment is completed (all n assignments are assigned
		#' in the w vector and all n responses in the y vector are recorded), i.e. throws 
		#' descriptive error if the experiment is incomplete.
		#' 
		#' @examples
		#' seq_des = SeqDesign$new(n = 6, p = 10, design = "CRD", response_type = "continuous")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' 
		#' #if run, it would throw an error since all of the covariate vectors are not yet recorded
		#' #seq_des$assert_experiment_completed() 
		#' 
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' 
		#' #if run, it would throw an error since the responses are not yet recorded
		#' #seq_des$assert_experiment_completed() 
		#' 
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des$assert_experiment_completed() #no response means the assert is true
		assert_experiment_completed = function(){
			#cat("SeqDesign assert_experiment_completed\n")
			if (private$fixed_sample & private$all_assignments_not_yet_allocated()){
				stop("This experiment is incomplete as all n assignments aren't administered yet.")
			}
			if (private$all_responses_not_yet_recorded()){
				stop("This experiment is incomplete as all responses aren't recorded yet.")
			}
		},
		
		#' @description
		#' Checks if the experiment is completed (all n assignments are assigned
		#' in the w vector and all n responses in the y vector are recorded).
		#' 
		#' @return	\code{TRUE} if experiment is complete, \code{FALSE} otherwise.
		#' 
		#' @examples
		#' seq_des = SeqDesign$new(n = 6, p = 10, design = "CRD", response_type = "continuous")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' 
		#' #returns FALSE since all of the covariate vectors are not yet recorded
		#' seq_des$check_experiment_completed() 
		#' 
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' 
		#' #returns FALSE since the responses are not yet recorded
		#' seq_des$check_experiment_completed() 
		#' 
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des$check_experiment_completed() #returns TRUE
		#' 
		check_experiment_completed = function(){
			if (private$fixed_sample & private$all_assignments_not_yet_allocated(self$t)){
				FALSE
			} else if (private$all_responses_not_yet_recorded()){
				FALSE
			} else {
				TRUE
			}
		},
		
		#' @description
		#' Checks if the experiment has a 50-50 allocation to treatment and control
		#' 
		#' @return	\code{TRUE} if 50-50, \code{FALSE} otherwise.
		#' 
		#' @examples
		#' seq_des = SeqDesignCRD$new(n = 6, p = 10, design = "CRD", response_type = "continuous")
		#' seq_des$assert_even_allocation() #returns TRUE
		#' seq_des = SeqDesignCRD$new(n = 6, p = 10, prob_T = 0.4, design = "CRD", response_type = "continuous")
		#' seq_des$assert_even_allocation() #returns FALSE
		assert_even_allocation = function(){
			if (private$prob_T != 0.5){
				stop("This type of design currently only works with even treatment allocation, i.e. you must set prob_T = 0.5 upon initialization")
			}
		},
		
		#' @description
		#' Checks if the experiment has a 50-50 allocation to treatment and control
		#' 
		#' @return	\code{TRUE} if 50-50, \code{FALSE} otherwise.
		#' 
		#' @examples
		#' seq_des = SeqDesignCRD$new(response_type = "continuous", n = 6)
		#' seq_des$assert_fixed_sample() #returns TRUE
		#' seq_des = SeqDesignCRD$new(response_type = "continuous")
		#' seq_des$assert_fixed_sample() #returns FALSE
		assert_fixed_sample = function(){
			if (!private$fixed_sample){
				stop("This type of design currently only works with fixed sample, i.e., you must specify n upon initialization")
			}
		},
		
		#' @description
		#' Checks if the experiment has any censored responses
		#' 
		#' @return	\code{TRUE} if there are any censored responses, \code{FALSE} otherwise.
		#' 
		#' @examples
		#' seq_des = SeqDesign$new(n = 6, p = 10, design = "CRD", response_type = "continuous")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' 
		#' #returns FALSE since all of the covariate vectors are not yet recorded
		#' seq_des$check_experiment_completed() 
		#' 
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' 
		#' #returns FALSE since the responses are not yet recorded
		#' seq_des$check_experiment_completed() 
		#' 
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des$any_censoring() #returns FALSE
		#' 			
		any_censoring = function(){
			sum(private$dead) < length(private$dead)
		},
		
		#' Get t
		#'
		#' @return 			The current number of subjects in this sequential experiment (begins at zero).					
		#' 			
		get_t = function(){
			private$t
		},
		
		#' Get raw X information
		#'
		#' @return 			A data frame (data.table object) of subject data with number of rows n (the number of subjects) and number of 
		#' 					columns p (the number of characteristics measured for each subject). This data frame is filled in
		#' 					sequentially by the experimenter and thus will have data present for rows 1...t (i.e. the number of subjects in the
		#' 					experiment currently) but otherwise will be missing.					
		get_X_raw = function(){
			private$Xraw
		},
		
		#' Get imputed X information
		#'
		#' @return 		Same as \code{Xraw} except with imputations for missing values (if necessary) and deletions of linearly dependent columns					
		get_X_imp = function(){
			private$Ximp
		},
		
		#' Get y
		#'
		#' @return 			A numeric vector of subject responses with number of entries n (the number of subjects). During
		#' 					the KK21 designs the experimenter fills these values in when they are measured.
		#' 					For non-KK21 designs, this vector can be set at anytime (but must be set before inference is desired).				
		get_y = function(){
			private$y
		},
		
		#' Get w
		#'
		#' @return 			A binary vector of subject assignments with number of entries n (the number of subjects). 
		#' 					This vector is filled in sequentially by this package (similar to X) and will have assignments present for
		#' 					entries 1...t (i.e. the number of subjects in the experiment currently) but otherwise will be missing.	
		get_w = function(){
			private$w
		},
		
		#' Get n, the sample size
		#'
		#' If n is fixed, it returns n, if n is not fixed, it returns the current number of subjects, t
		#'
		#' @return 			The number of subjects
		get_n = function(){
			ifelse(private$fixed_sample, private$n, private$t)
		},
		
		#' Get dead
		#'
		#' @return 			A binary vector of whether the subject is dead with number of entries n (the number of subjects). This 
		#' 					vector is filled in only for \code{response_type} values "survival". The value
		#' 					of 1 indicates uncensored (as the subject died) and a value 0 indicates the real survival value is censored 
		#' 					as the subject is still alive at the time of measurement. This follows the same convention as the \code{event} 
		#' 					argument in the canonical \code{survival} package in the constructor \code{survival::Surv}. During
		#' 					the KK21 designs the experimenter fills these values in when they are measured.
		#' 					For non-KK21 designs, this vector can be set at anytime (but must be set before inference is desired).				
		get_dead = function(){
			private$dead
		},
		
		#' Get probability of treatment
		#'
		#' @return 			The experimenter-specified probability a subject becomes wtated to the treatment arm.		
		get_prob_T = function(){
			private$prob_T
		},

		#' Get response type
		#'
		#' @return 			The experimenter-specified response type which is one of the following: 
		#' 							"continuous", 
		#' 							"incidence", 
		#' 							"proportion", 
		#' 							"count", 
		#' 							"survival"			
		get_response_type = function(){
			private$response_type
		},
		
		
		#' Get covariate weights
		#'
		#' @return 			For KK21 designs, the running values of the weights for each covariate.
		get_covariate_weights = function(){
			private$covariate_weights
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
		prob_T = NULL,
		response_type = NULL,
		fixed_sample = NULL,
		include_is_missing_as_a_new_feature = NULL,
		verbose = NULL,		
		y_i_t_i = list(),	 #at what point during the experiment are the subjects recorded?	
		
		thin_duplicate = function(){
			do.call(get(class(self)[1])$new, args = list(thin = TRUE))		
		},

		duplicate = function(){
			self$assert_experiment_completed() #can't duplicate without the experiment being done
			seq_design_class_constructor = get(class(self)[1])$new
			d = do.call(seq_design_class_constructor, args = list(
				response_type = 						private$response_type,  
				prob_T = 								private$prob_T, 
				include_is_missing_as_a_new_feature = 	private$include_is_missing_as_a_new_feature, 
				verbose = 								FALSE, #when we're duplicating, we want messages to be off
				n = 									private$n
			))
			#we are assuming the experiment is complete so we have private$t = 0 initialized
			d$.__enclos_env__$private$Xraw = 					private$Xraw
			d$.__enclos_env__$private$Ximp = 					private$Ximp
			d$.__enclos_env__$private$y = 						private$y
			d$.__enclos_env__$private$dead = 					private$dead
			d$.__enclos_env__$private$w = 						private$w		
			d$.__enclos_env__$private$t = 						private$t	
			d$.__enclos_env__$private$X = 						private$X	
			d$.__enclos_env__$private$fixed_sample = 			private$fixed_sample
			d$.__enclos_env__$private$y_i_t_i = 				private$y_i_t_i
			d
		},
		
		redraw_w_according_to_design = function(){
			for (t in 1 : private$t){
				private$w[t] = private$assign_wt()
			}
		},		
		
		all_assignments_allocated = function(){
			#cat("SeqDesign all_assignments_allocated\n")
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
			#cat("SeqDesign all_responses_not_yet_recorded\n")
			sum(!is.na(private$y)) != length(private$w)
		},
		
		covariate_impute_if_necessary_and_then_create_model_matrix = function(){
			#make a copy... sometimes the raw will be the same as the imputed if there are no imputations
			private$Ximp = copy(private$Xraw)
			
			column_has_missingness = private$Xraw[, lapply(.SD, function(xj) sum(is.na(xj)))] > 0
			if (any(column_has_missingness)){
				#deal with include_is_missing_as_a_new_feature here
				if (private$include_is_missing_as_a_new_feature){
					for (j in which(column_has_missingness)){
						private$Ximp = cbind(private$Ximp, ifelse(is.na(private$Ximp[, .SD, .SDcols = j]), 1, 0))
						names(private$Ximp)[ncol(private$Ximp)] = paste0(names(private$Ximp)[j], "_is_missing")
					}
				}
				
				#we need to convert characters into factor for the imputation to work
				col_types = private$Ximp[, lapply(.SD, function(xj){class(xj)})]
				idx_cols_to_convert_to_factor = which(col_types == "character")
				private$Ximp[, (idx_cols_to_convert_to_factor) := lapply(.SD, as.factor), .SDcols = idx_cols_to_convert_to_factor]
				
				#now do the imputation here by using missRanger (fast but fragile) and if that fails, use missForest (slow but more robust)
				private$Ximp = tryCatch({
										if (any(!is.na(private$y))){
											suppressWarnings(missRanger(cbind(private$Ximp, private$y[1 : nrow(private$Ximp)]), verbose = FALSE)[, 1 : ncol(private$Ximp)])
										} else {
											suppressWarnings(missRanger(private$Ximp, verbose = FALSE))
										}
									}, error = function(e){
										if (any(!is.na(private$y))){
											suppressWarnings(missForest(cbind(private$Ximp, private$y[1 : nrow(private$Ximp)]))$ximp[, 1 : ncol(private$Ximp)])
										} else {
											suppressWarnings(missForest(private$Ximp)$ximp)
										}
									}
								)
			}
			
			#now let's drop any columns that don't have any variation
			num_unique_values_per_column = private$Ximp[, lapply(.SD, function(xj){uniqueN(xj)})]
			private$Ximp = private$Ximp[, .SD, .SDcols = which(num_unique_values_per_column > 1)]
			
			#for nonblank data frames...
			if (ncol(private$Ximp) > 0){
				#now we need to convert character features into factors
				col_types = private$Ximp[, lapply(.SD, function(xj){class(xj)})]
				idx_cols_to_convert_to_factor = which(col_types == "character")
				private$Ximp[, (idx_cols_to_convert_to_factor) := lapply(.SD, as.factor), .SDcols = idx_cols_to_convert_to_factor]
				
				#now we need to update the numeric model matrix which may have expanded due to new factors, new missingness cols, etc
				private$X = model.matrix(~ ., private$Ximp)
				private$X = private$drop_linearly_dependent_cols(private$X)$M
				
				if (nrow(private$X) != nrow(private$Xraw) | nrow(private$X) != nrow(private$Ximp) | nrow(private$Ximp) != nrow(private$Xraw)){
					stop("boom")
				}
			} else {
				#blank covariate matrix
				private$X = matrix(NA, nrow = nrow(private$Xraw), ncol = 0)
			}			
		},
		
		compute_all_subject_data = function(){
			i_present_y = which(!is.na(private$y))
			i_all = 1 : private$t
			i_past = 	if (private$t == 1){
							c()
						} else {
							1 : (private$t - 1)
						}
			i_past_y_present = intersect(i_past, i_present_y)
			i_all_y_present =  intersect(i_all,  i_present_y)
			
			if (private$t == 1){
				xt = private$X[private$t, ]
				Xall = matrix(xt, nrow = 1)
				Xall_with_y = Xall[i_all_y_present, , drop = FALSE]
				list(
					Xprev = NA, 
					rank_previous = NA,
					xt_prev = NA,
					w_prev = NA,
					
					Xprev_with_y = NA,
					rank_previous_with_y = NA,
					xt_prev_with_y = NA,
					w_prev_with_y = NA,
					
					Xall = Xall,
					rank_all = length(xt),
					xt_all = xt,
					w_all = private$w[private$t],
					
					Xall_with_y = Xall_with_y,
					rank_all_with_y = ifelse(nrow(Xall_with_y) == 1, length(xt), NA),
					xt_all = ifelse(nrow(Xall_with_y) == 1, xt, NA),
					w_all = private$w[private$t],
					
					y_previous = NA,
					y_all = private$y[i_all_y_present],
					
					dead_previous = NA,
					dead_all = private$dead[i_all_y_present]
				)	
			} else {
				past_info = 				private$remove_linearly_dependent_covariates_and_compute_info(i_past)
				past_info_with_y = 			private$remove_linearly_dependent_covariates_and_compute_info(i_past_y_present)
				all_info = 					private$remove_linearly_dependent_covariates_and_compute_info(i_all)
				all_info_with_y = 			private$remove_linearly_dependent_covariates_and_compute_info(i_all_y_present)
				all_info_scaled = 			private$remove_linearly_dependent_covariates_and_compute_info(i_all, scaled = TRUE)
				all_info_with_y_scaled = 	private$remove_linearly_dependent_covariates_and_compute_info(i_all_y_present, scaled = TRUE)
		

				list(
					X_prev = past_info$Xint, 
					rank_prev = past_info$rank,
					xt_prev = past_info$xt,
					w_prev = private$w[i_past],
					
					X_prev_with_y = past_info_with_y$Xint,
					rank_prev_with_y = past_info_with_y$rank,
					xt_prev_with_y = past_info_with_y$xt,
					w_prev_with_y = private$w[i_past_y_present],
					
					X_all = all_info$Xint,
					rank_all = all_info$rank,
					xt_all = all_info$xt,
					w_all = private$w[i_all],
					
					X_all_with_y = all_info_with_y$Xint,
					rank_all_with_y = all_info_with_y$rank,
					xt_all_with_y = all_info_with_y$xt,
					w_all_with_y = private$w[i_all_y_present],
					
					X_all_scaled = all_info_scaled$Xint,
					rank_all_scaled = all_info_scaled$rank,
					xt_all_scaled = all_info_scaled$xt,
					w_all_scaled = private$w[i_all],
					
					X_all_with_y_scaled = all_info_with_y_scaled$Xint,
					rank_all_with_y_scaled = all_info_with_y_scaled$rank,
					xt_all_with_y_scaled = all_info_with_y_scaled$xt,
					w_all_with_y_scaled = private$w[i_all_y_present],
					
					y_prev = private$y[i_past_y_present],
					y_all = private$y[i_all_y_present],
					
					dead_prev = private$dead[i_past_y_present],
					dead_all = private$dead[i_all_y_present]
				)
			}		
		},
		
		remove_linearly_dependent_covariates_and_compute_info = function(is, scaled = FALSE){
			if (length(is) == 0){				
				list(Xint = NA, rank = NA, xt = NA)
			} else {
				Xint = private$X[is, , drop = FALSE]
				
				#we first kill columns that have no variation
				js = which_cols_vary_cpp(Xint)
				Xint = Xint[, js]
				xt = private$X[private$t, js]
				
				if (scaled & length(is) > 1){
					#we need to scale all data including the tth
					Xint = rbind(Xint, xt)
					Xint_colnames = colnames(Xint)
					Xint = scale_columns_cpp(Xint)
					colnames(Xint) = Xint_colnames
					#if the column is all one unique value, it goes NaN after scaling, so we have to patch that up
					Xint[is.nan(Xint)] = 0
				}
				
				drop_obj = private$drop_linearly_dependent_cols(Xint)
				Xint = drop_obj$M
				xt = xt[drop_obj$js] #make sure xt comports with Xint!
				rank = matrix_rank_cpp(Xint)
				
				if (scaled & length(is) > 1){
					list(Xint = Xint[1 : (nrow(Xint) - 1), , drop = FALSE], rank = rank, xt = Xint[nrow(Xint), ])
				} else {
					list(Xint = Xint, rank = rank, xt = xt)	
				}								
			}
		},
		
		drop_linearly_dependent_cols = function(M){
			rank = matrix_rank_cpp(M)
			js = 1 : ncol(M)
			#it's possible that there may be linearly dependent columns
			if (rank != ncol(M)){
				#kill linearly dependent column(s) via cool trick found at
				#https://stackoverflow.com/questions/19100600/extract-maximal-set-of-independent-columns-from-a-matrix
				qrX = qr(M)
				js = qrX$pivot[seq_len(qrX$rank)] #the true linearly independent column indicies
				M = M[, js, drop = FALSE]
			}
			list(M = M, js = js)
		},
		
		
		assign_wt_CRD = function(){
			rbinom(1, 1, private$prob_T)
		}
				
		

	)
)