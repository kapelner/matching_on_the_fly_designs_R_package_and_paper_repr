#' A Sequential Design
#' 
#' @description
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#' This class takes care of data intialization and sequential assignments. The class object
#' should be saved securely after each assignment e.g. on an encrypted cloud server.
#' 
#' @export
SeqDesign = R6::R6Class("SeqDesign",
	public = list(
			#' @field t			The current number of subjects in this sequential experiment (begins at zero).				
			#' @field design	The experimenter-specified type of sequential experimental design (see constructor's documentation).	
			#' @field Xraw		A data frame (data.table object) of subject data with number of rows n (the number of subjects) and number of 
			#' 					columns p (the number of characteristics measured for each subject). This data frame is filled in
			#' 					sequentially by the experimenter and thus will have data present for rows 1...t (i.e. the number of subjects in the
			#' 					experiment currently) but otherwise will be missing.					
			#' @field Ximp		Same as \code{Xraw} except with imputations for missing values (if necessary) and deletions of linearly dependent columns
			#' @field X			Same as \code{Ximp} except turned into a model matrix (i.e. all numeric with factors dummified) with no linearly dependent columns 
			#' 					(and it is also a matrix object, not a data.table object)
			#' @field y			A numeric vector of subject responses with number of entries n (the number of subjects). During
			#' 					the KK21 designs the experimenter fills these values in when they are measured.
			#' 					For non-KK21 designs, this vector can be set at anytime (but must be set before inference is desired).				
			#' @field dead		A binary vector of whether the subject is dead with number of entries n (the number of subjects). This 
			#' 					vector is filled in only for \code{response_type} values "survival". The value
			#' 					of 1 indicates uncensored (as the subject died) and a value 0 indicates the real survival value is censored 
			#' 					as the subject is still alive at the time of measurement. This follows the same convention as the \code{event} 
			#' 					argument in the canonical \code{survival} package in the constructor \code{survival::Surv}. During
			#' 					the KK21 designs the experimenter fills these values in when they are measured.
			#' 					For non-KK21 designs, this vector can be set at anytime (but must be set before inference is desired).				
			#' @field prob_T	The experimenter-specified probability a subject becomes wtated to the treatment arm.
			#' @field w			A binary vector of subject assignments with number of entries n (the number of subjects). 
			#' 					This vector is filled in sequentially by this package (similar to X) and will have assignments present for
			#' 					entries 1...t (i.e. the number of subjects in the experiment currently) but otherwise will be missing.	
			#' @field response_type		This is the experimenter-specified type of response value which is one of the following: 
			#' 							"continuous", 
			#' 							"incidence", 
			#' 							"proportion", 
			#' 							"count", 
			#' 							"survival"	
			#' @field covariate_weights		The running values of the weights for each covariate
			#' 
			t = 0,
			design = NULL,
			Xraw = data.table(),
			Ximp = data.table(),
			X = NULL,
			y = NULL,	
			dead = NULL,
			prob_T = NULL,
			w = NULL,
			response_type = NULL,
			covariate_weights = NULL,
			#' 				
			#' @description
			#' Initialize a sequential experimental design
			#' 
			#' @param design	The type of sequential experimental design. This must be one of the following
			#' 					"CRD" for the completely randomized design / Bernoulli design, 
			#' 					"iBCRD" for the incomplete / balanaced completely randomized design with appropriate permuted blocks based on \code{prob_T}
			#' 					(e.g., if \code{prob_T = 2}, then this design would enforce n/2 T's and n/2 C's),
			#' 					"Efron" for Efron's (1971) Biased Coin Design
			#' 					"Atkinson" for Atkinson's (1982) Covariate-Adjusted Biased Coin Design
			#' 					"KK14" for Kapelner and Krieger's (2014) Covariate-Adjusted Matching on the Fly Design
			#' 					"KK21" for Kapelner and Krieger's (2021) CARA Matching on the Fly with Differential Covariate Weights Design
			#' 					"KK21stepwise" for Kapelner and Krieger's (2021) CARA Matching on the Fly with Differential Covariate Weights Stepwise Design
			#' @param response_type 	The data type of response values which must be one of the following: 
			#' 							"continuous", 
			#' 							"incidence", 
			#' 							"proportion", 
			#' 							"count", 
			#' 							"survival".
			#' 							This package will enforce that all added responses via the \code{add_subject_response} method will be
			#' 							of the appropriate type.
			#' @param n 		Number of subjects fixed beforehand.
			#' @param prob_T	The probability of the treatment assignment. This defaults to \code{0.5}.
			#' @param include_is_missing_as_a_new_feature	If missing data is present in a variable, should we include another dummy variable for its
			#' 												missingness in addition to imputing its value? If the feature is type factor, instead of creating
			#' 												a new column, we allow missingness to be its own level. The default is \code{TRUE}.  
			#' @param verbose	A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}.
			#' @param ...		Design-specific parameters:
			#' 					"Efron" requires "weighted_coin_prob" which is the probability of the weighted coin for assignment. If unspecified, default is 2/3.
			#' 					All "KK" designs require "lambda", the quantile cutoff of the subject distance distribution for determining matches. If unspecified, default is 10%.
			#' 					All "KK" designs require "t_0_pct", the percentage of total sample size n where matching begins. If unspecified, default is 35%.
			#' 					All "KK" designs have optional flag KK_verbose with default \code{FALSE} which prints out debug messages about how the matching-on-the-fly is working.
			#' 					All "KK21" designs further require "num_boot" which is the number of bootstrap samples taken to approximate the subject-distance distribution. 
			#' 					If unspecified, default is 500. There is an optional flag "proportion_use_speedup = TRUE" which uses a continuous regression on log(y/(1-y))
			#' 					instead of a beta regression each time to generate the weights in KK21 designs. The default is this flag is on.
			#' @return A new `SeqDesign` object.
			#' 
			#' @examples
			#' seq_des = SeqDesign$new(design = "KK21stepwise", response_type = "continuous")
			#'  
			initialize = function(
					n, 
					design, 
					response_type, 
					prob_T = 0.5, 
					include_is_missing_as_a_new_feature = TRUE, 
					verbose = TRUE, 
					...
				) {
				assertChoice(design, c("CRD", "iBCRD", "Efron", "Atkinson", "KK14", "KK21", "KK21stepwise"))
				assertChoice(response_type, c("continuous", "incidence", "proportion", "count", "survival"))
				assertCount(n, positive = TRUE)
				assertNumeric(prob_T, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
				assertFlag(include_is_missing_as_a_new_feature)
				assertFlag(verbose)
				
				
				private$isKK = grepl("KK", design)				

				if (design == "iBCRD" & !all.equal(n * prob_T, as.integer(n * prob_T), check.attributes = FALSE)){
					stop("Design iBCRD requires that the fraction of treatments of the total sample size must be a natural number.")
				}
				
				self$prob_T = prob_T
				private$n = n				
				self$response_type = response_type	
				self$design = design
				private$include_is_missing_as_a_new_feature = include_is_missing_as_a_new_feature
				private$verbose = verbose
				
				self$y = 	array(NA, n)
				self$w = 	array(NA, n)
				self$dead = array(NA, n)
				
				#now deal with design-specific hyperparameters
				private$other_params = list(...)
				if (design == "Efron"){
					if (is.null(private$other_params$weighted_coin_prob)){
						private$other_params$weighted_coin_prob = 2 / 3 #default Efron coin
					} else {
						assertNumeric(private$other_params$weighted_coin_prob, lower = 0, upper = 1)
					}
				}
				if (private$isKK){
					private$match_indic = array(NA, n)
					
					if (is.null(private$other_params$lambda)){
						private$other_params$lambda = 0.1 #default
					} else {
						assertNumeric(private$other_params$lambda, lower = 0, upper = 1)
					}
					if (is.null(private$other_params$t_0_pct)){
						private$t_0 = round(0.35 * n) #default
					} else {
						assertNumeric(private$other_params$t_0_pct, lower = .Machine$double.eps, upper = 1)
						private$t_0 = round(private$other_params$t_0_pct * n)
					}
					if (is.null(private$other_params$KK_verbose)){
						private$other_params$KK_verbose = FALSE
					} else {
						assertFlag(private$other_params$KK_verbose)
					}
					if (grepl("KK21", design)){
						private$isKK21 = TRUE
						if (is.null(private$other_params$num_boot)){
							private$other_params$num_boot = 500
						} else {
							assertCount(private$other_params$num_boot, positive = TRUE)
						}
						if (is.null(private$other_params$proportion_use_speedup)){
							private$other_params$proportion_use_speedup = TRUE							
						} else {
							assertFlag(private$other_params$proportion_use_speedup)	
						}
					}
				}
				if (private$verbose){
					cat(paste0("Intialized a ", design, " experiment with response type ", response_type, ".\n"))
				}					
			},
			
			#' @description
			#' Add subject-specific measurements for the next subject entrant and return this new subject's treatment assignment
			#' 
			#' @param x_new 			A row of the data frame corresponding to the new subject to be added (must be type data.table).
			#' @param allow_new_cols	Should we allow new/different features than previously seen in previous subjects in the 
			#' 							new subject's covariates? Default is \code{TRUE}.
			#' @param KK_verbose		If \code{TRUE}, we will print out messages about the KK assignment. This is useful for understanding
			#' 							how the KK assignment is working
			#' 
			#' @examples
			#' seq_des = SeqDesign$new(n = 100, p = 10, design = "CRD", response_type = "continuous")
			#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
			#' 
			add_subject_to_experiment_and_assign = function(x_new, allow_new_cols = TRUE) {					
				assertDataFrame(x_new, nrows = 1)
				assertClass(x_new, "data.table")
				if (self$check_experiment_completed()){
					stop(paste("You cannot add any new subjects as all n =", private$n, "subjects have already been added."))
				}
				if (any(is.na(x_new)) & self$t == 0){
					x_new = x_new[which(!is.na(x_new))]
					if (!allow_new_cols){
						warning("There is missing data in the first subject's covariate value(s). Setting the flag allow_new_cols = FALSE will disallow additional subjects")
					}					
				}
				
				xnew_data_types = sapply(x_new, class)
				
				if ("ordered" %in% unlist(xnew_data_types)){
					stop("Ordered factor data type is not supported; please convert to either an unordered factor or numeric.")
				}
				
				if (self$t > 0){
					Xraw_data_types = sapply(self$Xraw, class)
					colnames_Xraw = names(self$Xraw)
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
									self$Xraw[, (new_Xraw_cols[j]) := switch(new_Xraw_col_types[j],
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
				
				#iterate t
				self$t = self$t + 1
				#add new subject's measurements to the raw data frame (there should be the same exact columns even if there are new ones introduced)
				self$Xraw = rbind(self$Xraw, x_new)
				
				#we only bother with imputation and model matrices if we have enough data otherwise it's a huge mess
				#thus, designs cannot utilize imputations nor model matrices until this condition is met
				#luckily, those are the designs implemented herein so we have complete control (if you are extending this package, you'll have to deal with this issue here)
				if (self$t > (ncol(self$Xraw) + 2) & !(self$design %in% c("CRD", "iBCRD", "Efron"))){ #we only need to impute if we need the X's to make the allocation decisions 
					private$covariate_impute_if_necessary_and_then_create_model_matrix()
				}
				
				#now make the assignment
				self$w[self$t] = private[[paste0("assign_wt_", self$design)]]() #equivalent to do.call(what = paste0("assign_wt_", self$design), args = list())
				#and return the new assignment
				self$w[self$t]
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
				cat("Subject number", self$t, "is assigned to", ifelse(self$w[self$t] == 1, "TREATMENT", "CONTROL"), "via design", self$design, "\n")
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
				if (t > self$t){
					stop(paste("You cannot add response for subject", t, "when the most recent subjects' record added is", self$t))	
				}
				
				if (!is.na(self$y[t])){
					warning(paste("Overwriting previous response for t =", t, "y[t] =", self$y[t]))
				}
				
				#deal with the myriad checks on the response value based on response_type
				if (self$response_type == "continuous"){
					assertNumeric(y, any.missing = FALSE)	
				} else if (self$response_type == "incidence"){
					assertChoice(y, c(0, 1))
				} else if (self$response_type == "proportion"){
					assertNumeric(y, any.missing = FALSE, lower = 0, upper = 1) #the betareg package can handle 0's and 1's exactly (see their documentation) this seems to be why the statmod / numDeriv packages is also required 
				} else if (self$response_type == "count"){
					assertCount(y, na.ok = FALSE)
				} else if (self$response_type == "survival"){
					assertNumeric(y, any.missing = FALSE, lower = 0)
					if (y == 0){
						warning("0 survival responses not allowed --- recording .Machine$double.eps as its value instead")
						y = .Machine$double.eps
					}						
				}
				#finally, record the response value and the time at which it was recorded
				self$y[t] = y
				self$dead[t] = dead
				private$y_i_t_i[[t]] = self$t
			},
			
			#' @description
			#' For non-CARA designs, add all subject responses
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
					deads = rep(1, self$t)
				}
				assertNumeric(ys, len = self$t)
				assertNumeric(deads, len = self$t)
				
				for (t in 1 : self$t){
					self$add_subject_response(t, ys[t], deads[t])
				}
			},
			
			#' @description
			#' For KK designs only, this returns a list with useful matching statistics.
			#' 
			#' @return 	A list with the following data: \code{num_matches}, \code{prop_subjects_matched}, 
			#' 			\code{num_subjects_remaining_in_reservoir}, \code{prop_subjects_remaining_in_reservoir}.
			#' 
			#' @examples
			#' seq_des = SeqDesign$new(n = 6, p = 10, design = "KK14", response_type = "continuous")
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
			#' seq_des$matching_statistics()
			#'
			matching_statistics = function(){
				if (!private$isKK){
					stop("Matching statistics are only available for KK designs")
				}
				if (self$t == 0){
					stop("The experiment has not begun yet")
				}
				num_subjects_matched = sum(private$match_indic != 0, na.rm = TRUE)
				num_subjects_remaining_in_reservoir = self$t - num_subjects_matched
				list(
					num_matches = length(unique(private$match_indic[private$match_indic != 0])) ,
					prop_subjects_matched = num_subjects_matched / self$t,
					num_subjects_remaining_in_reservoir = num_subjects_remaining_in_reservoir,
					prop_subjects_remaining_in_reservoir = num_subjects_remaining_in_reservoir / self$t
				)					
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
				if (private$all_assignments_not_yet_allocated()){
					stop("This experiment is incomplete as the assignments aren't all wtated yet.")
				}
				if (private$all_responses_not_yet_recorded()){
					stop("This experiment is incomplete as the responses aren't all recorded yet.")
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
				if (private$all_assignments_not_yet_allocated()){
					FALSE
				} else if (private$all_responses_not_yet_recorded()){
					FALSE
				} else {
					TRUE
				}
			}
	),
	
	
	private = list(
		n = NULL,	
		include_is_missing_as_a_new_feature = NULL,
		verbose = NULL,
		
		#when during the experiment the subjects recorded
		y_i_t_i = list(),			
		#design specific parameters
		isKK = FALSE,
		isKK21 = FALSE,		
		other_params = NULL,
		t_0 = NULL,
		
		match_indic = NULL,
		
		duplicate = function(){
			self$assert_experiment_completed() #can't duplicate without the experiment being done
			d = SeqDesign$new(
				design = self$design, 
				response_type = self$response_type, 
				n = private$n, 
				prob_T = self$prob_T, 
				include_is_missing_as_a_new_feature = private$include_is_missing_as_a_new_feature, 
				verbose = FALSE
			)
			#we are assuming the experiment is complete so we have self$t = 0 initialized
			d$Xraw = self$Xraw
			d$Ximp = self$Ximp
			d$y = self$y
			d$dead = self$dead
			d$w = self$w		
			d$X = self$X		
			d$.__enclos_env__$private$y_i_t_i = private$y_i_t_i
			d$.__enclos_env__$private$isKK = private$isKK
			d$.__enclos_env__$private$isKK21 = private$isKK21
			d$.__enclos_env__$private$other_params = private$other_params
			d$.__enclos_env__$private$t_0 = private$t_0
			d$.__enclos_env__$private$match_indic = private$match_indic
			d
		},	
		
		covariate_impute_if_necessary_and_then_create_model_matrix = function(){
			#make a copy... sometimes the raw will be the same as the imputed if there are no imputations
			self$Ximp = copy(self$Xraw)
			
			column_has_missingness = self$Xraw[, lapply(.SD, function(xj) sum(is.na(xj)))] > 0
			if (any(column_has_missingness)){
				#deal with include_is_missing_as_a_new_feature here
				if (private$include_is_missing_as_a_new_feature){
					for (j in which(column_has_missingness)){
						self$Ximp = cbind(self$Ximp, ifelse(is.na(self$Ximp[, .SD, .SDcols = j]), 1, 0))
						names(self$Ximp)[ncol(self$Ximp)] = paste0(names(self$Ximp)[j], "_is_missing")
					}
				}
				
				#we need to convert characters into factor for the imputation to work
				col_types = self$Ximp[, lapply(.SD, function(xj){class(xj)})]
				idx_cols_to_convert_to_factor = which(col_types == "character")
				self$Ximp[, (idx_cols_to_convert_to_factor) := lapply(.SD, as.factor), .SDcols = idx_cols_to_convert_to_factor]
				
				#now do the imputation here by using missRanger (fast but fragile) and if that fails, use missForest (slow but more robust)
				self$Ximp = tryCatch({
							if (any(!is.na(self$y))){
								suppressWarnings(missRanger(cbind(self$Ximp, self$y[1 : nrow(self$Ximp)]), verbose = FALSE)[, 1 : ncol(self$Ximp)])
							} else {
								suppressWarnings(missRanger(self$Ximp, verbose = FALSE))
							}
						}, error = function(e){
							if (any(!is.na(self$y))){
								suppressWarnings(missForest(cbind(self$Ximp, self$y[1 : nrow(self$Ximp)]))$ximp[, 1 : ncol(self$Ximp)])
							} else {
								suppressWarnings(missForest(self$Ximp)$ximp)
							}
						}
				)
			}
			
			#now let's drop any columns that don't have any variation
			num_unique_values_per_column = self$Ximp[, lapply(.SD, function(xj){uniqueN(xj)})]
			self$Ximp = self$Ximp[, .SD, .SDcols = which(num_unique_values_per_column > 1)]
			
			#for nonblank data frames...
			if (ncol(self$Ximp) > 0){
				#now we need to convert character features into factors
				col_types = self$Ximp[, lapply(.SD, function(xj){class(xj)})]
				idx_cols_to_convert_to_factor = which(col_types == "character")
				self$Ximp[, (idx_cols_to_convert_to_factor) := lapply(.SD, as.factor), .SDcols = idx_cols_to_convert_to_factor]
				
				#now we need to update the numeric model matrix which may have expanded due to new factors, new missingness cols, etc
				self$X = model.matrix(~ ., self$Ximp)
				self$X = private$drop_linearly_dependent_cols(self$X)$M
				
				if (nrow(self$X) != nrow(self$Xraw) | nrow(self$X) != nrow(self$Ximp) | nrow(self$Ximp) != nrow(self$Xraw)){
					stop("boom")
				}
			} else {
				#blank covariate matrix
				self$X = matrix(NA, nrow = nrow(self$Xraw), ncol = 0)
			}			
		},
		
		compute_all_subject_data = function(){
			i_present_y = which(!is.na(self$y))
			i_all = 1 : self$t
			i_past = 	if (self$t == 1){
							c()
						} else {
							1 : (self$t - 1)
						}
			i_past_y_present = intersect(i_past, i_present_y)
			i_all_y_present = intersect(i_all, i_present_y)
			
			if (self$t == 1){
				xt = self$X[self$t, ]
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
					w_all = self$w[self$t],
					
					Xall_with_y = Xall_with_y,
					rank_all_with_y = ifelse(nrow(Xall_with_y) == 1, length(xt), NA),
					xt_all = ifelse(nrow(Xall_with_y) == 1, xt, NA),
					w_all = self$w[self$t],
					
					y_previous = NA,
					y_all = self$y[i_all_y_present],
					
					dead_previous = NA,
					dead_all = self$dead[i_all_y_present]
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
					w_prev = self$w[i_past],
					
					X_prev_with_y = past_info_with_y$Xint,
					rank_prev_with_y = past_info_with_y$rank,
					xt_prev_with_y = past_info_with_y$xt,
					w_prev_with_y = self$w[i_past_y_present],
					
					X_all = all_info$Xint,
					rank_all = all_info$rank,
					xt_all = all_info$xt,
					w_all = self$w[i_all],
					
					X_all_with_y = all_info_with_y$Xint,
					rank_all_with_y = all_info_with_y$rank,
					xt_all_with_y = all_info_with_y$xt,
					w_all_with_y = self$w[i_all_y_present],
					
					X_all_scaled = all_info_scaled$Xint,
					rank_all_scaled = all_info_scaled$rank,
					xt_all_scaled = all_info_scaled$xt,
					w_all_scaled = self$w[i_all],
					
					X_all_with_y_scaled = all_info_with_y_scaled$Xint,
					rank_all_with_y_scaled = all_info_with_y_scaled$rank,
					xt_all_with_y_scaled = all_info_with_y_scaled$xt,
					w_all_with_y_scaled = self$w[i_all_y_present],
					
					y_prev = self$y[i_past_y_present],
					y_all = self$y[i_all_y_present],
					
					dead_prev = self$dead[i_past_y_present],
					dead_all = self$dead[i_all_y_present]
				)
			}		
		},
		
		remove_linearly_dependent_covariates_and_compute_info = function(is, scaled = FALSE){
			if (length(is) == 0){				
				list(Xint = NA, rank = NA, xt = NA)
			} else {
				Xint = self$X[is, , drop = FALSE]
				
				#we first kill columns that have no variation
				num_unique_values_per_column = apply(Xint, 2, function(xj){length(unique(xj))})
				js = which(num_unique_values_per_column > 1)
				Xint = Xint[, js]
				xt = self$X[self$t, js]
				
				if (scaled & length(is) > 1){
					#we need to scale all data including the tth
					Xint = rbind(Xint, xt)
					Xint = apply(Xint, 2, scale)
					#if the column is all one unique value, it goes NaN after scaling, so we have to patch that up
					Xint[is.nan(Xint)] = 0
				}
				
				drop_obj = private$drop_linearly_dependent_cols(Xint)
				Xint = drop_obj$M
				xt = xt[drop_obj$js] #make sure xt comports with Xint!
				rank = Matrix::rankMatrix(Xint)
				
				if (scaled & length(is) > 1){
					list(Xint = Xint[1 : (nrow(Xint) - 1), , drop = FALSE], rank = rank, xt = Xint[nrow(Xint), ])
				} else {
					list(Xint = Xint, rank = rank, xt = xt)	
				}								
			}
		},
		
		drop_linearly_dependent_cols = function(M){
			rank = Matrix::rankMatrix(M)
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
		
		redraw_w_according_to_design = function(){
			#the designs are implemented here for all n for speed
			if (self$design == "CRD"){
				self$w = rbinom(private$n, 1, self$prob_T)
			} else if (self$design == "iBCRD"){
				n_T_total = round(private$n * self$prob_T) #this quantity should never be a fraction anyway as it was checked during initialization
				self$w = sample(c(rep(1, n_T_total), rep(0, private$n - n_T_total)))
			} else if (self$design == "Efron"){
				for (t in 1 : private$n){
					self$w[t] = private$assign_wt_Efron()
				}
			} else if (self$design == "Atkinson"){
				self$t = 0 #reset
				for (t in 1 : private$n){
					self$t = self$t + 1
					self$w[t] = private$assign_wt_Atkinson()
				}				
			} else { #KK design
				m_vec = private$match_indic
				#we now rearrange within each match set (and the reservoir)
				for (m in 0 : max(m_vec)){
					self$w[m_vec == m] = sample(self$w[m_vec == m]) 
				}
			}				
			self$t = private$n #mark the experiment as complete
		},
		
		all_responses_not_yet_recorded = function(){
			sum(!is.na(self$y)) != private$n
		},
		
		all_assignments_not_yet_allocated = function(){
			self$t != private$n
		},
		
		assign_wt_CRD = function(){
			rbinom(1, 1, self$prob_T)
		},
		
		assign_wt_iBCRD = function(){
			n_T_total = round(private$n * self$prob_T) #this quantity should never be a fraction anyway as it was checked during initialization
			nT = sum(self$w == 1, na.rm = TRUE)
			nC = sum(self$w == 0, na.rm = TRUE)
			sample(c(rep(1, n_T_total - nT), rep(0, n_T_total - nC)), 1)
		},
		
		assign_wt_Efron = function(){
			n_T = sum(private$w, na.rm = TRUE)
			n_C = private$n - n_T
			if (n_T * self$prob_T > n_C * (1 - self$prob_T)){
				rbinom(1, 1, 1 - private$other_params$weighted_coin_prob)
			} else if (n_T * self$prob_T < n_C * (1 - self$prob_T)){
				rbinom(1, 1, private$other_params$weighted_coin_prob)
			} else {
				private$assign_wt_CRD()
			}				
		},
		
		assign_wt_Atkinson = function(){
			
			#if it's too early in the trial or if all the assignments are the same, then randomize
			if (self$t <= ncol(self$Xraw) + 2 + 1 | length(unique(self$t)) == 1){
				private$assign_wt_CRD()
			} else {
				all_subject_data = private$compute_all_subject_data()
				#this matrix is [w | 1 | X]
				Xprev_with_w = cbind(self$w[1 : (self$t - 1)], 1, all_subject_data$X_prev)
				XwtXw = t(Xprev_with_w) %*% Xprev_with_w	
				tryCatch({										
					M = (self$t - 1) * solve(XwtXw, tol = .Machine$double.xmin)
					A = M[1, 2 : (all_subject_data$rank_prev + 2)] %*% c(1, all_subject_data$xt_prev) 
					s_over_A_plus_one_sq = (M[1, 1] / A + 1)^2
					#assign via the Atkinson weighted biased coin
					rbinom(1, 1, s_over_A_plus_one_sq / (s_over_A_plus_one_sq + 1))	
				}, error = function(e){ #sometimes XtX is still not invertible... 
					#so in that case... just randomize
					private$assign_wt_CRD()
				})
			}
		},
		
		assign_wt_KK14 = function(){
			wt = 	if (sum(private$match_indic == 0, na.rm = TRUE) == 0 | self$t <= (ncol(self$Xraw) + 2) | self$t <= private$t_0){
						#we're early, so randomize
						private$match_indic[self$t] = 0
						private$assign_wt_CRD()
					} else {
						all_subject_data = private$compute_all_subject_data()
						# cat("else\n")
						#first calculate the threshold we're operating at	
						#when inverting, ensure full rank by adding eps * I			
						S_xs_inv = solve(var(all_subject_data$X_prev) + diag(.Machine$double.eps, all_subject_data$rank_prev), tol = .Machine$double.xmin)
						F_crit =  qf(private$other_params$lambda, all_subject_data$rank_prev, self$t - all_subject_data$rank_prev)
						T_cutoff_sq = all_subject_data$rank_prev * (private$n - 1) / (private$n - all_subject_data$rank_prev) * F_crit
						#now iterate over all items in reservoir and take the minimum distance x
						reservoir_indices = which(private$match_indic == 0)
						sqd_distances_times_two = array(NA, length(reservoir_indices))
						for (r in 1 : length(reservoir_indices)){
							x_r_x_new_delta = all_subject_data$xt_prev - all_subject_data$X_prev[reservoir_indices[r], ]
							sqd_distances_times_two[r] = t(x_r_x_new_delta) %*% S_xs_inv %*% x_r_x_new_delta		
						}					
						#find minimum distance index
						min_sqd_dist_index = which(sqd_distances_times_two == min(sqd_distances_times_two))
						if (length(sqd_distances_times_two[min_sqd_dist_index]) > 1 || length(T_cutoff_sq) > 1){
							min_sqd_dist_index = min_sqd_dist_index[1] #if there's a tie, just take the first one
						}
						#if it's smaller than the threshold, we're in business: match it
						if (sqd_distances_times_two[min_sqd_dist_index] < T_cutoff_sq){
							match_num = max(private$match_indic, na.rm = TRUE) + 1
							private$match_indic[reservoir_indices[min_sqd_dist_index]] = match_num
							private$match_indic[self$t] = match_num
							#assign opposite
							1 - self$w[reservoir_indices[min_sqd_dist_index]]
						} else { #otherwise, randomize and add it to the reservoir
							private$match_indic[self$t] = 0
							private$assign_wt_CRD()
						}
					}
			if (is.na(private$match_indic[self$t])){
				stop("no match data recorded")
			}
			wt
		},
		
		compute_weight_KK21_continuous = function(xs_to_date, ys_to_date, deaths_to_date, j){
			ols_mod = lm(ys_to_date ~ xs_to_date[, j])
			summary_ols_mod = suppressWarnings(coef(summary(ols_mod)))
			#1 - coef(summary(logistic_regr_mod))[2, 4]
			#if there was only one row, then this feature was all one unique value... so send back a weight of nada
			ifelse(nrow(summary_ols_mod) >= 2, abs(summary_ols_mod[2, 3]), .Machine$double.eps)
		},
		
		compute_weight_KK21_incidence = function(xs_to_date, ys_to_date, deaths_to_date, j){
			tryCatch({
				logistic_regr_mod = suppressWarnings(glm(ys_to_date ~ xs_to_date[, j], family = "binomial"))
			}, error = function(e){return(.Machine$double.eps)}) #sometimes these glm's blow up and we don't really care that much
			summary_logistic_regr_mod = coef(summary_glm_lean(logistic_regr_mod))
			#1 - coef(summary(logistic_regr_mod))[2, 4]
			#if there was only one row, then this feature was all one unique value... so send back a weight of nada
			ifelse(nrow(summary_logistic_regr_mod) >= 2, abs(summary_logistic_regr_mod[2, 3]), .Machine$double.eps)
		},
		
		compute_weight_KK21_count = function(xs_to_date, ys_to_date, deaths_to_date, j){
			tryCatch({
				negbin_regr_mod = suppressWarnings(MASS::glm.nb(y ~ x, data = data.frame(x = xs_to_date[, j], y = ys_to_date)))
				summary_negbin_regr_mod = coef(summary_glm_lean(negbin_regr_mod))
				#1 - coef(summary(negbin_regr_mod))[2, 4]
				#if there was only one row, then this feature was all one unique value... so send back a weight of nada
				return(ifelse(nrow(summary_negbin_regr_mod) >= 2, abs(summary_negbin_regr_mod[2, 3]), .Machine$double.eps))
			}, error = function(e){}) #sometimes these glm's blow up and we don't really care that much
			.Machine$double.eps #otherwise weight it nada
		},
		
		compute_weight_KK21_proportion = function(xs_to_date, ys_to_date, deaths_to_date, j){
			if (!private$other_params$proportion_use_speedup){
				tryCatch({
					beta_regr_mod = suppressWarnings(betareg::betareg(y ~ x, data = data.frame(x = xs_to_date[, j], y = ys_to_date)))
					summary_beta_regr_mod = coef(summary(beta_regr_mod)) 
					#1 - coef(summary(beta_regr_mod))$mean[2, 4]
					tab = 	if (!is.null(summary_beta_regr_mod$mean)){ #beta model
								summary_beta_regr_mod$mean 
							} else if (!is.null(summary_beta_regr_mod$mu)){ #extended-support xbetax model
								summary_beta_regr_mod$mu
							}
					if (nrow(tab) >= 2){
						return(abs(tab[2, 3]))
					}
				}, error = function(e){}) #sometimes these glm's blow up and we don't really care that much
			}
			#if that didn't work, let's just use the continuous weights on a transformed proportion
			ys_to_date[ys_to_date == 0] = .Machine$double.eps
			ys_to_date[ys_to_date == 1] = 1 - .Machine$double.eps
			private$compute_weight_KK21_continuous(xs_to_date, log(ys_to_date / (1 - ys_to_date)), deaths_to_date, j)
		},
		
		compute_weight_KK21_survival = function(xs_to_date, ys_to_date, deaths_to_date, j){
			surv_obj = survival::Surv(ys_to_date, deaths_to_date)
			#sometims the weibull is unstable... so try other distributions... this doesn't matter since we are just trying to get weights
			#and we are not relying on the model assumptions
			for (dist in c("weibull", "lognormal", "loglogistic")){
				surv_regr_mod = robust_survreg_with_surv_object(surv_obj, xs_to_date[, j], dist = dist)
				if (is.null(surv_regr_mod)){
					break
				}
				summary_surv_regr_mod = suppressWarnings(summary(surv_regr_mod)$table)
				if (any(is.nan(summary_surv_regr_mod))){
					break
				}
				weight = ifelse(nrow(summary_surv_regr_mod) >= 2, abs(summary_surv_regr_mod[2, 3]), NA)
				#1 - summary(weibull_regr_mod)$table[2, 4]
				if (!is.na(weight)){ #sometimes these glm's blow up and we don't really care that much
					return(weight)
				}
			}
			#if that didn't work, default to OLS and log the survival times... again... this doesn't matter since we are just trying to get weights
			#and we are not relying on the model assumptions being true
			private$compute_weight_KK21_continuous(xs_to_date, log(ys_to_date), deaths_to_date, j)
		},
				
		assign_wt_KK21 = function(){
			wt = 	if (sum(private$match_indic == 0, na.rm = TRUE) == 0 | (self$t <= ncol(self$Xraw) + 2) | self$t <= private$t_0){
						#we're early or the reservoir is empty, so randomize
						#cat("    assign_wt_KK21 CRD t", self$t, "\n")
						private$match_indic[self$t] = 0
						private$assign_wt_CRD()
					} else if (is.null(self$X) | (sum(!is.na(self$y)) < 2 * (ncol(self$X) + 2))){ 
						#This is the number of responses collected before
						#the algorithm begins estimating the covariate-specific weights. If left unspecified this defaults to \code{2 * (p + 2)} i.e. two data points
						#for every glm regression parameter estimated (p covariates, the intercept and the coefficient of the additive treatment effect). The minimum
						#value is p + 2 to allow OLS to estimate. If this number of not met, we default to KK14 matching (which is better than nothing as it matches
						#the covariate distributions at the very least).				
						private$assign_wt_KK14()
					}  else {
						all_subject_data = private$compute_all_subject_data()
						#1) need to calculate the weights - slightly different for each response type
						#we calculate with the scaled covariates putting them all on the same footing
						#and we only calculate with those subjects that actually have y values
						raw_weights = array(NA, all_subject_data$rank_all_with_y_scaled)
						for (j in 1 : all_subject_data$rank_all_with_y_scaled){
							raw_weights[j] = private[[paste0("compute_weight_KK21_", self$response_type)]](
								all_subject_data$X_all_with_y_scaled, 
								all_subject_data$y_all, 
								all_subject_data$dead_all, 
								j
							) 
						}
						if (any(is.na(raw_weights)) | any(is.infinite(raw_weights)) | any(is.nan(raw_weights)) | any(raw_weights < 0)){
							stop("raw weight values illegal KK21")
						}
						# (c) now we need to normalize the weights
						#cat("    assign_wt_KK21 using weights t", self$t, "raw weights", weights, "\n")
						self$covariate_weights = raw_weights / sum(raw_weights)
						names(self$covariate_weights) = colnames(all_subject_data$X_all_with_y_scaled)

						#cat("    assign_wt_KK21 using weights t", self$t, "sorted weights", sort(weights), "\n")
						
						#2) now iterate over all items in reservoir and calculate the weighted sqd distiance vs new guy 
						reservoir_indices = which(private$match_indic == 0)
						weighted_features = colnames(all_subject_data$X_all_with_y_scaled)
						x_new = all_subject_data$xt_all_scaled[weighted_features]
						X_all_scaled_col_subset = all_subject_data$X_all_scaled[, weighted_features, drop = FALSE]
						
						weighted_sqd_distances = array(NA, length(reservoir_indices))
						for (r in 1 : length(reservoir_indices)){
							x_r_x_new_delta = x_new - X_all_scaled_col_subset[reservoir_indices[r], ]
							weighted_sqd_distances[r] = x_r_x_new_delta^2 %*% self$covariate_weights
						}
						#3) find minimum weighted sqd distiance index
						min_weighted_sqd_dist_index = which(weighted_sqd_distances == min(weighted_sqd_distances))
						
						#generate a cutoff for the weighted minimum distance squared based on bootstrap
						bootstrapped_weighted_sqd_distances = array(NA, private$other_params$num_boot)
						for (b in 1 : private$other_params$num_boot){
							two_xs  = X_all_scaled_col_subset[sample.int(self$t, 2), ] #self$X[sample_int_ccrank(self$t, 2, rep(1, (self$t))), ] #
							delta_x = two_xs[1, ] - two_xs[2, ]
							bootstrapped_weighted_sqd_distances[b] = delta_x^2 %*% self$covariate_weights
						}
						
						min_weighted_dsqd_cutoff_sq = quantile(bootstrapped_weighted_sqd_distances, private$other_params$lambda)
						
						#5) Now, does the minimum make the cut?
						if (length(weighted_sqd_distances[min_weighted_sqd_dist_index]) > 1 || length(min_weighted_dsqd_cutoff_sq) > 1){
							min_weighted_sqd_dist_index = min_weighted_sqd_dist_index[1] #if there's a tie, just take the first one
						}
						#  (a) if it's smaller than the threshold, we're in business: match it
						if (weighted_sqd_distances[min_weighted_sqd_dist_index] < min_weighted_dsqd_cutoff_sq){
							match_num = max(private$match_indic, na.rm = TRUE) + 1
							private$match_indic[reservoir_indices[min_weighted_sqd_dist_index]] = match_num
							private$match_indic[self$t] = match_num
							#assign opposite
							1 - self$w[reservoir_indices[min_weighted_sqd_dist_index]]
						# (b) otherwise, randomize and add it to the reservoir
						} else { 
							private$match_indic[self$t] = 0	
							private$assign_wt_CRD()
						}
					}
			if (is.na(private$match_indic[self$t])){
				stop("no match data recorded")
			}
			wt
		},
		
		compute_weights_KK21stepwise = function(Xfull, response_obj, ws, abs_z_compute_fun){
			weights = array(NA, ncol(Xfull))
			j_droppeds = c()
			X_stepwise = matrix(NA, nrow = nrow(Xfull), ncol = 0)
			
			repeat {
				covs_to_try = setdiff(1 : ncol(Xfull), j_droppeds)				
				if (length(covs_to_try) == 0){ #if there's none left, we jet
					break
				}
				abs_approx_zs = array(NA, ncol(Xfull))
				for (j in covs_to_try){
					abs_approx_zs[j] = abs_z_compute_fun(response_obj, cbind(Xfull[, j], X_stepwise, ws))
				}
				j_max = which.max(abs_approx_zs)
				weights[j_max] = abs_approx_zs[j_max]
				j_droppeds = c(j_droppeds, j_max)
				X_stepwise = cbind(X_stepwise, Xfull[, j_max])
			}
			if (any(is.na(weights))){
				stop("boom")					
			}
			weights
		},
		
		compute_weights_KK21stepwise_continuous = function(xs, ys, ws, ...){
			private$compute_weights_KK21stepwise(xs, scale(ys), ws, function(response_obj, covariate_data_matrix){
				ols_mod = lm(response_obj ~ covariate_data_matrix)
				abs(coef(suppressWarnings(summary(ols_mod)))[2, 3])
			})
		},
		
		compute_weights_KK21stepwise_incidence = function(xs, ys, ws, ...){
			private$compute_weights_KK21stepwise(xs, ys, ws, function(response_obj, covariate_data_matrix){
				logistic_regr_mod = suppressWarnings(glm(response_obj ~ covariate_data_matrix, family = "binomial"))
				abs(coef(summary_glm_lean(logistic_regr_mod))[2, 3])
			})
		},
		
		compute_weights_KK21stepwise_count = function(xs, ys, ws, ...){	
			private$compute_weights_KK21stepwise(xs, ys, ws, function(response_obj, covariate_data_matrix){
				negbin_regr_mod = robust_negbinreg(response_obj ~ ., cbind(data.frame(response_obj = response_obj), covariate_data_matrix))
				abs(coef(summary_glm_lean(negbin_regr_mod))[2, 3])
			})	
		},
		
		compute_weights_KK21stepwise_proportion = function(xs, ys, ws, ...){
			if (!private$other_params$proportion_use_speedup){
				tryCatch({
					weight = 	private$compute_weights_KK21stepwise(xs, ys, ws, function(response_obj, covariate_data_matrix){					
									beta_regr_mod = suppressWarnings(betareg::betareg(response_obj ~ ., data = cbind(data.frame(response_obj = response_obj), covariate_data_matrix)))
									summary_beta_regr_mod = coef(summary(beta_regr_mod)) 
									tab = 	if (!is.null(summary_beta_regr_mod$mean)){ #beta model
												summary_beta_regr_mod$mean 
											} else if (!is.null(summary_beta_regr_mod$mu)){ #extended-support xbetax model
												summary_beta_regr_mod$mu
											}
									ifelse(nrow(tab) >= 2, abs(tab[2, 3]), NA)
								})
					if (!is.na(weight)){
						return(weight)
					}
				}, error = function(e){})
			}
			#if that didn't work, let's just use the continuous weights on a transformed proportion
			ys[ys == 0] = .Machine$double.eps
			ys[ys == 1] = 1 - .Machine$double.eps
			private$compute_weights_KK21stepwise_continuous(xs, log(ys / (1 - ys)), ws, ...)
		},
		
		compute_weights_KK21stepwise_survival = function(xs, ys, ws, deaths){		
			private$compute_weights_KK21stepwise(xs, survival::Surv(ys, deaths), ws, function(response_obj, covariate_data_matrix){
				#sometims the weibull is unstable... so try other distributions... this doesn't matter since we are just trying to get weights
				#and we are not relying on the model assumptions
				for (dist in c("weibull", "lognormal", "loglogistic")){
					surv_regr_mod = robust_survreg_with_surv_object(response_obj, covariate_data_matrix, dist = dist)
					if (is.null(surv_regr_mod)){
						break
					}
					summary_surv_regr_mod = suppressWarnings(summary(surv_regr_mod)$table)
					if (any(is.nan(summary_surv_regr_mod))){
						break
					}
					weight = ifelse(nrow(summary_surv_regr_mod) >= 2, abs(summary_surv_regr_mod[2, 3]), NA)
					#1 - summary(weibull_regr_mod)$table[2, 4]
					if (!is.na(weight)){
						return(weight)
					}
				}	
				#if that didn't work, default to OLS and log the survival times... again... this doesn't matter since we are just trying to get weights
				#and we are not relying on the model assumptions
				ols_mod = lm(log(as.numeric(response_obj)[1 : length(response_obj)]) ~ covariate_data_matrix)
				abs(coef(suppressWarnings(summary(ols_mod)))[2, 3])
			})	
		},
		
		assign_wt_KK21stepwise = function(){
			wt = 	if (sum(private$match_indic == 0, na.rm = TRUE) == 0 | self$t <= (ncol(self$Xraw) + 2) | self$t <= private$t_0){
						#we're early or the reservoir is empty, so randomize
						#cat("    assign_wt_KK21stepwise CRD t", self$t, "\n")
						private$match_indic[self$t] = 0
						private$assign_wt_CRD()
					} else if (is.null(self$X) | (sum(!is.na(self$y)) < 2 * (ncol(self$X) + 2))){ 
						#This is the number of responses collected before
						#the algorithm begins estimating the covariate-specific weights. If left unspecified this defaults to \code{2 * (p + 2)} i.e. two data points
						#for every glm regression parameter estimated (p covariates, the intercept and the coefficient of the additive treatment effect). The minimum
						#value is p + 2 to allow OLS to estimate. If this number of not met, we default to KK14 matching (which is better than nothing as it matches
						#the covariate distributions at the very least).	
						private$assign_wt_KK14()
					} else {		
						all_subject_data = private$compute_all_subject_data()		
						#1) need to calculate the weights...
						
						#we need to now run the appropriate (based on response_type) stepwise procedure to get the weights
						raw_weights = private[[paste0("compute_weights_KK21stepwise_", self$response_type)]](
							all_subject_data$X_all_with_y_scaled, #to calculate weights, we need to use only the data that has y's!
							all_subject_data$y_all,
							all_subject_data$w_all_with_y_scaled, 
							all_subject_data$dead_all
						)						
						if (any(is.na(raw_weights)) | any(is.infinite(raw_weights)) | any(is.nan(raw_weights)) | any(raw_weights < 0)){
							stop("raw weight values illegal KK21stepwise")
						}
						#ensure the weights are normalized
						self$covariate_weights = raw_weights / sum(raw_weights)
						names(self$covariate_weights) = colnames(all_subject_data$X_all_with_y_scaled)
						#cat("    assign_wt_KK21stepwise using weights t", self$t, "weights", weights, "\n")
						
						#2) now iterate over all items in reservoir and calculate the weighted sqd distiance vs new guy 
						reservoir_indices = which(private$match_indic == 0)
						weighted_features = colnames(all_subject_data$X_all_with_y_scaled)
						x_new = all_subject_data$xt_all_scaled[weighted_features]
						X_all_scaled_col_subset = all_subject_data$X_all_scaled[, weighted_features, drop = FALSE]
						
						weighted_sqd_distances = array(NA, length(reservoir_indices))
						for (r in 1 : length(reservoir_indices)){
							x_r_x_new_delta = x_new - X_all_scaled_col_subset[reservoir_indices[r], ] #we have to use all data here
							weighted_sqd_distances[r] = x_r_x_new_delta^2 %*% self$covariate_weights			
						}
						#3) find minimum weighted sqd distiance index
						min_weighted_sqd_dist_index = which(weighted_sqd_distances == min(weighted_sqd_distances))
						
						#generate a cutoff for the weighted minimum distance squared based on bootstrap
						bootstrapped_weighted_sqd_distances = array(NA, private$other_params$num_boot)
						for (b in 1 : private$other_params$num_boot){
							two_xs  = X_all_scaled_col_subset[sample.int(self$t, 2), ] #self$X[sample_int_ccrank(self$t, 2, rep(1, (self$t))), ] 
							delta_x = two_xs[1, ] - two_xs[2, ]
							bootstrapped_weighted_sqd_distances[b] = delta_x^2 %*% self$covariate_weights
						}
						
						min_weighted_dsqd_cutoff_sq = quantile(bootstrapped_weighted_sqd_distances, private$other_params$lambda)
						
						#5) Now, does the minimum make the cut?
						if (length(weighted_sqd_distances[min_weighted_sqd_dist_index]) > 1 || length(min_weighted_dsqd_cutoff_sq) > 1){
							min_weighted_sqd_dist_index = min_weighted_sqd_dist_index[1] #if there's a tie, just take the first one
						}
						#  (a) if it's smaller than the threshold, we're in business: match it
						if (weighted_sqd_distances[min_weighted_sqd_dist_index] < min_weighted_dsqd_cutoff_sq){
							new_match_id = max(private$match_indic, na.rm = TRUE) + 1 #the ID of a new match
							private$match_indic[reservoir_indices[min_weighted_sqd_dist_index]] = new_match_id
							private$match_indic[self$t] = new_match_id
							1 - self$w[reservoir_indices[min_weighted_sqd_dist_index]]
						} else { # (b) otherwise, randomize and add it to the reservoir
							private$match_indic[self$t] = 0	
							private$assign_wt_CRD()	
						}
					}
			if (is.na(private$match_indic[self$t])){ #this should never happen
				stop("no match data recorded")
			}
			wt
		}
	)
)