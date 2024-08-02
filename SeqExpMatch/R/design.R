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
			#' @field X			A numeric matrix of subject data with number of rows n (the number of subjects) and number of 
			#' 					columns p (the number of characteristics measured for each subject). This matrix is filled in
			#' 					sequentially by the experimenter and thus will have data present for rows 1...t (i.e. the number of subjects in the
			#' 					experiment currently) but otherwise will be missing.	
			#' @field Ximp		Same as \code{X} except			
			#' @field y			A numeric vector of subject responses with number of entries n (the number of subjects). During
			#' 					the KK21 designs the experimenter fills these values in when they are measured.
			#' 					For non-KK21 designs, this vector can be set at anytime (but must be set before inference is desired).				
			#' @field dead		A binary vector of whether the subject is dead with number of entries n (the number of subjects). This 
			#' 					vector is filled in only for \code{response_type} values "survival". The value
			#' 					of 1 indicates uncensored (as the subject died) and a value 0 indicates the real survival value is censored 
			#' 					as the subject is still alive at the time of measurement. This follows the same convention as the \code{event} 
			#' 					argument in the canonical \code{survival} package in the data constructor found here: \link{\code{survival::Surv}}. During
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
			t = 0,
			design = NULL,
			X = NULL,
			Ximp = NULL,
			y = NULL,	
			dead = NULL,
			prob_T = NULL,
			w = NULL,
			response_type = NULL,
			#' 				
			#' @description
			#' Initialize a sequential experimental design
			#' 
			#' @param n 		Number of subjects fixed beforehand. A future version of this software will allow
			#' 					for sequential stopping and thus n will not need to be prespecified.
			#' @param p 		Number of characteristics measured for each subject. If measurement j are 
			#' 					categorical with L_j levels, you must select a reference level and convert this information
			#' 					to L_j-1 dummies. Thus p := # of numeric variables + sum_j (L_j - 1).
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
			#' 							This package will enforce that all added responses via \link{add_subject_response} will be
			#' 							of the appropriate type.
			#' @param prob_T	The probability of the treatment assignment. This defaults to \code{0.5}.
			#' @param verbose	A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}.
			#' @param ...		Design-specific parameters:
			#' 					"Efron" requires "weighted_coin_prob" which is the probability of the weighted coin for assignment. If unspecified, default is 2/3.
			#' 					All "KK" designs require "lambda", the quantile cutoff of the subject distance distribution for determining matches. If unspecified, default is 10%.
			#' 					All "KK" designs require "t_0_pct", the percentage of total sample size n where matching begins. If unspecified, default is 35%.
			#' 					All "KK21" designs further require "num_boot" which is the number of bootstrap samples taken to approximate the subject-distance distribution. 
			#' 					If unspecified, default is 500. They also require "num_responses_to_begin_weight_estimation". This is the number of responses collected before
			#' 					the algorithm begins estimating the covariate-specific weights. If left unspecified this defaults to \code{2 * (p + 2)} i.e. two data points
			#' 					for every glm regression parameter estimated (p covariates, the intercept and the coefficient of the additive treatment effect). The minimum
			#' 					value is p + 2 to allow OLS to estimate. If this number of not met, we default to KK14 matching (which is better than nothing as it matches
			#' 					the covariate distributions at the very least).
			#' @return A new `SeqDesign` object.
			#' 
			#' @examples
			#' seq_des = SeqDesign$new(n = 100, p = 10, design = "KK21stepwise", response_type = "continuous")
			#'  
			initialize = function(n, p, design, response_type, prob_T = 0.5, verbose = TRUE, ...) {
				assertCount(n, positive = TRUE)
				assertCount(p, positive = TRUE)
				assertNumeric(prob_T, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
				assertChoice(design, c("CRD", "iBCRD", "Efron", "Atkinson", "KK14", "KK21", "KK21stepwise"))
				assertChoice(response_type, c("continuous", "incidence", "proportion", "count", "survival"))
				assertFlag(verbose)
				
				if (!all.equal(n * prob_T, as.integer(n * prob_T), check.attributes = FALSE) & design == "iBCRD"){
					stop("Design iBCRD requires that the fraction of treatments of the total sample size must be a natural number.")
				}
				self$prob_T = prob_T
				private$n = n
				private$p = p					
				self$response_type = response_type	
				self$design = design
				private$verbose = verbose
				
				#initialize data derived from n, p
				self$X = matrix(NA, nrow = n, ncol = p)
				self$y = array(NA, n)
				self$w = array(NA, n)
				
				#now deal with design-specific hyperparameters
				other_params = list(...)
				if (design == "Efron"){
					if (is.null(other_params$weighted_coin_prob)){
						private$weighted_coin_prob = 2 / 3 #default Efron coin
					} else {
						assertNumeric(other_params$weighted_coin_prob, lower = 0, upper = 1)
						private$weighted_coin_prob = other_params$weighted_coin_prob
					}
				}
				if (grepl("KK", design)){
					private$isKK = TRUE
					private$match_indic = array(NA, n)
					if (is.null(other_params$lambda)){
						private$lambda = 0.1 #default
					} else {
						assertNumeric(other_params$lambda, lower = 0, upper = 1)
						private$lambda = other_params$lambda
					}
					if (is.null(other_params$t_0_pct)){
						private$t_0 = round(0.35 * n) #default
					} else {
						assertNumeric(other_params$t_0_pct, lower = (p + 2) / n, upper = 1)
						private$t_0 = round(other_params$t_0_pct * n)
					}
					if (grepl("KK21", design)){
						private$isKK21 = TRUE
						self$y = array(NA, n)
						if (is.null(other_params$num_boot)){
							private$num_boot = 500
						} else {
							assertCount(other_params$num_boot, positive = TRUE)
							private$num_boot = other_params$num_boot
						}
						if (is.null(other_params$num_responses_to_begin_weight_estimation)){
							private$num_responses_to_begin_weight_estimation = 2 * (p + 2)
						} else {
							assertCount(other_params$num_responses_to_begin_weight_estimation)
							assertNumeric(other_params$num_responses_to_begin_weight_estimation, lower = p + 2)
							private$num_responses_to_begin_weight_estimation = other_params$num_responses_to_begin_weight_estimation
						}
					}
				}
				if (private$verbose){
					cat(paste0("Intialized a ", design, " design for ", n, " subjects each with ", p, " characteristics measured per subject.\n"))
				}					
			},
			
			#' @description
			#' Add subject-specific measurements for the next subject entrant and return this new subject's treatment assignment
			#' 
			#' @param x_vec A p-length numeric vector
			#' 
			#' @examples
			#' seq_des = SeqDesign$new(n = 100, p = 10, design = "CRD", response_type = "continuous")
			#' seq_des$add_subject_to_experiment_and_assign(c(1, 38, 142, 71, 5.3, 0, 0, 0, 1, 0))
			#' 
			add_subject_to_experiment_and_assign = function(x_vec) {					
				assertNumeric(x_vec, len = private$p, any.missing = FALSE)
				if (self$check_experiment_completed()){
					stop(paste("You cannot add any new subjects as all n =", private$n, "have already been added."))
				}
				
				#iterate t
				self$t = self$t + 1
				#add subject to data frame
				self$X[self$t, ] = x_vec
				#now make the assignment
				self$w[self$t] = private[[paste0("assign_wt_", self$design)]]() #equivalent to do.call(what = paste0("assign_wt_", self$design), args = list())
				#return the new assignment
				self$w[self$t]
			},
			
			#' @description
			#' Prints the current assignment to screen. Should be called after \code{add_subject_to_experiment_and_assign}.
			#' 
			#' @examples
			#' seq_des = SeqDesign$new(n = 100, p = 10, design = "CRD", response_type = "continuous")
			#' 
			#' seq_des$add_subject_to_experiment_and_assign(c(1, 38, 142, 71, 5.3, 0, 0, 0, 1, 0))
			#' seq_des$print_current_subject_assignment()
			#' 
			print_current_subject_assignment = function(){
				cat("Subject number", self$t, "is assigned to", ifelse(self$w[self$t] == 1, "TREATMENT", "CONTROL"), "via design", self$design, "\n")
			},
			
			#' @description
			#' For CARA designs, add subject response for the a subject
			#' 
			#' @param t 	 The subject index for which to attach a response (beginning with 1, ending with n). You cannot add responses
			#' 				 for subjects that have not yet been added to the experiment via \link{add_subject_to_experiment_and_assign}
			#' @param y 	 The response value which must be appropriate for the response_type. 
			#' @param dead	 If the response is censored, enter 0 for this value. This is only necessary to specify for response type
			#' 				 "survival" otherwise do not specify this argument (as it will default to 1).
			#' @examples
			#' seq_des = SeqDesign$new(n = 100, p = 10, design = "KK21", response_type = "continuous")
			#' 
			#' seq_des$add_subject_to_experiment_and_assign(c(1, 38, 142, 71, 5.3, 0, 0, 0, 1, 0))
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
				
				assertCount(t)
				if (t > self$t){
					stop(paste("You cannot add response for subject", t, "when the most recent subject added is", self$t))	
				}
				
				#deal with the myriad checks on the response value based on response_type
				if (self$response_type == "continuous"){
					assertNumeric(y, any.missing = FALSE)	
				} else if (self$response_type == "incidence"){
					assertChoice(y, c(0, 1))
				} else if (self$response_type == "proportion"){
					assertNumeric(y, any.missing = FALSE, lower = 0, upper = 1)
					if (y == 0){
						warning("0% proportion responses not allowed --- recording .Machine$double.eps as its value instead")
						y = .Machine$double.eps
					} else if (y == 1){
						warning("100% proportion responses not allowed --- recording 1 - .Machine$double.eps as its value instead")
						y = 1 - .Machine$double.eps							
					}
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
			#' seq_des$add_subject_to_experiment_and_assign(c(1, 38, 142, 71, 5.3, 0, 0, 0, 1, 0))
			#' seq_des$add_subject_to_experiment_and_assign(c(0, 27, 127, 60, 5.5, 0, 0, 0, 1, 0))
			#' seq_des$add_subject_to_experiment_and_assign(c(1, 42, 169, 74, 5.1, 0, 1, 0, 0, 0))
			#' seq_des$add_subject_to_experiment_and_assign(c(0, 59, 105, 62, 5.9, 0, 0, 0, 1, 0))
			#' seq_des$add_subject_to_experiment_and_assign(c(1, 32, 186, 66, 5.6, 1, 0, 0, 0, 0))
			#' seq_des$add_subject_to_experiment_and_assign(c(1, 37, 178, 75, 6.5, 0, 0, 0, 0, 1))
			#' 
			#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
			#' 				
			add_all_subject_responses = function(ys, deads = NULL) {
				if (is.null(deads)){
					deads = rep(1, private$n)
				}
				assertNumeric(ys, len = private$n)
				assertNumeric(deads, len = private$n)
				
				for (t in 1 : private$n){
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
			#' seq_des$add_subject_to_experiment_and_assign(c(1, 38, 142, 71, 5.3, 0, 0, 0, 1, 0))
			#' seq_des$add_subject_to_experiment_and_assign(c(0, 27, 127, 60, 5.5, 0, 0, 0, 1, 0))
			#' seq_des$add_subject_to_experiment_and_assign(c(1, 42, 169, 74, 5.1, 0, 1, 0, 0, 0))
			#' seq_des$add_subject_to_experiment_and_assign(c(0, 59, 105, 62, 5.9, 0, 0, 0, 1, 0))
			#' seq_des$add_subject_to_experiment_and_assign(c(1, 32, 186, 66, 5.6, 1, 0, 0, 0, 0))
			#' seq_des$add_subject_to_experiment_and_assign(c(1, 37, 178, 75, 6.5, 0, 0, 0, 0, 1))
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
			#' seq_des$add_subject_to_experiment_and_assign(c(1, 38, 142, 71, 5.3, 0, 0, 0, 1, 0))
			#' 
			#' #if run, it would throw an error since all of the covariate vectors are not yet recorded
			#' #seq_des$assert_experiment_completed() 
			#' 
			#' seq_des$add_subject_to_experiment_and_assign(c(0, 27, 127, 60, 5.5, 0, 0, 0, 1, 0))
			#' seq_des$add_subject_to_experiment_and_assign(c(1, 42, 169, 74, 5.1, 0, 1, 0, 0, 0))
			#' seq_des$add_subject_to_experiment_and_assign(c(0, 59, 105, 62, 5.9, 0, 0, 0, 1, 0))
			#' seq_des$add_subject_to_experiment_and_assign(c(1, 32, 186, 66, 5.6, 1, 0, 0, 0, 0))
			#' seq_des$add_subject_to_experiment_and_assign(c(1, 37, 178, 75, 6.5, 0, 0, 0, 0, 1))
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
			#' seq_des$add_subject_to_experiment_and_assign(c(1, 38, 142, 71, 5.3, 0, 0, 0, 1, 0))
			#' 
			#' #returns FALSE since all of the covariate vectors are not yet recorded
			#' seq_des$check_experiment_completed() 
			#' 
			#' seq_des$add_subject_to_experiment_and_assign(c(0, 27, 127, 60, 5.5, 0, 0, 0, 1, 0))
			#' seq_des$add_subject_to_experiment_and_assign(c(1, 42, 169, 74, 5.1, 0, 1, 0, 0, 0))
			#' seq_des$add_subject_to_experiment_and_assign(c(0, 59, 105, 62, 5.9, 0, 0, 0, 1, 0))
			#' seq_des$add_subject_to_experiment_and_assign(c(1, 32, 186, 66, 5.6, 1, 0, 0, 0, 0))
			#' seq_des$add_subject_to_experiment_and_assign(c(1, 37, 178, 75, 6.5, 0, 0, 0, 0, 1))
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
		p = NULL,		
		verbose = NULL,
		
		#when during the experiment the subjects recorded
		y_i_t_i = list(),			
		#design specific parameters
		isKK = FALSE,
		isKK21 = FALSE,
		#Efron
		weighted_coin_prob = NULL,
		#KK parameters
		lambda = NULL,
		t_0 = NULL,
		match_indic = NULL,
		#KK21 parameters
		num_boot = NULL,
		num_responses_to_begin_weight_estimation = NULL,
		
		duplicate = function(){
			self$assert_experiment_completed() #can't duplicate without the experiment being done
			d = SeqDesign$new(private$n, private$p, self$design, self$response_type, self$prob_T, verbose = FALSE)
			#we are assuming the experiment is complete so we have self$t = 0 initialized
			d$X = self$X
			d$y = self$y
			d$dead = self$dead
			d$w = self$w				
			d$.__enclos_env__$private$y_i_t_i = private$y_i_t_i
			d$.__enclos_env__$private$isKK = private$isKK
			d$.__enclos_env__$private$isKK21 = private$isKK21
			d$.__enclos_env__$private$weighted_coin_prob = private$weighted_coin_prob
			d$.__enclos_env__$private$lambda = private$lambda
			d$.__enclos_env__$private$t_0 = private$t_0
			d$.__enclos_env__$private$match_indic = private$match_indic
			d$.__enclos_env__$private$num_boot = private$num_boot
			d
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
				rbinom(1, 1, 1 - private$weighted_coin_prob)
			} else if (n_T * self$prob_T < n_C * (1 - self$prob_T)){
				rbinom(1, 1, private$weighted_coin_prob)
			} else {
				private$assign_wt_CRD()
			}				
		},
		
		assign_wt_Atkinson = function(){
			#if it's too early in the trial or if all the assignments are the same, then randomize
			if (self$t <= private$p + 2 + 1 | length(unique(self$t)) == 1){
				private$assign_wt_CRD()
			} else {
				tryCatch({
					#this matrix is [w | 1 | X]
					X_t_min_1 = cbind(self$w[1 : (self$t - 1)], 1, self$X[1 : (self$t - 1), ])
					XtX = t(X_t_min_1) %*% X_t_min_1						
					M = (self$t - 1) * solve(XtX)
					A = M[1, 2 : (private$p + 2)] %*% c(1, self$X[self$t, ]) 
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
			wt = NULL
			if (length(private$match_indic[private$match_indic == 0]) == 0 | self$t <= private$p | self$t <= private$t_0){
				#we're early, so randomize
				private$match_indic[self$t] = 0
				wt = private$assign_wt_CRD()
			} else {
				# cat("else\n")
				#first calculate the threshold we're operating at
				xs_to_date = self$X[1 : self$t, ]
				S_xs_inv = solve(var(xs_to_date))
				F_crit =  qf(private$lambda, p, self$t - private$p)
				T_cutoff_sq = private$p * (private$n - 1) / (private$n - private$p) * F_crit
				#now iterate over all items in reservoir and take the minimum distance x
				reservoir_indices = which(private$match_indic == 0)
				x_star = self$X[self$t, ]
				sqd_distances = array(NA, length(reservoir_indices))
				for (r in 1 : length(reservoir_indices)){
					sqd_distances[r] = 1 / 2 * 
						t(x_star - self$X[reservoir_indices[r], ]) %*%
						S_xs_inv %*%
						(x_star - self$X[reservoir_indices[r], ])			
				}					
				#find minimum distance index
				min_sqd_dist_index = which(sqd_distances == min(sqd_distances))
				if (length(sqd_distances[min_sqd_dist_index]) > 1 || length(T_cutoff_sq) > 1){
					min_sqd_dist_index = min_sqd_dist_index[1] #if there's a tie, just take the first one
				}
				#if it's smaller than the threshold, we're in business: match it
				if (sqd_distances[min_sqd_dist_index] < T_cutoff_sq){
					match_num = max(private$match_indic, na.rm = TRUE) + 1
					private$match_indic[reservoir_indices[min_sqd_dist_index]] = match_num
					private$match_indic[self$t] = match_num
					#assign opposite
					wt = 1 - self$w[reservoir_indices[min_sqd_dist_index]]
				} else { #otherwise, randomize and add it to the reservoir
					private$match_indic[self$t] = 0
					wt = private$assign_wt_CRD()
				}
			}
			if (is.na(private$match_indic[self$t])){
				stop("no match data recorded")
			}
			wt
		},
		
		compute_weight_KK21_continuous = function(xs_to_date, ys_to_date, deaths_to_date, j){
			ols_mod = lm(ys_to_date ~ xs_to_date[, j])
			#1 - coef(summary(logistic_regr_mod))[2, 4]
			abs(coef(summary(ols_mod))[2, 3])
		},
		
		compute_weight_KK21_incidence = function(xs_to_date, ys_to_date, deaths_to_date, j){
			logistic_regr_mod = suppressWarnings(glm(ys_to_date ~ xs_to_date[, j], family = "binomial"))
			#1 - coef(summary(logistic_regr_mod))[2, 4]
			abs(coef(summary(logistic_regr_mod))[2, 3])	
		},
		
		compute_weight_KK21_count = function(xs_to_date, ys_to_date, deaths_to_date, j){
			negbin_regr_mod = suppressWarnings(MASS::glm.nb(y ~ x, data = data.frame(x = xs_to_date[, j], y = ys_to_date)))
			#1 - coef(summary(negbin_regr_mod))[2, 4]
			abs(coef(summary(negbin_regr_mod))[2, 3])
		},
		
		compute_weight_KK21_proportion = function(xs_to_date, ys_to_date, deaths_to_date, j){
			beta_regr_mod = suppressWarnings(betareg::betareg(y ~ x, data = data.frame(x = xs_to_date[, j], y = ys_to_date)))
			#1 - coef(summary(beta_regr_mod))$mean[2, 4]
			abs(coef(summary(beta_regr_mod))$mean[2, 3])
		},
		
		compute_weight_KK21_survival = function(xs_to_date, ys_to_date, deaths_to_date, j){
			surv_obj = survival::Surv(ys_to_date, deaths_to_date)
			#sometims the weibull is unstable... so try other distributions... this doesn't matter since we are just trying to get weights
			#and we are not relying on the model assumptions
			for (dist in c("weibull", "lognormal", "loglogistic")){
				surv_regr_mod = robust_survreg_with_surv_object(surv_obj, xs_to_date[, j], dist = dist)
				weight = abs(summary(surv_regr_mod)$table[2, 3])
				#1 - summary(weibull_regr_mod)$table[2, 4]
				if (!is.na(weight)){
					return(weight)
				}
			}
			#if that didn't work, default to OLS and log the survival times... again... this doesn't matter since we are just trying to get weights
			#and we are not relying on the model assumptions
			private$compute_weight_KK21_continuous(xs_to_date, log(ys_to_date), deaths_to_date, j)
		},
		
		assign_wt_KK21 = function(){
			wt = NULL
			ys_to_date = self$y[1 : (self$t - 1)]
			
			if (length(private$match_indic[private$match_indic == 0]) == 0 | self$t <= private$p | self$t <= private$t_0){
				#we're early or the reservoir is empty, so randomize
				private$match_indic[self$t] = 0
				wt = private$assign_wt_CRD()
				#cat("    assign_wt_KK21 CRD t", self$t, "\n")
			} else if (sum(!is.na(ys_to_date)) < private$num_responses_to_begin_weight_estimation){ 
				wt = private$assign_wt_KK14()
			}  else {
				#1) need to calculate the weights
				# (a) get old x's and y's (we assume that x's are orthogonal for now)
				xs_to_date = self$X[1 : (self$t - 1), , drop = FALSE]
				deaths_to_date = self$dead[1 : (self$t - 1)]
				# (b) run simple correlations to get Rsq's and use them as the relative weights
				weights = array(NA, private$p)
				for (j in 1 : private$p){
					weights[j] = private[[paste0("compute_weight_KK21_", self$response_type)]](xs_to_date, ys_to_date, deaths_to_date, j) 
				}
				if (any(is.na(weights)) | any(is.infinite(weights)) | any(is.nan(weights))){
					stop("weight values illegal")
				}
				# (c) now we need to normalize the weights
				#cat("    assign_wt_KK21 using weights t", self$t, "raw weights", weights, "\n")
				weights = weights / sum(weights)
				#cat("    assign_wt_KK21 using weights t", self$t, "sorted weights", sort(weights), "\n")
				
				#2) now iterate over all items in reservoir and calculate the weighted sqd distiance vs new guy 
				reservoir_indices = which(private$match_indic == 0)
				x_star = self$X[self$t, ]
				weighted_sqd_distances = array(NA, length(reservoir_indices))
				for (r in 1 : length(reservoir_indices)){
					delta_x = x_star - self$X[reservoir_indices[r], ]
					weighted_sqd_distances[r] = delta_x^2 %*% weights			
				}
				#3) find minimum weighted sqd distiance index
				min_weighted_sqd_dist_index = which(weighted_sqd_distances == min(weighted_sqd_distances))
				
				#generate a cutoff for the weighted minimum distance squared based on bootstrap
				bootstrapped_weighted_sqd_distances = array(NA, private$num_boot)
				for (b in 1 : private$num_boot){
					two_xs  = self$X[sample.int(self$t, 2), ] #self$X[sample_int_ccrank(self$t, 2, rep(1, (self$t))), ] #
					delta_x = two_xs[1, ] - two_xs[2, ]
					bootstrapped_weighted_sqd_distances[b] = delta_x^2 %*% weights
				}
				
				min_weighted_dsqd_cutoff_sq = quantile(bootstrapped_weighted_sqd_distances, private$lambda)
				
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
					wt = 1 - self$w[reservoir_indices[min_weighted_sqd_dist_index]]
				# (b) otherwise, randomize and add it to the reservoir
				} else { 
					private$match_indic[self$t] = 0	
					wt = private$assign_wt_CRD()
				}
			}
			if (is.na(private$match_indic[self$t])){
				stop("no match data recorded")
			}
			wt
		},
		
		compute_weights_KK21stepwise = function(xs_to_date, response_obj, abs_z_compute_fun){
			X_full = data.matrix(xs_to_date[1 : (self$t - 1), ])
			w_to_date = self$w[1 : (self$t - 1)]
			
			#initialize values
			weights = array(NA, private$p)
			j_droppeds = c()
			X_stepwise = matrix(NA, nrow = self$t - 1, ncol = 0)
			
			repeat {
				covs_to_try = setdiff(1 : private$p, j_droppeds)
				#cat(paste("   num covs to try", length(covs_to_try)), "\n")
				#if there's none left, we jet
				if (length(covs_to_try) == 0){
					break
				}
				abs_approx_zs = array(NA, private$p)
				for (j in covs_to_try){
					abs_approx_zs[j] = abs_z_compute_fun(response_obj, cbind(X_full[, j], X_stepwise, w_to_date))
				}
				j_max = which.max(abs_approx_zs)
				weights[j_max] = abs_approx_zs[j_max]
				j_droppeds = c(j_droppeds, j_max)
				X_stepwise = cbind(X_stepwise, X_full[, j_max])
			}
			weights			
		},
		
		compute_weights_KK21stepwise_continuous = function(xs_to_date, ys_to_date, deaths_to_date){
			# we need to standardize the ys
			ys_to_date = scale(ys_to_date)			
			# extract out the relevant time portion and adjust for the effect of betaT for the eventual rand. test
			ys_seen_previously_adj_for_trt = lm.fit(as.matrix(self$w[1 : (self$t - 1)]), ys_to_date)$residuals
			
			private$compute_weights_KK21stepwise(xs_to_date, ys_seen_previously_adj_for_trt, function(response_obj, covariate_data_matrix){
				ols_mod = lm(response_obj ~ covariate_data_matrix)
				abs(coef(summary(ols_mod))[2, 3])
			})
		},
		
		compute_weights_KK21stepwise_incidence = function(xs_to_date, ys_to_date, deaths_to_date){
			private$compute_weights_KK21stepwise(xs_to_date, ys_to_date, function(response_obj, covariate_data_matrix){
				logistic_regr_mod = suppressWarnings(glm(response_obj ~ covariate_data_matrix, family = "binomial"))
				abs(coef(summary(logistic_regr_mod))[2, 3])
			})		
		},
		
		compute_weights_KK21stepwise_count = function(xs_to_date, ys_to_date, deaths_to_date){	
			private$compute_weights_KK21stepwise(xs_to_date, ys_to_date, function(response_obj, covariate_data_matrix){
				negbin_regr_mod = suppressWarnings(MASS::glm.nb(response_obj ~ ., data = cbind(data.frame(response_obj = response_obj), covariate_data_matrix)))
				abs(coef(summary(negbin_regr_mod))[2, 3])
			})	
		},
		
		compute_weights_KK21stepwise_proportion = function(xs_to_date, ys_to_date, deaths_to_date){	
			private$compute_weights_KK21stepwise(xs_to_date, ys_to_date, function(response_obj, covariate_data_matrix){
				beta_regr_mod = suppressWarnings(betareg::betareg(response_obj ~ ., data = cbind(data.frame(response_obj = response_obj), covariate_data_matrix)))
				abs(coef(summary(beta_regr_mod))$mean[2, 3])
			})
		},
		
		compute_weights_KK21stepwise_survival = function(xs_to_date, ys_to_date, deaths_to_date){		
			private$compute_weights_KK21stepwise(xs_to_date, survival::Surv(ys_to_date, deaths_to_date), function(response_obj, covariate_data_matrix){
				surv_regr_mod = robust_survreg_with_surv_object(response_obj, covariate_data_matrix)
				abs(summary(surv_regr_mod)$table[2, 3])
			})	
		},
		
		assign_wt_KK21stepwise = function(){
			wt = NULL
			ys_to_date = self$y[1 : (self$t - 1)]				
			
			if (length(private$match_indic[private$match_indic == 0]) == 0 | self$t <= private$p | self$t <= private$t_0){
				#we're early or the reservoir is empty, so randomize
				private$match_indic[self$t] = 0
				wt = private$assign_wt_CRD()
				#cat("    assign_wt_KK21stepwise CRD t", self$t, "\n")
			} else if (sum(!is.na(ys_to_date)) < private$num_responses_to_begin_weight_estimation){ 
				#we don't have enough data to estimate the weights
				wt = private$assign_wt_KK14()
			} else {				
				#1) need to calculate the weights.. 
				xs_to_date = self$X[1 : self$t,  , drop = FALSE]
				deaths_to_date = self$dead[1 : (self$t - 1)]
				
				# we need to standardize the xs - puts them all on equal footing
				xs_to_date = apply(xs_to_date, 2, scale)
				
				#we need to now run the appropriate (based on response_type) stepwise procedure to get the weights
				weights = private[[paste0("compute_weights_KK21stepwise_", self$response_type)]](xs_to_date, ys_to_date, deaths_to_date) 
				#ensure the weights are normalized
				weights = weights / sum(weights)
				#cat("    assign_wt_KK21stepwise using weights t", self$t, "weights", weights, "\n")
				
				#2) now iterate over all items in reservoir and calculate the weighted sqd distiance vs new guy 
				reservoir_indices = which(private$match_indic == 0)
				x_star = xs_to_date[self$t, ]
				weighted_sqd_distances = array(NA, length(reservoir_indices))
				for (r in 1 : length(reservoir_indices)){
					delta_x = x_star - xs_to_date[reservoir_indices[r], ]
					weighted_sqd_distances[r] = delta_x^2 %*% weights			
				}
				#3) find minimum weighted sqd distiance index
				min_weighted_sqd_dist_index = which(weighted_sqd_distances == min(weighted_sqd_distances))
				
				#generate a cutoff for the weighted minimum distance squared based on bootstrap
				bootstrapped_weighted_sqd_distances = array(NA, private$num_boot)
				for (b in 1 : private$num_boot){
					two_xs  = self$X[sample.int(self$t, 2), ] #self$X[sample_int_ccrank(self$t, 2, rep(1, (self$t))), ] 
					delta_x = two_xs[1, ] - two_xs[2, ]
					bootstrapped_weighted_sqd_distances[b] = delta_x^2 %*% weights
				}
				
				min_weighted_dsqd_cutoff_sq = quantile(bootstrapped_weighted_sqd_distances, private$lambda)
				
				#5) Now, does the minimum make the cut?
				if (length(weighted_sqd_distances[min_weighted_sqd_dist_index]) > 1 || length(min_weighted_dsqd_cutoff_sq) > 1){
					min_weighted_sqd_dist_index = min_weighted_sqd_dist_index[1] #if there's a tie, just take the first one
				}
				#  (a) if it's smaller than the threshold, we're in business: match it
				if (weighted_sqd_distances[min_weighted_sqd_dist_index] < min_weighted_dsqd_cutoff_sq){
					match_num = max(private$match_indic, na.rm = TRUE) + 1
					private$match_indic[reservoir_indices[min_weighted_sqd_dist_index]] = match_num
					private$match_indic[self$t] = match_num
					wt = 1 - self$w[reservoir_indices[min_weighted_sqd_dist_index]]
					# (b) otherwise, randomize and add it to the reservoir
				} else { 
					private$match_indic[self$t] = 0	
					wt = private$assign_wt_CRD()	
				}
			}
			if (is.na(private$match_indic[self$t])){
				stop("no match data recorded")
			}
			wt
		}
	)
)