#' Extended Robins Blocked Incidence Inference
#'
#' Unadjusted blocked-design incidence inference using the simple mean-difference
#' point estimate with a block-stratified standard error.
#'
#' @export
InferenceIncidExtendedRobins = R6::R6Class("InferenceIncidExtendedRobins",
	lock_objects = FALSE,
	inherit = InferenceIncidAzriel,

	private = list(
		get_standard_error = function(){
		  if (!is.null(private$cached_values$robins_s_beta_hat_T)) {
		    return(private$cached_values$robins_s_beta_hat_T)
		  }
           m = private$des_obj$get_block_ids()
           B = length(unique(m))
           n_B = sum(m == 1)
           n_B_over_two = n_B / 2
           var_robbins = 0
           for (b in 1 : B){
                   y_b = private$des_obj_priv_int$y[m == b]
                   w_b = private$des_obj_priv_int$w[m == b]
                   p_hat_T_b = sum(y_b[w_b == 1]) / n_B_over_two
                   p_hat_C_b = sum(y_b[w_b == 0]) / n_B_over_two
                   m_1_b = max(p_hat_T_b, p_hat_C_b)
                   m_0_b = min(p_hat_T_b, p_hat_C_b)
                   var_robbins = var_robbins + 
                           m_1_b * (1 - m_1_b) / n_B_over_two +
                           m_0_b * (1 - m_0_b) / n_B_over_two +
                           ((2 * m_0_b - m_1_b) * (1 - m_1_b) - m_0_b * (1 - m_0_b))  / n_B
           }

			p_hat_T = mean(private$des_obj_priv_int$y[private$des_obj_priv_int$w == 1])
			p_hat_C = mean(private$des_obj_priv_int$y[private$des_obj_priv_int$w == 0])
			var_robbins_ext = 1 / private$des_obj_priv_int$n * (
				p_hat_T * (1 - p_hat_T) + p_hat_C * (1 - p_hat_C)
			)
			private$cached_values$robins_s_beta_hat_T = sqrt(var_robbins + var_robbins_ext)
			private$cached_values$robins_s_beta_hat_T
		}
	)
)
