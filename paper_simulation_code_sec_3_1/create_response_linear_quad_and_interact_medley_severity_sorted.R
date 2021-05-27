z = betas[1] * x_s[, 1] +
		betas[2] * x_s[, 2] + 
		betas[3] * x_s[, 1]^2 +
		betas[4] * x_s[, 2]^2 +
		betas[5] * x_s[, 1] * x_s[, 2]
z_plus_errors = z + errors
z_plus_errors = sort(z_plus_errors)
y = beta_T * indic_T + z_plus_errors