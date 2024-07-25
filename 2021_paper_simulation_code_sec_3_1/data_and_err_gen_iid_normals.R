#build mvnp covariates
x_s = matrix(rnorm(n * p, mu_x, sigma_x), ncol = p)

#build errors independent of x (and thus a priori of treatment allocation)
errors = rnorm(n, 0, sigma_e)
