indic_T = sample(
	c(
		rep(0, n * (1 - prob_trt)), 
		rep(1, n * prob_trt)
	)
)
sys.source("common_crd_and_bcrd.R", envir = environment(), toplevel.env = environment())
