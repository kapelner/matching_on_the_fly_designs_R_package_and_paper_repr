library(testthat)
library(EDI)

if (identical(Sys.getenv("NOT_CRAN"), "true")) {
	test_check("EDI")
}
