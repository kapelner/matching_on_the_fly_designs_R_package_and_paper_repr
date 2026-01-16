library(testthat)
library(SeqExpMatch)

if (identical(Sys.getenv("NOT_CRAN"), "true")) {
	test_check("SeqExpMatch")
}
