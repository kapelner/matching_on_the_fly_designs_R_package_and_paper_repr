#/c/Program\ Files/R/R-devel/bin/R.exe CMD INSTALL -l ~/AppData/Local/R/win-library/4.5/ SeqExpMatch/
cat("DEBUG: Line 2\n"); flush.console()
rm(list = ls())
cat("DEBUG: Line 3 - After rm\n"); flush.console()
pacman::p_load(doParallel, PTE, datasets, qgam, mlbench, AppliedPredictiveModeling, dplyr, ggplot2, gridExtra, profvis, data.table, profvis)
cat("DEBUG: Line 4 - After pacman\n"); flush.console()
file_arg = "--file="
cat("DEBUG: Line 5\n"); flush.console()
script_path = sub(file_arg, "", commandArgs(trailingOnly = FALSE)[grep(file_arg, commandArgs(trailingOnly = FALSE))])
cat("DEBUG: Line 6 - script_path:", script_path, "\n"); flush.console()
test_dir = dirname(normalizePath(script_path))
cat("DEBUG: Line 7 - test_dir:", test_dir, "\n"); flush.console()
lib_path = normalizePath(file.path(test_dir, "..", "Rlib"), mustWork = FALSE)
cat("DEBUG: Line 8 - lib_path:", lib_path, "\n"); flush.console()
if (dir.exists(lib_path)){
  cat("DEBUG: Line 9 - lib_path exists, adding to .libPaths\n"); flush.console()
  .libPaths(c(lib_path, .libPaths()))
}
cat("DEBUG: Line 11 - Before library(SeqExpMatch)\n"); flush.console()
library(SeqExpMatch)
cat("DEBUG: Line 12 - After library(SeqExpMatch)\n"); flush.console()
debug_mode = nzchar(Sys.getenv("SEQEXPMATCH_DEBUG"))
cat("DEBUG: Line 13 - debug_mode:", debug_mode, "\n"); flush.console()
if (debug_mode){
  options(error = function(){
    traceback(2)
    q(status = 1)
  })
} else {
  options(error = NULL)
}
cat("DEBUG: Line 20 - After options\n"); flush.console()
# options(warn=2)
set.seed(1)
cat("DEBUG: Line 22 - After set.seed\n"); flush.console()

max_n_dataset = 500
cat("DEBUG: Line 24 - max_n_dataset set\n"); flush.console()
cat("DEBUG: Line 25 - About to source _dataset_load.R from:", test_dir, "\n"); flush.console()
source(file.path(test_dir, "_dataset_load.R"))
cat("DEBUG: Line 26 - After source\n"); flush.console()

cat("DEBUG: Reached end of script initialization\n"); flush.console()
