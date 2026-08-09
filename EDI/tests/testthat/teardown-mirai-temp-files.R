# Windows CI has been observed to leave stray 'RscriptXXXXXXXX' files (no
# extension, hex suffix) in R's session temp directory after the test
# suite runs -- most likely daemon/IPC artifacts from tests that spin up
# mirai daemons via set_num_cores()/SimulationFramework and tear them down
# with unset_num_cores()/mirai::daemons(0), which doesn't always remove
# every local transport file on that platform. None of this package's own
# code writes files matching this pattern, so sweeping it here at the end
# of the whole suite is a safe, root-cause-agnostic cleanup that keeps
# R CMD check's "detritus in the temp directory" check clean regardless of
# which specific test(s) created them.
#
# Must search tempdir() (and its parent), not the working directory: R CMD
# check's "detritus in the temp directory" check scans R's own session temp
# dir (confirmed via an actual CI run on 2026-08-09 -- this teardown ran
# but the NOTE still fired, since list.files() with no path defaults to
# getwd(), which during testthat execution is tests/testthat, not tempdir()).
# mirai launches daemons as literal background Rscript processes rather than
# via a filesystem-based IPC transport on Windows (named pipes, not files),
# so these are very likely artifacts of that launch mechanism rather than
# the transport itself -- exactly where such a file lands (R's per-session
# tempdir() vs. the shared TMPDIR root one level up, dirname(tempdir())) is
# not fully certain, so both are swept here to be safe.
stray_mirai_files = list.files(
	c(tempdir(), dirname(tempdir())),
	pattern = "^Rscript[0-9a-f]+$", all.files = TRUE, no.. = TRUE, full.names = TRUE
)
if (length(stray_mirai_files) > 0) {
	unlink(stray_mirai_files)
}
