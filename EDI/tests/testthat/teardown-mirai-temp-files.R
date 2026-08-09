# Windows CI has been observed to leave stray 'RscriptXXXXXXXX' files (no
# extension, hex suffix) in the check working directory after the test
# suite runs -- most likely daemon/IPC artifacts from tests that spin up
# mirai daemons via set_num_cores()/SimulationFramework and tear them down
# with unset_num_cores()/mirai::daemons(0), which doesn't always remove
# every local transport file on that platform. None of this package's own
# code writes files matching this pattern, so sweeping it here at the end
# of the whole suite is a safe, root-cause-agnostic cleanup that keeps
# R CMD check's "detritus in the temp directory" check clean regardless of
# which specific test(s) created them.
stray_mirai_files = list.files(pattern = "^Rscript[0-9a-f]+$", all.files = TRUE, no.. = TRUE)
if (length(stray_mirai_files) > 0) {
	unlink(stray_mirai_files)
}
