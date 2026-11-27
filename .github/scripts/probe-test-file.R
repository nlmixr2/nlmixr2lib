# Diagnostic only (windows-crash-probe.yaml). Runs one test file, or the whole
# directory with "all", against the installed package, serially, with
# testthat's location reporter writing to stderr so the last "Start test" line
# survives a crash of the R process.
part <- commandArgs(TRUE)[1]
library(nlmixr2lib)
for (p in c("rxode2", "rxode2lincmt", "nlmixr2est", "lotri", "dparser")) {
  d <- utils::packageDescription(p)
  message(p, " ", d$Version, " ", d$RemoteSha %||% "", " built ", d$Built)
}
message("DLLs loaded at start: ", length(getLoadedDLLs()), "; R_MAX_NUM_DLLS limit applies")
rep <- testthat::LocationReporter$new(file = stderr())
res <- if (identical(part, "all")) {
  testthat::test_dir(
    "tests/testthat",
    reporter = rep,
    package = "nlmixr2lib",
    load_package = "installed",
    stop_on_failure = FALSE
  )
} else {
  testthat::test_file(
    file.path("tests/testthat", paste0("test-", part, ".R")),
    reporter = rep,
    package = "nlmixr2lib",
    load_package = "installed",
    stop_on_failure = FALSE
  )
}
message("DLLs loaded at end: ", length(getLoadedDLLs()))
df <- as.data.frame(res)
bad <- df[df$failed > 0 | df$error, c("file", "test", "failed", "error"), drop = FALSE]
print(bad)
if (nrow(bad)) quit(status = 1)
