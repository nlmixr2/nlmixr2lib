# Diagnostic only (windows-crash-probe.yaml). On windows-latest, R dies with a
# segfault in test-checkModelConventions.R at the first rxModelVars() call of
# ".doseTimingCalls skips an empty first argument and a namespaced call", which
# hands rxode2's parser model text it must reject. The test before it parses
# every dose-history function and one rejected model. This replays that mix
# against rxode2 alone, in a loop, printing each parse before it starts.
#
# Argument: the number of rounds (default 40). Each round parses distinct valid
# models (a different constant per round, so rxModelVars() cannot answer from
# its cache) and then the four rejected spellings.
reps <- as.integer(commandArgs(TRUE)[1])
if (is.na(reps)) reps <- 40L
library(rxode2)
for (p in c("rxode2", "rxode2lincmt", "rxode2ll", "dparser", "lotri")) {
  d <- tryCatch(utils::packageDescription(p), warning = function(w) NULL)
  if (!is.null(d)) {
    message(p, " ", d$Version, " ", d$RemoteSha %||% "", " ", d$Repository %||% "", " built ", d$Built)
  }
}
good <- c("tad", "tafd", "tlast", "tfirst", "podo", "dosenum", "tlast0", "tfirst0")
bad <- c("x <- dosenum(depot)", "x <- rxode2::tad(depot)", "x <- tad(\"depot\")", "x <- tad(, depot)")
parseOne <- function(code) {
  m <- paste0("d/dt(depot) <- -depot; ", code)
  message("parse: ", code)
  res <- tryCatch(
    {
      invisible(utils::capture.output(suppressMessages(rxode2::rxModelVars(m))))
      "ok"
    },
    error = function(e) paste("error:", conditionMessage(e))
  )
  message("    -> ", substr(gsub("\n", " | ", res), 1, 200))
}
for (r in seq_len(reps)) {
  message("== round ", r)
  for (fn in good) parseOne(paste0("x <- ", fn, if (fn == "dosenum") "()" else "(depot)", " + ", r))
  for (code in bad) parseOne(code)
}
message("done: ", reps, " rounds without a crash")
