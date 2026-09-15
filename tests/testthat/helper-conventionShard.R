# Sharding for the whole-database convention sweep.
#
# Checking every model in one test_that() block put the entire suite behind a
# single ~1100s serial call: testthat parallelises across FILES, so one block
# cannot be shared between workers no matter how many are free.  The sweep is
# split across several files instead, and together they still cover every
# model, so the enumerating contract is unchanged.
#
# Shards are interleaved (every nth name) rather than grouped by modeldb
# subdirectory, because the subdirectories are wildly uneven -- specificDrugs
# holds most of the library -- and grouping by them would leave one shard
# carrying nearly all the work, which is the problem being solved.

.conventionShardCount <- 8L

.conventionShardModels <- function(shard) {
  checkmate::assertInt(shard, lower = 1L, upper = .conventionShardCount)
  nms <- get0("modeldb", envir = asNamespace("nlmixr2lib"))$name
  if (is.null(nms)) {
    return(character(0))
  }
  # method = "radix" sorts in the C locale. Plain sort() uses the collation
  # locale, and en_US.UTF-8 and C disagree about these names: 215 of shard 7's
  # 364 models move between the two. The union still covers every model either
  # way, but a failure would land in a different shard on a developer's machine
  # than in CI, which makes it unreproducible exactly when someone needs to
  # reproduce it.
  nms <- sort(nms, method = "radix")
  nms[seq_along(nms) %% .conventionShardCount == (shard - 1L)]
}

.conventionShardCheck <- function(shard) {
  nms <- .conventionShardModels(shard)
  skip_if(length(nms) == 0L, "modeldb is not available")
  # suppressMessages() as well as suppressWarnings(): readModelDb() reports
  # each model's full description through cli_alert_info(), and a shard walks
  # 364 of them. Left alone that is the overwhelming majority of this
  # package's CI log on a *successful* run, which buries anything worth
  # reading. Nothing here is diagnosed from the message stream -- a violation
  # is reported by the assertion, which names the model itself.
  do.call(rbind, lapply(nms, function(nm) {
    suppressMessages(suppressWarnings(nlmixr2lib:::.checkOneModel(nm, verbose = FALSE)))
  }))
}

# Stand-in for .checkOneModel() used to test the no-argument iteration in
# checkModelConventions() without rebuilding the whole library. The visited
# names accumulate in an environment rather than through `<<-`, so the
# assignment target is named at the assignment site.
.conventionStubAcc <- new.env(parent = emptyenv())
.conventionStubAcc$seen <- character(0)

.conventionStubCheckOneModel <- function(model, verbose) {
  .conventionStubAcc$seen <- c(.conventionStubAcc$seen, model)
  data.frame(
    model = model, category = "stub", severity = "info",
    name = "n", message = "m", suggestion = "s",
    stringsAsFactors = FALSE
  )
}
