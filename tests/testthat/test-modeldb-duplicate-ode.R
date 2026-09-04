# A model file that declares `d/dt(<state>)` more than once builds FEWER ODE
# states than the paper published, silently. rxode2 keeps one of the competing
# equations and drops the rest; nothing errors, nothing warns, and the model
# still solves -- it just solves the wrong system. Two shipped models carried
# this defect undetected:
#
#   * Ngo_2020_HL2351.R declared `d/dt(peripheral1)` twice -- once for the FcRn
#     absorption/distribution space and once for the true peripheral
#     compartment -- so a 4-compartment published model built 3 states and the
#     distribution space was merged into the peripheral. Fixed 2026-09-04 by
#     repointing the distribution space at the `abs_site` canonical.
#   * LeTilly_2021_trastuzumab.R declares `d/dt(lat1)` twice; see the
#     quarantine block below.
#
# Prose cannot catch this class of bug, so this file is the mechanical gate.

# Extract the `d/dt(...)` left-hand-side state names from a model file, in
# declaration order, by walking the parsed `model({})` block. Parsing rather
# than grepping means a state named inside a comment or a string cannot
# produce a false positive.
odeStateNames <- function(path) {
  env <- new.env(parent = globalenv())
  fname <- sub("[.]R$", "", basename(path))
  sys.source(path, envir = env)
  if (!exists(fname, envir = env, inherits = FALSE)) {
    stop("`", fname, "` not defined by ", path)
  }
  b <- body(get(fname, envir = env))
  modelBlock <- NULL
  for (i in seq_along(b)) {
    node <- b[[i]]
    if (is.call(node) && identical(node[[1]], as.name("model"))) {
      modelBlock <- node[[2]]
    }
  }
  if (is.null(modelBlock)) return(character(0))
  states <- character(0)
  for (i in seq_along(modelBlock)) {
    st <- modelBlock[[i]]
    isAssign <- is.call(st) && length(st) >= 3 &&
      (identical(st[[1]], as.name("<-")) || identical(st[[1]], as.name("=")))
    if (!isAssign) next
    target <- st[[2]]
    # `d/dt(x)` parses as `/`(d, dt(x))
    isDdt <- is.call(target) &&
      identical(target[[1]], as.name("/")) &&
      identical(target[[2]], as.name("d")) &&
      is.call(target[[3]]) &&
      identical(target[[3]][[1]], as.name("dt"))
    if (isDdt) states <- c(states, as.character(target[[3]][[2]]))
  }
  states
}

duplicatedOdeStates <- function(path) {
  states <- odeStateNames(path)
  counts <- table(states)
  names(counts)[counts > 1]
}

modelFiles <- function() {
  root <- system.file("modeldb", package = "nlmixr2lib")
  if (!nzchar(root)) skip("nlmixr2lib modeldb directory not found")
  list.files(root, pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
}

# Models known to carry a duplicate d/dt that CANNOT be repaired from sources
# on hand. Keep this list as close to empty as possible, and never add an entry
# without recording why the fix is blocked.
#
# LeTilly_2021_trastuzumab: `d/dt(lat1)` is declared at both the third and
# fourth position of the latent-HER2 transit chain, and `lat1(0)` is likewise
# set twice, so the intended final chain member (presumably `lat3`) is missing.
# Repairing it requires knowing the chain length Le Tilly 2021 actually fitted
# and which member drives the negative feedback (equation 8) and the
# binding-driven elimination -- none of which is recoverable from the model file
# alone. The source is Clin Pharmacol Ther 2021;110(1):210-219
# (doi:10.1002/cpt.2188), which is Wiley subscription-only with no PMC deposit
# and no open-access licence, so it has been queued for acquisition rather than
# guessed at. Remove this entry -- and this comment -- in the same commit that
# fixes the model.
knownDuplicateOdeModels <- c("LeTilly_2021_trastuzumab")

test_that("no model file declares the same d/dt state twice", {
  files <- modelFiles()
  expect_gt(length(files), 0)

  offenders <- list()
  for (path in files) {
    dups <- tryCatch(duplicatedOdeStates(path), error = function(e) character(0))
    if (length(dups)) offenders[[sub("[.]R$", "", basename(path))]] <- dups
  }

  unexpected <- setdiff(names(offenders), knownDuplicateOdeModels)
  expect_equal(
    unexpected, character(0),
    info = paste0(
      "Model(s) declaring a duplicate d/dt state: ",
      paste(vapply(unexpected, function(nm) {
        paste0(nm, " (", paste(offenders[[nm]], collapse = ", "), ")")
      }, character(1)), collapse = "; "),
      ". A duplicated d/dt silently builds fewer ODE states than declared -- ",
      "rename the state (see inst/references/compartment-names.md) rather than ",
      "adding it to knownDuplicateOdeModels."
    )
  )
})

test_that("the duplicate-d/dt quarantine list contains no already-fixed model", {
  # Without this, a quarantine entry outlives the bug it documents and the
  # gate above silently stops covering that model. Fixing a quarantined model
  # MUST fail this test until its entry is removed.
  files <- modelFiles()
  byName <- stats::setNames(files, sub("[.]R$", "", basename(files)))

  for (nm in knownDuplicateOdeModels) {
    expect_true(
      nm %in% names(byName),
      info = paste0(nm, " is quarantined but no such model file exists; ",
                    "remove it from knownDuplicateOdeModels.")
    )
    if (!nm %in% names(byName)) next
    dups <- duplicatedOdeStates(byName[[nm]])
    expect_gt(
      length(dups), 0,
      # nolint next: line_length_linter.
      label = paste0(nm, " no longer has a duplicate d/dt, so it must be removed from knownDuplicateOdeModels")
    )
  }
})
