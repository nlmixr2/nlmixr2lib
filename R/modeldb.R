#' Model library for nlmixr2
#'
#' This is a data frame of the available models in nlmixr2lib, it is generated
#' with the package.  Custom modeldb may be used.
#'
#' @eval buildModelDb()
"modeldb"

utils::globalVariables(c("modeldb"))

## nocov start

# Accumulator for error-severity convention issues found while loading models,
# and for model files that could not be loaded at all. An environment rather
# than `<<-`: the assignment target is named at the assignment site, and it
# survives being called from a nested loader.
.conventionErrorEnv <- new.env(parent = emptyenv())

.conventionErrorsReset <- function() {
  .conventionErrorEnv$byModel <- list()
  .conventionErrorEnv$failures <- list()
  invisible(NULL)
}

.conventionErrorsAdd <- function(modelName, issues) {
  if (is.null(.conventionErrorEnv$byModel)) {
    .conventionErrorsReset()
  }
  .conventionErrorEnv$byModel[[modelName]] <- issues
  invisible(NULL)
}

# Record a model file whose loading raised an R error (a parse error, a
# filename/function-name mismatch, a model nlmixr2 cannot read).
.modeldbFailureAdd <- function(file, message) {
  if (is.null(.conventionErrorEnv$failures)) {
    .conventionErrorsReset()
  }
  .conventionErrorEnv$failures[[file]] <- message
  invisible(NULL)
}

# The most load failures listed one by one in the stop message; a broken
# dependency fails every model, and thousands of identical lines bury the cause.
.modeldbFailureListMax <- 20L

.conventionErrorsStopIfAny <- function() {
  byModel <- .conventionErrorEnv$byModel
  failures <- .conventionErrorEnv$failures
  if (length(byModel) == 0 && length(failures) == 0) {
    return(invisible(NULL))
  }
  parts <- character()
  if (length(failures) > 0) {
    shown <- utils::head(names(failures), .modeldbFailureListMax)
    detail <- paste0("  ", shown, ": ", unlist(failures[shown], use.names = FALSE))
    if (length(failures) > length(shown)) {
      detail <- c(detail, sprintf("  ... and %d more", length(failures) - length(shown)))
    }
    parts <- c(
      parts,
      sprintf(
        "buildModelDb(): %d model file(s) could not be loaded.\n%s",
        length(failures),
        paste(detail, collapse = "\n")
      )
    )
  }
  if (length(byModel) > 0) {
    n <- sum(vapply(byModel, nrow, integer(1)))
    detail <- vapply(names(byModel), .conventionErrorsDetail, character(1), byModel = byModel)
    parts <- c(
      parts,
      sprintf(
        paste0(
          "buildModelDb(): %d convention error(s) in %d model(s).\n%s\n",
          "Run checkModelConventions(\"<model>\") for the full report and ",
          "suggested fix."
        ),
        n,
        length(byModel),
        paste(detail, collapse = "\n")
      )
    )
  }
  stop(paste(parts, collapse = "\n"), call. = FALSE)
}

.conventionErrorsDetail <- function(nm, byModel) {
  iss <- byModel[[nm]]
  paste0("  ", nm, ": ", paste(sprintf("%s (%s)", iss$name, iss$category), collapse = ", "))
}

#' Build the default model database used in nlmixr2lib
#'
#' The function is only used during the build process, and it is not needed for
#' general use.  See the documentation of "modeldb" for its output.
#'
#' Results are cached per model file in `data-raw/modeldb-cache.rds`, keyed on
#' the file's md5 hash. The cache is written whenever the function exits,
#' including when it stops on a convention error, a model that cannot be
#' loaded, or an interrupt: every model that loaded cleanly is kept, and the
#' models that failed are left out so the next build re-checks only them. A
#' parallel build records its models a chunk at a time, so an interrupt loses
#' the chunks still being parsed.
#'
#' @param workers Number of R processes that parse the models missing from the
#'   cache. `1` parses them in this R session. With more than one, the misses
#'   are split into chunks of at most 200 models (`.modeldbChunkMax`), each
#'   parsed by a fresh `callr` process, so a process never accumulates more
#'   than one chunk's worth of rxode2 state. 20 or fewer misses always run in
#'   this session. The rows, their order, the cache and the errors reported are
#'   the same for every worker count; the log comes out in the order the chunks
#'   finish. Workers run the code this session runs: under a pkgload
#'   development load (`devtools::document()`, `load_all()`) each worker re-runs
#'   `pkgload::load_all()` on the source tree as it starts, so a file edited
#'   during a build reaches the workers started after the edit; under an
#'   installed package they load it from the same library. When neither holds
#'   (the code was sourced, as `roxygenise(load_code = roxygen2::load_source)`
#'   does) or `callr` is not installed, the misses are parsed in this session
#'   and a message says why. The default is the option
#'   `nlmixr2lib.buildModelDb.workers` when set, otherwise
#'   `min(8, max(1, rxode2::getRxThreads() %/% 2))`: a quarter of the cores on
#'   rxode2's own default (half the cores), and smaller under `R CMD check`,
#'   `OMP_THREAD_LIMIT` or `rxode2::setRxThreads()`, leaving room for the other
#'   work on a shared machine. An override must be a single positive whole
#'   number.
#' @return The text needed for documenting the "modeldb" object (see the modeldb
#'   roxygen eval statement)
#' @noRd
buildModelDb <- function(workers = getOption("nlmixr2lib.buildModelDb.workers", .modeldbDefaultWorkers())) {
  workers <- .modeldbCheckWorkers(workers)
  # Find the package root directory
  packageDirectory <- normalizePath(".", winslash = "/")
  origPackageDirectory <- packageDirectory
  while (!file.exists(file.path(packageDirectory, "DESCRIPTION"))) {
    oldPackageDirectory <- packageDirectory
    packageDirectory <- normalizePath(file.path(oldPackageDirectory, ".."), winslash = "/")
    if (oldPackageDirectory == packageDirectory) {
      stop("Could not find root package directory from ", origPackageDirectory)
    }
  }
  message("Building the modeldb from ", packageDirectory)
  .conventionErrorsReset()
  cachePath <- file.path(packageDirectory, "data-raw/modeldb-cache.rds")
  cache <- .modeldbCacheRead(cachePath)
  sig <- .modeldbGlobalSig(packageDirectory)
  if (is.null(cache) || !identical(cache$globalSig, sig)) {
    if (!is.null(cache)) {
      message("modeldb cache invalidated (global signature changed)")
    }
    cache <- list(globalSig = sig, entries = list())
  }
  acc <- new.env(parent = emptyenv())
  acc$entries <- cache$entries
  # Written on every exit, so a build that stops part-way keeps what it did.
  on.exit(.modeldbCacheWrite(cachePath, list(globalSig = sig, entries = acc$entries)), add = TRUE)
  result <- .addDirToModelDbCached(
    dir = file.path(packageDirectory, "inst/modeldb"),
    entries = cache$entries,
    workers = workers,
    acc = acc,
    reloadRegisters = TRUE
  )
  modeldb <- result$modeldb
  # The canonical-name registers are checked against the tree being built, not
  # the installed copy, so a duplicate cannot ride in on a stale install.
  refIssues <- .referenceDuplicateIssues(
    .collateGlob(file.path(packageDirectory, "inst/references/*.md"))
  )
  if (nrow(refIssues)) {
    .conventionErrorsAdd("inst/references", refIssues)
  }
  # Abort before writing artifacts: a database with convention errors in it
  # must never be saved. The cache written on exit holds only the models that
  # loaded without error, so the next build re-checks the ones that did not.
  .conventionErrorsStopIfAny()
  # Drop the base package directory name so that will be installation-agnostic
  modeldb$filename <-
    gsub(
      x = modeldb$filename,
      pattern = "^.*modeldb/",
      replacement = ""
    )
  savefile <- file.path(packageDirectory, "data/modeldb.rda")
  message("Saving the modeldb to ", savefile)
  save(modeldb, file = savefile, compress = "bzip2", version = 2, ascii = FALSE)
  saveRDS(modeldb, file = file.path(packageDirectory, "inst/modeldb.rds"), compress = "bzip2", version = 2)
  .writePkgdownNavbar(modeldb, packageDirectory)
  message("Done saving the modeldb to ", savefile)

  colDesc <-
    list(
      name = "Model name that can be used to extract the model from the model library",
      description = "Model description in free from text; in model itself",
      parameters = paste(
        "A comma separated string listing either the parameter in the model",
        "defined by population/individual effects or a population effect parameter"
      ),
      DV = "The definition of the dependent variable(s)",
      linCmt = "Logical flag indicating if solved models are used (TRUE) or not (FALSE)",
      algebraic = paste(
        "Logical flag indicating if the model is purely algebraic:",
        "TRUE no linCmt() and no ODEs; FALSE otherwise"
      ),
      dosing = "A comma separated string of identified dosing compartments",
      depends = "A comma separated string of objects the model depends on",
      vignette = "Basename of the vignette associated with this model (without path or extension); NA if none",
      label = paste(
        "Human-readable navbar label derived from the filename: 'Drug (Author Year)'",
        "for the canonical '<Author>_<Year>_<drug>' form, 'DDMoRe: <drug>' for the",
        "parameterless NA_NA_<drug> ddmore entries, otherwise the basename with",
        "underscores replaced by spaces"
      ),
      category = paste(
        "Coarse-grained model bucket derived from the filename path:",
        "'specificDrugs', 'ddmore', or 'other' (anything not under those two",
        "directories)"
      ),
      filename = "Filename of the model.  By default these are installed in the model library and read on demand"
    )
  # The names must exactly match
  stopifnot(all(names(modeldb) %in% names(colDesc)))
  stopifnot(all(names(colDesc) %in% names(modeldb)))
  formatText <- sprintf("@format A data frame with %g rows and %g columns", nrow(modeldb), ncol(modeldb))
  describeText <-
    sprintf(
      "\\describe{\n%s}\n",
      paste(
        sprintf("  \\item{%s}{%s}\n", names(colDesc), unlist(colDesc)),
        collapse = ""
      )
    )
  paste(formatText, describeText, sep = "\n")
}

#' Add a directory to the modeldb
#'
#' @param dir Directory name containing model files
#' @param modeldb The starting modeldb data.frame
#' @return The updated modeldb data.frame
#' @export
addDirToModelDb <- function(dir, modeldb = data.frame()) {
  filesToLoad <-
    .collateListFiles(
      path = dir,
      pattern = "\\.R$",
      ignore.case = TRUE,
      recursive = TRUE
    )
  for (currentFile in filesToLoad) {
    message("parse currentFile")
    modeldb <- addFileToModelDb(dir = dir, file = currentFile, modeldb = modeldb)
  }
  modeldb
}

#' @describeIn addDirToModelDb Add a file to the modeldb
#' @param file The file name (without the directory name)
#' @return the model database
#' @export
addFileToModelDb <- function(dir, file, modeldb) {
  res <- .modeldbFileRow(dir = dir, file = file)
  if (!is.null(res$errors)) {
    .conventionErrorsAdd(res$row$name, res$errors)
  }
  modeldb <- rbind(modeldb, res$row)
  if (any(duplicated(modeldb$name))) {
    stop("Duplicated model name: ", modeldb$name[duplicated(modeldb$name)]) # nocov
  }
  modeldb
}

# Load one model file and return its modeldb row, with its error-severity
# convention issues (NULL when there are none). Touches no shared state, so it
# can run in a worker process; addFileToModelDb() and .modeldbProcessFiles()
# decide what to do with the errors.
.modeldbFileRow <- function(dir, file) {
  fileName <- file.path(dir, file)

  # Extract the model from the file
  parsedFile <- parse(file = fileName)
  stopifnot(identical(parsedFile[[1]][[1]], as.name("<-")))

  modelName <- as.character(parsedFile[[1]][[2]])
  packageStartupMessage("Loading ", modelName, " from ", fileName)
  if (modelName != tools::file_path_sans_ext(basename(file))) {
    stop("Loading model failed due to filename/modelName mismatch: ", fileName, call. = FALSE) # nocov
  }

  # Parse the model to get the fixed effects and DV parameters
  mod <- nlmixr2est::nlmixr(eval(parsedFile))
  # rxode2 names the ui after the deparsed argument ("eval(parsedFile)"). The
  # checker identifies the model by this name, and a covariate registered with
  # `Scope: specific` is allowed only in the models its entry lists.
  mod$modelName <- modelName

  # Convention check. Error-severity issues are accumulated and abort
  # `buildModelDb()` once every model has been visited, so one run reports
  # every offending model instead of one per rebuild.
  issues <- tryCatch(
    suppressWarnings(checkModelConventions(mod, verbose = FALSE)),
    error = function(e) NULL
  )
  errors <- NULL
  if (!is.null(issues) && nrow(issues) > 0) {
    n_err <- sum(issues$severity == "error")
    n_warn <- sum(issues$severity == "warning")
    if (n_err > 0) {
      errors <- issues[issues$severity == "error", , drop = FALSE]
    }
    if (n_err + n_warn > 0) {
      message(sprintf(
        "  %s: %d convention error(s), %d warning(s) - run checkModelConventions(\"%s\") for details",
        modelName,
        n_err,
        n_warn,
        modelName
      ))
    }
  }

  description <- mod$meta$description
  if (is.null(description)) {
    message("No description for model in ", fileName)
    description <- NA_character_
  }

  vignette <- mod$meta$vignette
  if (is.null(vignette)) {
    vignette <- NA_character_
  }

  # Finding dosing
  dosing <- NULL
  dosing_meta <- mod$meta$dosing
  if (!is.null(dosing_meta)) {
    dosing <- paste(dosing_meta, collapse = ",")
  } else {
    if ("depot" %in% mod$props$cmt) {
      dosing <- c(dosing, "depot")
    }
    if ("central" %in% mod$props$cmt) {
      dosing <- c(dosing, "central")
    }
    if (!is.null(dosing)) {
      dosing <- paste(dosing, collapse = ",")
    } else {
      dosing <- NA_character_
    }
  }

  # Finding depends
  depends <- NULL
  depends_meta <- mod$meta$depends
  if (!is.null(depends_meta)) {
    depends <- paste(depends_meta, collapse = ",")
  }
  if (is.null(depends)) {
    depends <- NA_character_
  }

  # Extract the parameter names
  modParam <- mod$iniDf
  # Fixed effects
  modParamFixed <- modParam$name[is.na(modParam$neta1) & is.na(modParam$err)]

  # swap modeled parameter names for the mu-ref parameter names, where
  # applicable
  .ref <- .getVarLhs(mod)
  for (nm in names(.ref)) {
    modParamFixed[modParamFixed %in% nm] <- .ref[nm]
  }

  # Error model
  paramErr <- mod$predDf$cond
  if ("rxLinCmt" %in% paramErr) {
    paramErr[paramErr %in% "rxLinCmt"] <- "linCmt()"
  }
  if (is.null(paramErr)) {
    paramErr <- ""
  }

  if (!mod$props$linCmt && (length(mod$props$cmt) == 0)) {
    algebraic <- TRUE
  } else {
    algebraic <- FALSE
  }

  ret <-
    data.frame(
      name = modelName,
      description = description,
      parameters = paste(modParamFixed, collapse = ","),
      DV = paste(paramErr, collapse = ","),
      linCmt = mod$props$linCmt,
      algebraic = algebraic,
      dosing = dosing,
      depends = depends,
      vignette = vignette,
      label = .parseModelLabel(fileName),
      category = .parseModelCategory(fileName),
      filename = fileName
    )
  list(row = ret, errors = errors)
}
## nocov end

# Walk a modeldb directory and return rows + updated cache entries. Reuses a
# cached row when the file's md5 matches the cached hash; otherwise re-parses
# via .modeldbFileRow(), in this session or in `workers` worker processes.
#
# `acc$entries` is the cache as it stands, so the caller can save it however
# the walk ends. It is updated after each model in this session and after each
# finished chunk in a parallel walk (an interrupt loses the chunks still
# running). It starts from the cached entries of the files still on disk; a
# model that loads cleanly replaces its entry and one that errors (convention
# or R error) loses it. Errors are recorded with .conventionErrorsAdd() /
# .modeldbFailureAdd() in file order, the same for every worker count, and
# raised by the caller.
#
# `chunkSize` overrides the chunking computed by .modeldbChunkPlan(); tests
# use it to put a handful of fixture models through the worker path.
#
# `reloadRegisters = TRUE` re-parses the canonical-name registers before the
# first miss (once per call, and not at all when every file is cached), so a
# build checks against the registers on disk rather than a session's memo.
.addDirToModelDbCached <- function(dir, entries, workers = 1L, acc = NULL, chunkSize = NULL, reloadRegisters = FALSE) {
  filesToLoad <-
    .collateListFiles(
      path = dir,
      pattern = "\\.R$",
      ignore.case = TRUE,
      recursive = TRUE
    )
  if (is.null(acc)) {
    acc <- new.env(parent = emptyenv())
  }
  acc$entries <- entries[names(entries) %in% filesToLoad]
  hashes <- unname(tools::md5sum(file.path(dir, filesToLoad)))
  rows <- vector("list", length(filesToLoad))
  miss <- logical(length(filesToLoad))
  for (i in seq_along(filesToLoad)) {
    cached <- entries[[filesToLoad[i]]]
    if (!is.null(cached) && identical(cached$hash, hashes[i])) {
      message("modeldb cache hit: ", filesToLoad[i])
      row <- cached$row
      row$filename <- file.path(dir, filesToLoad[i])
      rows[[i]] <- row
    } else {
      miss[i] <- TRUE
    }
  }
  if (any(miss)) {
    # Workers are handed this session's parse, so a build parses each register
    # once however many processes it uses.
    if (reloadRegisters) {
      .modeldbReloadRegisters()
    }
    missIdx <- which(miss)
    results <- .modeldbRunMisses(
      dir = dir,
      files = filesToLoad[missIdx],
      hashes = hashes[missIdx],
      workers = workers,
      acc = acc,
      chunkSize = chunkSize
    )
    for (j in seq_along(missIdx)) {
      res <- results[[j]]
      if (!is.null(res$failure)) {
        .modeldbFailureAdd(file.path(dir, res$file), res$failure)
        next
      }
      if (!is.null(res$errors)) {
        .conventionErrorsAdd(res$row$name, res$errors)
      }
      rows[[missIdx[j]]] <- res$row
    }
  }
  rows <- rows[lengths(rows) > 0L]
  modeldb <- if (length(rows) > 0L) do.call(rbind, rows) else data.frame()
  if (anyDuplicated(modeldb$name)) {
    stop(
      "Duplicated model name: ",
      paste(unique(modeldb$name[duplicated(modeldb$name)]), collapse = ", "),
      call. = FALSE
    )
  }
  acc$entries <- acc$entries[filesToLoad[filesToLoad %in% names(acc$entries)]]
  list(modeldb = modeldb, entries = acc$entries)
}

# Result for one parsed file: its cache key and hash, and either the row (with
# any error-severity convention issues) or `failure`, the message of the R
# error that stopped it loading.
.modeldbFileResult <- function(file, hash, row = NULL, errors = NULL, failure = NULL) {
  list(file = file, hash = hash, row = row, errors = errors, failure = failure)
}

# Keep the cache in step with one result: a clean model is (re)cached, and a
# model with an error loses any entry it had, so it is re-checked next build.
.modeldbAccRecord <- function(acc, res) {
  if (is.null(acc)) {
    return(invisible(NULL))
  }
  if (is.null(res$failure) && is.null(res$errors)) {
    acc$entries[[res$file]] <- list(hash = res$hash, row = res$row)
  } else {
    acc$entries[[res$file]] <- NULL
  }
  invisible(NULL)
}

# Parse `files` (relative to `dir`) in this R session. Runs in the build
# session for serial builds and inside each worker process for parallel ones.
.modeldbProcessFiles <- function(dir, files, hashes, acc = NULL) {
  out <- vector("list", length(files))
  for (i in seq_along(files)) {
    message("modeldb cache miss: ", files[i])
    loaded <- tryCatch(.modeldbFileRow(dir = dir, file = files[i]), error = identity)
    out[[i]] <-
      if (inherits(loaded, "error")) {
        .modeldbFileResult(files[i], hashes[i], failure = conditionMessage(loaded))
      } else {
        rownames(loaded$row) <- NULL
        .modeldbFileResult(files[i], hashes[i], row = loaded$row, errors = loaded$errors)
      }
    .modeldbAccRecord(acc, out[[i]])
  }
  out
}

# Largest number of models one worker process parses; a worker parses one
# chunk and exits. A guard, not a measured need: the parse-and-check path
# shows no per-model DLL leak (a session that parsed all 4,662 models ended
# with 38 DLLs), but this bounds whatever rxode2 state does accumulate.
.modeldbChunkMax <- 200L
# Smallest chunk worth a process: starting a worker and loading the package
# costs about as much as parsing 20 models.
.modeldbChunkMin <- 20L

# Split `n` cache misses into chunks for at most `workers` processes. Returns
# the chunk index of each miss (in file order) and the number of processes to
# start; one process means "parse in this session".
.modeldbChunkPlan <- function(n, workers, chunkSize = NULL) {
  if (is.null(chunkSize)) {
    chunkSize <- min(.modeldbChunkMax, max(.modeldbChunkMin, as.integer(ceiling(n / workers))))
  }
  chunk <- (seq_len(n) - 1L) %/% as.integer(chunkSize) + 1L
  list(chunk = chunk, workers = min(as.integer(workers), max(chunk, 0L)))
}

.modeldbRunMisses <- function(dir, files, hashes, workers, acc, chunkSize = NULL) {
  plan <- .modeldbChunkPlan(length(files), workers, chunkSize)
  if (plan$workers > 1L) {
    loader <- .modeldbWorkerLoader()
    if (!is.null(loader$serialReason)) {
      message(sprintf("modeldb: parsing %d model(s) in this session: %s", length(files), loader$serialReason))
      plan$workers <- 1L
    }
  }
  if (plan$workers <= 1L) {
    return(.modeldbProcessFiles(dir, files, hashes, acc))
  }
  .modeldbProcessParallel(dir, files, hashes, plan, acc, loader)
}

# Parse the misses chunk by chunk in up to `plan$workers` concurrent callr
# processes. Each process loads the package as `loader` says and is seeded
# with this session's parsed registers. Results come back in file order
# whatever order the chunks finish in. The cache in `acc` is updated as each
# chunk finishes, so an interrupt loses the chunks still running.
.modeldbProcessParallel <- function(dir, files, hashes, plan, acc, loader) {
  chunks <- split(seq_along(files), plan$chunk)
  registers <- .modeldbRegisterSnapshot()
  message(sprintf(
    "modeldb: parsing %d model(s) in %d chunk(s) on %d worker process(es)",
    length(files),
    length(chunks),
    plan$workers
  ))
  state <- new.env(parent = emptyenv())
  state$running <- list()
  # Stop the workers however this ends (an error here, an interrupt).
  on.exit(.modeldbKillWorkers(state), add = TRUE)
  results <- vector("list", length(chunks))
  nextChunk <- 1L
  while (nextChunk <= length(chunks) || length(state$running) > 0L) {
    while (length(state$running) < plan$workers && nextChunk <= length(chunks)) {
      idx <- chunks[[nextChunk]]
      state$running[[as.character(nextChunk)]] <- .modeldbStartWorker(
        loader = loader,
        dir = dir,
        files = files[idx],
        hashes = hashes[idx],
        registers = registers
      )
      nextChunk <- nextChunk + 1L
    }
    alive <- vapply(state$running, .modeldbWorkerAlive, logical(1))
    if (all(alive)) {
      Sys.sleep(0.2)
      next
    }
    for (k in names(state$running)[!alive]) {
      w <- state$running[[k]]
      state$running[[k]] <- NULL
      chunkResults <- .modeldbCollectWorker(w)
      for (res in chunkResults) {
        .modeldbAccRecord(acc, res)
      }
      results[[as.integer(k)]] <- chunkResults
    }
  }
  unlist(results, recursive = FALSE)
}

.modeldbStartWorker <- function(loader, dir, files, hashes, registers) {
  log <- tempfile("modeldb-worker-", fileext = ".log")
  proc <- callr::r_bg(
    func = .modeldbWorkerMain,
    args = list(loader = loader, dir = dir, files = files, hashes = hashes, registers = registers),
    stdout = log,
    stderr = "2>&1",
    package = FALSE
  )
  list(proc = proc, log = log, files = files, hashes = hashes)
}

.modeldbWorkerAlive <- function(w) {
  w$proc$is_alive()
}

.modeldbKillWorkers <- function(state) {
  for (w in state$running) {
    try(w$proc$kill(), silent = TRUE)
    unlink(w$log)
  }
  state$running <- list()
  invisible(NULL)
}

# Relay a finished worker's output and return its per-file results. A worker
# that died (a crash, the DLL limit) fails every file of its chunk, with the
# tail of its output, so the files are named and the build stops on them.
.modeldbCollectWorker <- function(w) {
  out <- tryCatch(w$proc$get_result(), error = identity)
  logLines <- if (file.exists(w$log)) readLines(w$log, warn = FALSE) else character()
  unlink(w$log)
  if (length(logLines) > 0L) {
    message(paste(logLines, collapse = "\n"))
  }
  if (!inherits(out, "error")) {
    return(out)
  }
  msg <- paste0(
    "worker process failed: ",
    conditionMessage(out),
    if (length(logLines) > 0L) {
      paste0("\nlast output:\n", paste(utils::tail(logLines, 10L), collapse = "\n"))
    }
  )
  mapply(
    .modeldbFileResult,
    file = w$files,
    hash = w$hashes,
    MoreArgs = list(failure = msg),
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )
}

# How a worker loads nlmixr2lib so that it runs the code this session runs.
# `codeEnv` is where the build's own functions live:
# - "source": a pkgload development namespace (devtools::document(),
#   load_all()). Each worker re-runs pkgload::load_all() on that tree as it
#   starts, so a file edited during a build reaches the workers started after
#   the edit.
# - "installed": the nlmixr2lib namespace loaded from an installed library;
#   workers load the package from that library.
# Anything else -- code sourced into a plain environment, as
# roxygenise(load_code = roxygen2::load_source) does, an installed copy no
# longer on disk, or no callr -- gives no way to run the same code in another
# process, and `serialReason` says why the misses stay in this session.
.modeldbWorkerLoader <- function(
  codeEnv = environment(.modeldbWorkerLoader),
  isDev = NULL,
  path = NULL,
  haveCallr = .modeldbHaveCallr()
) {
  if (!haveCallr) {
    return(list(serialReason = "the 'callr' package is not installed"))
  }
  if (!isNamespace(codeEnv) || !identical(unname(getNamespaceName(codeEnv)), "nlmixr2lib")) {
    return(list(
      serialReason = paste(
        "the build code is not a loaded nlmixr2lib namespace (it was sourced,",
        "e.g. by roxygenise(load_code = roxygen2::load_source)), so no worker",
        "process could run the same code"
      )
    ))
  }
  if (is.null(isDev)) {
    isDev <- .modeldbIsDevNamespace()
  }
  if (is.null(path)) {
    path <- getNamespaceInfo(codeEnv, "path")
  }
  if (!isDev && !file.exists(file.path(path, "Meta", "package.rds"))) {
    return(list(serialReason = paste0("the installed nlmixr2lib this session loaded is no longer at ", path)))
  }
  list(type = if (isDev) "source" else "installed", path = path)
}

.modeldbIsDevNamespace <- function() {
  requireNamespace("pkgload", quietly = TRUE) && pkgload::is_dev_package("nlmixr2lib")
}

.modeldbHaveCallr <- function() {
  requireNamespace("callr", quietly = TRUE)
}

# Runs in a callr process, whose copy of this function has the global
# environment as its enclosure, so it reaches the package only through `::`
# and asNamespace(). Coverage tools do not trace that process; the parallel
# test in test-modeldb.R checks its results against a serial build.
# nocov start
.modeldbWorkerMain <- function(loader, dir, files, hashes, registers) {
  if (identical(loader$type, "source")) {
    pkgload::load_all(loader$path, quiet = TRUE, helpers = FALSE, attach_testthat = FALSE)
  } else {
    loadNamespace("nlmixr2lib", lib.loc = dirname(loader$path))
  }
  ns <- asNamespace("nlmixr2lib")
  get(".modeldbSeedRegisters", envir = ns)(registers)
  get(".modeldbProcessFiles", envir = ns)(dir, files, hashes)
}
# nocov end

# Parse the three canonical-name registers afresh. The checker memoises them
# per session, which is right while checking models one at a time but would
# let a long-lived session build against a register edited since.
.modeldbReloadRegisters <- function() {
  .loadCanonicalCovariates(force = TRUE)
  .loadCanonicalParameters(force = TRUE)
  .loadCanonicalCompartments(force = TRUE)
  invisible(NULL)
}

.modeldbRegisterSnapshot <- function() {
  list(
    canonical = .loadCanonicalCovariates(),
    parameters = .loadCanonicalParameters(),
    compartments = .loadCanonicalCompartments()
  )
}

.modeldbSeedRegisters <- function(registers) {
  .covariateRegisterCache$canonical <- registers$canonical
  .covariateRegisterCache$aliasMap <- NULL
  .covariateRegisterCache$parameters <- registers$parameters
  .covariateRegisterCache$compartments <- registers$compartments
  invisible(NULL)
}

# rxode2::getRxThreads() is half the cores by default and already honours
# setRxThreads(), OMP_THREAD_LIMIT and R CMD check's two-core limit.
.modeldbDefaultWorkers <- function(threads = rxode2::getRxThreads()) {
  min(8L, max(1L, as.integer(threads) %/% 2L))
}

.modeldbCheckWorkers <- function(workers) {
  checkmate::assertCount(
    workers,
    positive = TRUE,
    .var.name = "workers (or option nlmixr2lib.buildModelDb.workers)"
  )
  as.integer(workers)
}

# Build the global signature used to invalidate all entries at once. Changes to
# the extraction code (R/modeldb.R) or to any upstream parser version force a
# full rebuild. So do changes to the convention checker and to the
# canonical-name registers: a cached model is not re-checked, so its cached
# verdict holds only for the checker and registers it was checked against.
# Silent upstream breakage without a version bump is not covered here; the
# escape hatch is `unlink("data-raw/modeldb-cache.rds")`.
.modeldbGlobalSig <- function(packageDirectory) {
  list(
    cacheVersion = 3L,
    rVersion = R.version.string,
    nlmixr2 = .modeldbSafeVersion("nlmixr2"),
    nlmixr2est = .modeldbSafeVersion("nlmixr2est"),
    rxode2 = .modeldbSafeVersion("rxode2"),
    modeldbRhash = unname(tools::md5sum(file.path(packageDirectory, "R/modeldb.R"))),
    checkerHash = .modeldbFileHashes(
      file.path(packageDirectory, c("R/checkModelConventions.R", "R/conventions.R"))
    ),
    registerHash = .modeldbFileHashes(
      .collateGlob(file.path(packageDirectory, "inst/references/*.md"))
    )
  )
}

.modeldbSafeVersion <- function(pkg) {
  tryCatch(
    as.character(utils::packageVersion(pkg)),
    error = function(e) NA_character_
  )
}

# md5 of each file, named by its basename; NA for a file that does not exist.
.modeldbFileHashes <- function(paths) {
  stats::setNames(unname(tools::md5sum(paths)), basename(paths))
}

.modeldbCacheRead <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  tryCatch(
    {
      cache <- readRDS(path)
      if (!is.list(cache) || !all(c("globalSig", "entries") %in% names(cache))) {
        return(NULL)
      }
      cache
    },
    error = function(e) NULL
  )
}

# Write through a temporary file in the cache's own directory and rename it
# into place, so a reader never sees a half-written cache. The temporary name
# is unique, so two builds of one tree cannot write the same temporary file,
# and it is removed whether or not the write succeeds.
.modeldbCacheWrite <- function(path, cache) {
  tmp <- tempfile("modeldb-cache-", tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  res <- tryCatch(
    {
      dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
      saveRDS(cache, file = tmp)
      file.rename(tmp, path)
    },
    warning = identity,
    error = identity
  )
  if (inherits(res, "condition")) {
    warning("Could not write modeldb cache to ", path, ": ", conditionMessage(res), call. = FALSE)
  } else if (!isTRUE(res)) {
    warning("Could not write modeldb cache to ", path, ": the rename into place failed", call. = FALSE)
  }
  invisible(NULL)
}

# Parse a model filename like "Bajaj_2017_nivolumab" into a friendlier
# label "Nivolumab (Bajaj 2017)". Year-letter collision suffixes
# ("Hansson_2013a_sunitinib") are passed through verbatim so each entry
# disambiguates from its siblings. Filenames whose first two
# underscore-separated tokens are not "<Author>_<4-digit year>[letter]"
# fall back to the bare name with underscores replaced by spaces. The
# "NA_NA_<drug>" convention used by some ddmore models is prefixed with
# "DDMoRe:" so the parameterless entries are clearly attributable
# wherever the label appears (vignette title, list-of-models name
# column, navbar dropdown).
.parseModelLabel <- function(filename) {
  base <- tools::file_path_sans_ext(basename(filename))
  parts <- strsplit(base, "_", fixed = TRUE)[[1]]
  if (
    length(parts) >= 3 &&
      grepl("^[0-9]{4}[a-z]?$", parts[2]) &&
      parts[1] != "NA" &&
      parts[2] != "NA"
  ) {
    drug <- paste(parts[-(1:2)], collapse = " ")
    drug <- paste0(toupper(substr(drug, 1, 1)), substr(drug, 2, nchar(drug)))
    return(sprintf("%s (%s %s)", drug, parts[1], parts[2]))
  }
  if (length(parts) >= 3 && parts[1] == "NA" && parts[2] == "NA") {
    return(sprintf("DDMoRe: %s", paste(parts[-(1:2)], collapse = " ")))
  }
  gsub("_", " ", base, fixed = TRUE)
}

# Coarse-grained category derived from a model file's path. Matches the
# inst/modeldb/ directory layout: "specificDrugs/<file>.R" -> "specificDrugs",
# "ddmore/<file>.R" -> "ddmore", everything else (top-level templates and the
# remaining therapeutic-area / endogenous / pharmacokinetics directories) ->
# "other". `filename` here may be the relative path stored on the modeldb row
# (e.g. "specificDrugs/Aguiar_2021_ustekinumab.R") or an absolute filesystem
# path -- the prefix match is anchored on the canonical directory names.
.parseModelCategory <- function(filename) {
  if (grepl("(^|/)specificDrugs/", filename)) {
    return("specificDrugs")
  }
  if (grepl("(^|/)ddmore/", filename)) {
    return("ddmore")
  }
  "other"
}

# Refresh the auto-generated regions of the package's `_pkgdown.yml` by
# rewriting only the lines between AUTOGEN signpost comments. Three
# regions are managed:
#
#   * `# AUTOGEN:specific_drugs:BEGIN` / `:END` -- the menu items under
#     `navbar.components.specific_drugs.menu`. One entry per specificDrugs/
#     model with a vignette.
#   * `# AUTOGEN:ddmore:BEGIN` / `:END` -- the menu items under
#     `navbar.components.ddmore.menu`. One entry per ddmore/ model with a
#     vignette.
#   * `# AUTOGEN:articles:BEGIN` / `:END` -- the top-level `articles:`
#     groups. A visible "General" group listing every top-level
#     `vignettes/*.Rmd`, and an `internal: true` group listing every
#     `vignettes/articles/*.Rmd` -- so the drug-specific pages remain
#     buildable and URL-addressable but are absent from the auto-generated
#     Articles index (they are reached from the navbar dropdowns instead).
#
# Anything outside the signposts -- template, url, navbar.structure, the
# `text:` lines of the dropdowns, comments, maintainer-added components,
# blank lines, formatting -- is preserved byte-for-byte. This avoids the
# yaml round-trip pitfalls of `yaml::read_yaml` + `yaml::write_yaml`
# (lost comments, reordered keys, boolean reformatting, broken style),
# and makes a maintainer-edited region trivially distinguishable from
# the auto-generated parts.
#
# If a signpost block is missing, the function logs a message and skips
# only that region; the others still refresh.
.writePkgdownNavbar <- function(modeldb, packageDirectory) {
  ymlPath <- file.path(packageDirectory, "_pkgdown.yml")
  if (!file.exists(ymlPath)) {
    return(invisible())
  }

  with_vignette <- !is.na(modeldb$vignette)
  spec <- modeldb[with_vignette & modeldb$category == "specificDrugs", , drop = FALSE]
  ddmo <- modeldb[with_vignette & modeldb$category == "ddmore", , drop = FALSE]

  # Each replacement line is emitted at column 0 (no leading whitespace).
  # .replaceAutogen() prefixes whatever indentation the BEGIN signpost
  # carries, so a marker at any nesting level renders correctly.
  buildMenuLines <- function(rows) {
    if (nrow(rows) == 0) {
      return(character())
    }
    ord <- .collateOrder(rows$label)
    rows <- rows[ord, , drop = FALSE]
    unlist(lapply(seq_len(nrow(rows)), function(i) {
      c(
        sprintf("- text: %s", .yamlQuote(rows$label[[i]])),
        sprintf("  href: articles/%s.html", rows$vignette[[i]])
      )
    }))
  }

  # Top-level (visible) vignettes and the drug-specific articles list.
  # Read straight from disk -- the modeldb only knows about model files,
  # not vignettes -- and stay consistent with what pkgdown will build.
  vignDir <- file.path(packageDirectory, "vignettes")
  topNames <- character()
  articleNames <- character()
  # The names are ordered after the extension is stripped, because stripping it
  # can reorder them ("a-b.Rmd" < "a.Rmd", but "a" < "a-b"), and that is the
  # order the file carries. The enumeration is collated too, so that no
  # ordering in this file is left to a judgement call about whether it matters.
  if (dir.exists(vignDir)) {
    topNames <- .collateSort(tools::file_path_sans_ext(
      .collateListFiles(vignDir, pattern = "\\.Rmd$")
    ))
    articleDir <- file.path(vignDir, "articles")
    if (dir.exists(articleDir)) {
      articleNames <- .collateSort(tools::file_path_sans_ext(
        .collateListFiles(articleDir, pattern = "\\.Rmd$")
      ))
    }
  }
  articlesLines <- character()
  if (length(topNames) > 0) {
    articlesLines <- c(
      articlesLines,
      "- title: General",
      "  desc: Cross-cutting guides and the list of models.",
      "  contents:",
      sprintf("  - %s", topNames)
    )
  }
  if (length(articleNames) > 0) {
    articlesLines <- c(
      articlesLines,
      "- title: internal",
      paste(
        "  desc: Drug-specific validation vignettes. They are reached from the",
        "Specific drug models and DDMoRe models navbar dropdowns; this group",
        "hides them from the main Articles index."
      ),
      "  contents:",
      sprintf("  - articles/%s", articleNames),
      "  internal: true"
    )
  }

  ymlLines <- readLines(ymlPath, encoding = "UTF-8", warn = FALSE)
  result <- .replaceAutogen(ymlLines, "specific_drugs", buildMenuLines(spec))
  result <- .replaceAutogen(result, "ddmore", buildMenuLines(ddmo))
  result <- .replaceAutogen(result, "articles", articlesLines)

  if (!identical(result, ymlLines)) {
    writeLines(result, ymlPath, useBytes = FALSE, sep = "\n")
  }
  message(
    "Refreshing navbar in ",
    ymlPath,
    " (specificDrugs=",
    nrow(spec),
    ", ddmore=",
    nrow(ddmo),
    ", internal-articles=",
    length(articleNames),
    ")"
  )
  invisible()
}

# Replace the block of lines between `# AUTOGEN:<name>:BEGIN` and
# `# AUTOGEN:<name>:END` markers with `replacement`. Each element of
# `replacement` is one output line, emitted at column 0; the leading
# whitespace of the BEGIN marker is prefixed onto every replacement
# line so a marker at any YAML nesting level renders with correct
# indentation. The marker lines themselves are preserved.
#
# If either marker is missing or malformed, the function leaves the
# lines untouched and warns once per region.
.replaceAutogen <- function(lines, name, replacement) {
  beginRx <- sprintf("^([[:space:]]*)#[[:space:]]*AUTOGEN:%s:BEGIN", name)
  endRx <- sprintf("^[[:space:]]*#[[:space:]]*AUTOGEN:%s:END", name)
  bi <- grep(beginRx, lines)
  ei <- grep(endRx, lines)
  if (length(bi) != 1L || length(ei) != 1L || ei <= bi) {
    warning(
      "AUTOGEN:",
      name,
      " signposts missing or malformed in _pkgdown.yml; ",
      "leaving that region unchanged.",
      call. = FALSE
    )
    return(lines)
  }
  indent <- sub(beginRx, "\\1", lines[bi])
  if (length(replacement) > 0) {
    replacement <- paste0(indent, replacement)
    # Drop the trailing whitespace from any line that ended up bare-indent,
    # so the file stays tidy without trailing spaces on otherwise-empty lines.
    replacement <- sub("[[:space:]]+$", "", replacement)
  }
  c(lines[seq_len(bi)], replacement, lines[ei:length(lines)])
}

# Minimal YAML-scalar quoter for the menu `text:` field. Adds double
# quotes when the label contains a leading colon, leading dash, hash,
# leading-or-trailing whitespace, or any of YAML's flow indicators that
# would otherwise need parsing care. Internal double quotes are escaped.
# The fields we emit are short navbar labels like
# `Ustekinumab (Aguiar 2021)` and `DDMoRe: lidocaine`; the latter has a
# leading-colon-like prefix that would parse as a mapping without
# quotes, so the conservative default is to always quote.
.yamlQuote <- function(x) {
  sprintf('"%s"', gsub('"', '\\\\"', x, fixed = TRUE))
}
