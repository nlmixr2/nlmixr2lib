#' Model library for nlmixr2
#'
#' This is a data frame of the available models in nlmixr2lib, it is generated
#' with the package.  Custom modeldb may be used.
#'
#' @eval buildModelDb()
"modeldb"

utils::globalVariables(c("modeldb"))

## nocov start

#' Build the default model database used in nlmixr2lib
#'
#' The function is only used during the build process, and it is not needed for
#' general use.  See the documentation of "modeldb" for its output.
#'
#' @return The text needed for documenting the "modeldb" object (see the modeldb
#'   roxygen eval statement)
#' @noRd
buildModelDb <- function() {
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
  cachePath <- file.path(packageDirectory, "data-raw/modeldb-cache.rds")
  cache <- .modeldbCacheRead(cachePath)
  sig <- .modeldbGlobalSig(packageDirectory)
  if (is.null(cache) || !identical(cache$globalSig, sig)) {
    if (!is.null(cache)) {
      message("modeldb cache invalidated (global signature changed)")
    }
    cache <- list(globalSig = sig, entries = list())
  }
  result <- .addDirToModelDbCached(
    dir = file.path(packageDirectory, "inst/modeldb"),
    entries = cache$entries
  )
  modeldb <- result$modeldb
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
  qs2::qs_save(modeldb, file = file.path(packageDirectory, "inst/modeldb.qs2"))
  .modeldbCacheWrite(cachePath, list(globalSig = sig, entries = result$entries))
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
    list.files(
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
  fileName <- file.path(dir, file)

  # Extract the model from the file
  parsedFile <- parse(file = fileName)
  stopifnot(identical(parsedFile[[1]][[1]], as.name("<-")))

  modelName <- as.character(parsedFile[[1]][[2]])
  packageStartupMessage("Loading ", modelName, " from ", fileName)
  if (modelName != tools::file_path_sans_ext(basename(file))) {
    stop("Loading model failed due to filename/modelName mismatch: ", fileName,
      call. = FALSE
    ) # nocov
  }

  # Parse the model to get the fixed effects and DV parameters
  mod <- nlmixr2est::nlmixr(eval(parsedFile))

  # Convention check. Reports deviations but does not halt the build so that
  # grandfathered models continue to regenerate while their issues surface.
  issues <- tryCatch(
    suppressWarnings(checkModelConventions(mod, verbose = FALSE)),
    error = function(e) NULL
  )
  if (!is.null(issues) && nrow(issues) > 0) {
    n_err <- sum(issues$severity == "error")
    n_warn <- sum(issues$severity == "warning")
    if (n_err + n_warn > 0) {
      message(sprintf(
        "  %s: %d convention error(s), %d warning(s) - run checkModelConventions(\"%s\") for details",
        modelName, n_err, n_warn, modelName
      ))
    }
  }

  description <- mod$meta$description
  if (is.null(description)) {
    message("No description for model in ", fileName)
    description <- NA_character_
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
      name        = modelName,
      description = description,
      parameters  = paste(modParamFixed, collapse = ","),
      DV          = paste(paramErr, collapse = ","),
      linCmt      = mod$props$linCmt,
      algebraic   = algebraic,
      dosing      = dosing,
      depends     = depends,
      filename    = fileName
    )
  modeldb <- rbind(modeldb, ret)
  if (any(duplicated(modeldb$name))) {
    stop("Duplicated model name: ", modeldb$name[duplicated(modeldb$name)]) # nocov
  }
  modeldb
}
## nocov end

# Walk a modeldb directory and return rows + updated cache entries. Reuses a
# cached row when the file's md5 matches the cached hash; otherwise re-parses
# via addFileToModelDb(). Files not present on disk are dropped from entries
# by virtue of iterating only over live files.
.addDirToModelDbCached <- function(dir, entries) {
  filesToLoad <-
    list.files(
      path = dir,
      pattern = "\\.R$",
      ignore.case = TRUE,
      recursive = TRUE
    )
  modeldb <- data.frame()
  newEntries <- list()
  for (currentFile in filesToLoad) {
    fullPath <- file.path(dir, currentFile)
    fileHash <- unname(tools::md5sum(fullPath))
    cached <- entries[[currentFile]]
    if (!is.null(cached) && identical(cached$hash, fileHash)) {
      message("modeldb cache hit: ", currentFile)
      row <- cached$row
      row$filename <- fullPath
      modeldb <- rbind(modeldb, row)
      newEntries[[currentFile]] <- cached
    } else {
      message("modeldb cache miss: ", currentFile)
      before <- nrow(modeldb)
      modeldb <- addFileToModelDb(dir = dir, file = currentFile, modeldb = modeldb)
      newRow <- modeldb[before + 1L, , drop = FALSE]
      rownames(newRow) <- NULL
      newEntries[[currentFile]] <- list(hash = fileHash, row = newRow)
    }
  }
  list(modeldb = modeldb, entries = newEntries)
}

# Build the global signature used to invalidate all entries at once. Changes to
# the extraction code (R/modeldb.R) or to any upstream parser version force a
# full rebuild. Silent upstream breakage without a version bump is not covered
# here; the escape hatch is `unlink("data-raw/modeldb-cache.rds")`.
.modeldbGlobalSig <- function(packageDirectory) {
  safeVersion <- function(pkg) {
    tryCatch(
      as.character(utils::packageVersion(pkg)),
      error = function(e) NA_character_
    )
  }
  list(
    cacheVersion = 1L,
    rVersion     = R.version.string,
    nlmixr2      = safeVersion("nlmixr2"),
    nlmixr2est   = safeVersion("nlmixr2est"),
    rxode2       = safeVersion("rxode2"),
    modeldbRhash = unname(tools::md5sum(file.path(packageDirectory, "R/modeldb.R")))
  )
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

.modeldbCacheWrite <- function(path, cache) {
  reportFailure <- function(c) {
    warning("Could not write modeldb cache to ", path, ": ", conditionMessage(c))
  }
  tryCatch(
    {
      dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
      tmp <- paste0(path, ".tmp")
      saveRDS(cache, file = tmp)
      file.rename(tmp, path)
    },
    warning = reportFailure,
    error = reportFailure
  )
  invisible(NULL)
}
