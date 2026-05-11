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
      name        = modelName,
      description = description,
      parameters  = paste(modParamFixed, collapse = ","),
      DV          = paste(paramErr, collapse = ","),
      linCmt      = mod$props$linCmt,
      algebraic   = algebraic,
      dosing      = dosing,
      depends     = depends,
      vignette    = vignette,
      label       = .parseModelLabel(fileName),
      category    = .parseModelCategory(fileName),
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
    cacheVersion = 3L,
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
  if (length(parts) >= 3 &&
      grepl("^[0-9]{4}[a-z]?$", parts[2]) &&
      parts[1] != "NA" && parts[2] != "NA") {
    drug <- paste(parts[-(1:2)], collapse = " ")
    drug <- paste0(toupper(substr(drug, 1, 1)),
                   substr(drug, 2, nchar(drug)))
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
  if (grepl("(^|/)specificDrugs/", filename)) return("specificDrugs")
  if (grepl("(^|/)ddmore/",        filename)) return("ddmore")
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
  if (!file.exists(ymlPath)) return(invisible())

  with_vignette <- !is.na(modeldb$vignette)
  spec <- modeldb[with_vignette & modeldb$category == "specificDrugs", , drop = FALSE]
  ddmo <- modeldb[with_vignette & modeldb$category == "ddmore",        , drop = FALSE]

  # Each replacement line is emitted at column 0 (no leading whitespace).
  # .replaceAutogen() prefixes whatever indentation the BEGIN signpost
  # carries, so a marker at any nesting level renders correctly.
  buildMenuLines <- function(rows) {
    if (nrow(rows) == 0) return(character())
    ord <- order(rows$label)
    rows <- rows[ord, , drop = FALSE]
    unlist(lapply(seq_len(nrow(rows)), function(i) {
      c(
        sprintf("- text: %s",            .yamlQuote(rows$label[[i]])),
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
  if (dir.exists(vignDir)) {
    topNames <- sort(tools::file_path_sans_ext(
      list.files(vignDir, pattern = "\\.Rmd$")
    ))
    articleDir <- file.path(vignDir, "articles")
    if (dir.exists(articleDir)) {
      articleNames <- sort(tools::file_path_sans_ext(
        list.files(articleDir, pattern = "\\.Rmd$")
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
  result <- .replaceAutogen(result,   "ddmore",         buildMenuLines(ddmo))
  result <- .replaceAutogen(result,   "articles",       articlesLines)

  if (!identical(result, ymlLines)) {
    writeLines(result, ymlPath, useBytes = FALSE, sep = "\n")
  }
  message("Refreshing navbar in ", ymlPath,
          " (specificDrugs=", nrow(spec),
          ", ddmore=", nrow(ddmo),
          ", internal-articles=", length(articleNames), ")")
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
  endRx   <- sprintf("^[[:space:]]*#[[:space:]]*AUTOGEN:%s:END",   name)
  bi <- grep(beginRx, lines)
  ei <- grep(endRx,   lines)
  if (length(bi) != 1L || length(ei) != 1L || ei <= bi) {
    warning("AUTOGEN:", name, " signposts missing or malformed in _pkgdown.yml; ",
            "leaving that region unchanged.", call. = FALSE)
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
