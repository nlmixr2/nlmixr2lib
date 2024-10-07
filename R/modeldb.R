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
  modeldb <- addDirToModelDb(file.path(packageDirectory, "inst/modeldb"))
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
  message("Done saving the modeldb to ", savefile)

  colDesc <-
    list(
      name = "Model name that can be used to extract the model from the model library",
      description = "Model description in free from text; in model itself",
      parameters  = "A comma separated string listing either the parameter in the model defined by population/individual effects or a population effect parameter",
      DV          = "The definition of the dependent variable(s)",
      linCmt      = "Logical flag indicating if solved models are used (TRUE) or not (FALSE)",
      algebraic   = "Logical flag indicating if the model is purely algebraic: TRUE no linCmt() and no ODEs; FALSE otherwise",
      dosing      = "A comma separated string of identified dosing compartments",
      depends     = "A comma separated string of objects the model depends on",
      filename    = "Filename of the model.  By default these are installed in the model library and read on demand"
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

  description <- mod$meta$description
  if (is.null(description)) {
    message("No description for model in ", fileName)
    description <- NA_character_
  }

  # Finding dosing
  dosing <- NULL
  dosing_meta <- mod$meta$dosing
  if(!is.null(dosing_meta)){
    dosing <- paste(dosing_meta, collapse=",")
  }else {
    if("depot" %in% mod$props$cmt) {
      dosing <- c(dosing, "depot")
    }
    if("central" %in% mod$props$cmt) {
      dosing <- c(dosing, "central")
    }
    if(!is.null(dosing)) {
      dosing <- paste(dosing, collapse=",")
    } else {
      dosing <- NA_character_
    }
  }

  # Finding depends
  depends <- NULL
  depends_meta <- mod$meta$depends
  if(!is.null(depends_meta)){
    depends = paste(depends_meta, collapse=",")
  }
  if(is.null(depends)){
    depends = NA_character_
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


  if(!mod$props$linCmt && (length(mod$props$cmt)  == 0)){
    algebraic = TRUE
  } else {
    algebraic = FALSE
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
