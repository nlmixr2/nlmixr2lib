#' Model library for nlmixr2
#'
#' This is a data frame of the availble models in nlmixr2lib, it is generated at compile time
#'
#'
#' @format A data frame with XXX rows and 5 columns
#' \describe{
#'   \item{name}{Model name that can be used to extract the model from the model library}
#'   \item{description}{Model description in free from text; in model itself}
#'   \item{parameters}{A comma separated string listing either the parameter in the model defined by population/individual effects or a population effect parameter}
#'   \item{DV}{The definition of the dependent variable(s)}
#'   \item{filename}{Filename of the model.  By default these are installed in the model library and read on demand}
#' }
"modeldb"

utils::globalVariables(c("modeldb"))

#' Add a directory to the modeldb
#'
#' @param dir Directory name containing model files
#' @param modeldb The starting modeldb data.frame
#' @return The updated modeldb data.frame
#' @export
addDirToModelDb <- function(dir, modeldb=data.frame()) {
  filesToLoad <-
    list.files(
      path = dir,
      pattern = "\\.R$",
      ignore.case = TRUE
    )
  for (currentFile in filesToLoad) {
    message("parse currentFile")
    addFileToModelDb(dir = dir, file = currentFile)
  }
  modeldb
}

#' @describeIn addDirToModelDb Add a file to the modeldb
#' @param file The file name (without the directory name)
#' @export
addFileToModelDb <- function(dir, file, modeldb) {

  fileName <- file.path(dir, file)
  # Extract the description from the first line of the file
  desc <- readLines(con = fileName, n = 1)
  descClean <- gsub(x = desc, pattern = "^# *Description: *", replacement = "")

  # Extract the model from the file
  parsedFile <- parse(file = fileName)
  stopifnot(identical(parsedFile[[1]][[1]], as.name("<-")))

  modelName <- as.character(parsedFile[[1]][[2]])
  packageStartupMessage("Loading ", modelName, " from ", fileName)
  if (modelName != tools::file_path_sans_ext(file)) {
    stop("Loading model failed due to filename/modelName mismatch: ", fileName,
         call.=FALSE) # nocov
  }

  # Parse the model to get the fixed effects and DV parameters
  mod <- nlmixr2est::nlmixr(eval(parsedFile[[1]][[3]]))

  description <- mod$meta$description
  if (is.null(description)) {
    message("No description for model in ", fileName)
    description <- NA_character_
  }

  # Extract the parameter names
  modParam <- mod$iniDf
  # Fixed effects
  modParamFixed <- modParam$name[is.na(modParam$neta1) & is.na(modParam$err)]

  # swap modeled parameter names for the mu-ref parameter names, where
  # applicable
  for (nm in names(mod$getSplitMuModel$pureMuRef)) {
    modParamFixed[modParamFixed %in% nm] <- mod$getSplitMuModel$pureMuRef[[nm]]
  }

  # Error model
  paramErr <- mod$predDf$cond
  if ("rxLinCmt" %in% paramErr) {
    paramErr[paramErr %in% "rxLinCmt"] <- "linCmt()"
  }

  ret <-
    data.frame(
      name=modelName,
      description=description,
      parameters=paste(modParamFixed, collapse = ","),
      DV=paramErr,
      filename = gsub("inst[/\\]", "", fileName)
    )
  modeldb <- rbind(modeldb, ret)
  if (any(duplicated(modeldb$name))) {
    stop("Duplicated model name: ", modeldb$name[duplicated(modeldb$name)]) # nocov
  }
  ret
}
