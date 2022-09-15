modeldb <- data.frame()

.onLoad <- function(libname, pkgname) {
  loadPath <- file.path(libname, pkgname, "inst")
  addDirToModelDb(dir = loadPath)
}

addDirToModelDb <- function(dir) {
  filesToLoad <-
    list.files(
      path = dir,
      pattern = "\\.R$",
      ignore.case = TRUE
    )
  for (currentFile in filesToLoad) {
    addFileToModelDb(dir = dir, file = currentFile)
  }
  modeldb
}

addFileToModelDb <- function(dir, file) {
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
    stop("Loading model failed due to filename/modelName mismatch: ", fileName) # nocov
  }

  # Parse the model to get the fixed effects and DV parameters
  mod <- rxode2::rxode(eval(parsedFile[[1]][[3]]))

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
      Parameters=paste(modParamFixed, collapse = ","),
      DV=paramErr,
      filename=fileName
    )
  modeldb <<- rbind(modeldb, ret)
  if (any(duplicated(modeldb$name))) {
    stop("Duplicated model name: ", modeldb$name[duplicated(modeldb$name)]) # nocov
  }
  ret
}
