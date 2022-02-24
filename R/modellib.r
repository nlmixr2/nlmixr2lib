#------------------------------------------ modellib ------------------------------------------
#' Get the model from the model library
#'
#' This function gets a model from the available model library
#'
#' @param model character with the name of the model to load (if NULL, lists all available base models)
#' @param iiv vector with the parameters to add IIV on
#' @param reserr character with the type of residual error (currently "add","prop" and "add+prop" are accepted)
#' @param rds character with the basename of the rds which include the model library
#' @param loc character with the location of the rds database and models
#'
#' @details This is a very first draft just to look at the proof of concept
#'
#' @return The function returns a character string with the model code
#'
#' @export
#' @examples
#'
#' \dontrun{
#'   modellib(model="PK_1cmt")
#'   modellib(model="PK_1cmt",iiv = c("ka","v"),reserr = "add")
#'   modellib(model="PK_1cmt",reserr = "add")
#' }
modellib <- function(model=NULL,iiv=NULL,reserr=NULL,rds="modellib.rds",loc=system.file(package = "nlmixr.lib")){
  modeldb <- try(readRDS(paste0(loc,"/",rds)))
  if("try-error"%in%class(modeldb)) stop("model database could not be loaded")
  if(is.null(model)){
    cat(paste0(modeldb$name," (",modeldb$description,")"),sep="\n")
  }else if(!model%in%modeldb$name){
    stop("selected model not in database")
  }else{
    modeldb <- modeldb[modeldb$name==model,]
    modr    <- readLines(paste0(loc,"/",model,".r"))
    inib <- modb <- inib2 <- modb2 <- NULL
    if(is.null(iiv)){
      modr <- modr
    }else if(!all(iiv%in%trimws(unlist(strsplit(modeldb$Parameters,split=","))))){
      warning("not all iiv parameters in model, iiv not added")
    }else{
      inib <- paste0("ini({",paste(paste0("eta.",iiv,"~",0.1),collapse = ";"),"})")
      modb <- paste0("model({",paste(paste0(iiv,"=","exp(",iiv," + eta.",iiv,")"),collapse = ";"),"})")
    }
    if(is.null(reserr) || reserr=="prop"){
      modr <- modr
    }else if(!reserr%in%c("add","prop","add+prop")){
      warning("unknown residual error, not added")
    }else if(reserr=="add"){
      inib2 <- "ini(add.err=0.1)"
      modb2 <- paste0("model({",modeldb$DV,"~add(add.err)})")
    }else if(reserr=="add+prop"){
      inib2 <- "ini(add.err=0.1)"
      modb2 <- paste0("model({",modeldb$DV,"~add(add.err)+prop(prop.err)})")
    }
    if(!is.null(inib) & is.null(inib2)){
      modr  <- c(modr,paste0(model," <- ",model," %>% ",modb," %>% ",inib))
    }else if(is.null(inib) & !is.null(inib2)){
      modr  <- c(modr,paste0(model," <- ",model," %>% ",modb2," %>% ",inib2))
    }else if(!is.null(inib) & !is.null(inib2)){
      modr  <- c(modr,paste0(model," <- ",model," %>% ",modb," %>% ",inib," %>% ",modb2," %>% ",inib2))
    }
    # currently a vector with the model code is returned
    # we could easily extend this, e.g.
    # to console: cat(modr,sep="\n)
    # to file: writeLines(modr,paste0(model,".r"))
    # directly evaluated: eval(parse(text=modr))
    # to Rstudio's current script: rstudioapi::insertText(paste(modr,"\n"))
    return(modr)
  }
}
