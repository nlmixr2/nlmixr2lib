#' To remove transit compartments from the model
#' @param model The model as a function
#' @param transit The number of transit compartments to remove
#' @export
#' @examples
#' 
#' # In this example the transit is added and then a few are removed
#' 
#' readModelDb("PK_1cmt_des") |>
#'    addTransit(4)  |>
#'    removeTransit(3)
removeTransit <- function(model,transit,central="central",depot="depot",transitComp ="transit", ktr="ktr"){
  checkmate::assertCharacter(central, pattern= "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",len=1,any.missing = FALSE,min.chars = 1)
  checkmate::assertCharacter(depot, pattern= "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",len=1,any.missing = FALSE,min.chars = 1)
  checkmate::assertCharacter(transitComp, pattern= "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",len=1,any.missing = FALSE,min.chars = 1)
  checkmate::assertCharacter(ktr, pattern= "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",len=1,any.missing = FALSE,min.chars = 1)
  
  if(!missing(transit)){
    checkmate::assertIntegerish(transit, lower=1,any.missing = FALSE,len = 1)
    
  }
  
  temp  <- rxode2::assertRxUi(model)
  
  #browser()
  mv <- rxode2::rxModelVars(temp)
  
  
  if (!(central %in% mv$state)){
    stop("'",central,"' needs to be in the model")
  }
  if (!(any(grepl("^transit",mv$state)))){
    stop("'",transitComp," need to be in the model")
  }
  if(!(depot %in% mv$state)){
    stop("'",depot,"' needs to be in the model")
  }
  #Extract model
  modelNew <- rxode2::modelExtract(temp,endpoint=NA)
  
  #modify ODE for central compartment
  center<- eval(str2lang(paste0("rxode2::modelExtract(temp,d/dt(",central,"),lines=TRUE)")))
  rhs <- sub(".*<-\\s*","",center)
  rhs <- sub(paste0("\\s*",ktr,"\\d*\\s*\\*\\s*",transitComp,"\\d*"),"",rhs)
  
  #Find total number of transit compartments in the model
  obj <- c(unlist(modelNew)[which(grepl("\\s*^ktr",mv$lhs))])
  totalTransit <- length(obj)
  
  if(missing(transit)){
    transit <- totalTransit
    line <- str2lang(paste0("d/dt(",central,") <- ",rhs))
  }else{
    line <- str2lang(paste0("d/dt(",central,") <- ",ktr,(totalTransit - transit),"*",transitComp,(totalTransit - transit),deparse(str2lang(rhs))))
  }
  
  
  #Modify ini{} 
  temp2<-temp$iniDf
  temp3 <- temp2$name
  temp3 <- rev(grep("^lktr([1-9]|[10-99]+)$",temp3,value=TRUE))
  temp3 <- temp3[1:transit]
  temp2 <- temp2[!temp2$name %in% temp3, ]
  rxode2::ini(temp)<-temp2
  
  #Modify model{}
  obj <- c()
  indices <- totalTransit:(totalTransit - transit + 1)
  obj <- unlist(lapply(indices, function(i) {
    obj1 <- eval(str2lang(paste0("rxode2::modelExtract(temp,",ktr,i,",lines = TRUE)")))
    obj2 <- eval(str2lang(paste0("rxode2::modelExtract(temp, d/dt(transit", i, "), lines = TRUE)")))
    obj <- c(obj1, obj2)
  }))
 
  for (i in obj){
    index <- which(modelNew==i)
    modelNew <- modelNew[-index]
  }
  
  rxode2::model(temp)<-modelNew
  
  #modify ODE for central compartment
  # center<- eval(str2lang(paste0("rxode2::modelExtract(temp,d/dt(",central,"),lines=TRUE)")))
  # rhs <- sub(".*<-\\s*","",center)
  # rhs <- sub(paste0("\\s*",ktr,"\\d*\\s*\\*\\s*",transitComp,"\\d*"),"",rhs)
  # line <- str2lang(paste0("d/dt(",central,") <- ",ktr,(totalTransit - transit),"*",transitComp,(totalTransit - transit),deparse(str2lang(rhs))))
  temp <- rxode2::model(temp, line)
  
  #return 
  temp
  
  
 }
