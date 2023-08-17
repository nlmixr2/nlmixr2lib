#' To remove peripheral compartments from the model
#' @param model The model as a function
#' @param peripheral The number of peripheral compartments to remove
#' @examples
#' library(rxode2)
#' readModelDb("PK_1cmt") %>% removeComp(.,3)
removeComp <- function(model,peripheral,central="central",depot="depot",peripheralComp ="peripheral",vp="vp",vc="vc",q="q"){
  
  checkmate::assertCharacter(central, pattern= "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",len=1,any.missing = FALSE,min.chars = 1)
  checkmate::assertCharacter(depot, pattern= "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",len=1,any.missing = FALSE,min.chars = 1)
  checkmate::assertCharacter(peripheralComp, pattern= "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",len=1,any.missing = FALSE,min.chars = 1)
  checkmate::assertCharacter(vp, pattern= "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",len=1,any.missing = FALSE,min.chars = 1)
  checkmate::assertCharacter(vc, pattern= "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",len=1,any.missing = FALSE,min.chars = 1)
  checkmate::assertCharacter(q, pattern= "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",len=1,any.missing = FALSE,min.chars = 1)
  
  if(!missing(peripheral)){
    checkmate::assertIntegerish(peripheral, lower=1,any.missing = FALSE,len = 1)
    
  }
  
  temp  <- rxode2::assertRxUi(model)
  
  #browser()
  mv <- rxode2::rxModelVars(temp)
  
  
  if (!(central %in% mv$state)){
    stop("'",central,"' needs to be in the model")
  }
  if (!(any(grepl("^peripheral",mv$state)))){
    stop("'",peripheralComp," need to be in the model")
  }
  
  #Extract model
  modelNew <- rxode2::modelExtract(temp,endpoint=NA)
  
  #modify ODE for central compartment to delete all elements related to peripheral compartments
  center<- eval(str2lang(paste0("rxode2::modelExtract(temp,d/dt(",central,"),lines=TRUE)")))
  centralLine <- attr(center,"lines")
  rhs <- sub(".*<-\\s*","",center)
  rhs <- gsub("\\s*[+-]?\\s*k(?:[0-9][0-9])\\s*\\*\\s*\\w+\\d*","",rhs)
  
  #Find total number of peripheral compartments in the model
  obj <- c(unlist(modelNew)[which(grepl("\\s*^peripheral",mv$state))])
  totalPeripheral <- length(obj)
  
  if(missing(peripheral)){
    peripheral <- totalPeripheral
    line <- str2lang(paste0("d/dt(",central,") <- ",rhs))
  }
  
  #Modify ini{}
  temp2<- temp$iniDf
  temp3<-temp2$name
  ini1 <- c(paste0("l",vp))
  ini2 <- c(paste0("l",q))
  if (totalPeripheral>1){
    for (i in totalPeripheral:(totalPeripheral-peripheral+1)){
      ini1 <- c(ini1,paste0("l",vp,i+1))
      ini2 <- c(ini2,paste0("l",q,i+1))
      
  }
  }
  temp4 <- temp3[!(temp3 %in% c(ini1,ini2))]
  temp2 <- temp2[temp2$name %in% temp4, ]
  rxode2::ini(temp) <-temp2
  #browser()
  #Locate the ODEs for peripheral compartments to be deleted
  obj=c()
  for (i in totalPeripheral:(totalPeripheral-peripheral+1)){
    obj1 <- eval(str2lang(paste0("rxode2::modelExtract(temp,d/dt(",peripheralComp,i,"), lines = TRUE)")))
    obj2 <- eval(str2lang(paste0("rxode2::modelExtract(temp, k1",i+1,",lines = TRUE)")))
    obj3 <- eval(str2lang(paste0("rxode2::modelExtract(temp, k",i+1,"1,lines = TRUE)")))
    obj4 <- eval(str2lang(paste0("rxode2::modelExtract(temp,",vp,i,",lines = TRUE)")))
    obj5 <- eval(str2lang(paste0("rxode2::modelExtract(temp,",q,i,",lines = TRUE)")))
    
    obj=c(obj,obj1,obj2,obj3,obj4,obj5)
  }
  obj6 <- rxode2::modelExtract(temp,vp,q,lines = TRUE)
  obj=c(obj,obj6)
  
  for (i in obj){
    index <- which(modelNew==i)
    modelNew <- modelNew[-index]
  }
  
  #Insert modified ODE for central compartment into the model and modify model{}
  rxode2::model(temp) <- modelNew
  temp2 <- temp %>%
    rxode2::model(line)
  temp2
  
  
  
  
  
  
  
  
  
  
  
  
  
}

