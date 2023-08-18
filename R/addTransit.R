#' To add transit compartments to the model
#' @param model The model as a function
#' @param transit the number of transit compartments to be added
#' @param transitComp the transit compartment prefix
#' @param ktr the parameter name for the transit compartment rate
#' @inheritParams addComp
#' @export
#' @examples
#' readModelDb("PK_1cmt_des") |>
#'   addTransit(3)
addTransit <- function(model,transit,central="central",depot="depot",transitComp ="transit",ktr="ktr"){
  checkmate::assertCharacter(central, pattern= "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",len=1,any.missing = FALSE,min.chars = 1)
  checkmate::assertCharacter(depot,pattern= "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",len=1,any.missing = FALSE,min.chars = 1)
  checkmate::assertCharacter(transitComp,pattern= "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",len=1,any.missing = FALSE,min.chars = 1)
  checkmate::assertNumeric(transit, lower=1, upper=10)
  checkmate::assertCharacter(ktr,pattern= "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",len=1,any.missing = FALSE,min.chars = 1)
  temp  <- rxode2::assertRxUi(model)
  mv <- rxode2::rxModelVars(temp)
  if (!(central %in% mv$state)){
    stop("'",central,"' needs to be in the model")
  }
  if(!(depot %in% mv$state)){
    stop("'",depot,"' needs to be in the model")
  }
  if ((any(grepl("^transit",mv$state)))){
    model <- removeTransit(model)
    temp <-rxode2::assertRxUi(model)
  }
  
  #Extract model and central ODE
  model <- rxode2::modelExtract(temp,endpoint=NA)
  center<- eval(str2lang(paste0("rxode2::modelExtract(temp,d/dt(",central,"),lines=TRUE)")))
  centralLine <- attr(center,"lines")
  
  #Modify ODE for central compartment
  rhs <- sub(".*<-\\s*","",center)
  rhs <- sub("\\s*ka\\s*\\*\\s*depot", "",rhs)
  line <- str2lang(paste0("d/dt(",central,") <- ",ktr,"*",transitComp,transit,deparse(str2lang(rhs))))
  
  
  #ODEs for transit compartments
  equation1 <- list(str2lang(paste0(ktr," <- exp(l",ktr,")")))
  equation2 <- list(str2lang(paste0("d/dt(",transitComp,"1) <- ka*",depot," - ",ktr,"*",transitComp,"1")))
  if(transit>1){
    for (i in 2:transit){
      equation2 <- c(equation2,str2lang(paste0("d/dt(",transitComp,i,")<- ",ktr,"*",transitComp,i-1,"-",ktr,"*",transitComp,i)))
    }
  }
  #modify ini block
  equationIni <- setNames(0.05, paste0("l",ktr))
  
  #modify model block 
  rxode2::model(temp) <- c(equation1, model[1:(centralLine-1)],equation2, line,model[(centralLine+1):length(model)])
  temp2 <- rxode2::ini(temp, equationIni)
  temp3 <- temp2$iniDf
  w <- which(temp3$name %in% names(equationIni))
  
  #Update labels in ini block
  for (i in 1:transit){
    suppressMessages(temp3$label[temp3$name==paste0("l",ktr,i)] <- paste0("First order transition rate (",ktr,i,")"))
  }
  rxode2::ini(temp2) <- temp3
  temp2
}