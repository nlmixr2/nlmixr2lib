#' To convert from Infusion/intravenous administration to first order oral absorption
#' @param model The model as a function
#' @param central central compartment name
#' @param depot depot name
#' @param absRate absorption rate
#' @export
#' @examples
#' # most of the examples in the model library already have a depot.
#' # for this example we will remove the depot and then add it back 
#' readModelDb("PK_1cmt_des") |>
#'   removeDepot() |>
#'   addDepot()
addDepot <- function(model,central="central",depot="depot",absRate="ka",lag=FALSE,tlag="lagD") {
  #browser()
  checkmate::assertCharacter(central,pattern= "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$", len=1,any.missing = FALSE,min.chars = 1)
  checkmate::assertCharacter(depot,pattern= "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$", len=1,any.missing = FALSE,min.chars = 1)
  checkmate::assertCharacter(absRate,pattern= "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",len=1,any.missing = FALSE,min.chars = 1)
  checkmate::assertLogical(lag,len=1,any.missing = FALSE)
  temp  <- rxode2::assertRxUi(model)
  mv <- rxode2::rxModelVars(temp)
  if (absRate %in% mv$params) {
    stop("'",absRate,"' cannot be in the model")
  }
  if (!(central %in% mv$state)) {
    stop("'",central,"' needs to be in the model")
  }
  if (depot %in% mv$state) {
    stop("'",depot,"' cannot be in the model")
  }
  if (any(grepl("^transit",mv$state))) {
    transit<- eval(str2lang(paste0("rxode2::modelExtract(temp,d/dt(transit1),lines=TRUE)")))
    transitLine <- attr(transit,"lines")
    transitRhs <- sub(".*<-\\s*","",transit)
    transitODE <- str2lang(paste0("d/dt(transit1) <- ",absRate,"*",depot,deparse1(str2lang(transitRhs))))
  }
  #Extract model
  model <- rxode2::modelExtract(temp,endpoint=NA)
  
  #Extract and modify ODE for central compartment
  center<- eval(str2lang(paste0("rxode2::modelExtract(temp,d/dt(",central,"),lines=TRUE)")))
  rhs <- sub(".*<-\\s*","",center)
  line <- str2lang(paste0("d/dt(",central,") <- ",absRate,"*",depot,"+",deparse1(str2lang(rhs))))
  lineNew <- str2lang(paste0("d/dt(",central,") <- ",deparse1(str2lang(rhs))))
  centralLine <- attr(center, "lines")
  
  #Additional equations to be added to model block
  absrateModel <- paste0(absRate," <- exp(l",absRate,")")
  fdepotModel <- paste0("f",depot," <- exp(lf",depot,")")
  lagModel <- paste0(tlag," <- exp(la",tlag,")")
  fdepotODE <- paste0("f(",depot,") <- f",depot)
  depotODE <-  paste0("d/dt(",depot,") <- -",absRate,"*",depot)
  lagODE <- paste0("alag(",depot,") <- ",tlag)
  
  #Modify model block
  if (lag==FALSE){
    rxode2::model(temp) <- c(absrateModel,fdepotModel,model[1:(centralLine-1)],depotODE, fdepotODE, line, model[(centralLine+1):length(model)])
  }else{
    rxode2::model(temp) <- c(absrateModel,fdepotModel,lagModel,model[1:(centralLine-1)],depotODE, fdepotODE,lagODE, line, model[(centralLine+1):length(model)])
  }
  
  if (any(grepl("^transit",mv$state))){
    rxode2::model(temp) <- c(absrateModel,fdepotModel,lagModel,model[1:(transitLine-1)],depotODE, fdepotODE,lagODE,transitODE,model[(transitLine+1):(centralLine-1)], lineNew, model[(centralLine+1):length(model)])
  }
  
  #Modify ini block
  rateIni <- str2lang(paste0("l",absRate," <-0.02"))
  lfdepotIni <- str2lang(paste0("lf",depot," <-0.04"))
  lalagIni <- str2lang(paste0("la",tlag," <-0.09"))
  if (lag==FALSE) {
    temp <- rxode2::ini(temp, rateIni,append=0) |>
      rxode2::ini(lfdepotIni)
    temp2 <- temp$iniDf
    suppressMessages(temp2$label[temp2$name==paste0("l",absRate)] <- paste0("First order absorption rate (",absRate,")"))
    suppressMessages(temp2$label[temp2$name==paste0("lf",depot)] <- "Bioavailability (F)")
  } else {
    temp <- temp |>
      rxode2::ini(rateIni,append=0) |>
      rxode2::ini(lfdepotIni) |>
      rxode2::ini(lalagIni) 
    temp2 <- temp$iniDf
    suppressMessages(temp2$label[temp2$name==paste0("l",absRate)] <- paste0("First order absorption rate (",absRate,")"))
    suppressMessages(temp2$label[temp2$name==paste0("lf",depot)] <- "Bioavailability (F)")
    suppressMessages(temp2$label[temp2$name==paste0("la",tlag)] <- paste0("Lag time (",tlag,")"))
  }
  rxode2::ini(temp) <- temp2
  
  #return
  temp
}




