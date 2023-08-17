#' To convert from first order oral absorption to IV/Intravenous
#' @param model The model as a function
#' @param equation The modified ODE for central compartment in the model
#' @export
#' @examples
removeDepot <- function(model,central="central",depot="depot"){
  checkmate::assertCharacter(central,len=1,any.missing = FALSE,min.chars = 1)
  checkmate::assertCharacter(depot,len=1,any.missing = FALSE,min.chars = 1)
  temp  <- rxode2::assertRxUi(model)
  mv <- rxode2::rxModelVars(temp)
  if (!(central %in% mv$state)){
    stop("'",central,"' needs to be in the model")
  }
  if (!(depot %in% mv$state)){
    stop("'",depot,"' needs to be in the model")
  }
  if (any(grepl("^transit",mv$state))){
    transit<- eval(str2lang(paste0("rxode2::modelExtract(temp,d/dt(transit1),lines=TRUE)")))
    transitLine <- attr(transit,"lines")
    transitNew <- str2lang(sub("\\s*ka\\s*\\*\\s*depot", "",transit))
  }

  
  model <- rxode2::modelExtract(temp,endpoint=NA)
  center<- eval(str2lang(paste0("rxode2::modelExtract(temp,d/dt(",central,"),lines=TRUE)")))
  centralLine <- attr(center,"lines")
  rhs <- str2lang(sub("\\s*ka\\s*\\*\\s*depot", "",center))
  if (any(grepl("^transit",mv$state))){
    rxode2::model(temp) <- c(model[1:(transitLine-1)],transitNew,model[(transitLine+1):(centralLine-1)], rhs, model[(centralLine+1):length(model)])
  }else{
    rxode2::model(temp) <- c(model[1:(centralLine-1)],rhs, model[(centralLine+1):length(model)]) 
  }

  ka <- fdepot <- depot <- NULL
  if ("fdepot" %in% mv$lhs){
    temp <- rxode2::model(temp, -ka,-fdepot,-f(depot),-d/dt(depot))
  }else{
    temp <- rxode2::model(temp, -ka,-d/dt(depot))
  }
  temp
}
