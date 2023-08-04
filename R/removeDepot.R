#' To convert from first order oral absorption to IV/Intravenous
#' @param model The model as a function
#' @param equation The modified ODE for central compartment in the model
#' @examples
#' library(rxode2)
#' 
removeDepot <- function(model,central="central",depot="depot"){
  #browser()
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
  #browser()
  
  model <- rxode2::modelExtract(temp,endpoint=NA)
  center<- eval(str2lang(paste0("rxode2::modelExtract(temp,d/dt(",central,"),lines=TRUE)")))
  centralLine <- attr(center,"lines")
  rhs <- str2lang(sub("\\s*ka\\s*\\*\\s*depot", "",center))
  rxode2::model(temp) <- c(model[1:(centralLine-1)],rhs, model[(centralLine+1):length(model)])
  if ("fdepot" %in% mv$lhs){
    temp <- temp%>%model(-ka,-fdepot,-f(depot),-d/dt(depot))
  }else{
    temp <- temp%>%model(-ka,-d/dt(depot))
  }
  temp
}