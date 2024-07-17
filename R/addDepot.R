#' To convert from Infusion/intravenous administration to first order
#' oral absorption
#' @param model The model as a function
#' @param central central compartment name
#' @param depot depot name
#' @param absRate absorption rate
#' @param lag A boolean representing if you are going to add a lag to
#'   this compartment
#' @param tlag a character vector representing the lag time
#' @param fdepot boolean that determines if the bioavailability of the
#'   depot compartment is included.
#' @param lagIni Initial value for the lag time
#' @param fdepotIni Initial value for the depot
#' @param absRateIni Initial value for the first order rate
#' @return a model with the depot added
#' @export
#' @examples
#' # most of the examples in the model library already have a depot.
#' # for this example we will remove the depot and then add it back
#' readModelDb("PK_1cmt_des") |>
#'   removeDepot() |>
#'   addDepot()
addDepot <- function(model, central = "central", depot = "depot", absRate = "ka", lag = FALSE, tlag = "lagD",
                     fdepot=FALSE,
                     lagIni=0.0, fdepotIni=1.0,
                     absRateIni=1.0) {
  assertCompartmentName(depot)
  checkmate::assertCharacter(absRate, pattern = "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$", len = 1, any.missing = FALSE, min.chars = 1)
  checkmate::assertLogical(lag, len = 1, any.missing = FALSE)
  checkmate::assertLogical(fdepot, len = 1, any.missing = FALSE)
  checkmate::assertCharacter(tlag, pattern = "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$", len = 1, any.missing = FALSE, min.chars = 1)
  checkmate::assertNumeric(lagIni, len=1, any.missing=FALSE, finite = TRUE)
  checkmate::assertNumeric(fdepotIni, len=1, any.missing=FALSE, finite = TRUE)
  checkmate::assertNumeric(absRateIni, len=1, any.missing=FALSE, finite = TRUE)
  temp <- rxode2::assertRxUi(model)
  mv <- rxode2::rxModelVars(temp)
  if (absRate %in% mv$params) {
    stop("'", absRate, "' cannot be in the model")
  }
  if (!(central %in% mv$state)) {
    stop("'", central, "' needs to be in the model")
  }
  if (depot %in% mv$state) {
    stop("'", depot, "' cannot be in the model")
  }
  if (any(grepl("^transit", mv$state))) {
    transit <- eval(str2lang(paste0("rxode2::modelExtract(temp,d/dt(transit1),lines=TRUE)")))
    transitLine <- attr(transit, "lines")
    transitRhs <- sub(".*<-\\s*", "", transit)
    transitODE <- str2lang(paste0("d/dt(transit1) <- ", absRate, "*", depot, deparse1(str2lang(transitRhs))))
  }
  # Extract model
  model <- rxode2::modelExtract(temp, endpoint = NA)

  # Extract and modify ODE for central compartment
  center <- eval(str2lang(paste0("rxode2::modelExtract(temp,d/dt(", central, "),lines=TRUE)")))
  rhs <- sub(".*<-\\s*", "", center)
  line <- str2lang(paste0("d/dt(", central, ") <- ", absRate, "*", depot, "+", deparse1(str2lang(rhs))))
  lineNew <- str2lang(paste0("d/dt(", central, ") <- ", deparse1(str2lang(rhs))))
  centralLine <- attr(center, "lines")

  # Additional equations to be added to model block
  absrateModel <- paste0(absRate, " <- exp(l", absRate, ")")
  lagModel <- paste0(tlag, " <- exp(la", tlag, ")")

  fdepotModel <- paste0("f", depot, " <- exp(lf", depot, ")")
  fdepotODE <- paste0("f(", depot, ") <- f", depot)
  if (!fdepot) {
    fdepotModel <- fdepotODE <- NULL
  }

  depotODE <- paste0("d/dt(", depot, ") <- -", absRate, "*", depot)
  lagODE <- paste0("alag(", depot, ") <- ", tlag)

  # Modify model block
  if (lag == FALSE) {
    rxode2::model(temp) <- c(absrateModel, fdepotModel, model[1:(centralLine - 1)], depotODE, fdepotODE, line, model[(centralLine + 1):length(model)])
  } else {
    rxode2::model(temp) <- c(absrateModel, fdepotModel, lagModel, model[1:(centralLine - 1)], depotODE, fdepotODE, lagODE, line, model[(centralLine + 1):length(model)])
  }

  if (any(grepl("^transit", mv$state))) {
    rxode2::model(temp) <- c(absrateModel, fdepotModel, lagModel, model[1:(transitLine - 1)], depotODE, fdepotODE, lagODE, transitODE, model[(transitLine + 1):(centralLine - 1)], lineNew, model[(centralLine + 1):length(model)])
  }

  # Modify ini block
  rateIni <- str2lang(paste0("l", absRate, " <-", absRateIni))
  lfdepotIni <- str2lang(paste0("lf", depot, " <-", fdepotIni))
  lalagIni <- str2lang(paste0("la", tlag, " <-", lagIni))
  if (lag == FALSE) {
    temp <- rxode2::ini(temp, rateIni, append = 0)
    if (fdepot) {
      temp <- rxode2::ini(temp, lfdepotIni)
    }
    temp2 <- temp$iniDf
    suppressMessages(temp2$label[temp2$name == paste0("l", absRate)] <- paste0("First order absorption rate (", absRate, ")"))
    suppressMessages(temp2$label[temp2$name == paste0("lf", depot)] <- "Bioavailability (F)")
  } else {
    temp <- temp |>
      rxode2::ini(rateIni, append = 0) |>
      rxode2::ini(lalagIni)
    if (fdepot) {
      temp <- rxode2::ini(temp, lfdepotIni)
    }
    temp2 <- temp$iniDf
    suppressMessages(temp2$label[temp2$name == paste0("l", absRate)] <- paste0("First order absorption rate (", absRate, ")"))
    if (fdepot) {
      suppressMessages(temp2$label[temp2$name == paste0("lf", depot)] <- "Bioavailability (F)")
    }
    suppressMessages(temp2$label[temp2$name == paste0("la", tlag)] <- paste0("Lag time (", tlag, ")"))
  }
  rxode2::ini(temp) <- temp2

  # return
  temp
}
