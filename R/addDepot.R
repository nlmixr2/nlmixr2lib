#' To convert from infusion/intravenous administration to first-order oral
#' absorption
#'
#' @param model The model as a function (or something convertible to an rxUi
#'   object)
#' @param central central compartment name
#' @param depot depot compartment name
#' @param ka absorption rate parameter name
#' @param lag lag parameter name
#' @param lagIni Initial value for the lag time (`NA` to omit)
#' @param fdepotIni Initial value for the depot bioavailability (`NA` to omit)
#' @param absRateIni Initial value for the first order rate
#' @return a model with the depot added
#' @export
#' @examples
#' # most of the examples in the model library already have a depot.
#' # for this example we will remove the depot and then add it back
#' readModelDb("PK_1cmt_des") |>
#'   removeDepot() |>
#'   addDepot()
addDepot <- function(model,
                     central = "central", depot = "depot",
                     ka = "ka", lag = paste0("lag", depot),
                     lagIni=NA, fdepotIni=NA,
                     absRateIni=1.0) {
  model <- rxode2::assertRxUi(model)
  assertCompartmentName(depot)
  assertCompartmentExists(model, central)
  assertVariableName(ka)
  assertVariableName(lag)
  if (!is.na(lagIni)) {
    assertParameterValue(lagIni)
  }
  if (!is.na(fdepotIni)) {
    assertParameterValue(fdepotIni)
  }
  assertParameterValue(kaIni)
  temp <- rxode2::assertRxUi(model)

  mv <- rxode2::rxModelVars(temp)
  if (ka %in% mv$params) {
    stop("'", ka, "' cannot be in the model")
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
    transitODE <- str2lang(paste0("d/dt(transit1) <- ", ka, "*", depot, deparse1(str2lang(transitRhs))))
  }
  # Extract model
  model <- rxode2::modelExtract(temp, endpoint = NA)

  # Extract and modify ODE for central compartment
  center <- eval(str2lang(paste0("rxode2::modelExtract(temp,d/dt(", central, "),lines=TRUE)")))
  rhs <- sub(".*<-\\s*", "", center)
  line <- str2lang(paste0("d/dt(", central, ") <- ", ka, "*", depot, "+", deparse1(str2lang(rhs))))
  lineNew <- str2lang(paste0("d/dt(", central, ") <- ", deparse1(str2lang(rhs))))
  centralLine <- attr(center, "lines")

  # Additional equations to be added to model block
  kaModel <- paste0(ka, " <- exp(l", ka, ")")

  # Lag equations
  if (is.na(lagIni)) {
    lagModel <- lagODE <- NULL
  } else {
    lagModel <- paste0(lag, " <- exp(la", lag, ")")
    lagODE <- paste0("alag(", depot, ") <- ", lag)
  }

  # Bioavailability equations
  if (is.na(fdepotIni)) {
    fdepotModel <- fdepotODE <- NULL
  } else {
    fdepotModel <- paste0("f", depot, " <- exp(lf", depot, ")")
    fdepotODE <- paste0("f(", depot, ") <- f", depot)
  }

  depotODE <- paste0("d/dt(", depot, ") <- -", ka, "*", depot)

  # Modify model block
  if (any(grepl("^transit", mv$state))) {
    rxode2::model(temp) <-
      c(
        kaModel,
        fdepotModel,
        lagModel,
        model[1:(transitLine - 1)],
        depotODE,
        fdepotODE,
        lagODE,
        transitODE,
        model[(transitLine + 1):(centralLine - 1)],
        lineNew,
        model[(centralLine + 1):length(model)]
      )
  } else {
    rxode2::model(temp) <-
      c(
        kaModel,
        fdepotModel,
        lagModel,
        model[1:(centralLine - 1)],
        depotODE,
        fdepotODE,
        lagODE,
        line,
        model[(centralLine + 1):length(model)]
      )
  }

  # Modify ini block
  rateIni <- str2lang(paste0("l", ka, " <-", kaIni))
  lfdepotIni <- str2lang(paste0("lf", depot, " <-", fdepotIni))
  if (is.na(lagIni)) {
    temp <- rxode2::ini(temp, rateIni, append = 0)
    if (!is.na(fdepotIni)) {
      temp <- rxode2::ini(temp, lfdepotIni)
    }
    temp2 <- temp$iniDf
    suppressMessages(temp2$label[temp2$name == paste0("l", ka)] <- paste0("First order absorption rate (", ka, ")"))
    suppressMessages(temp2$label[temp2$name == paste0("lf", depot)] <- "Bioavailability (F)")
  } else {
    lalagIni <- str2lang(paste0("la", lag, " <-", lagIni))
    temp <- temp |>
      rxode2::ini(rateIni, append = 0) |>
      rxode2::ini(lalagIni)
    if (!is.na(fdepotIni)) {
      temp <- rxode2::ini(temp, lfdepotIni)
    }
    temp2 <- temp$iniDf
    suppressMessages(temp2$label[temp2$name == paste0("l", ka)] <- paste0("First order absorption rate (", ka, ")"))
    if (!is.na(fdepotIni)) {
      suppressMessages(temp2$label[temp2$name == paste0("lf", depot)] <- "Bioavailability (F)")
    }
    suppressMessages(temp2$label[temp2$name == paste0("la", lag)] <- paste0("Lag time (", lag, ")"))
  }
  rxode2::ini(temp) <- temp2

  # return
  temp
}

#' To convert from first order oral absorption to IV/Intravenous
#' @param model The model as a function
#' @param central This is a character vector that represents the central compartment in the model
#' @param depot  This is a character vector that represents the depot in the model
#' @return Returns a model with the depot from a first order absorption model removed
#' @export
#' @examples
#' readModelDb("PK_1cmt_des") |>
#'   removeDepot()
removeDepot <- function(model, central = "central", depot = "depot") {
  assertCompartmentName(central)
  assertCompartmentName(depot)

  temp <- rxode2::assertRxUi(model)
  mv <- rxode2::rxModelVars(temp)
  if (!(central %in% mv$state)) {
    stop("'", central, "' needs to be in the model")
  }
  if (!(depot %in% mv$state)) {
    stop("'", depot, "' needs to be in the model")
  }
  if (any(grepl("^transit", mv$state))) {
    transit <- eval(str2lang(paste0("rxode2::modelExtract(temp,d/dt(transit1),lines=TRUE)")))
    transitLine <- attr(transit, "lines")
    transitNew <- str2lang(sub("\\s*ka\\s*\\*\\s*depot", "", transit))
  }


  model <- rxode2::modelExtract(temp, endpoint = NA)
  center <- eval(str2lang(paste0("rxode2::modelExtract(temp,d/dt(", central, "),lines=TRUE)")))
  centralLine <- attr(center, "lines")
  rhs <- str2lang(sub("\\s*ka\\s*\\*\\s*depot", "", center))
  if (any(grepl("^transit", mv$state))) {
    rxode2::model(temp) <- c(model[1:(transitLine - 1)], transitNew, model[(transitLine + 1):(centralLine - 1)], rhs, model[(centralLine + 1):length(model)])
  } else {
    rxode2::model(temp) <- c(model[1:(centralLine - 1)], rhs, model[(centralLine + 1):length(model)])
  }

  ka <- fdepot <- depot <- d <- dt <- f <- NULL
  if ("fdepot" %in% mv$lhs) {
    temp <- rxode2::model(temp, -ka, -fdepot, -f(depot), -d / dt(depot))
  } else {
    temp <- rxode2::model(temp, -ka, -d / dt(depot))
  }
  temp
}
