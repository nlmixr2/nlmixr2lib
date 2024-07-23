#' To add additional compartments to the model
#' @param model The model as a function
#' @param numPeripheral number of peripheral compartments to be added
#'   to the model
#' @param central a character vector representing the central
#'   compartment
#' @param peripheralComp A character vector representing the prefix of
#'   peripheral compartments
#' @param vp parameter representing the peripheral volume of the first
#'   (central) compartment and the prefix of the other compartments
#' @param vc parameter representing the central volume of the first
#'   (central) compartment and the prefix of the other compartment's
#'   volume
#' @param q inter-compartmental clearance parameter or prefix
#'   (depending on the model)
#' @return A rxode2 model function with an additional compartment added
#' @export
#' @examples
#' readModelDb("PK_1cmt_des") |>
#'   addComp(1)
addComp <- function(model, numPeripheral, central = "central",
                    peripheralComp = "peripheral", vp = "vp", vc = "vc", q = "q") {
  assertCompartmentName(central)
  assertCompartmentName(peripheralComp)
  assertVariableName(vp)
  assertVariableName(vc)
  assertVariableName(q)

  checkmate::assertIntegerish(numPeripheral, lower = 1, upper = 4)
  temp <- rxode2::assertRxUi(model)
  mv <- rxode2::rxModelVars(temp)
  if (!(central %in% mv$state)) {
    stop("'", central, "' needs to be in the model")
  }
  if ((any(grepl("^peripheral", mv$state)))) {
    model <- removeComp(model)
    temp <- rxode2::assertRxUi(model)
  }

  # Find the location of central compartment
  model <- rxode2::modelExtract(temp, endpoint = NA)
  center <- eval(str2lang(paste0("rxode2::modelExtract(temp,d/dt(", central, "),lines=TRUE)")))
  centralLine <- attr(center, "lines")

  # Generate additional rate constants,clearances,volumes for additional peripheral compartments
  periVP <- list(paste0(vp, " <- exp(l", vp, ")"))
  periq <- list(paste0(q, " <- exp(l", q, ")"))
  rates <- list("k12 <- q/vc", "k21 <- q/vp")
  if (numPeripheral > 1) {
    for (i in 2:numPeripheral) {
      periVP <- c(periVP, paste0(vp, i, "<- exp(l", vp, i, ")"))
      periq <- c(periq, paste0(q, i, "<- exp(l", q, i, ")"))
      eq1 <- paste0("k1", i + 1, " <- ", q, i - 1, "/", vc)
      eq2 <- paste0("k", i + 1, "1 <- ", q, i - 1, "/", vp, i)
      rates <- c(rates, eq1, eq2)
    }
  }

  # Generate new ODEs for additional peripheral compartments
  periODE <- list()
  for (i in 1:numPeripheral) {
    periODE <- c(periODE, str2lang(paste0("d/dt(", peripheralComp, i, ") <- k1", i + 1, "*", central, "- k", i + 1, "1 *", peripheralComp, i)))
  }

  # Extract ODE of central compartment without peripheral component
  # currently code k12*abcd#
  # fix to allow central/periph# * k12 and vice versa
  rhs <- gsub("\\s*[+-]?\\s*k(?:[0-9][0-9])\\s*\\*\\s*\\w+\\d*", "", center)

  # Modify ODE of central compartment with the number of peripheral components specified by user
  newElements <- ""
  for (i in 1:numPeripheral) {
    toElement <- paste0(" + k", i + 1, "1 * ", peripheralComp, i)
    fromElement <- paste0(" - k1", i + 1, " * ", central)
    newElements <- paste0(newElements, toElement, fromElement)
  }
  center <- paste0(rhs, newElements)

  vcLine <- eval(str2lang(paste0("rxode2::modelExtract(temp,", vc, ",lines=TRUE)")))
  vcLine <- attr(vcLine, "lines")
  rxode2::model(temp) <- c(model[1:vcLine], periVP, periq, rates, center, periODE, model[(length(model) - 1):length(model)])

  temp

  # Modify ini{}

  ini <- list("lvp", "lq")
  equationIni <- rep(0.05, 2)
  if (i > 1) {
    for (i in 2:numPeripheral) {
      ini1 <- paste0("lvp", i)
      ini2 <- paste0("lq", i)
      ini <- c(ini, ini1, ini2)
      equationIni <- c(equationIni, 0.05, 0.05)
    }
  }
  names(equationIni) <- ini
  temp2 <- rxode2::ini(temp, equationIni)
  temp3 <- temp2$iniDf
  rxode2::ini(temp) <- temp3
  temp
}

#' To remove peripheral compartments from the model
#' @param model The model as a function
#' @param peripheral The number of peripheral compartments to remove
#' @inheritParams addComp
#' @return rxode2 model function/ui with a compartment removed
#' @export
#' @examples
#' library(rxode2)
#' readModelDb("PK_2cmt_des") |> removeComp(1)
removeComp <- function(model, peripheral, central = "central", peripheralComp = "peripheral", vp = "vp", vc = "vc", q = "q") {
  assertCompartmentName(central)
  assertCompartmentName(peripheralComp)
  assertVariableName(vp)
  assertVariableName(vc)
  assertVariableName(q)

  if (!missing(peripheral)) {
    checkmate::assertIntegerish(peripheral, lower = 1, any.missing = FALSE, len = 1)
  }

  temp <- rxode2::assertRxUi(model)

  mv <- rxode2::rxModelVars(temp)


  if (!(central %in% mv$state)) {
    stop("'", central, "' needs to be in the model")
  }
  if (!(any(grepl("^peripheral", mv$state)))) {
    stop("'", peripheralComp, " need to be in the model")
  }

  # Extract model
  modelNew <- rxode2::modelExtract(temp, endpoint = NA)

  # modify ODE for central compartment to delete all elements related to peripheral compartments
  center <- eval(str2lang(paste0("rxode2::modelExtract(temp,d/dt(", central, "),lines=TRUE)")))
  centralLine <- attr(center, "lines")
  rhs <- sub(".*<-\\s*", "", center)
  rhs <- gsub("\\s*[+-]?\\s*k(?:[0-9][0-9])\\s*\\*\\s*\\w+\\d*", "", rhs)

  # Find total number of peripheral compartments in the model
  obj <- c(unlist(modelNew)[which(grepl("\\s*^peripheral", mv$state))])
  totalPeripheral <- length(obj)

  if (missing(peripheral)) {
    peripheral <- totalPeripheral
  }
  line <- str2lang(paste0("d/dt(", central, ") <- ", rhs))

  # Modify ini{}
  temp2 <- temp$iniDf
  temp3 <- temp2$name
  ini1 <- c(paste0("l", vp))
  ini2 <- c(paste0("l", q))
  if (totalPeripheral > 1) {
    for (i in totalPeripheral:(totalPeripheral - peripheral + 1)) {
      ini1 <- c(ini1, paste0("l", vp, i + 1))
      ini2 <- c(ini2, paste0("l", q, i + 1))
    }
  }
  temp4 <- temp3[!(temp3 %in% c(ini1, ini2))]
  temp2 <- temp2[temp2$name %in% temp4, ]
  rxode2::ini(temp) <- temp2

  # Locate the ODEs for peripheral compartments to be deleted
  obj <- NULL
  for (i in totalPeripheral:(totalPeripheral - peripheral + 1)) {
    obj1 <- eval(str2lang(paste0("rxode2::modelExtract(temp,d/dt(", peripheralComp, i, "), lines = TRUE)")))
    obj2 <- eval(str2lang(paste0("rxode2::modelExtract(temp, k1", i + 1, ",lines = TRUE)")))
    obj3 <- eval(str2lang(paste0("rxode2::modelExtract(temp, k", i + 1, "1,lines = TRUE)")))
    obj4 <- eval(str2lang(paste0("rxode2::modelExtract(temp,", vp, i, ",lines = TRUE)")))
    obj5 <- eval(str2lang(paste0("rxode2::modelExtract(temp,", q, i, ",lines = TRUE)")))

    obj <- c(obj, obj1, obj2, obj3, obj4, obj5)
  }
  obj6 <- rxode2::modelExtract(temp, vp, q, lines = TRUE)
  obj <- c(obj, obj6)

  for (i in obj) {
    index <- which(modelNew == i)
    modelNew <- modelNew[-index]
  }

  # Insert modified ODE for central compartment into the model and modify model{}
  rxode2::model(temp) <- modelNew

  temp <- rxode2::model(temp, line)

  temp
}
