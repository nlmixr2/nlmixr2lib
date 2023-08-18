#' To add additional compartments to the model
#' @param model The model as a function
#' @param numPeripheral number of peripheral compartments to be added
#'   to the model
#' @param central a character vector representing the central
#'   compartment
#' @param depot a character vector representing the depot compartment
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
addComp <- function(model, numPeripheral, central = "central", depot = "depot", peripheralComp = "peripheral", vp = "vp", vc = "vc", q = "q") {
  checkmate::assertCharacter(central, len = 1, any.missing = FALSE, min.chars = 1)
  checkmate::assertCharacter(depot, len = 1, any.missing = FALSE, min.chars = 1)
  checkmate::assertCharacter(peripheralComp, pattern = "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$", len = 1, any.missing = FALSE, min.chars = 1)
  checkmate::assertCharacter(vp, pattern = "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$", len = 1, any.missing = FALSE, min.chars = 1)
  checkmate::assertCharacter(vc, pattern = "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$", len = 1, any.missing = FALSE, min.chars = 1)
  checkmate::assertCharacter(q, pattern = "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$", len = 1, any.missing = FALSE, min.chars = 1)

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
