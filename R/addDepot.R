#' To convert from infusion/intravenous administration to first-order oral
#' absorption
#'
#' @param ui The model as a function (or something convertible to an rxUi
#'   object)
#' @param central central compartment name
#' @param depot depot compartment name
#' @param ka absorption rate parameter name
#' @return a model with the depot added
#' @export
#' @examples
#' # most of the examples in the model library already have a depot
#' # the PK_2cmt_no_depot is an exception
#' readModelDb("PK_2cmt_no_depot")  |> addDepot()
addDepot <- function(ui,
                     central = "central", depot = "depot",
                     ka="ka") {
  .ui <- rxode2::assertRxUi(ui)
  assertCompartmentName(depot)
  assertCompartmentExists(.ui, central)
  assertVariableName(ka)
  .mv <- rxode2::rxModelVars(.ui)
  # Get the central ODE and add depot to it
  .modelLines <- .ui$lstExpr
  .w <- .whichDdt(.modelLines, central)
  .tmp <- .extractModelLinesAtW(.modelLines, .w)
  .tmp$w <- str2lang(paste0(deparse1(.tmp$w), "+", ka, "*", depot))
  .modelLines <- c(list(str2lang(paste0(ka, "<- exp(l", ka, ")"))),
                   .tmp$pre,
                   list(str2lang(paste0("d/dt(", depot, ") <- -", ka, "*", depot))),
                   list(.tmp$w),
                   .tmp$post)

  .tmp <- .getEtaThetaTheta1(.ui)
  .iniDf <- .tmp$iniDf
  .theta <- .tmp$theta
  .theta1 <- .tmp$theta1
  .eta <- .tmp$eta
  if (length(.iniDf$name) == 0L)  {
    .ntheta <- 0
  } else {
    .ntheta <- max(.iniDf$ntheta)
  }

  .thetaka <- .get1theta(ka, .theta1, .ntheta,
                         label=paste0("First order absorption rate (", ka, ")"))
  .ntheta <- .ntheta + 1

  .ui <- rxode2::rxUiDecompress(.ui)
  .ui$iniDf <- rbind(.theta,
                     .thetaka,
                     .eta)
  if (exists("description", envir=.ui$meta)) {
    rm("description", envir=.ui$meta)
  }
  rxode2::model(.ui) <- .modelLines
  rxode2::rxUiCompress(.ui)
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
