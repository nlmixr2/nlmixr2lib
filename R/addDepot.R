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

#' @inheritParams addDepot
#' @return Returns a model with the depot from a first order absorption model removed
#' @export
#' @examples
#' readModelDb("PK_1cmt_des") |> removeDepot()
removeDepot <- function(ui, central = "central", depot = "depot",
                        ka="ka") {
  .ui <- rxode2::assertRxUi(ui)
  assertCompartmentExists(.ui, central)
  assertCompartmentExists(.ui, depot)
  assertVariableName(ka)
  .modelLines <- .rmDdt(.ui$lstExpr, depot)
  .w <- .whichDdt(.modelLines, central)
  .tmp <- .extractModelLinesAtW(.modelLines, .w)
  .tmp$w <- .dropDotAddExpr(.replaceMult(.tmp$w, ka, depot, "."))
  .modelLines <- c(.tmp$pre,
                   .tmp$w,
                   .tmp$post)
  .tmp <- .getEtaThetaTheta1(.ui)
  .iniDf <- .tmp$iniDf
  .eta <- .tmp$eta
  .theta <- .tmp$theta
  .theta1 <- .tmp$theta1
  .theta <- .dropTheta(.theta, ka)
  .eta <- .dropEta(.eta, ka)
  .tmp <- .dropLines(.ui, .modelLines, .theta, .eta, ka)
  .modelLines <- .tmp$modelLines
  .theta <- .tmp$theta
  .eta <- .tmp$eta
  .ui <- rxode2::rxUiDecompress(.ui)
  .ui$iniDf <- rbind(.theta,
                     .eta)
  if (exists("description", envir=.ui$meta)) {
    rm("description", envir=.ui$meta)
  }
  rxode2::model(.ui) <- .modelLines
  rxode2::rxUiCompress(.ui)
}
