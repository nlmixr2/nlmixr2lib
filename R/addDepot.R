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
#' Drop additive expression that is a .
#'
#' @param x expression
#' @return expression with dot removed
#' @author Matthew L. Fidler
#' @noRd
#' @examples
#' .dropDotAddExpr(str2lang("x+.+4"))
#' .dropDotAddExpr(str2lang("+.+4"))
#' .dropDotAddExpr(str2lang(".+4"))
#' .dropDotAddExpr(str2lang("4+x+."))
#' .dropDotAddExpr(str2lang("x-.+4"))
#' .dropDotAddExpr(str2lang("-.+4"))
#' .dropDotAddExpr(str2lang("+.-4"))
#' .dropDotAddExpr(str2lang("+4-."))
#' .dropDotAddExpr(str2lang("d/dt(central) <- . - kel * central"))
#' .dropDotAddExpr(list(str2lang("d/dt(central) <- . - kel * central")))
.dropDotAddExpr <- function(x) {
  if (is.list(x)) {
    return(lapply(seq_along(x),
                  function(i) {
                    .dropDotAddExpr(x[[i]])
                  }))
  }
  if (is.call(x)) {
    if (length(x) == 2) {
      if (length(x[[1]]) != 1) {
        .x1 <- .dropDotAddExpr(x[[1]])
      } else {
        .x1 <- x[[1]]
      }
      if (length(x[[2]]) != 1) {
        .x2 <- .dropDotAddExpr(x[[2]])
      } else {
        .x2 <- x[[2]]
      }
      if (identical(.x1, quote(`+`)) ||
            identical(.x1, quote(`-`))) {
        return(.x2)
      }
    }
    if (length(x) == 3) {
      if (length(x[[1]]) != 1) {
        .x1 <- .dropDotAddExpr(x[[1]])
      } else {
        .x1 <- x[[1]]
      }
      if (length(x[[2]]) != 1) {
        .x2 <- .dropDotAddExpr(x[[2]])
      } else {
        .x2 <- x[[2]]
      }
      if (length(x[[3]]) != 1) {
        .x3 <- .dropDotAddExpr(x[[3]])
      } else {
        .x3 <- x[[3]]
      }
      if (identical(.x1, quote(`=`)) ||
            identical(.x1, quote(`<-`))) {
        return(as.call(list(.x1, .x2, .x3)))
      }
      if (identical(.x1, quote(`-`))) {
        if (identical(.x2, quote(`.`))) {
          return(str2lang(paste0("-", deparse1(.x3))))
        }
        if (identical(.x3, quote(`.`))) {
          return(.x2)
        }
      }
      if (identical(.x1, quote(`+`))) {
        if (identical(.x2, quote(`.`))) {
          return(.x3)
        }
        if (identical(.x3, quote(`.`))) {
          return(.x2)
        }
      }
    }
    as.call(lapply(x, .dropDotAddExpr))
  } else {
    x
  }
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
