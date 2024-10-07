#' To add transit compartments to the model
#'
#' @param ui The model as a function
#' @param ntransit the number of transit compartments to be added
#' @param transit the transit compartment prefix
#' @param ktr the parameter name for the transit compartment rate
#' @inheritParams addDepot
#' @family absorption
#' @return a model with transit compartment added
#'
#' This matches
#'
#' `dose->a0->a1->abs cmt->central`
#'
#' But `a0` is depot so dosing records labeled depot do not need to be
#' changed
#'
#' The abs cmt becomes the last "transit" compartment
#'
#' This is simply for convenience
#'
#' @export
#' @examples
#' readModelDb("PK_1cmt_des") |> addTransit(3)
addTransit <- function(ui, ntransit, central = "central",
                       depot = "depot",
                       transit = "transit",
                       ktr = "ktr",
                       ka="ka") {
  checkmate::assertIntegerish(ntransit, lower = 1)
  rxode2::assertCompartmentName(transit)
  .ui <- rxode2::assertRxUi(ui)
  rxode2::assertCompartmentExists(.ui, central)
  .mv <- rxode2::rxModelVars(.ui)
  if (!rxode2::testCompartmentExists(.ui, depot)) {
    .ui <- addDepot(.ui, central=central, depot=depot, ka=ka)
    .mv <- rxode2::rxModelVars(.ui)
    warning("'", depot, "' added to model for transit model", call.=FALSE)
  } else if (rxode2::testCompartmentExists(.ui, paste0(transit, "1"))) {
    .ui <- removeTransit(ui,
                         central = central,
                         depot = depot, transit=transit,
                         ktr = ktr,
                         ka=ka)
  }
  rxode2::assertCompartmentExists(.ui, depot)

  # Extract model and central ODE
  .tmp <- .getEtaThetaTheta1(.ui)
  .iniDf <- .tmp$iniDf
  .theta <- .tmp$theta
  .theta1 <- .tmp$theta1
  .eta <- .tmp$eta
  .modelLines <- .ui$lstExpr
  # Get the central ODE and modify the depot expression to a transit
  # expression
  .w <- .whichDdt(.modelLines, central)
  .tmp <- .extractModelLinesAtW(.modelLines, .w)
  .pre <- .tmp$pre
  .central <- .replaceMult(.tmp$w,
                           v1=depot, v2=ka,
                           ret=paste0(ka, "*", transit, ntransit))
  .post <- .tmp$post
  .v <- seq_len(ntransit)
  # ODEs for the transit compartment (except the one from the depot)
  .transMid <- lapply(.v,
                      function(i) {
                        if (i == 1) {
                          str2lang(paste0("d/dt(", transit, i, ")<- ",
                                          ktr, "*", depot, "-",
                                          ifelse(ntransit == 1, ka, ktr),
                                          "*", transit, i))
                        } else if (i == ntransit) {
                          str2lang(paste0("d/dt(", transit, i, ")<- ",
                                          ktr, "*", transit, i - 1, "-", ka,
                                          "*", transit, i))
                        } else {
                          str2lang(paste0("d/dt(", transit, i, ")<- ",
                                          ktr, "*", transit, i - 1, "-", ktr,
                                          "*", transit, i))
                        }
                      })
  # combine the lines for now
  .modelLines <- c(.pre,
                   .transMid,
                   .central,
                   .post)

  # Now extract the depot and split the model based on the depot cmt
  .w <- .whichDdt(.modelLines, depot)
  .tmp <- .extractModelLinesAtW(.modelLines, .w)
  .modelLines <- c(list(str2lang(paste0(ktr, " <- exp(l", ktr, ")"))),
                   .tmp$pre,
                   .replaceMult(.tmp$w,
                                v1=ka, v2=depot,
                                ret=paste0(ktr, "*", depot)),
                   .tmp$post)
  if (length(.theta$name) == 0L) {
    .ntheta <- 0
  } else {
    .ntheta <- max(.theta$ntheta)
  }
  .thetaktr <- .get1theta(ktr, .theta1, .ntheta,
                          label=paste0("First order transition rate (", ktr, ")"))
  .ntheta <- .ntheta + 1

  .ui <- rxode2::rxUiDecompress(.ui)
  .ui$iniDf <- rbind(.theta,
                     .thetaktr,
                     .eta)
  if (exists("description", envir=.ui$meta)) {
    rm("description", envir=.ui$meta)
  }

  # modify model block
  rxode2::model(.ui) <- .modelLines
  .ui
}

#' To remove transit compartments from the model
#' @param ui The model as a function
#' @param transit The number of transit compartments to remove
#' @inheritParams addDepot
#' @inheritParams addTransit
#' @return rxode2 model with transit compartment removed
#' @family absorption
#' @export
#' @examples
#'
#' # In this example the transit is added and then a few are removed
#'
#' readModelDb("PK_1cmt_des") |> addTransit(4) |> removeTransit(3)
#'
#' readModelDb("PK_1cmt_des") |> addTransit(4) |> removeTransit()
removeTransit <- function(ui, ntransit, central = "central",
                          depot = "depot", transit = "transit",
                          ktr = "ktr",
                          ka="ka") {
  if (!missing(ntransit)) {
    checkmate::assertIntegerish(ntransit, lower = 1, any.missing = FALSE)
  }
  .ui <- rxode2::assertRxUi(ui)
  rxode2::assertCompartmentExists(.ui, central)
  rxode2::assertCompartmentExists(.ui, depot)
  rxode2::assertCompartmentExists(.ui, paste0(transit, "1"))
  rxode2::assertVariableName(ktr)
  .mv <- rxode2::rxModelVars(.ui)
  .transitCmts <- .mv$state
  .transitCmts <- .transitCmts[grepl(paste0("^", transit), .transitCmts)]
  .nc <- nchar(transit) + 1
  .totTransit <- max(vapply(.transitCmts,
                            function(n) {
                              as.integer(substr(n, .nc, nchar(n)))
                            }, integer(1), USE.NAMES = FALSE))
  if (!missing(ntransit)) {
    checkmate::assertIntegerish(ntransit, lower = 1, any.missing = FALSE, len = 1)
  } else {
    ntransit <- .totTransit
  }
  if (ntransit > .totTransit) {
    warning("reset ntransit to ", .totTransit, call.=FALSE)
    ntransit <- .totTransit
  }
  .ui <- rxode2::rxUiDecompress(.ui)
  if (ntransit == .totTransit) {
    # remove all
    .tmp <- .getEtaThetaTheta1(.ui)
    .iniDf <- .tmp$iniDf
    .eta <- .tmp$eta
    .theta <- .tmp$theta
    .theta1 <- .tmp$theta1
    .theta <- .dropTheta(.theta, ktr)
    .eta <- .dropEta(.eta, ktr)

    .rm <- seq_len(.totTransit)
    .transit <- paste0(transit, .rm)
    .modelLines <- .rmDdt(.ui$lstExpr, .transit)
    .w <- .whichDdt(.modelLines, central)
    .tmp <- .extractModelLinesAtW(.modelLines, .w)
    .tmp$w <- .replaceMult(.tmp$w,
                           v1=paste0(transit, .totTransit), v2=ka,
                           ret=paste0(ka, "*", depot))
    .tmp$pre <- .replaceMult(.tmp$pre,
                             v1=depot, v2=ktr,
                             ret=paste0(ka, "*", depot))
    .modelLines <- c(.tmp$pre,
                     .tmp$w,
                     .tmp$post)
    .tmp <- .dropLines(.ui, .modelLines, .theta, .eta, ktr)
    .modelLines <- .tmp$modelLines
    .theta <- .tmp$theta
    .eta <- .tmp$eta
    .ui$iniDf <- rbind(.theta,
                       .eta)
  } else {
    # remove some, but not all
    .ftransit <- .totTransit - ntransit
    .transit <- paste0(transit, seq_len(.totTransit)[-seq_len(.ftransit)])
    .modelLines <- .rmDdt(.ui$lstExpr, .transit)
    .w <- .whichDdt(.modelLines, central)
    .tmp <- .extractModelLinesAtW(.modelLines, .w)
    .tmp$w <- .replaceMult(.tmp$w,
                           v1=paste0(transit, .totTransit), v2=ka,
                           ret=paste0(ka, "*", transit, .ftransit))
    .tmp$pre <- .replaceMult(.tmp$pre,
                           v1=paste0(transit, ntransit), v2=ktr,
                           ret=paste0(ka, "*", transit, .ftransit))
    .modelLines <- c(.tmp$pre,
                     .tmp$w,
                     .tmp$post)

  }
  if (exists("description", envir=.ui$meta)) {
    rm("description", envir=.ui$meta)
  }
  rxode2::model(.ui) <- .modelLines
  return(rxode2::rxUiCompress(.ui))
}
