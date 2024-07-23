#' Extract model lines at W
#'
#' @param modelLines modelLines list
#' @param w where to split the modelLines
#' @return model lines list with:
#'   list(pre=pre model lines, w=model lines, post=post model lines)
#' @noRd
#' @author Matthew L. Fidler
.extractModelLinesAtW <- function(modelLines, w) {
  checkmate::assertInteger(w, lower=1, len=1)
  .w <- w
  #.pre will be list() if .w is at 1
  .pre <- lapply(seq(1, .w)[-.w],
                 function(i) {
                   modelLines[[i]]
                 })
  # .post will be list() if .w is at the end of the line
  .post <- lapply(seq(.w, length(modelLines))[-1],
                  function(i) {
                    modelLines[[i]]
                  })
  list(pre=.pre, w=modelLines[[.w]], post=.post)

}
#' Figures out what line d/dt(central) is on
#'
#'
#' @param modelLines modelLines expression list
#' @param central name of central compartment
#' @return which item in modelLines is the central compartment (or
#'   error if there is multiple lines)
#' @noRd
#' @author Matthew L. Fidler
.whichDdt <- function(modelLines, central) {
  .ddtCentral1 <- str2lang(paste0("d/dt(", central, ") <- ."))
  .ddtCentral2 <- str2lang(paste0("d/dt(", central, ") = ."))
  .w <- which(vapply(seq_along(modelLines),
                     function(i) {
                       .cur <- modelLines[[i]]
                       rxode2::.matchesLangTemplate(.cur, .ddtCentral1) ||
                         rxode2::.matchesLangTemplate(.cur, .ddtCentral2)
                     }, logical(1), USE.NAMES = FALSE))
  # Modify ODE for central compartment
  if (length(.w) != 1) stop("'d/dt(", central, ")' must be on a single line")
  .w
}

#' To add transit compartments to the model
#' @param ui The model as a function
#' @param ntransit the number of transit compartments to be added
#' @param transit the transit compartment prefix
#' @param ktr the parameter name for the transit compartment rate
#' @inheritParams addComp
#' @return a model with transit compartment added
#' @export
#' @examples
#' readModelDb("PK_1cmt_des") |> addTransit(3)
addTransit <- function(ui, ntransit, central = "central",
                       depot = "depot",
                       transit = "transit",
                       ktr = "ktr",
                       ka="ka") {
  .ui <- rxode2::assertRxUi(ui)
  rxode2::assertCompartmentExists(.ui, central)
  rxode2::assertCompartmentExists(.ui, depot)
  rxode2::assertCompartmentName(transit)
  rxode2::assertVariableName(ka)
  checkmate::assertIntegerish(ntransit, lower = 1)
  rxode2::assertVariableName(ktr)
  .mv <- rxode2::rxModelVars(.ui)

  # Extract model and central ODE
  .iniDf <- .ui$iniDf
  .eta <- .iniDf[!is.na(.iniDf$neta1),, drop = FALSE]
  .theta <- .iniDf[is.na(.iniDf$neta1),, drop = FALSE]
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
                           ret=paste0(ktr, "*", transit, ntransit))
  .post <- .tmp$post
  .v <- seq_len(ntransit)
  .v <- .v[-1]
  # ODEs for the transit compartment (except the one from the depot)
  .transMid <- list()
  if (length(.v) > 0) {
    .transMid <- lapply(.v,
                        function(i) {
                          str2lang(paste0("d/dt(", transit, i, ")<- ",
                                          ktr, "*", transit, i - 1, "-", ktr,
                                          "*", transit, i))
                        })
  }

  # combine the lines for now
  .modelLines <- c(.pre,
                   list(str2lang(paste0("d/dt(", transit,
                                        "1) <- ", ka, "*", depot, "-",
                                        ktr, "*", transit, "1"))),
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
                                ret=paste0(ka, "*", depot, " - ",
                                           ktr, "*", transit, "1")),
                   .tmp$post)
  if (length(.theta$name == 0L)) {
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
#' @param model The model as a function
#' @param transit The number of transit compartments to remove
#' @inheritParams addTransit
#' @inheritParams addComp
#' @return rxode2 model with transit compartment removed
#' @export
#' @examples
#'
#' # In this example the transit is added and then a few are removed
#'
#' readModelDb("PK_1cmt_des") |>
#'   addTransit(4) |>
#'   removeTransit(3)
removeTransit <- function(model, transit, central = "central", depot = "depot", transitComp = "transit", ktr = "ktr") {
  assertCompartmentName(central)
  assertCompartmentName(depot)
  assertCompartmentName(transitComp)
  assertVariableName(ktr)

  if (!missing(transit)) {
    checkmate::assertIntegerish(transit, lower = 1, any.missing = FALSE, len = 1)
  }
  temp <- rxode2::assertRxUi(model)

  mv <- rxode2::rxModelVars(temp)

  if (!(central %in% mv$state)) {
    stop("'", central, "' needs to be in the model")
  }
  if (!(any(grepl("^transit", mv$state)))) {
    stop("'", transitComp, " need to be in the model")
  }
  if (!(depot %in% mv$state)) {
    stop("'", depot, "' needs to be in the model")
  }
  # Extract model
  modelNew <- rxode2::modelExtract(temp, endpoint = NA)

  # modify ODE for central compartment
  center <- eval(str2lang(paste0("rxode2::modelExtract(temp,d/dt(", central, "),lines=TRUE)")))
  rhs <- sub(".*<-\\s*", "", center)
  # could be made less fragile by using transitive property
  rhs <- sub(paste0("\\s*", ktr, "\\s*\\*\\s*", transitComp, "\\d*"), "", rhs)

  # Find total number of transit compartments in the model
  totalTransit <- sum(grepl(paste0("\\s*^", transitComp, "[1-9][0-9]*"), mv$state))

  if (missing(transit)) {
    transit <- totalTransit
    line <- str2lang(paste0("d/dt(", central, ") <- ", rhs))
  } else {
    line <- str2lang(paste0("d/dt(", central, ") <- ", ktr, "*", transitComp, (totalTransit - transit), deparse(str2lang(rhs))))
  }


  # Modify ini{}
  if (transit == totalTransit) {
    # remove parameter
    rxode2::ini(temp) <- temp$iniDf[which(temp$iniDf$name != "lktr"), ]
  }

  # Modify model{}
  obj <- NULL
  indices <- totalTransit:(totalTransit - transit + 1)
  obj <- unlist(lapply(indices, function(i) {
    obj1 <- eval(str2lang(paste0("rxode2::modelExtract(temp, d/dt(transit", i, "), lines = TRUE)")))
    obj <- c(obj, obj1)
  }))


  obj2 <- eval(str2lang(paste0("rxode2::modelExtract(temp,", ktr, ",lines = TRUE)")))
  if (transit == totalTransit) {
    obj <- c(obj, obj2)
  }


  for (i in obj) {
    index <- which(modelNew == i)
    modelNew <- modelNew[-index]
  }

  rxode2::model(temp) <- modelNew

  # modify ODE for central compartment
  temp <- rxode2::model(temp, line)

  # return
  temp
}
