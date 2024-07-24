#' Remove model lines
#'
#' @param modelLines model line expression
#' @param lhs left handed variable
#' @returns model with lhs assignments dropped without checking
#'
#' @noRd
#' @author Matthew L. Fidler
.removeLines <- function(modelLines, lhs) {
  # Look for assignment of elimination in model
  .w <- vapply(seq_along(modelLines),
               function(i) {
                 rxode2::.matchesLangTemplate(modelLines[[i]],
                                              str2lang(paste0(lhs,
                                                              " <- ."))) ||
                   rxode2::.matchesLangTemplate(modelLines[[i]],
                                                str2lang(paste0(lhs,
                                                                "=.")))
               }, logical(1), USE.NAMES = FALSE)
  .w <- which(!.w)
  .seq <- seq_along(modelLines)[.w]
  lapply(.seq,
         function(i) {modelLines[[i]]})
}

# .replaceMultC(tmp[[5]], str2lang("kel"), str2lang("central"), str2lang("(vm*central/vc)/(km+central/vc)"))

#' Replace multiplication expressions
#'
#' @param x expression to modify
#' @param v1 expression for the first multiplication term
#' @param v2 expression for the second multiplication term
#' @param ret return expression that v1*v2 will be replaced
#' @returns expression with multiplications replaced with ret
#'
#' @noRd
#' @author Matthew L. Fidler
.replaceMultC <- function(x, v1, v2, ret) {
  if (is.call(x)) {
    if (length(x) == 3 &&
          identical(x[[1]], quote(`*`))) {
      .neg <- FALSE
      if (length(x[[2]]) == 2 &&
            identical(x[[2]][[1]], quote(`+`))) {
        x[[2]] <- x[[2]][[2]]
      }
      if (length(x[[2]]) == 2 &&
            identical(x[[2]][[1]], quote(`-`))) {
        .neg <- TRUE
        x[[2]] <- x[[2]][[2]]
      }
      if ((identical(x[[2]], v1) &&
             identical(x[[3]], v2)) ||
            (identical(x[[3]], v1) &&
               identical(x[[2]], v2))) {
        if (.neg) {
          return(str2lang(paste0("-", deparse1(ret))))
        } else {
          return(ret)
        }
      }
      as.call(lapply(x, .replaceMultC, v1=v1, v2=v2, ret=ret))
    } else {
      as.call(lapply(x, .replaceMultC, v1=v1, v2=v2, ret=ret))
    }
  } else {
    x
  }
}

#' Replace multiplication expressions
#'
#' @param modelLines model lines to replace multiplication expressions
#' @param v1 variable one to replace
#' @param v2 variable two to replace
#' @param ret new expression with multiplication replaced
#' @returns modelLines with multiplication expressions removedes
#'
#' @noRd
#' @author Matthew L. Fidler
.replaceMult <- function(modelLines, v1, v2, ret) {
  if (!is.list(modelLines)) modelLines <- list(modelLines)
  .v1 <- str2lang(v1)
  .v2 <- str2lang(v2)
  .ret <- str2lang(ret)
  lapply(seq_along(modelLines),
         function(i) {
           as.call(.replaceMultC(modelLines[[i]], v1=.v1, v2=.v2, ret=.ret))
         })
}
#' Drop thetas in the theta section of an iniDf
#'
#' @param theta theta iniDf data.frame
#' @param pars  parameters to drop
#' @returns theta data frame
#'
#' @noRd
#' @author Matthew L. Fidler
.dropTheta <- function(theta, pars) {
  .theta <- theta
  .w <- which(.theta$name %in% pars)
  if (length(.w) > 0) {
    # These are directly estimated, drop
    .theta <- .theta[-.w,, drop = FALSE]
    .theta$ntheta <- seq_along(.theta$name)
  }
  .theta
}

#' Drop the etas in the eta iniDf data frame
#'
#' @param eta eta iniDf data.frame
#' @param pars parameters
#' @returns eta data frame
#'
#' @noRd
#' @author Matthew L. Fidler
.dropEta <- function(eta, pars) {
  .eta <- eta
  .w <- which(.eta$name %in% pars)
  if (length(.w) > 0) {
    # Convert to matrix and drop columns
    .e <- vapply(.w, function(i) {
      .eta$neta1[i]
    }, double(1), USE.NAMES=FALSE)
    .eta <- .eta[!(.eta$neta1 %in% .e),, drop = FALSE]
    .eta <- .eta[!(.eta$neta2 %in% .e),, drop = FALSE]
    if (length(.eta$name) > 0) {
      .eta <- .eta[order(.eta$neta1, .eta$neta2),]
      .eta$neta1 <- as.integer(factor(.eta$neta1))
      .eta$neta2 <- as.integer(factor(.eta$neta2))
    }
  }
  .eta
}

#' Drops a single line from a model and recursively removes items
#'
#' @param ui original ui
#' @param modelLines model lines
#' @param theta iniDf theta name
#' @param eta iniDf eta name
#' @param par1 a single parameter to remove from the model
#' @returns list of modelLines, theta, and eta
#'
#' @noRd
#' @author Matthew L. Fidler
.dropLine1 <- function(ui, modelLines, theta, eta, par1) {
  .line <- rxode2::modelExtract(ui, par1)
  if (length(.line) == 0) {
    return(list(modelLines=modelLines, theta=theta, eta=eta))
  }
  .modelLines <- .removeLines(modelLines, par1)
  .vars <- rxode2::rxModelVars(.line)$params
  .theta <- .dropTheta(theta, .vars)
  .eta <- .dropEta(eta, .vars)
  for (.v in .vars) {
    .ret <- .dropLine1(ui, .modelLines, .theta, .eta, .v)
    .theta <- .ret$theta
    .eta <- .ret$eta
    .modelLines <- .ret$modelLines
  }
  list(modelLines=.modelLines, theta=.theta, eta=.eta)
}
#' Drop the lines from the model
#'
#'
#' @param ui original ui
#' @param modelLines model lines where this will be dropped
#' @param theta theta section of iniDf
#' @param eta eta section of iniDf
#' @param vars variables to drop
#' @returns list of modelLines, theta, and eta
#' @noRd
#' @author Matthew L. Fidler
.dropLines <- function(ui, modelLines, theta, eta, vars) {
  .modelLines <- modelLines
  .theta <- theta
  .eta <- eta
  for (.v in vars) {
    .ret <- .dropLine1(ui, .modelLines, .theta, .eta, .v)
    .theta <- .ret$theta
    .eta <- .ret$eta
    .modelLines <- .ret$modelLines
  }
  list(modelLines=.modelLines, theta=.theta, eta=.eta)
}
#' Returns the iniDf, theta, eta and theta1 data frames
#'
#' @param ui rxode2 ui function
#' @return A list of iniDf, theta, eta, and theta1. If there is no
#'   population information the theta1 will be generated
#' @noRd
#' @author Matthew L. Fidler
.getEtaThetaTheta1 <- function(ui) {
  .ui <- ui
  .iniDf <- .ui$iniDf
  .eta <- .iniDf[!is.na(.iniDf$neta1),, drop = FALSE]
  .theta <- .iniDf[is.na(.iniDf$neta1),, drop = FALSE]
  if (length(.theta$ntheta) == 0) {
    .theta1 <- lapply(names(.theta),
                      function(n) {
                        switch(n,
                               ntheta=1L,
                               neta1=NA,
                               neta2=NA,
                               name="_dummy",
                               lower= -Inf,
                               est=0,
                               upper=Inf,
                               fix=FALSE,
                               label=NA_character_,
                               backTransform=NA_character_,
                               condition=NA_character_,
                               err=NA_character_,
                               NA)
                      })
    names(.theta1) <- names(.theta)
    .theta1 <- as.data.frame(.theta1)
  } else {
    .theta1 <- .theta[1, ]
  }
  list(iniDf=.iniDf, theta=.theta, theta1=.theta1, eta=.eta)
}
#' Get a single theta estimate
#'
#'
#' @param vm name of the estimate; will pre-pend l to this
#' @param theta1 theta1 dataset
#' @param ntheta number of thetas, will increment to ntheta+1
#' @param lower lower estimate, default -Inf
#' @param est estimate, default 0.1
#' @param upper upper estimate, default Inf
#' @param fix fixed default FALSE
#' @param label default NA_character_
#' @return a single theta for integration with iniDf
#' @noRd
#' @author Matthew L. Fidler
.get1theta <- function(vm, theta1, ntheta,
                       lower= -Inf, est=0.1, upper=Inf,
                       fix=FALSE, label=NA_character_) {
  .thetaVm <- theta1
  .thetaVm$ntheta <- ntheta + 1
  .thetaVm$name <- paste0("l", vm)
  .thetaVm$lower <-lower
  .thetaVm$est <- est
  .thetaVm$upper <- upper
  .thetaVm$fix <- fix
  .thetaVm$label <- label
  .thetaVm
}

#' Convert models from linear elimination to Michaelis-Menten elimination
#'
#' @param ui model to convert
#' @param central the central compartment where the elimination is present
#' @param elimination variable for the elimination constant in the
#'   model
#' @param vm variable name for Vmax in the model
#' @param km variable name for Km in the model
#' @param vc variable name for Vc in the model
#' @returns new model changing linear elimination to Michaelis-Menten elimination
#'
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' readModelDb("PK_1cmt_des") |> convertMM()
#'
#' readModelDb("PK_2cmt_des") |> convertMM()
#'
#' readModelDb("PK_3cmt_des") |> convertMM()
#'
#' readModelDb("PK_3cmt_des") |> removeDepot() |> convertMM()
convertMM <- function(ui, central="central",
                      elimination="kel",
                      vm="vm", km="km", vc="vc") {
  rxode2::assertVariableName(elimination)
  rxode2::assertVariableName(vm)
  rxode2::assertVariableName(km)
  rxode2::assertVariableName(vc)
  .ui <- rxode2::assertRxUi(ui)
  rxode2::assertCompartmentExists(.ui, central)
  .tmp <- .getEtaThetaTheta1(.ui)
  .iniDf <- .tmp$iniDf
  .theta <- .tmp$theta
  .theta1 <- .tmp$theta1
  .eta <- .tmp$eta
  .line <- rxode2::modelExtract(.ui, elimination)
  .modelLines <- .ui$lstExpr
  if (identical(.line, character(0))) {
    # line not in model; kel estimated directly
    .theta <- .dropTheta(.theta, elimination)
    .eta <- .dropEta(.eta, elimination)
  } else {
    .modelLines <- .removeLines(.modelLines, elimination)
    .vars <- rxode2::rxModelVars(.line)$params
    .vars <- .vars[!(.vars %in% c(vm, km, vc))]
    .theta <- .dropTheta(.theta, .vars)
    .eta <- .dropEta(.eta, .vars)
    .ret <- .dropLines(.ui, .modelLines, .theta, .eta, .vars)
    .modelLines <- .ret$modelLines
    .theta <- .ret$theta
    .eta <- .ret$eta
  }
  if(length(.theta$ntheta) == 0) {
    .ntheta <- 0
  } else {
    .ntheta <- max(.theta$ntheta)
  }
  .thetaVm <- .get1theta(vm, .theta1, .ntheta)
  .ntheta <- .ntheta + 1

  .thetakm <- .get1theta(km, .theta1, .ntheta)
  .ntheta <- .ntheta + 1

  .ui <- rxode2::rxUiDecompress(.ui)
  .ui$iniDf <- rbind(.theta,
        .thetaVm,
        .thetakm,
        .eta)
  .model <- c(list(str2lang(paste0(vm, " <- log(l", vm, ")")),
         str2lang(paste0(km, " <- log(l", km, ")"))),
    .replaceMult(.modelLines, elimination, central,
               paste0("(", vm, "*", central, "/", vc, ")/(", km,
                      "+", central, "/", vc, ")")))
  if (exists("description", envir=.ui$meta)) {
    rm("description", envir=.ui$meta)
  }
  rxode2::model(.ui) <- .model
  .ui
}
