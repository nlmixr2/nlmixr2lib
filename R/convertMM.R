#' Remove model lines
#'
#' @param modelLines model line expression
#'
#' @param lhs left handed variable
#'
#' @return model with lhs assignments dropped without checking
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
#' @return expression with multiplications replaced with ret
#' @noRd
#' @author Matthew L. Fidler
.replaceMultC <- function(x, v1, v2, ret) {
  if (is.call(x)) {
    if (length(x) == 3 &&
          identical(x[[1]], quote(`*`)) &&
          ((identical(x[[2]], v1) &&
              identical(x[[3]], v2)) ||
             (identical(x[[3]], v1) &&
                identical(x[[2]], v2)))) {
      return(ret)
    } else {
      return(as.call(lapply(x, .replaceMultC, v1=v1, v2=v2, ret=ret)))
    }
  } else {
    return(x)
  }
}
#' replace multiplication expressions
#'
#' @param modelLines model lines to replace multiplication expressions
#' @param v1 variable one to replace
#' @param v2 variable two to replace
#' @param ret new expression with multiplication replaced
#' @return modelLines with multiplication expressions removedes
#' @noRd
#' @author Matthew L. Fidler
.replaceMult <- function(modelLines, v1, v2, ret) {
  .v1 <- str2lang(v1)
  .v2 <- str2lang(v2)
  .ret <- str2lang(ret)
  lapply(seq_along(modelLines),
         function(i) {
           .replaceMultC(modelLines[[i]], v1=.v1, v2=.v2, ret=.ret)
         })
}

#' Convert models from linear elimination to MM elimination
#'
#' @param ui model to convert
#'
#' @param central the central compartment where the elimination is present
#'
#' @param elimination variable for the elimination constant in the
#'   model
#'
#' @param vm variable name for Vmax in the model
#'
#' @param km variable name for Km in the model
#'
#' @param vc variable name for Vc in the model
#'
#' @return new model changing linear elimination to MM elimination
#'
#' @export
#'
#' @author Matthew L. Fidler
#'
#' @examples
#'
#' readModelDb("PK_1cmt_des") |> convertMM()
#'
#' readModelDb("PK_2cmt_des") |> convertMM()
#'
#' readModelDb("PK_3cmt_des") |> convertMM()
#'
#' readModelDb("PK_3cmt_des") |> removeDepot() |> convertMM()
#'
convertMM <- function(ui, central="central",
                      elimination="kel",
                      vm="vm", km="km", vc="vc",
                      cl="cl") {
  ui <- rxode2::assertRxUi(ui)
  rxode2::assertCompartmentExists(ui, central)
  ui <- eval(str2lang(paste0("rxode2::model(ui, -", cl, ")")))
  .modelLines <- ui$lstExpr
  .l <- length(.modelLines)
  .modelLines <- .removeLines(.modelLines, elimination)
  if (length(.modelLines) == .l) {
    stop("assumes '", elimination, "' is a derived variable",
         call.=FALSE)
  }
  .iniDf <- ui$iniDf
  .eta <- .iniDf[!is.na(.iniDf$neta1),, drop = FALSE]
  .theta <- .iniDf[is.na(.iniDf$neta1),, drop = FALSE]
  if (length(.theta$ntheta) == 0) {
    stop("need to have at least one population/residual parameter in the model", call.=FALSE)
  }
  .ntheta <- max(.theta$ntheta)
  .theta1 <- .theta[1, ]
  .thetaVm <- .theta1
  .thetaVm$name <- paste0("l", vm)
  .thetaVm$lower <- -Inf
  .thetaVm$est <- 0.1
  .thetaVm$upper <- Inf
  .thetaVm$fix <- FALSE
  .thetaVm$label <- NA_character_
  .thetakm <- .theta1
  .thetakm$name <- paste0("l", km)
  .thetakm$lower <- -Inf
  .thetakm$est <- 0.1
  .thetakm$upper <- Inf
  .thetakm$fix <- FALSE
  .thetakm$label <- NA_character_
  ui <- rxode2::rxUiDecompress(ui)
  ui$iniDf <- rbind(.theta,
        .thetaVm,
        .thetakm,
        .eta)
  .model <- c(list(str2lang(paste0(vm, " <- log(l", vm, ")")),
         str2lang(paste0(km, " <- log(l", km, ")"))),
    .replaceMult(.modelLines, elimination, central,
               paste0("(", vm, "*", central, "/", vc, ")/(", km,
                      "+", central, "/", vc, ")")))
  if (exists("description", envir=ui$meta)) {
    rm("description", envir=ui$meta)
  }
  rxode2::model(ui) <- .model
  ui
}
