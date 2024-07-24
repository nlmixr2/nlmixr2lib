#' Add linear indirect response model
#'
#' @param ui rxode2 model
#'
#' @param stim what type of stimulation indirect response model:
#'
#' - `in`: stimulation of input
#'
#' - `out`: stimulation of output
#'
#' @param inhib what type of inhibition indirect response model:
#'
#' - `in`: inhibition of input
#'
#' - `out`: inhibition of output
#'
#' @param ek simulation linear constant
#'
#' @param ik inhibition linear constant
#'
#' @param kin this is the kin parameter name
#'
#' @param kout this is the kout parameter name
#'
#' @param cc the concentration value
#'
#' @param effect the effect compartment that will be modeled
#'
#' @return model with linear indirect effect added;
#'
#' note that while linear indirect effects are not common, it allows
#' an easier hook to produce all sorts of standard effect curves, like
#' Emax/Imax, Hill, etc.
#'
#' @export
#' @author Matthew L. Fidler
#' @family idr, PK
#' @examples
#'
#' readModelDb("PK_2cmt_no_depot") |> addIndirectLin(stim="in")
#'
#' readModelDb("PK_2cmt_no_depot") |> addIndirectLin(stim="out")
#'
#' readModelDb("PK_2cmt_no_depot") |> addIndirectLin(inhib="in")
#'
#' readModelDb("PK_2cmt_no_depot") |> addIndirectLin(inhib="out")
#'
addIndirectLin <- function(ui,
                           stim=c("in", "out"),
                           inhib=c("in", "out"),
                           ek="Ek",
                           ik="Ik",
                           kin="kin", kout="kout",
                           cc="Cc",
                           R="R",
                           R0="R0",
                           effect="effect") {
  if ((missing(stim) && missing(inhib)) ||
        (!missing(stim) && !missing(inhib))) {
    stop("need to either 'stim' or 'inhib'",
         call.=FALSE)
  }
  .doStim <- FALSE
  if (!missing(stim)) {
    .doStim <- TRUE
    stim <- match.arg(stim)
  } else {
    inhib <- match.arg(inhib)
  }
  .effectSd <- paste0(effect, "Sd")
  .ui <- rxode2::assertRxUi(ui)
  rxode2::assertVariableExists(.ui, cc)
  rxode2::assertVariableNew(.ui, kin)
  rxode2::assertVariableNew(.ui, kout)
  rxode2::assertCompartmentNew(.ui, R)
  rxode2::assertVariableNew(.ui, effect)
  rxode2::assertVariableNew(.ui, .effectSd)
  .modelLines <- .ui$lstExpr
  if (.doStim) {
    rxode2::assertVariableNew(.ui, ek)
    print(stim)
    if (stim == "in") {
      .eff <- str2lang(paste0("d/dt(", R, ") <- ", kin,
                              "*(1+", ek, "*", cc, ") - ", kout, "*", R))
    } else {
      .eff <- str2lang(paste0("d/dt(", R, ") <- ", kin,
                              " - ", kout, "*", R, "*(1+", ek, "*", cc, ")"))
    }
  } else {
    rxode2::assertVariableNew(.ui,ik)
    if (inhib == "in") {
      .eff <- str2lang(paste0("d/dt(", R, ") <- ", kin,
                              "*(1-", ik, "*", cc, ") - ", kout, "*", R))
    } else {
      .eff <- str2lang(paste0("d/dt(", R, ") <- ", kin,
                              " - ", kout, "*", R, "*(1-", ik, "*", cc, ")"))
    }
  }
  .modelLines <- c(list(str2lang(paste0(kin, " <- exp(l", kin,")")),
                        str2lang(paste0(kout, " <- exp(l", kout,")")),
                        str2lang(ifelse(.doStim,
                                        paste0(ek, " <- exp(l", ek,")"),
                                        paste0(ik, " <- exp(l", ik,")")))),
                   .modelLines,
                   list(str2lang(paste0(R, "(0) <- ", kin, "/", kout)),
                        .eff,
                        str2lang(paste0(effect, " <- ", R)),
                        str2lang(paste0(effect, " ~ add(", .effectSd, ")"))
                        ))

  .tmp <- .getEtaThetaTheta1(.ui)
  .iniDf <- .tmp$iniDf
  .theta <- .tmp$theta
  .theta1 <- .tmp$theta1
  .eta <- .tmp$eta
  if(length(.theta$ntheta) == 0) {
    .ntheta <- 0
  } else {
    .ntheta <- max(.theta$ntheta)
  }

  .thetaKin <- .get1theta(kin, .theta1, .ntheta,
                          label=paste0("zero order response production(", kin, ")"))
  .ntheta <- .ntheta + 1

  .thetaKout <- .get1theta(kout, .theta1, .ntheta,
                           label=paste0("first order rate response loss (", kout, ")"))
  .ntheta <- .ntheta + 1

  if (.doStim) {
    .thetaK <-  .get1theta(ek, .theta1, .ntheta,
                           label=paste0("linear effect constant (", ek, ")"))
  } else {
    .thetaK <-  .get1theta(ik, .theta1, .ntheta,
                           lower= -Inf, upper=1,
                           label=paste0("linear inhibition constant (", ik, ")"))
  }

  .ntheta <- .ntheta + 1

  .thetaErr <-  .get1theta(.effectSd, .theta1, .ntheta,
                           lower=0,
                           label=paste0("additive error for ", effect),
                           name=.effectSd)
  .thetaErr$condition <- effect
  .thetaErr$err <- "add"

  .ui <- rxode2::rxUiDecompress(.ui)
  .ui$iniDf <- rbind(.theta,
                     .thetaKin,
                     .thetaKout,
                     .thetaK,
                     .thetaErr,
                     .eta)
  if (exists("description", envir=.ui$meta)) {
    rm("description", envir=.ui$meta)
  }
  rxode2::model(.ui) <- .modelLines
  .ui

}
