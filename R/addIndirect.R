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
#' @param R drug response compartment
#'
#' @param cc the concentration value
#'
#' @param effect the effect variable that will be modeled
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
  rxode2::assertExists(.ui, cc)
  rxode2::assertVariableNew(.ui, kin)
  rxode2::assertVariableNew(.ui, kout)
  rxode2::assertCompartmentNew(.ui, R)
  rxode2::assertVariableNew(.ui, effect)
  rxode2::assertVariableNew(.ui, .effectSd)
  .modelLines <- .ui$lstExpr
  if (.doStim) {
    rxode2::assertVariableNew(.ui, ek)
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
                        str2lang(paste0(effect, " ~ add(", .effectSd, ")"))))

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

#' Convert linear effect to Emax effect
#'
#' @param ui rxode2 model
#'
#' @param emax Emax parameter
#' @param ec50 EC50 parameter
#' @inheritParams  addIndirectLin
#'
#' @family PD
#'
#' @return Model with the linear effect converted to an Emax effect
#'
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' readModelDb("PK_2cmt_no_depot") |>
#'   addIndirectLin(stim="in") |>
#'   convertEmax()
#'
#' # When emax=1
#' readModelDb("PK_2cmt_no_depot") |>
#'   addIndirectLin(stim="in") |>
#'   convertEmax(emax=1)
#'
convertEmax <- function(ui, emax="Emax", ec50="EC50",
                        ek="Ek", cc="Cc") {
  .ui <- rxode2::assertRxUi(ui)
  rxode2::assertExists(.ui, cc)
  .emaxMult <- paste0(emax, "*")
  if (inherits(emax, "character")) {
    rxode2::assertVariableNew(.ui, emax)
  } else if (is.numeric(emax) && emax == 1.0) {
    .emaxMult <- ""
  } else {
    stop("'", emax, "' not specified correctly",
         call.=FALSE)
  }
  rxode2::assertVariableNew(.ui, ec50)
  rxode2::assertVariableExists(.ui, ek)
  .modelLines <- .replaceMult(.ui$lstExpr,
                              v1=ek, v2=cc,
                              ret=paste0(.emaxMult,
                                         cc, "/(", cc, "+", ec50, ")"))
  .tmp <- .getEtaThetaTheta1(.ui)
  .iniDf <- .tmp$iniDf
  .theta <- .tmp$theta
  .theta1 <- .tmp$theta1
  .eta <- .tmp$eta
  .tmp <- .dropLines(.ui, .modelLines, .theta, .eta, ek)
  .modelLines <- .tmp$modelLines
  .theta <- .tmp$theta
  .eta <- .tmp$eta

  if (length(.theta$ntheta) == 0) {
    .ntheta <- 0
  } else {
    .ntheta <- max(.theta$ntheta)
  }

  .emaxLine <- NULL
  if (.emaxMult != "") {
    .emaxLine <- list(str2lang(paste0(emax, "<- exp(l", emax, ")")))
    .thetaEmax <- .get1theta(emax, .theta1, .ntheta,
                             label=paste0("Maximum effect (", emax, ")"))
    .ntheta <- .ntheta + 1
  } else {
    .thetaEmax <- NULL
  }

  .thetaEc50 <- .get1theta(ec50, .theta1, .ntheta,
                           label=paste0("Concentration of 50% ", emax,
                                        " (", emax, ")"))
  .ntheta <- .ntheta + 1

  .ui <- rxode2::rxUiDecompress(.ui)
  .ui$iniDf <- rbind(.theta,
                     .thetaEmax,
                     .thetaEc50,
                     .eta)
  if (exists("description", envir=.ui$meta)) {
    rm("description", envir=.ui$meta)
  }
  rxode2::model(.ui) <- c(.emaxLine,
                          list(str2lang(paste0(ec50, "<- exp(l", ec50, ")"))),
                          .modelLines)
  .ui
}

#'  Convert linear effect to Emax-Hill effect
#'
#' @inheritParams convertEmax
#' @param g hill coefficient
#' @return Model with the linear effect converted to an Emax effect
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' readModelDb("PK_2cmt_no_depot") |>
#'   addIndirectLin(stim="in") |>
#'   convertEmaxHill()
#'
#' # can also specify as emax=1
#'
#' readModelDb("PK_2cmt_no_depot") |>
#'   addIndirectLin(stim="in") |>
#'   convertEmaxHill(emax=1)
#'
convertEmaxHill <- function(ui, emax="Emax", ec50="EC50", g="g",
                            ek="Ek", cc="Cc") {
  .ui <- rxode2::assertRxUi(ui)
  rxode2::assertExists(.ui, cc)
  .emaxMult <- paste0(emax, "*")
  if (inherits(emax, "character")) {
    rxode2::assertVariableNew(.ui, emax)
  } else if (is.numeric(emax) && emax == 1.0) {
    .emaxMult <- ""
  } else {
    stop("'", emax, "' not specified correctly",
         call.=FALSE)
  }
  rxode2::assertVariableNew(.ui, ec50)
  rxode2::assertVariableExists(.ui, ek)
  .modelLines <- .replaceMult(.ui$lstExpr,
                              v1=ek, v2=cc,
                              ret=paste0(.emaxMult, cc, "^", g,
                                         "/(", cc, "^", g,
                                         "+", ec50, "^", g, ")"))
  .tmp <- .getEtaThetaTheta1(.ui)
  .iniDf <- .tmp$iniDf
  .theta <- .tmp$theta
  .theta1 <- .tmp$theta1
  .eta <- .tmp$eta
  .tmp <- .dropLines(.ui, .modelLines, .theta, .eta, ek)
  .modelLines <- .tmp$modelLines
  .theta <- .tmp$theta
  .eta <- .tmp$eta

  if (length(.theta$ntheta) == 0) {
    .ntheta <- 0
  } else {
    .ntheta <- max(.theta$ntheta)
  }

  .emaxLine <- NULL
  if (.emaxMult != "") {
    .emaxLine <- list(str2lang(paste0(emax, "<- exp(l", emax, ")")))
    .thetaEmax <- .get1theta(emax, .theta1, .ntheta,
                             label=paste0("Maximum effect (", emax, ")"))
    .ntheta <- .ntheta + 1
  } else {
    .thetaEmax <- NULL
  }

  .thetaEc50 <- .get1theta(ec50, .theta1, .ntheta,
                           label=paste0("Concentration of 50% ", emax,
                                        " (", emax, ")"))
  .ntheta <- .ntheta + 1

  .thetaEc50 <- .get1theta(ec50, .theta1, .ntheta,
                           label=paste0("Concentration of 50% ", emax,
                                        " (", emax, ")"))
  .ntheta <- .ntheta + 1

  .thetaG <- .get1theta(g, .theta1, .ntheta, name=paste0("lg", g),
                        est=logit(1, 0.1, 10),
                        label=paste0("logit-constrained Hill coefficient ", g))
  .ntheta <- .ntheta + 1


  .ui <- rxode2::rxUiDecompress(.ui)
  .ui$iniDf <- rbind(.theta,
                     .thetaEmax,
                     .thetaEc50,
                     .thetaG,
                     .eta)
  if (exists("description", envir=.ui$meta)) {
    rm("description", envir=.ui$meta)
  }
  rxode2::model(.ui) <- c(.emaxLine,
                          list(str2lang(paste0(ec50, "<- exp(l", ec50, ")")),
                               str2lang(paste0(g, "<- expit(lg", g, ", 0.1, 10)"))),
                          .modelLines)
  .ui
}
#' Add an indirect response model to a PK model
#'
#' @inheritParams addIndirectLin
#' @inheritParams convertEmaxHill
#' @param hill boolean stating if a hill sigmoid cofficient will be added
#' @param imax maximum inhibitory concentration
#' @param ic50 concentration where half of the imax occurs
#' @return pk model with indirect response model added
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' readModelDb("PK_2cmt_no_depot") |>
#'   addIndirect(stim="in",hill=TRUE)
#'
#' readModelDb("PK_2cmt_no_depot") |>
#'   addIndirect(inhib="out", imax=1)
addIndirect <- function(ui,
                        stim=c("in", "out"),
                        inhib=c("in", "out"),
                        hill=FALSE,
                        ek="Ek",
                        ik="Ik",
                        emax="Emax",
                        ec50="EC50",
                        imax="Imax",
                        ic50="IC50",
                        kin="kin", kout="kout",
                        g="g",
                        cc="Cc",
                        R="R",
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
  if (.doStim) {
    .mod1 <- addIndirectLin(ui, stim=stim,
                            ek=ek,
                            ik=ik,
                            kin=kin, kout=kout,
                            cc=cc,
                            R=R,
                            effect=effect)
    if (hill) {
      convertEmaxHill(.mod1, emax=emax, ec50=ec50, g=g,
                      ek=ek, cc=cc)
    } else {
      convertEmax(.mod1, emax=emax, ec50=ec50,
                  ek=ek, cc=cc)
    }
  } else {
    .mod1 <- addIndirectLin(ui, inhib=inhib,
                            ek=ek,
                            ik=ik,
                            kin=kin, kout=kout,
                            cc=cc,
                            R=R,
                            effect=effect)
    if (hill) {
      convertEmaxHill(.mod1, emax=imax, ec50=ic50, g=g,
                      ek=ik, cc=cc)
    } else {
      convertEmax(.mod1, emax=imax, ec50=ic50,
                  ek=ik, cc=cc)
    }
  }
}