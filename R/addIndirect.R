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
#' @return model with linear indirect effect added
#'
#' Note that while linear indirect effects are not common, it allows an easier
#' hook to produce other standard effect curves like Emax/Imax, Hill, etc.
#'
#' @export
#' @author Matthew L. Fidler
#' @family Indirect response
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
    stop("need either 'stim' or 'inhib' specified",
         call.=FALSE)
  }
  .doStim <- FALSE
  if (!missing(stim)) {
    if (missing(ui)) {
      return(fakeCc(addIndirectLin, stim=stim,
                    ek=ek, ik=ik,
                    kin=kin, kout=kout,
                    cc=cc,
                    R=R,
                    effect=effect))
    }
    .doStim <- TRUE
    stim <- match.arg(stim)
  } else {
    if (missing(ui)) {
      return(fakeCc(addIndirectLin, stim=stim,
                    ek=ek, ik=ik,
                    kin=kin, kout=kout,
                    cc=cc,
                    R=R,
                    effect=effect))
    }
    inhib <- match.arg(inhib)
  }
  .effectSd <- defaultCombine(effect, "sd")
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

#' Add an indirect response model to a PK model
#'
#' @inheritParams addIndirectLin
#' @inheritParams convertEmaxHill
#' @param hill Boolean stating if a hill sigmoid coefficient will be added
#' @param imax maximum inhibitory concentration
#' @param ic50 concentration where half of the Imax occurs
#' @return pk model with indirect response model added
#' @export
#' @author Matthew L. Fidler
#' @family Indirect response
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
    stop("need either 'stim' or 'inhib' specified",
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
                            ik=ek,
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
                            ek=ik,
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
