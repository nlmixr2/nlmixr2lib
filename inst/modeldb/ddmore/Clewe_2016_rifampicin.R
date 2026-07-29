Clewe_2016_rifampicin <- function() {
  description <- "Multistate Tuberculosis Pharmacometric (MTP) model for in vitro M. tuberculosis H37Rv natural growth (no drug effect): three bacterial states (fast-multiplying, slow-multiplying, non-multiplying) with Gompertz growth on F and time-varying F->S transfer (Clewe 2016; DDMODEL00000240, scenario 4 -- natural-growth backbone of the framework that the publication then couples with rifampicin exposure-response; the bundled .mod ships only the natural-growth scaffold)"
  reference <- paste(
    "Clewe O, Aulin L, Hu Y, Coates AR, Simonsson US. (2016).",
    "A multistate tuberculosis pharmacometric model: a framework for studying",
    "anti-tubercular drug effects in vitro.",
    "J Antimicrob Chemother 71(4):964-974.",
    "doi:10.1093/jac/dkv416.",
    "DDMORE Foundation Model Repository: DDMODEL00000240 (scenario 4)."
  )
  vignette <- "Clewe_2016_rifampicin"
  units <- list(
    time = "day",
    dosing = "(none; in vitro natural-growth scaffold, no drug effect encoded)",
    concentration = "CFU/mL (colony-forming units per mL of bacterial culture)"
  )
  ddmore_id <- "DDMODEL00000240"
  replicate_of <- NULL

  covariateData <- list()

  population <- list(
    n_subjects     = 12,
    n_studies      = 1,
    age_range      = "(not applicable; in vitro experiment)",
    weight_range   = "(not applicable; in vitro experiment)",
    sex_female_pct = NULL,
    disease_state  = "(in vitro M. tuberculosis H37Rv culture)",
    dose_range     = "(no drug; natural-growth scaffold only)",
    notes          = "In vitro time-kill experiments on Mycobacterium tuberculosis H37Rv (St George's University strain). The bundled DDMODEL00000240 .mod fits 12 replicate cultures (per the Output_real_MTP.lst ETABAR N = 12) followed for ~200 days under untreated growth conditions. Initial F (fast-multiplying), S (slow-multiplying) and N (non-multiplying) bacterial counts and the inter-state transfer rate constants are estimated from the CFU/mL trajectory; only IIV on initial F is identifiable. The Clewe 2016 publication then uses this scaffold + a separate rifampicin exposure-response layer (scenario 4 in Model_Accomodations.txt), but the rifampicin layer is NOT shipped in the DDMORE bundle; the model file therefore contains only the natural-growth scaffold."
  )

  ini({
    # All values are FINAL PARAMETER ESTIMATES from
    # Output_real_MTP.lst lines 386-411 (post `MINIMIZATION SUCCESSFUL`,
    # FOCE; see ddmore-source.md "Reading final estimates from .lst").
    # The Clewe 2016 publication is not on disk under
    # /home/bill/github/mab_human_consensus/literature; values therefore
    # have not been cross-checked against the paper's printed tables.

    lkg <- log(0.206)
    label("Growth rate of fast-multiplying bacteria (1/day)")
    # Output_real_MTP.lst TH 1 = 2.06E-01; .mod $THETA(1) = 0.206361 (kG).

    lkfslin <- log(0.166 / 100)
    label("Slope of time-varying F->S transfer (1/day^2; kFSLIN, scaled /100 in .mod)")
    # Output_real_MTP.lst TH 2 = 1.66E-01; .mod $THETA(2) = 0.1657 with KFSLIN = THETA(2)/100.

    lkfn <- log(0.897 / 1e6)
    label("Transfer rate F->N (1/day; kFN, scaled /1e6 in .mod)")
    # Output_real_MTP.lst TH 3 = 8.97E-01; .mod $THETA(3) = 0.9 with KFN = THETA(3)/1e6.

    lksf <- log(0.145 / 10)
    label("Transfer rate S->F (1/day; kSF, scaled /10 in .mod)")
    # Output_real_MTP.lst TH 4 = 1.45E-01; .mod $THETA(4) = 0.14478 with KSF = THETA(4)/10.

    lksn <- log(0.186)
    label("Transfer rate S->N (1/day; kSN)")
    # Output_real_MTP.lst TH 5 = 1.86E-01; .mod $THETA(5) = 0.185568.

    lkns <- log(0.123 / 100)
    label("Transfer rate N->S (1/day; kNS, scaled /100 in .mod)")
    # Output_real_MTP.lst TH 6 = 1.23E-01; .mod $THETA(6) = 0.1227 with KNS = THETA(6)/100.

    lbmax <- log(242 * 1e6)
    label("System carrying capacity (CFU/mL; Bmax, scaled *1e6 in .mod)")
    # Output_real_MTP.lst TH 7 = 2.42E+02; .mod $THETA(7) = 241.6170 with BMAX = THETA(7)*1e6.

    lf0 <- log(4.10)
    label("Typical initial F (fast-multiplying) bacterial count (CFU/mL)")
    # Output_real_MTP.lst TH 8 = 4.10E+00; .mod $THETA(8) = 4.109880.

    lrbase <- log(9770)
    label("Typical initial S (slow-multiplying) bacterial count (CFU/mL)")
    # Output_real_MTP.lst TH 9 = 9.77E+03; .mod $THETA(9) = 9770.730.

    # IIV on log(F0). NONMEM `F0 = TVF0 * EXP(ETA(1))` puts ETA on log scale,
    # so etalf0 ~ var translates the diagonal $OMEGA verbatim.
    etalf0 ~ 22.4
    # Output_real_MTP.lst OMEGA ETA1 = 2.24E+01 (variance, log-scale);
    # very wide IIV -- etabar p-value 0.06, 41% etashrink -- see vignette Errata.

    # Residual error: NONMEM `IPRED = LOG(F+S); Y = IPRED + EPS(1)` with
    # SIGMA = 0.160 (variance on log-scale residual). Per
    # naming-conventions.md Section "$ERROR block patterns" this maps to a
    # proportional error in linear space: cfu ~ prop(propSd) with
    # propSd = sqrt(SIGMA) = sqrt(0.160) = 0.400.
    propSd <- sqrt(0.160)
    label("Proportional residual error on culturable CFU (F+S) -- fraction; sqrt of SIGMA(1) on log-scale")
    # Output_real_MTP.lst SIGMA EPS1 = 1.60E-01 (variance);
    # propSd = sqrt(0.160) = 0.4 ~ 40% CV in linear space.
  })

  model({
    # Derived rate constants in absolute units (back out the .mod's
    # /100, /10, /1e6, *1e6 NONMEM-bookkeeping scalings).
    kg     <- exp(lkg)
    kfslin <- exp(lkfslin)
    kfn    <- exp(lkfn)
    ksf    <- exp(lksf)
    ksn    <- exp(lksn)
    kns    <- exp(lkns)
    bmax   <- exp(lbmax)
    f0     <- exp(lf0 + etalf0)
    rbase     <- exp(lrbase)

    # Initial bacterial population (CFU/mL).
    fbugs(0) <- f0
    sbugs(0) <- rbase
    nbugs(0) <- 0.00001

    # Gompertz growth of F bacteria; clamp to non-negative when the
    # combined population exceeds Bmax (population at carrying capacity
    # cannot grow further). Matches .mod $DES `IF(GROWTHFUNC.LT.0) GROWTHFUNC=0`.
    total      <- fbugs + sbugs + nbugs
    growthfunc <- kg * log(bmax / total)
    growthfunc <- max(growthfunc, 0)

    # Time-varying F->S transfer rate (kFS = kFSLIN * t).
    kfs <- kfslin * t

    # Three-state ODE system from .mod $DES.
    d/dt(fbugs) <- fbugs * growthfunc + ksf * sbugs - kfs * fbugs - kfn * fbugs
    d/dt(sbugs) <- kfs * fbugs + kns * nbugs - ksn * sbugs - ksf * sbugs
    d/dt(nbugs) <- ksn * sbugs + kfn * fbugs - kns * nbugs

    # Observation: culturable CFU = F + S (matches IPRED = LOG(A(1)+A(2))
    # in .mod $ERROR; non-multiplying N bacteria do not form colonies).
    cfu <- fbugs + sbugs
    cfu ~ prop(propSd)
  })
}
