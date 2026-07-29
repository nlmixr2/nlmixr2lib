Wilkins_2008_rifampicin <- function() {
  description <- "One-compartment population PK model for oral rifampicin in adult South African pulmonary tuberculosis patients (Wilkins 2008), with an analytical transit-compartment chain (Savic 2007 form) preceding first-order absorption, multiplicative formulation effects of single-drug-combination (vs fixed-dose-combination reference) on apparent oral clearance and on mean transit time, IIV on CL/V (correlated)/Ka/MTT/NN, and 6-occasion inter-occasion variability on log-CL and log-MTT."
  reference <- paste(
    "Wilkins JJ, Savic RM, Karlsson MO, Langdon G, McIlleron H, Pillai G, Smith PJ,",
    "Simonsson US. (2008). Population pharmacokinetics of rifampin in pulmonary",
    "tuberculosis patients, including a semimechanistic model to describe variable",
    "absorption. Antimicrob Agents Chemother 52(6):2138-2148.",
    "doi:10.1128/AAC.00461-07.",
    "DDMORE Foundation Model Repository: DDMODEL00000280.",
    sep = " "
  )
  vignette <- "Wilkins_2008_rifampicin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")
  ddmore_id    <- "DDMODEL00000280"
  replicate_of <- NULL

  covariateData <- list(
    FORM_FDC = list(
      description        = "Fixed-dose-combination antitubercular formulation indicator (1 = FDC tablet, 0 = single-drug-combination, separate tablets).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (FDC; the most-common formulation in the Wilkins 2008 cohort and the typical-value reference for `lmtt` / `lcl`)",
      notes              = "Source `.mod` $PK block uses `IF(FDC.EQ.1) MTTFDC = 0` and `IF(FDC.EQ.0) MTTFDC = THETA(8)` (and an analogous block for CLLOC vs THETA(9)), so the FDC = 1 subgroup carries the typical-value `lmtt` / `lcl` and the SDC subgroup (FDC = 0) gets a `(1 + theta) *` multiplicative shift. The .mod's `CLLOC` variable name is a vestigial label from earlier model-development; the actual covariate driving CLCOV is `FDC` (not the data-only `LOC` column kept in $INPUT but unused in $PK). Multiplicative shifts: SDC subjects get 1 + 1.04 = 2.04 x MTT (longer absorption) and 1 + 0.236 = 1.236 x CL (faster clearance) than the FDC reference. The `LOC` data column is therefore not declared as a covariate in this model file.",
      source_name        = "FDC"
    ),
    OCC = list(
      description        = "Integer-valued occasion / period indicator for inter-occasion-variability multiplexing.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Values 1, 2, 3, 4, 5, 6 identify the dosing / sampling occasion within subject (the Wilkins 2008 cohort had up to 6 sampling occasions per subject across the antitubercular treatment course). Decomposed inside `model()` into binary indicators `oc1` .. `oc6` that multiplex the 6 IOV etas on log-CL and the 6 IOV etas on log-MTT.",
      source_name        = "OCC"
    )
  )

  population <- list(
    n_subjects     = 263L,
    n_studies      = 1L,
    age_range      = "Adults with pulmonary tuberculosis; demographic detail not extracted (Wilkins 2008 publication PDF not on disk in this worktree).",
    weight_range   = "(not extracted; Wilkins 2008 publication not on disk for cross-check)",
    sex_female_pct = "(not extracted)",
    disease_state  = "Adults with newly-diagnosed pulmonary tuberculosis treated with a standard antitubercular backbone (rifampicin + isoniazid + pyrazinamide +/- ethambutol).",
    dose_range     = "Oral rifampicin 450, 480, or 600 mg daily, multiple-dose at steady-state. The DDMORE bundle's Simulated_TB_Rifampicin_PK_Wilkins_2008.csv carries 250 simulated subjects across these three dose levels with sampling on six occasions (day 1, week 1, 2, 4, 8, 24-ish).",
    regions        = "South Africa (Wilkins 2008 cohort).",
    n_observations = 2913L,
    notes          = "Population descriptors are inferred from the bundle metadata (DDMODEL00000280.rdf model-has-description-long, the .lst's `TOT. NO. OF INDIVIDUALS: 263 / TOT. NO. OF OBS RECS: 2913` header, and the .mod $INPUT column list) because the Wilkins 2008 publication PDF is not on disk in /home/bill/github/mab_human_consensus/literature/ at extraction time. n_subjects = 263 is taken directly from the .lst header. The cohort's age / weight / sex breakdown is not derivable from the bundle alone; the validation vignette's Errata section documents this caveat."
  )

  ini({
    # Structural PK parameters. Final estimates from the DDMODEL00000280 bundle's
    # `.mod` $THETA initial values (Executable_real_TB_Rifampicin_PK_Wilkins_2008.mod
    # lines 127-135), which equal the FINAL THETAs in TB_Rifampicin_PK_Wilkins_2008_real.lst
    # ("FINAL PARAMETER ESTIMATE" block, lines 244-246). The .mod's $PROBLEM line names
    # this run "RIF Run71 newinits" -- the THETAs are FINAL estimates from a previous
    # successful upstream run carried forward as initial values; the .lst's MINIMIZATION
    # TERMINATED status reflects a smoke-test re-run that did not progress past iteration 4
    # on a slightly different real dataset, NOT a non-converged fit. The DDMORE-curated
    # `Output_real_TB_Rifampicin_PK_Wilkins_2008` summary table reproduces the same values
    # with sensible RSEs (e.g. CL/F = 19.2 L/h with 1.29% RSE), confirming these are the
    # published Wilkins 2008 final estimates. See the validation vignette's Errata section
    # for the full caveat list (incl. the summary's likely typo on additive residual error).
    lcl  <- log(19.2);  label("Apparent oral clearance CL/F at FDC reference (L/h)")            # .mod $THETA(1) FINAL = 1.92E+01; matches Output_real_* summary CL/F = 19.2 (RSE 1.29%)
    lvc  <- log(53.2);  label("Apparent central volume of distribution V/F (L)")                # .mod $THETA(2) FINAL = 5.32E+01; matches Output_real_* summary V/F = 53.2 (RSE 1.16%)
    lka  <- log(1.15);  label("Absorption rate constant from depot to central, ka (1/h)")       # .mod $THETA(3) FINAL = 1.15E+00; matches Output_real_* summary ka = 1.15 (RSE 3.91%)
    lmtt <- log(0.424); label("Mean transit time through the absorption transit chain, MTT (h, FDC reference)") # .mod $THETA(6) FINAL = 4.24E-01; matches Output_real_* summary MTT = 0.424 (RSE 3.82%)
    lnn  <- log(7.13);  label("Number of transit compartments NN (continuous, dimensionless)")  # .mod $THETA(7) FINAL = 7.13E+00; matches Output_real_* summary n = 7.13 (RSE 8.42%)

    # Formulation covariate effects. NONMEM source uses multiplicative `(1 + theta * (1 - FORM_FDC))`
    # shifts on MTT (THETA(8)) and on CL (THETA(9)); both effects vanish for the FDC = 1
    # typical-value reference. The .mod's CLLOC variable name is vestigial -- the conditional
    # is `IF(FDC.EQ.0) CLLOC = THETA(9)`, so the effect is FDC-driven, not LOC-driven.
    e_fdc0_mtt <- 1.04;  label("Multiplicative MTT increase for SDC (FDC = 0) vs FDC = 1 reference (unitless)")  # .mod $THETA(8) FINAL = 1.04E+00; Output_real_* summary 'SDC on MTT' = 1.04 (RSE 8.78%)
    e_fdc0_cl  <- 0.236; label("Multiplicative CL increase for SDC (FDC = 0) vs FDC = 1 reference (unitless)")   # .mod $THETA(9) FINAL = 2.36E-01; Output_real_* summary 'SDC on CL'  = 0.236 (RSE 9.63%)

    # Inter-individual variability (IIV). NONMEM `$OMEGA BLOCK(2)` correlates ETA(1)/ETA(2)
    # (CL/V); diagonal $OMEGA on ETA(3) (KA), ETA(4) (MTT), ETA(17) (NN). nlmixr2 BLOCK(2)
    # form expects lower-triangle order in c(): (var_CL, cov, var_V).
    etalcl + etalvc ~ c(0.279, 0.217, 0.188)  # .mod $OMEGA BLOCK(2) on ETA(1)/ETA(2); .lst FINAL OMEGA(1,1) = 2.79E-01, OMEGA(1,2) = 2.17E-01, OMEGA(2,2) = 1.88E-01
    etalka  ~ 0.439  # .mod $OMEGA on ETA(3); .lst FINAL OMEGA(3,3) = 4.39E-01 -- IIV Ka  (log-normal variance)
    etalmtt ~ 0.361  # .mod $OMEGA on ETA(4); .lst FINAL OMEGA(4,4) = 3.61E-01 -- IIV MTT (log-normal variance)
    etalnn  ~ 2.44   # .mod $OMEGA on ETA(17); .lst FINAL OMEGA(17,17) = 2.44E+00 -- IIV NN (log-normal variance)

    # Inter-occasion variability (IOV) on log-CL across 6 occasions. NONMEM source declares
    # `$OMEGA BLOCK(1) 0.0508` followed by five `BLOCK(1) SAME` re-uses the same single-element
    # variance across the six occasions; the FINAL estimate of that shared variance is
    # OMEGA(5,5) = 5.08E-02. nlmixr2 has no `SAME` shortcut, so each occasion gets its own
    # eta with the variance fixed to the shared value after the first (matching the
    # Jonsson_2011_ethambutol pattern).
    etaiov_cl_1 ~ 0.0508         # .mod $OMEGA BLOCK(1) on ETA(5); .lst FINAL OMEGA(5,5) = 5.08E-02 -- estimated occasion-1 IOV variance
    etaiov_cl_2 ~ fix(0.0508)    # $OMEGA BLOCK(1) SAME on ETA(6); fixed equal to OMEGA(5,5)
    etaiov_cl_3 ~ fix(0.0508)    # $OMEGA BLOCK(1) SAME on ETA(7); fixed equal to OMEGA(5,5)
    etaiov_cl_4 ~ fix(0.0508)    # $OMEGA BLOCK(1) SAME on ETA(8); fixed equal to OMEGA(5,5)
    etaiov_cl_5 ~ fix(0.0508)    # $OMEGA BLOCK(1) SAME on ETA(9); fixed equal to OMEGA(5,5)
    etaiov_cl_6 ~ fix(0.0508)    # $OMEGA BLOCK(1) SAME on ETA(10); fixed equal to OMEGA(5,5)

    # IOV on log-MTT across 6 occasions. Same `BLOCK(1) + 5 SAME` pattern; FINAL shared
    # variance is OMEGA(11,11) = 4.61E-01.
    etaiov_mtt_1 ~ 0.461         # .mod $OMEGA BLOCK(1) on ETA(11); .lst FINAL OMEGA(11,11) = 4.61E-01 -- estimated occasion-1 IOV variance
    etaiov_mtt_2 ~ fix(0.461)    # $OMEGA BLOCK(1) SAME on ETA(12)
    etaiov_mtt_3 ~ fix(0.461)    # $OMEGA BLOCK(1) SAME on ETA(13)
    etaiov_mtt_4 ~ fix(0.461)    # $OMEGA BLOCK(1) SAME on ETA(14)
    etaiov_mtt_5 ~ fix(0.461)    # $OMEGA BLOCK(1) SAME on ETA(15)
    etaiov_mtt_6 ~ fix(0.461)    # $OMEGA BLOCK(1) SAME on ETA(16)

    # Combined add+prop residual error on the linear (mg/L) scale. The .mod $ERROR block
    # uses `W = SQRT(THETA(4)**2 + THETA(5)**2 * F * F)` and `Y = IPRED + W * EPS(1)` with
    # `$SIGMA 1 FIX`; this is nlmixr2's default Pythagorean-SD `combined2` form with
    # additive SD = THETA(4) and proportional SD = THETA(5).
    addSd  <- 0.0923;  label("Additive residual error (mg/L)")             # .mod $THETA(4) FINAL = 9.23E-02; .lst FINAL TH 4 = 9.23E-02. The DDMORE Output_real_* summary file lists `Additive residual error 0.0508 (mg/L)` for the same parameter -- likely a transcription typo (the row immediately above reads `IOV on CL 0.0508 (variance)`); see vignette Errata.
    propSd <- 0.222;   label("Proportional residual error (fraction)")     # .mod $THETA(5) FINAL = 2.22E-01; matches Output_real_* summary 'Proportional residual error 0.222' (RSE 2.89%)
  })

  model({
    # Decompose the integer-valued occasion column into binary indicators for IOV
    # multiplexing on log-CL and log-MTT (matches the .mod `IOVCL` / `IOVMT` `IF(OCC.EQ.k)`
    # ladder in $PK lines 40-61).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)

    iov_cl  <- oc1 * etaiov_cl_1  + oc2 * etaiov_cl_2  + oc3 * etaiov_cl_3 +
               oc4 * etaiov_cl_4  + oc5 * etaiov_cl_5  + oc6 * etaiov_cl_6
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 + oc3 * etaiov_mtt_3 +
               oc4 * etaiov_mtt_4 + oc5 * etaiov_mtt_5 + oc6 * etaiov_mtt_6

    # Formulation effect: SDC subjects (FDC = 0) get a (1 + theta) multiplicative shift
    # on MTT and on CL, vanishing for the FDC = 1 typical-value reference (see .mod $PK
    # lines 18-34: `IF(FDC.EQ.1) MTTFDC = 0 ; Most common ... MTTCOV = (1 + MTTFDC)` and
    # the analogous `CLLOC` / `CLCOV` block).
    fdc0    <- 1 - FORM_FDC
    cl_cov  <- 1 + e_fdc0_cl  * fdc0
    mtt_cov <- 1 + e_fdc0_mtt * fdc0

    # Individual PK parameters. CL and MTT carry both IIV and the per-occasion IOV eta;
    # V, Ka, NN carry only IIV. NN is positive-constrained so the `lnn` log-transform
    # tracks the .mod's `THETA(7)` lower bound at 1 (`$THETA (1, 7.13, 80)`).
    cl  <- exp(lcl  + etalcl  + iov_cl)  * cl_cov
    vc  <- exp(lvc  + etalvc)
    ka  <- exp(lka  + etalka)
    mtt <- exp(lmtt + etalmtt + iov_mtt) * mtt_cov
    nn  <- exp(lnn  + etalnn)

    kel <- cl / vc

    # 1-compartment PK with an analytical transit-compartment chain (Savic JG, Jonsson EN,
    # Wahlby U, Karlsson MO. J Pharmacokinet Pharmacodyn 2007;34:711-726). The .mod $DES
    # implements the closed-form transit-chain absorption rate using a Stirling expansion
    # of `log(NN!)` to permit non-integer NN; rxode2's built-in `transit(n, mtt, bio)`
    # function emits the same gamma-PDF input rate (verified algebraically against the
    # .mod's `EXP(LOG(PD) + LOG(KTR) + NN*LOG(KTR*(T-TDOS)) - KTR*(T-TDOS) - L)` form
    # with `KTR = (NN + 1)/MTT` and `L = log(2.5066) + (NN+0.5)*log(NN) - NN +
    # log(1 + 1/(12*NN))`). The dose lands in `depot` with `f(depot) = 0` to suppress the
    # bolus content (matching the .mod's `F1 = 0`); `transit()` reads the raw dose amount
    # via `podo(depot)` regardless of `f(depot)` to drive the input rate.
    d/dt(depot)   <- transit(nn, mtt) - ka * depot
    d/dt(central) <-                    ka * depot - kel * central

    # Bioavailability fraction on the dosing compartment. Source $PK declares F1 = 0 so
    # no bolus enters depot; the entire dose is delivered via the analytical transit-chain
    # input rate above.
    f(depot) <- 0

    # Concentration in plasma. Dose units mg, V units L -> Cc units mg/L (= ug/mL).
    Cc <- central / vc

    # Combined additive + proportional residual error (Pythagorean / combined2 form,
    # matching the source $ERROR `W = SQRT(addSd^2 + propSd^2 * F^2)` linearization).
    Cc ~ add(addSd) + prop(propSd)
  })
}
