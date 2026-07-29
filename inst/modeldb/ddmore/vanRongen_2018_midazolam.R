vanRongen_2018_midazolam <- function() {
  description <- "Two-compartment population PK model for midazolam with a five-transit-compartment first-order oral absorption chain (KA = KTR), supporting oral and intravenous dosing, in 19 obese adolescents (median total body weight 102.7 kg) and 20 morbidly obese adults (median 144 kg). Adult and adolescent typical clearances are estimated as separate intercepts; total body weight enters as a power covariate on clearance only in adolescents (reference 104.7 kg) and on peripheral volume only in adults (reference 141.8 kg)."
  reference <- paste(
    "van Rongen A, Brill MJE, Vaughns JD, Valitalo PAJ, van Dongen EPA,",
    "van Ramshorst B, Barrett JS, van den Anker JN, Knibbe CAJ (2018).",
    "Higher Midazolam Clearance in Obese Adolescents Compared with",
    "Morbidly Obese Adults.",
    "Clin Pharmacokinet 57(5):601-611.",
    "doi:10.1007/s40262-017-0579-4.",
    "DDMORE Foundation Model Repository: DDMODEL00000250.",
    sep = " "
  )
  vignette  <- "vanRongen_2018_midazolam"
  units     <- list(time = "minute", dosing = "microgram", concentration = "microgram/L")
  ddmore_id    <- "DDMODEL00000250"
  replicate_of <- NULL

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column 'TBW' (total body weight in kg) maps to the canonical WT.",
        "Time-fixed at baseline. Power-form covariate with two distinct reference",
        "weights, applied conditionally on the ADOLESCENT indicator: in adolescents",
        "(ADOLESCENT == 1) WT enters CL as `(WT/104.7)^e_wt_cl`; in adults",
        "(ADOLESCENT == 0) WT enters peripheral volume as `(WT/141.8)^e_wt_vp`.",
        "The 104.7 kg and 141.8 kg reference values come from the .mod source",
        "(Executable_FinalModelCode.mod) and are close to but not identical to",
        "the published cohort medians (102.7 kg adolescents, 144 kg adults; van Rongen",
        "2018 Abstract)."
      ),
      source_name        = "TBW"
    ),
    ADOLESCENT = list(
      description        = "Adolescent indicator (1 = obese adolescent, 0 = morbidly obese adult)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (morbidly obese adult)",
      notes              = paste(
        "The source .mod encodes the indicator implicitly through the subject ID:",
        "ID <= 30 are adults, ID > 30 are adolescents (.mod $INPUT comment line and",
        "$PK 'IF (ID.LE.30) ... IF (ID.GT.30) ...' branches). Users supplying their",
        "own dataset must pre-compute the canonical 0/1 ADOLESCENT column from the",
        "subject's age category (no implicit derivation from ID is performed here)."
      ),
      source_name        = "ID"
    )
  )

  population <- list(
    n_subjects     = 39L,
    n_studies      = 1L,
    age_range      = "Obese adolescents and morbidly obese adults; numeric age range not reproduced in the DDMORE bundle and not in the publication abstract",
    weight_range   = "Adolescents: median 102.7 kg (62-149.5 kg); morbidly obese adults: median 144 kg (112-186 kg)",
    weight_median  = "Adolescents 102.7 kg; adults 144 kg",
    sex_female_pct = NA_real_,
    disease_state  = "Combined cohort of obese adolescents and morbidly obese adults (CYP3A activity study); no specific concurrent disease named in the publication abstract.",
    dose_range     = "Combined oral and intravenous midazolam dosing. The DDMORE bundle's Simulated_DatafileMidaObesity.csv ships subjects receiving an oral dose of 7500 microgram followed approximately 2.5 hours later by an intravenous bolus of 5000 microgram (RATE 29997 microgram/min, ~10 s); the publication's exact dose levels are not reproduced in the bundle.",
    regions        = NA_character_,
    notes          = paste(
      "Demographic detail summarized from the linked publication's PubMed abstract",
      "(PMID 28785981, doi:10.1007/s40262-017-0579-4) and the DDMORE Foundation Model",
      "Repository bundle (DDMODEL00000250). The full publication PDF is not on disk",
      "in this worktree; numeric age range, sex distribution, and per-region",
      "breakdown are not in the abstract and are recorded as NA. The simulated dataset",
      "shipped with the bundle has 9 representative subjects (4 adults, 5 adolescents)",
      "and is intended as a regression-test cohort, not a representative population."
    )
  )

  ini({
    # Final estimates from Output_real_FinalModelCode.lst FINAL PARAMETER ESTIMATE
    # block (after MINIMIZATION SUCCESSFUL, line 353; objective function 524.451,
    # line 407). The DDMORE bundle re-fits the model on the shipped
    # Simulated_DatafileMidaObesity.csv (9 subjects, 138 records), so these final
    # estimates differ slightly from the publication's reported point values
    # (van Rongen 2018 abstract: CL adults 0.44 L/min, CL adolescents 0.71 L/min,
    # which match the .mod $THETA initial values). See vignette Errata for the
    # full bundle-vs-publication comparison and the near-boundary OMEGA(1,1)
    # caveat.

    # Structural PK parameters
    lcl     <- log(0.540)  ; label("Typical clearance in morbidly obese adults (L/min)")                                # FINAL THETA TH 1
    lvc     <- log(57.3)   ; label("Central volume of distribution V2 (L)")                                             # FINAL THETA TH 3
    lq      <- log(1.27)   ; label("Inter-compartmental clearance Q (L/min)")                                           # FINAL THETA TH 4
    lvp     <- log(166)    ; label("Peripheral volume of distribution V3 at the adult reference WT = 141.8 kg (L)")     # FINAL THETA TH 5
    lka     <- log(0.111)  ; label("First-order absorption / transit rate KA = KTR (1/min)")                            # FINAL THETA TH 6
    lfdepot <- log(0.684)  ; label("Oral bioavailability F1 (unitless)")                                                # FINAL THETA TH 2

    # Sub-population shift on log-CL: adolescent typical CL at WT = 104.7 kg
    # is THETA(7) = 0.793 L/min vs. adult typical CL THETA(1) = 0.540 L/min; the
    # log-ratio log(0.793 / 0.540) is the additive shift on log(CL) when ADOLESCENT == 1.
    e_adolescent_cl <- log(0.793 / 0.540) ; label("Log-additive ADOLESCENT effect on CL: log(THETA(7) / THETA(1)) at WT = 104.7 kg (unitless)") # derived from FINAL THETA TH 7 / TH 1

    # Power-form weight covariate effects, applied conditionally
    e_wt_cl <- 1.05        ; label("Power exponent of (WT / 104.7) on CL in obese adolescents only (unitless)")        # FINAL THETA TH 8
    e_wt_vp <- 3.36        ; label("Power exponent of (WT / 141.8) on peripheral volume in morbidly obese adults only (unitless)") # FINAL THETA TH 9

    # IIV (variances on the log-normal eta scale; OMEGA diagonal in the .lst).
    # OMEGA(1,1) = 4.33e-06 is at the lower boundary in the bundle's re-fit
    # (NONMEM "PARAMETER ESTIMATE IS NEAR ITS BOUNDARY" warning, .lst line 356;
    # ETA shrinkage 98.97% per .lst line 368). The eta on CL is retained for
    # structural fidelity to the source .mod but contributes negligibly to
    # simulation variability.
    etalcl     ~ 4.33e-06    # FINAL OMEGA(1,1); near-boundary
    etalfdepot ~ 9.22e-02    # FINAL OMEGA(2,2)
    etalq      ~ 4.32e-01    # FINAL OMEGA(3,3)
    etalka     ~ 8.61e-02    # FINAL OMEGA(4,4)
    etalvc     ~ 1.86e-02    # FINAL OMEGA(5,5)
    etalvp     ~ 7.96e-02    # FINAL OMEGA(6,6)

    # Residual error: NONMEM Y = F * (1 + ERR(1)) with var(ERR) = SIGMA(1,1) = 9.88e-02
    # is proportional in linear space; propSd = sqrt(SIGMA(1,1)) is the linear-scale SD.
    propSd <- sqrt(9.88e-02) ; label("Proportional residual error (fraction of IPRED)")                                 # FINAL SIGMA(1,1)
  })

  model({
    # Sub-group switch: ADOLESCENT is supplied per subject as 0 (morbidly obese adult)
    # or 1 (obese adolescent). The .mod uses ID <= 30 vs ID > 30 as the implicit gate;
    # users must materialize the ADOLESCENT column themselves (see covariateData).

    # Individual PK parameters with sub-group-conditional weight covariate effects
    # (.mod $PK lines 26-35: TVCL by group, TVV3 by group, KA = absorption = transit).
    cl <- exp(lcl + etalcl + e_adolescent_cl * ADOLESCENT) *
          (WT / 104.7)^(e_wt_cl * ADOLESCENT)
    vp <- exp(lvp + etalvp) *
          (WT / 141.8)^(e_wt_vp * (1 - ADOLESCENT))
    vc  <- exp(lvc + etalvc)
    q   <- exp(lq  + etalq)
    ka  <- exp(lka + etalka)
    ktr <- ka                                         # .mod sets KTR = KA explicitly

    # Micro-constants for the central / peripheral compartments. Two-compartment
    # disposition: K23 = Q/V2, K32 = Q/V3, K20 = CL/V2 in the .mod's $PK block.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Compartment ordering matches the .mod $MODEL block so user data CMT codes
    # carry over: CMT 1 = depot (PODOSE oral input), CMT 2 = central (IV input
    # and the observed-concentration compartment), CMT 3 = peripheral1, and
    # CMT 4..8 = transit1..transit5. Oral doses load depot, IV bolus / infusion
    # loads central directly.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ktr * transit5 - kel * central -
                          k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(transit1)    <-  ka  * depot     - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1  - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2  - ktr * transit3
    d/dt(transit4)    <-  ktr * transit3  - ktr * transit4
    d/dt(transit5)    <-  ktr * transit4  - ktr * transit5

    # Oral bioavailability F1 applied to the depot only; IV doses (CMT == 2)
    # bypass this scaling, matching the .mod's `F1 = THETA(2) * EXP(ETA(2))`.
    f(depot) <- exp(lfdepot + etalfdepot)

    # Observed plasma midazolam concentration (microgram/L = ng/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
