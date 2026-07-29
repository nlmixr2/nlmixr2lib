Voller_2017_phenobarbital <- function() {
  description <- "One-compartment first-order-absorption population PK model for phenobarbital in preterm and term newborns (Voller 2017), as packaged in DDMORE Foundation Model Repository entry DDMODEL00000256."
  reference <- paste(
    "Voller S, Pichlmeier U, Bauer-Brandl A, Kloft C (2017).",
    "Pharmacokinetics of phenobarbital in newborns:",
    "Towards model-based optimisation of the loading dose.",
    "European Journal of Pharmaceutical Sciences 109S:S90-S97.",
    "doi:10.1016/j.ejps.2017.05.026.",
    "DDMORE Foundation Model Repository: DDMODEL00000256.",
    sep = " "
  )
  vignette <- "Voller_2017_phenobarbital"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  ddmore_id    <- "DDMODEL00000256"
  replicate_of <- NULL

  covariateData <- list(
    WT = list(
      description        = "Current body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Drives a linear-deviation effect on V relative to a 2.70 kg",
        "newborn reference: VWEIGHT = 1 + 0.309 * (WT - 2.70).",
        "Distinct from WT_BIRTH (birth weight, time-fixed).",
        "Source data column WEIGHT carries the same kg unit; no conversion required."
      ),
      source_name        = "WEIGHT"
    ),
    WT_BIRTH = list(
      description        = "Birth weight (time-fixed per subject)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at birth. Drives a linear-deviation effect on CL relative to a",
        "2.59 kg newborn reference: CLBW = 1 + 0.369 * (WT_BIRTH - 2.59).",
        "Distinct from WT (current body weight, time-varying).",
        "Source data column BWEIGHT carries the same kg unit; no conversion required."
      ),
      source_name        = "BWEIGHT"
    ),
    PNA = list(
      description        = "Postnatal age (chronological time since birth)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Drives a linear-deviation effect on CL relative to a 4.50-day",
        "postnatal-age reference: CLAGE = 1 + 0.0533 * (PNA_days - 4.50).",
        "The canonical PNA is in months (per inst/references/covariate-columns.md);",
        "the source paper / DDMORE bundle reported AGE in days, so the model converts",
        "internally as PNA_days = PNA_months * 30.4375 and keeps the per-day slope and",
        "the 4.50-day reference offset on the source-paper days scale for traceability",
        "against the .lst final estimates.",
        "When postnatal age is < 1 day (immediately post-partum), supply a small",
        "positive value (e.g., the time fraction in days / 30.4375) rather than 0;",
        "the linear-deviation form does not divide by PNA, so PNA = 0 is mathematically",
        "well-defined but may sit below the model's calibration range."
      ),
      source_name        = "AGE (days; multiply by 1/30.4375 to obtain canonical PNA in months)"
    )
  )

  population <- list(
    n_subjects     = 53,
    n_studies      = 1,
    age_range      = "Preterm and term newborns; postnatal age (PNA) range not extractable from the DDMORE bundle (Voller 2017 PDF not on disk). The bundle's simulated dataset spans PNA 0-58 days across 5 representative subjects.",
    weight_range   = "Birth weight (BWEIGHT) range not extractable from the DDMORE bundle. The bundle's simulated dataset includes subjects from 0.8 kg (extreme preterm) to 4.2 kg (term).",
    sex_female_pct = "Not extractable from DDMORE bundle.",
    race_ethnicity = "Not extractable from DDMORE bundle.",
    disease_state  = "Preterm and term newborns receiving phenobarbital (typical clinical indication: prevention or treatment of neonatal seizures). The DDMORE bundle does not specify the indication or NICU setting.",
    dose_range     = paste(
      "Phenobarbital given as an IV loading dose (typically a short infusion to the",
      "central compartment) followed by oral maintenance doses to the gastrointestinal",
      "depot. Doses in the bundle's simulated dataset range from approximately 4-17 mg",
      "per dose. Real clinical regimens are typically a 15-20 mg/kg IV loading dose",
      "followed by 3-5 mg/kg/day oral maintenance, but per-subject mg/kg dosing is not",
      "extractable from the bundle without the publication."
    ),
    regions        = "Not extractable from DDMORE bundle.",
    notes          = paste(
      "Population description is reconstructed from the .mod / .lst $PROBLEM line",
      "('Phenobarbital PK in newborns'), the 53-subject / 229-observation totals from",
      "the .lst data-summary block, and the simulated event-table demographics. The",
      "full Voller 2017 publication PDF is not on disk under",
      "/home/bill/github/mab_human_consensus/literature/; detailed demographics",
      "(age range, weight range, sex distribution, race, indication, regional setting)",
      "could not be cross-checked. The DDMORE entry's RDF metadata describes the",
      "purpose as 'The PK of phenobarbital was quantified in preterm and term newborns,",
      "to optimize drug dosing.'"
    )
  )

  ini({
    # Structural typical values: from FINAL PARAMETER ESTIMATE block of
    # Output_real_run522.lst (DDMODEL00000256) after MINIMIZATION SUCCESSFUL
    # (OBJV 1129.151, 53 subjects / 229 observations). Reference subject:
    # PNA = 4.50 days, WT_BIRTH = 2.59 kg, WT = 2.70 kg.
    lcl <- log(0.00909); label("Typical CL at PNA = 4.50 d, WT_BIRTH = 2.59 kg (L/h)")  # .lst TH 1 = 9.09E-03
    lvc <- log(2.38);    label("Typical V at WT = 2.70 kg (L)")                          # .lst TH 2 = 2.38E+00

    # KA was FIXED at 50 1/h in the source $THETA (`50 FIX`); kept fixed here.
    lka <- fixed(log(50)); label("First-order absorption rate constant (1/h; FIXED)")    # .lst TH 7 = 50 (FIXED)

    # Oral bioavailability of the depot (gastrointestinal) compartment.
    # Source $THETA bound (0, 0.594, 1) is a bounded parameter; the typical
    # value is taken verbatim from the .lst final estimate. For simulation
    # within nlmixr2lib the value is held as exp(lfdepot); users re-fitting
    # on real data should re-introduce the (0, 1) bound via a logit
    # transformation in their fitting workflow.
    lfdepot <- log(0.594); label("Oral bioavailability of the depot compartment (fraction)")  # .lst TH 8 = 5.94E-01

    # Linear-deviation covariate slopes. Slopes are kept on the source-paper
    # units (per day for PNA, per kg for WT_BIRTH and WT) so the .lst final
    # estimates carry across without rescaling. The model converts the
    # canonical PNA (months) back to days at use site so the 0.0533/day
    # slope and the 4.50-day reference offset apply unchanged.
    e_pna_cl       <- 0.0533;  label("PNA slope on CL (per day; linear-deviation about 4.50 d)")  # .lst TH 4 = 5.33E-02
    e_wtbirth_cl   <- 0.369;   label("WT_BIRTH slope on CL (per kg; linear-deviation about 2.59 kg)")  # .lst TH 5 = 3.69E-01
    e_wt_vc        <- 0.309;   label("WT slope on V (per kg; linear-deviation about 2.70 kg)")  # .lst TH 6 = 3.09E-01

    # Inter-individual variability. NONMEM $OMEGA is a diagonal matrix in
    # this model (no BLOCK); ETA1 acts on CL, ETA2 acts on V.
    # Variances on the internal log-normal scale.
    etalcl ~ 0.0898   # .lst OMEGA(1,1) = 8.98E-02 (CL)
    etalvc ~ 0.0504   # .lst OMEGA(2,2) = 5.04E-02 (V)

    # Residual error: NONMEM $SIGMA is `1 FIX` and W = sqrt(THETA(3) * IPRED^2),
    # so Y = IPRED * (1 + EPS(1) * sqrt(THETA(3))) is a pure proportional
    # error with proportional SD = sqrt(THETA(3)). nlmixr2 expects the SD
    # directly via `prop(propSd)`.
    propSd <- sqrt(0.0258); label("Proportional residual error (fraction; SD scale)")  # .lst TH 3 = 2.58E-02 (variance)
  })
  model({
    # 1. Convert canonical PNA (months) back to source-paper days so that
    #    the 0.0533/day slope and the 4.50-day reference offset apply on
    #    the same numerical scale used in the .lst.
    pna_days <- PNA * 30.4375

    # 2. Linear-deviation covariate factors (.mod $PK; CLAGE / CLBW on CL,
    #    VWEIGHT on V). Centred at the source-paper reference values.
    clage <- 1 + e_pna_cl     * (pna_days - 4.50)
    clbw  <- 1 + e_wtbirth_cl * (WT_BIRTH - 2.59)
    vwt   <- 1 + e_wt_vc      * (WT       - 2.70)

    # 3. Individual PK parameters with covariate effects and IIV
    cl <- exp(lcl + etalcl) * clage * clbw
    vc <- exp(lvc + etalvc) * vwt
    ka <- exp(lka)

    # 4. Micro-constant for the 1-compartment ODE
    kel <- cl / vc

    # 5. ODE system: oral depot (CMT 1) -> central (CMT 2). IV doses go
    #    directly into central via cmt = 2 in the event table.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 6. Bioavailability applies only to the depot (oral) compartment;
    #    IV doses to central are not multiplied by F.
    f(depot) <- exp(lfdepot)

    # 7. Observation (mg/L; dose in mg, volume in L) and proportional residual error
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
