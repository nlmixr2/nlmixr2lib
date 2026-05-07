Vet_2016_midazolam <- function() {
  description <- "Two-compartment population PK model for IV midazolam in critically ill children (Vet 2016) with body-weight allometric scaling on CL and V1 (reference 5 kg), a CRP power effect on CL (reference 32 mg/L), per-stratum typical CL values for the number of failing organs (ORG_FAIL_COUNT strata 0 / 1 / 2 / 3 / >=4; the source NMTRAN dataset names this column ORGF), inter-individual variability on CL and V1, and inter-occasion variability on CL across six daily occasions. Packaged in DDMORE Foundation Model Repository entry DDMODEL00000249."
  reference <- paste(
    "Vet NJ, Brussee JM, de Hoog M, Mooij MG, Verlaat CWM, Jerchel IS,",
    "van Schaik RHN, Koch BCP, Tibboel D, Knibbe CAJ, de Wildt SN; SKIC",
    "(Dutch collaborative PICU research network) (2016).",
    "Inflammation and Organ Failure Severely Affect Midazolam Clearance",
    "in Critically Ill Children.",
    "Am J Respir Crit Care Med 194(1):58-66.",
    "doi:10.1164/rccm.201510-2114OC.",
    "DDMORE Foundation Model Repository: DDMODEL00000249.",
    sep = " "
  )
  vignette <- "Vet_2016_midazolam"
  units <- list(time = "hour", dosing = "ug", concentration = "ug/L")

  ddmore_id    <- "DDMODEL00000249"
  replicate_of <- NULL

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL and V1 with paper-estimated exponents (THETA(5) = 1.02 on CL, THETA(6) = 1.34 on V1) and a 5 kg reference weight (Vet 2016 .mod $PK lines 51-58). Time-varying within subject is permitted by the source code, but the bundled simulated dataset uses a single per-subject baseline value.",
      source_name        = "WT"
    ),
    CRP = list(
      description        = "C-reactive protein concentration; standard (non-hs) assay used in Vet 2016 to quantify systemic inflammation in critically ill children.",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL with paper-estimated exponent THETA(11) = -0.312 and reference 32 mg/L (Vet 2016 .mod $PK line 51). Time-varying within subject (re-evaluated each ICU day in the source dataset).",
      source_name        = "CRP"
    ),
    ORG_FAIL_COUNT = list(
      description        = "Number of organs failing in the critically ill child on the observation day, ascertained per-day and reported as the worst-of-day count.",
      units              = "(count)",
      type               = "categorical",
      reference_category = "0 (no organs failing)",
      notes              = "Decomposed inside `model()` into mutually exclusive binary indicators `orgf1`, `orgf2`, `orgf3`, `orgf_ge4` that select per-stratum typical CL values (Vet 2016 .mod $PK lines 51-55: `IF (ORGF.EQ.0) TVCL = THETA(1) * ...`, `IF (ORGF.EQ.1) TVCL = THETA(7) * ...`, etc., with a final `IF (ORGF.GT.3.5) TVCL = THETA(10) * ...` collapsing the 4-and-5-organ strata). The Vet 2016 source dataset records this column as `ORGF` with values in {0, 1, 2, 3, 4, 5}; the canonical column name in nlmixr2lib is `ORG_FAIL_COUNT` and the canonical decomposition treats >=4 as a single combined stratum to match the source's parameterization.",
      source_name        = "ORGF"
    ),
    OCC = list(
      description        = "Integer-valued ICU-day occasion indicator for inter-occasion-variability multiplexing (1 = day 1, 2 = day 2, ..., 6 = day 6 or later).",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "The Vet 2016 .mod derives OCC from cumulative TIME (hours since first dose) within subject: OCC = 1 + floor(TIME / 24) capped at 6, matching .mod $PK lines 24-29 (`OCC=1; IF(TIME.GE.24)OCC=2; ... ; IF(TIME.GE.120)OCC=6`). The DDMORE-bundled simulated dataset (Simulated_MidaCriticallyIll.csv) does not carry an OCC column; it is derived from the cumulative TIME by the same rule when assembling the simulation event table for the vignette. Decomposed inside `model()` into binary indicators `oc1` .. `oc6` that multiplex the six IOV etas on log-CL.",
      source_name        = "OCC"
    )
  )

  population <- list(
    n_subjects     = 83L,
    n_studies      = 1L,
    age_range      = "Not extractable from DDMORE bundle (Vet 2016 PDF not on disk).",
    weight_range   = "Not extractable from DDMORE bundle (Vet 2016 PDF not on disk).",
    weight_reference = "5 kg (allometric reference per .mod $PK)",
    sex_female_pct = "Not extractable from DDMORE bundle (Vet 2016 PDF not on disk).",
    race_ethnicity = "Not extractable from DDMORE bundle (Vet 2016 PDF not on disk).",
    disease_state  = "Critically ill paediatric patients receiving continuous IV midazolam in the paediatric intensive care unit (PICU). Inflammation (CRP) and number of failing organs (ORG_FAIL_COUNT, source-data column `ORGF`, 0..>=4) are the two retained covariates.",
    dose_range     = "Continuous IV infusion at clinically titrated rates; the bundled simulated dataset spans 300-22500 ug/h infusions over 1-3 day infusion durations following an initial bolus.",
    crp_reference  = "32 mg/L (CRP power-effect reference per .mod $PK)",
    regions        = "Netherlands (SKIC paediatric ICU research network).",
    notes          = "Population descriptors are derived from the DDMODEL00000249 RDF `model-has-description` (`Midazolam PK in critically ill pediatric patients, using inflammation (quantified as CRP concentrations) and number of organs failing are most important covariates`) and the .mod / .lst FINAL PARAMETER ESTIMATE block. The Vet 2016 publication itself (Am J Respir Crit Care Med 194(1):58-66, doi:10.1164/rccm.201510-2114OC) is not on disk in this worktree, so demographics here come from the DDMORE bundle metadata rather than the paper's Table 1; the absence of a paper cross-check is documented in the validation vignette's Errata."
  )

  ini({
    # Structural typical values from DDMODEL00000249 Output_real_OriginalModelCode.lst
    # FINAL PARAMETER ESTIMATE block (THETA vector) captured after `MINIMIZATION SUCCESSFUL`
    # at .lst line 386 (NO. OF FUNCTION EVALUATIONS USED: 773; OBJV 6301.530 at .lst line 439).
    # Reference weight is 5 kg; reference CRP is 32 mg/L (Vet 2016 .mod $PK lines 51-58).

    # ORG_FAIL_COUNT=0 typical CL (FIXED at .mod $THETA line 78 "1.6 FIX") for a 5 kg child with
    # CRP = 32 mg/L. Per-stratum CL values for ORG_FAIL_COUNT >= 1 are encoded as additive shifts
    # `e_orgf<k>_cl` on the log scale below; the reference category is ORG_FAIL_COUNT = 0.
    lcl <- fixed(log(1.60))   ; label("Typical CL for ORG_FAIL_COUNT = 0, WT = 5 kg, CRP = 32 mg/L (L/h, FIXED)") # .lst line 455 TH 1 = 1.60E+00 (THETA(1) FIX in .mod)
    lvc <- log(3.28)          ; label("Central volume of distribution V1 at WT = 5 kg (L)")             # .lst line 455 TH 2 = 3.28E+00
    lq  <- log(1.52)          ; label("Inter-compartmental clearance Q (L/h)")                          # .lst line 455 TH 3 = 1.52E+00
    lvp <- log(5.44)          ; label("Peripheral volume of distribution V2 (L)")                       # .lst line 455 TH 4 = 5.44E+00

    # Allometric exponents on body weight (estimated, not theory-fixed) for CL and V1.
    e_wt_cl  <- 1.02          ; label("WT exponent on CL (unitless), reference 5 kg")  # .lst line 455 TH 5 = 1.02E+00
    e_wt_vc  <- 1.34          ; label("WT exponent on V1 (unitless), reference 5 kg")  # .lst line 455 TH 6 = 1.34E+00

    # Per-stratum CL shifts on the log scale relative to the ORG_FAIL_COUNT = 0 reference. Each
    # value is `log(THETA(k) / THETA(1))` for k = 7, 8, 9, 10, with THETA(1) = 1.6 (the
    # ORG_FAIL_COUNT = 0 typical CL). The numerical shifts are pre-computed below to avoid an
    # arithmetic step in `ini()`; the .lst line 455 source THETAs are quoted in each
    # comment for traceability.
    e_orgf1_cl    <- -0.2154  ; label("Additive shift in log(CL) when ORG_FAIL_COUNT = 1 (unitless)")     # log(1.29  / 1.60) ; .lst TH 7  = 1.29E+00
    e_orgf2_cl    <- -0.5141  ; label("Additive shift in log(CL) when ORG_FAIL_COUNT = 2 (unitless)")     # log(0.957 / 1.60) ; .lst TH 8  = 9.57E-01
    e_orgf3_cl    <- -0.6420  ; label("Additive shift in log(CL) when ORG_FAIL_COUNT = 3 (unitless)")     # log(0.842 / 1.60) ; .lst TH 9  = 8.42E-01
    e_orgf_ge4_cl <- -0.8585  ; label("Additive shift in log(CL) when ORG_FAIL_COUNT >= 4 (unitless)")    # log(0.678 / 1.60) ; .lst TH 10 = 6.78E-01

    # CRP power effect on CL with reference 32 mg/L. Negative exponent reflects reduced
    # midazolam clearance under elevated systemic inflammation (Vet 2016 paper finding).
    e_crp_cl <- -0.312        ; label("CRP exponent on CL (unitless), reference 32 mg/L")  # .lst line 455 TH 11 = -3.12E-01

    # Inter-individual variability (IIV). NONMEM `$OMEGA` diagonals on CL and V1 only;
    # Q and V2 carry no IIV in the .mod. Variances reported on the log-normal internal
    # scale per NONMEM convention.
    etalcl ~ 0.345            # OMEGA(1,1) FINAL = 3.45E-01 ; .lst line 465
    etalvc ~ 1.19             # OMEGA(2,2) FINAL = 1.19E+00 ; .lst line 468

    # Inter-occasion variability (IOV) on log-CL across six ICU-day occasions. NONMEM
    # `$OMEGA BLOCK(1) 0.05` followed by five `$OMEGA BLOCK(1) SAME` re-uses the same
    # single-element variance across all six occasions; the FINAL estimate of that
    # shared variance is OMEGA(3,3) = 1.97E-01 (.lst lines 471-486 on the diagonal).
    # nlmixr2 has no `SAME` shortcut so each occasion gets its own eta with the
    # variance fixed to the shared estimate after the first (the Jonsson 2011 pattern).
    etaiov_cl_1 ~ 0.197       # OMEGA(3,3) FINAL = 1.97E-01 ; estimated occasion-1 IOV variance
    etaiov_cl_2 ~ fix(0.197)  # OMEGA(4,4) fixed equal to OMEGA(3,3) per `$OMEGA BLOCK(1) SAME`
    etaiov_cl_3 ~ fix(0.197)  # OMEGA(5,5) fixed equal to OMEGA(3,3) per `$OMEGA BLOCK(1) SAME`
    etaiov_cl_4 ~ fix(0.197)  # OMEGA(6,6) fixed equal to OMEGA(3,3) per `$OMEGA BLOCK(1) SAME`
    etaiov_cl_5 ~ fix(0.197)  # OMEGA(7,7) fixed equal to OMEGA(3,3) per `$OMEGA BLOCK(1) SAME`
    etaiov_cl_6 ~ fix(0.197)  # OMEGA(8,8) fixed equal to OMEGA(3,3) per `$OMEGA BLOCK(1) SAME`

    # Combined proportional + additive residual error on the linear (ug/L) scale, per
    # the .mod $ERROR block `Y = F * (1 + ERR(1)) + ERR(2)`. The .lst SIGMA block
    # reports variances; the SDs below are sqrt of those variances.
    propSd <- 0.313           ; label("Proportional residual error (fraction)")  # sqrt(SIGMA(1,1) = 9.77E-02) ; .lst line 496
    addSd  <- 0.371           ; label("Additive residual error (ug/L)")          # sqrt(SIGMA(2,2) = 1.38E-01) ; .lst line 499
  })

  model({
    # 1. Decompose the integer-valued ORG_FAIL_COUNT column (source-data alias `ORGF`
    # in the .mod $INPUT) into mutually exclusive binary indicators that select the
    # per-stratum CL shift. Reference category is ORG_FAIL_COUNT = 0 (all four
    # indicators evaluate to 0; no shift applied).
    orgf1    <- (ORG_FAIL_COUNT == 1)
    orgf2    <- (ORG_FAIL_COUNT == 2)
    orgf3    <- (ORG_FAIL_COUNT == 3)
    orgf_ge4 <- (ORG_FAIL_COUNT >= 4)

    cl_orgf_shift <- orgf1 * e_orgf1_cl + orgf2 * e_orgf2_cl + orgf3 * e_orgf3_cl + orgf_ge4 * e_orgf_ge4_cl

    # 2. Decompose the integer-valued OCC column into binary indicators for IOV
    # multiplexing on log-CL. Vet 2016 derives OCC from cumulative TIME (hours since
    # first dose); the vignette pre-computes the column when assembling events.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)

    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
              oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5 + oc6 * etaiov_cl_6

    # 3. Individual PK parameters with allometric WT scaling, the per-stratum
    # ORG_FAIL_COUNT CL shift, the CRP power effect on CL, and the IOV multiplexing
    # on log-CL.
    cl <- exp(lcl + etalcl + cl_orgf_shift + iov_cl) * (WT / 5)^e_wt_cl * (CRP / 32)^e_crp_cl
    vc <- exp(lvc + etalvc) * (WT / 5)^e_wt_vc
    q  <- exp(lq)
    vp <- exp(lvp)

    # 4. Micro-constants for the 2-compartment IV ODE (matches .mod $PK lines 65-67).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 5. ODE system (matches .mod $DES lines 70-71, with `central` <-> A(1) and
    # `peripheral1` <-> A(2)).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 6. Plasma midazolam concentration (ug/L; numerically equal to ng/mL because
    # AMT is in ug and Vc is in L, per the .mod $INPUT comments at lines 4-13).
    Cc <- central / vc

    # 7. Combined residual error model (matches .mod $ERROR `Y = F*(1+ERR(1)) + ERR(2)`).
    Cc ~ prop(propSd) + add(addSd)
  })
}
