Jonsson_2005_disufenton <- function() {
  description <- "Two-compartment intravenous PK model for disufenton sodium (NXY-059) in adult patients with acute ischaemic or haemorrhagic stroke (Jonsson 2005), as packaged in DDMORE Foundation Model Repository entry DDMODEL00000245. Continuous IV infusion (1-h loading + 71-h maintenance) with a piecewise-linear creatinine-clearance effect on CL (no effect at CLCR <= 40 mL/min, linear above 40) and a linear weight effect on the central volume of distribution (centered at 76 kg). Correlated inter-individual variability on CL and Vc and a log-transform-both-sides residual error model."
  reference <- paste(
    "Jonsson S, Cheng Y-F, Edenius C, Lees KR, Odergren T, Karlsson MO. (2005).",
    "Population pharmacokinetic modelling and estimation of dosing strategy for",
    "NXY-059, a nitrone being developed for stroke.",
    "Clin Pharmacokinet 44(8):863-878.",
    "doi:10.2165/00003088-200544080-00007.",
    "DDMORE Foundation Model Repository: DDMODEL00000245.",
    sep = " "
  )
  vignette <- "Jonsson_2005_disufenton"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")
  ddmore_id    <- "DDMODEL00000245"
  replicate_of <- NULL

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Linear-deviation effect on the central volume of distribution centered at 76 kg: vc = vc * (1 + e_wt_vc * (WT - 76)). The 76 kg reference is encoded directly in the DDMODEL00000245 .mod $PK block (line 37: `V1WT = THETA(7)*(WT-76.00)`); the rerun comment `;Rerun to estimate parameters at CLCR 70 and WT 75` describes the pop-typical patient discussed in the publication, not the parameterisation reference. Source data column WT is missing-coded as -99 in the bundle's NMTRAN dataset; for those rows the .mod sets V1WT = 0 (i.e., V1 falls back to the 76 kg typical value).",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Creatinine clearance (raw, measured; not BSA-normalised)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying (CRCL is updated within-subject across the maintenance",
        "infusion in the bundle's NMTRAN dataset). Piecewise-linear effect on",
        "CL with breakpoint at 40 mL/min: cl = cl * (1 + e_crcl_cl * max(0, CRCL - 40)).",
        "Below 40 mL/min CL takes its CRCL=40 typical value. Reference 40 is the",
        "breakpoint; the publication discusses parameter values at the population",
        "centre CRCL = 70 mL/min (TVCL = 4.54 L/h there). Deviation from canonical",
        "CRCL: the canonical register entry is BSA-normalised (mL/min/1.73 m^2),",
        "but Jonsson 2005 uses raw measured creatinine clearance in mL/min, and",
        "the slope (0.0187 / mL/min) was estimated under that raw-mL/min",
        "parameterisation. The bundle's data column is named CLCR with sentinel",
        "-99 for missing rows; the .mod imputes 61.48 mL/min for most missing",
        "rows (and 34.55 mL/min for ID 10, 124.22 mL/min for ID 21). When",
        "supplying CRCL on input, do not include the -99 sentinel."
      ),
      source_name        = "CLCR"
    )
  )

  population <- list(
    n_subjects     = 179L,
    n_studies      = 2L,
    age_range      = "34-92 years",
    weight_range   = NA_character_,
    weight_median  = "76 kg (parameterisation reference; close to the population mean)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = "Adults with acute ischaemic or haemorrhagic stroke. Renal function ranges from severe impairment to normal (estimated CRCL 20-143 mL/min).",
    dose_range     = "Continuous intravenous infusion of NXY-059 over 72 hours, comprising a 1-hour loading infusion followed by a 71-hour maintenance infusion. Maintenance infusion rate was individualised to the subject's baseline creatinine clearance.",
    crcl_range     = "20-143 mL/min (raw, measured)",
    regions        = NA_character_,
    notes          = "Demographics summarised from the DDMODEL00000245 RDF model-has-description-long abstract, which mirrors Jonsson 2005's Methods. Pooled across the SA-NXY-0003 and SA-NXY-0004 stroke trials; .lst reports 177 individuals contributing observations after EVID filtering of the bundle's 179-subject simulated dataset. The Jonsson 2005 paper itself is not on disk in this worktree, so weight, sex, and race breakdowns could not be cross-checked against the publication's Table 1; see the validation vignette's Errata for the full caveat list."
  )

  ini({
    # Structural PK parameters - DDMODEL00000245 Output_real_run111.lst FINAL
    # PARAMETER ESTIMATE block (line 243), captured after `MINIMIZATION SUCCESSFUL`
    # on line 186. The .mod (Executable_run111.mod $PK lines 31-46) parameterises
    # CL with a piecewise CLCR effect (slope above 40 mL/min only) and Vc with a
    # linear deviation centered at 76 kg, so THETA(5) is the CL "intercept"
    # at CLCR <= 40 mL/min and THETA(2) is the Vc at WT = 76 kg.
    lcl <- log(2.91)  ; label("Clearance for CRCL <= 40 mL/min, typical individual (CL, L/h)")              # THETA(5) FINAL = 2.91E+00
    lvc <- log(7.91)  ; label("Central volume of distribution at WT = 76 kg, typical individual (V1, L)")    # THETA(2) FINAL = 7.91E+00
    lq  <- log(13.1)  ; label("Inter-compartmental clearance, typical individual (Q, L/h)")                  # THETA(3) FINAL = 1.31E+01
    lvp <- log(7.17)  ; label("Peripheral volume of distribution, typical individual (V2, L)")               # THETA(4) FINAL = 7.17E+00

    # Covariate effects (linear-deviation form per .mod $PK lines 31-41)
    e_crcl_cl <- 0.0187 ; label("Slope of (CRCL - 40) on CL above the 40 mL/min breakpoint: cl = cl * (1 + e_crcl_cl * max(0, CRCL - 40)) (1/(mL/min))") # THETA(6) FINAL = 1.87E-02
    e_wt_vc   <- 0.0194 ; label("Linear WT effect on Vc: vc = vc * (1 + e_wt_vc * (WT - 76)) (1/kg)")         # THETA(7) FINAL = 1.94E-02

    # Inter-individual variability - $OMEGA BLOCK(2) on (ETA1 = log-CL, ETA2 = log-Vc).
    # Final estimates from Output_real_run111.lst lines 250-256:
    #   OM(1,1) = 5.43E-02, OM(2,1) = 2.55E-02, OM(2,2) = 1.62E-01.
    # On the log scale, var(eta) = log(1 + CV^2) implies CV(CL) = sqrt(exp(0.0543) - 1) = 0.236
    # and CV(Vc) = sqrt(exp(0.162) - 1) = 0.421, matching the publication's
    # 23% CV (CL) and 40% CV (Vc) point estimates.
    etalcl + etalvc ~ c(0.0543, 0.0255, 0.162)  # NONMEM $OMEGA BLOCK(2): var(CL), cov(CL,Vc), var(Vc)

    # Residual error - $ERROR block uses log-transform-both-sides
    # (Y = LOG(F) + EPS(1) * THETA(1) with $SIGMA 1 FIX), so the SD on the
    # log-concentration scale equals THETA(1). NONMEM "additive on log scale"
    # is equivalent to nlmixr2's lnorm() residual; see naming-conventions.md
    # section  "$ERROR block patterns" and the existing Netterberg_2017_docetaxel /
    # Wu_2024_inotuzumab models for the same translation.
    propSd <- 0.165 ; label("Log-normal residual error SD on log-concentration (unitless; ~16.5% CV on linear scale)")  # THETA(1) FINAL = 1.65E-01
  })

  model({
    # Individual PK parameters with the linear-deviation covariate effects from
    # the .mod $PK block. The CLCR effect is a hockey-stick: zero contribution
    # below the 40 mL/min breakpoint, slope * (CRCL - 40) above it.
    cl <- exp(lcl + etalcl) * (1 + e_crcl_cl * max(0, CRCL - 40))
    vc <- exp(lvc + etalvc) * (1 + e_wt_vc   * (WT   - 76))
    q  <- exp(lq)
    vp <- exp(lvp)

    # Micro-constants for the 2-compartment IV ODE
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system (ADVAN3 TRANS4 in the .mod, mapped to canonical compartment names)
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Plasma concentration (dose mg, vc L -> mg/L = ug/mL)
    Cc <- central / vc

    # Log-normal residual error: Y = log(Cc) + eps with eps ~ N(0, propSd^2),
    # equivalent to NONMEM's Y = LOG(F) + EPS(1) * THETA(1) under SIGMA = 1 fix.
    Cc ~ lnorm(propSd)
  })
}
