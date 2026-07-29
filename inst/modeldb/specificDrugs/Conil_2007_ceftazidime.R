Conil_2007_ceftazidime <- function() {
  description <- "Two-compartment IV population PK model for ceftazidime in adult burn-ICU patients, with creatinine clearance on CL and sex / mechanical ventilation / creatinine clearance on the peripheral volume V2 (Conil 2007)"
  reference <- "Conil JM, Georges B, Lavit M, Laguerre J, Samii K, Houin G, Saivin S. A population pharmacokinetic approach to ceftazidime use in burn patients: influence of glomerular filtration, gender and mechanical ventilation. Br J Clin Pharmacol. 2007;64(1):27-35. doi:10.1111/j.1365-2125.2007.02857.x"
  vignette <- "Conil_2007_ceftazidime"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Source column SEX with men = 0, women = 1; orientation matches canonical SEXF directly, no transformation required. 38 male and 12 female subjects in the cohort (Conil 2007 Table 1).",
      source_name        = "SEX"
    ),
    MECH_VENT = list(
      description        = "Invasive mechanical-ventilation treatment-status indicator at study entry, 1 = on mechanical ventilation, 0 = not ventilated",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no mechanical ventilation)",
      notes              = "Source column VENT with values 0 (without) and 1 (with); orientation matches canonical MECH_VENT directly, no transformation required. Time-fixed at the subject level for the PK sampling window. 16 of 50 patients ventilated (Conil 2007 Table 1). Hypothesised mechanism: positive end-expiratory pressure (PEEP) reduces glomerular filtration rate and stimulates antidiuretic-hormone secretion, increasing extracellular fluid and the peripheral volume of distribution V2.",
      source_name        = "VENT"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CLCR. Computed by the Cockcroft-Gault equation in raw mL/min (NOT BSA-normalized to mL/min/1.73 m^2); a correction factor of 0.85 is applied to women per Cockcroft-Gault convention (Conil 2007 Methods, Estimation of population parameters). Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization, with the per-model description recording the assay form; precedent Delattre_2010_amikacin). Cohort mean 105.3 +/- 39.3 mL/min, range 33-191 (Conil 2007 Table 1). Enters CL as an additive linear effect and V2 as a chained multiplicative fractional-change term.",
      source_name        = "CLCR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 50L,
    n_studies      = 1L,
    age_range      = "15-90 years",
    age_median     = "52.3 years (mean +/- 20.7 SD)",
    weight_range   = "50-108 kg",
    weight_median  = "71.3 kg (mean +/- 13.5 SD)",
    sex_female_pct = 24,
    race_ethnicity = "Not reported (French university-hospital burn-ICU population, Toulouse)",
    disease_state  = "Burn injury during the secondary phase; treated for local infection or sepsis with ceftazidime; mean burned surface area 23 +/- 13.5% of total body surface, UBS index 64.8 +/- 50.0, Baux index 75.6 +/- 22.4, Tobiasen index 7.0 +/- 1.8. Renal function in the normal range (serum creatinine 76.5 +/- 21.8 umol/L; creatinine clearance 105.3 +/- 39.3 mL/min).",
    dose_range     = "Ceftazidime 6 g/24 h administered as either three 2 g doses every 8 h or six 1 g doses every 4 h; each dose given as a 20-minute IV infusion via electric syringe (Conil 2007 Methods, Subjects and sampling).",
    regions        = "France (Burns Unit, University Hospital of Toulouse-Rangueil)",
    mech_vent_pct  = 32,
    notes          = "Single-center observational study over 4 years, 50 patients with 237 serum ceftazidime concentrations (mean 4.7 samples per patient). Concentrations measured by reversed-phase HPLC with UV detection at 260 nm; LLOQ 1 mg/L, calibration linear 1-100 mg/L with CV < 15%. Population PK by NONMEM with Visual-NM; one-compartment (ADVAN1 TRANS2) vs two-compartment (ADVAN3 TRANS4) models compared by likelihood ratio and AIC, two-compartment with proportional error retained as the final structural model (OFV 1749.280 for the basic two-compartment vs 1824.628 for one-compartment). Final covariate model OFV 1618.37 after retaining CRCL on CL and SEX + VENT + CRCL on V2; covariate retention threshold delta-OFV > 10.83 (p < 0.0001, chi-squared 1 df) on backward elimination."
  )

  ini({
    # Structural parameters at the typical baseline (no covariates active);
    # Conil 2007 Results final-model equations (p. 31) and Table 3 values.
    lcl <- log(1.08);  label("Non-renal / baseline component of CL (intercept of the additive linear CL ~ CRCL model) (L/h)")  # Conil 2007 final-model eqn: CL = 1.08 + 0.0536 * CLCR (intercept term 1.08 L/h, Results p. 31 and Abstract)
    lvc <- log(18.81); label("Central volume of distribution V1 (L)")                                                          # Conil 2007 final-model eqn: V1 = 18.81 L (Results p. 31 and Abstract; Table 3 reports 18.8 in all four sex / ventilation strata)
    lvp <- log(2.69);  label("Peripheral volume of distribution V2 baseline scale factor (L)")                                  # Conil 2007 final-model eqn: TVV2 = 2.69 * (1 + 1.43*SEX) * (1 + 2.44*VENT) * (1 + 0.00414*CLCR) (Results p. 31 baseline scale 2.69 L)
    lq  <- log(6.881); label("Intercompartmental clearance Q (L/h)")                                                            # Conil 2007 final-model eqn: Q = 6.881 L/h (Results p. 31 and Abstract; Table 3 reports 6.9 in all four strata)

    # Covariate effects. Conil 2007 reports the CL covariate as an additive
    # linear effect on the linear scale (CL_typ = intercept + slope * CRCL),
    # and the V2 covariates as a chained product of (1 + e * COV) fractional
    # changes, NOT the more common (COV / ref)^exponent or exp(e * (COV - ref))
    # forms. The encoded values therefore reproduce the published equations
    # verbatim (see vignette Source trace for the four-strata reproduction
    # check against Table 3).
    e_crcl_cl       <- 0.0536;  label("CRCL slope on CL (L/h per mL/min)")          # Conil 2007 final-model eqn: 0.0536 (slope of additive linear CL ~ CRCL covariate term)
    e_sexf_vp       <- 1.43;    label("SEXF fractional effect on V2 (unitless)")    # Conil 2007 final-model eqn: 1.43 (multiplicative factor (1 + 1.43*SEXF) on V2 baseline; females have 2.43x the male contribution)
    e_mech_vent_vp  <- 2.44;    label("MECH_VENT fractional effect on V2 (unitless)") # Conil 2007 final-model eqn: 2.44 (multiplicative factor (1 + 2.44*MECH_VENT) on V2 baseline; ventilated patients have 3.44x the non-ventilated contribution)
    e_crcl_vp       <- 0.00414; label("CRCL fractional effect on V2 (per mL/min)") # Conil 2007 final-model eqn: 0.00414 (slope of (1 + 0.00414*CLCR) multiplicative term on V2 baseline)

    # Inter-individual variability. Final-model CL and V1 CV% reported on
    # Results p. 31 ("Interindividual variability in the clearance was
    # decreased to 16% and that of the central volume of distribution to
    # 13%"). The paper is silent on Q and V2 IIV for the final model and on
    # the final-model residual error; per the standing instruction (see
    # operator sidecar response 2026-06-17), the basic-model values for Q
    # IIV, V2 IIV, and the proportional residual are carried forward
    # unchanged. omega^2 = log(CV^2 + 1) for log-normal etas.
    etalcl ~ 0.02528 # log(0.16^2 + 1); 16% CV on CL (Conil 2007 Results p. 31, final-model value)
    etalvc ~ 0.01676 # log(0.13^2 + 1); 13% CV on V1 (Conil 2007 Results p. 31, final-model value)
    etalq  ~ 2.6520  # log(3.63^2 + 1); 363% CV on Q (Conil 2007 Results p. 31 basic-model value carried forward; final model silent on Q IIV)
    etalvp ~ 1.5280  # log(1.90^2 + 1); 190% CV on V2 (Conil 2007 Results p. 31 basic-model value carried forward; final model silent on V2 IIV)

    # Proportional residual error (Conil 2007 Results p. 31 basic-model
    # value 38%, carried forward to the final model per the operator's
    # standing instruction; the final-model section does not re-report the
    # residual error).
    propSd <- 0.38; label("Proportional residual error (fraction)") # Conil 2007 Results p. 31 basic-model proportional CV = 38%; final-model residual not separately reported
  })
  model({
    # Individual PK parameters. CL has an additive linear covariate term on
    # CRCL (intercept + slope*CRCL, on the linear scale per the paper's
    # equation), wrapped in exp(etalcl) to give log-normal IIV. V2 is a
    # chained product of fractional-change factors on a baseline scale of
    # 2.69 L, wrapped in exp(etalvp).
    cl <- (exp(lcl) + e_crcl_cl * CRCL) * exp(etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp) *
      (1 + e_sexf_vp      * SEXF)      *
      (1 + e_mech_vent_vp * MECH_VENT) *
      (1 + e_crcl_vp      * CRCL)      *
      exp(etalvp)
    q  <- exp(lq + etalq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg administered as a 20-minute IV infusion into central
    # (rate-encoded on the event row); volumes in L -> central/vc gives
    # mg/L (= ug/mL), matching the paper's reported concentration units
    # (range 1-214.7 mg/L, Discussion).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
