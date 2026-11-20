Sinnollareddy_2026_epcoritamab <- function() {
  description <- "Updated two-compartment quasi-steady-state (QSS) target-mediated drug disposition population PK model with first-order subcutaneous absorption and constant total target for epcoritamab (CD3xCD20 bispecific antibody) in adults with relapsed or refractory B cell non-Hodgkin lymphoma (Sinnollareddy 2026), re-estimated from the Li 2025 model with EPCORE NHL-1 expansion and optimization data. Body weight (reference 75 kg) scales CL/F, Q/F, Vc/F and Vp/F by power laws; age (reference 65 years) and body mass index (reference 25 kg/m^2) lower ka; tumour size (sum of the products of perpendicular diameters, reference 30 cm^2) raises and indolent NHL subtype lowers the drug-target complex internalization rate kint. Fixed effects and residual error are the published final estimates; the random effects are carried from the Li 2025 parent model because the update does not report them."
  reference <- paste(
    "Sinnollareddy M, Sanghavi K, Gibiansky L, Li T, Putnins M, Mohamed MEF,",
    "Patah P, Favaro E, Gupta M, Parikh A, Xu S. Exposure-Response Analyses to",
    "Inform the Optimal Epcoritamab Monotherapy Dosing Regimen in Relapsed or",
    "Refractory Follicular Lymphoma. Clin Pharmacokinet. 2026;65(9):1415-1428.",
    "doi:10.1007/s40262-026-01661-1. Final estimates in supplementary Table S2.",
    "Parent model (structure, covariate reference values and random effects):",
    "Li T, Gibiansky L, Parikh A, et al. Population Pharmacokinetics of",
    "Epcoritamab Following Subcutaneous Administration in Relapsed or Refractory",
    "B Cell Non-Hodgkin Lymphoma. Clin Pharmacokinet. 2025;64(1):127-141.",
    "doi:10.1007/s40262-024-01464-2 (modellib('Li_2025_epcoritamab'))."
  )
  vignette <- "Sinnollareddy_2026_epcoritamab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # etaRUV is the random effect on the magnitude of the residual error carried
  # from the Li 2025 parent model; it has no canonical IIV name.
  paper_specific_etas <- c("etaRUV")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power covariate on CL/F (estimated exponent 0.949), Q/F (fixed 0.75), Vc/F (estimated 0.641) and Vp/F (fixed 1) (supplementary Table S2). Reference 75 kg from the Li 2025 parent model (supplementary Fig. S1 there); Sinnollareddy 2026 does not restate it.",
      source_name = "WT"
    ),
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power covariate on ka (exponent -0.381, supplementary Table S2). Reference 65 years from the Li 2025 parent model.",
      source_name = "AGE"
    ),
    BMI = list(
      description = "Baseline body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power covariate on ka (exponent -0.689, supplementary Table S2), added in this update. The paper does not print the reference value; 25 kg/m^2 is a rounded value near the cohort medians (25.6-26.4 kg/m^2 in the EPCORE NHL-1 cohorts, supplementary Table S1). The vignette shows the geometric-mean ka this reference reproduces against Table 1 and Table S3.",
      source_name = "BMI"
    ),
    TUMSZ = list(
      description = "Baseline tumour size, sum of the products of perpendicular diameters (SUMPPD) of target lesions per investigator",
      units = "cm^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power covariate on kint (exponent 0.312, supplementary Table S2), added in this update. Units are cm^2 as in supplementary Table S1 (multiply by 100 for mm^2). The paper does not print the reference value; 30 cm^2 is a rounded value near the cohort medians (20.7-37.9 cm^2; natural-log median 3.48 in the EPCORE NHL-1 expansion cohort, i.e. 32.5 cm^2).",
      source_name = "SUMPPD"
    ),
    TUMTP_INHL = list(
      description = "Indolent non-Hodgkin lymphoma subtype indicator (1 = indolent B-NHL, including follicular lymphoma; 0 = large B cell lymphoma or mantle cell lymphoma)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (aggressive B-NHL: LBCL or MCL)",
      notes = "Multiplicative effect on kint: kint is multiplied by 0.484 for indolent NHL (supplementary Table S2, 'iNHL subtype effect on kint' = 0.484). The ratio reading is the maintainers' interpretation of an unprinted functional form; it is the only reading under which the published geometric-mean kint of follicular lymphoma patients (Table 1, Table S3) is reproduced with a plausible tumour-size reference -- see the vignette.",
      source_name = "iNHL"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "epcoritamab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(
      analyte = "epcoritamab (total: free + target-bound)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(analyte = "epcoritamab", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 466,
    n_studies = 2,
    age_range = "20-89 years",
    weight_range = "39-172 kg",
    bmi_range = "15.7-59.3 kg/m^2",
    sex_female_pct = 39.1,
    race_ethnicity = c(White = 51.9, Asian = 24.0, Other = 24.0),
    disease_state = "Relapsed or refractory B cell non-Hodgkin lymphoma: indolent B-NHL 39.9% (follicular lymphoma n = 158), large B cell lymphoma 47.9%, mantle cell lymphoma 12.2%",
    dose_range = "Subcutaneous 2-step step-up dosing to full doses of 0.0128-60 mg (424 of 466 patients on 48 mg: 0.16/0.8/48 mg, QW cycles 1-3, Q2W cycles 4-9, Q4W from cycle 10); 3-step step-up dosing 0.16/0.8/3 (or 6)/48 mg in the follicular lymphoma optimization cohort (n = 81)",
    regions = "Europe, Asia (all EPCORE NHL-3 patients Japanese), North America, Australia",
    notes = "EPCORE NHL-1 (NCT03625037; escalation n = 48, expansion n = 354, optimization n = 81) and EPCORE NHL-3 (NCT04542824; n = 64). The Results describe a 466-patient 2-step dataset; Methods state the model was updated with expansion and optimization (3-step) data and externally validated in both regimens. Demographics are pooled from supplementary Table S1 over the 2-step cohorts."
  )

  ini({
    # Structural parameters for a 75 kg, 65-year-old, BMI 25 kg/m^2 patient
    # with aggressive B-NHL and SUMPPD 30 cm^2 (supplementary Table S2)
    lcl <- log(0.526); label("Apparent nonspecific (linear) clearance CL/F (L/day)") # Table S2: CL/F = 0.526 L/day (RSE 2.15%)
    lq <- log(0.527); label("Apparent intercompartmental clearance Q/F (L/day)") # Table S2: Q/F = 0.527 L/day (RSE 4.91%)
    lvc <- log(9.42); label("Apparent central volume of distribution Vc/F (L)") # Table S2: Vc/F = 9.42 L (RSE 2.51%)
    lvp <- log(10.7); label("Apparent peripheral volume of distribution Vp/F (L)") # Table S2: Vp/F = 10.7 L (RSE 7.68%)
    lka <- log(0.534); label("First-order subcutaneous absorption rate constant ka (1/day)") # Table S2: ka = 0.534 1/day (RSE 3.21%)

    # QSS target-mediated disposition (constant total target; Table S2)
    lrbase_target <- log(2.37); label("Baseline (constant) total target concentration BASE (ug/mL)") # Table S2: BASE = 2.37 ug/mL (RSE 4.83%)
    lkss <- log(0.277); label("Quasi-steady-state constant KSS (ug/mL)") # Table S2: KSS = 0.277 ug/mL (RSE 6.16%)
    lkint <- log(0.0253); label("Internalization rate constant of the drug-target complex kint (1/day)") # Table S2: kint = 0.0253 1/day (RSE 9.13%)

    # Covariate effects (Table S2)
    e_wt_cl <- 0.949; label("Power exponent of (WT/75) on CL/F (unitless)") # Table S2: weight effect on CL/F = 0.949 (RSE 8.81%)
    e_wt_q <- fixed(0.75); label("Power exponent of (WT/75) on Q/F (unitless)") # Table S2: weight effect on Q/F = 0.75 (Fixed)
    e_wt_vc <- 0.641; label("Power exponent of (WT/75) on Vc/F (unitless)") # Table S2: weight effect on Vc/F = 0.641 (RSE 14%)
    e_wt_vp <- fixed(1); label("Power exponent of (WT/75) on Vp/F (unitless)") # Table S2: weight effect on Vp/F = 1 (Fixed)
    e_age_ka <- -0.381; label("Power exponent of (AGE/65) on ka (unitless)") # Table S2: age effect on ka = -0.381 (RSE 36.3%)
    e_bmi_ka <- -0.689; label("Power exponent of (BMI/25) on ka (unitless)") # Table S2: BMI effect on ka = -0.689 (RSE 22.5%); reference 25 kg/m^2 assumed
    e_tumsz_kint <- 0.312; label("Power exponent of (SUMPPD/30) on kint (unitless)") # Table S2: tumor size effect on kint = 0.312 (RSE 20.2%); reference 30 cm^2 assumed
    e_tumtp_inhl_kint <- 0.484; label("Multiplicative factor on kint for indolent NHL subtype (unitless ratio)") # Table S2: iNHL subtype effect on kint = 0.484 (RSE 13.4%)

    # IIV: not reported by Sinnollareddy 2026; carried from the Li 2025 parent
    # model (Li 2025 Table 2 CV%, omega^2 = log(1 + CV^2))
    etalcl ~ 0.063959 # Li 2025 Table 2: IIV CL/F CV 25.7% (carried forward)
    etalq ~ 0.56355 # Li 2025 Table 2: IIV Q/F CV 87.0% (carried forward)
    etalvc ~ 0.092893 # Li 2025 Table 2: IIV Vc/F CV 31.2% (carried forward)
    etalvp ~ 1.0615 # Li 2025 Table 2: IIV Vp/F CV 137.5% (carried forward)
    etalka ~ 0.26176 # Li 2025 Table 2: IIV ka CV 54.7% (carried forward)
    etalrbase_target ~ 0.32885 # Li 2025 Table 2: IIV BASE CV 62.4% (carried forward)
    etalkss ~ 0.54082 # Li 2025 Table 2: IIV KSS CV 84.7% (carried forward)
    etalkint ~ 0.45109 # Li 2025 Table 2: IIV kint CV 75.5% (carried forward)
    etaRUV ~ 0.043963 # Li 2025 Table 2: random effect on residual-error magnitude CV 21.2% (carried forward)

    # Residual error (Table S2): proportional part as CV, additive part as SD
    propSd <- 0.164; label("Proportional residual error (fraction)") # Table S2: residual error proportional part (CV) = 0.164
    addSd <- 0.0151; label("Additive residual error (ug/mL)") # Table S2: residual error additive part (SD) = 0.0151 ug/mL
  })
  model({
    # Individual parameters. Weight and age forms and references follow the
    # Li 2025 parent model (supplementary Fig. S1 there); BMI and SUMPPD use
    # the same power form, and indolent NHL multiplies kint by a ratio.
    cl <- exp(lcl + etalcl) * (WT / 75)^e_wt_cl
    q <- exp(lq + etalq) * (WT / 75)^e_wt_q
    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / 75)^e_wt_vp
    ka <- exp(lka + etalka) * (AGE / 65)^e_age_ka * (BMI / 25)^e_bmi_ka
    rbase_target <- exp(lrbase_target + etalrbase_target)
    kss <- exp(lkss + etalkss)
    kint <- exp(lkint + etalkint) * (TUMSZ / 30)^e_tumsz_kint *
      e_tumtp_inhl_kint^TUMTP_INHL

    # QSS free concentration, written in the rationalised form
    # C = 2 * KSS * Ctot / (sqrt(a^2 + 4 * KSS * Ctot) - a), a = Ctot - BASE - KSS,
    # which is algebraically identical to 0.5 * (a + sqrt(a^2 + 4 * KSS * Ctot))
    # but does not cancel to zero when Ctot is far below BASE + KSS.
    ctot <- central / vc
    qss_a <- ctot - rbase_target - kss
    cfree <- 2 * kss * ctot / (sqrt(qss_a^2 + 4 * kss * ctot) - qss_a)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl * cfree - q * cfree + q / vp * peripheral1 -
      kint * rbase_target * cfree * vc / (kss + cfree)
    d/dt(peripheral1) <- q * cfree - q / vp * peripheral1

    # Overall target saturation C / (KSS + C)
    target_saturation <- cfree / (kss + cfree)

    # Observation: free epcoritamab concentration (ug/mL = mg/L). The residual
    # SD sqrt(C^2 * propSd^2 + addSd^2) is scaled by exp(etaRUV) as in the
    # Li 2025 parent model.
    Cc <- cfree
    propSdi <- propSd * exp(etaRUV)
    addSdi <- addSd * exp(etaRUV)
    Cc ~ add(addSdi) + prop(propSdi) + combined2()
  })
}
