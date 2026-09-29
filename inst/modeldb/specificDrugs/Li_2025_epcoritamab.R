Li_2025_epcoritamab <- function() {
  description <- "Two-compartment quasi-steady-state (QSS) target-mediated drug disposition population PK model with first-order subcutaneous absorption and constant total target for epcoritamab (CD3xCD20 bispecific antibody) in adults with relapsed or refractory B cell non-Hodgkin lymphoma (Li 2025). Body weight (reference 75 kg) scales CL/F, Q/F, Vc/F and Vp/F by power laws, age (reference 65 years) lowers ka, and the combined additive + proportional residual error carries a per-subject random effect on its magnitude."
  reference <- "Li T, Gibiansky L, Parikh A, van der Linden M, Sanghavi K, Putnins M, Sacchi M, Feng H, Ahmadi T, Gupta M, Xu S. Population Pharmacokinetics of Epcoritamab Following Subcutaneous Administration in Relapsed or Refractory B Cell Non-Hodgkin Lymphoma. Clin Pharmacokinet. 2025;64(1):127-141. doi:10.1007/s40262-024-01464-2"
  vignette <- "Li_2025_epcoritamab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # etaRUV is the random effect on the magnitude of the residual error
  # (supplementary Fig. S1, eta9), which has no canonical IIV name.
  paper_specific_etas <- c("etaRUV")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline (time-fixed) in the source analysis (supplementary Fig. S1 legend: 'WT baseline weight (kg)'). Power covariate on CL/F (estimated exponent), Q/F (fixed 0.75), Vc/F (estimated exponent) and Vp/F (fixed 1), all normalised to 75 kg.",
      source_name = "WT"
    ),
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline in the source analysis. Power covariate on ka normalised to 65 years (supplementary Fig. S1: COVka = (AGE/65)^theta15).",
      source_name = "AGE"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex, 1 = female",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened (Section 2.4.2) but not retained in the final model (Section 3.4)."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = "Screened (Section 2.4.2) but not retained in the final model (Section 3.4); lower apparent clearance and volumes in EPCORE NHL-3 were attributed to lower body weight."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Section 2.4.2) but not retained in the final model (Section 3.4)."
    ),
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Renal function screened (Section 2.4.2; Cockcroft-Gault per supplementary Fig. S1) but not retained in the final model (Section 3.4)."
    ),
    ADA_POS = list(
      description = "Antidrug-antibody positive status, 1 = positive",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ADA-negative)",
      notes = "Screened (Section 2.4.2) but not retained; Section 3.4 reports no apparent effect of ADA on exposure (14 ADA-positive of 327)."
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
    n_subjects = 327,
    n_studies = 2,
    n_observations = 6819,
    age_range = "20-89 years",
    age_median = "67 years",
    weight_range = "39-144 kg",
    weight_median = "70 kg",
    sex_female_pct = 40.4,
    race_ethnicity = c(White = 56.6, Asian = 30.3, `Native American` = 0.3, Other = 12.8),
    disease_state = "Relapsed or refractory B cell non-Hodgkin lymphoma: large B cell lymphoma (64.8%), indolent B-NHL (29.7%), mantle cell lymphoma (5.5%)",
    dose_range = "Subcutaneous step-up dosing then full doses of 0.0128-60 mg (overall 0.004-60 mg); 298 patients on the approved 0.16/0.8/48 mg regimen (48 mg QW cycles 1-3, Q2W cycles 4-9, Q4W from cycle 10; 28-day cycles)",
    regions = "Europe (52.0%), Asia (29.4%), Australia (9.5%), North America (9.2%)",
    renal_function = "CrCl >= 90 mL/min 33.6%, 60-<90 mL/min 45.0%, 30-<60 mL/min 20.8%, missing 0.6%",
    hepatic_function = "Normal 82.6%, mild 16.2%, moderate 0.3%, missing 0.9%",
    notes = "EPCORE NHL-1 (NCT03625037; dose escalation n = 35, expansion n = 232) and EPCORE NHL-3 (NCT04542824; Japan, n = 60). Baseline characteristics from Table 1; trial design from supplementary Table S1."
  )

  ini({
    # Structural parameters for a 75 kg, 65-year-old patient (Table 2)
    lcl <- log(0.481); label("Apparent nonspecific (linear) clearance CL/F (L/day)") # Table 2: CL/F = 0.481 L/day (RSE 2.66%)
    lq <- log(0.488); label("Apparent intercompartmental clearance Q/F (L/day)") # Table 2: Q/F = 0.488 L/day (RSE 7.88%)
    lvc <- log(9.33); label("Apparent central volume of distribution Vc/F (L)") # Table 2: Vc/F = 9.33 L (RSE 3.18%)
    lvp <- log(14.1); label("Apparent peripheral volume of distribution Vp/F (L)") # Table 2: Vp/F = 14.1 L (RSE 10.9%)
    lka <- log(0.584); label("First-order subcutaneous absorption rate constant ka (1/day)") # Table 2: ka = 0.584 1/day (RSE 4.94%)

    # QSS target-mediated disposition (constant total target; Table 2)
    lrbase_target <- log(2.03); label("Baseline (constant) total target concentration BASE (ug/mL)") # Table 2: BASE = 2.03 ug/mL (RSE 5.84%)
    lkss <- log(0.214); label("Quasi-steady-state constant KSS (ug/mL)") # Table 2: KSS = 0.214 ug/mL (RSE 7.27%)
    lkint <- log(0.0278); label("Internalization rate constant of the drug-target complex kint (1/day)") # Table 2: kint = 0.0278 1/day (RSE 9.31%)

    # Covariate effects (Table 2; supplementary Fig. S1 covariate model)
    e_wt_cl <- 0.875; label("Power exponent of (WT/75) on CL/F (unitless)") # Table 2: effect of weight on CL/F = 0.875 (RSE 10.7%)
    e_wt_q <- fixed(0.75); label("Power exponent of (WT/75) on Q/F (unitless)") # Table 2: effect of weight on Q/F = 0.75 (fixed)
    e_wt_vc <- 0.603; label("Power exponent of (WT/75) on Vc/F (unitless)") # Table 2: effect of weight on Vc/F = 0.603 (RSE 16.6%)
    e_wt_vp <- fixed(1); label("Power exponent of (WT/75) on Vp/F (unitless)") # Table 2: effect of weight on Vp/F = 1 (fixed)
    e_age_ka <- -0.503; label("Power exponent of (AGE/65) on ka (unitless)") # Table 2: effect of age on ka = -0.503 (RSE 33.0%)

    # IIV: Table 2 reports CV%; omega^2 = log(1 + CV^2)
    etalcl ~ 0.063959 # Table 2: IIV CL/F CV 25.7%
    etalq ~ 0.56355 # Table 2: IIV Q/F CV 87.0%
    etalvc ~ 0.092893 # Table 2: IIV Vc/F CV 31.2%
    etalvp ~ 1.0615 # Table 2: IIV Vp/F CV 137.5%
    etalka ~ 0.26176 # Table 2: IIV ka CV 54.7%
    etalrbase_target ~ 0.32885 # Table 2: IIV BASE CV 62.4%
    etalkss ~ 0.54082 # Table 2: IIV KSS CV 84.7%
    etalkint ~ 0.45109 # Table 2: IIV kint CV 75.5%
    etaRUV ~ 0.043963 # Table 2: random effect on residual-error magnitude CV 21.2% (footnote a)

    # Residual error (Table 2; supplementary Fig. S1). sigma^2 = 1 is fixed,
    # so the THETAs below are the standard deviations themselves.
    propSd <- 0.189; label("Proportional residual error (fraction)") # Table 2: residual error proportional (CV) = 0.189
    addSd <- 0.0133; label("Additive residual error (ug/mL)") # Table 2: residual error additive (SD) = 0.0133 ug/mL
  })
  model({
    # Individual parameters (supplementary Fig. S1). The printed Fig. S1
    # equation for Q omits exp(eta2), but Table 2 reports an estimated IIV on
    # Q/F (CV 87.0%, shrinkage 30.0%) and Fig. S1 / S3 define eta2 as the
    # random effect of Q, so it is included here.
    cl <- exp(lcl + etalcl) * (WT / 75)^e_wt_cl
    q <- exp(lq + etalq) * (WT / 75)^e_wt_q
    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / 75)^e_wt_vp
    ka <- exp(lka + etalka) * (AGE / 65)^e_age_ka
    rbase_target <- exp(lrbase_target + etalrbase_target)
    kss <- exp(lkss + etalkss)
    kint <- exp(lkint + etalkint)

    # QSS free concentration (supplementary Fig. S1). central holds the total
    # (free + target-bound) drug amount; the total target is constant at BASE.
    # Fig. S1 prints C = 0.5 * (a + sqrt(a^2 + 4 * KSS * Ctot)) with
    # a = Ctot - BASE - KSS. That form cancels catastrophically when Ctot is
    # far below BASE + KSS (C rounds to exactly 0 and elimination stops), so
    # the algebraically identical rationalised form is used:
    # C = 2 * KSS * Ctot / (sqrt(a^2 + 4 * KSS * Ctot) - a).
    ctot <- central / vc
    qss_a <- ctot - rbase_target - kss
    cfree <- 2 * kss * ctot / (sqrt(qss_a^2 + 4 * kss * ctot) - qss_a)

    # ODEs (supplementary Fig. S1: A1 = depot, A2 = central, A3 = peripheral1)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl * cfree - q * cfree + q / vp * peripheral1 -
      kint * rbase_target * cfree * vc / (kss + cfree)
    d/dt(peripheral1) <- q * cfree - q / vp * peripheral1

    # Derived quantities shown in Figure 3: overall target saturation
    # C / (KSS + C) and concentration-dependent apparent total clearance
    # CLtot/F = CL/F + kint * BASE * (Vc/F) / (KSS + C).
    target_saturation <- cfree / (kss + cfree)
    cltot <- cl + kint * rbase_target * vc / (kss + cfree)

    # Observation: free epcoritamab concentration (ug/mL = mg/L). Residual
    # SD = sqrt(C^2 * propSd^2 + addSd^2) * exp(etaRUV) (supplementary Fig. S1);
    # scaling both components by exp(etaRUV) is exactly that form.
    Cc <- cfree
    propSdi <- propSd * exp(etaRUV)
    addSdi <- addSd * exp(etaRUV)
    Cc ~ add(addSdi) + prop(propSdi) + combined2()
  })
}
