Maggs_2026_fusidicAcid <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption for",
    "oral fusidic acid (500 mg three times daily) given with oral rifampicin",
    "(300 mg twice daily) to adults with staphylococcal bone, joint or",
    "endovascular infection (Maggs 2026). Apparent clearance and volume are",
    "allometrically scaled to 70 kg body weight (exponents fixed at 0.75 and",
    "1). Apparent clearance carries two time-dependent step factors: an",
    "autoinhibition step from 24 h after the first fusidic acid dose",
    "(CL/F x (1 - 0.756)) and a rifampicin-mediated induction step from 7",
    "days after the first dose (CL/F x (1 + 1.33)), which act together, so",
    "the net steady-state clearance is 0.57 times the Day-1 value.",
    "Inter-individual variability on ka and CL/F; proportional residual",
    "error on log-transformed concentrations."
  )
  reference <- paste(
    "Maggs C, Manning L, Yoo O, Salman S, Page-Sharp M, Aboltins C. (2026).",
    "A prospective pharmacokinetic interaction study between rifampicin and",
    "fusidic acid for the treatment of staphylococcal infections.",
    "J Antimicrob Chemother 81(10):959-967.",
    "doi:10.1093/jac/dkag293."
  )
  vignette <- "Maggs_2026_rifampicin_fusidicAcid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "fusidic acid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "fusidic acid", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight, used for a priori allometric scaling of CL/F and V/F.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference 70 kg (Table 2 units 'L h-1.70 kg-1' and footnote 'for an",
        "individual with body weight of 70 kg'); exponents 0.75 on clearance",
        "and 1 on volume fixed a priori (Methods, Pharmacokinetic analysis).",
        "The cohort's body weights are not reported (median BMI 31 kg/m2)."
      ),
      source_name = "body weight"
    ),
    T_FIRSTDOSE = list(
      description = "Time since the first fusidic acid dose of the treatment course (h).",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives two step factors on CL/F: autoinhibition '(1 + INT) x CL/F",
        "after 24 hours from the first dose' (INT = -0.756) and rifampicin",
        "induction '(1 + rIND) x CL/F after 7 days from the first dose' (rIND",
        "= 1.33) (Table 2). The two factors multiply: 1 / (1 - 0.756) = 4.10",
        "and 1 / ((1 - 0.756) x 2.33) = 1.76 reproduce the Day-2 (4.09-fold)",
        "and steady-state (1.76-fold) AUCtau ratios of the Results. The paper",
        "does not say whether 'the first dose' is the first fusidic acid dose",
        "or the first study dose (which was rifampicin for 6 of 10",
        "participants, 24 h earlier); the first fusidic acid dose is used",
        "here because the autoinhibition is driven by fusidic acid itself,",
        "and Figure 1 agrees: the rifampicin-index participants' Day-2",
        "profiles (their first fusidic acid dose) are low, as expected",
        "before autoinhibition starts.",
        "Rifampicin co-administration is not itself a covariate: every",
        "participant was on rifampicin from the second study day, so the",
        "induction step is purely time-based, as in the source."
      ),
      source_name = "time after the first dose"
    )
  )

  covariatesDataExcluded <- list(
    CRP = list(
      description = "C-reactive protein.",
      units = "mg/L",
      type = "continuous",
      notes = "Screened on ka, CL and V (Methods; Table S1 runs 3250-3264); not retained. The paper prints the baseline median as 44 (IQR 22-112) mg/dL, which is implausible for CRP and is most likely mg/L."
    ),
    ALB = list(
      description = "Serum albumin.",
      units = "g/L",
      type = "continuous",
      notes = "Screened on ka and CL (Methods); not retained."
    ),
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      notes = "Screened on ka, CL and V (Methods; Table S1); not retained."
    ),
    HCT = list(
      description = "Haematocrit.",
      units = "(fraction)",
      type = "continuous",
      notes = "Screened on ka (Methods); not retained."
    ),
    CREAT = list(
      description = "Serum creatinine.",
      units = "umol/L",
      type = "continuous",
      notes = "Screened on CL (Table S1 runs 3250-3256); not retained."
    ),
    HGB = list(
      description = "Haemoglobin.",
      units = "g/L",
      type = "continuous",
      notes = "Screened on CL (Table S1 runs 3250-3256); not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 10L,
    n_studies = 1L,
    age_range = "median 59 years (IQR 48-68)",
    weight_range = "not reported; median BMI 31 kg/m2 (IQR 29-35)",
    sex_female_pct = 30,
    race_ethnicity = "not reported",
    disease_state = paste(
      "Presumed or confirmed staphylococcal infection: osteoarticular infection",
      "of the axial skeleton (4), large-joint septic arthritis (4, including 2",
      "prosthetic joint infections), endovascular infection (1) and small-joint",
      "infection (1). Baseline median creatinine 78 umol/L (IQR 63-89), albumin",
      "27 g/L (21-35), CRP 44 (22-112)."
    ),
    dose_range = paste(
      "Fusidic acid 500 mg orally three times daily with rifampicin 300 mg",
      "orally twice daily. Randomised index drug for the first 24 h (fusidic",
      "acid n = 4, rifampicin n = 6), then the combination for the treatment",
      "course."
    ),
    regions = "Melbourne, Victoria, Australia (St Vincent's Hospital Melbourne and Northern Health), 2019-2022.",
    n_observations = "96 fusidic acid concentrations.",
    notes = paste(
      "Maggs 2026 Results paragraph 1. Concentrations were measured by",
      "LC-MS/MS in capillary dried blood spots (LLOQ 0.15 mg/L) on Day 1, Day",
      "2 and Day 28 +/- 7 (pre-dose, 30-40 min, 1.5-2 h and 3.5-4.5 h, plus",
      "7-8 h on Days 1 and 2; Table 1). NONMEM 7.5.1, FOCE-I."
    )
  )

  ini({
    # Maggs 2026 Table 2, final fusidic acid model (Table S1 run 1654, OFV
    # 63.799). Typical values are for a 70 kg individual.
    lka <- log(0.419)
    label("Absorption rate constant ka (1/h)") # Table 2 'Ka, h-1' = 0.419 (RSE 47%)
    lcl <- log(6.9)
    label("Apparent clearance CL/F at 70 kg before autoinhibition (L/h)") # Table 2 'CL/F, L h-1.70 kg-1' = 6.9 (RSE 32%)
    lvc <- log(45.3)
    label("Apparent volume V/F at 70 kg (L)") # Table 2 'V/F, l.70 kg-1' = 45.3 (RSE 25%)

    e_autoinh_cl <- -0.756
    label("Fractional change in CL/F from 24 h after the first dose, autoinhibition (fraction)") # Table 2 'Inhibition (1 + INT) x CL/F after 24 hours from the first dose' = -0.756 (RSE 17%)
    e_rif_ind_cl <- 1.33
    label("Fractional change in CL/F from 7 days after the first dose, rifampicin induction (fraction)") # Table 2 'Rifampicin induction (1 + rIND) x CL/F after 7 days from the first dose' = 1.33 (RSE 81%)

    # IIV and residual error are printed as 100% x sqrt(variance estimate)
    # (Table 2 footnote), so each variance is the squared fraction.
    etalka ~ 0.962361 # Table 2 'IIV in Ka' = 98.1% (shrinkage 21%); 0.981^2
    etalcl ~ 0.605284 # Table 2 'IIV in CL/F' = 77.8% (shrinkage 9%); 0.778^2

    propSd <- 0.722
    label("Proportional residual error SD (fraction)") # Table 2 'Proportional error (%)' = 72.2 (shrinkage 8%)
  })
  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 70)^0.75 *
      (1 + e_autoinh_cl * (T_FIRSTDOSE >= 24)) *
      (1 + e_rif_ind_cl * (T_FIRSTDOSE >= 168))
    vc <- exp(lvc) * (WT / 70)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
