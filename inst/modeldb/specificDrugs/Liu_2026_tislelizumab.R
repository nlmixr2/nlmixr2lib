Liu_2026_tislelizumab <- function() {
  description <- "Three-compartment population PK model for subcutaneous (thigh or abdomen) and intravenous tislelizumab (anti-PD-1 IgG4) in treatment-naive adults with advanced NSCLC (Liu 2026): two parallel first-order SC absorption depots with site-specific ka and absolute bioavailability, disposition and covariate effects fixed to the Budha 2023 intravenous model"
  reference <- paste(
    "Liu X, Guan Y, Liu Y, Yu T, Abdrashitov R, Niu Z, Lin X, Ye X, Wang Y,",
    "Sahasranaman S, Hanley WD, Budha N. Model-Informed Development of a",
    "Subcutaneous Formulation of Tislelizumab: Phase 1 Pharmacokinetics in",
    "Patients With Cancer. Clin Transl Sci. 2026;19(7):e70649.",
    "doi:10.1111/cts.70649.",
    "Disposition parameters, covariate effects, their IIV and residual error",
    "fixed from Budha et al. (2023) CPT Pharmacometrics Syst Pharmacol",
    "12(1):95-109, doi:10.1002/psp4.12880; see modellib('Budha_2023_tislelizumab').",
    sep = " "
  )
  vignette <- "Liu_2026_tislelizumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  # depot: subcutaneous abdomen injection (paper 'Depot 1');
  # depot2: subcutaneous thigh injection (paper 'Depot 2');
  # central: intravenous infusion.
  dosing <- c("depot", "depot2", "central")

  compartmentData <- list(
    depot = list(analyte = "tislelizumab", units = "mg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "tislelizumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tislelizumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "tislelizumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral2 = list(analyte = "tislelizumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL and Vc, reference 65 kg, fixed from the intravenous model (Liu 2026 Table 1; reference value from Budha 2023 Equations 5 and 6).",
      source_name = "WT"
    ),
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on Vc, reference 60 years, fixed from the intravenous model (Liu 2026 Table 1; reference value from Budha 2023 Equation 6).",
      source_name = "AGE"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL, reference 41 g/L, fixed from the intravenous model (Liu 2026 Table 1; reference value from Budha 2023 Equation 5). Liu 2026 Table S2 also reports albumin in g/L.",
      source_name = "ALB"
    ),
    TUMSZ = list(
      description = "Baseline tumor size (sum of target-lesion diameters)",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = "Liu 2026 calls it 'baseline tumor burden' (Table 1, Table S2, mm); Budha 2023 calls the same covariate TUMSZ. Power effect on CL, reference 63 mm (Budha 2023 Equation 5).",
      source_name = "baseline tumor burden"
    ),
    SEXF = list(
      description = "Biological sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Exponential effect on Vc for females, fixed from the intravenous model (Liu 2026 Table 1 'Effect of sex on Vc/F'; Budha 2023 Equation 6, SEX = Female).",
      source_name = "SEX"
    ),
    ADA_POS = list(
      description = "Treatment-emergent antidrug-antibody status",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ADA-negative)",
      notes = "Exponential effect on CL, fixed from the intravenous model (Liu 2026 Table 1 'Effect of treatment-emergent ADA status on CL/F').",
      source_name = "ADA"
    ),
    TUMTP_GASTRIC = list(
      description = "Tumor-type indicator for gastric cancer",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all other tumor types)",
      notes = "Exponential effect on CL, fixed from the intravenous model (Liu 2026 Table 1 'Effect of gastric cancer status on CL/F'). 0 for the NSCLC patients the absorption model was fitted to; 1 for the first-line gastric-cancer population of the Liu 2026 phase 3 dose-selection simulations.",
      source_name = "TUMTP"
    ),
    TUMTP_HODGKIN_CLASSICAL = list(
      description = "Tumor-type indicator for classical Hodgkin lymphoma",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all other tumor types)",
      notes = "Exponential effect on CL, fixed from the intravenous model (Liu 2026 Table 1 'Effect of classic Hodgkin lymphoma status on CL/F'). 0 for every patient in Liu 2026.",
      source_name = "TUMTP"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 60,
    n_studies = 1,
    n_observations = 922,
    disease_state = "Treatment-naive locally advanced or metastatic non-small cell lung cancer (ECOG PS 0-1)",
    dose_range = "Tislelizumab 300 mg SC (thigh or abdomen) or 200 mg IV in Cycles 1-3 of Part 1 (six randomized site sequences, Table S1), then 200 mg IV Q3W from Cycle 4; Part 2 300 mg SC thigh Q3W throughout",
    regions = "China and Europe (15 sites)",
    notes = paste(
      "Phase 1 study BGB-A317-103 (NCT06091943): 39 Part 1 and 21 Part 2 patients",
      "with evaluable PK (data cutoff 31 July 2025). Demographics are published",
      "in the primary clinical report (Liu 2026 reference 25), not in this paper.",
      "No covariate was retained on the absorption parameters; the disposition",
      "covariates were carried from the Budha 2023 intravenous model (2596",
      "patients, 12 studies)."
    )
  )

  ini({
    # Disposition parameters fixed to the intravenous model estimates, at the
    # reference patient (65 kg, 60 y, ALB 41 g/L, tumor size 63 mm, male,
    # ADA-negative, not gastric cancer or cHL).
    lcl <- fixed(log(0.153)); label("Clearance (L/day)") # Liu 2026 Table 1 'CL/F, L/d' 0.153 Fixed
    lvc <- fixed(log(3.05)); label("Central volume of distribution (L)") # Liu 2026 Table 1 'Vc/F, L' 3.05 Fixed
    lq <- fixed(log(0.74)); label("Intercompartmental clearance, central-peripheral1 (L/day)") # Liu 2026 Table 1 'Q2/F, L/d' 0.74 Fixed
    lvp <- fixed(log(1.27)); label("Peripheral1 volume of distribution (L)") # Liu 2026 Table 1 'V2/F, L' 1.27 Fixed
    lq2 <- fixed(log(0.092)); label("Intercompartmental clearance, central-peripheral2 (L/day)") # Liu 2026 Table 1 'Q3/F, L/d' 0.092 Fixed
    lvp2 <- fixed(log(2.1)); label("Peripheral2 volume of distribution (L)") # Liu 2026 Table 1 'V3/F, L' 2.1 Fixed

    e_wt_cl <- fixed(0.565); label("Power exponent of WT on CL (unitless)") # Liu 2026 Table 1 'Effect of baseline weight on CL/F' 0.565 Fixed
    e_alb_cl <- fixed(-0.457); label("Power exponent of ALB on CL (unitless)") # Liu 2026 Table 1 'Effect of baseline albumin on CL/F' -0.457 Fixed
    e_tumsz_cl <- fixed(0.0735); label("Power exponent of TUMSZ on CL (unitless)") # Liu 2026 Table 1 'Effect of baseline tumor burden on CL/F' 0.0735 Fixed
    e_ada_cl <- fixed(0.111); label("Exponential coefficient of ADA-positive status on CL (unitless)") # Liu 2026 Table 1 'Effect of treatment-emergent ADA status on CL/F' 0.111 Fixed
    e_tumtp_gc_cl <- fixed(0.069); label("Exponential coefficient of gastric cancer on CL (unitless)") # Liu 2026 Table 1 'Effect of gastric cancer status on CL/F' 0.069 Fixed
    e_tumtp_chl_cl <- fixed(-0.216); label("Exponential coefficient of classical Hodgkin lymphoma on CL (unitless)") # Liu 2026 Table 1 'Effect of classic Hodgkin lymphoma status on CL/F' -0.216 Fixed
    e_wt_vc <- fixed(0.397); label("Power exponent of WT on Vc (unitless)") # Liu 2026 Table 1 'Effect of baseline body weight on Vc/F' 0.397 Fixed
    e_sexf_vc <- fixed(-0.116); label("Exponential coefficient of female sex on Vc (unitless)") # Liu 2026 Table 1 'Effect of sex on Vc/F' -0.116 Fixed
    e_age_vc <- fixed(0.0966); label("Power exponent of AGE on Vc (unitless)") # Liu 2026 Table 1 'Effect of age on Vc/F' 0.0966 Fixed

    # Subcutaneous absorption, estimated by injection site.
    lka_abdomen <- log(0.247); label("Absorption rate constant, abdomen SC depot (1/day)") # Liu 2026 Table 1 'Ka abdomen, 1/d' 0.247 (RSE 10.9%)
    lka_thigh <- log(0.189); label("Absorption rate constant, thigh SC depot (1/day)") # Liu 2026 Table 1 'Ka thigh, 1/d' 0.189 (RSE 6.1%)
    lfdepot_abdomen <- log(0.731); label("Absolute bioavailability, abdomen SC injection (fraction)") # Liu 2026 Table 1 'Bioavailability, abdomen' 0.731 (RSE 5.9%)
    lfdepot_thigh <- log(0.866); label("Absolute bioavailability, thigh SC injection (fraction)") # Liu 2026 Table 1 'Bioavailability, thigh' 0.866 (RSE 3.2%)

    # IIV. Liu 2026 Table 1 reports omega on the SD scale in percent, so each
    # variance is (value / 100)^2. The CL-Vc covariance is not printed in
    # Liu 2026; it is the intravenous-model value the disposition block was
    # fixed to.
    etalcl + etalvc ~ fixed(c(0.069169, 0.020, 0.027889)) # Liu 2026 Table 1 omega CL/F 26.3% and Vc/F 16.7% (intravenous-model values); covariance 0.020 from Budha 2023 Table 3
    etalvp ~ fixed(0.558009) # Liu 2026 Table 1 omega V2/F 74.7% (intravenous-model value)
    etalvp2 ~ fixed(0.998001) # Liu 2026 Table 1 omega V3/F 99.9% (intravenous-model value)
    etalka_abdomen ~ 0.322624 # Liu 2026 Table 1 omega Ka abdomen 56.8% (shrinkage 34.5%)
    etalka_thigh ~ 0.145924 # Liu 2026 Table 1 omega Ka thigh 38.2% (shrinkage 17.1%)
    etalfdepot_abdomen ~ 0.069696 # Liu 2026 Table 1 omega bioavailability abdomen 26.4% (shrinkage 36.3%)
    etalfdepot_thigh ~ 0.018769 # Liu 2026 Table 1 omega bioavailability thigh 13.7% (shrinkage 38.9%)

    propSd <- fixed(0.126); label("Proportional residual error (fraction)") # Liu 2026 Table 1 'Proportional error' 12.6% Fixed
    addSd <- fixed(2.09); label("Additive residual error (mg/L)") # Liu 2026 Table 1 'Additive error' 2.09 Fixed
  })
  model({
    # Disposition: power effects of continuous covariates about their
    # reference values and exponential effects of indicators (Budha 2023
    # Equations 5 and 6), all fixed in Liu 2026.
    cl <- exp(lcl + etalcl) *
      (WT / 65)^e_wt_cl *
      (ALB / 41)^e_alb_cl *
      (TUMSZ / 63)^e_tumsz_cl *
      exp(e_ada_cl * ADA_POS) *
      exp(e_tumtp_gc_cl * TUMTP_GASTRIC) *
      exp(e_tumtp_chl_cl * TUMTP_HODGKIN_CLASSICAL)
    vc <- exp(lvc + etalvc) *
      (WT / 65)^e_wt_vc *
      (AGE / 60)^e_age_vc *
      exp(e_sexf_vc * SEXF)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)
    q2 <- exp(lq2)
    vp2 <- exp(lvp2 + etalvp2)

    # Site-specific SC absorption (Liu 2026 Figure S1: abdomen = Depot 1,
    # thigh = Depot 2, each with its own ka and F). Exponential IIV on F as on
    # every other parameter; F is not bounded at 1.
    ka_abdomen <- exp(lka_abdomen + etalka_abdomen)
    ka_thigh <- exp(lka_thigh + etalka_thigh)
    fdepot_abdomen <- exp(lfdepot_abdomen + etalfdepot_abdomen)
    fdepot_thigh <- exp(lfdepot_thigh + etalfdepot_thigh)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(depot) <- -ka_abdomen * depot
    d/dt(depot2) <- -ka_thigh * depot2
    d/dt(central) <- ka_abdomen * depot + ka_thigh * depot2 -
      kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    f(depot) <- fdepot_abdomen
    f(depot2) <- fdepot_thigh

    # mg / L
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
