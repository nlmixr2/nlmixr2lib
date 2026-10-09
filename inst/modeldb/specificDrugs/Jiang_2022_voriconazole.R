Jiang_2022_voriconazole <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption for intravenous and oral voriconazole in Chinese adults with talaromycosis (Talaromyces marneffei infection) (Jiang 2022); C-reactive protein enters clearance as an exponential inflammation effect."
  reference <- "Jiang Z, Wei Y, Huang W, Li B, Zhou S, Liao L, Li T, Liang T, Yu X, Li X, Zhou C, Cao C, Liu T. Population pharmacokinetics of voriconazole and initial dosage optimization in patients with talaromycosis. Front Pharmacol. 2022;13:982981. doi:10.3389/fphar.2022.982981"
  vignette <- "Jiang_2022_voriconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "voriconazole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "voriconazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRP = list(
      description = "C-reactive protein concentration (clinical laboratory assay; standard vs high-sensitivity not stated), time-varying: Jiang 2022 related clearance to CRP measured within the same period as each concentration",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters CL as an exponential effect scaled by 43.6 mg/L (Jiang 2022 Eq. 6):",
        "exp(e_crp_cl * CRP / 43.6). This is a SCALED but not CENTERED form, so the",
        "tabulated typical CL of 4.34 L/h is the value at CRP = 0; at CRP = 43.6 mg/L the",
        "typical CL is 4.34 * exp(-0.135) = 3.79 L/h. The paper does not say what 43.6 mg/L is",
        "(Table 1 gives only per-site medians of 70.5 and 59.1 mg/L; the Methods say covariates",
        "were centered by their medians, so it is presumably the median over the analysed records).",
        "The scaled reading is the printed equation and reproduces the paper's own Monte Carlo",
        "target-attainment tables; see the model vignette.",
        "Cohort CRP: site 1 mean 93.9, SD 71.0, median 70.5, range 1.6-207.7 mg/L; site 2 mean",
        "61.3, SD 41.8, median 59.1, range 0.9-202 mg/L (Jiang 2022 Table 1)."
      ),
      source_name = "CRP"
    )
  )

  # Covariates Jiang 2022 screened during stepwise covariate modelling but did
  # not retain in the final model. Documentation only -- these names are
  # deliberately absent from model(). Height, hemoglobin, neutrophils, AST,
  # total protein, total and direct bilirubin and urea were also screened
  # (Methods 'Clinical data collection') and are recorded in population$notes.
  covariatesDataExcluded <- list(
    CYP2C19_IM = list(
      description = "CYP2C19 intermediate-metabolizer phenotype indicator",
      units = "(binary)",
      type = "binary",
      notes = "31 of 69 patients (*1/*2, *1/*3; Jiang 2022 Results). CYP2C19 phenotype had no significant effect on voriconazole PK; extensive-metabolizer status on V (dOFV 5.235) passed forward inclusion only and was removed in backward elimination."
    ),
    CYP2C19_PM = list(
      description = "CYP2C19 poor-metabolizer phenotype indicator",
      units = "(binary)",
      type = "binary",
      notes = "8 of 69 patients (*2/*2, *2/*3; Jiang 2022 Results). Not retained."
    ),
    ALB = list(
      description = "Serum albumin concentration",
      units = "g/L",
      type = "continuous",
      notes = "Site medians 26 and 24.7 g/L (Jiang 2022 Table 1). Correlated with the interindividual variability of CL but its forward-inclusion dOFV was below 3.84, so it was excluded (Discussion)."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Site medians 61 and 50 kg, range 38-87 kg (Jiang 2022 Table 1). Screened and not retained."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Site medians 57 and 30 years, range 20-69 years (Jiang 2022 Table 1). Screened and not retained."
    ),
    SEXF = list(
      description = "Sex indicator, 1 = female",
      units = "(binary)",
      type = "binary",
      notes = "14 of 69 patients (20.3%) female (Jiang 2022 Table 1). Screened and not retained."
    ),
    HIV_POS = list(
      description = "HIV-positive comorbidity indicator",
      units = "(binary)",
      type = "binary",
      notes = "34 of 69 patients (all from the Baise site) were newly diagnosed HIV positive with no antiretroviral history (Jiang 2022 Results). Screened and not retained."
    ),
    CONMED_PPI = list(
      description = "Concomitant proton-pump inhibitor use",
      units = "(binary)",
      type = "binary",
      notes = "32 of 69 patients received a PPI (omeprazole, pantoprazole, lansoprazole or rabeprazole), maximum 40 mg/day (Jiang 2022 Results, Table 1). No significant effect was found."
    ),
    CONMED_STEROID = list(
      description = "Concomitant systemic corticosteroid use",
      units = "(binary)",
      type = "binary",
      notes = "7 of 69 patients received a glucocorticoid (Jiang 2022 Table 1). Screened and not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Site medians 27.2 and 29 U/L, range 4.0-236 U/L (Jiang 2022 Table 1). Screened and not retained."
    ),
    GGT = list(
      description = "Gamma-glutamyltransferase",
      units = "U/L",
      type = "continuous",
      notes = "Site medians 290.1 and 95 U/L, range 19-1154 U/L (Jiang 2022 Table 1). Screened and not retained."
    ),
    PLT = list(
      description = "Platelet count",
      units = "10^9/L",
      type = "continuous",
      notes = "Site medians 417 and 102 x 10^9/L, range 6-625.7 (Jiang 2022 Table 1). Screened and not retained."
    ),
    WBC = list(
      description = "White blood cell count",
      units = "10^9/L",
      type = "continuous",
      notes = "Site medians 15.7 and 3.6 x 10^9/L, range 1.0-27.81 (Jiang 2022 Table 1). Screened and not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 69L,
    n_studies = 1L,
    n_observations = 233L,
    age_range = "20-69 years",
    age_median = "57 years (site 1), 30 years (site 2)",
    weight_range = "38-87 kg",
    weight_median = "61 kg (site 1), 50 kg (site 2)",
    sex_female_pct = 20.3,
    race_ethnicity = c(Chinese = 100),
    cyp2c19_phenotype = c(EM_pct = 43.5, IM_pct = 44.9, PM_pct = 11.6, RM_pct = 0, UM_pct = 0),
    co_medication = c(ProtonPumpInhibitor_pct = 46.4, Glucocorticoid_pct = 10.1),
    disease_state = "Adults (>= 18 years) with confirmed talaromycosis (Talaromyces marneffei infection) treated with voriconazole as initial therapy; 34 of 69 (49.3%) newly diagnosed HIV positive with no antiretroviral history. Excluded: Child-Pugh C hepatic dysfunction, creatinine > 3x upper limit of normal, other antifungals or interacting drugs, tuberculosis, chemotherapy, pregnancy.",
    dose_range = "Standard voriconazole dosing: IV 6 mg/kg q12h on day 1 (or oral 400 mg q12h) then 4 mg/kg q12h IV (or oral 200 mg q12h); non-loading regimen 4 mg/kg q12h IV or 200 mg q12h oral; oral dose halved below 40 kg. 50 patients IV loading + IV maintenance, 4 oral loading + oral maintenance, 1 IV loading + oral maintenance, 14 without loading dose (3 IV, 11 oral).",
    regions = "Guangxi, China: The First Affiliated Hospital of Guangxi Medical University, Nanning (n = 35) and People's Hospital of Baise (n = 34).",
    notes = paste(
      "Prospective observational study, February 2019 - November 2021. Sparse sampling: at least one",
      "sample within 30 min pre-dose and at 0.5, 1, 2, 4, 6, 8, 10 or 12 h post-dose; 233 concentrations",
      "(median 4 per patient, range 1-9) including 75 steady-state troughs. Concentrations by",
      "two-dimensional HPLC, LLOQ 0.2 mg/L. Observed troughs ranged 0.23-16.95 mg/L; 47.7% of 65",
      "troughs lay in the 1.0-5.5 mg/L target range, 12.3% below and 40.0% above.",
      "Covariates screened and not retained: sex, age, weight, height, HIV, PPIs, glucocorticoids,",
      "WBC, hemoglobin, platelets, neutrophils, ALT, AST, albumin, total protein, total bilirubin,",
      "GGT, urea and CYP2C19 phenotype. Demographics per Jiang 2022 Table 1; final-model estimates",
      "and 1000-replicate bootstrap (986 successful) per Jiang 2022 Table 2 and Eqs. 6-9."
    )
  )

  ini({
    # Absorption: ka fixed at 1.1/h from Pascual 2012 because few samples
    # fell in the absorption phase.
    lka <- fixed(log(1.1)); label("Absorption rate constant (1/h)")  # Jiang 2022 Methods 'Population pharmacokinetic model' ('Ka was fixed at 1.1 h-1 as reported in the literature'); Table 2 ka = 1.1 (Fixed); Eq. 8

    lcl <- log(4.34); label("Clearance at CRP = 0 mg/L (L/h)")  # Jiang 2022 Table 2 final model: CL = 4.34 L/h (RSE 18.6%, bootstrap median 4.39, 95% CI 2.86-4.46); Eq. 6

    lvc <- log(97.4); label("Volume of distribution (L)")  # Jiang 2022 Table 2 final model: V = 97.4 L (RSE 7.1%, bootstrap median 97.3, 95% CI 84.5-111.9); Eq. 7

    lfdepot <- log(0.951); label("Oral bioavailability (fraction)")  # Jiang 2022 Table 2 final model: F1 = 95.1% (RSE 20.5%, bootstrap median 93.7, 95% CI 46.4-134); Eq. 9

    e_crp_cl <- -0.135; label("Exponential CRP coefficient on CL, per unit of CRP/43.6 (unitless)")  # Jiang 2022 Table 2 final model: 'CRP on CL' = -0.135 (RSE 65.1%, bootstrap median -0.151, 95% CI -0.367 to 0.098); Eq. 6

    # IIV. Eqs. 6-7 print the omega^2 values in the exponent (e^1.01 and
    # e^0.0973) and Table 2 reports sqrt(omega^2) * 100: sqrt(1.01) = 100.5%,
    # sqrt(0.0973) = 31.2%.
    etalcl ~ 1.01  # Jiang 2022 Eq. 6 omega^2 = 1.01; Table 2 IIV_CL = 100.5% (RSE 11.7%, shrinkage 6.50%)
    etalvc ~ 0.0973  # Jiang 2022 Eq. 7 omega^2 = 0.0973; Table 2 IIV_V = 31.2% (RSE 14.7%, shrinkage 42.5%)

    # Residual error: combined model Y = F * (1 + eps1) + eps2 (Methods Eq. 5;
    # Results 'a combined model was used').
    propSd <- 0.071; label("Proportional residual error (fraction)")  # Jiang 2022 Table 2 final model: RSV_CV = 7.1% (RSE 16.5%, shrinkage 19.8%)
    addSd <- 0.373; label("Additive residual error (mg/L)")  # Jiang 2022 Table 2 final model: RSV_SD = 0.373 mg/L (RSE 13.5%, shrinkage 19.8%)
  })

  model({
    # Jiang 2022 Eqs. 6-9:
    #   CL (L/h) = 4.34 * exp(-0.135 * CRP(mg/L) / 43.6) * exp(eta_CL)
    #   V (L)    = 97.4 * exp(eta_V)
    #   Ka       = 1.1 /h (fixed);  F = 95.1%
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * exp(e_crp_cl * CRP / 43.6)
    vc <- exp(lvc + etalvc)

    # One-compartment disposition with first-order oral absorption.
    # Intravenous doses bypass the depot and enter central directly.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (cl / vc) * central

    # Oral bioavailability applies only to doses entering via the depot.
    f(depot) <- exp(lfdepot)

    # Amounts in mg and volume in L give mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
