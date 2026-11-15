Jovanovic_2026_adalimumab <- function() {
  description <- paste(
    "One-compartment population PK model with first-order subcutaneous",
    "absorption and linear elimination for adalimumab in adults with Crohn's",
    "disease, developed on sparse therapeutic-drug-monitoring trough",
    "concentrations from a Serbian real-world cohort with NONMEM $PRIOR",
    "informative priors from the Ternant 2015 model (on V/F, IIV of V/F and",
    "the anti-drug-antibody effect). Apparent clearance carries a",
    "proportional anti-drug-antibody effect, a log-linear C-reactive-protein",
    "effect and a proportional increase for weekly dosing (Jovanovic 2026)."
  )
  reference <- paste(
    "Jovanovic M, Topic Vucenovic V, Roganovic M, Pavlovic G, Kralj D,",
    "Markovic S, Svorcan P, Vucicevic K. Unlocking the Potential of",
    "Population Pharmacokinetic Models of Adalimumab in Patients with",
    "Crohn's Disease. Pharmaceutics. 2026;18(7):788.",
    "doi:10.3390/pharmaceutics18070788.",
    "Final-model parameter values from Table 4; clearance covariate model",
    "from Equation (6)."
  )
  vignette <- "Jovanovic_2026_adalimumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(
      analyte = "adalimumab",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "adalimumab",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    )
  )

  covariateData <- list(
    ADA_POS = list(
      description = "Anti-adalimumab antibody status, 1 = present, 0 = absent",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ADA absent)",
      notes = paste(
        "Proportional effect on CL/F, (1 + 6.9 * ADA_POS), i.e. a 7.9-fold",
        "higher CL/F when ADA are present (Equation (6); Table 4 row 'ADA on",
        "CL/F' = 6.90). Retained a priori and estimated with an informative",
        "prior from Ternant 2015 because only 9 patients were ADA-positive.",
        "ADA were measured on 294 occasions with the R-Biopharm ELISA and",
        "categorised as present or absent; the status is time-varying.",
        "Missing ADA status was imputed as absent (0) in the source analysis."
      ),
      source_name = "ADA"
    ),
    CRP = list(
      description = "Serum C-reactive protein, standard assay, time-varying",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Log-linear effect on CL/F, (1 + 0.075 * log(CRP)), with CRP in mg/L",
        "and NO centring value (Equation (6); Table 4 row 'ln(CRP) on CL/F' =",
        "0.075). The unit multiplier is therefore reached at CRP = 1 mg/L, not",
        "at the cohort mean (8.46 mg/L, Table 2). The form is undefined at",
        "CRP = 0 and the paper does not state how zero values were floored;",
        "Table 2 reports a minimum of 0.00 mg/L while the lowest ln(CRP) in",
        "Figure S1 is about -4.6 (about 0.01 mg/L). Supply CRP > 0. Missing",
        "values were imputed with the patient's own median, else the",
        "population median, in the source analysis."
      ),
      source_name = "CRP"
    ),
    REGI_QW = list(
      description = "1 = adalimumab dosed every week, 0 = any other interval",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (standard every-other-week dosing)",
      notes = paste(
        "Proportional effect on CL/F, (1 + 0.218 * REGI_QW): CL/F is 21.8%",
        "higher on weekly administration than on the standard every-other-week",
        "schedule (Results 3.5; Equation (6) term 'REGIMEN'; Table 4 row",
        "'Dosing regimen on CL/F' = 0.218). The covariate screened was",
        "'dosing interval (one week or otherwise)' (Methods 2.4). The authors",
        "interpret it as a marker of treatment intensification in response to",
        "disease activity rather than a causal effect of dosing frequency",
        "(Discussion), so it should not be used to predict the effect of",
        "switching an individual patient to weekly dosing."
      ),
      source_name = "REGIMEN"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 195L,
    n_studies = 1L,
    age_range = "17-67 years",
    age_mean = "38.05 years (SD 11.83)",
    weight_range = "38-105 kg",
    weight_mean = "69.44 kg (SD 13.47)",
    sex_female_pct = 50.8,
    race_ethnicity = "Not reported (single-centre Serbian cohort).",
    disease_state = "Adults with Crohn's disease on subcutaneous adalimumab undergoing therapeutic drug monitoring during induction or maintenance.",
    dose_range = paste(
      "Adalimumab (Humira 40 mg pre-filled pen) subcutaneously per the SPC:",
      "induction 160/80 mg or 80/40 mg at weeks 0/2, then 40 mg every other",
      "week, with optimisation (weekly dosing or 80 mg every other week)",
      "when clinically indicated."
    ),
    regions = "University Medical Center 'Zvezdara', Belgrade, Serbia",
    crp_mean = "8.46 mg/L (SD 13.47, range 0.00-72.70)",
    ada_positive = "9 patients (3.06% of 294 ADA measurements)",
    azathioprine_pct = 60.5,
    prior_anti_tnf_pct = 33.8,
    notes = paste(
      "Baseline demographics from Jovanovic 2026 Table 2. 942 trough",
      "concentrations (median time after dose 6.75 or 13.75 days; 4.46 +/-",
      "2.92 per patient, range 1-14), of which 15 were below the limit of",
      "quantification (M5 method) and 70 were reported qualitatively as",
      "> 12 mg/L (handled as right-censored, M3 method). Concentrations and",
      "ADA measured with R-Biopharm ELISA. NONMEM 7.6, Laplace conditional",
      "estimation with interaction, with $PRIOR NWPRI informative priors on",
      "V/F, IIV of V/F and the ADA effect from Ternant 2015 (CL/F and its",
      "IIV estimated without priors)."
    )
  )

  ini({
    lka <- fixed(log(0.15)); label("First-order SC absorption rate constant ka (1/day), Ternant 2015 value") # Table 4 row 'Ka (day-1)' = 0.15 (fix); Results 3.5
    lcl <- log(0.334); label("Apparent clearance CL/F for an ADA-negative patient on every-other-week dosing with CRP = 1 mg/L (L/day)") # Table 4 row 'CL/F (L/day)' Final model = 0.334 (RSE 3.4); Equation (6)
    lvc <- log(13.5); label("Apparent volume of distribution V/F (L)") # Table 4 row 'V/F (L)' Final model = 13.5 (RSE 6.6), informative prior from Ternant 2015

    e_ada_cl <- 6.90; label("Proportional ADA effect on CL/F, applied as (1 + e_ada_cl * ADA_POS) (fraction)") # Table 4 row 'ADA on CL/F' Final model = 6.90 (RSE 28.4); Equation (6) 6.9
    e_crp_cl <- 0.075; label("Log-linear CRP effect on CL/F, applied as (1 + e_crp_cl * log(CRP)) with CRP in mg/L (per log unit)") # Table 4 row 'ln(CRP) on CL/F' Final model = 0.075 (RSE 20.1); Equation (6)
    e_regi_qw_cl <- 0.218; label("Proportional weekly-dosing effect on CL/F, applied as (1 + e_regi_qw_cl * REGI_QW) (fraction)") # Table 4 row 'Dosing regimen on CL/F' Final model = 0.218 (RSE 20.7); Equation (6); Results 3.5 '21.8% higher'

    # IIV reported as %CV of exponential (log-normal) random effects and read
    # on the SD scale: omega = CV/100, variance = (CV/100)^2.
    etalcl ~ 0.1109 # Table 4 row 'IIVCL (%)' Final model = 33.3 (RSE 8.3) -> 0.333^2
    etalvc ~ 0.2894 # Table 4 row 'IIVV (%)' Final model = 53.8 (RSE 10) -> 0.538^2

    propSd <- 0.337; label("Proportional residual error (fraction)") # Table 4 row 'Proportional residual error' Final model = 0.337 (RSE 6.9)
    addSd <- 2.71; label("Additive residual error (mg/L)") # Table 4 row 'Additive residual error (mg/L)' Final model = 2.71 (RSE 7.6)
  })

  model({
    # Equation (6): CL/F = 0.334 * (1 + 6.9*ADA) * (1 + 0.075*ln(CRP)) * (1 + 0.218*REGIMEN)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (1 + e_ada_cl * ADA_POS) *
      (1 + e_crp_cl * log(CRP)) *
      (1 + e_regi_qw_cl * REGI_QW)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg and V/F in L give mg/L
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
