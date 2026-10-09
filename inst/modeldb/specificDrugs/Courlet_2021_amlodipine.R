Courlet_2021_amlodipine <- function() {
  description <- "One-compartment population PK model for oral amlodipine in adults living with HIV (Courlet 2021), with first-order absorption after a lag time, between-subject variability on apparent clearance only, an additive residual error, and two antiretroviral drug-drug-interaction effects on CL/F: a -49% linear effect of strong CYP3A4-inhibiting ARVs (boosted darunavir, boosted atazanavir, cobicistat-boosted elvitegravir) and a +140% linear effect of efavirenz."
  reference <- paste(
    "Courlet P, Guidi M, Alves Saldanha S, Cavassini M, Stoeckle M, Buclin T,",
    "Marzolini C, Decosterd LA, Csajka C. Population pharmacokinetic modelling",
    "to quantify the magnitude of drug-drug interactions between amlodipine",
    "and antiretroviral drugs. Eur J Clin Pharmacol. 2021;77(7):979-987.",
    "doi:10.1007/s00228-020-03060-2.",
    sep = " "
  )
  vignette <- "Courlet_2021_amlodipine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "amlodipine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "amlodipine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CONMED_CYP3A4_INH_STRONG = list(
      description = "Concomitant antiretroviral strong CYP3A4 inhibitor (1 = ritonavir-boosted darunavir, cobicistat-boosted darunavir, ritonavir-boosted atazanavir or cobicistat-boosted elvitegravir; 0 = none)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no strong CYP3A4-inhibiting ARV)",
      notes = "Linear multiplicative effect on CL/F: (1 - 0.49 * CONMED_CYP3A4_INH_STRONG), Courlet 2021 Table 2 footnote 'Final model'. The paper classifies the four boosted ARVs as strong CYP3A4 inhibitors per the FDA list (Methods, Covariate model) and pools them into one indicator; ritonavir-boosted darunavir dominates (27 of the 36 inhibitor-exposed samples, Table 1).",
      source_name = "CYP3A4 inhibitors"
    ),
    CONMED_EFV = list(
      description = "Concomitant efavirenz (1 = on efavirenz, 0 = not on efavirenz)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no efavirenz; ARVs without CYP3A4 interaction potential)",
      notes = "Linear multiplicative effect on CL/F: (1 + 1.40 * CONMED_EFV), Courlet 2021 Table 2 footnote 'Final model', i.e. a 2.4-fold CL/F. No subject received a CYP3A4 inhibitor and efavirenz together (Results, Model-based simulations). Etravirine (also a moderate CYP3A4 inducer) was not retained because it was co-prescribed with ritonavir-boosted darunavir in both etravirine patients; do not code etravirine as CONMED_EFV.",
      source_name = "efavirenz"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested and not retained (Delta OFV > -1.1, p > 0.29; Courlet 2021 Results)."
    ),
    SEXF = list(
      description = "Biological sex, 1 = female",
      units = "(binary)",
      type = "binary",
      notes = "Tested and not retained (Courlet 2021 Results)."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested and not retained (Courlet 2021 Results)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Tested and not retained (Courlet 2021 Results)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = "Tested and not retained (Courlet 2021 Results)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = "Tested and not retained (Courlet 2021 Results)."
    ),
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault)",
      units = "mL/min",
      type = "continuous",
      notes = "Tested and not retained (Courlet 2021 Results)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 55L,
    n_studies = 2L,
    n_observations = 163L,
    age_median = "61 years (IQR 53-70)",
    weight_median = "79 kg (IQR 71-91)",
    sex_female_pct = 25,
    disease_state = "People living with HIV on antiretroviral therapy and taking amlodipine for hypertension",
    dose_range = "Amlodipine 2.5-10 mg orally once daily (three subjects 5 mg twice daily); steady state assumed for all subjects",
    regions = "Switzerland (Lausanne and Basel)",
    co_medication = "Strong CYP3A4-inhibiting ARVs in 36 of 163 samples (ritonavir-boosted darunavir 27); efavirenz in 7; etravirine 18; nevirapine 8; dolutegravir 103 (Table 1)",
    notes = "Swiss HIV Cohort Study project 815 (sparse, median 1 sample per subject) plus the rich-sampling PK study NCT03515772 (8 subjects, 84 concentrations, sampling to 28 h post-dose). Baseline demographics from Courlet 2021 Table 1."
  )

  ini({
    lka <- log(0.69); label("First-order absorption rate constant ka (1/h)") # Table 2, 'ka (h-1)' = 0.69 (RSE 25%)
    ltlag <- log(0.87); label("Absorption lag time ALAG (h)") # Table 2, 'ALAG (h)' = 0.87 (RSE 28%)
    lvc <- log(1000); label("Apparent volume of distribution V/F (L)") # Table 2, 'V/F (L)' = 1000 (RSE 16%)
    lcl <- log(17.0); label("Apparent clearance CL/F without interacting ARVs (L/h)") # Table 2, 'CL/F (L/h)' = 17.0 (RSE 9%)

    e_conmed_cyp3a4_inh_strong_cl <- -0.49; label("Fractional change in CL/F with a strong CYP3A4-inhibiting ARV (unitless)") # Table 2, 'theta CYP3A4 inhibitors' = -0.49 (RSE 12%)
    e_conmed_efv_cl <- 1.40; label("Fractional change in CL/F with efavirenz (unitless)") # Table 2, 'theta efavirenz' = 1.40 (RSE 37%)

    etalcl ~ 0.16246 # Table 2, 'BSV CL (CV%)' = 42; omega^2 = log(1 + 0.42^2)

    addSd <- 2.85; label("Additive residual error (ng/mL)") # Table 2, 'sigma add (ng/mL)' = 2.85 (RSE 11%)
  })

  model({
    ka <- exp(lka)
    tlag <- exp(ltlag)
    vc <- exp(lvc)
    # Table 2 footnote: CL/F = 17.0 x (1 - 0.49 x CYP3A4 inhibitors) x (1 + 1.40 x efavirenz)
    cl <- exp(lcl + etalcl) *
      (1 + e_conmed_cyp3a4_inh_strong_cl * CONMED_CYP3A4_INH_STRONG) *
      (1 + e_conmed_efv_cl * CONMED_EFV)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    # central in mg, vc in L -> mg/L; x 1000 -> ng/mL
    Cc <- 1000 * central / vc
    Cc ~ add(addSd)
  })
}
