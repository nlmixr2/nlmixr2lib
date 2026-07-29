`Drweesh_2026_adalimumab` <- function() {
  description <- "One-compartment population PK model with first-order subcutaneous absorption and linear elimination for adalimumab originator (Humira) and biosimilars (Amgevita, Hyrimoz) in adults with inflammatory bowel disease and other autoimmune disorders, fit to multicenter therapeutic-drug-monitoring trough data from Saudi Arabia and Qatar (Drweesh 2026). Structural backbone (V/F, IIV variances, residual error) inherited from Marquez-Megias 2023 because Drweesh 2026 reports only ka (fixed) and the typical clearance value."
  reference <- paste(
    "Drweesh H, Alotaibi A, Alqassim HA, Eid MM, Emsaad SA, Saad MO,",
    "Higazy AS, Alturaiki A, Almutairi F, Alhussain A, Assiri A,",
    "Alkofide H, Alruthia Y, Albassam A, Alsultan A.",
    "Comparing the pharmacokinetics of adalimumab originator and biosimilar",
    "product in patients with Inflammatory bowel disease or autoimmune disease.",
    "Saudi Pharm J. 2026;34:14. doi:10.1007/s44446-026-00063-5.",
    "PK structure (V/F, IIV, residual error) adapted from",
    "Marquez-Megias S, Nalda-Molina R, Mas-Serrano P, Ramon-Lopez A.",
    "Population Pharmacokinetic Model of Adalimumab Based on Prior Information",
    "Using Real World Data. Biomedicines. 2023;11(10):2822.",
    "doi:10.3390/biomedicines11102822;",
    "see modellib('Marquez-Megias_2023_adalimumab').",
    sep = " "
  )
  vignette <- "Drweesh_2026_adalimumab"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    species          = "human",
    n_subjects       = 99L,
    n_studies        = 1L,
    n_observations   = 121L,
    age_range        = "adults over 18 years (mean 36.8 +/- 10.9 years)",
    age_mean_sd      = "36.8 +/- 10.9 years",
    weight_mean_sd   = "69.11 +/- 18.7 kg",
    adjusted_bw_mean_sd = "61.68 +/- 12.0 kg",
    lean_bw_mean_sd  = "49.7 +/- 9.9 kg",
    sex_female_pct   = 51.5,
    race_ethnicity   = "Not reported (Saudi Arabian and Qatari multicenter cohort).",
    disease_state    = "Adults with inflammatory bowel disease (Crohn's disease 59.6%, ulcerative colitis 21.2%) or other autoimmune disease (rheumatoid arthritis 13.1%, psoriasis 2.0%).",
    dose_range       = "Subcutaneous adalimumab (mean dose 41.6 +/- 12.7 mg; mean dosing interval 276.6 +/- 80.7 h).",
    regions          = "Multicenter retrospective: Saudi Arabia (Armed Forces Hospital Southern Region, Khamis Mushait; King Abdulaziz Medical City, Riyadh; King Fahd Specialist Hospital, Dammam; King Fahd University Hospital, Khobar) and Qatar (three hospitals in Hamad Medical Corporation).",
    products         = "Humira (originator) 70.7%, Amgevita 16.2%, Hyrimoz 13.1%.",
    aaa_positive_pct = 40.4,
    aaa_negative_pct = 37.4,
    aaa_unknown_pct  = 22.2,
    albumin_mean_sd  = "36.16 +/- 8.5 g/L",
    scr_mean_sd      = "61.95 +/- 20.4 umol/L",
    crp_mean_sd      = "13.8 +/- 26.9 mg/L",
    azathioprine_pct = 43.4,
    trough_mean_sd   = "7.8 +/- 5.9 ug/mL (mg/L)",
    notes            = "Multicenter retrospective TDM analysis. Concentrations measured by ELISA at all centers. Eight patients with all-BLQ samples were removed from the analysis. Covariate effects were assessed by stepwise linear regression on individual Bayesian-estimated clearances (NOT within the popPK model); only age and anti-adalimumab antibody status remained significant in the multivariable regression (R = 0.45, P = 0.001). Slope coefficients were not reported in the source, so the within-model covariate structure is empty here; see vignette Assumptions and deviations."
  )

  ini({
    # Structural PK parameters.
    # ka was fixed by Drweesh 2026 to 0.01 1/h based on prior literature
    # (Marquez-Megias 2021; Kang 2020). Note this differs from the
    # Marquez-Megias 2023 value of 0.00625 1/h (Ternant 2015 reference).
    lka <- fixed(log(0.01))
    label("First-order SC absorption rate constant ka (1/h); fixed")  # Drweesh 2026 Methods Section 2.3

    # Typical CL/F is the mean of the individual empirical-Bayes clearance
    # estimates reported in Drweesh 2026 Section 3.2 (0.018 +/- 0.012 L/hr).
    lcl <- log(0.018)
    label("Apparent clearance CL/F (L/h)")  # Drweesh 2026 Section 3.2

    # Drweesh 2026 does not report V/F. Inherited from the Marquez-Megias
    # 2023 final model (Table 3), which Drweesh 2026 explicitly cites as the
    # structural-form precedent ("a 1 compartment model with linear
    # elimination and first order absorption consistent with prior studies
    # (Marquez-Megias et al. 2023; Marquez-Megias et al. 2021)").
    lvc <- log(7.76)
    label("Apparent volume of distribution V/F (L)")  # Marquez-Megias 2023 Table 3 (Final model); inherited because Drweesh 2026 does not report V/F

    # Inter-individual variability. Drweesh 2026 does not report omega
    # values; values inherited from Marquez-Megias 2023 Table 3 (Final
    # model: omega_CL = 0.667, omega_V = 0.477, interpreted as Monolix
    # log-scale standard deviations and squared to variances).
    etalcl ~ 0.667^2  # Marquez-Megias 2023 Table 3 (Final model); inherited because Drweesh 2026 does not report IIV magnitudes
    etalvc ~ 0.477^2  # Marquez-Megias 2023 Table 3 (Final model); inherited because Drweesh 2026 does not report IIV magnitudes

    # Residual error. Drweesh 2026 Section 3.2 states the proportional
    # error model was best, but does not report the magnitude. Magnitude
    # inherited from Marquez-Megias 2023 Table 3 (Final model: propSd =
    # 0.547, the proportional-only residual after dropping the high-RSE
    # additive term).
    propSd <- 0.547
    label("Proportional residual error (fraction)")  # Form: Drweesh 2026 Section 3.2; magnitude: Marquez-Megias 2023 Table 3 inherited
  })

  model({
    # Individual PK parameters for a typical patient (Drweesh 2026 did not
    # report a within-model covariate structure; covariate associations
    # were assessed by post-hoc linear regression on EB clearance estimates
    # without published slopes -- see population$notes and the vignette
    # Assumptions and deviations section).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Concentration: dose in mg, V/F in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
