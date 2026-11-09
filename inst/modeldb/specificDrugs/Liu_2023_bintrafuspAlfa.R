Liu_2023_bintrafuspAlfa <- function() {
  description <- paste(
    "Two-compartment population PK model for bintrafusp alfa (a bifunctional anti-PD-L1 / TGF-beta trap",
    "fusion protein) given as a 1-hour IV infusion of 500 or 1200 mg every 2 weeks to adults with",
    "second-line non-small cell lung cancer (phase I trial NCT02517398, 2L NSCLC expansion cohort).",
    "Time-constant linear clearance; peripheral volume and intercompartmental clearance fixed to the",
    "Wilkins 2019 pan-tumour model; correlated IIV on CL and Vc; combined additive + proportional",
    "residual error whose magnitude carries its own between-subject random effect."
  )
  reference <- paste(
    "Liu H, Milenkovic-Grisic AM, Krishnan SM, Jonsson S, Friberg LE, Girard P, Venkatakrishnan K,",
    "Vugmeyster Y, Khandelwal A, Karlsson MO.",
    "A multistate modeling and simulation framework to learn dose-response of oncology drugs:",
    "Application to bintrafusp alfa in non-small cell lung cancer.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(11):1738-1750. doi:10.1002/psp4.12976.",
    "PMCID: PMC10681430. Parameter values from Supplementary Table S2; structure from Appendix S1",
    "'Population pharmacokinetic model'. Vp and Q fixed to Wilkins JJ et al. Adv Ther.",
    "2019;36(9):2414-2433."
  )
  vignette <- "Liu_2023_bintrafuspAlfa"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # etaRUV is the random effect on the magnitude of the residual error
  # (Supplementary Table S2 row 'sigma^2 EPS-pk'); see the model() comment.
  paper_specific_etas <- c("etaRUV")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "bintrafusp alfa", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "bintrafusp alfa", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 80L,
    n_studies = 1L,
    age_range = "38-85 years",
    age_mean = "63.9 years (500 mg arm), 62.1 years (1200 mg arm)",
    weight_range = "41.9-109.0 kg",
    weight_mean = "68.0 kg (500 mg arm), 70.8 kg (1200 mg arm)",
    sex_female_pct = 28.75,
    race_ethnicity = c(White = 48.75, Black = 1.25, Asian = 50),
    disease_state = "Advanced non-small cell lung cancer, second line (2L NSCLC); ECOG 0-1; 81% with metastases at baseline",
    dose_range = "500 or 1200 mg IV every 2 weeks (1:1 randomized, 40 patients per arm)",
    regions = "USA, United Kingdom, Belgium, Italy, Spain, South Korea, Taiwan (25 sites)",
    notes = paste(
      "Baseline characteristics from Liu 2023 Table 1 (n = 40 per arm) and Supplementary Table S4.",
      "998 plasma concentrations; PK sampling schedule in Supplementary Table S1 (day 1 pre-dose, end of",
      "infusion and 4 h after start of infusion; days 15, 43 and 85 pre-dose, end of infusion and 2-8 h",
      "after end of infusion; pre-dose troughs thereafter)."
    )
  )

  ini({
    lcl <- log(0.0143)
    label("Clearance CL (L/h)") # Supplementary Table S2: CL = 0.0143 L/h (RSE 3.5%)
    lvc <- log(3.02)
    label("Central volume of distribution Vc (L)") # Supplementary Table S2: Vc = 3.02 L (RSE 2.6%)
    lvp <- fixed(log(0.483))
    label("Peripheral volume of distribution Vp (L)") # Supplementary Table S2: Vp = 0.483 L, Fixed (Wilkins 2019)
    lq <- fixed(log(0.00512))
    label("Intercompartmental clearance Q (L/h)") # Supplementary Table S2: Q = 0.00512 L/h, Fixed (Wilkins 2019)

    etalcl + etalvc ~ c(
      0.0795,
      0.0472,
      0.0424
    ) # Supplementary Table S2: omega^2 CL = 0.0795, cov(CL, Vc) = 0.0472, omega^2 Vc = 0.0424

    etaRUV ~ 0.260 # Supplementary Table S2: 'sigma^2 EPS-pk (Variance)' = 0.260 (RSE 27%), random effect on the residual-error magnitude

    propSd <- 0.17
    label("Proportional residual error (fraction)") # Supplementary Table S2: PK-proportional residual error (CV) = 0.17
    addSd <- 7.57
    label("Additive residual error (ug/mL)") # Supplementary Table S2: PK-additive residual error (SD) = 7.57 ug/mL
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L = ug/mL
    Cc <- central / vc

    # Residual SD = sqrt(propSd^2 * Cc^2 + addSd^2) * exp(etaRUV): the
    # combined error of Appendix S1 with a subject-level random effect on its
    # magnitude. Scaling both components by exp(etaRUV) is exactly that form.
    propSdi <- propSd * exp(etaRUV)
    addSdi <- addSd * exp(etaRUV)
    Cc ~ add(addSdi) + prop(propSdi) + combined2()
  })
}
