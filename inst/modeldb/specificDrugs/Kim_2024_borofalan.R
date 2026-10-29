Kim_2024_borofalan <- function() {
  description <- paste(
    "Two-compartment population PK model for whole-blood boron-10 after",
    "intravenous infusion of borofalan (10B) (L-4-boronophenylalanine, BPA,",
    "given as the BPA-fructose complex) in cancer patients undergoing boron",
    "neutron capture therapy. Fitted in NONMEM 7.4.3 (ADVAN13, FOCE-I) to",
    "377 blood boron concentrations from 27 patients, digitised by the",
    "authors from seven published clinical reports. Parameterised in",
    "micro-constant form exactly as the authors report it (V1, k10, k12,",
    "k21) with log-normal IIV on all four and a combined additive +",
    "proportional residual error. Dosing and volume are per kilogram of",
    "body weight: the dose is mg/kg of BPA, V1 is L/kg, and the model",
    "converts the BPA amount to boron-10 with the authors' factor",
    "10/208.21, so the output is the boron-10 concentration in ug/g (ppm).",
    "The paper uses this model as the prior for Bayesian (MAP) prediction",
    "of the blood boron level during neutron irradiation.",
    sep = " "
  )
  reference <- paste(
    "Kim W, Won JY, Yi J, Choi SC, Lee SM, Mun K, Lim HS.",
    "PK Modeling of L-4-Boronophenylalanine and Development of Bayesian",
    "Predictive Platform for L-4-Boronophenylalanine PKs for Boron Neutron",
    "Capture Therapy.",
    "Pharmaceuticals (Basel). 2024;17(3):301.",
    "doi:10.3390/ph17030301.",
    sep = " "
  )
  vignette <- "Kim_2024_borofalan"
  units <- list(
    time = "min",
    dosing = "mg/kg",
    concentration = "ug/g"
  )

  # Per-kg normalisation: the dose is mg/kg of BPA and V1 is L/kg, so the
  # states hold mg/kg of BPA. The boron-10 concentration is
  # central * (10/208.21) / vc in mg/L, read as ug/g (ppm) as the source does
  # (blood density taken as 1 g/mL).
  compartmentData <- list(
    central = list(analyte = "boron-10 (borofalan)", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "boron-10 (borofalan)", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 27L,
    n_studies = 7L,
    n_observations = 377L,
    disease_state = paste(
      "Cancer patients treated with boron neutron capture therapy",
      "(glioma / glioblastoma, malignant melanoma, and head and neck",
      "squamous cell carcinoma in the pooled source reports)."
    ),
    dose_range = paste(
      "Single intravenous infusion of 10B-enriched BPA (BPA-fructose",
      "complex) dosed per kg; the pooled reports used different doses and",
      "infusion durations (Figure 2 shows 250 mg/kg over 65 min and",
      "290 mg/kg over 120 min). The paper's simulations use 500 mg/kg over",
      "180 min."
    ),
    regions = "Pooled published reports from Finland, the USA, China, Europe (EORTC) and Japan.",
    notes = paste(
      "Individual whole-blood boron-10 concentrations were graphically",
      "extracted from seven published reports (Kim 2024 references 8 and",
      "10-15) and pooled; no individual demographics or covariates were",
      "available and none were tested. See Kim 2024 Methods 2.1."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Kim 2024 Table 2 'Estimates' column.
    # Paper symbols: K10 is the canonical kel; V1 is vc. V1 is printed as
    # 'L' but is per kilogram: Figure 4A prints an individual Bayesian
    # V1 of '0.22902 l/kg', and only L/kg reproduces the mg/kg-dosed
    # boron levels of Figures 2 and 4A.
    # ------------------------------------------------------------------
    lkel <- log(0.006); label("Elimination rate constant from central (1/min)") # Table 2, K10 = 0.006 1/min (95% CI 0.005-0.007)
    lk12 <- log(0.023); label("Rate constant central -> peripheral1 (1/min)") # Table 2, K12 = 0.023 1/min (95% CI 0.018-0.028)
    lk21 <- log(0.012); label("Rate constant peripheral1 -> central (1/min)") # Table 2, K21 = 0.012 1/min (95% CI 0.009-0.014)
    lvc <- log(0.252); label("Central volume of distribution (L/kg)") # Table 2, V1 = 0.252 (95% CI 0.216-0.288)

    # ------------------------------------------------------------------
    # IIV -- Table 2 prints omega^2 (variance of eta) with the CV% in
    # parentheses; every CV% equals sqrt(exp(omega^2) - 1), e.g.
    # sqrt(exp(0.071) - 1) = 27.1%. Methods 2.1: P_i = P_TV * exp(eta).
    # No off-diagonal elements are reported (diagonal OMEGA).
    # ------------------------------------------------------------------
    etalkel ~ 0.050 # Table 2, IIV K10 = 0.050 (CV 22.6%)
    etalk12 ~ 0.071 # Table 2, IIV K12 = 0.071 (CV 27.2%)
    etalk21 ~ 0.037 # Table 2, IIV K21 = 0.037 (CV 19.4%)
    etalvc ~ 0.061 # Table 2, IIV V1 = 0.061 (CV 25.0%)

    # ------------------------------------------------------------------
    # Residual error -- Results 3.1: Cobs = Cpred + eps1 + eps2 * Cpred.
    # Table 2 prints eps1 = 0.540 (additive) and eps2 = 0.001
    # (proportional). These are read as the NONMEM SIGMA variances, the
    # same raw-estimate scale as the omega^2 rows of the same table, so
    # the SDs are their square roots. See the vignette 'Assumptions and
    # deviations' section for the SD reading and its (small) effect.
    # ------------------------------------------------------------------
    addSd <- sqrt(0.540); label("Additive residual error SD (ug/g)") # Table 2, eps1 additive = 0.540 (variance)
    propSd <- sqrt(0.001); label("Proportional residual error SD (fraction)") # Table 2, eps2 proportional = 0.001 (variance)
  })

  model({
    # Individual parameters (log-normal IIV, Methods 2.1)
    kel <- exp(lkel + etalkel)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    vc <- exp(lvc + etalvc)

    # Two-compartment linear disposition (Results 3.1), amounts in mg/kg BPA
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Boron-10 concentration: Methods 2.1 applies the factor 10/208.21
    # (boron-10 mass / BPA molecular weight) to the BPA dose. mg/kg over
    # L/kg gives mg/L, reported as ug/g (ppm).
    Cc <- central * (10 / 208.21) / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
