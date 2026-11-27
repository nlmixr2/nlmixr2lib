Choi_2018_kaurenoic_acid <- function() {
  description <- paste(
    "Two-compartment PK model with first-order absorption for kaurenoic",
    "acid (KAU), the representative diterpene of Araliae continentalis",
    "Radix (Aralia continentalis Kitagawa root), after a single oral dose",
    "of 1.0054 g of the herbal extract powder (1.147 mg KAU) to healthy",
    "Korean men. Fitted with the naive-pooled engine of Phoenix NLME 8.0,",
    "so the model carries no between-subject variability; residual error",
    "is proportional. CL, Q, Vc and Vp are apparent (/F) values. Dose is",
    "entered in micrograms of KAU so that central/vc is in ug/L (ng/mL).",
    sep = " "
  )
  reference <- paste(
    "Choi EJ, Choi GW, Yang SJ, Lee YB, Cho HY.",
    "Pharmacokinetic Profile of Kaurenoic Acid after Oral Administration",
    "of Araliae Continentalis Radix Extract Powder to Humans.",
    "Pharmaceutics. 2018;10(4):253.",
    "doi:10.3390/pharmaceutics10040253.",
    sep = " "
  )
  vignette <- "Choi_2018_kaurenoic_acid"
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "ug/L"
  )

  compartmentData <- list(
    depot = list(analyte = "kaurenoic acid", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "kaurenoic acid", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "kaurenoic acid", units = "ug", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  # Covariates screened in the exploratory mixed-effects analysis
  # (Choi 2018 Results 3.4, final paragraph) and not retained: only
  # CrCl on CL/F showed a relationship and it did not significantly
  # improve the basic model, so the basic-model estimates were published.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on the primary parameters; not retained (Choi 2018 Results 3.4)."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened on the primary parameters; not retained (Choi 2018 Results 3.4)."
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      notes = "Screened on the primary parameters; not retained (Choi 2018 Results 3.4)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened (paper's SCr); not retained (Choi 2018 Results 3.4)."
    ),
    CRCL = list(
      description = "Creatinine clearance, absolute (not BSA-normalised)",
      units = "mL/min",
      type = "continuous",
      notes = paste(
        "Paper's CrCl (Table 3, mL/min). The only covariate relationship",
        "identified (on CL/F), but it did not significantly improve the basic",
        "model and no coefficient was published (Choi 2018 Results 3.4)."
      )
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained (Choi 2018 Results 3.4)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained (Choi 2018 Results 3.4)."
    ),
    GGT = list(
      description = "Gamma-glutamyltransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened (paper's GGTP); not retained (Choi 2018 Results 3.4)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 10L,
    n_studies = 1L,
    age_range = "21-33 years (mean 25.8, median 24.5)",
    weight_range = "60.4-86.7 kg (mean 73.3, median 73.65)",
    height_range = "157.0-183.0 cm in Table 3 (Methods 2.5 prints 167-183 cm)",
    sex_female_pct = 0,
    race_ethnicity = "Korean (100%)",
    disease_state = "Healthy volunteers",
    renal_function = "CrCl 98.0-139.1 mL/min (mean 115.0); SCr 0.80-1.10 mg/dL",
    dose_range = paste(
      "Single oral dose of 1.0054 g Araliae continentalis Radix extract",
      "powder (the usual clinical dose approved by the Korean MFDS)",
      "containing 1.147 mg kaurenoic acid, taken with 240 mL water after",
      "a >= 10 h fast."
    ),
    regions = "Korea",
    notes = paste(
      "Open-label investigator-initiated study (KCT0001211). Plasma",
      "sampled pre-dose and at 0.25, 0.5, 0.75, 1, 2, 3, 4, 6, 8, 10 and",
      "12 h; KAU quantified by UPLC-MS/MS (LLOQ 0.2 ng/mL). Demographics",
      "in Choi 2018 Table 3."
    )
  )

  ini({
    # Structural parameters -- Choi 2018 Table 4 'Primary parameters'
    # (naive-pooled Phoenix NLME estimates; CV% is the parameter RSE).
    # All are apparent values that include the unknown oral fraction F
    # (Methods 2.6).
    lka <- log(1.72)  ; label("First-order absorption rate constant ka (1/h)")          # Table 4: Ka = 1.72 h-1 (CV 28.79%)
    lcl <- log(23.89) ; label("Apparent clearance CL/F (L/h)")                            # Table 4: CL/F = 23.89 L/h (CV 9.30%)
    lvc <- log(24.44) ; label("Apparent central volume V/F (L)")                         # Table 4: V/F = 24.44 L (CV 29.16%)
    lq  <- log(15.55) ; label("Apparent inter-compartmental clearance CLD/F (L/h)")     # Table 4: CLD/F = 15.55 L/h (CV 13.77%)
    lvp <- log(64.05) ; label("Apparent peripheral volume V2/F (L)")                    # Table 4: V2/F = 64.05 L (CV 21.05%)

    # No between-subject variability: the published estimates come from the
    # naive-pooled engine (Methods 2.6), and the exploratory mixed-effects
    # analysis (Results 3.4) reports no omega values.

    # Residual error -- Methods eq. 3, Cobs = Cpred * (1 + eps); Table 4
    # 'Proportional residual error' = 0.48. Phoenix NLME reports the
    # residual-error parameter as a standard deviation.
    propSd <- 0.48 ; label("Proportional residual error SD (fraction)") # Table 4: epsilon = 0.48 (CV 11.09%)
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Choi 2018 Methods eq. 2 (absorption compartment Aa, central A1,
    # peripheral A2).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in ug and volume in L, so Cc is in ug/L (= ng/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
