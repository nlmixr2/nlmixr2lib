Li_2026_midazolam <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption for",
    "oral microdose midazolam (CYP3A probe) in Chinese healthy volunteers and",
    "type 2 diabetes patients (Li 2026). CL/F rises with the stool relative",
    "abundance of the gut genus Salmonella (power model centred on 0.187%);",
    "type 2 diabetes itself is not a covariate."
  )
  reference <- paste(
    "Li Y, Yang J, Zhou E, Geng K, Lei Z, Yao X, Feng F, Jin Y, Kong W, Hu Y,",
    "Liang H, Li H, Yang G-P, Hong T, Liu D.",
    "A microdose cocktail and population pharmacokinetic study suggests",
    "potentially reduced CYP3A and P-gp activities in type 2 diabetes.",
    "Front Pharmacol. 2026;17:1924307. doi:10.3389/fphar.2026.1924307.",
    sep = " "
  )
  vignette <- "Li_2026_microdose_cocktail_t2d"
  units <- list(time = "h", dosing = "ug", concentration = "pg/mL")

  covariateData <- list(
    GUT_SALMONELLA_PCT = list(
      description = "Stool relative abundance of the gut genus Salmonella",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "16S rRNA V3-V4 sequencing (primers 338F/806R, MGI-2000) of a",
        "pre-dose stool sample. Enters as (Salmonella / 0.187)^0.124 on CL/F",
        "(Li 2026 Eq. 1), so the effect is 1 at 0.187%. Reported means 0.619",
        "+/- 0.505% in healthy volunteers and 0.0880 +/- 0.124% in type 2",
        "diabetes (Results, 'PopPK analysis'); Figure 3a/b spans about",
        "0-1.6%. The power form gives CL/F = 0 at an abundance of exactly 0,",
        "so supply a positive value (e.g. the sequencing detection floor)."
      ),
      source_name = "Salmonella"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "midazolam", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "midazolam", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "midazolam", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 27L,
    n_studies = 2L,
    age_range = "18-65 years (inclusion criterion)",
    age_mean = "35.6 years (healthy volunteers); 49.2 years (type 2 diabetes)",
    weight_mean = "58.8 kg (healthy volunteers); 67.1 kg (type 2 diabetes)",
    sex_female_pct = 66.7,
    race_ethnicity = c(Chinese = 100),
    disease_state = paste(
      "14 healthy volunteers and 13 type 2 diabetes patients (HbA1c 7.62 +/-",
      "0.77%, disease duration 1.5-8.5 years, normal hepatic and renal",
      "function; 7 on metformin, 2 acarbose, 2 glyburide, 1 dapagliflozin)."
    ),
    dose_range = paste(
      "Single oral dose of midazolam 10 ug in a fasted microdose cocktail",
      "with dabigatran etexilate 375 ug, pitavastatin 10 ug, rosuvastatin",
      "50 ug and atorvastatin 100 ug."
    ),
    regions = "China (Peking University Third Hospital, Beijing)",
    renal_function = "Mean eGFR 106.6 (healthy volunteers) and 103.6 (type 2 diabetes) mL/min/1.73 m^2",
    notes = paste(
      "Li 2026 Table 1. Healthy-volunteer data come from the earlier Kong",
      "2025 cocktail study and were sampled at 0.5-24 h; type 2 diabetes",
      "patients were sampled at 1, 4 and 24 h. Microbiome data were",
      "available for 9 healthy volunteers and 8 patients (Figure 3).",
      "NCT06273215."
    )
  )

  ini({
    lka <- log(0.773); label("Absorption rate constant (1/h)") # Table S4 MDZ Ka 0.773 (RSE 12.0%)
    lcl <- log(57.2); label("Apparent clearance CL/F at Salmonella 0.187% (L/h)") # Table S4 MDZ CL/F 57.2 (RSE 9.40%); Eq. 1
    lvc <- log(8.32); label("Apparent central volume Vc/F (L)") # Table S4 MDZ VC/F 8.32 (RSE 116%)
    lvp <- log(326); label("Apparent peripheral volume Vp/F (L)") # Table S4 MDZ VT/F 326 (RSE 19.9%)
    lq <- log(44.1); label("Apparent intercompartmental clearance Q/F (L/h)") # Table S4 MDZ Q/F 44.1 (RSE 16.2%)

    e_gut_salmonella_pct_cl <- 0.124; label("Power exponent of Salmonella/0.187 on CL/F (unitless)") # Table S4 'Salmonella on CL/F' 0.124 (RSE 32.7%); Eq. 1

    etalcl ~ 0.102 # Table S4 MDZ omega^2 (CL/F) 0.102 (RSE 44.3%)

    propSd <- sqrt(0.193); label("Proportional residual error (fraction)") # Table S4 MDZ sigma^2 0.193 (RSE 17.4%); SD = sqrt(variance)
  })
  model({
    # Eq. 1: CL/F = 57.2 x (Salmonella / 0.187)^0.124
    salm_cl <- (GUT_SALMONELLA_PCT / 0.187)^e_gut_salmonella_pct_cl

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * salm_cl
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ug / L = ng/mL; x 1000 -> pg/mL (the paper's concentration unit)
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
