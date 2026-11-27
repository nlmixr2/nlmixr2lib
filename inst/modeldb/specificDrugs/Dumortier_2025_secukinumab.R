Dumortier_2025_secukinumab <- function() {
  description <- paste(
    "Two-compartment population PK model for secukinumab (anti-IL-17A IgG1)",
    "after subcutaneous and intravenous dosing in adults with psoriatic",
    "arthritis, ankylosing spondylitis or non-radiographic axial",
    "spondyloarthritis (Dumortier 2025; 14 phase II/III trials, 5270",
    "patients). First-order SC absorption with ka fixed at 0.18/day and",
    "bioavailability 0.783; zero-order IV infusion into the central",
    "compartment; linear elimination. Baseline body weight (reference 84 kg)",
    "scales CL, Vc and Vp by power functions, and CL is 1.13-fold higher in",
    "Asian patients. Correlated IIV on CL and Vc, independent IIV on Vp,",
    "proportional residual error. Used to select the approved 1.75 mg/kg IV",
    "q4w regimen (with or without a 6 mg/kg IV loading dose)."
  )
  reference <- paste(
    "Dumortier T, Valenzuela G, Churchill M, Mijatovic J, Bruin G,",
    "Pricop L, Richards H, Renard D, Singhal A, Marathe A.",
    "Model-Informed Drug Development-Based Bridging from Subcutaneous to",
    "Intravenous Secukinumab Dosing: Approval in Psoriatic Arthritis and",
    "Axial Spondyloarthritis. Clin Pharmacol Ther. 2025;118(2):480-488.",
    "doi:10.1002/cpt.3716."
  )
  vignette <- "Dumortier_2025_secukinumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")
  dosing <- c("depot", "central")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline value (NONMEM WEIGHT0), not time-varying. Power effects",
        "on CL, Vc and Vp normalised to 84 kg (Table S3 footnote 'The",
        "reference is an 84 kg non-Asian patient'). Q carries no weight",
        "effect in the control stream, although the main text lists one."
      ),
      source_name = "WEIGHT0"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = non-Asian",
      notes = paste(
        "Control stream: ASIAN = 1 when RACEN = 30, else 0; CL is",
        "multiplied by THETA(10)**ASIAN. Asian patients were 15% of the",
        "model-development population (765/5270, Table 1)."
      ),
      source_name = "ASIAN (derived from RACEN = 30)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "secukinumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "secukinumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "secukinumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 5270L,
    n_studies = 14L,
    n_observations = "24149 serum secukinumab concentrations (4396 after IV and 19753 after SC dosing)",
    age_range = "mean 45.1 years (SD 12.9)",
    weight_range = "mean 81.6 kg (SD 19.4)",
    sex_female_pct = 45,
    race_ethnicity = c(
      White = 80,
      Asian = 15,
      Other = 3,
      Black = 1,
      NativeAmerican = 1,
      PacificIslander = 0,
      Unknown = 0
    ),
    disease_state = "Psoriatic arthritis (58%), ankylosing spondylitis (31%) and non-radiographic axial spondyloarthritis (10%)",
    dose_range = paste(
      "0.1, 1 or 10 mg/kg IV at weeks 0 and 3 (phase II); 75, 150 or 300 mg SC",
      "q4w with or without weekly SC loading or 10 mg/kg IV loading at weeks",
      "0, 2 and 4; 6 mg/kg IV loading then 3 mg/kg IV q4w (INVIGORATE-2)"
    ),
    regions = "Multinational",
    notes = paste(
      "Model-development population of Dumortier 2025 Table 1 (n = 5270).",
      "INVIGORATE-1 (n = 516, axSpA) was held out for external validation",
      "(Figure S1). Trial list in Table S1."
    )
  )

  ini({
    # Table S3 (final NONMEM estimates; reference = 84 kg non-Asian patient).
    # The $THETA block of the supplementary control stream holds the
    # INITIAL values (e.g. CL 0.181, V2 3.74) and is not used here.
    lcl <- log(0.171); label("Clearance for an 84 kg non-Asian patient (L/day)") # Table S3 'CL, L/day' 0.171 (RSE 0.836%)
    lvc <- log(3.32); label("Central volume for an 84 kg patient (L)") # Table S3 'Vc, L' 3.32 (RSE 1.44%)
    lq <- log(0.35); label("Intercompartmental clearance (L/day)") # Table S3 'Q, L/day' 0.35 (RSE 4.37%)
    lvp <- log(2.63); label("Peripheral volume for an 84 kg patient (L)") # Table S3 'Vp, L' 2.63 (RSE 1.81%)
    lka <- fixed(log(0.18)); label("SC absorption rate constant (1/day)") # Table S3 'KA, 1/day' 0.18 (fixed); control stream '(0.18 FIX)'; value from Bruin 2017 psoriasis model
    lfdepot <- log(0.783); label("SC bioavailability (fraction)") # Table S3 'F1, bioavailability' 0.783 (RSE 0.789%)

    e_wt_cl <- 0.892; label("Power exponent of body weight on CL (unitless)") # Table S3 'CL ~weight' 0.892 (RSE 2.32%)
    e_wt_vc <- 0.651; label("Power exponent of body weight on Vc (unitless)") # Table S3 'Vc ~weight' 0.651 (RSE 7.51%)
    e_wt_vp <- 0.712; label("Power exponent of body weight on Vp (unitless)") # Table S3 'Vp ~weight' 0.712 (RSE 8.02%)
    e_race_asian_cl <- 1.13; label("Multiplicative factor on CL for Asian patients (unitless)") # Table S3 'CL ~Asian' 1.13 (RSE 1.18%)

    # Table S3 inter-individual variances (CV% = sqrt(variance) x 100).
    etalcl + etalvc ~ c(0.0912, 0.0489, 0.0691) # Table S3 omega2 CL 0.0912, cov 0.0489 (R = 0.62), omega2 Vc 0.0691
    etalvp ~ 0.113 # Table S3 omega2 Vp 0.113 (33.6% CV)

    propSd <- sqrt(0.0544); label("Proportional residual error (fraction)") # Table S3 'delta2 prop' variance 0.0544 (23.3% CV)
  })

  model({
    asian_cl <- e_race_asian_cl^RACE_ASIAN

    cl <- exp(lcl + etalcl) * (WT / 84)^e_wt_cl * asian_cl
    vc <- exp(lvc + etalvc) * (WT / 84)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp + etalvp) * (WT / 84)^e_wt_vp
    ka <- exp(lka)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot)

    # Dose in mg and volume in L give mg/L = ug/mL, the unit of every
    # concentration reported in the paper.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
