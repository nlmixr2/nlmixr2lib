Yarur_2025_ustekinumab <- function() {
  description <- paste0(
    "Two-compartment population PK model with first-order subcutaneous ",
    "absorption for ustekinumab in adults and adolescents with ",
    "inflammatory bowel disease, with concentrations measured by the ",
    "drug-tolerant homogeneous mobility shift assay (HMSA) (Yarur 2025). ",
    "This is the 'Final Bayesian Priors' model of Table 1: the Wang 2022 ",
    "Crohn's disease model structure re-fitted to HMSA re-assayed samples ",
    "(albumin power effect on CL), plus an antibodies-to-ustekinumab ",
    "(ATU) effect on CL estimated from a clinical PK-laboratory database, ",
    "with between-subject variability on CL only. The paper used it, with ",
    "every parameter fixed, to obtain empirical-Bayes clearances for the ",
    "Yarur_2025_ustekinumab_* exposure-response models."
  )
  reference <- paste(
    "Yarur AJ, Dervieux T, Ungaro R, Spencer EA, Bruss A, Nunez L,",
    "Berens B, Vermeire S, Wang Z, Panetta JC, Dreesen E, Dubinsky MC.",
    "Ustekinumab Drug Clearance Is Better Associated with Disease Control",
    "than Serum Trough Concentrations in a Prospective Cohort of",
    "Inflammatory Bowel Disease. Pharmaceutics. 2025;17(2):187.",
    "doi:10.3390/pharmaceutics17020187.",
    "Model structure from Wang Z, Verstockt B, Sabino J, Vermeire S,",
    "Ferrante M, Declerck P, Dreesen E. Br J Clin Pharmacol.",
    "2022;88(1):323-335. doi:10.1111/bcp.14971.",
    sep = " "
  )
  vignette <- "Yarur_2025_ustekinumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")
  dosing <- c("depot", "central")

  covariateData <- list(
    ALB = list(
      description = "Serum albumin concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Power effect on CL centred at 4.0 g/dL (Table 1 footnote, ",
        "covariate model 2: CL = CLpop x (Albumin/4.0)^beta). The paper ",
        "reports albumin in g/dL; the canonical column is g/L, so ",
        "model() converts with ALB / 10 before applying the published ",
        "4.0 g/dL reference. Measured at each trough sample by ",
        "immunoassay (IMMAGE 800), so it may be time-varying. Study ",
        "median 3.8 g/dL, IQR 3.5-4.1 (Table 4)."
      ),
      source_name = "ALB"
    ),
    ADA_POS = list(
      description = "Antibodies to ustekinumab (ATU) positive status",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ATU-negative)",
      notes = paste0(
        "ATU measured by the drug-tolerant HMSA; > 1.6 U/mL is positive ",
        "(Methods 2.2). Exponential effect on CL (Table 1 footnote, ",
        "covariate model 1: CL = CLpop x exp(beta x (ATU = Positive))), ",
        "exp(0.566) = 1.76-fold higher CL when positive (Results 3.1). ",
        "Estimated from the Prometheus PK database (57 ATU-positive of ",
        "11,372 specimens); none of the 312 study samples were ",
        "ATU-positive (Table 4)."
      ),
      source_name = "ATU"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "ustekinumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ustekinumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "ustekinumab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 65L,
    n_studies = 1L,
    n_observations = "308 HMSA re-assayed serum concentrations (training cohort); ATU effect from 11,372 PK-laboratory specimens",
    disease_state = "Crohn's disease (training cohort); the Bayesian application cohort had Crohn's disease (78/83) or ulcerative colitis (5/83)",
    dose_range = "IV weight-based loading dose then 90 mg SC every 8 weeks (training cohort); 90 mg SC every 4-8 weeks maintenance (application cohort)",
    regions = "Belgium (training cohort, Wang 2022); United States (application cohort and PK-laboratory database)",
    notes = paste0(
      "The structural and albumin parameters were estimated on 308 ",
      "samples from 65 Crohn's disease patients of the Wang 2022 ",
      "training cohort, re-assayed with the HMSA (Methods 2.3). The ATU ",
      "effect was then estimated with all other parameters fixed on a ",
      "clinical PK-laboratory database of 11,372 specimens (mean ",
      "concentration 6.7 ug/mL, 57 ATU-positive). The resulting model ",
      "was applied, with all parameters fixed, to a prospective cohort ",
      "of 83 patients (median age 42 years, IQR 28-57; 52% female; ",
      "312 dose cycles; Table 4). Demographics of the training cohort ",
      "are not reprinted in Yarur 2025."
    )
  )

  ini({
    # Table 1, 'Final Bayesian Priors' column. Every entry is fixed in that
    # column because the paper used the model only as a Bayesian prior for
    # post hoc clearance estimation.
    lka <- fixed(log(0.142)); label("First-order SC absorption rate constant (1/day)")  # Table 1 'Ka (day-1)' = 0.142 (fixed) in all four columns; Methods 2.3 prose says 0.2, see vignette
    lfdepot <- fixed(log(0.735)); label("SC bioavailability (fraction)")  # Table 1 'Bioavailability (F)' = 0.735 (fixed), the HMSA training-cohort estimate
    lcl <- fixed(log(0.161)); label("Clearance at albumin 4.0 g/dL, ATU-negative (L/day)")  # Table 1 'CL_pop (L/day)' = 0.161 (fixed)
    lvc <- fixed(log(3.72)); label("Central volume of distribution (L)")  # Table 1 'V1_pop (L)' = 3.72 (fixed)
    lq <- fixed(log(0.0322)); label("Intercompartmental clearance (L/day)")  # Table 1 'Q (L/day)' = 0.0322 (fixed)
    lvp <- fixed(log(0.393)); label("Peripheral volume of distribution (L)")  # Table 1 'V2 (L)' = 0.393 (fixed)

    e_ada_cl <- fixed(0.566); label("Exponential effect of ATU positivity on CL (unitless)")  # Table 1 'effect of ATU on CL' = 0.566 (fixed); estimated 0.566 (RSE 8.65%) on the PK database
    e_alb_cl <- fixed(-0.887); label("Power exponent of albumin (per 4.0 g/dL) on CL (unitless)")  # Table 1 'effect of ALB on CL' = -0.887 (fixed)

    etalcl ~ fixed(0.0924)  # Table 1 'IIV of CL' = 30.4%, Final Bayesian Priors column; variance = 0.304^2

    propSd <- fixed(0.11); label("Proportional residual error (fraction)")  # Table 1 'Proportional error (%)' = 11.0% (fixed)
    addSd <- fixed(0.495); label("Additive residual error (ug/mL)")  # Table 1 'Additive error (ug/mL)' = 0.495 (fixed)
  })

  model({
    # Albumin enters in g/dL, centred at 4.0 g/dL (Table 1 footnote)
    alb_gdl <- ALB / 10

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (alb_gdl / 4.0)^e_alb_cl * exp(e_ada_cl * ADA_POS)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Bioavailability applies to SC doses into depot; IV loading doses go
    # to central
    f(depot) <- exp(lfdepot)

    # mg / L = ug/mL
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
