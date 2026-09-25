Shemesh_2019_atezolizumab <- function() {
  description <- paste(
    "Two-compartment population PK model with intravenous infusion input",
    "for atezolizumab (anti-PD-L1 IgG1) in pediatric and young adult",
    "patients (7 months to 29 years, 8.7-154 kg) with relapsed or refractory",
    "solid tumors or lymphoma from the iMATRIX-atezolizumab phase I/II study",
    "(Shemesh 2019, n = 87, 431 serum concentrations). The adult",
    "atezolizumab popPK structure was refitted to the pediatric data with",
    "every parameter re-estimated: power effects of body weight and albumin",
    "on CL and V1, of baseline tumor burden on CL, and a multiplicative",
    "factor for treatment-emergent ADA on CL. Sex effects were not evaluated",
    "and V2 and Q carry no covariates."
  )
  reference <- paste(
    "Shemesh CS, Chanu P, Jamsen K, Wada R, Rossato G, Donaldson F, Garg A,",
    "Winter H, Ruppel J, Wang X, Bruno R, Jin J, Girish S.",
    "Population pharmacokinetics, exposure-safety, and immunogenicity of",
    "atezolizumab in pediatric and young adult patients with cancer.",
    "J Immunother Cancer. 2019;7:314. doi:10.1186/s40425-019-0791-x."
  )
  vignette <- "Shemesh_2019_atezolizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "atezolizumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "atezolizumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effects on CL (exponent 0.795) and V1 (exponent 0.766),",
        "normalized to 77 kg. The normalization reference is printed only",
        "in the adult CL and V1 equations of Shemesh 2019 Methods ('PopPK",
        "model'); the pediatric model 'utilized the same structure' with",
        "every parameter re-estimated, so the adult reference is retained.",
        "Pediatric cohort range 8.7-154 kg (Table 1)."
      ),
      source_name = "BWT"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effects on CL (exponent -1.18) and V1 (exponent -0.566),",
        "normalized to 40 g/L (reference printed in the Shemesh 2019 adult",
        "equations). Missing values imputed to the median (Table 1",
        "footnote)."
      ),
      source_name = "ALBU"
    ),
    TUMSZ = list(
      description = "Baseline tumor burden (sum of longest diameters of target lesions)",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL only (exponent 0.122), normalized to 63 mm",
        "(reference printed in the Shemesh 2019 adult CL equation).",
        "Missing values (up to 10 of 38 adolescents) imputed to the median",
        "(Table 1 footnote)."
      ),
      source_name = "TUM"
    ),
    ADA_POS = list(
      description = "Treatment-emergent anti-drug-antibody status (1 = ADA-positive, 0 = ADA-negative)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ADA-negative)",
      notes = paste(
        "Multiplicative factor 1.23^ADA_POS on CL. Shemesh 2019 Methods:",
        "'Proportional changes (i.e., for ADA and sex) in the pediatric and",
        "young adult model were parameterized theta^cov ... (both coded 0",
        "or 1)'; Discussion confirms '23% increase in pediatric patients'.",
        "Post-baseline ADA status, treated as a per-subject indicator;",
        "missing ADA records (10 patients) were imputed as ADA-negative."
      ),
      source_name = "ADA"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Present in the adult model (0.871 on V1, 0.728 on V2) but 'NE'",
        "(not evaluated) in the pediatric and young adult model (Shemesh",
        "2019 Table 2); Results: 'Sex effects had minimal impact on the",
        "objective function.' Not part of this model."
      )
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = paste(
        "Not a model covariate. Examined graphically against eta.CL and",
        "eta.V1 (Fig. 1d, no trend); a sensitivity analysis adding weight",
        "and age on V2 and Q was not retained."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 87,
    n_studies = 1,
    n_observations = 431,
    age_range = "7 months to 29 years",
    age_median = "12 years (69 patients < 18 years); 22 years (18 young adults)",
    weight_range = "8.7-154 kg",
    weight_median = "38.9 kg (< 18 years); 61.0 kg (>= 18 years)",
    sex_female_pct = 46,
    disease_state = paste(
      "Relapsed or refractory pediatric solid tumors (Ewing sarcoma,",
      "neuroblastoma, osteosarcoma, rhabdomyosarcoma, non-rhabdomyosarcoma",
      "soft tissue sarcoma, Wilms tumor, rhabdoid tumors), Hodgkin and",
      "non-Hodgkin lymphoma, and other rare tumors"
    ),
    dose_range = paste(
      "15 mg/kg IV q3w (maximum 1200 mg) for patients < 18 years; 1200 mg IV",
      "q3w flat for patients >= 18 years; 60-min infusion in cycle 1 and",
      "30-min infusions thereafter"
    ),
    age_groups = "2 infants < 2 y, 29 children 2 to < 12 y, 38 adolescents 12 to < 18 y, 18 young adults >= 18 y",
    ada_positive_pct = 13,
    regions = "Multinational (iMATRIX-atezolizumab, NCT02541604)",
    notes = "Shemesh 2019 Table 1 (baseline demographics by age group); 40 of 87 female."
  )

  ini({
    # Structural parameters at WT 77 kg, ALB 40 g/L, TUMSZ 63 mm, ADA-negative.
    lcl <- log(0.217); label("Clearance at reference covariates (L/day)")                  # Shemesh 2019 Table 2, CL (L/day) 0.217
    lvc <- log(3.01);  label("Central volume at reference covariates (L)")                 # Shemesh 2019 Table 2, V1 (L) 3.01
    lvp <- log(1.36);  label("Peripheral volume (L)")                                      # Shemesh 2019 Table 2, V2 (L) 1.36
    lq  <- log(0.183); label("Intercompartmental clearance (L/day)")                       # Shemesh 2019 Table 2, Q (L/day) 0.183

    e_wt_cl    <- 0.795; label("Power exponent of WT on CL (unitless)")                    # Shemesh 2019 Table 2, Weight on CL
    e_alb_cl   <- -1.18; label("Power exponent of ALB on CL (unitless)")                   # Shemesh 2019 Table 2, Albumin on CL
    e_tumsz_cl <- 0.122; label("Power exponent of TUMSZ on CL (unitless)")                 # Shemesh 2019 Table 2, Tumor burden on CL
    e_ada_cl   <- 1.23;  label("Multiplicative factor on CL for ADA-positive (unitless)")  # Shemesh 2019 Table 2, Positive ADA on CL; theta^cov per Methods
    e_wt_vc    <- 0.766; label("Power exponent of WT on V1 (unitless)")                    # Shemesh 2019 Table 2, Weight on V1
    e_alb_vc   <- -0.566; label("Power exponent of ALB on V1 (unitless)")                  # Shemesh 2019 Table 2, Albumin on V1

    # BSV variances (log scale); CL-V1 covariance = 0.510 * sqrt(0.0458 * 0.0140)
    etalcl + etalvc ~ c(0.0458, 0.012914, 0.0140) # Shemesh 2019 Table 2, BSV CL 0.0458, BSV V1 0.0140, correlation 0.510
    etalvp ~ 0.311                                # Shemesh 2019 Table 2, BSV for V2 0.311

    # Residual variances 0.051 and 68.9 (ug/mL)^2 reported; SDs are square roots
    propSd <- 0.2258; label("Proportional residual error (fraction)")                      # Shemesh 2019 Table 2, proportional residual variance 0.051 -> sqrt = 0.2258
    addSd  <- 8.301;  label("Additive residual error (ug/mL)")                             # Shemesh 2019 Table 2, additive residual variance 68.9 -> sqrt = 8.301
  })
  model({
    cl <- exp(lcl + etalcl) *
      (WT / 77)^e_wt_cl *
      (ALB / 40)^e_alb_cl *
      (TUMSZ / 63)^e_tumsz_cl *
      e_ada_cl^ADA_POS
    vc <- exp(lvc + etalvc) *
      (WT / 77)^e_wt_vc *
      (ALB / 40)^e_alb_vc
    vp <- exp(lvp + etalvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L = ug/mL
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
