Saito_2019_ipragliflozin <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption for",
    "oral ipragliflozin in Japanese adults with type 2 diabetes mellitus",
    "(Saito 2019). Apparent clearance CL/F is a product of power functions",
    "of absolute (BSA-de-normalised) GFR, total serum protein, total",
    "bilirubin and body surface area; IIV on CL/F only. Companion",
    "exposure-response model for 24 h urinary glucose excretion:",
    "Saito_2019_ipragliflozin_uge."
  )
  reference <- paste(
    "Saito M, Kaibara A, Kadokura T, Toyoshima J, Yoshida S, Kazuta K,",
    "Ueyama E. Pharmacokinetic and pharmacodynamic modelling for renal",
    "function dependent urinary glucose excretion effect of ipragliflozin,",
    "a selective sodium-glucose cotransporter 2 inhibitor, both in healthy",
    "subjects and patients with type 2 diabetes mellitus.",
    "Br J Clin Pharmacol. 2019;85(8):1808-1819. doi:10.1111/bcp.13972"
  )
  vignette <- "Saito_2019_ipragliflozin"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description = "Absolute glomerular filtration rate: Japanese-modified MDRD eGFR de-normalised by the individual Du Bois body surface area",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Saito 2019 Results 'Demographics': eGFR (mL/min/1.73 m^2) from the",
        "MDRD equation modified for Japanese CKD patients, then 'GFR (mL/min)",
        "corrected by individual BSA was used for modelling' -- i.e. raw",
        "mL/min, NOT BSA-normalised (Table 1 footnote; Table 1 GFR row is",
        "in mL/min). Supply eGFR x BSA / 1.73. Reference 90 mL/min",
        "(Equation 4). Power term (CRCL / 90)^0.233 on CL/F."
      ),
      source_name = "GFR"
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source reports g/dL (Table 1; reference 7.0 g/dL in Equation 4).",
        "The model converts the canonical g/L column inline (TPRO / 10)",
        "before applying (TPRO_gdL / 7.0)^-0.417. Baseline value per the",
        "covariate screen, which also used values at each assessment visit."
      ),
      source_name = "TPRO"
    ),
    TBILI = list(
      description = "Total serum bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source reports mg/dL (Table 1; reference 0.8 mg/dL in Equation 4).",
        "The model converts the canonical umol/L column inline (TBILI / 17.1)",
        "before applying (TBIL_mgdL / 0.8)^-0.0681."
      ),
      source_name = "TBIL"
    ),
    BSA = list(
      description = "Body surface area (Du Bois)",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Du Bois equation (Saito 2019 Results 'Demographics'). Reference 1.73 m^2 (Equation 4). Power term (BSA / 1.73)^0.610 on CL/F.",
      source_name = "BSA"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "ipragliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ipragliflozin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ipragliflozin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 673L,
    n_studies = 6L,
    studies = paste(
      "Prior (base) model: rich sampling in clinical pharmacology studies",
      "B (CL-0070, NCT01023945) and C (CL-0073, NCT01097681; 534",
      "concentrations from 43 patients). Final model: trough",
      "concentrations from phase II/III studies D (CL-0103), E (CL-0105),",
      "F (CL-0121) and G (CL-0072; 3714 concentrations from 630 patients)",
      "fitted with NONMEM $PRIOR on the base model."
    ),
    age_range = "26-86 years (Table 1, T2DM total)",
    weight_range = "41.5-128.0 kg (Table 1, T2DM total; mean 68.19 kg)",
    sex_female_pct = 31.7,
    race_ethnicity = "Japanese",
    disease_state = "Type 2 diabetes mellitus; normal renal function to moderate renal impairment (eGFR >= 30 mL/min/1.73 m^2)",
    dose_range = "Ipragliflozin 12.5-100 mg once daily orally",
    regions = "Japan",
    renal_function = "GFR 24.1-181.5 mL/min (Table 1, T2DM total; mean 84.46); 35.9% normal, 52.5% mild, 11.6% moderate impairment by eGFR class",
    notes = "Demographics from Saito 2019 Table 1 (T2DM total, n = 887 including 192 placebo patients; 673 contributed PK data)."
  )

  ini({
    lka <- log(6.38)
    label("First-order absorption rate constant ka (1/h)") # Table 3 'Ka (h-1)' = 6.38
    lcl <- log(9.47)
    label("Apparent clearance CL/F at the reference covariates (L/h)") # Table 3 'CL (L/h)' = 9.47; Equation 4
    lvc <- log(39.4)
    label("Apparent central volume Vc/F (L)") # Table 3 'Vc (L)' = 39.4
    lq <- log(6.63)
    label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 'Q/F (L/h)' = 6.63
    lvp <- log(68.1)
    label("Apparent peripheral volume Vp/F (L)") # Table 3 'Vp (L)' = 68.1

    e_crcl_cl <- 0.233
    label("Power exponent of GFR (mL/min, ref 90) on CL/F (unitless)") # Table 3 'GFR effect on CL' = 0.233; Equation 4
    e_tpro_cl <- -0.417
    label("Power exponent of total protein (g/dL, ref 7.0) on CL/F (unitless)") # Table 3 'TPRO effect on CL' = -0.417; Equation 4
    e_tbili_cl <- -0.0681
    label("Power exponent of total bilirubin (mg/dL, ref 0.8) on CL/F (unitless)") # Table 3 'TBIL effect on CL' = -0.0681; Equation 4
    e_bsa_cl <- 0.610
    label("Power exponent of body surface area (m^2, ref 1.73) on CL/F (unitless)") # Table 3 'BSA effect on CL' = 0.610; Equation 4

    etalcl ~ 0.0533 # Table 3 'omega2: CL' = 0.0533 (CV 23.4%)

    propSd <- 0.2441
    label("Proportional residual error (fraction)") # Table 3 sigma2 = 0.0596 -> sqrt(0.0596) = 0.2441 (paper CV 24.8% uses sqrt(exp(sigma2)-1))
  })
  model({
    # Canonical SI covariate columns converted to the paper's US units
    tpro_gdl <- TPRO / 10
    tbili_mgdl <- TBILI / 17.1

    # Equation 4
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (CRCL / 90)^e_crcl_cl *
      (tpro_gdl / 7.0)^e_tpro_cl *
      (tbili_mgdl / 0.8)^e_tbili_cl *
      (BSA / 1.73)^e_bsa_cl
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L = ug/mL; x 1000 -> ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
