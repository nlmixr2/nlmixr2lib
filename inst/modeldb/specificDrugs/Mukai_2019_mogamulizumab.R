Mukai_2019_mogamulizumab <- function() {
  description <- "Two-compartment population PK model for mogamulizumab in adults with cutaneous T-cell lymphoma or adult T-cell lymphoma (Mukai 2019)"
  reference <- "Mukai M, Mould DR, Nishimura K, Gallerani E, Grimwood D. Population Pharmacokinetic Modeling of Mogamulizumab in Adults With Cutaneous T-Cell Lymphoma or Adult T-Cell Lymphoma. J Clin Pharmacol. 2020;60(1):58-66. doi:10.1002/jcph.1564"
  vignette <- "Mukai_2019_mogamulizumab"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    ALB = list(
      description        = "Serum albumin concentration",
      units = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline or time-varying. Inversely related to CL and V2. Reference value 4.0 g/dL.",
      source_name        = "ALB"
    ),
    AST = list(
      description        = "Aspartate aminotransferase activity",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline or time-varying. Positively related to CL. Reference value 24 U/L.",
      source_name        = "AST"
    ),
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Computed from height and weight (method unspecified in source). Positively related to V1. Reference value 1.82 m^2.",
      source_name        = "BSA"
    ),
    HEPIMP = list(
      description        = "Hepatic impairment indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function)",
      notes              = "Mild and moderate hepatic impairment combined (83 patients: 80 mild + 3 moderate) vs normal (361 patients). No severe cases. Mild-to-moderate impairment associated with 14% increase in CL.",
      source_name        = "HI"
    ),
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Female subjects (204, 45.9%) had approximately 23% lower CL than males (240, 54.1%).",
      source_name        = "SEX"
    )
  )

  population <- list(
    n_subjects     = 444,
    n_studies      = 6,
    age_range      = "22-101 years",
    age_mean       = "61.4 years (SD 13.0)",
    weight_range   = "36.0-149.7 kg",
    weight_mean    = "73.6 kg (SD 18.5)",
    sex_female_pct = 45.9,
    race_ethnicity = c(
      White = 49.5,
      Asian = 28.6,
      Black = 18.0,
      "Native American or Alaskan" = 0.2,
      Other = 1.6,
      Unknown = 11.9
    ),
    disease_state  = "Relapsed/refractory or previously untreated adult T-cell lymphoma (ATL, 29.1%) or cutaneous T-cell lymphoma (CTCL, 70.9%). ATL subtypes: acute (56.6%), lymphoma (28.7%), chronic (14.7%). CTCL subtypes: mycosis fungoides (54.9%), Sézary syndrome (45.1%).",
    dose_range     = "0.01-1.0 mg/kg IV; 98% received 1.0 mg/kg",
    regions        = "Studies 0761-0501, 0761-002, 0761-003, 0761-004 in Japan; 0761-009 (phase 2, global); 0761-010 (phase 3, global including US, Europe, Japan)",
    notes          = "Pooled from 6 phase 1-3 clinical trials. Dosing regimens varied: once weekly × 4-8 weeks, then every 2 weeks until progression for studies 0761-009 and 0761-010. ECOG PS 0-2. Concomitant mLSG15 chemotherapy in 29 patients (6.5%). Renal impairment: normal 46.2%, mild 35.4%, moderate 18.0%, severe 0.5%. Hepatic impairment: normal 81.3%, mild/moderate 18.7%, severe 0%."
  )

  ini({
    # Structural PK parameters (log-transformed)
    lcl <- log(0.0138); label("Clearance at reference covariates (L/hr)") # Table 5
    lvc <- log(3.65);   label("Central volume of distribution at reference covariates (L)") # Table 5, V1
    lvp <- log(2.48);   label("Peripheral volume of distribution at reference covariates (L)") # Table 5, V2
    lq  <- log(0.0532); label("Intercompartmental clearance (L/hr)") # Table 5, Q

    # Covariate effects on CL
    e_alb_cl  <- -1.81; label("ALB effect on CL (power exponent, reference 4.0 g/dL)") # Table 5
    e_ast_cl  <-  0.282; label("AST effect on CL (power exponent, reference 24 U/L)") # Table 5
    e_hepimp  <-  1.14; label("Hepatic impairment effect on CL (multiplicative factor for mild/moderate vs normal)") # Table 5, HI on CL
    e_sexf    <-  0.768; label("Female sex effect on CL (multiplicative factor vs male reference)") # Table 5, SEX on CL

    # Covariate effect on V1 (central volume)
    e_bsa_vc  <-  0.884; label("BSA effect on V1 (power exponent, reference 1.82 m^2)") # Table 5, BSA on V1

    # Covariate effect on V2 (peripheral volume)
    e_alb_vp  <- -1.47; label("ALB effect on V2 (power exponent, reference 4.0 g/dL)") # Table 5, ALB on V2

    # Inter-individual variability (log-normal; omega^2 = log(CV^2 + 1))
    # Table 5 reports SD%: CL 50.5%, V1 20.5%, V2 72.7%
    # omega^2: CL = log(1 + 0.505^2) = 0.227, V1 = log(1 + 0.205^2) = 0.0411, V2 = log(1 + 0.727^2) = 0.424
    etalcl ~ 0.227  # 50.5% CV
    etalvc ~ 0.0411 # 20.5% CV
    etalvp ~ 0.424  # 72.7% CV

    # Residual error (log-additive, equivalent to proportional on back-transform)
    propSd <- 0.261; label("Proportional residual error (CV)") # Table 5, residual error
  })

  model({
    # SI -> US-convention unit conversion (canonical ALB is in SI g/L per the
    # 2026-06-19 register standardization audit; the original calibration
    # used the g/dL reference value, so convert inline here).
    alb_gdL <- ALB * 0.1  # SI g/L -> US-convention g/dL (factor 0.1)

    # Reference covariate values (Table 6)
    alb_ref <- 4.0   # g/dL
    ast_ref <- 24    # U/L
    bsa_ref <- 1.82  # m^2

    # Covariate effects
    alb_cl_effect    <- (alb_gdL / alb_ref)^e_alb_cl
    ast_cl_effect    <- (AST / ast_ref)^e_ast_cl
    hepimp_cl_effect <- 1 + (e_hepimp - 1) * HEPIMP  # 1 when HEPIMP=0, 1.14 when HEPIMP=1
    sex_cl_effect    <- 1 + (e_sexf - 1) * SEXF      # 1 when SEXF=0 (male), 0.768 when SEXF=1 (female)
    bsa_vc_effect    <- (BSA / bsa_ref)^e_bsa_vc
    alb_vp_effect    <- (alb_gdL / alb_ref)^e_alb_vp

    # Individual PK parameters
    cl <- exp(lcl + etalcl) * alb_cl_effect * ast_cl_effect * hepimp_cl_effect * sex_cl_effect
    vc <- exp(lvc + etalvc) * bsa_vc_effect
    vp <- exp(lvp + etalvp) * alb_vp_effect
    q  <- exp(lq)

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment model with IV dosing to central
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Concentration (dose in mg, volume in L -> mg/L = ug/mL)
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
