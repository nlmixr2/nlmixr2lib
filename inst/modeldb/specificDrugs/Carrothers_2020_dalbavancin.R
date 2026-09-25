Carrothers_2020_dalbavancin <- function() {
  description <- paste(
    "Three-compartment intravenous population PK model for dalbavancin in adults with",
    "acute bacterial skin and skin structure infections or catheter-related bloodstream",
    "infections (pooled phase 2/3 studies VER001-4, VER001-5, VER001-9 and DUR001-303).",
    "Zero-order infusion into the central compartment and first-order elimination.",
    "Power covariate effects: albumin, creatinine clearance and body weight on CL; albumin",
    "and weight on V1; age, albumin and weight on V2; albumin and weight on V3. IIV on CL,",
    "V1, V2 and V3 as two correlated blocks (CL with V2, V1 with V3); no IIV on Q2 or Q3.",
    "Proportional residual error.",
    sep = " "
  )
  reference <- paste(
    "Carrothers TJ, Chittenden JT, Critchley I. Dalbavancin Population Pharmacokinetic",
    "Modeling and Target Attainment Analysis. Clin Pharmacol Drug Dev. 2020;9(1):21-31.",
    "doi:10.1002/cpdd.695",
    sep = " "
  )
  vignette <- "Carrothers_2020_dalbavancin"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "dalbavancin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dalbavancin", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "dalbavancin", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL, V1, V2 and V3, each normalised to 85.5 kg (the printed reference",
        "in the four covariate equations under Table 3). Table 2 pooled median 85.0 kg, range",
        "43-320 kg. Weight replaced the body surface area covariate of the earlier",
        "(2-compartment) dalbavancin popPK model (Discussion)."
      ),
      source_name = "WT"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Supply in canonical SI g/L. The source reports albumin in US-convention g/dL (Table 2",
        "pooled median 3.7 g/dL, range 1.1-5.1) and every covariate equation normalises to",
        "ALB/3.7, so model() converts inline with alb_gdL <- ALB * 0.1 to keep the published",
        "exponents on their original calibration. Power effect (negative exponent) on CL, V1,",
        "V2 and V3."
      ),
      source_name = "ALB"
    ),
    CRCL = list(
      description = "Baseline creatinine clearance, raw mL/min (not BSA-normalised)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL only, normalised to 100 mL/min (printed in the CL equation under",
        "Table 3). Creatinine clearance was tested only on CL 'based on biological",
        "plausibility' (Methods). The paper does not state the estimating equation (none is",
        "named in the main text or the supplement) and does not say the value was",
        "BSA-normalised; Table 2 summarises it in a raw clearance unit (pooled median 113",
        "mL/min, range 22-440; the table header misprints the unit as 'mg/mL'). The 30 mL/min",
        "dose-reduction threshold used throughout the paper is likewise a raw mL/min value.",
        "Few subjects had CLCR < 30 mL/min and the authors caution against simulating below",
        "30 mL/min (Discussion)."
      ),
      source_name = "CLCR"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on V2 only, normalised to 47 years (printed in the V2 equation under",
        "Table 3; equals the Table 2 pooled median of 47.0 years, range 18-93). An age effect",
        "on V3 selected by the stepwise search was removed during refinement (Supplementary",
        "Table S3, Run 1003) and replaced by weight and albumin."
      ),
      source_name = "AGE"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tested on each parameter with IIV in the stepwise covariate search (Methods) but not retained (Supplementary Table S2). Table 2: 59.2% male pooled."
    ),
    RACE_BLACK = list(
      description = "Black race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Race was tested on each parameter with IIV (Methods) but not retained (Supplementary Table S2). Table 2 pooled: Caucasian 71.7%, Hispanic 14.7%, Black 11.1%, Asian 1.0%, Other 1.6%."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 703L,
    n_studies = 4L,
    n_observations = "2310 dalbavancin plasma concentrations (Methods, Data Inclusion)",
    age_range = "18-93 years",
    age_median = "47 years",
    weight_range = "43-320 kg",
    weight_median = "85 kg",
    sex_female_pct = 40.8,
    race_ethnicity = "Caucasian 71.7%, Hispanic 14.7%, Black 11.1%, Asian 1.0%, Other 1.6% (Table 2, pooled)",
    disease_state = paste(
      "Adults with catheter-related bloodstream infection (VER001-4), skin and soft tissue",
      "infection (VER001-5), complicated skin and soft tissue infection (VER001-9) or acute",
      "bacterial skin and skin structure infection (DUR001-303) caused by suspected or",
      "confirmed gram-positive pathogens."
    ),
    renal_function = "Creatinine clearance median 113 mL/min, range 22-440 (Table 2, pooled); albumin median 3.7 g/dL, range 1.1-5.1.",
    dose_range = paste(
      "Intravenous dalbavancin: 1000 mg day 1 + 500 mg day 8 (VER001-4, VER001-9, part of",
      "VER001-5 and DUR001-303); 1100 mg day 1 (VER001-5); 1500 mg single dose on day 1",
      "(DUR001-303). DUR001-303 reduced doses to 1000 mg (single dose) or 750 + 375 mg",
      "(two-dose) for CLCR < 30 mL/min without dialysis. The target-attainment simulations",
      "used a 30-minute infusion."
    ),
    regions = "Multinational phase 2/3 trials (US-led; study sites listed in the online supplemental table).",
    notes = paste(
      "NONMEM 7.3, FOCE. Final model = supplement Run 1005. Plasma dalbavancin measured by",
      "LC-MS/MS (LLOQ 0.5 ug/mL); 2 BLQ/missing samples and 1 outlier were excluded.",
      "Covariate search via PsN stepwise covariate modelling (forward p < 0.01, backward",
      "p < 0.001). Bootstrap uncertainty. Target attainment used fu = 0.07 and daily",
      "average fAUC = fAUC(0-120 h)/5 against a murine-thigh stasis target of fAUC/MIC 27.1."
    )
  )

  ini({
    # Structural parameters: Carrothers 2020 Table 3 (final model), typical values at the
    # reference covariates ALB 3.7 g/dL, CLCR 100 mL/min, WT 85.5 kg, AGE 47 y.
    lcl <- log(0.0531); label("Clearance CL at reference covariates (L/h)") # Table 3 theta1 CL = 0.0531 L/h (RSE 1.1%)
    lvc <- log(3.04); label("Central volume V1 at reference covariates (L)") # Table 3 theta2 V1 = 3.04 L (RSE 4.1%)
    lvp <- log(8.78); label("First peripheral volume V2 at reference covariates (L)") # Table 3 theta3 V2 = 8.78 L (RSE 3.9%)
    lvp2 <- log(3.28); label("Second peripheral volume V3 at reference covariates (L)") # Table 3 theta4 V3 = 3.28 L (RSE 9.6%)
    lq <- log(0.288); label("Intercompartmental clearance Q2, central-V2 (L/h)") # Table 3 theta5 Q2 = 0.288 L/h (RSE 13.2%)
    lq2 <- log(2.11); label("Intercompartmental clearance Q3, central-V3 (L/h)") # Table 3 theta6 Q3 = 2.11 L/h (RSE 10.8%)

    # Covariate power exponents: Table 3 theta7-theta17 (theta15 is not listed; see model() note).
    e_alb_cl <- -0.477; label("Power exponent of albumin on CL (unitless)") # Table 3 theta7 CL-ALB = -0.477
    e_crcl_cl <- 0.273; label("Power exponent of creatinine clearance on CL (unitless)") # Table 3 theta8 CL-CLCR = 0.273
    e_wt_cl <- 0.391; label("Power exponent of body weight on CL (unitless)") # Table 3 theta9 CL-WT = 0.391
    e_alb_vc <- -0.340; label("Power exponent of albumin on V1 (unitless)") # Table 3 theta10 V1-ALB = -0.340
    e_wt_vc <- 0.683; label("Power exponent of body weight on V1 (unitless)") # Table 3 theta11 V1-WT = 0.683
    e_age_vp <- 0.486; label("Power exponent of age on V2 (unitless)") # Table 3 theta12 V2-AGE = 0.486
    e_alb_vp <- -0.413; label("Power exponent of albumin on V2 (unitless)") # Table 3 theta13 V2-ALB = -0.413
    e_wt_vp <- 0.365; label("Power exponent of body weight on V2 (unitless)") # Table 3 theta14 V2-WT = 0.365
    e_alb_vp2 <- -0.551; label("Power exponent of albumin on V3 (unitless)") # Table 3 theta16 V3-ALB = -0.551
    e_wt_vp2 <- 0.518; label("Power exponent of body weight on V3 (unitless)") # Table 3 theta17 V3-WT = 0.518

    # IIV: Table 3 omegas are variances on the log scale (sqrt(0.0489) = 22% CV as printed).
    # Two blocks per supplement Table S3 Run 1005, 'Omega(CL + V2, V1 + V3)'.
    etalcl + etalvp ~ c(0.0489, 0.0823, 0.153) # Table 3 omega1.1 CL = 0.0489, omega2.1 CL,V2 = 0.0823, omega2.2 V2 = 0.153
    etalvc + etalvp2 ~ c(0.0566, 0.111, 0.437) # Table 3 omega3.3 V1 = 0.0566, omega4.3 V1,V3 = 0.111, omega4.4 V3 = 0.437

    # Residual error: Table 3 sigma1.1 proportional variance 0.0362 -> SD sqrt(0.0362) = 0.190.
    propSd <- 0.190; label("Proportional residual error (fraction)") # Table 3 sigma1.1 = 0.0362 (variance)
  })

  model({
    # Albumin is supplied in canonical g/L; the source equations use g/dL (reference 3.7).
    alb_gdL <- ALB * 0.1

    # Covariate equations printed under Table 3. The V3 equation prints its exponents as
    # theta15 (ALB) and theta16 (WT), while Table 3 lists V3-ALB as theta16 and V3-WT as
    # theta17 with no theta15; the values are mapped by NAME (V3-ALB, V3-WT), which is
    # unambiguous. The missing theta15 is consistent with the age-on-V3 effect removed in
    # supplement Run 1003.
    cl <- exp(lcl + etalcl) * (alb_gdL / 3.7)^e_alb_cl * (CRCL / 100)^e_crcl_cl * (WT / 85.5)^e_wt_cl
    vc <- exp(lvc + etalvc) * (alb_gdL / 3.7)^e_alb_vc * (WT / 85.5)^e_wt_vc
    vp <- exp(lvp + etalvp) * (AGE / 47)^e_age_vp * (alb_gdL / 3.7)^e_alb_vp * (WT / 85.5)^e_wt_vp
    vp2 <- exp(lvp2 + etalvp2) * (alb_gdL / 3.7)^e_alb_vp2 * (WT / 85.5)^e_wt_vp2
    q <- exp(lq)
    q2 <- exp(lq2)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(central) <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Total plasma dalbavancin (dose mg / volume L = mg/L = ug/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
