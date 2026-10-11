Carrothers_2023_dalbavancin <- function() {
  description <- paste(
    "Three-compartment intravenous population PK model for dalbavancin in pediatric",
    "patients from birth to <18 years (pooled studies A8841004, DUR001-106, DUR001-107",
    "and DUR001-306). Zero-order infusion into the central compartment and first-order",
    "elimination. Fixed allometric scaling of all clearances (exponent 0.75) and volumes",
    "(exponent 1) to 70 kg; serum albumin acts on every PK parameter through a power",
    "effect on the relative bioavailability of the central compartment; an age-switched",
    "renal-function power effect on CL (maturation eGFR below 2 years, bedside Schwartz",
    "creatinine clearance from 2 years). Correlated IIV on CL, V1 and V2; proportional",
    "residual error with a separate magnitude for the phase 3 study.",
    sep = " "
  )
  reference <- paste(
    "Carrothers TJ, Lagraauw HM, Lindbom L, Riccobene TA. Population Pharmacokinetic and",
    "Pharmacokinetic/Pharmacodynamic Target Attainment Analyses for Dalbavancin in",
    "Pediatric Patients. Pediatr Infect Dis J. 2023;42(2):99-105.",
    "doi:10.1097/INF.0000000000003764.",
    "Covariate reference values and baseline demographics from the FDA Multi-disciplinary",
    "Review and Evaluation, NDA 021883/S-010 (Dalvance pediatric efficacy supplement,",
    "Reference ID 4829661, 2021).",
    sep = " "
  )
  vignette <- "Carrothers_2023_dalbavancin"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "dalbavancin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dalbavancin", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "dalbavancin", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Fixed allometric scaling to a 70 kg reference: exponent 0.75 on CL, Q and Q2,",
        "exponent 1 on V1, V2 and V3 (Results; Supplemental Digital Content 2 footnote).",
        "Pooled median 26.4 kg, range 2.6-105 kg (FDA review Table 29)."
      ),
      source_name = "WT"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Supply in canonical SI g/L. The source reports albumin in g/dL; model() converts",
        "inline (ALB * 0.1) so the published exponent keeps its calibration. Power effect",
        "on the central-compartment relative bioavailability F1 (Supplemental Digital",
        "Content 2 theta8), which scales all apparent clearances and volumes together.",
        "Centred on 4.4 g/dL, the pooled median of the analysis population (FDA review",
        "Table 29; range 1.9-5.3 g/dL). Neither the paper nor the review prints the",
        "functional form or the centring value; see the vignette Assumptions."
      ),
      source_name = "ALB"
    ),
    CRCL = list(
      description = "Age-appropriate size-normalised renal function (mL/min/1.73 m^2)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Piece-wise by AGE, as in the source. For AGE >= 2 years supply the bedside",
        "Schwartz creatinine clearance, 0.413 * height (cm) / serum creatinine (mg/dL)",
        "(source column CLCRN; Methods ref 25). For AGE < 2 years supply the Rhodin 2009",
        "postmenstrual-age maturation eGFR, 121.2 * PMA^3.4 / (47.7^3.4 + PMA^3.4) with PMA",
        "in weeks (source column EGFR; Methods ref 26; the size term is carried by the",
        "allometric WT scaling of CL). Centring values are the age-stratum medians the FDA",
        "reviewer reported: 74.06 (Rhodin, <2 y) and 100.88 (Schwartz, >=2 y) mL/min/1.73 m^2.",
        "Patients with Schwartz CLcr < 30 mL/min/1.73 m^2 were excluded from all studies."
      ),
      source_name = "CLCRN / EGFR"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Selects the renal-function descriptor and exponent on CL: AGE < 2 uses the",
        "maturation eGFR exponent (theta9), AGE >= 2 the bedside Schwartz exponent",
        "(theta10). No other age effect. Pooled range 0.011-17.9 years (FDA review",
        "Table 29)."
      ),
      source_name = "AGE"
    ),
    STUDY_PHASE3 = list(
      description = "Phase 3 study DUR001-306 indicator (1 = DUR001-306, 0 = phase 1 studies)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase 1 studies A8841004, DUR001-106, DUR001-107)",
      notes = paste(
        "Selects the residual-error magnitude only: proportional SD 0.123 in the phase 1",
        "studies and 0.173 in DUR001-306 (Supplemental Digital Content 2 sigma1.1 and",
        "sigma3.3). It touches no structural or covariate parameter; simulate with 0 for",
        "rich-profile predictions."
      ),
      source_name = "STUDY"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL, V1, V2, V3, Q2 and Q3 in the stepwise covariate search but not retained (FDA review Table 31; Results)."
    ),
    RACE_BLACK = list(
      description = "Black race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Race was tested on CL, V1, V2, V3, Q2 and Q3 in the stepwise covariate search but not retained (FDA review Table 31; Results). Pooled: White 84.8%, Black 10.0%."
    )
  )

  # The paper estimated one residual SD for the phase 1 studies and another for
  # the phase 3 study; nlmixr2 takes a single residual-SD symbol per output, so
  # the two are combined in model() with the STUDY_PHASE3 indicator.
  paper_specific_residual_sds <- c("propSdPh1", "propSdPh3")

  population <- list(
    species = "human",
    n_subjects = 211L,
    n_studies = 4L,
    n_observations = "1124 dalbavancin plasma concentrations (Methods)",
    age_range = "4 days to 17.9 years (0.011-17.9 years)",
    age_median = "8.5 years",
    weight_range = "2.6-105 kg",
    weight_median = "26.4 kg",
    sex_female_pct = 36.5,
    race_ethnicity = "White 84.8%, Black 10.0%, Other 2.4%, American Indian or Alaska Native 1.9%, Asian 0.9% (FDA review Table 29)",
    disease_state = paste(
      "Pediatric patients with acute bacterial skin and skin structure infection",
      "(DUR001-306), hospitalized children on standard anti-infective treatment for",
      "bacterial infection (A8841004, DUR001-106), and preterm or term neonates and",
      "infants < 3 months with suspected or confirmed bacterial infection (DUR001-107)."
    ),
    renal_function = "Serum creatinine median 0.52 mg/dL (range 0.13-1.29); Schwartz CLcr < 30 mL/min/1.73 m^2 excluded.",
    albumin = "Serum albumin median 4.4 g/dL (range 1.9-5.3) (FDA review Table 29).",
    dose_range = paste(
      "Single 30-minute intravenous infusions of 10-25 mg/kg (capped at 1000 or 1500 mg)",
      "or 1000 mg; DUR001-306 also gave a two-dose regimen (day 1 and day 8, 15 + 7.5 or",
      "12 + 6 mg/kg). Supplemental Digital Content 1."
    ),
    regions = "Multinational (US, Europe).",
    notes = paste(
      "NONMEM 7.4.0, FOCE with interaction; PsN 4.2.0 stepwise covariate search. Model",
      "updates the earlier pediatric model of Gonzalez et al. (2017) with studies",
      "DUR001-107 and DUR001-306. Target attainment used fu = 0.07 and daily average",
      "fAUC = fAUC(0-120 h)/5."
    )
  )

  ini({
    # Structural parameters: Supplemental Digital Content 2 (identical to FDA review
    # Table 32), typical values for a 70 kg patient at the reference covariates.
    lcl <- log(0.0578); label("Clearance CL for a 70 kg patient (L/h)") # SDC 2 theta1 CL = 0.0578 L/h
    lvc <- log(4.58); label("Central volume V1 for a 70 kg patient (L)") # SDC 2 theta2 V = 4.58 L
    lvp <- log(6.1); label("First peripheral volume V2 for a 70 kg patient (L)") # SDC 2 theta3 V2 = 6.1 L
    lq <- log(0.794); label("Intercompartmental clearance Q, central-V2, for a 70 kg patient (L/h)") # SDC 2 theta4 Q = 0.794 L/h
    lq2 <- log(0.00996); label("Intercompartmental clearance Q2, central-V3, for a 70 kg patient (L/h)") # SDC 2 theta5 Q2 = 0.00996 L/h
    lvp2 <- log(5.57); label("Second peripheral volume V3 for a 70 kg patient (L)") # SDC 2 theta6 V3 = 5.57 L

    # Typical relative bioavailability of the IV dose, the anchor for the albumin effect.
    lfcentral <- fixed(log(1)); label("Typical relative bioavailability F1 of the IV dose (fraction)") # FDA review p.92 'ALB was parameterized via F (fixed to 1)'; the theta7 gap in SDC 2

    # Allometric exponents fixed a priori.
    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT/70) on CL, Q and Q2 (unitless)") # Methods / SDC 2 footnote 'fixed exponents of 0.75 for all clearances'
    e_wt_vc <- fixed(1); label("Allometric exponent of (WT/70) on V1, V2 and V3 (unitless)") # Methods / SDC 2 footnote 'and 1 for all volumes'

    # Covariate exponents.
    e_alb_fcentral <- 0.385; label("Power exponent of albumin (g/dL / 4.4) on F1 (unitless)") # SDC 2 theta8 F1,ALB = 0.385
    e_crcl_cl_lt2y <- 0.167; label("Power exponent of maturation eGFR (/74.06) on CL, age < 2 years (unitless)") # SDC 2 theta9 CL,eGFR = 0.167
    e_crcl_cl <- 0.0681; label("Power exponent of Schwartz CLcr (/100.88) on CL, age >= 2 years (unitless)") # SDC 2 theta10 CL,CrCLN = 0.0681

    # IIV: SDC 2 reports SDs (CV% = sqrt(exp(SD^2) - 1)) and correlations; converted
    # here to a variance-covariance block. Variances 0.319^2, 0.454^2, 0.321^2;
    # covariances corr * SD_i * SD_j.
    etalcl + etalvc + etalvp ~ c(
      0.101761,
      0.118902, 0.206116,
      0.087449, 0.099536, 0.103041
    ) # SDC 2 omega1.1 0.319, omega2.2 0.454, omega3.3 0.321 (SDs); corr CL-V 0.821, CL-V2 0.854, V-V2 0.683

    # Residual error: proportional, separate magnitude for the phase 3 study.
    propSdPh1 <- 0.123; label("Proportional residual SD, phase 1 studies (fraction)") # SDC 2 sigma1.1 PropErr = 0.123 (SD)
    propSdPh3 <- 0.173; label("Proportional residual SD, phase 3 study DUR001-306 (fraction)") # SDC 2 sigma3.3 PropErr-PhIII = 0.173 (SD)
  })

  model({
    # Albumin is supplied in canonical g/L; the source works in g/dL.
    alb_gdL <- ALB * 0.1

    # Age-switched renal function on CL: maturation eGFR below 2 years, bedside
    # Schwartz CLcr from 2 years, each centred on its stratum median.
    lt2y <- (AGE < 2)
    renal_cl <- lt2y * (CRCL / 74.06)^e_crcl_cl_lt2y + (1 - lt2y) * (CRCL / 100.88)^e_crcl_cl

    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * renal_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc
    vp2 <- exp(lvp2) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_cl
    q2 <- exp(lq2) * (WT / 70)^e_wt_cl

    # Albumin enters as a relative bioavailability of the IV dose, which divides
    # every apparent clearance and volume by the same factor.
    fcentral <- exp(lfcentral) * (alb_gdL / 4.4)^e_alb_fcentral

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(central) <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    f(central) <- fcentral

    # Total plasma dalbavancin (dose mg / volume L = mg/L = ug/mL).
    Cc <- central / vc

    propSdCc <- propSdPh3 * STUDY_PHASE3 + propSdPh1 * (1 - STUDY_PHASE3)
    Cc ~ prop(propSdCc)
  })
}
