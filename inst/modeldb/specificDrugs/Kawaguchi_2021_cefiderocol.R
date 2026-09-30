Kawaguchi_2021_cefiderocol <- function() {
  description <- "Three-compartment population PK model for intravenous cefiderocol in uninfected subjects (healthy volunteers and subjects spanning normal renal function to end-stage renal disease) and patients with pneumonia, bloodstream infection/sepsis, or complicated urinary tract infection, with Cockcroft-Gault creatinine clearance (capped at 150 mL/min) and infection site on CL, body weight on V1 and V2, and serum albumin and any active infection on V1"
  reference <- paste(
    "Kawaguchi N, Katsube T, Echols R, Wajima T.",
    "Population pharmacokinetic and pharmacokinetic/pharmacodynamic analyses",
    "of cefiderocol, a parenteral siderophore cephalosporin, in patients with",
    "pneumonia, bloodstream infection/sepsis, or complicated urinary tract",
    "infection.",
    "Antimicrob Agents Chemother. 2021;65(3):e01437-20.",
    "doi:10.1128/AAC.01437-20",
    sep = " "
  )
  vignette <- "Kawaguchi_2021_cefiderocol"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance calculated by the Cockcroft-Gault equation.",
        "RAW mL/min and NOT BSA-normalized."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters CL as a power function (CRCL / 83.0)^0.682 up to 150 mL/min",
        "and is held constant at (150 / 83.0)^0.682 above it (Table 2",
        "footnote b; Table S2 $PK). The 150 mL/min cutoff was selected from",
        "a visual inspection of CL against CRCL and then confirmed against",
        "120 and 180 mL/min by objective function; a power-plus-linear model",
        "had given a slope above 150 mL/min of less than 0.0001.",
        "The reference 83.0 mL/min is the centring constant printed in the",
        "footnote equation and in the control stream; it is not the median of",
        "any Table 1 column (the per-study medians are 121.0, 83.0, 69.0 and",
        "73.0 mL/min), and it happens to equal the phase 2 APEKS-cUTI study",
        "median.",
        "Time-varying in the source data set for the phase 3 patients:",
        "Table S2 reads the baseline column CLCR for subjects with ID < 2000",
        "(phase 1 and phase 2) and the time-varying column TCLCR for subjects",
        "with ID >= 2000 (phase 3 CREDIBLE-CR and APEKS-NP). Supply a",
        "per-record value when renal function changes over the course of",
        "treatment.",
        "Observed range 4-306 mL/min (Table 1). Patients on hemodialysis in",
        "the phase 3 studies (n = 32) were excluded from the analysis."
      ),
      source_name = "CLCR / TCLCR"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters V1 and V2 with ONE shared power exponent, 0.580, centred on",
        "72.6 kg (Table 2, 'Effect of body weight on V 1 and V 2'; Table S2",
        "WT1 = (WT/72.6)**THETA(9) multiplies both TVV1 and TVV2). The 72.6",
        "kg centre is printed in the footnote equation; the per-study medians",
        "in Table 1 run 68.4-76.4 kg. Observed range 25.0-156.0 kg."
      ),
      source_name = "WT"
    ),
    ALB = list(
      description = "Serum albumin concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The paper reports albumin in g/dL and the model was calibrated in",
        "g/dL with a reference of 3.9 g/dL: V1 is multiplied by",
        "(ALB / 3.9)^-0.617 (Table 2 footnote b; Table S2",
        "ALB2 = (ALB/3.9)**THETA(14)). The canonical column is in g/L, so",
        "model() converts it with alb_gdL <- ALB * 0.1 before applying the",
        "published equation. Supply 39 g/L for a patient at the reference.",
        "Lower albumin gives a LARGER V1. Observed range 1.2-5.3 g/dL; the",
        "phase 3 medians were 3.0 (APEKS-NP) and 2.7 (CREDIBLE-CR) g/dL",
        "against 4.2 g/dL in the phase 1 and phase 2 studies (Table 1)."
      ),
      source_name = "ALB"
    ),
    DIS_PNEUMONIA = list(
      description = paste(
        "Index infection is pneumonia: 1 = patient enrolled with pneumonia in",
        "the phase 3 APEKS-NP or CREDIBLE-CR study; 0 = otherwise."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all four infection indicators 0 = subject without infection, the phase 1 cohort)",
      notes = paste(
        "Source coding PT = 2 in Table S2. The enrolled pneumonia was",
        "hospital-acquired, ventilator-associated or health care-associated",
        "(Materials and Methods, Data for analyses), but the model carries",
        "ONE pneumonia level with no acquisition-setting split and the",
        "registry has no health care-associated pneumonia column, so the",
        "general DIS_PNEUMONIA flag is used rather than DIS_HABP + DIS_VABP.",
        "Mutually exclusive with DIS_BACTEREMIA and DIS_CUTI.",
        "Multiplies CL by 0.981 (Table 2) and, as one of the infection",
        "indicators, V1 by 1.39. Mechanical ventilation was screened and was",
        "not retained."
      ),
      source_name = "PT = 2"
    ),
    DIS_BACTEREMIA = list(
      description = paste(
        "Index infection is bloodstream infection/sepsis: 1 = patient",
        "enrolled with BSI/sepsis in the phase 3 CREDIBLE-CR study; 0 =",
        "otherwise."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all four infection indicators 0 = subject without infection, the phase 1 cohort)",
      notes = paste(
        "Source coding PT = 3 in Table S2 ('BSI in Ph3 CR'). The paper's",
        "category is 'bloodstream infection/sepsis', i.e. the protocol",
        "pooled sepsis with bacteremia; there is no separate sepsis",
        "coefficient, so DIS_SEPSIS is not used. Only CREDIBLE-CR enrolled",
        "this category (n = 20, Table 1). Mutually exclusive with",
        "DIS_PNEUMONIA and DIS_CUTI. Multiplies CL by 1.08 (Table 2) and, as",
        "one of the infection indicators, V1 by 1.39."
      ),
      source_name = "PT = 3"
    ),
    DIS_CUTI = list(
      description = paste(
        "Index infection is complicated urinary tract infection (including",
        "acute uncomplicated pyelonephritis in the phase 2 study): 1 =",
        "patient enrolled with cUTI or AUP; 0 = otherwise."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all four infection indicators 0 = subject without infection, the phase 1 cohort)",
      notes = paste(
        "Source coding PT = 1 in Table S2. The model carries a single",
        "pooled coefficient for cUTI and AUP, so DIS_AP is not used. The",
        "cUTI effect on CL is STUDY-SPECIFIC: 1.27 for the phase 2",
        "APEKS-cUTI study (STUDY_CEFIDEROCOL_PHASE2 = 1) and 0.872 for the",
        "phase 3 CREDIBLE-CR study (STUDY_CEFIDEROCOL_PHASE2 = 0). The",
        "paper's own Monte-Carlo simulations used the CREDIBLE-CR value for",
        "cUTI, to represent critically ill patients with carbapenem-resistant",
        "infections. Mutually exclusive with DIS_PNEUMONIA and",
        "DIS_BACTEREMIA. As one of the infection indicators it also",
        "multiplies V1 by 1.39."
      ),
      source_name = "PT = 1"
    ),
    STUDY_CEFIDEROCOL_PHASE2 = list(
      description = paste(
        "1 = patient enrolled in the phase 2 APEKS-cUTI study",
        "(NCT02321800); 0 = otherwise."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase 3 CREDIBLE-CR cUTI coefficient)",
      notes = paste(
        "Only acts in combination with DIS_CUTI = 1, selecting the phase 2",
        "cUTI/AUP coefficient 1.27 on CL instead of the phase 3 CREDIBLE-CR",
        "cUTI coefficient 0.872. Table S2 derives it from the subject ID",
        "(ID < 2000 for phase 1 and phase 2, ID >= 2000 for phase 3); since",
        "every phase 1 subject has PT = 0 it has no effect outside cUTI.",
        "Leave at 0 to reproduce the paper's simulations."
      ),
      source_name = "ID < 2000 (with PT = 1)"
    )
  )

  # Covariates Kawaguchi 2021 screened but did not retain in the final model.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL and V1 (Materials and Methods) and not retained. Range 18-93 years (Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened on CL and V1 and not retained. Sex still enters indirectly through the Cockcroft-Gault CRCL column."
    ),
    AST = list(
      description = "Aspartate aminotransferase concentration",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL and not retained (Table 1 ranges 3-367 U/L)."
    ),
    ALT = list(
      description = "Alanine aminotransferase concentration",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL and not retained (Table 1 ranges 4-153 U/L)."
    ),
    BILI = list(
      description = "Total bilirubin concentration",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL and not retained (Table 1 ranges 0.10-15.20 mg/dL)."
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-White)",
      notes = "Race was screened on CL and V1 and not retained."
    ),
    VENT = list(
      description = "Mechanical ventilation during PK sampling",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not ventilated)",
      notes = paste(
        "Screened on CL and V1 and not retained (Discussion: 'mechanical",
        "ventilation was not a significant covariate on CL or V 1');",
        "post hoc Cmax and daily AUC were similar between ventilated and",
        "non-ventilated pneumonia patients (Fig. 2B)."
      )
    )
  )

  compartmentData <- list(
    central = list(analyte = "cefiderocol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefiderocol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "cefiderocol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 516L,
    n_studies = 5L,
    n_observations = 3427L,
    age_range = "18-93 years (per-study medians 36.0 phase 1, 65.0 APEKS-cUTI, 68.0 APEKS-NP, 67.5 CREDIBLE-CR)",
    weight_range = "25.0-156.0 kg (per-study medians 68.4-76.4 kg)",
    sex_female_pct = 40.1,
    race_ethnicity = c(White = 72.3, Asian = 22.1, Black = 3.3, `Native American or Alaska Native` = 0.2, Other = 2.1),
    disease_state = paste(
      "Pooled: 91 subjects without infection (two phase 1 studies: healthy",
      "volunteers and subjects spanning normal renal function to end-stage",
      "renal disease), 238 patients with cUTI or acute uncomplicated",
      "pyelonephritis (phase 2 APEKS-cUTI), 115 patients with hospital-",
      "acquired, ventilator-associated or health care-associated pneumonia",
      "(phase 3 APEKS-NP), and 72 patients with pneumonia (31), bloodstream",
      "infection/sepsis (20) or cUTI (21) caused by carbapenem-resistant",
      "Gram-negative pathogens (phase 3 CREDIBLE-CR)."
    ),
    renal_function = "Cockcroft-Gault CRCL 4-306 mL/min (per-study medians 69.0-121.0 mL/min); hemodialysis patients in the phase 3 studies were excluded.",
    dose_range = paste(
      "Phase 1: single 0.1-2 g doses or 1-2 g q8h for 10 days over 1 h;",
      "renal-impairment study single 1 g over 1 h. Phase 2: 2 g q8h over",
      "1 h with renal adjustment. Phase 3: 2 g q8h over 3 h with renal",
      "adjustment (2 g q6h for CRCL >= 120 mL/min) (Table S1)."
    ),
    regions = "Japan (phase 1 ascending-dose study), United States (phase 1 renal-impairment study), multinational (phase 2 and phase 3 studies).",
    bioanalytic_methods = "Total plasma cefiderocol by validated LC-MS/MS after 1:1 dilution with 0.2 mol/L ammonium acetate (pH 5); linear 0.1-100 ug/mL, LLOQ 0.1 ug/mL.",
    protein_binding = "In vitro plasma protein binding 57.8%; the paper computes free concentrations with a fixed unbound fraction of 0.422.",
    notes = paste(
      "Baseline characteristics in Table 1; study designs in Table S1; the",
      "final NONMEM control stream (ADVAN11 TRANS4, FOCE-I) in Table S2.",
      "363 below-limit-of-quantification concentrations and 6 anomalous",
      "concentrations were excluded, as were all data from 32 hemodialysis",
      "patients in the phase 3 studies."
    )
  )

  ini({
    # Structural parameters. Table 2, 'Final model' column; THETA(1)-THETA(6)
    # of Table S2. CL, V1 and V2 are the leading constants of the Table 2
    # footnote b equations, i.e. the typical values for a subject without
    # infection at CRCL 83.0 mL/min, 72.6 kg and albumin 3.9 g/dL.
    lcl <- log(4.04); label("Typical total clearance in an uninfected subject at CRCL 83.0 mL/min (L/h)") # Table 2, 'CL (liter/h)' = 4.04 (RSE 1.8%)
    lvc <- log(7.78); label("Typical central volume V1 in an uninfected subject at 72.6 kg and albumin 3.9 g/dL (L)") # Table 2, 'V 1 (liter)' = 7.78 (RSE 5.2%)
    lq <- log(6.19); label("Intercompartmental clearance Q2 between central and peripheral1 (L/h)") # Table 2, 'Q 2 (liter/h)' = 6.19 (RSE 5.7%)
    lvp <- log(5.77); label("Typical peripheral volume V2 at 72.6 kg (L)") # Table 2, 'V 2 (liter)' = 5.77 (RSE 3.2%)
    lq2 <- log(0.127); label("Intercompartmental clearance Q3 between central and peripheral2 (L/h)") # Table 2, 'Q 3 (liter/h)' = 0.127 (RSE 14.1%)
    lvp2 <- log(0.798); label("Peripheral volume V3 (L)") # Table 2, 'V 3 (liter)' = 0.798 (RSE 6.4%)

    # Covariate effects. Table 2 rows and footnote b; Table S2 THETA(8)-THETA(15).
    e_crcl_cl <- 0.682; label("Power exponent on (min(CRCL, 150) / 83.0) for CL (unitless)") # Table 2, 'Effect of CrCL on CL (CrCL cutoff value of 150ml/min)' = 0.682 (RSE 4.0%)
    e_wt_vc_vp <- 0.580; label("Power exponent on (WT / 72.6) shared by V1 and V2 (unitless)") # Table 2, 'Effect of body weight on V 1 and V 2' = 0.580 (RSE 12.2%)
    e_pneumonia_cl <- 0.981; label("Multiplicative factor on CL for pneumonia (unitless)") # Table 2, 'Effect of infection with pneumonia on CL' = 0.981 (RSE 4.1%)
    e_bacteremia_cl <- 1.08; label("Multiplicative factor on CL for bloodstream infection/sepsis (unitless)") # Table 2, 'Effect of infection with BSI/sepsis on CL' = 1.08 (RSE 10.4%)
    e_cuti_cl <- 0.872; label("Multiplicative factor on CL for cUTI in the phase 3 CREDIBLE-CR study (unitless)") # Table 2, 'Effect of infection with cUTI in CREDIBLE-CR study on CL' = 0.872 (RSE 6.4%)
    e_cuti_phase2_cl <- 1.27; label("Multiplicative factor on CL for cUTI/AUP in the phase 2 APEKS-cUTI study (unitless)") # Table 2, 'Effect of infection with cUTI/AUP in APEKS-cUTI study on CL' = 1.27 (RSE 3.1%)
    e_alb_vc <- -0.617; label("Power exponent on (albumin in g/dL / 3.9) for V1 (unitless)") # Table 2, 'Effect of albumin concentration on V 1' = -0.617 (RSE 10.9%); sign from footnote b and the Table S2 THETA(14) bounds (-2.0, -0.5, 0.5)
    e_infect_vc <- 1.39; label("Multiplicative factor on V1 for any active infection (unitless)") # Table 2, 'Effect of infection on V 1' = 1.39 (RSE 6.7%)

    # Interindividual variability. Table 2 reports IIV as CV%; the variances
    # are (CV/100)^2. The table also prints the covariances and their
    # correlations, which over-determines the scale: e.g. for CL-V1,
    # 0.415 * 0.375 * 0.569 = 0.0886, exactly the printed covariance, while
    # the log-normal reading sqrt(log(1 + CV^2)) would give 0.0797. The
    # printed covariances are used as-is. Table S2 $OMEGA BLOCK(3) on
    # ETA(1) CL, ETA(2) V1, ETA(3) V2.
    etalcl + etalvc + etalvp ~ c(
      0.140625,
      0.0886, 0.323761,
      0.0792, 0.150, 0.112896
    ) # Table 2: CL 37.5% (0.375^2), V1 56.9% (0.569^2), V2 33.6% (0.336^2); covariances CL-V1 0.0886 (R 0.415), CL-V2 0.0792 (R 0.629), V1-V2 0.150 (R 0.784)

    # Residual error. Table S2: W = IPRED * THETA(7), Y = IPRED + W * EPS(1)
    # with $SIGMA 1 FIX, so THETA(7) is the proportional SD directly.
    propSd <- 0.205; label("Proportional residual error (fraction)") # Table 2, 'Proportional residual error' = 20.5% (RSE 5.1%)
  })

  model({
    # Covariate terms, Table 2 footnote b and Table S2 $PK. CL rises as a
    # power of CRCL up to 150 mL/min and is constant above it.
    crcl_cl <- min(CRCL, 150)
    alb_gdL <- ALB * 0.1 # SI g/L -> US-convention g/dL, the units of calibration
    infect <- DIS_PNEUMONIA + DIS_BACTEREMIA + DIS_CUTI # mutually exclusive; Table S2 NPT = 1 for PT >= 1
    cuti_phase2 <- DIS_CUTI * STUDY_CEFIDEROCOL_PHASE2
    cuti_phase3 <- DIS_CUTI * (1 - STUDY_CEFIDEROCOL_PHASE2)

    cl <- exp(lcl + etalcl) * (crcl_cl / 83.0)^e_crcl_cl *
      e_pneumonia_cl^DIS_PNEUMONIA * e_bacteremia_cl^DIS_BACTEREMIA *
      e_cuti_cl^cuti_phase3 * e_cuti_phase2_cl^cuti_phase2
    vc <- exp(lvc + etalvc) * (WT / 72.6)^e_wt_vc_vp * (alb_gdL / 3.9)^e_alb_vc *
      e_infect_vc^infect
    q <- exp(lq)
    vp <- exp(lvp + etalvp) * (WT / 72.6)^e_wt_vc_vp
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)

    # Three-compartment intravenous PK (NONMEM ADVAN11 / TRANS4). Doses go
    # into `central` as infusions (3 h for the phase 3 regimens).
    d/dt(central) <- q / vp * peripheral1 + q2 / vp2 * peripheral2 -
      (cl + q + q2) / vc * central
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1
    d/dt(peripheral2) <- q2 / vc * central - q2 / vp2 * peripheral2

    # Total plasma cefiderocol (mg / L = ug/mL). Free concentrations for
    # fT>MIC are 0.422 times this; that scaling is left to post-processing.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
