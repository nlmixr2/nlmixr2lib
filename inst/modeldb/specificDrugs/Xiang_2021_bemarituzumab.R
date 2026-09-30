Xiang_2021_bemarituzumab <- function() {
  description <- "Two-compartment population PK model with parallel linear and Michaelis-Menten elimination for bemarituzumab (anti-FGFR2b antibody) in adults with advanced solid tumours, mainly gastric and gastroesophageal junction adenocarcinoma, given as monotherapy or with mFOLFOX6 (Xiang 2021)"
  reference <- "Xiang H, Liu L, Gao Y, Ahene A, Collins H. Covariate effects and population pharmacokinetic analysis of the anti-FGFR2b antibody bemarituzumab in patients from phase 1 to phase 2 trials. Cancer Chemother Pharmacol. 2021;88(5):899-910. doi:10.1007/s00280-021-04333-y"
  vignette <- "Xiang_2021_bemarituzumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # Bemarituzumab was measured in SERUM by a validated ELISA (Methods,
  # 'Bemarituzumab serum concentration assay in humans').
  compartmentData <- list(
    central = list(analyte = "bemarituzumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "bemarituzumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL (exponent 0.695) and Vc (exponent 0.369), centred on 64 kg",
        "(Equations 7-8). Cohort median 63.9 kg, range 35.5-148 kg (Table 1). Also sets the",
        "mg/kg dose.",
        sep = " "
      ),
      source_name = "WT"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL (exponent -0.657), centred on 38 g/L (Equation 7). The paper reports",
        "albumin in g/L, the canonical unit, so no conversion is applied. Cohort median 38.0 g/L,",
        "range 19.0-50.2 g/L (Table 1).",
        sep = " "
      ),
      source_name = "ALB"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "male (SEXF = 0)",
      notes = paste(
        "Equation 8 enters '(Gender = female)' as an exponential shift on Vc: exp(-0.164) = 0.849,",
        "i.e. 15.1% smaller Vc in females (Results). 60/173 (34.7%) female (Table 1).",
        sep = " "
      ),
      source_name = "Gender"
    ),
    CONMED_CHEMO = list(
      description = "Bemarituzumab given in combination with mFOLFOX6 chemotherapy (1) vs monotherapy (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "monotherapy (CONMED_CHEMO = 0)",
      notes = paste(
        "The paper's 'combotherapy/study' indicator (Equation 7: -0.200 x (Therapy =",
        "combotherapy/study)) on CL: exp(-0.200) = 0.819, i.e. 18.1% lower CL (Results). It is",
        "1 for the 88 patients of study FPA144-004 (FIGHT phase 1, n = 12, and phase 2, n = 76),",
        "all of whom received bemarituzumab with mFOLFOX6 (5-fluorouracil, leucovorin,",
        "oxaliplatin), and 0 for the 85 monotherapy patients of FPA144-001 and FPA144-002. The",
        "authors attribute the effect to the different (first-line) patient population rather",
        "than to a drug interaction, since therapy and study are fully confounded (Discussion).",
        sep = " "
      ),
      source_name = "Therapy"
    )
  )

  # Covariates in Table 1 that were screened on the CL and Vc etas (or tested in
  # the forward search) and not retained (Results, 'Base model development and
  # covariate assessment' and 'Final population PK model').
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age; screened, not retained.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Median 60 (range 23-86) years; Table 1."
    ),
    CRCL = list(
      description = "CKD-EPI estimated glomerular filtration rate; screened, not retained.",
      units = "mL/min/1.73m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Median 95.1 (range 29.1-145) mL/min/1.73 m^2; Table 1."
    ),
    AST = list(
      description = "Aspartate aminotransferase; significant in the eta screen, not retained in the forward search.",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Median 22.0 (range 7.00-113) U/L; Table 1."
    ),
    TBILI = list(
      description = "Total bilirubin; significant in the eta screen, not retained in the forward search.",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Median 8.04 (range 1.71-34.2) umol/L; Table 1."
    ),
    LDH = list(
      description = "Lactate dehydrogenase; significant in the eta screen, not retained in the forward search.",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Median 215 (range 91.0-1041) U/L; Table 1."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 173L,
    n_studies = 3L,
    n_observations = 1552L,
    age_range = "23-86 years",
    age_median = "60 years",
    weight_range = "35.5-148 kg",
    weight_median = "63.9 kg",
    sex_female_pct = 34.7,
    race_ethnicity = c(White = 40.5, Asian = 56.1, Other = 3.47),
    disease_state = paste(
      "Advanced solid tumours: gastric cancer (119/173, 68.8%), gastroesophageal junction",
      "adenocarcinoma (16/173, 9.25%) and other gastrointestinal or solid tumours (38/173, 22.0%).",
      sep = " "
    ),
    dose_range = paste(
      "0.3-15 mg/kg IV Q2W as a 30-min infusion; most patients received 15 mg/kg Q2W, in",
      "FPA144-004 (FIGHT) with an extra 7.5 mg/kg dose on Cycle 1 Day 8.",
      sep = " "
    ),
    regions = "US/Europe/Australia 45.1%, mainland China 8.09%, rest of Asia 46.8% (Table 1).",
    albumin_range = "19.0-50.2 g/L (median 38.0)",
    co_medication = "mFOLFOX6 in all 88 FPA144-004 (FIGHT) patients; the other 85 received bemarituzumab monotherapy.",
    notes = paste(
      "Studies FPA144-001 (phase 1, n = 79), FPA144-002 (phase 1 in Japanese patients, n = 6)",
      "and FPA144-004 (FIGHT; phase 1, n = 12, and phase 2, n = 76); Table 1 and Supplementary",
      "Table 1. No patient developed anti-drug antibodies. Assay LLOQ 0.125 ug/mL; 19 of 1552",
      "samples (1.22%) were below it and omitted.",
      sep = " "
    )
  )

  ini({
    # Structural parameters: typical values for a 64 kg male on monotherapy with
    # albumin 38 g/L (Results, 'For a typical male patient on monotherapy ...').
    lcl <- log(0.311); label("Linear clearance (L/day)") # Table 2, row 'Linear clearance, CL (L/day)' final = 0.311 (3.68%); Equation 7 intercept exp(-4.35) L/h
    lvc <- log(3.58); label("Central volume of distribution (L)") # Table 2, row 'Volume of central compartment, V c (L)' final = 3.58 (1.64%)
    lq <- log(0.952); label("Intercompartmental clearance (L/day)") # Table 2, row 'Inter-compartmental clearance, Q (L/day)' final = 0.952 (12.5%)
    lvp <- log(2.71); label("Peripheral volume of distribution (L)") # Table 2, row 'Volume of peripheral compartment, V p (L)' final = 2.71 (4.92%)
    # Table 2 labels Vmax as 2.80 ug/day. With mg doses, volumes in L and Km in
    # ug/mL (= mg/L), the value is read as 2.80 mg/day; at 2.80 ug/day the
    # Michaelis-Menten pathway would be about 1000-fold smaller than the linear
    # one at every studied dose. See the vignette section 'Units of Vmax'.
    lvmax <- log(2.80); label("Maximum Michaelis-Menten elimination rate (mg/day)") # Table 2, row 'V max (ug/day)' final = 2.80 (4.13%)
    lkm <- log(4.45); label("Michaelis-Menten constant (ug/mL)") # Table 2, row 'K m (ug/mL)' final = 4.45 (8.32%)

    # Covariate effects (Equations 5-8).
    e_wt_cl <- 0.695; label("Power exponent of WT/64 kg on CL (unitless)") # Table 2, row 'Influence of body weight on CL' (theta7) = 0.695 (6.46%)
    e_alb_cl <- -0.657; label("Power exponent of ALB/38 g/L on CL (unitless)") # Table 2, row 'Influence of albumin on CL' (theta9) = -0.657 (15.2%)
    e_conmed_chemo_cl <- -0.200; label("Log-scale shift in CL with mFOLFOX6 combination therapy (unitless)") # Table 2, row 'Influence of combotherapy/study on CL' (theta10) = -0.200 (23.2%)
    e_wt_vc <- 0.369; label("Power exponent of WT/64 kg on Vc (unitless)") # Table 2, row 'Influence of body weight on V c' (theta8) = 0.369 (13.7%)
    e_sexf_vc <- -0.164; label("Log-scale shift in Vc for females (unitless)") # Table 2, row 'Influence of gender on V c' (theta11) = -0.164 (15.9%)

    # IIV. Table 2 gives IIV as CV% = 100 * sqrt(OMEGA): the Results give
    # OMEGA(CL) = 0.0854 and OMEGA(Vc) = 0.0221, whose square roots are the
    # tabulated 29.2% and 14.9%. Vp and Vmax variances are (CV/100)^2.
    etalcl + etalvc ~ c(
      0.0854, # Results text, final-model OMEGA(CL) = 0.0854; Table 2 IIV CL = 29.2%
      0.0128, # Table 2, row 'Covariance between CL and V c' final = 0.0128 (39.8%)
      0.0221 # Results text, final-model OMEGA(Vc) = 0.0221; Table 2 IIV Vc = 14.9%
    )
    etalvp ~ 0.364816 # 0.604^2; Table 2, IIV row 'V p' = 60.4 (10.7%)
    etalvmax ~ 0.948676 # 0.974^2; Table 2, IIV row 'V max' = 97.4 (6.84%)

    # Equation 4: additive error on log-transformed concentrations.
    expSd <- 0.146; label("Additive residual error on the log scale (SD)") # Table 2, row 'Residual variability (%)' = 14.6 (5.19%)
  })
  model({
    # Equations 7-8 (power form of Equation 5; indicator form of Equation 6).
    cl <- exp(lcl + e_wt_cl * log(WT / 64) + e_alb_cl * log(ALB / 38) + e_conmed_chemo_cl * CONMED_CHEMO + etalcl)
    vc <- exp(lvc + e_wt_vc * log(WT / 64) + e_sexf_vc * SEXF + etalvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)
    vmax <- exp(lvmax + etalvmax)
    km <- exp(lkm)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Equations 1-2: parallel linear and Michaelis-Menten elimination from the
    # central compartment; Km is a concentration (ug/mL = mg/L).
    d/dt(central) <- -(vmax / (km + central / vc)) / vc * central - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
