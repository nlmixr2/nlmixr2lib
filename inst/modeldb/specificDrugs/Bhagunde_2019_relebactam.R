Bhagunde_2019_relebactam <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for relebactam, a",
    "small-molecule class A / class C beta-lactamase inhibitor developed as a",
    "fixed-dose combination with imipenem/cilastatin, in adult healthy",
    "volunteers and patients with complicated intra-abdominal infection,",
    "complicated urinary tract infection, or hospital-acquired /",
    "ventilator-associated bacterial pneumonia (Bhagunde 2019; 649",
    "relebactam participants and 4,814 quantifiable plasma concentrations",
    "pooled across 10 phase I-III studies). Zero-order infusion into a",
    "central compartment with first-order linear elimination and linear",
    "distribution to one peripheral compartment. Cockcroft-Gault creatinine",
    "clearance enters clearance and body weight enters the central volume, as",
    "power functions centred on the cohort medians 109 mL/min and 76 kg.",
    "Inter-individual variability is log-normal on CL, V1 and V2 with an",
    "estimated CL-V1 correlation, and residual error is proportional. The",
    "companion imipenem model fitted in the same NONMEM run is",
    "Bhagunde_2019_imipenem."
  )
  reference <- paste(
    "Bhagunde P, Patel P, Lala M, Watson K, Copalu W, Xu M, Kulkarni P,",
    "Young K, Rizk ML (2019).",
    "Population pharmacokinetic analysis for imipenem-relebactam in healthy",
    "volunteers and patients with bacterial infections.",
    "CPT Pharmacometrics Syst Pharmacol 8(10):748-758.",
    "doi:10.1002/psp4.12462. PMCID PMC6813166.",
    sep = " "
  )
  vignette <- "Bhagunde_2019_imipenem_relebactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482. Relebactam (development code MK-7655) was co-administered
  # with imipenem/cilastatin; the assayed analyte here is relebactam in
  # plasma (Results, 'Data analysis': 649 participants provided 4,814
  # quantifiable relebactam plasma concentrations).
  #
  # STATE UNITS. The source NONMEM dataset carried AMT in nmol and DV in
  # nmol/L (Supplementary Material S2 $INPUT comment: "AMT = nmol; DV =
  # nmol/L; TIME = h; CL and Q = L/hr; V1 & V2 = L"). Because the model is
  # linear and CL / Q / V1 / V2 were estimated in L/h and L, the parameters
  # are unchanged by the choice of amount unit: dosing in mg yields Cc in
  # mg/L (= ug/mL, the clinical reporting unit), which is what this file
  # declares. To reproduce the paper's molar exposure metrics divide by the
  # relebactam molar mass 348.37 g/mol and multiply by 1000 to obtain uM.
  compartmentData <- list(
    central = list(analyte = "relebactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "relebactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance calculated with the Cockcroft-Gault equation,",
        "reported as raw mL/min and NOT normalised to 1.73 m^2 body surface",
        "area (Table S1 abbreviation list: 'CrCL, creatinine clearance",
        "(Cockcroft-Gault)'; Table 2 footnote a: 'Calculated using the",
        "Cockroft-Gault equation')."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject. Power effect on CL, centred on the cohort",
        "median: CL = 7.02 * (CrCL / 109)^0.75 * exp(eta). The centring",
        "constant is confirmed twice: Table 3 gives the cohort median CrCL as",
        "109 mL/min, and the final-model control stream (Supplementary",
        "Material S2) codes CL_RLCRCL = ((CRCL/109)**THETA(11)). Relebactam",
        "is predominantly renally excreted and CrCL on relebactam CL was 'the",
        "highest apparent correlation of all covariate-parameter",
        "relationships' in the whole analysis (dMVOF -242.22; Results,",
        "'Covariate model'). The exponent 0.75 is markedly steeper than the",
        "imipenem value of 0.46, which is why relebactam exposure is the more",
        "renally sensitive of the two components. Fitted over an observed",
        "range of 8-406 mL/min (Table 3). Under 15% of participants had",
        "missing CrCL and were imputed with the population median (Methods,",
        "'Data analysis'); the control stream guards the missing code with",
        "IF(CRCL.EQ.-99) CL_RLCRCL = 1. The paper's predicted fold changes in",
        "steady-state AUC0-24 relative to normal renal function are 1.38,",
        "1.89 and 3.05 for mild, moderate and severe renal impairment",
        "(Results, 'Assessing clinical relevance of covariates'; Figure 1)."
      ),
      source_name = "CRCL"
    ),
    WT = list(
      description = "Total body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject. A single estimated power effect on V1 only,",
        "centred on the cohort median 76 kg (Table 3 median weight; control",
        "stream V1_RLWT = ((WT/76)**THETA(14))): exponent 0.70. Weight was",
        "identified as a significant covariate of relebactam V1 (dMVOF",
        "-37.44) but NOT of relebactam CL -- unlike imipenem, which carries",
        "an estimated weight exponent on both. Standard allometric exponents",
        "were deliberately not used anywhere in this analysis; see the",
        "Discussion ('The impact of weight on estimated CrCL was also why",
        "standard allometric exponents were not applied to weight but were",
        "instead estimated directly'). Observed range 39-180 kg (Table 3).",
        "Under 15% of participants had missing weight and were imputed with",
        "the population median; the control stream guards the missing code",
        "with IF(WT.EQ.-99) V1_RLWT = 1. Because weight enters V1 and not CL,",
        "its effect on steady-state AUC0-24 is negligible: the paper reports",
        "predicted fold changes of 1.03 (40-50 kg), 1.02 (50-60 kg) and 1.02",
        "(60-70 kg) versus the 70-90 kg reference, and states that 'The",
        "impact of weight on relebactam exposure was not considered to be of",
        "significance.'"
      ),
      source_name = "WT"
    )
  )

  # Screened as candidate covariates (Table S1: "Covariates investigated for
  # their potential impact on the pharmacokinetic parameters of relebactam
  # and imipenem" -- CrCL, WT, HLTH, Age, Sex, Race on CL and WT, Age, Sex,
  # Race on V1) but not retained in the final relebactam model. No point
  # estimate is published for any of them, so nothing is encoded.
  covariatesDataExcluded <- list(
    DIS_HEALTHY = list(
      description = paste(
        "Health status: 1 = healthy volunteer, 0 = patient with a bacterial",
        "infection. Table S1 abbreviation list: 'HLTH, health status (healthy",
        "or with infection)'."
      ),
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Entered the forward-selection phase on relebactam CL and was then",
        "REMOVED at backward deletion: 'The least significant",
        "covariate-parameter pair was health status on relebactam CL, which",
        "was dropped during the backward deletion phase (dMVOF, -9.11)'",
        "(Results, 'Covariate model'). The backward-deletion criterion was",
        "P < 0.001, i.e. dMVOF 10.83 for 1 degree of freedom, which -9.11",
        "does not meet. No point estimate is published and the final-model",
        "control stream contains no relebactam health-status term (only",
        "V1_IPHLTH, for imipenem). Retained as a covariate effect in the",
        "companion Bhagunde_2019_imipenem model, on V1 rather than CL. Table",
        "3: healthy volunteers 231 (27.0%), patients 624 (73.0%)."
      ),
      source_name = "HLTH"
    ),
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      notes = paste(
        "Screened on both CL and V1 (Table S1, 'Age and sex are of",
        "exploratory interest'). Methods, 'Covariate analysis': 'The effect",
        "of age on the PK of each compound was initially tested with a",
        "linear relationship, with subsequent testing of the power form if",
        "merited.' Not retained. Table 3 gives 18-90 years, median 51. Age is",
        "itself an input to the Cockcroft-Gault CrCL that was retained."
      ),
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Sex, female indicator.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened on both CL and V1 (Table S1) and not retained. The source",
        "dataset column is MALE (Supplementary Material S2 $INPUT), i.e. the",
        "complementary orientation to the canonical SEXF; the transformation",
        "would be SEXF = 1 - MALE. Table 3 gives 519 male (60.7%) and 336",
        "female (39.3%). Sex is also an input to the Cockcroft-Gault CrCL",
        "that was retained."
      ),
      source_name = "MALE"
    ),
    RACE_BLACK = list(
      description = "Race, Black indicator.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened as one level of the RACE covariate (Table S1, 'Race is of",
        "exploratory interest'). Results, 'Covariate model': 'Meaningful",
        "evaluation of race as a covariate was not possible because of an",
        "insufficient number of nonwhite participants (Table 3); a post hoc",
        "analysis showed no obvious trends between unexplained BSV and",
        "race.' Table 3: White 739 (86.4%), Black 32 (3.7%), Asian 44",
        "(5.2%), Other 40 (4.7%)."
      ),
      source_name = "RACE"
    ),
    RACE_ASIAN = list(
      description = "Race, Asian indicator.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened as one level of the RACE covariate and not evaluable for",
        "the same reason as RACE_BLACK; see that entry. Table 3: Asian 44",
        "(5.2%). Two of the ten pooled studies (PN012 and PN019) enrolled",
        "healthy Japanese participants (Table 1)."
      ),
      source_name = "RACE"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 649,
    n_subjects_pooled = 855,
    n_observations = 4814,
    n_studies = 10,
    age_range = "18-90 years",
    age_median = "51 years",
    weight_range = "39-180 kg",
    weight_median = "76 kg",
    sex_female_pct = 39.3,
    race_ethnicity = c(White = 86.4, Black = 3.7, Asian = 5.2, Other = 4.7),
    disease_state = paste(
      "Pooled healthy volunteers (231, 27.0%) and patients (624, 73.0%) with",
      "complicated intra-abdominal infection, complicated urinary tract",
      "infection, or hospital-acquired / ventilator-associated bacterial",
      "pneumonia"
    ),
    renal_function = paste(
      "Creatinine clearance (Cockcroft-Gault) 8-406 mL/min, median 109",
      "mL/min. Distribution across the 852 participants with a recorded",
      "value: <15 mL/min 0.6%, 15 to <30 mL/min 1.1%, 30 to <60 mL/min 8.2%,",
      "60 to <90 mL/min 23.4%, 90 to <150 mL/min 51.5%, 150 to <180 mL/min",
      "12.1%, 180 to <210 mL/min 2.1%, 210 to <250 mL/min 0.5%, >=250",
      "mL/min 0.6%"
    ),
    dose_range = paste(
      "Relebactam 25-1,150 mg as an intravenous infusion, single dose or",
      "every 6 hours; the proposed fixed-dose combination is",
      "imipenem/relebactam 500/250 mg every 6 hours as a 30-minute infusion",
      "in normal renal function"
    ),
    notes = paste(
      "Table 3 ('Summary of the clinical and demographic data for study",
      "participants') and Table 1 ('Summary of studies included in the",
      "analysis'). The 855 pooled participants split into 815 with imipenem",
      "measurements and 649 with relebactam measurements; the demographic",
      "summary in Table 3 is reported for the pooled 855. Three participants",
      "had missing CrCL, so the renal-function percentages are over 852.",
      "Studies PN001, PN002, PN005, PN007, PN009, PN012 and PN019 were phase",
      "I; PN003 and PN004 were phase II; PN013 was phase III."
    )
  )

  ini({
    # =========================================================================
    # Structural disposition parameters, Table 4 ('Final imipenem and
    # relebactam model parameter estimates'), Relebactam 'Final model /
    # Estimate (RSE%)' column. Typical values refer to the covariate
    # reference subject: CrCL 109 mL/min and weight 76 kg. Bootstrap medians
    # over 1,000 replicates agree with the point estimates to within
    # rounding, and all parameters were estimated with RSE < 33%.
    #
    # STRUCTURE. Methods, 'Modeling approach': "two-compartment model of
    # disposition with zero-order i.v. infusion and first-order linear
    # elimination". The control stream (Supplementary Material S2) uses
    # $SUBROUTINE ADVAN3 TRANS4 with S1 = V1_RL, i.e. the CL / V1 / Q / V2
    # parameterisation encoded below.
    # =========================================================================
    lcl <- log(7.02)
    label("Clearance at CrCL 109 mL/min (L/h)") # Table 4, relebactam CL row: 7.02 (RSE 2.0%); 95% CI 6.75-7.29; bootstrap 7.02, 95% CI 6.74-7.31
    lvc <- log(11.08)
    label("Central volume of distribution at weight 76 kg (L)") # Table 4, relebactam V1 row: 11.08 (RSE 2.9%); 95% CI 10.45-11.71; bootstrap 11.06, 95% CI 10.46-11.68
    lvp <- log(6.41)
    label("Peripheral volume of distribution (L)") # Table 4, relebactam V2 row: 6.41 (RSE 3.8%); 95% CI 5.94-6.89; bootstrap 6.43, 95% CI 5.95-6.89
    lq <- log(10.45)
    label("Intercompartmental clearance (L/h)") # Table 4, relebactam Q row: 10.45 (RSE 6.8%); 95% CI 9.04-11.85; bootstrap 10.48, 95% CI 9.09-11.87

    # =========================================================================
    # Covariate effects. Methods, 'Covariate analysis': "The continuous
    # covariates CrCL and body weight were investigated via power
    # relationships centered on the median value." Both centring constants
    # are confirmed against the final-model control stream (Supplementary
    # Material S2): (CRCL/109) and (WT/76). The relebactam model carries no
    # weight effect on CL and no health-status effect on any parameter; the
    # Table 4 relebactam columns print '-' for both rows.
    # =========================================================================
    e_crcl_cl <- 0.75
    label("Power exponent of (CRCL / 109) on CL (unitless)") # Table 4, relebactam 'Covariates on CL / CrCL (power)': 0.75 (RSE 8.3%); 95% CI 0.62-0.87; bootstrap 0.75, 95% CI 0.63-0.86
    e_wt_vc <- 0.70
    label("Power exponent of (WT / 76) on V1 (unitless)") # Table 4, relebactam 'Covariates on V1 WT (power)': 0.70 (RSE 15.8%); 95% CI 0.48-0.92; bootstrap 0.72, 95% CI 0.46-0.96

    # =========================================================================
    # Inter-individual variability. Results, 'Base model': "the data
    # supported similar BSV for relebactam PK", i.e. the same log-normal BSV
    # on CL, V1 and V2 that was selected for imipenem. The control stream
    # confirms the exponential form (CL_RL = TVCL_RL*EXP(ETA(5)), etc.) and
    # confirms that Q carries no BSV ($OMEGA 0 FIX ; [BSV_Q_RL]).
    #
    # SCALE. The Table 4 footnote d is explicit and is NOT the usual
    # log-normal conversion: "Obtained according to the following equation:
    # CV% = sqrt(omega^2) x 100". The tabulated CV% is therefore 100 x the
    # omega standard deviation on the log scale, so omega^2 = (CV%/100)^2
    # directly -- do NOT apply omega^2 = log(CV^2 + 1) to these numbers.
    #
    # The reading is confirmed arithmetically against the control stream's
    # $OMEGA initial estimates, which were seeded from a near-final run:
    #   BSV_V2_RL init 0.168 -> sqrt = 0.410 -> 41.0% vs Table 4 41.1%
    #   BSV_V1_RL init 0.328 -> sqrt = 0.573 -> 57.3% vs Table 4 59.5%
    # and against the $SIGMA initial estimate for the proportional residual
    # (0.0223 -> sqrt = 0.149 -> 14.9% vs Table 4 15.3%). Under the
    # log-normal formula the same inits would imply CV% of 42.9 and 62.5, and
    # a residual CV of 15.0%, which track the table less well; the V2 pair
    # 41.0 / 41.1 is decisive.
    #
    # CORRELATION. Table 4 footnote e: "Correlation between variance
    # parameters calculated as omega^2_ij / sqrt(omega^2_ii * omega^2_jj)",
    # the ordinary correlation coefficient, so the covariance element is
    # corr x sd_CL x sd_V1 = 0.63 x 0.450 x 0.595 = 0.1686825.
    # =========================================================================
    etalcl + etalvc ~ c(
      0.2025000, # CV 45.0% (Table 4, relebactam 'BSV in CL'; RSE 11.6%, shrinkage 16.4%): 0.450^2
      0.1686825, # corr 0.63 (Table 4, relebactam 'Corr CL ~ V1'; RSE 14.2%): 0.63 x 0.450 x 0.595
      0.3540250 # CV 59.5% (Table 4, relebactam 'BSV in V1'; RSE 11.9%, shrinkage 18.8%): 0.595^2
    )
    etalvp ~ 0.1689210 # CV 41.1% (Table 4, relebactam 'BSV in V2'; RSE 30.3%, shrinkage 49.8%): 0.411^2

    # =========================================================================
    # Residual error. Results, 'Base model': "A proportional error model was
    # selected to describe the residual variability." The control stream's
    # $ERROR block writes RV_MK = X2*(F*(1+ERR(3)) + ERR(4)) but fixes the
    # additive term to zero ($SIGMA 0 FIX ; [RES_add_RL]), so the final model
    # is purely proportional.
    # =========================================================================
    propSd <- 0.153
    label("Proportional residual error (fraction)") # Table 4, relebactam 'Proportional Error': 15.3% (RSE 5.6%, shrinkage 14.1%); bootstrap 15.3%
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Individual disposition parameters. Control stream (Supplementary
    #    Material S2), $PK block:
    #      CL_RLCOV = CL_RLCRCL = (CRCL/109)^T11
    #      V1_RLCOV = V1_RLWT   = (WT/76)^T14
    #      TVCL_RL = CL_RLCOV * THETA(5);  CL_RL = TVCL_RL * EXP(ETA(5))
    #      TVV1_RL = V1_RLCOV * THETA(6);  V1_RL = TVV1_RL * EXP(ETA(6))
    #      V2_RL   = THETA(7) * EXP(ETA(7))
    #      Q_RL    = THETA(8) * EXP(ETA(8))   with ETA(8) variance FIXED to 0
    #
    #    Supply CRCL in raw mL/min (Cockcroft-Gault, not BSA-normalised) and
    #    WT in kg.
    # -----------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (CRCL / 109)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 76)^e_wt_vc
    vp <- exp(lvp + etalvp)
    q <- exp(lq)

    # -----------------------------------------------------------------------
    # 2. Micro-rate constants (ADVAN3 TRANS4 uses CL / V1 / Q / V2).
    # -----------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # -----------------------------------------------------------------------
    # 3. Two-compartment intravenous disposition. The zero-order infusion
    #    rate is supplied by the event table (relebactam was given as a
    #    30-minute intravenous infusion in the proposed regimen), so it does
    #    not appear as a model term. Doses in mg with volumes in L give Cc
    #    directly in mg/L (= ug/mL).
    # -----------------------------------------------------------------------
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # -----------------------------------------------------------------------
    # 4. Observation: total plasma relebactam concentration. The control
    #    stream sets S1 = V1_RL, so the predicted quantity is central / vc.
    #    The unbound fraction of relebactam in human plasma is 0.78 (Methods,
    #    'Probability of target attainment simulations'); apply it downstream
    #    when computing free-drug PK/PD metrics such as fAUC0-24/MIC.
    # -----------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
