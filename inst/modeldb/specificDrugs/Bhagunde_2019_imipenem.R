Bhagunde_2019_imipenem <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for imipenem in adult",
    "healthy volunteers and patients with complicated intra-abdominal",
    "infection, complicated urinary tract infection, or hospital-acquired /",
    "ventilator-associated bacterial pneumonia (Bhagunde 2019; 815 imipenem",
    "participants and 4,454 quantifiable plasma concentrations pooled across",
    "10 phase I-III studies of the imipenem/relebactam fixed-dose",
    "combination). Zero-order infusion into a central compartment with",
    "first-order linear elimination and linear distribution to one peripheral",
    "compartment. Cockcroft-Gault creatinine clearance and body weight enter",
    "clearance as power functions centred on the cohort medians 109 mL/min",
    "and 76 kg; body weight and healthy-volunteer status enter the central",
    "volume. Inter-individual variability is log-normal on CL, V1 and V2 with",
    "an estimated CL-V1 correlation, and residual error is proportional. The",
    "companion relebactam model fitted in the same NONMEM run is",
    "Bhagunde_2019_relebactam."
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

  # Issue #482. Imipenem is administered as the imipenem/cilastatin
  # fixed-dose combination (PRIMAXIN) with relebactam; cilastatin has no
  # antibacterial activity and was not modelled. The assayed analyte is
  # imipenem in plasma (Table 1; Results, 'Data analysis': "Quantifiable
  # plasma concentrations from 855 participants (815 imipenem, 649
  # relebactam)").
  #
  # STATE UNITS. The source NONMEM dataset carried AMT in nmol and DV in
  # nmol/L (Supplementary Material S2 $INPUT comment: "AMT = nmol; DV =
  # nmol/L; TIME = h; CL and Q = L/hr; V1 & V2 = L"). Because the model is
  # linear and CL / Q / V1 / V2 were estimated in L/h and L, the parameters
  # are unchanged by the choice of amount unit: dosing in mg yields Cc in
  # mg/L (= ug/mL, the clinical reporting unit), which is what this file
  # declares. To reproduce the paper's molar exposure metrics divide by the
  # imipenem molar mass 299.35 g/mol and multiply by 1000 to obtain uM.
  compartmentData <- list(
    central = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = TRUE)
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
        "median: CL = 12.53 * (CrCL / 109)^0.46 * (WT / 76)^0.33 * exp(eta).",
        "The centring constant is confirmed twice: Table 3 gives the cohort",
        "median CrCL as 109 mL/min, and the final-model control stream",
        "(Supplementary Material S2) codes CL_IPCRCL = ((CRCL/109)**THETA(9)).",
        "CrCL was the most influential covariate on imipenem CL (dMVOF",
        "-134.07) and the only covariate judged to warrant dose adjustment.",
        "Fitted over an observed range of 8-406 mL/min (Table 3); the",
        "covariate distribution is heavily weighted to normal and augmented",
        "renal function (51.5% of participants 90 to <150 mL/min, 15.3%",
        ">=150 mL/min, 9.9% <60 mL/min). Under 15% of participants had",
        "missing CrCL and were imputed with the population median (Methods,",
        "'Data analysis'); the control stream guards the missing code with",
        "IF(CRCL.EQ.-99) CL_IPCRCL = 1. The paper's predicted fold changes in",
        "steady-state AUC0-24 relative to normal renal function are 1.22,",
        "1.50 and 2.01 for mild, moderate and severe renal impairment",
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
        "Time-fixed per subject. Two separate estimated power effects,",
        "both centred on the cohort median 76 kg (Table 3 median weight;",
        "control stream CL_IPWT = ((WT/76)**THETA(10)) and V1_IPWT =",
        "((WT/76)**THETA(13))): exponent 0.33 on CL and 0.74 on V1.",
        "Standard allometric exponents were deliberately NOT used. The",
        "Discussion explains why: 'Given that weight is an input variable in",
        "the Cockcroft-Gault calculation of estimated CrCL, this would",
        "suggest that there are additional weight effects not accounted for",
        "by CrCL alone. The impact of weight on estimated CrCL was also why",
        "standard allometric exponents were not applied to weight but were",
        "instead estimated directly.' Both exponents are therefore",
        "conditional on CrCL already being in the model and must not be",
        "reused as free-standing allometric terms. Observed range 39-180 kg",
        "(Table 3). Under 15% of participants had missing weight and were",
        "imputed with the population median; the control stream guards the",
        "missing code with IF(WT.EQ.-99) ... = 1. Predicted fold changes in",
        "steady-state AUC0-24 relative to the 70-90 kg reference are 1.22",
        "(40-50 kg), 1.14 (50-60 kg) and 1.08 (60-70 kg)."
      ),
      source_name = "WT"
    ),
    DIS_HEALTHY = list(
      description = paste(
        "Health status: 1 = healthy volunteer, 0 = patient with a bacterial",
        "infection (complicated intra-abdominal infection, complicated",
        "urinary tract infection, or hospital-acquired / ventilator-associated",
        "bacterial pneumonia). Table S1 abbreviation list: 'HLTH, health",
        "status (healthy or with infection)'."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (patient with bacterial infection)",
      notes = paste(
        "Time-fixed per subject. Linear proportional effect on V1 only:",
        "V1 = 15.83 * (WT / 76)^0.74 * (1 - 0.29 * DIS_HEALTHY) * exp(eta).",
        "Methods, 'Covariate analysis': 'Categorical covariates were modeled",
        "linearly as a proportional change.' The ORIENTATION is fixed by the",
        "final-model control stream (Supplementary Material S2), which codes",
        "IF(HLTH.EQ.0) V1_IPHLTH = 1 ; Most common / IF(HLTH.EQ.1) V1_IPHLTH",
        "= ( 1 + THETA(12)) -- so the tabulated typical V1 of 15.83 L is the",
        "PATIENT value and healthy volunteers carry the 0.71 multiplier.",
        "Patients are indeed the more common group (624 of 855, 73.0%;",
        "Table 3), consistent with the control stream's 'Most common'",
        "comment. Health status was the second-largest covariate effect on",
        "imipenem V1 (dMVOF -56.10) and was retained through backward",
        "deletion. Imipenem V1 is the only parameter this covariate enters;",
        "health status on relebactam CL was dropped at backward deletion",
        "(dMVOF -9.11) and no health-status effect is in the relebactam",
        "model. The Discussion offers the mechanism: 'patients hospitalized",
        "for bacterial infections, including those in the intensive care",
        "unit, often receive i.v. fluids that increase the volume of",
        "distribution.' Health status did not affect steady-state AUC0-24 of",
        "either compound, because it enters V1 and not CL."
      ),
      source_name = "HLTH"
    )
  )

  # Screened as candidate covariates (Table S1: "Covariates investigated for
  # their potential impact on the pharmacokinetic parameters of relebactam
  # and imipenem" -- CrCL, WT, HLTH, Age, Sex, Race on CL and WT, Age, Sex,
  # Race on V1) but not retained in the final model. No point estimate is
  # published for any of them, so nothing is encoded.
  covariatesDataExcluded <- list(
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
    n_subjects = 815,
    n_subjects_pooled = 855,
    n_observations = 4454,
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
      "Imipenem 200-1,000 mg as a 30-minute intravenous infusion, single",
      "dose or every 6 hours; the proposed fixed-dose combination is",
      "imipenem/relebactam 500/250 mg every 6 hours in normal renal function"
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
    # relebactam model parameter estimates'), Imipenem 'Final model /
    # Estimate (RSE%)' column. Typical values refer to the covariate
    # reference subject: a PATIENT (DIS_HEALTHY = 0, see covariateData) with
    # CrCL 109 mL/min and weight 76 kg. Bootstrap medians over 1,000
    # replicates agree with the point estimates to within rounding, and all
    # parameters were estimated with RSE < 33%.
    #
    # STRUCTURE. Methods, 'Modeling approach': "two-compartment model of
    # disposition with zero-order i.v. infusion and first-order linear
    # elimination". The control stream (Supplementary Material S2) uses
    # $SUBROUTINE ADVAN3 TRANS4 with S1 = V1, i.e. the CL / V1 / Q / V2
    # parameterisation encoded below.
    # =========================================================================
    lcl <- log(12.53)
    label("Clearance at CrCL 109 mL/min and weight 76 kg (L/h)") # Table 4, imipenem CL row: 12.53 (RSE 2.0%); 95% CI 12.04-13.02; bootstrap 12.53, 95% CI 12.08-13.04
    lvc <- log(15.83)
    label("Central volume of distribution in a patient at weight 76 kg (L)") # Table 4, imipenem V1 row: 15.83 (RSE 3.2%); 95% CI 14.82-16.83; bootstrap 15.76, 95% CI 14.78-16.85
    lvp <- log(5.84)
    label("Peripheral volume of distribution (L)") # Table 4, imipenem V2 row: 5.84 (RSE 4.0%); 95% CI 5.39-6.29; bootstrap 5.86, 95% CI 5.44-6.32
    lq <- log(11.09)
    label("Intercompartmental clearance (L/h)") # Table 4, imipenem Q row: 11.09 (RSE 6.5%); 95% CI 9.68-12.49; bootstrap 11.15, 95% CI 9.75-12.69

    # =========================================================================
    # Covariate effects. Methods, 'Covariate analysis': "The continuous
    # covariates CrCL and body weight were investigated via power
    # relationships centered on the median value" and "Categorical covariates
    # were modeled linearly as a proportional change." Both centring
    # constants are confirmed against the final-model control stream
    # (Supplementary Material S2): (CRCL/109) and (WT/76).
    # =========================================================================
    e_crcl_cl <- 0.46
    label("Power exponent of (CRCL / 109) on CL (unitless)") # Table 4, imipenem 'Covariates on CL / CrCL (power)': 0.46 (RSE 8.1%); 95% CI 0.39-0.53; bootstrap 0.46, 95% CI 0.38-0.53
    e_wt_cl <- 0.33
    label("Power exponent of (WT / 76) on CL (unitless)") # Table 4, imipenem 'Covariates on CL / WT (power)': 0.33 (RSE 30.5%); 95% CI 0.13-0.53; bootstrap 0.34, 95% CI 0.13-0.55
    e_wt_vc <- 0.74
    label("Power exponent of (WT / 76) on V1 (unitless)") # Table 4, imipenem 'Covariates on V1 WT (power)': 0.74 (RSE 19.1%); 95% CI 0.46-1.01; bootstrap 0.76, 95% CI 0.42-1.07

    # The 'Healthy' row is the linear proportional shift applied to V1 when
    # DIS_HEALTHY = 1, i.e. V1_healthy = V1_patient * (1 + e_healthy_vc).
    # The printed 95% CI reads '-0.34 to 0.23' in the Table 4 body; the upper
    # bound is a lost minus sign. RSE 9.5% on an estimate of -0.29 gives a
    # standard error of 0.0276 and a Wald interval of -0.34 to -0.24, and a
    # CI spanning zero would contradict the parameter's retention through
    # backward deletion at P < 0.001 (dMVOF -56.10).
    e_healthy_vc <- -0.29
    label("Linear proportional effect of healthy-volunteer status on V1 (unitless)") # Table 4, imipenem 'Healthy' row: -0.29 (RSE 9.5%); 95% CI -0.34 to -0.23; bootstrap -0.28, 95% CI -0.34 to -0.23

    # =========================================================================
    # Inter-individual variability. Results, 'Base model': "The most
    # appropriate stochastic model to describe imipenem PK incorporated a
    # log-normal distribution for BSV in CL, V1 and V2." The control stream
    # confirms the exponential form (CL_IP = TVCL_IP*EXP(ETA(1)), etc.) and
    # confirms that Q carries no BSV ($OMEGA 0 FIX ; [BSV_Q_IP]).
    #
    # SCALE. The Table 4 footnote d is explicit and is NOT the usual
    # log-normal conversion: "Obtained according to the following equation:
    # CV% = sqrt(omega^2) x 100". The tabulated CV% is therefore 100 x the
    # omega standard deviation on the log scale, so omega^2 = (CV%/100)^2
    # directly -- do NOT apply omega^2 = log(CV^2 + 1) to these numbers.
    #
    # The reading is confirmed arithmetically against the control stream's
    # $OMEGA initial estimates, which were seeded from a near-final run:
    #   BSV_V2_IP init 0.126  -> sqrt = 0.355 -> 35.5% vs Table 4  35.0%
    #   BSV_V1_IP init 0.599  -> sqrt = 0.774 -> 77.4% vs Table 4  74.4%
    #   BSV_CL_IP init 0.306  -> sqrt = 0.553 -> 55.3% vs Table 4  51.8%
    # and against the $SIGMA initial estimate for the proportional residual
    # (0.025 -> sqrt = 0.158 -> 15.8% vs Table 4 16.1%). Under the log-normal
    # formula the same inits would imply CV% of 36.7, 89.9 and 61.4, which do
    # not track the table.
    #
    # CORRELATION. Table 4 footnote e: "Correlation between variance
    # parameters calculated as omega^2_ij / sqrt(omega^2_ii * omega^2_jj)",
    # the ordinary correlation coefficient, so the covariance element is
    # corr x sd_CL x sd_V1 = 0.77 x 0.518 x 0.744 = 0.2967518.
    # =========================================================================
    etalcl + etalvc ~ c(
      0.2683240, # CV 51.8% (Table 4, imipenem 'BSV in CL'; RSE 9.3%, shrinkage 5.1%): 0.518^2
      0.2967518, # corr 0.77 (Table 4, imipenem 'Corr CL ~ V1'; RSE 12.0%): 0.77 x 0.518 x 0.744
      0.5535360 # CV 74.4% (Table 4, imipenem 'BSV in V1'; RSE 11.6%, shrinkage 7.4%): 0.744^2
    )
    etalvp ~ 0.1225000 # CV 35.0% (Table 4, imipenem 'BSV in V2'; RSE 32.7%, shrinkage 52.9%): 0.350^2

    # =========================================================================
    # Residual error. Results, 'Base model': "A proportional error model was
    # selected to describe the residual variability." The control stream's
    # $ERROR block writes RV_IP = X1*(F*(1+ERR(1)) + ERR(2)) but fixes the
    # additive term to zero ($SIGMA 0 FIX ; [RES_add_IP]), so the final model
    # is purely proportional.
    # =========================================================================
    propSd <- 0.161
    label("Proportional residual error (fraction)") # Table 4, imipenem 'Proportional Error': 16.1% (RSE 6.7%, shrinkage 19.0%); bootstrap 16.1%
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Individual disposition parameters. Control stream (Supplementary
    #    Material S2), $PK block:
    #      CL_IPCOV = CL_IPCRCL * CL_IPWT = (CRCL/109)^T9 * (WT/76)^T10
    #      V1_IPCOV = V1_IPHLTH * V1_IPWT = (1 + T12*HLTH) * (WT/76)^T13
    #      TVCL_IP = CL_IPCOV * THETA(1);  CL_IP = TVCL_IP * EXP(ETA(1))
    #      TVV1_IP = V1_IPCOV * THETA(2);  V1_IP = TVV1_IP * EXP(ETA(2))
    #      V2_IP   = THETA(3) * EXP(ETA(3))
    #      Q_IP    = THETA(4) * EXP(ETA(4))   with ETA(4) variance FIXED to 0
    #
    #    Supply CRCL in raw mL/min (Cockcroft-Gault, not BSA-normalised) and
    #    WT in kg. DIS_HEALTHY is 1 for a healthy volunteer and 0 for a
    #    patient, so the typical values above are patient values.
    # -----------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (CRCL / 109)^e_crcl_cl * (WT / 76)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 76)^e_wt_vc * (1 + e_healthy_vc * DIS_HEALTHY)
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
    #    rate is supplied by the event table (imipenem was given as a
    #    30-minute intravenous infusion), so it does not appear as a model
    #    term. Doses in mg with volumes in L give Cc directly in mg/L
    #    (= ug/mL).
    # -----------------------------------------------------------------------
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # -----------------------------------------------------------------------
    # 4. Observation: total plasma imipenem concentration. The control stream
    #    sets S1 = V1_IP, so the predicted quantity is central / vc. The
    #    unbound fraction of imipenem in human plasma is 0.80 (Methods,
    #    'Probability of target attainment simulations'); apply it downstream
    #    when computing free-drug PK/PD metrics such as %fT>MIC.
    # -----------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
