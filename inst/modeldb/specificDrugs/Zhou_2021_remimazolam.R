Zhou_2021_remimazolam <- function() {
  description <- paste(
    "Three-compartment population pharmacokinetic model for intravenous",
    "remimazolam pooled across 11 phase I-III studies (359 subjects: healthy",
    "volunteers, procedural-sedation patients and general-anaesthesia",
    "patients) with arterial and venous plasma sampling (Zhou 2021). Clearance",
    "and volumes scale allometrically on body weight (fixed exponents 0.75 and",
    "1, reference 70 kg); clearance is 10% higher in women and 13% lower in",
    "African Americans, and all three volumes are 16% lower in African",
    "Americans. The central volume is fixed at 4.83 L/70 kg from a pilot fit",
    "to two infusion studies, and its inter-individual variability applies",
    "only to subjects from studies with early (< 2 min) sampling. Cc is the",
    "arterial concentration; venous samples are predicted as Cvenous = Cc",
    "times a venous:arterial ratio that rises as an Emax function of time",
    "since the start of an infusion and is a constant 1.28 after the",
    "infusion or bolus ends. One proportional residual error is shared by",
    "both sampling sites.",
    sep = " "
  )
  reference <- paste(
    "Zhou J, Curd L, Lohmer LL, Ossig J, Schippers F, Stoehr T, Schmith V.",
    "Population Pharmacokinetics of Remimazolam in Procedural Sedation With",
    "Nonhomogeneously Mixed Arterial and Venous Concentrations.",
    "Clin Transl Sci. 2021;14(1):326-334. doi:10.1111/cts.12875. PMC7877848.",
    sep = " "
  )
  vignette <- "Zhou_2021_remimazolam"
  units <- list(time = "min", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds. Verified against the Figure S3
  # control stream ($MODEL COMP=(C1) (P1) (P2); ADVAN13 with a three-state
  # $DES), which models remimazolam amounts in mg with S1 = V1/1000 so that
  # concentrations are in ng/mL.
  compartmentData <- list(
    central = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling on every clearance and volume, reference 70 kg (Eq. 1 and Figure S3 'ASCL = (WT/70)**0.75', 'ASV = (WT/70)**1'). The paper states there is no relationship between body weight and remimazolam PK and that the fixed allometry was added to allow potential predictions in paediatric subjects (Methods, Model development). Cohort mean 76.1 kg (SD 17.5) (Table 1).",
      source_name = "WT"
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Figure S3 'SEXEFF = THETA(10)**SEX' with Table 2 'Sex effect on CL, female:male ratio' = 1.1, so the source column SEX is 1 for women. 133 of 359 subjects were female (Table 1).",
      source_name = "SEX"
    ),
    RACE_BLACK = list(
      description = "African American race indicator (1 = African American, 0 = White, Asian or Other).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White, Asian and Other pooled)",
      notes = "Figure S3 'RACE1 = 0; IF (RACE.EQ.1) RACE1 = 1'; Table 2 labels the effects 'Race, African Americans vs. Asians and whites'. Derive RACE_BLACK = as.integer(RACE == 1) from the source column. 82 of 359 subjects were African American (Table 1; Results, Demographics 22.8%).",
      source_name = "RACE"
    ),
    STUDY_NOEARLYPK = list(
      description = "Subject is from a study without early (< 2 min post-dose) sampling (1 = yes; central-volume IIV switched off).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (study with early sampling; full inter-individual variability on V1)",
      notes = "Figure S3 '$PK IF (STDY.GT.11.AND.STDY.NE.17.OR.STDY.EQ.3) THEN V1=TVV1 ... ELSE V1=TVV1*EXP(ETA(4))'. Set to 1 for the procedural-sedation studies CNS7056-002, -004, -006, -008 and -015 and the general-anaesthesia study ONO-2745-03, whose first concentration generally came > 2 minutes after a bolus or rate change; 0 for CNS7056-001, CNS7056-017, ONO-2745-01, ONO-2745-02 and ONO-2745-IVU007 (Methods, Pilot vs. full dataset modeling). The switch is an estimation device -- V1 could not be estimated from sparse early data -- so for a simulation of new subjects set it to 0 to carry the full V1 variability.",
      source_name = "STDY"
    ),
    TINF = list(
      description = "Duration of the current (most recent) intravenous infusion, in hours; 0 for a bolus.",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters only the venous:arterial ratio of the Cvenous output. The source flags a sample as during an infusion when its time since the end of infusion is zero (Figure S3 'IF (TSEOI.EQ.0) TFLG=1'); rxode2 does not expose the infusion duration to model(), so the flag is reconstructed as tad(central) < TINF * 60 (model time is in minutes). For a multi-step infusion give each rate step its own dose record and set TINF to that step's duration. Arterial Cc does not depend on it.",
      source_name = "TSEOI"
    )
  )

  # Screened but not retained in the final model (Methods, Addition of
  # covariates to PopPK Model; Results, Addition of covariates).
  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index.",
      units = "kg/m^2",
      type = "continuous",
      notes = "Significant on Vss in forward addition as a linear effect relative to BMI 25 kg/m^2, then removed in backward elimination (Results). No estimate reported."
    ),
    AGE = list(
      description = "Subject age.",
      units = "years",
      type = "continuous",
      notes = "Prespecified; ONO-2745-03 was included only to span the age range. Not retained. Cohort mean 46.1 years (SD 16.2); 52 subjects aged 65 or older (Table 1)."
    ),
    CRCL = list(
      description = "Creatinine clearance / estimated glomerular filtration rate.",
      units = "mL/min",
      type = "continuous",
      notes = "Treated as correlated with age and ASA class (r^2 > 0.7 rule) and not retained. Cohort mean CrCL 107 mL/min (SD 31.8), eGFR 91.3 mL/min/1.73 m^2 (SD 23.7) (Table 1)."
    ),
    ASA_CLASS = list(
      description = "American Society of Anesthesiologists physical status class (1-4).",
      units = "(category)",
      type = "categorical",
      notes = "Screened; not retained. 234/87/23/15 subjects in class 1/2/3/4 (Table 1)."
    ),
    CONMED_CES1_INHIBITOR = list(
      description = "Chronic co-medication inhibiting carboxylesterase 1 (lovastatin, simvastatin, clopidogrel or telmisartan), pooled (1 = any).",
      units = "(binary)",
      type = "binary",
      notes = "Evaluated as 0 = none vs 1 = any of the four because only 22-24 subjects received one (Methods; Table 1 335/24). Not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 359L,
    n_studies = 11L,
    n_observations = "3642 plasma concentrations (2168 arterial, 1474 venous)",
    age_range = "mean 46.1 years (SD 16.2); study means 25.0-79.2 years",
    weight_range = "mean 76.1 kg (SD 17.5); study means 56.1-91.0 kg",
    sex_female_pct = 100 * 133 / 359,
    race_ethnicity = c(White = 51.3, Black = 22.8, Asian = 25.3, Other = 0.6),
    disease_state = "126 healthy volunteers, 193 procedural-sedation patients (colonoscopy, bronchoscopy; ASA class 1-4) and 40 surgical patients undergoing induction and maintenance of general anaesthesia.",
    dose_range = "Procedural sedation: 5-8 mg over 1 minute with 2, 2.5 or 3 mg top-ups. Healthy volunteers: single IV bolus 0.01-0.5 mg/kg; 1 mg/kg/h infusion for up to 1 h (ONO-2745-02); 5 mg/min for 5 min, 3 mg/min for 15 min, 1 mg/min for 15 min (85 mg, CNS7056-017). General anaesthesia: 4-30 mg/kg/h induction then 1 mg/kg/h maintenance titrated to BIS 40-60.",
    regions = "Japan, United States, European Union",
    notes = "Studies CNS7056-001, -002, -004, -006, -008, -015, -017 and ONO-2745-01, -02, -03, -IVU007 (Table S1). Two studies sampled arterial plasma only, five venous only and three both; simultaneous arterial-venous pairs came from ONO-2745-02 (during and after a 1-h infusion) and from CNS7056-001 and ONO-2745-01 (2-4 h after a bolus). Fitted with NONMEM 7.3 FOCE-I. Demographics from Table 1."
  )

  ini({
    # ---------------------------------------------------------------------
    # Disposition typical values for a 70 kg subject, Zhou 2021 Table 2
    # ('NONMEM estimate (%RSE)' column; identical to Table S2 'Final Model').
    # Q2/V2 map to the canonical peripheral1 pair and Q3/V3 to peripheral2,
    # following the Figure S3 $PK ('peripheral 1' = Q2/V2, 'peripheral 2' =
    # Q3/V3).
    # ---------------------------------------------------------------------
    lcl <- log(1.18)
    label("Clearance, CL (L/min/70 kg)") # Table 2: CL = 1.18 L/minute/70 kg (2%); THETA(1)
    lvc <- fixed(log(4.83))
    label("Central volume, V1 (L/70 kg)") # Table 2: V1 = 4.83 L/70 kg Fixed (pilot-model estimate, Table S2 4.83 (7.9%)); Figure S3 '4.83 FIX'
    lq <- log(0.284)
    label("Intercompartmental clearance to peripheral1, Q2 (L/min/70 kg)") # Table 2: Q2 = 0.284 L/minute/70 kg (2.7%); THETA(3)
    lvp <- log(18.7)
    label("Peripheral1 volume, V2 (L/70 kg)") # Table 2: V2 = 18.7 L/70 kg (2.8%); THETA(4)
    lq2 <- log(1.92)
    label("Intercompartmental clearance to peripheral2, Q3 (L/min/70 kg)") # Table 2: Q3 = 1.92 L/minute/70 kg (6%); THETA(5)
    lvp2 <- log(18)
    label("Peripheral2 volume, V3 (L/70 kg)") # Table 2: V3 = 18 L/70 kg (5.3%); THETA(6)

    # Allometric exponents, fixed (Methods Eq. 1: 'fixed to 0.75 for CL and
    # intercompartmental clearances, Q2 and Q3, and 1.0 for V1 and volumes of
    # the peripheral compartments (V2 and V3)'; Figure S3 ASCL / ASV).
    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent on CL, Q2 and Q3 (unitless)") # Methods Eq. 1; Figure S3 'ASCL = (WT/70)**0.75'
    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent on V1, V2 and V3 (unitless)") # Methods Eq. 1; Figure S3 'ASV = (WT/70)**1'

    # Covariate effects, Table 2. The source enters each as THETA**indicator
    # (Figure S3), i.e. a multiplicative ratio; carried here as its log so
    # that exp(e * indicator) reproduces the ratio exactly.
    e_sexf_cl <- log(1.1)
    label("Log female:male ratio of CL (unitless)") # Table 2: 'Sex effect on CL, female:male ratio' = 1.1 (2.5%); Figure S3 'SEXEFF = THETA(10)**SEX'
    e_race_black_cl <- log(0.87)
    label("Log African-American:other ratio of CL (unitless)") # Table 2: 'Race, African Americans vs. Asians and whites, effect on CL' = 0.87 (2.9%); Figure S3 'RACEEFF1 = THETA(11)**RACE1'
    e_race_black_vc_vp <- log(0.839)
    label("Log African-American:other ratio of V1, V2 and V3 (unitless)") # Table 2: '... effect on Vss' = 0.839 (2.1%); Figure S3 'ASV = (WT/70)**1*COEFFV ;COEFFV is applicable to all volume items'

    # Venous:arterial concentration ratio (Methods Eqs. 3-4; Figure S3
    # $ERROR). All three were estimated in the pilot model and fixed in the
    # full-dataset fits.
    cfven_max <- fixed(1)
    label("Maximum venous:arterial ratio during an infusion, Rmax (unitless)") # Table 2: 'Rmax, max venous:arterial ratio' = 1 Fixed; Figure S3 '1 FIX ; RMAX'
    cfven_t50 <- fixed(1.63)
    label("Time since infusion start to half of Rmax, T50 (min)") # Table 2: 'T50, minutes' = 1.63 Fixed (pilot Table S2 1.63 (15.6%)); Figure S3 '1.63 FIX ;T50 min'
    cfven <- fixed(1.28)
    label("Venous:arterial ratio after the infusion or bolus ends, Ratio2 (unitless)") # Table 2: 'Ratio2, venous:arterial ratio after infusion' = 1.28 Fixed (pilot Table S2 1.28 (1.6%)); Figure S3 '1.28 FIX ;RATIO2'

    # ---------------------------------------------------------------------
    # Inter-individual variability, Table 2. The paper prints each IIV as a
    # bare percentage without stating the transform; converted with the
    # log-normal CV relation omega^2 = log(1 + CV^2). Covariances are
    # r * sqrt(omega_i^2 * omega_j^2) with the printed correlations. The
    # block (CL, Q3, V3) and diagonal (V1, V2) structure follows the Figure S3
    # $OMEGA BLOCK(3) + $OMEGA; there is no IIV on Q2 (Results).
    # ---------------------------------------------------------------------
    etalcl + etalq2 + etalvp2 ~ c(
      0.0511122,
      0.0909497, 0.622210,
      0.0822611, 0.469657, 0.437662
    ) # Table 2: CL IIV 22.9%, Q3 IIV 92.9%, V3 IIV 74.1%; correlations CL/Q3 0.51, CL/V3 0.55, Q3/V3 0.9
    etalvc ~ 0.322583 # Table 2: V1 IIV 61.7%; log(1 + 0.617^2); estimated only in studies with early sampling
    etalvp ~ 0.0596868 # Table 2: V2 IIV 24.8%; log(1 + 0.248^2)

    # Proportional residual error on the site-specific prediction (Eq. 4
    # 'Yobs = Ypred * (1 + eps) * VtoA Ratio'; Figure S3 'Y = RATIO*(F +
    # F*EPS(1))'), one sigma for arterial and venous samples alike.
    propSd <- 0.207
    label("Proportional residual error, arterial and venous (fraction)") # Table 2: 'Residual error' = 20.7 (0.7%)
  })

  model({
    # Amounts in mg, volumes in L; Figure S3 'S1=V1/1000 ; dose in mg and DV
    # in ng/mL'.
    mgL_to_ngmL <- 1000

    # Covariate multipliers (Figure S3 COEFFCL and COEFFV).
    cov_cl <- exp(e_sexf_cl * SEXF + e_race_black_cl * RACE_BLACK)
    cov_v <- exp(e_race_black_vc_vp * RACE_BLACK)
    as_cl <- (WT / 70)^e_wt_cl_q
    as_v <- (WT / 70)^e_wt_vc_vp * cov_v

    # Individual parameters. V1 carries its eta only for subjects from a
    # study with early sampling (Figure S3 IF (STDY...) V1=TVV1 ELSE
    # V1=TVV1*EXP(ETA(4))).
    cl <- exp(lcl + etalcl) * as_cl * cov_cl
    vc <- exp(lvc + etalvc) * as_v
    if (STUDY_NOEARLYPK == 1) vc <- exp(lvc) * as_v
    q <- exp(lq) * as_cl
    vp <- exp(lvp + etalvp) * as_v
    q2 <- exp(lq2 + etalq2) * as_cl
    vp2 <- exp(lvp2 + etalvp2) * as_v

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Arterial plasma concentration (ng/mL); the source sets the ratio to 1
    # for arterial samples.
    Cc <- central / vc * mgL_to_ngmL

    # Venous:arterial ratio (Eq. 3). TSLC, the time since the start of the
    # infusion or bolus, is the time since the last dose into central; the
    # sample is during an infusion while that time is shorter than the
    # infusion duration.
    tslc <- tad(central)
    infusing <- 0
    if (tslc < TINF * 60) infusing <- 1
    ratio_ven <- infusing * cfven_max * tslc / (tslc + cfven_t50) +
      (1 - infusing) * cfven
    Cvenous <- Cc * ratio_ven

    # One sigma shared by both sampling sites (Figure S3 $SIGMA has a single
    # EPS(1)); nlmixr2 requires a distinct endpoint parameter per output, so
    # the venous output reads the same estimate through an alias.
    propSd_Cvenous <- propSd

    Cc ~ prop(propSd)
    Cvenous ~ prop(propSd_Cvenous)
  })
}
