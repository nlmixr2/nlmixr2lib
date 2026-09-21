Cirincione_2018_apixaban <- function() {
  description <- "Population PK model for apixaban in subjects with nonvalvular atrial fibrillation (Cirincione 2018, stage 2 final model): two-compartment with first-order absorption and elimination, clearance split into renal and nonrenal arms, and dose-dependent relative bioavailability."
  reference <- "Cirincione B, Kowalski K, Nielsen J, Roy A, Thanneer N, Byon W, Boyd R, Wang X, Leil T, LaCreta F, Ueno T, Oishi M, Frost C. Population Pharmacokinetics of Apixaban in Subjects With Nonvalvular Atrial Fibrillation. CPT: Pharmacometrics & Systems Pharmacology. 2018;7(11):728-738. doi:10.1002/psp4.12347"
  vignette <- "Cirincione_2018_apixaban"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot = list(analyte = "apixaban", units = "mg", specimen = "administration site", verified = FALSE),
    central = list(analyte = "apixaban", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "apixaban", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the nonrenal clearance arm; reference age 65 years (Table S2 Eq. [1c]), which is also the age of the paper's typical reference NVAF subject.",
      source_name = "AGE"
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the apparent central volume; reference weight 70 kg (Table S2 Eq. [1e]), the weight of the paper's typical reference NVAF subject. Time-fixed at baseline; mapped to the general-scope WT canonical rather than WT_BASE because the source analysis carries no paired time-varying weight column.",
      source_name = "D_WTB"
    ),
    CRCL = list(
      description = "Baseline creatinine clearance calculated with the Cockcroft-Gault equation (raw, NOT body-surface-area normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Proportional (exponent fixed to 1) effect on the apparent renal clearance arm; reference 80 mL/min. The source NONMEM control stream (Supplemental File S2) caps the covariate at 150 mL/min before forming the ratio ('FLG=0; IF(D_CCRCLB.GE.150) FLG=1; CRCL=(1-FLG)*D_CCRCLB+FLG*150'), which the paper introduces in Eq. 1 to avoid the Cockcroft-Gault formula's limitations at extremes of body weight. Raw mL/min rather than the register's mL/min/1.73 m^2 default, following the Delattre 2010 amikacin precedent for raw Cockcroft-Gault clearance.",
      source_name = "D_CCRCLB"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Fractional-change effect on the nonrenal clearance arm; female subjects have 21.6 percent lower CL_NR/F than male. Source column SEX is coded 1 = male, 2 = female and the control stream derives 'FEM=0; IF (SEX.EQ.2) FEM=1', matching the canonical orientation with no sign inversion.",
      source_name = "SEX"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = "Fractional-change effect on total apparent clearance. Asians were 14.9 percent of the analysis dataset (Table 2). The companion model Cirincione_2018_apixaban_asian_subgroups.R replaces this single indicator with the paper's ad hoc Japanese / Korean / other-Asian stratification.",
      source_name = "RACE"
    ),
    DIS_NVAF = list(
      description = "Nonvalvular atrial fibrillation subject-status indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy subject or, when DIS_ACS is 1, acute coronary syndrome patient)",
      notes = "Fractional-change effects on BOTH total apparent clearance and the apparent central volume. DIS_NVAF and DIS_ACS are mutually exclusive; both zero selects the healthy-volunteer reference. Derived in the control stream from the subject-type column as 'AF=0; IF (STYP.EQ.4) AF=1'.",
      source_name = "STYP"
    ),
    DIS_ACS = list(
      description = "Acute coronary syndrome subject-status indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy subject or, when DIS_NVAF is 1, nonvalvular atrial fibrillation patient)",
      notes = "Fractional-change effects on BOTH total apparent clearance and the apparent central volume. Mutually exclusive with DIS_NVAF; both zero selects the healthy-volunteer reference. Derived in the control stream as 'ACS=0; IF (STYP.EQ.5) ACS=1'. Contributed by the APPRAISE-1 and Japan ACS phase II studies.",
      source_name = "STYP"
    ),
    CONMED_CYP3A4_PGP_INH = list(
      description = "Concomitant strong or moderate CYP3A4 / P-glycoprotein inhibitor at the time of the pharmacokinetic sample",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant strong or moderate CYP3A4 / P-gp inhibitor)",
      notes = "Time-varying at the observation record. The strong and moderate categories were deliberately collapsed by the authors because only three subjects contributed a sample while on a strong inhibitor (Table 2: 3 of 4,385 strong, 718 of 4,385 moderate). Derived in the control stream as 'INH=0; IF(INHS.EQ.1.OR.INHM.EQ.1)INH=1'. The paper does not enumerate the individual agents. A companion strong-CYP3A4/P-gp-inducer indicator was tested and dropped (Table 3: theta20 not evaluable; control stream '(0 FIX) ;20 INDS on TCL'), so it is not carried here.",
      source_name = "INHS / INHM"
    ),
    DOSETIME_EVENING = list(
      description = "Indicator that the dose was administered in the evening clock-time window",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (morning dose, 12:01 am to 11:00 am, or afternoon dose, 11:01 am to 5:00 pm)",
      notes = "Time-varying at the dose record. Fractional-change effect on the absorption rate constant: evening dosing reduced ka by 43 percent. The source column AMPM has three levels (1 = morning 12:01 am to 11:00 am, 2 = evening 5:01 pm to 12:00 am, 3 = afternoon 11:01 am to 5:00 pm); only level 2 carries an estimated effect. The afternoon effect (theta27) was not evaluable in the final models (Table 3) and is fixed to zero in the control stream, so morning and afternoon are pooled into the reference category.",
      source_name = "AMPM"
    ),
    DOSE_APIXABAN_MG = list(
      description = "Administered apixaban dose amount",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the dose-dependent reduction in relative bioavailability (Eq. 4 / Table S2 Eq. [1n]) attributed to decreased dissolution at higher doses. Anchored at 2.5 mg, where relative bioavailability is fixed to 1; the 47.5 mg denominator makes I50 the fractional reduction at a 50 mg dose. Dose range in the analysis dataset was 2.5 to 50 mg.",
      source_name = "DOSE"
    ),
    STUDY_APPRAISE1 = list(
      description = "Indicator that the observation record comes from the global phase II acute coronary syndrome study APPRAISE-1",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase I healthy-volunteer studies, Japan NVAF phase II, or Japan ACS phase II)",
      notes = "Selects the APPRAISE-1 log-scale residual error magnitude. Source NONMEM control stream: 'W=THETA(21); IF (STUDYID.EQ.23) W=THETA(22)'. Paired with STUDY_ARISTOTLE; both zero selects the pooled healthy-volunteer plus Japan phase II residual magnitude.",
      source_name = "STUDYID"
    ),
    STUDY_ARISTOTLE = list(
      description = "Indicator that the observation record comes from the global phase III nonvalvular atrial fibrillation study ARISTOTLE",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase I healthy-volunteer studies, Japan NVAF phase II, or Japan ACS phase II)",
      notes = "Selects the ARISTOTLE log-scale residual error magnitude. Source NONMEM control stream: 'IF (STUDYID.EQ.30) W=THETA(23)'. Paired with STUDY_APPRAISE1; both zero selects the pooled healthy-volunteer plus Japan phase II residual magnitude.",
      source_name = "STUDYID"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 4385,
    n_studies = 12,
    age_range = "18-94 years",
    age_mean = "65.44 years (SD 13.0); median 68",
    weight_range = "32-198.2 kg",
    weight_mean = "83.44 kg (SD 19.6); median 81.4",
    crcl_range = "11.9-319.7 mL/min",
    crcl_mean = "85.2 mL/min (SD 34.5); median 79.3",
    sex_female_pct = 29.76,
    race_ethnicity = c(White = 83.22, Asian = 14.87, `Black/African American` = 1.39, Other = 0.52),
    disease_state = "Pooled healthy volunteers (270), subjects with nonvalvular atrial fibrillation (3,071) and subjects with acute coronary syndrome (1,044)",
    dose_range = "2.5-50 mg, single dose and once- or twice-daily multiple dose",
    regions = "Global (phase I in the United States; phase II and III including Japan, Korea and China)",
    n_observations = 11968,
    notes = "Stage 2 final model of a two-stage analysis. Stage 1 used 9,036 concentrations from 270 phase I subjects and 1,183 phase II subjects; stage 2 added 2,932 concentrations from 2,932 ARISTOTLE phase III NVAF subjects. Typical reference subject used for the paper's covariate displays: non-Asian male with NVAF, aged 65 years, 70 kg, cCrCL 80 mL/min, no concomitant CYP3A4/P-gp inhibitor."
  )

  ini({
    # Structural parameters: Table 3, 'stage 2 final model' column.
    lka <- log(0.473)
    label("Absorption rate constant (ka, 1/h)")
    lcl_renal <- log(1.57)
    label("Apparent renal clearance at the reference cCrCL of 80 mL/min (CLR/F, L/h)")
    lcl_nonren <- log(2.02)
    label("Apparent nonrenal clearance for the reference 65-year-old male (CLNR/F, L/h)")
    lvc <- log(30.0)
    label("Apparent central volume of distribution for the reference 70-kg healthy subject (Vc/F, L)")
    lq <- log(1.91)
    label("Apparent intercompartmental clearance (Q/F, L/h)")
    lvp <- log(27.0)
    label("Apparent peripheral volume of distribution (Vp/F, L)")

    # Covariate effects: Table 3, 'stage 2 final model' column; functional
    # forms from Table S2 Eqs. [1a]-[1e] and Supplemental File S2.
    e_dosetime_evening_ka <- -0.433
    label("Fractional change in ka for an evening dose vs a morning or afternoon dose (theta10, unitless)")
    e_crcl_cl_renal <- fixed(1)
    label("Power exponent on (cCrCL capped at 150 / 80) for CLR/F (theta7, unitless)")
    e_age_cl_nonren <- -0.429
    label("Power exponent on (AGE / 65) for CLNR/F (theta14, unitless)")
    e_sexf_cl_nonren <- -0.216
    label("Fractional change in CLNR/F for female vs male (theta15, unitless)")
    e_race_asian_cl <- -0.119
    label("Fractional change in total CL/F for Asian vs non-Asian (theta16, unitless)")
    e_dis_nvaf_cl <- -0.139
    label("Fractional change in total CL/F for NVAF subjects vs healthy subjects (theta17, unitless)")
    e_dis_acs_cl <- -0.215
    label("Fractional change in total CL/F for ACS subjects vs healthy subjects (theta18, unitless)")
    e_conmed_cyp3a4_pgp_inh_cl <- -0.146
    label("Fractional change in total CL/F with a concomitant strong or moderate CYP3A4/P-gp inhibitor (theta19, unitless)")
    e_wt_vc <- 0.790
    label("Power exponent on (WT / 70) for Vc/F (theta11, unitless)")
    e_dis_nvaf_vc <- -0.0405
    label("Fractional change in Vc/F for NVAF subjects vs healthy subjects (theta12, unitless)")
    e_dis_acs_vc <- -0.180
    label("Fractional change in Vc/F for ACS subjects vs healthy subjects (theta13, unitless)")

    # Dose-dependent relative bioavailability: Eq. 4 and Table S2 Eqs.
    # [1l]-[1n]. The stage 1 Imax/ED50 form of Eq. 2 was replaced by this
    # power form in the updated stage 1 and stage 2 models because gamma,
    # logit-Imax and ED50 were correlated above rho = 0.95.
    lgamma <- log(0.857)
    label("Shape parameter of the dose-relative-bioavailability relationship (gamma, theta8, unitless)")
    logiti50 <- -0.322
    label("Logit of the fractional reduction in relative bioavailability at a 50-mg dose (I50, theta9, unitless)")

    # IIV: Table 3, 'stage 2 final model' column. Exponential on the
    # microconstants of the ADVAN4/TRANS1 parameterisation, diagonal OMEGA
    # ('$OMEGA DIAGONAL(5)' in Supplemental File S2). The authors moved the
    # random effects onto the microconstants because random effects on CL/F
    # and Vc/F were estimated with a correlation of ~1.
    etalka ~ 0.263 # Table 3, stage 2 final, 'omega2-ka' = 0.263 (51.3 %CV)
    etalkel ~ 0.0954 # Table 3, stage 2 final, 'omega2-k' = 0.0954 (30.9 %CV)
    etalvc ~ 0.0294 # Table 3, stage 2 final, 'omega2-Vc/F' = 0.0294 (17.1 %CV)
    etalk21 ~ 0.240 # Table 3, stage 2 final, 'omega2-k21' = 0.240 (49.0 %CV)
    etalk12 ~ 1.55 # Table 3, stage 2 final, 'omega2-k12' = 1.55 (124 %CV)

    # Residual error: Table 3, 'stage 2 final model' column. The source fits
    # log-transformed observations with an additive error on the log scale
    # and SIGMA fixed to 1, so these thetas ARE the log-scale residual SDs.
    expSdHvJapan <- 0.31
    label("Log-scale residual SD for phase I healthy volunteers and the Japan NVAF / Japan ACS phase II studies (theta21)")
    expSdAppraise1 <- 0.668
    label("Log-scale residual SD for the APPRAISE-1 phase II study (theta22)")
    expSdAristotle <- 0.460
    label("Log-scale residual SD for the ARISTOTLE phase III study (theta23)")
  })
  model({
    # Absorption rate, Table S2 Eq. [1a]: only the evening window shifts ka;
    # the afternoon effect was not evaluable and is pooled with morning.
    ka <- exp(lka + etalka) * (1 + e_dosetime_evening_ka * DOSETIME_EVENING)

    # Renal clearance arm, Eq. 1 / Table S2 Eq. [1b]. cCrCL is capped at
    # 150 mL/min before forming the ratio, reproducing the control stream's
    # FLG branch rather than extrapolating Cockcroft-Gault past its range.
    crclFlag <- (CRCL >= 150)
    crclEff <- (1 - crclFlag) * CRCL + crclFlag * 150
    cl_renal <- exp(lcl_renal) * (crclEff / 80)^e_crcl_cl_renal

    # Nonrenal clearance arm, Table S2 Eq. [1c].
    cl_nonren <- exp(lcl_nonren) * (AGE / 65)^e_age_cl_nonren *
      (1 + e_sexf_cl_nonren * SEXF)

    # Total apparent clearance, Table S2 Eq. [1d]: the two arms add, then
    # the race, subject-status and comedication effects multiply the sum.
    cl <- (cl_renal + cl_nonren) *
      (1 + e_race_asian_cl * RACE_ASIAN) *
      (1 + e_dis_nvaf_cl * DIS_NVAF) *
      (1 + e_dis_acs_cl * DIS_ACS) *
      (1 + e_conmed_cyp3a4_pgp_inh_cl * CONMED_CYP3A4_PGP_INH)

    # Typical apparent volumes, Table S2 Eqs. [1e], [1h], [1i].
    vcTypical <- exp(lvc) * (WT / 70)^e_wt_vc *
      (1 + e_dis_nvaf_vc * DIS_NVAF) *
      (1 + e_dis_acs_vc * DIS_ACS)
    q <- exp(lq)
    vp <- exp(lvp)

    # Microconstants, Table S2 Eqs. [1f], [1g], [1j], [1k]. Note that kel,
    # k12 and k21 are formed from the TYPICAL volumes -- the eta on vc
    # scales the concentration only, it does not feed back into the rate
    # constants. This reproduces the control stream verbatim:
    #   K = (TVCL/TVV2)*EXP(ETA(2)); V2 = TVV2*EXP(ETA(3))
    #   K32 = (TVQ/TVV3)*EXP(ETA(4)); K23 = (TVQ/TVV2)*EXP(ETA(5))
    vc <- vcTypical * exp(etalvc)
    kel <- (cl / vcTypical) * exp(etalkel)
    k12 <- (q / vcTypical) * exp(etalk12)
    k21 <- (q / vp) * exp(etalk21)

    # Dose-dependent relative bioavailability, Eq. 4 / Table S2 Eqs.
    # [1l]-[1n]. Frel is anchored at 1 for the 2.5-mg dose; the
    # (DOSE - 2.5) excess is clamped at zero so doses at or below the
    # anchor cannot produce a fractional power of a negative number.
    gamma <- exp(lgamma)
    i50 <- expit(logiti50)
    doseExcess <- (DOSE_APIXABAN_MG > 2.5) * (DOSE_APIXABAN_MG - 2.5)
    frel <- 1 - i50 * (doseExcess / 47.5)^gamma

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- frel

    # Amounts are in mg and vc is in L, so central/vc is mg/L; multiply by
    # 1000 to report ng/mL. This is the control stream's S2 = V2/1000.
    Cc <- 1000 * central / vc

    # Study-specific log-scale residual SD, Supplemental File S2 $ERROR:
    #   W = THETA(21); IF (STUDYID.EQ.23) W = THETA(22)
    #                  IF (STUDYID.EQ.30) W = THETA(23)
    expSd <- expSdAppraise1 * STUDY_APPRAISE1 +
      expSdAristotle * STUDY_ARISTOTLE +
      expSdHvJapan * (1 - STUDY_APPRAISE1 - STUDY_ARISTOTLE)
    Cc ~ lnorm(expSd)
  })
}
