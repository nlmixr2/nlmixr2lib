Nassar_2022_midazolam <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for intravenous microdosed",
    "midazolam and 1'-hydroxymidazolam in healthy adults, with time-resolved",
    "CYP3A perpetrator effects on midazolam clearance (Nassar 2022). One",
    "compartment for each analyte with linear elimination; a fixed 92% of",
    "midazolam clearance forms 1'-hydroxymidazolam (mass-converted by the",
    "molar-mass ratio 341.77/325.78). Oral voriconazole and intravenous",
    "voriconazole (inhibition), oral efavirenz (activation) and oral",
    "rifampicin (induction) each shift midazolam clearance by a separately",
    "estimated fraction in each of a set of discrete time intervals after",
    "the first perpetrator dose, CL = CLpop * (1 + theta_ij). Doses are",
    "micrograms and concentrations pg/mL."
  )
  reference <- paste(
    "Nassar YM, Hohmann N, Michelet R, Gottwalt K, Meid AD, Burhenne J,",
    "Huisinga W, Haefeli WE, Mikus G, Kloft C. Quantification of the Time",
    "Course of CYP3A Inhibition, Activation, and Induction Using a",
    "Population Pharmacokinetic Model of Microdosed Midazolam Continuous",
    "Infusion. Clin Pharmacokinet. 2022;61(11):1595-1607.",
    "doi:10.1007/s40262-022-01175-6. Final estimates from Table 1; model",
    "structure and perpetrator time windows from the NONMEM control stream",
    "in the Electronic Supplementary Material (ESM 4)."
  )
  vignette <- "Nassar_2022_midazolam"
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "pg/mL"
  )

  covariateData <- list(
    CONMED_VORICONAZOLE_ORAL = list(
      description = "Oral voriconazole perpetrator arm (1 = received a single 400 mg oral voriconazole dose, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no voriconazole; the pooled placebo arm)",
      notes = paste(
        "Time-fixed per subject. Control stream FLAGMED = 4. The effect",
        "on midazolam clearance switches between eight separately",
        "estimated fractions according to the time since the voriconazole",
        "dose, t - T_CONMED: (0, 1], (1, 2], ..., (6, 7] h and (7, 8] h.",
        "The last fraction applies to every time after 7 h, as in the",
        "control stream; sampling ended at 8 h, so later predictions",
        "extrapolate it. Mutually exclusive with the other three",
        "perpetrator indicators."
      ),
      source_name = "FLAGMED (= 4, voriconazole oral)"
    ),
    CONMED_VORICONAZOLE_IV = list(
      description = "Intravenous voriconazole perpetrator arm (1 = received 400 mg voriconazole as a 2-hour intravenous infusion, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no voriconazole; the pooled placebo arm)",
      notes = paste(
        "Time-fixed per subject. Control stream FLAGMED = 5. Seven",
        "time-interval fractions on midazolam clearance: (0, 1], (1, 3],",
        "(3, 4], (4, 5], (5, 6], (6, 7] and (7, 8] h after the start of",
        "the voriconazole infusion. The last applies to every time after",
        "7 h; sampling ended at 8 h."
      ),
      source_name = "FLAGMED (= 5, voriconazole i.v.)"
    ),
    CONMED_EFV = list(
      description = "Oral efavirenz perpetrator arm (1 = received a single 400 mg oral efavirenz dose, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no efavirenz; the pooled placebo arm)",
      notes = paste(
        "Time-fixed per subject. Control stream FLAGMED = 2. Five",
        "time-interval fractions on midazolam clearance: (0, 2], (2, 3],",
        "(3, 4], (4, 5] and (5, 6] h after the efavirenz dose. After 6 h",
        "clearance returns to the placebo value: the authors judged the",
        "(6, 8] h interval uninformative and pooled it with placebo",
        "(Results 3.2; Table 1 CL row label). Efavirenz here is a single-dose",
        "CYP3A activator, not the chronic inducer of most other uses of",
        "this column."
      ),
      source_name = "FLAGMED (= 2, efavirenz)"
    ),
    CONMED_RIFAMPICIN = list(
      description = "Oral rifampicin perpetrator arm (1 = received 600 mg oral rifampicin every 24 h for 2 days, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no rifampicin; the pooled placebo arm)",
      notes = paste(
        "Time-fixed per subject. Control stream FLAGMED = 3. No effect up",
        "to 22 h after the first rifampicin dose (pooled with placebo,",
        "Results 3.2), then five time-interval fractions on midazolam",
        "clearance: (22, 24], (24, 26], (26, 28], (28, 30] and (30, 34] h",
        "after the first dose (the second dose is at 24 h). The last",
        "applies to every time after 30 h; sampling ended at 34 h."
      ),
      source_name = "FLAGMED (= 3, rifampicin)"
    ),
    T_CONMED = list(
      description = "Time of the first perpetrator-drug administration relative to model time zero",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject. Model time zero is the start of the",
        "midazolam bolus plus infusion; the protocol gave the perpetrator",
        "2 h later, so T_CONMED = 2 for every subject in the source trial",
        "(Methods 2.1; the control stream hard-codes the perpetrator",
        "windows on its midazolam-time axis, e.g. TIME.GT.2.AND.TIME.LE.3",
        "for the first voriconazole interval). Inside model() the",
        "time since the perpetrator dose is t - T_CONMED. Irrelevant when",
        "all four perpetrator indicators are 0; supply any finite value."
      ),
      source_name = "derived (perpetrator administered at TIME = 2 h)"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened (Methods 2.3.2) and not retained: 'the demographic characteristics of the healthy population were insignificant as covariates on any PK model parameter' (Results 3.2). 12 of 24 female.",
      source_name = "SEX"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained (Results 3.2). Mean 29.6, range 22-54 years (ESM Table S3).",
      source_name = "AGE"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained (Results 3.2). Mean 71.3, range 55.3-90.5 kg (ESM Table S3).",
      source_name = "WT"
    ),
    HT = list(
      description = "Height",
      units = "m",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained (Results 3.2). Mean 1.73, range 1.59-1.90 m (ESM Table S3).",
      source_name = "HT"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained (Results 3.2). Mean 23.9, range 19.6-29.9 kg/m^2 (ESM Table S3).",
      source_name = "BMI"
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "midazolam",
      units = "ug",
      specimen = "plasma",
      verified = TRUE
    ),
    central_1ohm = list(
      analyte = "1'-hydroxymidazolam",
      units = "ug",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    age_range = "22-54 years (mean 29.6, median 27.0)",
    weight_range = "55.3-90.5 kg (mean 71.3, median 69.3)",
    sex_female_pct = 50,
    disease_state = "healthy volunteers",
    dose_range = paste(
      "Midazolam intravenous bolus 2.70-6.10 ug followed by a continuous",
      "infusion of 2.00-4.40 ug/h for 10 h (36 h in the rifampicin arm),",
      "individualised to a target of about 100 pg/mL (ESM Table S2).",
      "Perpetrators 2 h after the start of midazolam: voriconazole 400 mg",
      "once orally or as a 2-hour intravenous infusion, efavirenz 400 mg",
      "once orally, or rifampicin 600 mg orally every 24 h for 2 days."
    ),
    regions = "Germany (Heidelberg University Hospital)",
    notes = paste(
      "Open-label, fixed-sequence, randomised four-arm phase I trial",
      "(EudraCT 2013-004869-14). Each arm had four perpetrator and two",
      "placebo subjects; the eight placebo subjects were pooled into a",
      "fifth arm. 1858 concentrations, split evenly between midazolam and",
      "1'-hydroxymidazolam, none below the limit of quantification.",
      "Sampling every 15 min for 10 h (hourly for 36 h in the rifampicin",
      "arm). NONMEM 7.4.3, ADVAN6, FOCE-I; uncertainty from sampling",
      "importance resampling."
    )
  )

  ini({
    # Structural parameters (Table 1; control stream THETA(1)-THETA(5))
    lcl <- log(43.9)
    label("Midazolam clearance, placebo / no perpetrator effect (L/h)")      # Table 1 'CL MDZ' 43.9 L/h (RSE 4.37%); THETA(1)
    lvc <- log(56.7)
    label("Midazolam volume of distribution (L)")                            # Table 1 'V MDZ' 56.7 L (RSE 6.84%); THETA(2)
    lcl_1ohm <- log(264)
    label("1'-hydroxymidazolam clearance (L/h)")                             # Table 1 'CL 1'-OH-MDZ' 264 L/h (RSE 7.68%); THETA(3)
    lvc_1ohm <- log(300)
    label("1'-hydroxymidazolam volume of distribution (L)")                  # Table 1 'V 1'-OH-MDZ' 300 L (RSE 8.87%); THETA(4)
    fm <- fixed(0.92)
    label("Fraction of midazolam clearance forming 1'-hydroxymidazolam (fraction)") # Table 1 'Fm' 0.92, footnote b 'Fixed value'; $THETA 0.92 FIX; Methods 2.3.1

    # Oral voriconazole: fractional change in midazolam CL per interval of
    # time since the dose (Table 1; control stream THETA(16)-THETA(23))
    e_conmed_voriconazole_oral_cl_t0to1 <- -0.399
    label("Oral voriconazole effect on midazolam CL, (0, 1] h after dose (fraction)") # Table 1 '[0, 1]' -39.9%
    e_conmed_voriconazole_oral_cl_t1to2 <- -0.290
    label("Oral voriconazole effect on midazolam CL, (1, 2] h after dose (fraction)") # Table 1 '(1, 2]' -29.0%
    e_conmed_voriconazole_oral_cl_t2to3 <- -0.238
    label("Oral voriconazole effect on midazolam CL, (2, 3] h after dose (fraction)") # Table 1 '(2, 3]' -23.8%
    e_conmed_voriconazole_oral_cl_t3to4 <- -0.319
    label("Oral voriconazole effect on midazolam CL, (3, 4] h after dose (fraction)") # Table 1 '(3, 4]' -31.9%
    e_conmed_voriconazole_oral_cl_t4to5 <- -0.639
    label("Oral voriconazole effect on midazolam CL, (4, 5] h after dose (fraction)") # Table 1 '(4, 5]' -63.9%
    e_conmed_voriconazole_oral_cl_t5to6 <- -0.540
    label("Oral voriconazole effect on midazolam CL, (5, 6] h after dose (fraction)") # Table 1 '(5, 6]' -54.0%
    e_conmed_voriconazole_oral_cl_t6to7 <- -0.706
    label("Oral voriconazole effect on midazolam CL, (6, 7] h after dose (fraction)") # Table 1 '(6, 7]' -70.6%
    e_conmed_voriconazole_oral_cl_t7to8 <- -0.694
    label("Oral voriconazole effect on midazolam CL, after 7 h (fraction)")  # Table 1 '(7, 8]' printed '--69.4' (double minus is a typesetting slip); ESM Table S4 relative change -69.4%

    # Intravenous voriconazole (Table 1; control stream THETA(24)-THETA(30))
    e_conmed_voriconazole_iv_cl_t0to1 <- -0.166
    label("IV voriconazole effect on midazolam CL, (0, 1] h after infusion start (fraction)") # Table 1 '[0, 1]' -16.6%
    e_conmed_voriconazole_iv_cl_t1to3 <- -0.111
    label("IV voriconazole effect on midazolam CL, (1, 3] h after infusion start (fraction)") # Table 1 '(1, 3]' -11.1%
    e_conmed_voriconazole_iv_cl_t3to4 <- -0.152
    label("IV voriconazole effect on midazolam CL, (3, 4] h after infusion start (fraction)") # Table 1 '(3, 4]' -15.2%
    e_conmed_voriconazole_iv_cl_t4to5 <- -0.364
    label("IV voriconazole effect on midazolam CL, (4, 5] h after infusion start (fraction)") # Table 1 '(4, 5]' -36.4%
    e_conmed_voriconazole_iv_cl_t5to6 <- -0.393
    label("IV voriconazole effect on midazolam CL, (5, 6] h after infusion start (fraction)") # Table 1 '(5, 6]' -39.3%
    e_conmed_voriconazole_iv_cl_t6to7 <- -0.611
    label("IV voriconazole effect on midazolam CL, (6, 7] h after infusion start (fraction)") # Table 1 '(6, 7]' -61.1%
    e_conmed_voriconazole_iv_cl_t7to8 <- -0.583
    label("IV voriconazole effect on midazolam CL, after 7 h (fraction)")    # Table 1 '(7, 8]' -58.3%

    # Oral efavirenz (Table 1; control stream THETA(6)-THETA(10))
    e_conmed_efv_cl_t0to2 <- 0.150
    label("Efavirenz effect on midazolam CL, (0, 2] h after dose (fraction)") # Table 1 '[0, 2]' 15.0%
    e_conmed_efv_cl_t2to3 <- 0.591
    label("Efavirenz effect on midazolam CL, (2, 3] h after dose (fraction)") # Table 1 '(2, 3]' 59.1%
    e_conmed_efv_cl_t3to4 <- 0.560
    label("Efavirenz effect on midazolam CL, (3, 4] h after dose (fraction)") # Table 1 '(3, 4]' 56.0%
    e_conmed_efv_cl_t4to5 <- 0.285
    label("Efavirenz effect on midazolam CL, (4, 5] h after dose (fraction)") # Table 1 '(4, 5]' 28.5%
    e_conmed_efv_cl_t5to6 <- 0.333
    label("Efavirenz effect on midazolam CL, (5, 6] h after dose (fraction)") # Table 1 '(5, 6]' 33.3%

    # Oral rifampicin (Table 1; control stream THETA(11)-THETA(15))
    e_conmed_rifampicin_cl_t22to24 <- 0.279
    label("Rifampicin effect on midazolam CL, (22, 24] h after first dose (fraction)") # Table 1 '[22, 24]' 27.9%
    e_conmed_rifampicin_cl_t24to26 <- 0.176
    label("Rifampicin effect on midazolam CL, (24, 26] h after first dose (fraction)") # Table 1 '(24, 26]' 17.6%
    e_conmed_rifampicin_cl_t26to28 <- 0.237
    label("Rifampicin effect on midazolam CL, (26, 28] h after first dose (fraction)") # Table 1 '(26, 28]' 23.7%
    e_conmed_rifampicin_cl_t28to30 <- 0.467
    label("Rifampicin effect on midazolam CL, (28, 30] h after first dose (fraction)") # Table 1 '(28, 30]' 46.7%
    e_conmed_rifampicin_cl_t30to34 <- 0.102
    label("Rifampicin effect on midazolam CL, after 30 h (fraction)")        # Table 1 '(30, 34]' 10.2%

    # Between-subject variability. Table 1 reports %CV; variances are
    # omega^2 = log(CV^2 + 1) for the exponential (log-normal) IIV of the
    # control stream (CL = TVCL * EXP(ETA(1)) etc.).
    etalcl ~ 0.04685      # Table 1 'omega2 CL MDZ' 21.9% CV; log(0.219^2 + 1)
    etalvc ~ 0.08563      # Table 1 'omega2 V MDZ' 29.9% CV; log(0.299^2 + 1)
    etalcl_1ohm ~ 0.16246 # Table 1 'omega2 CL 1'-OH-MDZ' 42.0% CV; log(0.420^2 + 1)
    etalvc_1ohm ~ 0.14773 # Table 1 'omega2 V 1'-OH-MDZ' 39.9% CV; log(0.399^2 + 1)

    # Residual error: proportional for each analyte (control stream
    # Y = IPRED + IPRED * EPS(n)); Table 1 %CV taken as the proportional SD.
    propSd <- 0.126
    label("Proportional residual SD, midazolam (fraction)")                  # Table 1 'Residual variability MDZ' 12.6% CV
    propSd_1ohm <- 0.226
    label("Proportional residual SD, 1'-hydroxymidazolam (fraction)")        # Table 1 'Residual variability 1'-OH-MDZ' 22.6% CV
  })

  model({
    # Time since the first perpetrator dose. The control stream codes each
    # window as TIME.GT.a.AND.TIME.LE.b on the midazolam-time axis with the
    # perpetrator at TIME = 2 h, i.e. left-open, right-closed intervals.
    tconmed <- t - T_CONMED

    # Oral voriconazole: eight windows, the last open-ended (TIME.GT.9).
    eff_vori_oral <- e_conmed_voriconazole_oral_cl_t0to1 * (tconmed > 0) * (tconmed <= 1) +
      e_conmed_voriconazole_oral_cl_t1to2 * (tconmed > 1) * (tconmed <= 2) +
      e_conmed_voriconazole_oral_cl_t2to3 * (tconmed > 2) * (tconmed <= 3) +
      e_conmed_voriconazole_oral_cl_t3to4 * (tconmed > 3) * (tconmed <= 4) +
      e_conmed_voriconazole_oral_cl_t4to5 * (tconmed > 4) * (tconmed <= 5) +
      e_conmed_voriconazole_oral_cl_t5to6 * (tconmed > 5) * (tconmed <= 6) +
      e_conmed_voriconazole_oral_cl_t6to7 * (tconmed > 6) * (tconmed <= 7) +
      e_conmed_voriconazole_oral_cl_t7to8 * (tconmed > 7)

    # Intravenous voriconazole: seven windows, the last open-ended.
    eff_vori_iv <- e_conmed_voriconazole_iv_cl_t0to1 * (tconmed > 0) * (tconmed <= 1) +
      e_conmed_voriconazole_iv_cl_t1to3 * (tconmed > 1) * (tconmed <= 3) +
      e_conmed_voriconazole_iv_cl_t3to4 * (tconmed > 3) * (tconmed <= 4) +
      e_conmed_voriconazole_iv_cl_t4to5 * (tconmed > 4) * (tconmed <= 5) +
      e_conmed_voriconazole_iv_cl_t5to6 * (tconmed > 5) * (tconmed <= 6) +
      e_conmed_voriconazole_iv_cl_t6to7 * (tconmed > 6) * (tconmed <= 7) +
      e_conmed_voriconazole_iv_cl_t7to8 * (tconmed > 7)

    # Efavirenz: five windows, back to placebo after 6 h (TIME.GT.8 -> 1).
    eff_efv <- e_conmed_efv_cl_t0to2 * (tconmed > 0) * (tconmed <= 2) +
      e_conmed_efv_cl_t2to3 * (tconmed > 2) * (tconmed <= 3) +
      e_conmed_efv_cl_t3to4 * (tconmed > 3) * (tconmed <= 4) +
      e_conmed_efv_cl_t4to5 * (tconmed > 4) * (tconmed <= 5) +
      e_conmed_efv_cl_t5to6 * (tconmed > 5) * (tconmed <= 6)

    # Rifampicin: placebo up to 22 h (TIME.LE.24 -> 1), then five windows,
    # the last open-ended (TIME.GT.32).
    eff_rif <- e_conmed_rifampicin_cl_t22to24 * (tconmed > 22) * (tconmed <= 24) +
      e_conmed_rifampicin_cl_t24to26 * (tconmed > 24) * (tconmed <= 26) +
      e_conmed_rifampicin_cl_t26to28 * (tconmed > 26) * (tconmed <= 28) +
      e_conmed_rifampicin_cl_t28to30 * (tconmed > 28) * (tconmed <= 30) +
      e_conmed_rifampicin_cl_t30to34 * (tconmed > 30)

    # Eq. 1: CL_ij = CLpop * (1 + theta_ij) * exp(eta_CL)
    eff_conmed_cl <- CONMED_VORICONAZOLE_ORAL * eff_vori_oral +
      CONMED_VORICONAZOLE_IV * eff_vori_iv +
      CONMED_EFV * eff_efv +
      CONMED_RIFAMPICIN * eff_rif

    cl <- exp(lcl + etalcl) * (1 + eff_conmed_cl)
    vc <- exp(lvc + etalvc)
    cl_1ohm <- exp(lcl_1ohm + etalcl_1ohm)
    vc_1ohm <- exp(lvc_1ohm + etalvc_1ohm)

    # Control stream: K10 = (1 - FMET) * CL / V1, K12 = FMET * CL / V1,
    # K20 = CLM / V2. Metabolite amounts are mass of 1'-hydroxymidazolam:
    # formation is scaled by the molar-mass ratio 341.77 / 325.78 = 1.049
    # (Methods 2.3.1; control stream $DES).
    kel <- (1 - fm) * cl / vc
    kmet <- fm * cl / vc
    kel_1ohm <- cl_1ohm / vc_1ohm
    mw_ratio_1ohm <- 341.77 / 325.78

    d/dt(central) <- -kel * central - kmet * central
    d/dt(central_1ohm) <- kmet * central * mw_ratio_1ohm - kel_1ohm * central_1ohm

    # ug / L * 1000 = pg/mL (control stream S1 = V1 / 1000)
    Cc <- 1000 * central / vc
    Cc_1ohm <- 1000 * central_1ohm / vc_1ohm

    Cc ~ prop(propSd)
    Cc_1ohm ~ prop(propSd_1ohm)
  })
}
