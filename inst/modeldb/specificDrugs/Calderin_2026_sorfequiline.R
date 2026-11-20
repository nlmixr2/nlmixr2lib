Calderin_2026_sorfequiline <- function() {
  description <- paste(
    "Semi-mechanistic joint parent-metabolite population PK model for oral",
    "sorfequiline (TBAJ-876, a diarylquinoline antituberculosis drug) and",
    "its active metabolite M3 in healthy adult volunteers (Calderin 2026).",
    "Absorption runs through a chain of three transit compartments (depot ->",
    "transit1 -> transit2 -> transit3, shared rate ktr = 3 / MTT) into an",
    "absorption compartment (transit3) that empties at rate ka: a fraction",
    "1 - FAM reaches the sorfequiline central compartment and the fraction",
    "FAM (gut first-pass metabolism, logit-normal between-subject",
    "variability) reaches the M3 central compartment directly. Both analytes",
    "have three-compartment disposition; sorfequiline is cleared entirely to",
    "M3 (complete conversion, molecular-weight corrected 643.55 / 657.55).",
    "All clearances and volumes scale allometrically with fat-free mass",
    "(exponents 0.75 and 1, reference 53 kg). Food raises ka by 13.6% and",
    "lowers FAM from 0.464 (fasted) to 0.157 (fed); food lengthens the",
    "tablet's mean transit time by 188%. Relative bioavailability is 47.2%",
    "higher in study CL-002 and 46.7% lower for the 400 and 800 mg fasted",
    "suspension doses. Correlated between-subject variability on the two",
    "clearances, between-occasion variability on F, ka and MTT (three",
    "occasions), and combined additive + proportional residual error on",
    "each analyte.",
    sep = " "
  )
  reference <- paste(
    "Calderin JM, Nedelman J, Salinger DH, Vongjarudech T, Denti P,",
    "Svensson EM. Semi-Mechanistic Joint Population Pharmacokinetic Modeling",
    "of the Novel Antituberculosis Drug Sorfequiline and Its Metabolite M3 in",
    "Healthy Volunteers. CPT Pharmacometrics Syst Pharmacol.",
    "2026;15(9):e70331.",
    "doi:10.1002/psp4.70331.",
    sep = " "
  )
  vignette <- "Calderin_2026_sorfequiline"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "sorfequiline", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "sorfequiline", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "sorfequiline", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "sorfequiline", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sorfequiline", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "sorfequiline", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "sorfequiline", units = "mg", specimen = "tissue", verified = TRUE),
    central_m3 = list(analyte = "sorfequiline metabolite M3", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_m3 = list(analyte = "sorfequiline metabolite M3", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2_m3 = list(analyte = "sorfequiline metabolite M3", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass (Janmahasatian formula)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling of every clearance (CL, Q, Q2, CL_m3, Q_m3, Q2_m3; exponent 0.75) and volume (exponent 1), referenced to the cohort median 53 kg (control stream 'TVFFM = 53'; Table 2 footnote b). Computed by the authors from weight, height and sex with the Janmahasatian 2005 formula (Methods, Pharmacokinetic Modeling).",
      source_name = "FFM"
    ),
    FED = list(
      description = "Fed vs fasted dose-record indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = "Per-dose record. 1 = dose given after a high-calorie, high-fat breakfast (Discussion). Raises ka by 13.6% ('IF (FED.EQ.1) FOOD_KA = 1 + THETA(21)'), replaces the fasted FAM typical value 0.464 by the fed value 0.157 ('IF (FED.EQ.1) TVFAM = THETA(24)'), and, for tablets only, lengthens the mean transit time by 188% (see FORM_TABLET).",
      source_name = "FED"
    ),
    FORM_TABLET = list(
      description = "Tablet vs oral suspension formulation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (oral suspension)",
      notes = "Per-dose record. The control-stream column is SOLUTION (1 = suspension, 0 = tablet), so FORM_TABLET = 1 - SOLUTION. The only formulation effect is the food effect on the tablet's mean transit time ('IF (SOLUTION.EQ.0.AND.FED.EQ.1) FOOD_MTT = 1 + THETA(22)'); tablet and suspension otherwise share all parameters, and tablet number (4 x 25 mg vs 1 x 100 mg) had no effect (Results).",
      source_name = "SOLUTION"
    ),
    STUDY_CL002 = list(
      description = "Study CL-002 (NCT05526911) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (study CL-001)",
      notes = "Time-fixed per subject. Relative bioavailability is 47.2% higher in CL-002 than in CL-001 ('IF (PART.EQ.4) BIO_PART4 = 1 + THETA(25)'; PART 4 = CL-002). The authors did not apply this effect in their regimen simulations because its source could not be established; set STUDY_CL002 = 0 to reproduce them.",
      source_name = "PART"
    ),
    DOSE_HIGH = list(
      description = "400 mg or 800 mg dose indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (dose of 200 mg or less)",
      notes = "Per-dose record: 1 when the administered dose is 400 or 800 mg ('IF (DOSES.EQ.400.OR.DOSES.EQ.800) HIGH_DOSE_BIO = 1 + THETA(26)'). Both doses were given only as single fasted oral suspension doses in the SAD part of CL-001, so the 46.7% reduction in F is characterised only there; the authors left it out of their tablet-with-food simulations.",
      source_name = "DOSES"
    ),
    OCC = list(
      description = "Occasion index for between-occasion variability on F, ka and MTT",
      units = "(count)",
      type = "categorical",
      reference_category = "n/a -- selects the per-occasion BOV etas",
      notes = "An occasion is a dosing event and its subsequent observations (Methods). The control stream carries three occasion slots ('IF(OCC==1)', 'IF (OCC==2)', 'IF(OCC==3)'); records with OCC outside 1-3 receive no BOV, exactly as in the control stream. Decomposed in model() into occ1..occ3.",
      source_name = "OCC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 138L,
    n_studies = 2L,
    n_observations = 4638L,
    age_range = "19-55 years (median 35)",
    weight_range = "49-153 kg (median 77)",
    ffm_range = "33-82 kg (median 53)",
    sex_female_pct = 52,
    race_ethnicity = c(
      White = 57,
      Black = 38,
      Asian = 1,
      AmericanIndian = 1,
      Multiple = 3,
      HispanicLatino = 51
    ),
    disease_state = "Healthy adult volunteers",
    dose_range = paste(
      "CL-001 SAD: single 10-800 mg oral suspension fasted plus 100 mg fed;",
      "MAD: 25, 75 or 200 mg once daily for 14 days, fed; RBA: single 100 mg",
      "tablet fed or fasted, or 4 x 25 mg tablets fasted. CL-002: 200 mg once",
      "daily suspension (fed) days 6-13, 165 mg days 14-19, then 200 mg fasted",
      "on days 20-21 and 150 mg fed on days 22-24.",
      sep = " "
    ),
    regions = "United States (San Antonio, Texas; Fair Lawn, New Jersey)",
    notes = paste(
      "Pooled phase 1 studies CL-001 (NCT04493671; SAD n = 55, MAD n = 27,",
      "RBA n = 30) and CL-002 (NCT05526911, n = 26). Demographics from",
      "Calderin 2026 Table 1. Sorfequiline and M3 quantified in the same",
      "plasma samples by HPLC-MS/MS (LLOQ 1 ng/mL each). NONMEM 7.5.1,",
      "FOCE-I; BLQ data handled with the M7+ method.",
      sep = " "
    )
  )

  ini({
    # ---- Sorfequiline disposition (Table 2; control stream $THETA 1-6), typical FFM 53 kg ----
    lcl <- log(4.24); label("Sorfequiline clearance CL/F (L/h)") # Table 2 'Clearance, CLP' 4.24 L/h; $THETA 1
    lvc <- log(17.1); label("Sorfequiline central volume Vc/F (L)") # Table 2 'VC,p' 17.1 L; $THETA 2
    lq <- log(6.64); label("Sorfequiline intercompartmental clearance Q1/F (L/h)") # Table 2 'Q1,p' 6.64 L/h; $THETA 3
    lvp <- log(214); label("Sorfequiline first peripheral volume V2/F (L)") # Table 2 'V2,p' 214 L; $THETA 4
    lq2 <- log(10.4); label("Sorfequiline second intercompartmental clearance Q2/F (L/h)") # Table 2 'Q2,p' 10.4 L/h; $THETA 5
    lvp2 <- log(9180); label("Sorfequiline second peripheral volume V3/F (L)") # Table 2 'V3,p' 9180 L; $THETA 6

    # ---- Absorption (Table 2; $THETA 7-10) ----
    lka <- log(0.127); label("Absorption rate constant from the absorption compartment, fasted (1/h)") # Table 2 'Ka' 0.127 1/h; $THETA 7
    lfdepot <- fixed(log(1)); label("Relative bioavailability F, reference (unitless)") # Table 2 'Bioavailability, F' 1 Fixed; $THETA 8 '(1) FIX'
    lmtt <- log(1.13); label("Mean absorption transit time MTT, ktr = 3 / MTT (h)") # Table 2 'MTT' 1.13 h; $THETA 9

    # ---- M3 disposition (Table 2; $THETA 13-18) ----
    lcl_m3 <- log(38.4); label("M3 clearance CLM (L/h)") # Table 2 'CLM' 38.4 L/h; $THETA 14
    lvc_m3 <- log(3330); label("M3 central volume VC,M (L)") # Table 2 'VC,M' 3330 L; $THETA 13
    lq_m3 <- log(43.5); label("M3 intercompartmental clearance Q1,M (L/h)") # Table 2 'Q1,M' 43.5 L/h; $THETA 16
    lvp_m3 <- log(6460); label("M3 first peripheral volume V2,M (L)") # Table 2 'V2,M' 6460 L; $THETA 15
    lq2_m3 <- log(192); label("M3 second intercompartmental clearance Q2,M (L/h)") # Table 2 'Q2,M' 192 L/h; $THETA 18
    lvp2_m3 <- log(2870); label("M3 second peripheral volume V3,M (L)") # Table 2 'V3,M' 2870 L; $THETA 17

    # ---- Fraction of the absorbed dose entering as M3 (Table 2; $THETA 23-24), logit scale ----
    logitfm <- log(0.464 / (1 - 0.464)); label("Logit of FAM, fraction of absorbed dose entering as M3, fasted (FAM = 0.464)") # Table 2 'FAM' 0.464; $THETA 23
    logitfm_fed <- log(0.157 / (1 - 0.157)); label("Logit of FAM in fed individuals (FAM = 0.157)") # Table 2 'FAM in fed individuals' 0.157; $THETA 24

    # ---- Allometric exponents on fat-free mass (Methods; Table 2 footnote b) ----
    e_ffm_cl <- fixed(0.75); label("Allometric exponent on all clearances (unitless)") # Methods 'fixed to 0.75' for clearances
    e_ffm_vc <- fixed(1); label("Allometric exponent on all volumes (unitless)") # Methods 'fixed to ... 1' for volumes

    # ---- Covariate effects (Table 2 'Covariates'; $THETA 21, 22, 25, 26) ----
    e_fed_ka <- 0.136; label("Fractional change in ka when fed (unitless)") # Table 2 'Food effect on Ka' +13.6%; $THETA 21
    e_fed_tablet_mtt <- 1.88; label("Fractional change in MTT for a tablet taken fed (unitless)") # Table 2 'Food effect on tablet MTT' +188%; $THETA 22
    e_study_cl002_fdepot <- 0.472; label("Fractional change in F in study CL-002 (unitless)") # Table 2 'Change in F in Study CL-002' +47.2%; $THETA 25
    e_dose_high_fdepot <- -0.467; label("Fractional change in F for 400 and 800 mg doses (unitless)") # Table 2 'Change in F at doses >= 400 mg' -46.7%; $THETA 26

    # ---- Between-subject variability ($OMEGA BLOCK(2) and $OMEGA 16) ----
    # Table 2 %CV = sqrt(omega^2): sqrt(0.289) = 53.8%, sqrt(0.147) = 38.3%;
    # correlation 0.147 / sqrt(0.289 * 0.147) = 0.713 (Table 2 'BSVCL correlation' 71.3%).
    etalcl + etalcl_m3 ~ c(
      0.289,
      0.147, 0.147
    ) # $OMEGA BLOCK(2) 0.289 / 0.147 0.147; Table 2 BSV CLP 53.8%, CLM 38.3%, correlation 71.3%
    etalogitfm ~ 0.429 # $OMEGA 16 BSV FAM; Table 2 footnote f 'omega2 = 0.429' on the logit scale (31.8% CV after back-transformation)

    # ---- Between-occasion variability ($OMEGA 20-28; BLOCK(1) SAME for occasions 2-3) ----
    # Three occasion slots, as in the control stream; occasions 2 and 3 repeat
    # the occasion-1 variance (NONMEM SAME) and are therefore fixed().
    # Table 2 %CV = sqrt(omega^2): sqrt(0.282) = 53.1%, sqrt(0.069) = 26.3%, sqrt(0.236) = 48.6%.
    etaiov_lfdepot_1 ~ 0.282 # $OMEGA 20 BOVBIO 0.282; Table 2 BOV F 53.1%
    etaiov_lfdepot_2 ~ fixed(0.282) # $OMEGA 21 BLOCK(1) SAME
    etaiov_lfdepot_3 ~ fixed(0.282) # $OMEGA 22 BLOCK(1) SAME
    etaiov_lka_1 ~ 0.069 # $OMEGA 23 BOVKA 0.069; Table 2 BOV Ka 26.3%
    etaiov_lka_2 ~ fixed(0.069) # $OMEGA 24 BLOCK(1) SAME
    etaiov_lka_3 ~ fixed(0.069) # $OMEGA 25 BLOCK(1) SAME
    etaiov_lmtt_1 ~ 0.236 # $OMEGA 26 BOVMTT 0.236; Table 2 BOV MTT 48.6%
    etaiov_lmtt_2 ~ fixed(0.236) # $OMEGA 27 BLOCK(1) SAME
    etaiov_lmtt_3 ~ fixed(0.236) # $OMEGA 28 BLOCK(1) SAME

    # ---- Residual error (Table 2 'Residual unexplained variability'; $THETA 11, 12, 19, 20) ----
    # Additive SD = $THETA 12 (0 FIX) + 0.2 * LLOQ (0.001 mg/L) = 0.0002 mg/L.
    propSd <- 0.217; label("Sorfequiline proportional residual SD (fraction)") # Table 2 'Sorfequiline proportional error' 21.7%; $THETA 11
    addSd <- fixed(0.0002); label("Sorfequiline additive residual SD (mg/L)") # Table 2 'Sorfequiline additive error' 2 x 10-4 mg/L Fixed
    propSd_m3 <- 0.165; label("M3 proportional residual SD (fraction)") # Table 2 'M3 proportional error' 16.5%; $THETA 19
    addSd_m3 <- fixed(0.0002); label("M3 additive residual SD (mg/L)") # Table 2 'M3 additive error' 2 x 10-4 mg/L Fixed
  })

  model({
    # Molecular weights (Methods): sorfequiline 657.55 g/mol, M3 643.55 g/mol.
    # The control stream tracks moles (F1 = BIO / (657.55 * 1000); CP and CM
    # multiply back by the analyte's MW); tracking mg here and converting
    # parent mass to M3 mass at formation is the same model.
    mw_ratio_m3 <- 643.55 / 657.55

    # Occasion indicators for the BOV etas (control stream IF (OCC==k) blocks).
    occ1 <- (OCC == 1)
    occ2 <- (OCC == 2)
    occ3 <- (OCC == 3)
    iov_fdepot <- occ1 * etaiov_lfdepot_1 + occ2 * etaiov_lfdepot_2 + occ3 * etaiov_lfdepot_3
    iov_ka <- occ1 * etaiov_lka_1 + occ2 * etaiov_lka_2 + occ3 * etaiov_lka_3
    iov_mtt <- occ1 * etaiov_lmtt_1 + occ2 * etaiov_lmtt_2 + occ3 * etaiov_lmtt_3

    # Allometric scaling on fat-free mass, reference 53 kg.
    ffm_factor_cl <- (FFM / 53)^e_ffm_cl
    ffm_factor_v <- (FFM / 53)^e_ffm_vc

    # Sorfequiline disposition.
    cl <- exp(lcl + etalcl) * ffm_factor_cl
    vc <- exp(lvc) * ffm_factor_v
    q <- exp(lq) * ffm_factor_cl
    vp <- exp(lvp) * ffm_factor_v
    q2 <- exp(lq2) * ffm_factor_cl
    vp2 <- exp(lvp2) * ffm_factor_v

    # M3 disposition.
    cl_m3 <- exp(lcl_m3 + etalcl_m3) * ffm_factor_cl
    vc_m3 <- exp(lvc_m3) * ffm_factor_v
    q_m3 <- exp(lq_m3) * ffm_factor_cl
    vp_m3 <- exp(lvp_m3) * ffm_factor_v
    q2_m3 <- exp(lq2_m3) * ffm_factor_cl
    vp2_m3 <- exp(lvp2_m3) * ffm_factor_v

    # Absorption. The food effect on MTT applies to tablets only
    # ('IF (SOLUTION.EQ.0.AND.FED.EQ.1)').
    ka <- exp(lka + iov_ka) * (1 + e_fed_ka * FED)
    mtt <- exp(lmtt + iov_mtt) * (1 + e_fed_tablet_mtt * FED * FORM_TABLET)
    # Three ktr transfers (depot -> transit1 -> transit2 -> transit3); the
    # control stream writes KTR = (NN + 1) / MTT with NN = 2 FIX.
    ktr <- 3 / mtt

    # Fraction of the absorbed dose entering the M3 central compartment
    # (separate fasted and fed typical values on the logit scale).
    fm <- expit(logitfm * (1 - FED) + logitfm_fed * FED + etalogitfm)

    # Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2
    kel_m3 <- cl_m3 / vc_m3
    k12_m3 <- q_m3 / vc_m3
    k21_m3 <- q_m3 / vp_m3
    k13_m3 <- q2_m3 / vc_m3
    k31_m3 <- q2_m3 / vp2_m3

    d / dt(depot) <- -ktr * depot
    d / dt(transit1) <- ktr * depot - ktr * transit1
    d / dt(transit2) <- ktr * transit1 - ktr * transit2
    d / dt(transit3) <- ktr * transit2 - ka * transit3
    d / dt(central) <- ka * (1 - fm) * transit3 - kel * central -
      k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1
    d / dt(peripheral2) <- k13 * central - k31 * peripheral2
    # Sorfequiline is cleared entirely to M3; first-pass M3 enters from transit3.
    d / dt(central_m3) <- mw_ratio_m3 * (ka * fm * transit3 + kel * central) -
      kel_m3 * central_m3 - k12_m3 * central_m3 + k21_m3 * peripheral1_m3 -
      k13_m3 * central_m3 + k31_m3 * peripheral2_m3
    d / dt(peripheral1_m3) <- k12_m3 * central_m3 - k21_m3 * peripheral1_m3
    d / dt(peripheral2_m3) <- k13_m3 * central_m3 - k31_m3 * peripheral2_m3

    f(depot) <- exp(lfdepot + iov_fdepot) * (1 + e_study_cl002_fdepot * STUDY_CL002) *
      (1 + e_dose_high_fdepot * DOSE_HIGH)

    Cc <- central / vc
    Cc_m3 <- central_m3 / vc_m3

    Cc ~ add(addSd) + prop(propSd)
    Cc_m3 ~ add(addSd_m3) + prop(propSd_m3)
  })
}
