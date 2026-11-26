Mohamed_2020_upadacitinib <- function() {
  description <- paste(
    "Two-compartment population PK model for oral upadacitinib in healthy",
    "adults and adults with rheumatoid arthritis (RA) or Crohn's disease (CD),",
    "with first-order absorption after an absorption lag time (Mohamed 2020,",
    "CELEST exposure-response analysis). Covariates: RA on CL/F, RA or CD on",
    "the absorption rate, sex and creatinine clearance on CL/F, and sex and",
    "body weight on Vc/F. Interindividual variability on CL/F differs between",
    "healthy, RA and CD subjects, and on Vc/F between healthy subjects and",
    "patients. The lag time carries interindividual variability in phase 1",
    "studies and Box-Cox-transformed inter-occasion variability at phase 2",
    "study visits. The immediate-release (IR) formulation is the reference;",
    "the extended-release (ER) branch carries the ER absorption and relative",
    "bioavailability estimates that the paper fitted to a phase 1 ER study and",
    "used to simulate ER regimens in Crohn's disease (Table 2)."
  )
  reference <- paste(
    "Mohamed M-EF, Klunder B, Lacerda AP, Othman AA. Exposure-Response",
    "Analyses for Upadacitinib Efficacy and Safety in the Crohn's Disease",
    "CELEST Study and Bridging to the Extended-Release Formulation. Clin",
    "Pharmacol Ther. 2020;107(3):639-649. doi:10.1002/cpt.1668"
  )
  vignette <- "Mohamed_2020_upadacitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "upadacitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "upadacitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "upadacitinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline total body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on Vc/F, exponent 0.490, reference 74 kg (Table S3 footnote c; control stream '(WTKG/74)**THETA(11)').",
      source_name = "WTKG"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male). Table S3 footnotes b and c set the sex factors to 1 for males.",
      notes = "Female CL/F ratio 0.850 and Vc/F ratio 0.748 (Table S3). The control stream applies both when SEX.EQ.1, and Table S3 labels them 'Females Compared to Males', so the source SEX column is already female = 1.",
      source_name = "SEX"
    ),
    CRCL = list(
      description = "Baseline creatinine clearance by the Cockcroft-Gault equation (raw mL/min, not BSA-normalized).",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL/F, exponent 0.313. The control stream divides by 106.9 mL/min; Table S3 footnote b prints the rounded 107 mL/min. Table S2 reports 'Creatinine Clearance (mL/min)' separately from BSA-normalized eGFR.",
      source_name = "CRCL"
    ),
    DIS_RA = list(
      description = "Rheumatoid arthritis patient indicator (1 = RA, 0 = healthy subject or Crohn's disease patient).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-RA: healthy subjects and CD patients share the CL/F reference).",
      notes = "RA patients have 0.779 times the CL/F of non-RA subjects (Table S3 'Ratio of CL/F in RA Patients Compared to non-RA subjects (healthy subjects or CD patients)'). Also selects the RA interindividual variance on CL/F, and contributes to the patient indicator for the Ka ratio and the Vc/F variance.",
      source_name = "RA (control stream: RA.EQ.1 for RA; RA.GT.0 for any patient)"
    ),
    DIS_CD = list(
      description = "Crohn's disease patient indicator (1 = CD, 0 = healthy subject or RA patient).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy subject or RA patient).",
      notes = "No effect on typical CL/F (the paper reports CD CL/F comparable with healthy subjects). Selects the CD interindividual variance on CL/F (41.2%), and with DIS_RA forms the patient indicator for the Ka ratio of 1.63 and the patient Vc/F variance. The control stream identifies CD by STDY.EQ.13740 (CELEST) and the patient Ka term by RA.GT.0.",
      source_name = "STDY = 13740 (CELEST)"
    ),
    STUDY_UPA_PHASE1 = list(
      description = "Phase 1 study indicator (1 = one of the four phase 1 studies, healthy volunteers plus 14 RA patients in Study 3 substudy 2; 0 = phase 2 RA studies M13-537 and M13-550 or CELEST).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase 2 visit-based sampling).",
      notes = "Selects how the lag time varies: an interindividual random effect in phase 1 (36.6%), and a Box-Cox-transformed inter-occasion random effect per study visit in phase 2. The control stream branches on STDY.EQ.13537 / 13550 / 13740.",
      source_name = "STDY"
    ),
    OCC = list(
      description = "Occasion (study visit) index for the inter-occasion variability on the lag time.",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Consumed by 'etaiov_ltlag ~ ... | OCC' in ini(). The source carries five occasion etas (ETA(4)-ETA(8), '$ABBR REPLACE ETA(OCC_ALAG)'), one per phase 2 study visit. Needed in every event table; for phase 1 records and ER doses the occasion eta does not enter the lag time.",
      source_name = "OCC_ALAG"
    ),
    FORM_UPA_ER = list(
      description = "Upadacitinib extended-release formulation indicator (1 = extended-release tablet, 0 = immediate-release reference).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (immediate release).",
      notes = "Switches the absorption branch to the Table S4 extended-release estimates: first-order Ka 0.100 1/h for 44.6% of the absorbed dose, a 6.11 h zero-order input into the central compartment for the remainder, a common 0.225 h lag time, and 74.2% bioavailability relative to IR. Each ER administration needs two dose records at the same time and amount: one to 'depot' and one to 'central' with rate = -2. The ER control stream uses FORM.EQ.2 for ER, so FORM_UPA_ER = FORM - 1.",
      source_name = "FORM (1 = IR, 2 = ER)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 753L,
    n_studies = 6L,
    age_range = "19-85 years",
    weight_range = "39-152 kg",
    sex_female_pct = 63,
    race_ethnicity = c(White = 84, Black = 7, Asian = 5, Hispanic = 3, Other = 1),
    disease_state = "Healthy adults (phase 1), adults with rheumatoid arthritis (one phase 1 substudy and phase 2 studies M13-537 and M13-550), and adults with moderately-to-severely active Crohn's disease (phase 2 CELEST, N = 180 with PK).",
    dose_range = "Immediate release: single doses 1-48 mg, 3-24 mg twice daily, 24 mg once daily. Extended-release branch from a separate fit to phase 1 study M14-680 in healthy subjects: 15 and 30 mg once daily.",
    regions = "Multinational",
    renal_function = "Creatinine clearance 41-253 mL/min",
    notes = "Table S2 lists 824 subjects, including the 71 of the ER formulation study M14-680; the IR control stream excludes that study (IGNORE=(STDY.EQ.14680)), so the counts and percentages here are Table S2 totals minus the 'Study 4' column. CELEST (NCT02365649) PK samples at weeks 2, 4, 8, 12 and 16. NONMEM 7.3, FOCEI then SAEM with an importance-sampling objective function. The ER branch comes from the separate M14-680 fit (Table S4; packaged on its own as Mohamed_2020_upadacitinib_hv_er), in which all disposition parameters were fixed to Table S3."
  )

  ini({
    # ---- Disposition (Table S3, IR formulation, healthy male reference) -----
    lcl <- log(39.6)
    label("Apparent clearance CL/F, healthy or CD male at CrCL 106.9 mL/min (L/h)") # Table S3 'CL/F (L/h)' = 39.6 (RSE 1.80%)
    lvc <- log(141)
    label("Apparent central volume Vc/F, male at 74 kg (L)") # Table S3 'Vc/F (L)' = 141 (RSE 2.14%)
    lvp <- log(66.3)
    label("Apparent peripheral volume Vp/F (L)") # Table S3 'Vp/F (L)' = 66.3 (RSE 1.10%)
    lq <- log(3.39)
    label("Apparent intercompartmental clearance Q/F (L/h)") # Table S3 'Q/F (L/h)' = 3.39 (RSE 2.42%)

    # ---- IR absorption (Table S3) ----------------------------------------------
    lka <- log(3.60)
    label("Absorption rate constant Ka, IR formulation, healthy subject (1/h)") # Table S3 'KA (1/h)' = 3.60 (RSE 17.8%); footnote a: modeled as exp(theta)
    ltlag <- log(0.469)
    label("Absorption lag time, IR formulation (h)") # Table S3 'ALAG1 (h)' = 0.4690 (RSE 0.0147%)

    # ---- Covariate effects (Table S3 footnotes b-d) -----------------------------
    e_dis_ra_cl <- log(0.779)
    label("Log ratio of CL/F, RA patients vs non-RA subjects (unitless)") # Table S3 'Ratio of CL/F in RA Patients Compared to non-RA subjects' = 0.779
    e_sexf_cl <- log(0.850)
    label("Log ratio of CL/F, female vs male (unitless)") # Table S3 'Ratio of CL/F in Females Compared to Males' = 0.850
    e_crcl_cl <- 0.313
    label("Power exponent of creatinine clearance on CL/F (unitless)") # Table S3 'Exponent for the Effect of CrCL on CL/F' = 0.313
    e_sexf_vc <- log(0.748)
    label("Log ratio of Vc/F, female vs male (unitless)") # Table S3 'Ratio of Vc/F in Females Compared to Males' = 0.748
    e_wt_vc <- 0.490
    label("Power exponent of body weight on Vc/F (unitless)") # Table S3 'Exponent for the Effect of Bodyweight on Vc/F' = 0.490
    e_pat_ka <- log(1.63)
    label("Log ratio of Ka, RA or CD patients vs healthy subjects (unitless)") # Table S3 'Ratio of KA in RA or CD Patients Compared to healthy subjects' = 1.63

    # ---- Lag-time inter-occasion variability shape (Table S3) --------------------
    bc_tlag <- fixed(8)
    label("Box-Cox shape of the inter-occasion variability on the lag time (unitless)") # Table S3 'Box-cox Skewness Parameter for the Distribution of IOV on ALAG1' = 8.00 fix

    # ---- ER absorption (Table S4; fitted to study M14-680) ----------------------
    lka_er <- log(0.100)
    label("First-order absorption rate constant Ka, ER formulation (1/h)") # Table S4 'KA ER (1/h)' = 0.100
    ltlag_er <- log(0.225)
    label("Absorption lag time of both ER input routes (h)") # Table S4 'ALAG1 ER (h)' = 0.225; control stream ALAG1 = ALAG2 = LAGT2 for FORM.EQ.2
    logitffo <- log(0.446 / (1 - 0.446))
    label("Logit of the fraction of the ER dose absorbed by the first-order route (unitless)") # Table S4 'Fraction of Extended-release Dose Absorbed through First-Order Process (%)' = 44.6; footnote d: Frac = exp(theta + eta)/(exp(theta + eta) + 1)
    ld2 <- log(6.11)
    label("Duration of the zero-order ER input into the central compartment (h)") # Table S4 'Zero-Order Infusion duration D2 (h)' = 6.11
    lfrel_er <- log(0.742)
    label("Bioavailability of the ER formulation relative to IR (fraction)") # Table S4 'Bioavailability of the Extended-release Formulation Relative to the Immediate-Release Formulation (%)' = 74.2

    # ---- Interindividual variability ----------------------------------------------
    # The source writes CL = TVCL * EXP(THETA * ETA) with OMEGA = 1 FIX, so each
    # THETA is a standard deviation; Table S3 reports it as %ISV =
    # SQRT(omega^2) x 100. Variances below are (%ISV / 100)^2.
    etalcl_hv ~ 0.026244 # Table S3 'ISV CL in healthy subjects (%)' = 16.2; 0.162^2
    etalcl_ra ~ 0.073984 # Table S3 'ISV CL in subjects with RA (%)' = 27.2; 0.272^2
    etalcl_cd ~ 0.169744 # Table S3 'ISV CL in subjects with CD (%)' = 41.2; 0.412^2
    etalvc_hv ~ 0.017956 # Table S3 'ISV Vc in healthy subjects (%)' = 13.4; 0.134^2
    etalvc_pat ~ 0.081796 # Table S3 'ISV Vc in subjects with RA and CD (%)' = 28.6; 0.286^2
    etalka ~ 1.2769 # Table S3 'ISV KA (%)' = 113; 1.13^2
    etaltlag ~ 0.133956 # Table S3 'ISV ALAG1 (%) for Phase 1 Studies' = 36.6; 0.366^2
    etaiov_ltlag ~ 0.0936 | OCC # Table S3 'omega^2 IOV on ALAG1 (Phase 2 study visits only)' = 0.0936
    etalka_er ~ 1.2321 # Table S4 '%ISV' of 'KA ER' = 111; 1.11^2
    etaltlag_er ~ 0.218089 # Table S4 '%ISV' of 'ALAG1 ER' = 46.7; 0.467^2
    etalogitffo ~ 1.49 # Table S4 'Variance of the ISV for the Fraction of Extended-release Dose Absorbed through First-Order Process' = 1.49
    etald2 ~ 1.4641 # Table S4 '%ISV' of 'Zero-Order Infusion duration D2' = 121; 1.21^2
    etalfrel_er ~ 0.0289 # Table S4 '%ISV' of the relative bioavailability = 17.0; 0.170^2

    # ---- Residual error (Table S3; control stream Y = IPRED*(1+EPS(1)) + EPS(2)) ---
    propSd <- 0.305
    label("Proportional residual error SD (fraction)") # Table S3 'Proportional Residual Error SD (%)' = 30.5
    addSd <- 0.152
    label("Additive residual error SD (ng/mL)") # Table S3 'Additive Residual Error SD (ng/mL)' = 0.1520
  })

  model({
    # Reference covariate values (control stream).
    ref_crcl <- 106.9
    ref_wt <- 74

    # Patient indicator for the Ka ratio and the patient Vc/F variance.
    pat <- DIS_RA + DIS_CD
    is_hv <- 1 - pat

    # ---- Disease-stratified interindividual variability --------------------
    # Only the variance differs between the strata; the typical CL/F is shared
    # by healthy and CD subjects and scaled by e_dis_ra_cl for RA.
    etalcl_eff <- etalcl_hv * is_hv + etalcl_ra * DIS_RA + etalcl_cd * DIS_CD
    etalvc_eff <- etalvc_hv * is_hv + etalvc_pat * pat

    cl <- exp(
      lcl + etalcl_eff + e_dis_ra_cl * DIS_RA + e_sexf_cl * SEXF +
        e_crcl_cl * log(CRCL / ref_crcl)
    )
    vc <- exp(lvc + etalvc_eff + e_sexf_vc * SEXF + e_wt_vc * log(WT / ref_wt))
    vp <- exp(lvp)
    q <- exp(lq)

    # ---- Absorption ----------------------------------------------------------
    # IR lag time: phase 1 subjects carry an interindividual eta; phase 2
    # visits carry a Box-Cox-transformed occasion eta,
    #   TETA = (exp(eta)^SHP - 1) / SHP,  ALAG1 = THETA(6) * exp(TETA).
    iov_bc <- (exp(etaiov_ltlag)^bc_tlag - 1) / bc_tlag
    tlag_ir <- exp(
      ltlag + etaltlag * STUDY_UPA_PHASE1 + iov_bc * (1 - STUDY_UPA_PHASE1)
    )
    ka_ir <- exp(lka + e_pat_ka * pat + etalka)

    # ER branch (Table S4 control stream, FORM.EQ.2): KA, ALAG1 = ALAG2,
    # F1 = Frel * Frac, F2 = Frel * (1 - Frac), D2.
    ka_er <- exp(lka_er + etalka_er)
    tlag_er <- exp(ltlag_er + etaltlag_er)
    ffo_er <- expit(logitffo + etalogitffo)
    frel_er <- exp(lfrel_er + etalfrel_er)

    ka <- ka_er * FORM_UPA_ER + ka_ir * (1 - FORM_UPA_ER)
    tlag <- tlag_er * FORM_UPA_ER + tlag_ir * (1 - FORM_UPA_ER)
    frel <- frel_er * FORM_UPA_ER + (1 - FORM_UPA_ER)
    ffo <- ffo_er * FORM_UPA_ER + (1 - FORM_UPA_ER)
    d2 <- exp(ld2 + etald2)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- ODEs (ADVAN4) --------------------------------------------------------
    # IR doses: one record to depot. ER doses: one record to depot plus one to
    # central with rate = -2 at the same time and amount.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- frel * ffo
    f(central) <- frel * (1 - ffo)
    dur(central) <- d2
    alag(depot) <- tlag
    alag(central) <- tlag

    # mg / L x 1000 = ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
