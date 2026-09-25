Shah_2019_cefotiam <- function() {
  description <- paste(
    "Three-compartment population PK model for IV cefotiam with a",
    "simultaneous fit of total plasma concentrations and the fraction of",
    "dose excreted unchanged in urine. Built from 14 Caucasian adults",
    "(8 cystic fibrosis [CF] patients, 6 healthy volunteers [HVs]) each",
    "given a single 1027.5 mg cefotiam dose as a 3 min IV infusion.",
    "The novel feature of this model is that every disposition parameter",
    "is referenced to the UNBOUND cefotiam concentration, and the observed",
    "TOTAL plasma concentration is reconstructed algebraically as",
    "Cc = Cunbound / fu. Body size and body composition are captured by",
    "allometric scaling on lean body mass (LBM) with fixed exponents 0.75",
    "on all clearance terms and 1.0 on all volumes (reference LBM = 53 kg,",
    "equivalent to a standard 70 kg total body weight). Unbound total",
    "clearance is split into a renal arm (CL_R,u, whose mass is tracked in",
    "the canonical urine compartment) and a non-renal arm (CL_NR,u).",
    "Because the unbound parameters are shared across all four cohorts,",
    "the entire remaining CF-vs-HV and female-vs-male difference is",
    "carried by the unbound fraction fu, which is estimated separately for",
    "each combination of the DIS_CF and SEXF indicators: fu = 0.500",
    "(fixed, HV males), 0.545 (HV females), 0.563 (CF males) and 0.744",
    "(CF females). This replaces the disease-specific scale factors (FCYF)",
    "used by the same group in Bulitta_2011_cefpirome.R and",
    "Bulitta_2007_piperacillin.R; the authors removed all FCYF terms when",
    "the unbound fractions were estimated."
  )
  reference <- paste(
    "Shah NR, Bulitta JB, Kinzig M, Landersdorfer CB, Jiao Y, Sutaria DS,",
    "Tao X, Hohl R, Holzgrabe U, Kees F, Stephan U, Sorgel F. Novel",
    "Population Pharmacokinetic Approach to Explain the Differences",
    "between Cystic Fibrosis Patients and Healthy Volunteers via Protein",
    "Binding. Pharmaceutics. 2019 Jun 18;11(6):286.",
    "doi:10.3390/pharmaceutics11060286"
  )
  vignette <- "Shah_2019_cefotiam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the authors' own Berkeley Madonna
  # Monte Carlo source listing in the Supplementary Materials (states Cent,
  # Shal, Deep and Urin), which is the authoritative implementation of the
  # final model. The three disposition states carry the TOTAL amount of
  # cefotiam; it is the VOLUMES that are unbound-referenced, so that
  # central / vc is the UNBOUND plasma concentration.
  compartmentData <- list(
    central = list(analyte = "cefotiam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefotiam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "cefotiam", units = "mg", specimen = "plasma", verified = TRUE),
    urine = list(analyte = "cefotiam", units = "mg", specimen = "urine", verified = TRUE)
  )

  covariateData <- list(
    LBM = list(
      description = "Lean body mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling on CL_R,u, CL_NR,u, CLd_shallow,u and",
        "CLd_deep,u (fixed exponent 0.75, Shah 2019 Equation 2) and on V1u,",
        "V2u and V3u (fixed exponent 1.0, Shah 2019 Equation 1), with a",
        "reference LBM_STD of 53 kg stated to be equivalent to a standard",
        "total body weight of 70 kg (Shah 2019 Section 2.6.3 and the",
        "Table 4 caption). LBM was computed by the formula of Cheymol and",
        "James (Shah 2019 Table 1 footnote a); the cohort medians were",
        "40.3 kg (CF, range 28.8-46.2) and 50.6 kg (HV, range 44.6-65.4).",
        "The authors compared five body-size models (none, linear WT,",
        "allometric WT, linear LBM, allometric LBM) and retained",
        "allometric LBM because it gave disease factors closest to 1.0",
        "(Shah 2019 Table 3)."
      ),
      source_name = "LBM"
    ),
    DIS_CF = list(
      description = paste(
        "Cystic-fibrosis cohort indicator: 1 = adult CF patient, 0 = healthy",
        "adult volunteer (reference). Time-fixed per subject."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "The paper does not assign a column name to the CF/HV cohort flag;",
        "in the authors' Berkeley Madonna listing it is the variable GRP",
        "(1 = CF, otherwise HV). In this final model the cohort indicator",
        "does NOT act on any clearance or volume: all disposition",
        "parameters are shared on the unbound scale. DIS_CF acts ONLY in",
        "combination with SEXF to select which of the four estimated",
        "unbound fractions applies (Shah 2019 Section 2.6.5 and Table 4).",
        "This is the deliberate contrast with Bulitta_2011_cefpirome.R,",
        "where the same group carried the CF effect as FCYF scale factors",
        "on CL_R, CL_NR and the volumes."
      ),
      source_name = NA_character_
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "In the authors' Berkeley Madonna listing this is the variable",
        "SEX1F_0M with exactly the canonical orientation (1 = female,",
        "0 = male), so no transformation is needed. Sex enters only",
        "through the unbound fraction: within each cohort the female fu is",
        "higher than the male fu (HV 0.545 vs 0.500; CF 0.744 vs 0.563).",
        "The supplement rationalises this as lower albumin concentrations",
        "in females and in CF patients. The cohort was 7 female / 7 male",
        "(4/4 CF, 3/3 HV; Shah 2019 Table 1)."
      ),
      source_name = "SEX1F_0M"
    ),
    DOSE_CEFOTIAM_MG = list(
      description = "Administered cefotiam dose in milligrams for the current dosing interval.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Required only to express the second fitted endpoint, urinePct,",
        "which Shah 2019 Section 2.6.5 defines as the FRACTION of the dose",
        "excreted into urine as unchanged cefotiam (reported in percent",
        "throughout the paper, e.g. Table 2 medians of 70.3% in CF and",
        "66.3% in HV). Carrying the dose as an explicit covariate keeps the",
        "published additive residual SD of 0.384% attached to the paper's",
        "own observation scale and keeps the endpoint dose-general. Set",
        "DOSE_CEFOTIAM_MG = 1027.5 to reproduce the study; the urine",
        "compartment itself still holds the cumulative amount in mg, so a",
        "user who prefers an amount endpoint can read the state directly."
      ),
      source_name = NA_character_
    )
  )

  population <- list(
    species = "human",
    n_subjects = 14L,
    n_studies = 1L,
    age_range = "17-26 years (CF median 19, range 17-24; HV median 23.5, range 21-26)",
    weight_range = "33.0-80.0 kg total body weight (CF median 45.5, range 33.0-59.0; HV median 68.5, range 58.0-80.0)",
    lbm_range = "28.8-65.4 kg lean body mass (CF median 40.3, range 28.8-46.2; HV median 50.6, range 44.6-65.4)",
    height_range = "157-190 cm (CF median 167, range 157-173; HV median 169, range 164-190)",
    bmi_range = "13.4-27.9 kg/m^2 (CF median 17.0, range 13.4-19.9; HV median 22.5, range 20.3-27.9)",
    sex_female_pct = 50,
    race_ethnicity = "100% Caucasian (Shah 2019 Section 2.1).",
    disease_state = paste(
      "Two parallel groups: 8 patients with cystic fibrosis (4 female,",
      "4 male) and 6 healthy volunteers (3 female, 3 male). The CF",
      "patients were markedly smaller and leaner than the healthy",
      "volunteers (median BMI 17.0 vs 22.5 kg/m^2). One CF patient was 17",
      "years old and was enrolled with consent from a legal representative."
    ),
    dose_range = "Single 1027.5 mg cefotiam dose given as a 3 min IV infusion via a syringe-driver perfusor.",
    regions = "Germany (single-centre study; ethics approval University Hospital Essen, 1984).",
    notes = paste(
      "Baseline demographics from Shah 2019 Table 1. Plasma was sampled",
      "pre-dose, at the end of the 3 min infusion, and at 5, 10, 15, 20,",
      "30, 45, 60 and 90 min plus 2, 3, 4, 5, 6, 8, 12 and 24 h after the",
      "end of infusion; urine was collected over 0-1, 1-2, 2-3, 3-4, 4-5,",
      "5-6, 6-8, 8-12 and 12-24 h from the start of infusion. Plasma and",
      "urine data were fitted SIMULTANEOUSLY, which is what identifies the",
      "renal and non-renal clearance arms separately. Estimation used the",
      "importance-sampling algorithm (pmethod = 4) of S-ADAPT 1.57 driven",
      "by SADAPT-TRAN. The parameter values encoded here are the final",
      "model of Shah 2019 Table 4; where the authors' Berkeley Madonna",
      "listing in the Supplementary Materials prints more significant",
      "figures than Table 4, the listing's value is used and the Table 4",
      "rounding is given alongside it. Two alternative parameterisations",
      "reported in the Supplementary Materials are NOT encoded here",
      "because the paper does not present them as final: Table S1 (disease",
      "specific FCYF scale factors with fu fixed to 0.5 in every group) and",
      "Table S2 (renal clearance split into a fixed 7.2 L/h glomerular",
      "filtration arm plus an unbound-fraction-independent tubular",
      "secretion arm, which had a -2x log-likelihood worse by 7.1).",
      "Simulate the study dose with rate = amt / (3/60) or dur = 3/60 to",
      "encode the 3 min zero-order infusion."
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # UNBOUND DISPOSITION PARAMETERS -- Shah 2019 Table 4. Every clearance and
    # volume below refers to the UNBOUND cefotiam concentration and applies to
    # ALL four subject groups at the reference LBM of 53 kg; the final model
    # carries no disease- or sex-specific scale factor on any of them. Typical
    # values are taken from the authors' Berkeley Madonna listing in the
    # Supplementary Materials, which prints the estimates to six significant
    # figures; the Table 4 rounding is quoted in each trailing comment.
    # ------------------------------------------------------------------------
    lcl_renal <- log(23.7911)
    label("Unbound renal clearance CL_R,u at LBM = 53 kg (L/h)")
    # Supplement Madonna listing Mean_CLR = 23.7911; Table 4 row 'Unbound renal clearance' = 23.8 L/h (SE 6.9%)
    lcl_nonren <- log(10.9943)
    label("Unbound non-renal clearance CL_NR,u at LBM = 53 kg (L/h)")
    # Supplement Madonna listing Mean_CLNR = 10.9943; Table 4 row 'Unbound nonrenal clearance' = 11.0 L/h (SE 7.0%)
    lvc <- log(15.5927)
    label("Unbound central volume V1u at LBM = 53 kg (L)")
    # Supplement Madonna listing Mean_V1 = 15.5927; Table 4 row 'Unbound volume ... central compartment' = 15.6 L (SE 6.5%)
    lvp <- log(6.90932)
    label("Unbound shallow peripheral volume V2u at LBM = 53 kg (L)")
    # Supplement Madonna listing Mean_V2 = 6.90932; Table 4 row 'Unbound volume ... shallow peripheral compartment' = 6.91 L (SE 14.1%)
    lvp2 <- log(4.5646)
    label("Unbound deep peripheral volume V3u at LBM = 53 kg (L)")
    # Supplement Madonna listing Mean_V3 = 4.5646; Table 4 row 'Unbound volume ... deep peripheral compartment' = 4.56 L (SE 16.4%)
    lq <- log(13.8043)
    label("Unbound distribution CL central <-> shallow peripheral, CLd_shallow,u (L/h)")
    # Supplement Madonna listing Mean_CLD = 13.8043; Table 4 row 'Unbound distribution clearance for shallow peripheral compartment' = 13.8 L/h (SE 15.0%)
    lq2 <- log(1.83622)
    label("Unbound distribution CL central <-> deep peripheral, CLd_deep,u (L/h)")
    # Supplement Madonna listing Mean_CLD3 = 1.83622; Table 4 row 'Unbound distribution clearance for deep peripheral compartment' = 1.84 L/h (SE 26.1%)

    # ------------------------------------------------------------------------
    # ALLOMETRIC SCALING ON LEAN BODY MASS -- Shah 2019 Section 2.6.3:
    # 'The allometric body size models used a fixed exponent of 1.0 (i.e.,
    # linear scaling) for volume of distribution and a fixed exponent of 0.75
    # ... for clearances.' Equation 1 gives FSize,V = LBM / LBM_STD and
    # Equation 2 gives FSize,CL = (LBM / LBM_STD)^0.75, with LBM_STD = 53 kg.
    # Reproduced verbatim by the supplement Madonna listing as
    # FWTCL = (LBM/53)**0.75 and FWTV = (LBM/53).
    # ------------------------------------------------------------------------
    e_lbm_cl_q <- fixed(0.75)
    label("Allometric (LBM) exponent on CL_R,u, CL_NR,u, CLd_shallow,u, CLd_deep,u (unitless)")
    # Shah 2019 Equation 2 and Section 2.6.3 (exponent fixed, not estimated)
    e_lbm_vc_vp <- fixed(1.00)
    label("Allometric (LBM) exponent on V1u, V2u, V3u (unitless)")
    # Shah 2019 Equation 1 and Section 2.6.3 (exponent fixed, not estimated)

    # ------------------------------------------------------------------------
    # UNBOUND FRACTION IN PLASMA -- Shah 2019 Table 4, four separately
    # estimated population means selected by cohort (DIS_CF) and sex (SEXF).
    # Table 4 footnote d: 'Unbound fraction was fixed to 0.5 in male healthy
    # volunteers based on literature data. The population means of the
    # remaining three unbound fractions were estimated separately for males
    # and females with a small fixed between subject variability (5%
    # coefficient of variation).' The four values are the ONLY thing that
    # distinguishes the four subject groups in this model. Separate stratified
    # typical values follow the shipped lvc_male / lvc_female pattern of
    # Kim_2025_infliximab_ternant.R rather than a reference-plus-effect
    # reparameterisation, so each published estimate appears verbatim.
    # ------------------------------------------------------------------------
    lfu_hv_male <- fixed(log(0.5))
    label("Unbound fraction in plasma, male healthy volunteers (unitless)")
    # Shah 2019 Table 4 row 'Unbound fraction in plasma for male healthy volunteers' = 0.50 (fixed); supplement Madonna listing FU_HVM = 0.5
    lfu_hv_female <- log(0.544842)
    label("Unbound fraction in plasma, female healthy volunteers (unitless)")
    # Supplement Madonna listing Mean_FU_HVF = 0.544842; Table 4 row 'Unbound fraction in plasma for female healthy volunteers' = 0.545 (SE 13.6%)
    lfu_cf_male <- log(0.562527)
    label("Unbound fraction in plasma, male cystic fibrosis patients (unitless)")
    # Supplement Madonna listing Mean_FU_CFM = 0.562527; Table 4 row 'Unbound fraction in plasma for males with CF' = 0.563 (SE 13.5%)
    lfu_cf_female <- log(0.74359)
    label("Unbound fraction in plasma, female cystic fibrosis patients (unitless)")
    # Supplement Madonna listing Mean_FU_CFF = 0.74359; Table 4 row 'Unbound fraction in plasma for females with CF' = 0.744 (SE 4.5%)

    # ------------------------------------------------------------------------
    # BETWEEN-SUBJECT VARIABILITY -- Shah 2019 Table 4, BSV column. Section
    # 2.6.4 states that the BSV for clearances and volumes was log-normal and
    # that eta_BSV was 'a normally distributed random variable with mean zero
    # and standard deviation BSV', and Table 4 footnote a calls the reported
    # figure an 'apparent coefficient of variation of a normal distribution on
    # natural logarithmic scale'. Both readings agree that the tabulated
    # number IS the log-scale SD, so omega^2 = BSV^2 (NOT log(1 + CV^2)); the
    # supplement Madonna listing confirms it by drawing normal(0, CV_<param>).
    # The SDs below are the Madonna listing's six-figure values; the Table 4
    # rounding is quoted in each trailing comment.
    # ------------------------------------------------------------------------
    etalcl_renal ~ 0.0559564 # Madonna CV_CLR = 0.236551; Table 4 BSV(CLr,u) = 0.237 (SE 52.7%) -> 0.236551^2
    etalcl_nonren ~ 0.0561837 # Madonna CV_CLNR = 0.237031; Table 4 BSV(CLnr,u) = 0.237 (SE 50.2%) -> 0.237031^2
    etalvc ~ 0.0358917 # Madonna CV_V1 = 0.189451; Table 4 BSV(V1u) = 0.189 (SE 74.0%) -> 0.189451^2
    etalvp ~ 0.0656651 # Madonna CV_V2 = 0.256252; Table 4 BSV(V2u) = 0.256 (SE 88.2%) -> 0.256252^2
    etalvp2 ~ 0.2033234 # Madonna CV_V3 = 0.450914; Table 4 BSV(V3u) = 0.451 (SE 131%) -> 0.450914^2
    etalq ~ 0.1729113 # Madonna CV_CLD = 0.415826; Table 4 BSV(CLd shallow,u) = 0.416 (SE 183%) -> 0.415826^2
    etalq2 ~ 0.0952290 # Madonna CV_CLD3 = 0.308592; Table 4 BSV(CLd deep,u) = 0.309 (SE 83.8%) -> 0.308592^2

    # The three ESTIMATED unbound fractions each carry a BSV that the authors
    # FIXED at a 5% coefficient of variation (Table 4 footnote d), so
    # omega^2 = 0.05^2 = 0.0025. The unbound fraction of male healthy
    # volunteers was itself fixed and carries no random effect at all -- the
    # supplement Madonna listing draws ETA_FU_CFF, ETA_FU_CFM and ETA_FU_HVF
    # but has no ETA_FU_HVM.
    etalfu_hv_female ~ fixed(0.0025) # Table 4 footnote d: BSV held at a 5% CV -> 0.05^2
    etalfu_cf_male ~ fixed(0.0025) # Table 4 footnote d: BSV held at a 5% CV -> 0.05^2
    etalfu_cf_female ~ fixed(0.0025) # Table 4 footnote d: BSV held at a 5% CV -> 0.05^2

    # ------------------------------------------------------------------------
    # RESIDUAL ERROR -- Shah 2019 Table 4. Section 2.6.5: 'The residual
    # unidentified variability was described by a combined additive plus
    # proportional residual error model for plasma concentrations. The
    # fractions of dose excreted into urine as unchanged cefotiam were fit
    # using an additive residual error model.' Note that the Table 4 caption
    # exempts the additive residual errors from the 'refers to unbound
    # cefotiam' statement: both plasma residual terms act on the TOTAL plasma
    # concentration Cc, which is the quantity that was actually assayed.
    # ------------------------------------------------------------------------
    addSd <- 0.0186
    label("Additive residual SD, total plasma cefotiam (mg/L)")
    # Shah 2019 Table 4 row 'SD of additive residual error for plasma concentrations' SDin = 0.0186 mg/L (SE 53.7%)
    propSd <- 0.166
    label("Proportional residual SD, total plasma cefotiam (fraction)")
    # Shah 2019 Table 4 row 'Proportional residual error for plasma concentrations' SDsl = 0.166 (SE 7.8%)
    addSd_urinePct <- 0.384
    label("Additive residual SD, cumulative percent of dose excreted in urine (%)")
    # Shah 2019 Table 4 row 'SD of additive residual error for fraction of dose in urine' UDin = 0.384% (SE 76.1%)
  })

  model({
    # ---------------------------------------------------------------------
    # Reference body size. Shah 2019 Section 2.6.3 and Table 4 caption:
    # LBM_STD = 53 kg, 'equivalent to a standard weight of 70 kg'.
    # ---------------------------------------------------------------------
    lbm_ref <- 53 # kg

    # ---------------------------------------------------------------------
    # Individual unbound disposition parameters. Shah 2019 Equation 4:
    #   CLru,i = CLru * FSize,CL,i * exp(eta_BSV_CLr,u,i)
    # i.e. allometric LBM scaling and a log-normal random effect, with NO
    # disease-specific scale factor -- Section 2.6.5 states that 'all disease
    # specific scale factors FCYF were removed from the model' once the
    # unbound fractions were estimated. The same form applies to every
    # clearance and, with the volume exponent, to every volume.
    # ---------------------------------------------------------------------
    cl_renal <- exp(lcl_renal + etalcl_renal) * (LBM / lbm_ref)^e_lbm_cl_q
    cl_nonren <- exp(lcl_nonren + etalcl_nonren) * (LBM / lbm_ref)^e_lbm_cl_q
    q <- exp(lq + etalq) * (LBM / lbm_ref)^e_lbm_cl_q
    q2 <- exp(lq2 + etalq2) * (LBM / lbm_ref)^e_lbm_cl_q
    vc <- exp(lvc + etalvc) * (LBM / lbm_ref)^e_lbm_vc_vp
    vp <- exp(lvp + etalvp) * (LBM / lbm_ref)^e_lbm_vc_vp
    vp2 <- exp(lvp2 + etalvp2) * (LBM / lbm_ref)^e_lbm_vc_vp

    # ---------------------------------------------------------------------
    # Unbound fraction: the four-way cohort-by-sex selector. This reproduces
    # the nested IF cascade of the supplement Madonna listing
    #   FU = IF (GRP = 1) THEN IF (SEX1F_0M = 0) THEN FU_CFM ELSE FU_CFF
    #                     ELSE IF (SEX1F_0M = 0) THEN FU_HVM ELSE FU_HVF
    # written as an indicator-weighted sum so it is differentiable for the
    # solver. Exactly one of the four products is non-zero for any subject.
    # Each group's unbound fraction gets its own simple line so that rxode2
    # can mu-reference the three random effects.
    # ---------------------------------------------------------------------
    fu_hv_male <- exp(lfu_hv_male)
    fu_hv_female <- exp(lfu_hv_female + etalfu_hv_female)
    fu_cf_male <- exp(lfu_cf_male + etalfu_cf_male)
    fu_cf_female <- exp(lfu_cf_female + etalfu_cf_female)
    fu <- (1 - DIS_CF) * ((1 - SEXF) * fu_hv_male + SEXF * fu_hv_female) +
      DIS_CF * ((1 - SEXF) * fu_cf_male + SEXF * fu_cf_female)

    # ---------------------------------------------------------------------
    # Micro-constants.
    # ---------------------------------------------------------------------
    kel_renal <- cl_renal / vc
    kel_nonren <- cl_nonren / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ---------------------------------------------------------------------
    # Three-compartment disposition with parallel renal plus non-renal
    # elimination and cumulative urinary excretion, transcribed from the
    # supplement Madonna listing states Cent / Shal / Deep / Urin. IV doses
    # target central directly; supply rate = amt / (3/60) or dur = 3/60 to
    # encode the study's 3 min zero-order infusion.
    # ---------------------------------------------------------------------
    d / dt(central) <- -(kel_renal + kel_nonren) * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1
    d / dt(peripheral2) <- k13 * central - k31 * peripheral2
    d / dt(urine) <- kel_renal * central

    # ---------------------------------------------------------------------
    # Observations.
    # Cunbound: unbound plasma concentration (mg/L). The volumes are
    #           unbound-referenced, so central / vc is already unbound.
    # Cc      : TOTAL plasma concentration (mg/L), the assayed quantity.
    #           Shah 2019 Section 2.6.5: 'The observed plasma concentration
    #           of total cefotiam was calculated as the modelled unbound
    #           cefotiam concentration divided by the unbound fraction.'
    #           Supplement Madonna listing: C1 = Cent / V1cov / FU.
    # urinePct: cumulative percent of the administered dose recovered in
    #           urine as unchanged cefotiam (Section 2.6.5; Table 2 reports
    #           medians of 70.3% in CF and 66.3% in HV).
    # ---------------------------------------------------------------------
    Cunbound <- central / vc
    Cc <- Cunbound / fu
    urinePct <- 100 * urine / DOSE_CEFOTIAM_MG

    Cc ~ add(addSd) + prop(propSd)
    urinePct ~ add(addSd_urinePct)
  })
}
