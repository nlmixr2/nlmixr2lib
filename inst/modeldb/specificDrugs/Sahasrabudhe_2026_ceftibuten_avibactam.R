Sahasrabudhe_2026_ceftibuten_avibactam <- function() {
  description <- paste(
    "Joint (simultaneously fitted) oral population PK model for cis-ceftibuten",
    "and avibactam after oral ceftibuten and the avibactam prodrug avibactam",
    "tomilopil, given alone or together, in healthy adults from three Phase 1",
    "studies (Sahasrabudhe 2026, final model). Ceftibuten: one-compartment",
    "disposition with a depot plus five transit compartments sharing one",
    "transit rate constant. Avibactam: two-compartment disposition with a depot",
    "plus three transit compartments. Body weight enters allometrically with",
    "fixed exponents (0.75 on CL and Q, 1 on Vc and Vp) and BSA-normalised",
    "creatinine clearance as a power function on both clearances. Fasting,",
    "high-fat meal, suspension formulation, single-agent dosing and dose",
    "amount shift the transit rate constants, relative bioavailability and",
    "avibactam clearance relative to co-administration with a regular meal.",
    "The avibactam tomilopil dose is entered in mg of prodrug and converted to",
    "avibactam equivalents inside the model. A four-way OMEGA block links",
    "ceftibuten and avibactam CL and Vc. Ceftibuten is the unsuffixed parent;",
    "avibactam carries the sibling-drug suffix _avi throughout. See",
    "modellib('Sahasrabudhe_2026_ceftibuten_avibactam_crcl') for the",
    "simulation model the authors used to select the Phase 3 dose."
  )
  reference <- paste(
    "Sahasrabudhe SA, Shoji S, Jakate AS, Abraham BK, Toussi SS, Maligalig BJ,",
    "Darekar AS, Tripathy S, Kantaridis C, Mfopou JK, Selivanova O, Bergman A,",
    "Alami NN, Singh RSP. Dose Selection for Ceftibuten Combined with",
    "Avibactam Tomilopil, a Potential Oral Alternative to Carbapenems for",
    "Patients with Complicated Urinary Tract Infections. J Clin Pharmacol.",
    "2026;66(9):e70263. doi:10.1002/jcph.70263.",
    "All parameter estimates are from Supplemental Appendix Table S12",
    "'Parameter Estimates for CTB/AVI-PO Joint Pharmacokinetic Model (Final",
    "Model)'; covariate forms from the Table S12 footnote and Methods, checked",
    "against the Figure S4 (top) forest plot.",
    sep = " "
  )
  vignette <- "Sahasrabudhe_2026_ceftibuten_avibactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "ceftibuten", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "ceftibuten", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "ceftibuten", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "ceftibuten", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "ceftibuten", units = "mg", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "ceftibuten", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cis-ceftibuten", units = "mg", specimen = "plasma", verified = TRUE),
    depot_avi = list(analyte = "avibactam", units = "mg", specimen = "administration site", verified = TRUE),
    transit1_avi = list(analyte = "avibactam", units = "mg", specimen = "administration site", verified = TRUE),
    transit2_avi = list(analyte = "avibactam", units = "mg", specimen = "administration site", verified = TRUE),
    transit3_avi = list(analyte = "avibactam", units = "mg", specimen = "administration site", verified = TRUE),
    central_avi = list(analyte = "avibactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_avi = list(analyte = "avibactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  # Ceftibuten doses go to depot (mg ceftibuten); avibactam tomilopil doses go
  # to depot_avi in mg of PRODRUG (converted to avibactam equivalents by
  # f(depot_avi)).
  dosing <- c("depot", "depot_avi")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling on CL and Q (fixed exponent 0.75) and on Vc and",
        "Vp (fixed exponent 1) of both drugs, reference 70 kg (Methods,",
        "'Base Population Models'; Table S12 footnote). Analysis-set median",
        "76.6 kg, range 50.1-101 kg (Table 2)."
      ),
      source_name = "WT"
    ),
    CRCL = list(
      description = paste(
        "BSA-normalised Cockcroft-Gault creatinine clearance (nCrCL) in",
        "mL/min/1.73 m^2"
      ),
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on ceftibuten and avibactam CL/F,",
        "(nCrCL/112)^theta (Table S12 footnote). CrCL is the Cockcroft-Gault",
        "estimate normalised to 1.73 m^2 body surface area; the authors",
        "normalised it so that renal function is separated from body size,",
        "which is already carried by the allometric weight terms. The",
        "centring value 112 is the value printed in Methods and Table S12;",
        "Table 2 reports an analysis-set median nCrCL of 109 (range",
        "68.3-160) mL/min/1.73 m^2."
      ),
      source_name = "nCrCL"
    ),
    FED = list(
      description = "Fed state at dosing (1 = fed, 0 = fasted)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (fed; the model reference is co-administration with a regular meal)",
      notes = paste(
        "The model reference is FED = 1 with FED_HIGHFAT = 0 (a regular",
        "meal). The paper's 'fasted' effects enter as (1 - FED): on the",
        "transit rate constant and relative bioavailability of both drugs and",
        "on avibactam CL/F (Table S12). Fasted records: Study C4691001 Part 2",
        "Day 7 and Study C4691003 treatments A and B. May vary between",
        "dosing occasions of the same subject."
      ),
      source_name = "fasted"
    ),
    FED_HIGHFAT = list(
      description = "Dosed after a high-fat, high-calorie meal (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (regular meal or fasted)",
      notes = paste(
        "High-fat-meal records also have FED = 1. In this final model the",
        "high-fat meal shifts the ceftibuten transit rate constant and",
        "relative bioavailability and the avibactam transit rate constant",
        "(Table S12). Only Study C4691003 treatment C used a high-fat meal;",
        "the other fed records used a regular (standard) meal."
      ),
      source_name = "High-fat meal"
    ),
    FORM_AVIBACTAMTOMILOPIL_SUSP = list(
      description = "Avibactam tomilopil given as the oral suspension (1) rather than a capsule or tablet (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (capsule or fixed-dose combination tablet)",
      notes = paste(
        "Raises the avibactam transit rate constant and relative",
        "bioavailability (Table S12 rows 'SUSP on Ktr-AVI' and 'SUSP on",
        "rBA-AVI'). The suspension was used only in Study C4691002, where",
        "avibactam tomilopil was also always given alone, so this effect and",
        "the single-agent effect on avibactam CL/F cannot be separated",
        "(Results, 'Final Population Pharmacokinetic Model'). Tablet and",
        "capsule were not different (formulation was not retained)."
      ),
      source_name = "SUSP"
    ),
    CONMED_CEFTIBUTEN = list(
      description = "Ceftibuten co-administered with the avibactam tomilopil dose (1 = yes, 0 = avibactam tomilopil alone)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (co-administered)",
      notes = paste(
        "The paper's 'AVP alone' effect on avibactam CL/F enters as",
        "(1 - CONMED_CEFTIBUTEN): avibactam CL/F is 27% higher when avibactam",
        "tomilopil is given without ceftibuten (Results; Table S12). Set to 1",
        "for subjects who receive both drugs."
      ),
      source_name = "AVP alone"
    ),
    CONMED_AVIBACTAMTOMILOPIL = list(
      description = "Avibactam tomilopil co-administered with the ceftibuten dose (1 = yes, 0 = ceftibuten alone)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (co-administered)",
      notes = paste(
        "The paper's 'CTB alone' effect on ceftibuten relative",
        "bioavailability enters as (1 - CONMED_AVIBACTAMTOMILOPIL):",
        "ceftibuten relative bioavailability is 34% higher when ceftibuten is",
        "given without avibactam tomilopil (Results; Table S12). Set to 1 for",
        "subjects who receive both drugs."
      ),
      source_name = "CTB alone"
    ),
    DOSE_CEFTIBUTEN_MG = list(
      description = "Ceftibuten dose of the current administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on ceftibuten relative bioavailability,",
        "(DOSE/400)^theta, reference 400 mg (Methods; Table S12 footnote).",
        "Doses of 400-1600 mg were studied."
      ),
      source_name = "CTB DOSE"
    ),
    DOSE_AVIBACTAMTOMILOPIL_MG = list(
      description = "Avibactam tomilopil (prodrug) dose of the current administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on avibactam relative bioavailability,",
        "(DOSE/1350)^theta, reference 1350 mg of PRODRUG (Methods; Table S12",
        "footnote). Doses of 300-1350 mg were studied. Should equal the amount",
        "dosed into depot_avi."
      ),
      source_name = "AVP DOSE"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 71L,
    n_studies = 3L,
    age_range = "18-59 years",
    age_median = "32 years",
    weight_range = "50.1-101 kg",
    weight_median = "76.6 kg",
    sex_female_pct = 19.7,
    race_ethnicity = c(White = 62, Asian = 22.5, Black = 5.6, Other = 7, `Not reported` = 2.8),
    disease_state = "Healthy adult volunteers",
    dose_range = paste(
      "Ceftibuten 400-1600 mg single doses and 400-800 mg q8h for 7 days;",
      "avibactam tomilopil 300-1350 mg single doses and 900-1350 mg q8h for",
      "7 days; oral capsules, oral suspension (avibactam tomilopil only) or a",
      "ceftibuten 402 mg + avibactam tomilopil 1200 mg fixed-dose tablet"
    ),
    regions = "Australia, Belgium",
    renal_function = "CrCL 58.5-179 mL/min (median 124); nCrCL 68.3-160 mL/min/1.73 m^2 (median 109)",
    notes = paste(
      "Studies C4691002 (single ascending avibactam tomilopil suspension,",
      "fed), C4691001 (single and multiple ascending ceftibuten +",
      "avibactam tomilopil capsules, regular meal; Part 2 Day 7 fasted) and",
      "C4691003 (relative bioavailability of capsule vs tablet, fasted and",
      "high-fat fed). 1467 cis-ceftibuten and 1707 avibactam plasma",
      "observations above the LLOQ. Demographics from Table 2; study designs",
      "from Table 1."
    )
  )

  ini({
    # ---- Ceftibuten (unsuffixed parent), Table S12 ----
    lcl <- log(4.24); label("Ceftibuten apparent clearance CL/F at WT 70 kg, nCrCL 112, co-administration with a regular meal (L/h)") # Table S12: 'theta CL/F-CTB (L/h)' = 4.24 (RSE 4.3%)
    lvc <- log(17.9); label("Ceftibuten apparent central volume Vc/F at WT 70 kg (L)") # Table S12: 'theta Vc/F-CTB (L)' = 17.9 (RSE 4.17%)
    lktr <- log(2.55); label("Ceftibuten transit rate constant, regular meal (1/h)") # Table S12: 'theta ktr-CTB (1/h)' = 2.55 (RSE 7.32%)

    e_crcl_cl <- 0.272; label("Power exponent of nCrCL/112 on ceftibuten CL/F (unitless)") # Table S12: 'theta nCRCL on CL/F-CTB' = 0.272 (RSE 37.1%)
    e_fasted_ktr <- 1.62; label("Proportional shift in ceftibuten ktr when fasted (fraction)") # Table S12: 'theta fasted on Ktr-CTB' = 1.62 (RSE 16.4%)
    e_fasted_fdepot <- -0.187; label("Proportional shift in ceftibuten relative bioavailability when fasted (fraction)") # Table S12: 'theta fasted on rBA-CTB' = -0.187 (RSE 21.9%)
    e_mono_fdepot <- 0.343; label("Proportional shift in ceftibuten relative bioavailability when given alone (fraction)") # Table S12: 'theta CTB alone on rBA-CTB' = 0.343 (RSE 38.6%)
    e_dose_fdepot <- -0.211; label("Power exponent of ceftibuten dose/400 mg on ceftibuten relative bioavailability (unitless)") # Table S12: 'theta CTB DOSE on rBA-CTB' = -0.211 (RSE 33.1%)
    e_highfat_ktr <- 1.02; label("Proportional shift in ceftibuten ktr after a high-fat meal (fraction)") # Table S12: 'theta High-fat meal on Ktr-CTB' = 1.02 (RSE 23.6%)
    e_highfat_fdepot <- -0.255; label("Proportional shift in ceftibuten relative bioavailability after a high-fat meal (fraction)") # Table S12: 'theta High-fat meal on rBA-CTB' = -0.255 (RSE 17.4%)

    # ---- Avibactam (suffix _avi), Table S12 ----
    lcl_avi <- log(21.1); label("Avibactam apparent clearance CL/F at WT 70 kg, nCrCL 112, co-administration with a regular meal (L/h)") # Table S12: 'theta CL/F-AVI (L/h)' = 21.1 (RSE 4.16%)
    lvc_avi <- log(73.2); label("Avibactam apparent central volume Vc/F at WT 70 kg (L)") # Table S12: 'theta Vc/F-AVI (L)' = 73.2 (RSE 7.3%)
    lq_avi <- log(0.34); label("Avibactam apparent inter-compartmental clearance Q/F at WT 70 kg (L/h)") # Table S12: 'theta Q/F-AVI (L/h)' = 0.34 (RSE 11.3%)
    lvp_avi <- log(12.1); label("Avibactam apparent peripheral volume Vp/F at WT 70 kg (L)") # Table S12: 'theta Vp/F-AVI (L)' = 12.1 (RSE 28.9%)
    lktr_avi <- log(1.65); label("Avibactam transit rate constant, capsule with a regular meal (1/h)") # Table S12: 'theta ktr-AVI (1/h)' = 1.65 (RSE 7.75%)

    e_crcl_cl_avi <- 0.316; label("Power exponent of nCrCL/112 on avibactam CL/F (unitless)") # Table S12: 'theta nCrCL on CL/F-AVI' = 0.316 (RSE 34%)
    e_fasted_ktr_avi <- 2.32; label("Proportional shift in avibactam ktr when fasted (fraction)") # Table S12: 'theta fasted on Ktr-AVI' = 2.32 (RSE 13.4%)
    e_fasted_fdepot_avi <- -0.556; label("Proportional shift in avibactam relative bioavailability when fasted (fraction)") # Table S12: 'theta fasted on rBA-AVI' = -0.556 (RSE 6.52%)
    e_susp_ktr_avi <- 2.37; label("Proportional shift in avibactam ktr for the avibactam tomilopil suspension (fraction)") # Table S12: 'theta SUSP on Ktr-AVI' = 2.37 (RSE 16.9%)
    e_susp_fdepot_avi <- 0.384; label("Proportional shift in avibactam relative bioavailability for the suspension (fraction)") # Table S12: 'theta SUSP on rBA-AVI' = 0.384 (RSE 41.3%)
    e_mono_cl_avi <- 0.272; label("Proportional shift in avibactam CL/F when avibactam tomilopil is given alone (fraction)") # Table S12: 'theta AVP alone on CL/F-AVI' = 0.272 (RSE 30.5%)
    e_dose_fdepot_avi <- 0.184; label("Power exponent of avibactam tomilopil dose/1350 mg on avibactam relative bioavailability (unitless)") # Table S12: 'theta AVP DOSE on rBA-AVI' = 0.184 (RSE 72.4%)
    e_highfat_ktr_avi <- 0.555; label("Proportional shift in avibactam ktr after a high-fat meal (fraction)") # Table S12: 'theta High-fat meal on Ktr-AVI' = 0.555 (RSE 30.8%)
    e_fasted_cl_avi <- -0.423; label("Proportional shift in avibactam CL/F when fasted (fraction)") # Table S12: 'theta fasted on CL/F-AVI' = -0.423 (RSE 8.63%)

    # ---- Allometry, shared by both drugs (fixed) ----
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of WT/70 on CL/F and Q/F of both drugs (unitless)") # Methods: 'fixed exponent ... 0.75 ... on clearance'; Table S12 footnote
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of WT/70 on Vc/F and Vp/F of both drugs (unitless)") # Methods: '... and 1.0 on ... volume of distribution terms'; Table S12 footnote

    # ---- Inter-individual variability, Table S12 ----
    # The diagonals are printed as CV% and are sqrt(omega^2) x 100: omega^2 =
    # (CV/100)^2 reproduces the correlations printed for the simulation model
    # (Table S14 vs Table S15) exactly, so the same convention is applied here.
    # Off-diagonals are the printed covariances. Order: CL-CTB, CL-AVI, Vc-CTB,
    # Vc-AVI, as listed in Table S12.
    etalcl + etalcl_avi + etalvc + etalvc_avi ~ c(
      0.038809,
      0.0267, 0.047524,
      0.0331, 0.0122, 0.036864,
      0.0251, 0.0501, 0.0105, 0.085849
    ) # Table S12: CL-CTB 19.7% -> 0.197^2; CL-AVI:CL-CTB 0.0267; CL-AVI 21.8% -> 0.218^2; Vc-CTB:CL-CTB 0.0331; Vc-CTB:CL-AVI 0.0122; Vc-CTB 19.2% -> 0.192^2; Vc-AVI:CL-CTB 0.0251; Vc-AVI:CL-AVI 0.0501; Vc-AVI:Vc-CTB 0.0105; Vc-AVI 29.3% -> 0.293^2
    etalktr ~ 0.120409 # Table S12: 'omega ktr-CTB:ktr-CTB (CV%)' = 34.7 -> 0.347^2
    etalktr_avi ~ 0.1444 # Table S12: 'omega ktr-AVI:ktr-AVI (CV%)' = 38 -> 0.38^2

    # ---- Residual error: additive on log-transformed concentrations ----
    expSd <- 0.463; label("Ceftibuten residual SD on the log scale") # Table S12: 'sigma CTB:CTB (CV%)' = 46.3 -> SD 0.463
    expSd_avi <- 0.42; label("Avibactam residual SD on the log scale") # Table S12: 'sigma AVI:AVI (CV%)' = 42 -> SD 0.42; the CTB-AVI residual covariance 0.117 (r = 0.60) cannot be expressed and is omitted
  })

  model({
    # Indicators relative to the reference condition (both drugs given
    # together with a regular meal; capsule or tablet).
    fasted <- 1 - FED
    mono <- 1 - CONMED_AVIBACTAMTOMILOPIL
    mono_avi <- 1 - CONMED_CEFTIBUTEN

    # Ceftibuten individual parameters. Covariate shifts are (1 + theta x
    # indicator); dose and nCrCL are power functions (Table S12 footnote;
    # the Figure S4 forest-plot ratios reproduce these forms).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q * (CRCL / 112)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    ktr <- exp(lktr + etalktr) * (1 + e_fasted_ktr * fasted) * (1 + e_highfat_ktr * FED_HIGHFAT)
    fdepot <- (DOSE_CEFTIBUTEN_MG / 400)^e_dose_fdepot *
      (1 + e_fasted_fdepot * fasted) *
      (1 + e_mono_fdepot * mono) *
      (1 + e_highfat_fdepot * FED_HIGHFAT)

    # Avibactam individual parameters.
    cl_avi <- exp(lcl_avi + etalcl_avi) * (WT / 70)^e_wt_cl_q * (CRCL / 112)^e_crcl_cl_avi *
      (1 + e_mono_cl_avi * mono_avi) * (1 + e_fasted_cl_avi * fasted)
    vc_avi <- exp(lvc_avi + etalvc_avi) * (WT / 70)^e_wt_vc_vp
    q_avi <- exp(lq_avi) * (WT / 70)^e_wt_cl_q
    vp_avi <- exp(lvp_avi) * (WT / 70)^e_wt_vc_vp
    ktr_avi <- exp(lktr_avi + etalktr_avi) *
      (1 + e_fasted_ktr_avi * fasted) *
      (1 + e_susp_ktr_avi * FORM_AVIBACTAMTOMILOPIL_SUSP) *
      (1 + e_highfat_ktr_avi * FED_HIGHFAT)
    fdepot_avi <- (DOSE_AVIBACTAMTOMILOPIL_MG / 1350)^e_dose_fdepot_avi *
      (1 + e_fasted_fdepot_avi * fasted) *
      (1 + e_susp_fdepot_avi * FORM_AVIBACTAMTOMILOPIL_SUSP)

    # Avibactam-equivalent fraction of an avibactam tomilopil dose (molecular
    # weight ratio). Not printed; back-calculated from Table S11, whose
    # fraction-excreted column divides the amount excreted by the
    # avibactam-equivalent dose: 92.0 mg / 0.455 / 300 mg = 0.674,
    # 273.5 / 0.451 / 900 = 0.674, 510.8 / 0.561 / 1350 = 0.674.
    mw_ratio_avi <- 0.674

    kel <- cl / vc
    kel_avi <- cl_avi / vc_avi
    k12_avi <- q_avi / vc_avi
    k21_avi <- q_avi / vp_avi

    # Ceftibuten: depot -> 5 transit compartments -> central, every step at ktr.
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(central) <- ktr * transit5 - kel * central

    # Avibactam: depot -> 3 transit compartments -> central <-> peripheral1.
    d/dt(depot_avi) <- -ktr_avi * depot_avi
    d/dt(transit1_avi) <- ktr_avi * depot_avi - ktr_avi * transit1_avi
    d/dt(transit2_avi) <- ktr_avi * transit1_avi - ktr_avi * transit2_avi
    d/dt(transit3_avi) <- ktr_avi * transit2_avi - ktr_avi * transit3_avi
    d/dt(central_avi) <- ktr_avi * transit3_avi - (kel_avi + k12_avi) * central_avi + k21_avi * peripheral1_avi
    d/dt(peripheral1_avi) <- k12_avi * central_avi - k21_avi * peripheral1_avi

    f(depot) <- fdepot
    f(depot_avi) <- fdepot_avi * mw_ratio_avi

    # Total plasma concentrations in mg/L (dose mg, volume L). Unbound
    # fractions used by the paper for its PK/PD targets: 0.61 (cis-ceftibuten)
    # and 0.92 (avibactam).
    Cc <- central / vc
    Cc_avi <- central_avi / vc_avi

    Cc ~ lnorm(expSd)
    Cc_avi ~ lnorm(expSd_avi)
  })
}
