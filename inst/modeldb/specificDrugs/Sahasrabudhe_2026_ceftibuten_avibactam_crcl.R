Sahasrabudhe_2026_ceftibuten_avibactam_crcl <- function() {
  description <- paste(
    "Simulation version of the joint oral population PK model for",
    "cis-ceftibuten and avibactam after oral ceftibuten and the avibactam",
    "prodrug avibactam tomilopil (Sahasrabudhe 2026), re-estimated on the",
    "healthy-adult Phase 1 data for the joint probability of target",
    "attainment simulations that selected the Phase 3 dose for complicated",
    "urinary tract infection. Same structure as the final model (ceftibuten:",
    "one compartment, depot plus five transit compartments; avibactam: two",
    "compartments, depot plus three transit compartments; fixed allometric",
    "weight exponents), but the power effect of BSA-normalised creatinine",
    "clearance is replaced by a linear effect of absolute Cockcroft-Gault",
    "creatinine clearance on both clearances, with slopes fixed from published",
    "renal-impairment studies, and the covariate set is pared down: the",
    "ceftibuten food and single-agent effects on bioavailability and the",
    "avibactam tomilopil dose effect are removed, the ceftibuten dose effect",
    "becomes a step for doses above 800 mg, and a high-fat-meal effect on",
    "avibactam bioavailability is added. The avibactam tomilopil dose is",
    "entered in mg of prodrug and converted to avibactam equivalents inside",
    "the model. Ceftibuten is the unsuffixed parent; avibactam carries the",
    "sibling-drug suffix _avi throughout. See",
    "modellib('Sahasrabudhe_2026_ceftibuten_avibactam') for the final model."
  )
  reference <- paste(
    "Sahasrabudhe SA, Shoji S, Jakate AS, Abraham BK, Toussi SS, Maligalig BJ,",
    "Darekar AS, Tripathy S, Kantaridis C, Mfopou JK, Selivanova O, Bergman A,",
    "Alami NN, Singh RSP. Dose Selection for Ceftibuten Combined with",
    "Avibactam Tomilopil, a Potential Oral Alternative to Carbapenems for",
    "Patients with Complicated Urinary Tract Infections. J Clin Pharmacol.",
    "2026;66(9):e70263. doi:10.1002/jcph.70263.",
    "All parameter estimates are from Supplemental Appendix Table S15",
    "'Parameter Estimates for CTB/AVI-Oral Joint Pharmacokinetic Model for",
    "Simulations (Simulation Model)'; the changes from the final model are",
    "listed in Table S13 and Table S1 (step S6), and the covariate forms are",
    "checked against the Figure S4 (bottom) forest plot.",
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
        "Vp (fixed exponent 1) of both drugs, reference 70 kg (Table S15",
        "footnote). Kept alongside the absolute CrCL term, as in the paper."
      ),
      source_name = "WT"
    ),
    CRCL = list(
      description = paste(
        "ABSOLUTE (not BSA-normalised) Cockcroft-Gault creatinine clearance",
        "in mL/min"
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Raw mL/min, NOT mL/min/1.73 m^2: the simulation model replaces the",
        "final model's BSA-normalised power term with 'a linear effect of",
        "CrCL on CL/F' based on 'absolute CrCL' (Results; Discussion), so",
        "supplying a BSA-normalised value would mis-scale both clearances.",
        "Enters as 1 + slope x (CrCL - 112) on ceftibuten and avibactam CL/F",
        "(Table S15 footnote), with slopes 0.00794 and 0.00931 per mL/min",
        "fixed from published renal-impairment studies. The linear form makes",
        "avibactam CL/F non-positive below about 4.6 mL/min, outside the",
        "range the published slopes support."
      ),
      source_name = "CrCL"
    ),
    FED = list(
      description = "Fed state at dosing (1 = fed, 0 = fasted)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (fed; the model reference is co-administration with a regular meal)",
      notes = paste(
        "The paper's 'fasted' effects enter as (1 - FED): on the ceftibuten",
        "and avibactam transit rate constants, avibactam relative",
        "bioavailability and avibactam CL/F (Table S15). The ceftibuten",
        "fasted effect on bioavailability was removed for simulation (Table",
        "S13). The jPTA simulations were run separately with all patients",
        "fed or all fasted."
      ),
      source_name = "fasted"
    ),
    FED_HIGHFAT = list(
      description = "Dosed after a high-fat, high-calorie meal (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (regular meal or fasted)",
      notes = paste(
        "High-fat-meal records also have FED = 1. In the simulation model the",
        "only high-fat-meal effect is a 19.2% increase in avibactam relative",
        "bioavailability (Table S13, Table S15); the ceftibuten and avibactam",
        "transit-rate effects were removed."
      ),
      source_name = "high-fat meal"
    ),
    FORM_AVIBACTAMTOMILOPIL_SUSP = list(
      description = "Avibactam tomilopil given as the oral suspension (1) rather than a capsule or tablet (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (capsule or fixed-dose combination tablet)",
      notes = paste(
        "Raises the avibactam transit rate constant and relative",
        "bioavailability (Table S15). Used only in Study C4691002, where",
        "avibactam tomilopil was always given alone; set to 0 for capsule or",
        "tablet simulations."
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
        "(1 - CONMED_CEFTIBUTEN) (Table S15). Set to 1 for subjects who",
        "receive both drugs."
      ),
      source_name = "AVP alone"
    ),
    DOSE_CEFTIBUTEN_MG = list(
      description = "Ceftibuten dose of the current administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters only through the step indicator DOSE_CEFTIBUTEN_MG > 800:",
        "ceftibuten relative bioavailability is lower by 12.6% for doses",
        "above 800 mg and unchanged at or below 800 mg (Table S1 footnote e;",
        "Table S13; Table S15 footnote)."
      ),
      source_name = "CTB dose"
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
    disease_state = paste(
      "Healthy adult volunteers (estimation); intended for simulation of",
      "adults with complicated urinary tract infection across renal function"
    ),
    dose_range = paste(
      "Ceftibuten 400-1600 mg single doses and 400-800 mg q8h; avibactam",
      "tomilopil 300-1350 mg single doses and 900-1350 mg q8h; simulated",
      "regimens ceftibuten 400 mg + avibactam tomilopil 800-1350 mg q8h"
    ),
    regions = "Australia, Belgium",
    renal_function = "CrCL 58.5-179 mL/min (median 124) in the estimation data set",
    notes = paste(
      "Same 71 healthy adults and three Phase 1 studies as the final model",
      "(Table 2). The CrCL slopes were fixed from published studies that",
      "included mild, moderate and severe renal impairment, so the model is",
      "meant for simulation beyond the healthy-adult CrCL range (Discussion)."
    )
  )

  ini({
    # ---- Ceftibuten (unsuffixed parent), Table S15 ----
    lcl <- log(4.4); label("Ceftibuten apparent clearance CL/F at WT 70 kg, CrCL 112 mL/min, fed (L/h)") # Table S15: 'theta CL/F-CTB (L/h)' = 4.4 (RSE 4.31%)
    lvc <- log(19.4); label("Ceftibuten apparent central volume Vc/F at WT 70 kg (L)") # Table S15: 'theta Vc/F-CTB (L)' = 19.4 (RSE 3.78%)
    lktr <- log(3.02); label("Ceftibuten transit rate constant, fed (1/h)") # Table S15: 'theta Ktr-CTB (1/h)' = 3.02 (RSE 6.91%)

    e_crcl_cl <- fixed(0.00794); label("Linear slope of (CrCL - 112) on ceftibuten CL/F (per mL/min)") # Table S15: 'theta CrCL on CL/F-CTB' = 0.00794, RSE 'Fixed'
    e_fasted_ktr <- 0.65; label("Proportional shift in ceftibuten ktr when fasted (fraction)") # Table S15: 'theta fasted on Ktr-CTB' = 0.65 (RSE 22.2%)
    e_dose_fdepot <- -0.126; label("Proportional shift in ceftibuten relative bioavailability for doses above 800 mg (fraction)") # Table S15: 'theta CTB DOSE on rBA-CTB' = -0.126 (RSE 54.7%)

    # ---- Avibactam (suffix _avi), Table S15 ----
    lcl_avi <- log(20.8); label("Avibactam apparent clearance CL/F at WT 70 kg, CrCL 112 mL/min, fed, with ceftibuten (L/h)") # Table S15: 'theta CL/F-AVI (L/h)' = 20.8 (RSE 4.52%)
    lvc_avi <- log(76.6); label("Avibactam apparent central volume Vc/F at WT 70 kg (L)") # Table S15: 'theta Vc/F-AVI (L)' = 76.6 (RSE 6.75%)
    lq_avi <- log(0.356); label("Avibactam apparent inter-compartmental clearance Q/F at WT 70 kg (L/h)") # Table S15: 'theta Q/F-AVI (L/h)' = 0.356 (RSE 10.3%)
    lvp_avi <- log(12.4); label("Avibactam apparent peripheral volume Vp/F at WT 70 kg (L)") # Table S15: 'theta Vp/F-AVI (L)' = 12.4 (RSE 30%)
    lktr_avi <- log(1.83); label("Avibactam transit rate constant, capsule or tablet, fed (1/h)") # Table S15: 'theta Ktr-AVI (1/h)' = 1.83 (RSE 6.41%)

    e_crcl_cl_avi <- fixed(0.00931); label("Linear slope of (CrCL - 112) on avibactam CL/F (per mL/min)") # Table S15: 'theta CrCL on CL/F-AVI' = 0.00931, RSE 'Fixed'
    e_fasted_ktr_avi <- 1.51; label("Proportional shift in avibactam ktr when fasted (fraction)") # Table S15: 'theta fasted on Ktr-AVI' = 1.51 (RSE 16.3%)
    e_fasted_fdepot_avi <- -0.523; label("Proportional shift in avibactam relative bioavailability when fasted (fraction)") # Table S15: 'theta fasted on rBA-AVI' = -0.523 (RSE 7.55%)
    e_susp_ktr_avi <- 2.02; label("Proportional shift in avibactam ktr for the avibactam tomilopil suspension (fraction)") # Table S15: 'theta SUSP on Ktr-AVI' = 2.02 (RSE 16.4%)
    e_susp_fdepot_avi <- 0.317; label("Proportional shift in avibactam relative bioavailability for the suspension (fraction)") # Table S15: 'theta SUSP on rBA-AVI' = 0.317 (RSE 40.4%)
    e_mono_cl_avi <- 0.237; label("Proportional shift in avibactam CL/F when avibactam tomilopil is given alone (fraction)") # Table S15: 'theta AVP alone on CL/F-AVI' = 0.237 (RSE 40.9%)
    e_highfat_fdepot_avi <- 0.192; label("Proportional shift in avibactam relative bioavailability after a high-fat meal (fraction)") # Table S15: 'theta high-fat meal on rBA-AVI' = 0.192 (RSE 33.2%)
    e_fasted_cl_avi <- -0.439; label("Proportional shift in avibactam CL/F when fasted (fraction)") # Table S15: 'theta fasted on CL/F-AVI' = -0.439 (RSE 8%)

    # ---- Allometry, shared by both drugs (fixed) ----
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of WT/70 on CL/F and Q/F of both drugs (unitless)") # Table S15 footnote: '(WT/70)^0.75 for CL/F and Q/F'
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of WT/70 on Vc/F and Vp/F of both drugs (unitless)") # Table S15 footnote: '(WT/70)^1.00 for Vc/F and Vp/F'

    # ---- Inter-individual variability, Table S15 ----
    # Diagonals printed as CV% = sqrt(omega^2) x 100; off-diagonals are the
    # printed covariances. omega^2 = (CV/100)^2 reproduces the Table S14
    # correlations (e.g. 0.0582 / (0.270 x 0.274) = 0.787 vs 0.788).
    etalcl + etalcl_avi + etalvc + etalvc_avi ~ c(
      0.0729,
      0.0582, 0.075076,
      0.0436, 0.0174, 0.046225,
      0.0127, 0.0332, 0.00419, 0.068644
    ) # Table S15: CL-CTB 27% -> 0.27^2; CL-AVI:CL-CTB 0.0582; CL-AVI 27.4% -> 0.274^2; Vc-CTB:CL-CTB 0.0436; Vc-CTB:CL-AVI 0.0174; Vc-CTB 21.5% -> 0.215^2; Vc-AVI:CL-CTB 0.0127; Vc-AVI:CL-AVI 0.0332; Vc-AVI:Vc-CTB 0.00419; Vc-AVI 26.2% -> 0.262^2
    etalktr ~ 0.174724 # Table S15: 'omega ktr-CTB:ktr-CTB (CV%)' = 41.8 -> 0.418^2
    etalktr_avi ~ 0.142884 # Table S15: 'omega ktr-AVI:ktr-AVI (CV%)' = 37.8 -> 0.378^2

    # ---- Residual error: additive on log-transformed concentrations ----
    expSd <- 0.484; label("Ceftibuten residual SD on the log scale") # Table S15: 'sigma CTB:CTB (CV%)' = 48.4 -> SD 0.484
    expSd_avi <- 0.422; label("Avibactam residual SD on the log scale") # Table S15: 'sigma AVI:AVI (CV%)' = 42.2 -> SD 0.422; the CTB-AVI residual covariance 0.123 (Table S14 r = 0.601) cannot be expressed and is omitted
  })

  model({
    # Indicators relative to the reference condition (both drugs given
    # together, fed; capsule or tablet).
    fasted <- 1 - FED
    mono_avi <- 1 - CONMED_CEFTIBUTEN
    dose_gt800 <- DOSE_CEFTIBUTEN_MG > 800

    # Renal function: 1 + slope x (CrCL - 112) (Table S15 footnote; the
    # Figure S4 bottom forest plot gives 1.30 / 0.59 for ceftibuten and
    # 1.35 / 0.52 for avibactam at CrCL 150 / 60 mL/min, matching this form).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q * (1 + e_crcl_cl * (CRCL - 112))
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    ktr <- exp(lktr + etalktr) * (1 + e_fasted_ktr * fasted)
    fdepot <- 1 + e_dose_fdepot * dose_gt800

    cl_avi <- exp(lcl_avi + etalcl_avi) * (WT / 70)^e_wt_cl_q * (1 + e_crcl_cl_avi * (CRCL - 112)) *
      (1 + e_mono_cl_avi * mono_avi) * (1 + e_fasted_cl_avi * fasted)
    vc_avi <- exp(lvc_avi + etalvc_avi) * (WT / 70)^e_wt_vc_vp
    q_avi <- exp(lq_avi) * (WT / 70)^e_wt_cl_q
    vp_avi <- exp(lvp_avi) * (WT / 70)^e_wt_vc_vp
    ktr_avi <- exp(lktr_avi + etalktr_avi) *
      (1 + e_fasted_ktr_avi * fasted) *
      (1 + e_susp_ktr_avi * FORM_AVIBACTAMTOMILOPIL_SUSP)
    fdepot_avi <- (1 + e_fasted_fdepot_avi * fasted) *
      (1 + e_susp_fdepot_avi * FORM_AVIBACTAMTOMILOPIL_SUSP) *
      (1 + e_highfat_fdepot_avi * FED_HIGHFAT)

    # Avibactam-equivalent fraction of an avibactam tomilopil dose (molecular
    # weight ratio). Not printed; back-calculated from Table S11 (amount
    # excreted / fraction excreted / prodrug dose = 0.674 at 300, 900 and
    # 1350 mg).
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

    # Total plasma concentrations in mg/L. The jPTA targets are on unbound
    # concentrations: unbound fractions 0.61 (cis-ceftibuten) and 0.92
    # (avibactam) (Methods).
    Cc <- central / vc
    Cc_avi <- central_avi / vc_avi

    Cc ~ lnorm(expSd)
    Cc_avi ~ lnorm(expSd_avi)
  })
}
