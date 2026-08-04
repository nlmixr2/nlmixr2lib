Mitra_2026_ziftomenib <- function() {
  description <- "Sequential two-stage population PK model for oral ziftomenib (a potent, selective, oral menin inhibitor for R/R NPM1-mutated acute myeloid leukemia) and its two active metabolites KO-739 and KO-516 (Mitra 2026 Kura Oncology KOMET-001 + KO-MEN-003). Parent PK is a 2-compartment model with first-order absorption, absorption lag time, and linear elimination from the central compartment; oral bioavailability F1 is fixed at 0.129 (identifiability constraint from the human ADME + absolute-BA study KO-MEN-005). Each metabolite is 2-compartment with linear elimination; the metabolic clearance is split between KO-739 and KO-516 by a fixed 1:1 in-vitro-anchored biotransformation ratio (FM_KO516 = 0.5), with the total metabolized fraction FM held fixed at 0.535 after an initial identifiability-limited estimation. Covariate effects retained in the final model: FED and PPI on parent F1 (logit-scale shifts +3.21 fed; -0.520 PPI = 6.09x and 0.627x multipliers on F1), PPI on parent Ka (log-scale shift -0.485 = 0.616x), FED on parent absorption lag time (log-scale shift +0.322 = 1.38x), strong CYP3A4 inhibitor on parent CL/F (log-scale -0.778 = 0.459x), healthy-volunteer status on parent CL/F (log-scale +0.950 = 2.59x), healthy-volunteer status on FM (logit-scale -1.62 = 0.348x multiplier on FM), strong CYP3A4 inhibitor on KO-739 CL (log-scale -1.64 = 0.195x), strong CYP3A4 inhibitor on KO-516 CL (log-scale -0.802 = 0.449x), healthy-volunteer status on KO-739 Vc (log-scale -1.62 = 0.197x), and healthy-volunteer status on KO-516 Vc (log-scale -1.87 = 0.154x). No effect of NPM1-m vs KMT2A-r mutational status, body weight, sex, race, age, mild/moderate renal or hepatic impairment, or P-gp inhibitor coadministration on ziftomenib PK. IIV: parent 47.3% CV on CL and 120% CV on Vc; metabolites 74.7% (KO-739 CL), 110% (KO-739 Vc), 162% (KO-739 Q), 31.2% (KO-516 CL), 191% (KO-516 Vc), 118% (KO-516 Q), and 56.8% CV on FM (all independent diagonals). Inter-occasion variability on F1 (Omega 1.06 corresponding to 137.3% CV) reported in the parent NONMEM run across 3 occasions is not encoded structurally here (no operational occasion column is defined for the model-library use case; see vignette Assumptions and deviations). Residual error: proportional 43.7% CV on parent Cc; proportional 45.2% CV plus additive 0.128 ng/mL on Cc_ko739; proportional 36.4% CV on Cc_ko516."
  reference   <- "Mitra A, Yang X, Ortiz RH, Jomphe C, Leoni M, Gosselin NH. Population Pharmacokinetics and Exposure-Response Analysis of Ziftomenib in Relapsed or Refractory Acute Myeloid Leukemia Patients With NPM1 Mutation. CPT Pharmacometrics Syst Pharmacol. 2026. doi:10.1002/psp4.70244."
  vignette    <- "Mitra_2026_ziftomenib"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot             = list(analyte = "ziftomenib", units = "mg", specimen = "administration site", verified = FALSE),
    central           = list(analyte = "ziftomenib", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1       = list(analyte = "ziftomenib", units = "mg", specimen = "plasma", verified = FALSE),
    central_ko739     = list(analyte = "KO-739", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1_ko739 = list(analyte = "KO-739", units = "mg", specimen = "plasma", verified = FALSE),
    central_ko516     = list(analyte = "KO-516", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1_ko516 = list(analyte = "KO-516", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator, 1 = healthy subject (KO-MEN-003 crossover food-effect / PPI-effect study), 0 = patient with R/R NPM1-m or KMT2A-r AML (KOMET-001).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (R/R AML patient; the KOMET-001 reference cohort)",
      notes              = "Time-fixed per subject. Multiplicative effects (via the paper's log/logit-scale shifts): 2.59x on parent CL/F when DIS_HEALTHY = 1 (Table 1 'Healthy volunteers on CL'; corresponds to log-scale THETA(12) = 0.950); 0.348x on FM (Table 1 'Healthy volunteers on FM'; logit-scale THETA(17) = -1.62 -> multiplicative reduction of parent-to-metabolite fraction by 65.2% in healthy subjects); 0.197x on KO-739 Vc (Table 1 'Healthy volunteers on Vc of KO-739'; log-scale THETA(21) = -1.62); 0.154x on KO-516 Vc (Table 1 'Healthy volunteers on Vc of KO-516'; log-scale THETA(22) = -1.87). Source paper column 'HV' (0/1). Paper notes the elevated ziftomenib exposure in AML patients relative to healthy volunteers (GMR ~2.55 for AUCss, ~2.18 for Cmax,ss) is attributed to disease-state effects on CL/F and (for metabolite Vc) distribution.",
      source_name        = "HV"
    ),
    CONMED_CYP3A4_INH_STRONG = list(
      description        = "Concomitant strong CYP3A4 inhibitor coadministration indicator (per-dose-record), 1 = ziftomenib dose administered during a period of concomitant strong CYP3A4 inhibitor (e.g., posaconazole, voriconazole, itraconazole, ketoconazole) coadministration, 0 = no strong CYP3A4 inhibitor coadministration during the dose period.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no strong CYP3A4 inhibitor coadministration)",
      notes              = "Time-varying per record. Multiplicative effects (via the paper's log-scale shifts): 0.459x on parent CL/F when CONMED_CYP3A4_INH_STRONG = 1 (Table 1 'Strong CYP3A4 inhibitors on CL'; log-scale THETA(11) = -0.778, corresponding to a 54.1% reduction in parent CL/F, consistent with ziftomenib being primarily metabolized by CYP3A4); 0.195x on KO-739 CL (Table 1 'Strong CYP3A4 inhibitors on CL of KO-739'; log-scale THETA(18) = -1.64, 80.5% reduction); 0.449x on KO-516 CL (Table 1 'Strong CYP3A4 inhibitors on CL of KO-516'; log-scale THETA(19) = -0.802, 55.1% reduction). Source paper column 'DDICYPSH'. Weak- and moderate-CYP3A4-inhibitor effects on parent CL were tested (paper's THETA(10)) but fixed at 0 in the final model; the CONMED_CYP3A4_INH_STRONG canonical captures only the strong-inhibitor stratum. Concomitant antifungal azoles are the dominant driver of this indicator in the KOMET-001 AML cohort (used for prophylaxis against invasive fungal infections).",
      source_name        = "DDICYPSH"
    ),
    CONMED_PPI = list(
      description        = "Concomitant proton-pump inhibitor coadministration indicator (per-dose-record), 1 = ziftomenib dose administered during a period of concomitant PPI (e.g., rabeprazole in study KO-MEN-003; or an oncology-standard PPI in KOMET-001) coadministration, 0 = no PPI coadministration during the dose period.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no PPI coadministration)",
      notes              = "Time-varying per record. Multiplicative effects (via the paper's log/logit-scale shifts): 0.627x on parent F1 when CONMED_PPI = 1 (Table 1 'Effect of PPI on F1'; logit-scale THETA(15) = -0.520, 37.3% reduction in oral bioavailability); 0.616x on parent Ka (Table 1 'Effect of PPI on Ka'; log-scale THETA(13) = -0.485, 38.4% slower first-order absorption). The dual effect on F1 and Ka is consistent with ziftomenib's pH-solubility profile of reduced solubility at elevated gastric pH. Source paper column 'DDIPPI'.",
      source_name        = "DDIPPI"
    ),
    FED = list(
      description        = "Fed-vs-fasted dose-record indicator (per-dose-record), 1 = ziftomenib dose administered in the fed state (any meal), 0 = ziftomenib dose administered fasted (recommended clinical administration state).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "Time-varying per record. Multiplicative effects (via the paper's log/logit-scale shifts): 6.09x on parent F1 when FED = 1 (Table 1 'Effect of FED on F1'; logit-scale THETA(14) = 3.21, corresponding to a 6.1-fold increase in F1 due to food, consistent with ziftomenib's low solubility); 1.38x on parent absorption lag time (Table 1 'Effect of FED on Lag'; log-scale THETA(16) = 0.322). Source paper column 'FED'. Paper recommends fasted administration in clinical use (Ziftomenib label recommendation; Discussion).",
      source_name        = "FED"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 188L,
    n_studies       = 2L,
    n_patients      = 174L,
    n_healthy       = 14L,
    n_observations_parent    = 2436L,
    n_observations_ko739     = 2376L,
    n_observations_ko516     = 2299L,
    age_range       = "18-86 years (overall median 64.5); patients 18-86 (median 66.0); healthy volunteers 22-55 (median 39.0)",
    weight_range    = "41.0-135 kg (overall median 71.5); patients 41.0-135 (median 71.1); healthy volunteers 61.5-93.0 (median 76.7)",
    sex_female_pct  = 51.1,
    race_ethnicity  = c(White = 75.0, `Black or African American` = 5.3, Asian = 3.2, Other = 2.1, Missing = 14.4),
    disease_state   = "Pooled cohort of 174 adult patients with relapsed or refractory acute myeloid leukemia (R/R AML; predominantly NPM1-m in Phase 1b/2 and KMT2A-r in Phase 1a of KOMET-001) enrolled in KOMET-001 (NCT04067336; a Phase 1/2 first-in-human study of oral ziftomenib) plus 14 healthy volunteers enrolled in KO-MEN-003 (a Phase 1, open-label, crossover food-effect / rabeprazole DDI study of ziftomenib). At baseline in KOMET-001: 137/174 (78.7%) with normal hepatic function, 32/174 (18.4%) mild hepatic impairment, 5/174 (2.9%) moderate hepatic impairment; 73/174 (42.0%) normal renal function, 66/174 (37.9%) mild, 21/174 (12.1%) moderate renal impairment. No patients with severe renal or hepatic impairment were enrolled.",
    dose_range      = "KOMET-001 Phase 1a: 50 mg QD dose escalation up to 1000 mg QD (7 cohorts); Phase 1b: single-agent ziftomenib PO QD; Phase 2 (registrational-enabling): 600 mg PO QD. KO-MEN-003: single 400 mg oral dose in each period (fed or fasted) with or without a 20 mg PO rabeprazole pretreatment.",
    regions         = "International (KOMET-001 conducted at 21 sites)",
    assay           = "Validated bioanalytical assay for ziftomenib, KO-739, and KO-516 (LLOQ = 0.2 ng/mL for each analyte per Supplementary Table S3 footnote).",
    iov_structure   = "Paper reports inter-occasion variability (IOV) on ziftomenib F1 across 3 occasions in the final parent NONMEM run (Table 1 IOV Omega = 1.06, corresponding to 137.3% CV; Occasion #1 and #2 = rich sample profiles in each study, Occasion #3 = sparse samples in KOMET-001 or third rich profile in KO-MEN-003). This model file does NOT encode the IOV structurally -- no operational occasion column is defined for the model-library use case, and the nlmixr2lib convention (Yin_2020_pexidartinib precedent) is to omit IOV when no occasion mapping is defined; see vignette Assumptions and deviations.",
    notes           = "Sequential two-stage popPK: parent (ziftomenib) fit first, then individual empirical Bayes estimates (EBEs) of parent PK parameters were fixed and used as inputs to the metabolite (KO-739 + KO-516) model. Final metabolite NONMEM run also has fixed F1 = 0.129, FM = 0.535, FM_ko516 = 0.5 (see reference Discussion). Reference subject for typical-value simulations: R/R AML patient, fasted, no PPI, no strong CYP3A4 inhibitor coadministration (Table 1 base). Parent Phase 1 NONMEM run: OBJV = 23621.837 (2436 obs, 188 individuals). Metabolite Phase 4 NONMEM run: OBJV = 28349.343 (metabolites only after fixing parent EBEs)."
  )

  ini({
    # ---------------- Parent (ziftomenib) structural fixed effects ----------------
    # Paper Table 1 reports linear-scale values with 95% CIs. The Supplement's
    # NONMEM control stream stores THETAs in the log domain and (for F1) on the
    # logit scale; the log-scale THETA values quoted below in each source-trace
    # comment are copied verbatim from the Supplement's Final Parameter Estimate
    # block ("Final Population PK Model of Ziftomenib - NONMEM Output Code").
    lka       <- log(0.0928); label("First-order absorption rate KA (1/h)")                                            # Mitra 2026 Table 1 KA = 0.0928 1/h; Supp NONMEM TH2 = -2.38 -> exp(-2.38) = 0.0927
    lcl       <- log(11.6);   label("Apparent parent clearance CL/F (L/h) in R/R AML patients, no strong CYP3A4 inhibitor")    # Table 1 CL = 11.6 L/h; Supp TH4 = 2.45 -> exp(2.45) = 11.6
    lvc       <- log(54.6);   label("Apparent parent central volume Vc/F (L)")                                          # Table 1 Vc = 54.6 L; Supp TH3 = 4.00 -> exp(4.00) = 54.6
    lq        <- log(27.7);   label("Apparent parent inter-compartmental clearance Q/F (L/h)")                         # Table 1 Q = 27.7 L/h; Supp TH5 = 3.32 -> exp(3.32) = 27.7
    lvp       <- log(1106);   label("Apparent parent peripheral volume Vp/F (L)")                                       # Table 1 Vp = 1106 L; Supp TH6 = 7.01 -> exp(7.01) = 1108
    ltlag     <- log(0.325);  label("Parent absorption lag time ALAG1 (hr)")                                            # Table 1 Lag = 0.325 h; Supp TH7 = 0.325 (linear scale in this THETA position)

    # F1 is on logit scale in the source NONMEM code (THETA(1) FIXED at -1.91,
    # yielding logit^-1(-1.91) = 0.129 = 12.9%, matching the KO-MEN-005 absolute
    # bioavailability point estimate of 12.9% used to break the CL vs F
    # identifiability constraint). Covariate effects on F1 (FED, PPI) enter
    # additively on the logit scale.
    logitfdepot <- fixed(-1.91); label("Logit-transformed base fraction of oral absorption F1 (unitless); logit^-1(1.91) = 0.129 = paper's F1")  # Supp NONMEM TH1 = -1.91 FIX; Table 1 F1 = 0.129 Fixed

    # ---------------- Covariate effects on parent parameters ----------------
    e_fed_logitfdepot          <- 3.21    ; label("Additive shift on logit(F1) for fed vs fasted state (unitless; +3.21 -> F1 x 6.09)")                # Supp TH14 = 3.21 -> logit^-1(-1.91 + 3.21) / logit^-1(-1.91) = 0.786 / 0.129 = 6.09; Table 1 'Effect of FED on F1' = x6.09
    e_conmed_ppi_logitfdepot   <- -0.520  ; label("Additive shift on logit(F1) for concomitant PPI (unitless; -0.520 -> F1 x 0.627)")                  # Supp TH15 = -0.520 -> logit^-1(-1.91 - 0.520) / logit^-1(-1.91) = 0.0809 / 0.129 = 0.627; Table 1 'Effect of PPI on F1' = x0.627
    e_conmed_ppi_ka            <- -0.485  ; label("Additive shift on log(KA) for concomitant PPI (unitless; -0.485 -> KA x 0.616)")                    # Supp TH13 = -0.485 -> exp(-0.485) = 0.616; Table 1 'Effect of PPI on Ka' = x0.616
    e_fed_ltlag                <- 0.322   ; label("Additive shift on log(ALAG1) for fed vs fasted state (unitless; +0.322 -> ALAG1 x 1.38)")           # Supp TH16 = 0.322 -> exp(0.322) = 1.38; Table 1 'Effect of FED on Lag' = x1.38
    e_dis_healthy_cl           <- 0.950   ; label("Additive shift on log(CL) for healthy-volunteer cohort (unitless; +0.950 -> CL x 2.59)")            # Supp TH12 = 0.950 -> exp(0.950) = 2.586; Table 1 'Healthy volunteers on CL' = x2.59
    e_conmed_cyp3a4_inh_strong_cl <- -0.778 ; label("Additive shift on log(CL) for concomitant strong CYP3A4 inhibitor (unitless; -0.778 -> CL x 0.459)")  # Supp TH11 = -0.778 -> exp(-0.778) = 0.459; Table 1 'Strong CYP3A4 inhibitors on CL' = x0.459

    # ---------------- Parent IIV ----------------
    # Only CL and Vc carry IIV in the paper's final parent model (Supp OMEGA
    # block; all other diagonals are FIX 0). Variances taken from the OMEGA
    # matrix's final estimates. %CV verified: sqrt(exp(omega^2) - 1) * 100.
    #   OMEGA(3,3) = 0.893 (V2) -> sqrt(exp(0.893) - 1) = 1.20 = 120% CV
    #   OMEGA(4,4) = 0.202 (CL) -> sqrt(exp(0.202) - 1) = 0.473 = 47.3% CV
    etalcl ~ 0.202
    etalvc ~ 0.893

    # ---------------- Parent residual error ----------------
    # Table 1 residual error = proportional 43.7% CV = 0.437 fraction.
    # Supp NONMEM ERRPROP TH8 = 0.437 (linear-scale proportional SD).
    propSd <- 0.437; label("Parent proportional residual SD (fraction)")

    # ---------------- Metabolite structural fixed effects ----------------
    # Metabolite Phase 4 NONMEM control stream (Supp "Final Population PK
    # Model of Metabolites - NONMEM Output Code") stores THETAs 1-8 in the
    # log domain; the linear-scale metabolite values in Table 1 back out to
    # those log-scale THETAs to within rounding.
    lq_ko739  <- log(4.13);  label("Apparent KO-739 inter-compartmental clearance Q_KO-739 (L/h)")           # Table 1 Q_KO-739 = 4.13 L/h; Supp meta TH1 = 1.42 -> exp(1.42) = 4.14
    lvc_ko739 <- log(8.20);  label("Apparent KO-739 central volume Vc_KO-739 (L) in R/R AML patients")        # Table 1 Vc of KO-739 = 8.20 L; Supp meta TH2 = 2.10 -> exp(2.10) = 8.17
    lvp_ko739 <- log(240);   label("Apparent KO-739 peripheral volume Vp_KO-739 (L)")                          # Table 1 Vp of KO-739 = 240 L; Supp meta TH3 = 5.48 -> exp(5.48) = 240
    lcl_ko739 <- log(8.50);  label("Apparent KO-739 clearance CL_KO-739 (L/h) in R/R AML patients, no strong CYP3A4 inhibitor")  # Table 1 CL of KO-739 = 8.50 L/h; Supp meta TH4 = 2.14 -> exp(2.14) = 8.50
    lq_ko516  <- log(9.55);  label("Apparent KO-516 inter-compartmental clearance Q_KO-516 (L/h)")           # Table 1 Q of KO-516 = 9.55 L/h; Supp meta TH5 = 2.26 -> exp(2.26) = 9.58
    lvc_ko516 <- log(11.8);  label("Apparent KO-516 central volume Vc_KO-516 (L) in R/R AML patients")        # Table 1 Vc of KO-516 = 11.8 L; Supp meta TH6 = 2.47 -> exp(2.47) = 11.8
    lvp_ko516 <- log(604);   label("Apparent KO-516 peripheral volume Vp_KO-516 (L)")                          # Table 1 Vp of KO-516 = 604 L; Supp meta TH7 = 6.40 -> exp(6.40) = 601
    lcl_ko516 <- log(21.7);  label("Apparent KO-516 clearance CL_KO-516 (L/h) in R/R AML patients, no strong CYP3A4 inhibitor")  # Table 1 CL of KO-516 = 21.7 L/h; Supp meta TH8 = 3.08 -> exp(3.08) = 21.8

    # Fraction of parent ziftomenib metabolized (FM) is on logit scale in the
    # source NONMEM code (THETA(9) FIXED at 0.14, yielding logit^-1(0.14) =
    # 0.535). The parent-to-metabolite split FM_KO516 = 0.5 (1:1 KO-739 :
    # KO-516) is fixed based on in-vitro human-liver-microsome relative
    # abundance (Kura Oncology Inc., data on file per Table 1 footnote d)
    # and is hardcoded at 0.5 inside model() below since it has no covariate
    # effect and no IIV.
    logitfm <- fixed(0.14); label("Logit-transformed base fraction of parent metabolized to KO-739 + KO-516 (unitless); logit^-1(0.14) = 0.535")  # Supp meta TH9 = 0.14 FIX; Table 1 FM = 0.535 Fixed

    # ---------------- Covariate effects on metabolite parameters ----------------
    e_dis_healthy_logitfm            <- -1.62 ; label("Additive shift on logit(FM) for healthy-volunteer cohort (unitless; -1.62 -> FM x 0.348)")        # Supp meta TH17 = -1.62 -> logit^-1(0.14 - 1.62) / logit^-1(0.14) = 0.186 / 0.535 = 0.348; Table 1 'Healthy volunteers on FM' = x0.348
    e_conmed_cyp3a4_inh_strong_cl_ko739 <- -1.64 ; label("Additive shift on log(CL_KO-739) for concomitant strong CYP3A4 inhibitor (unitless; -1.64 -> CL_KO-739 x 0.195)")  # Supp meta TH18 = -1.64 -> exp(-1.64) = 0.194; Table 1 'Strong CYP3A4 inhibitors on CL of KO-739' = x0.195
    e_conmed_cyp3a4_inh_strong_cl_ko516 <- -0.802; label("Additive shift on log(CL_KO-516) for concomitant strong CYP3A4 inhibitor (unitless; -0.802 -> CL_KO-516 x 0.449)")  # Supp meta TH19 = -0.802 -> exp(-0.802) = 0.448; Table 1 'Strong CYP3A4 inhibitors on CL of KO-516' = x0.449
    e_dis_healthy_vc_ko739           <- -1.62 ; label("Additive shift on log(Vc_KO-739) for healthy-volunteer cohort (unitless; -1.62 -> Vc_KO-739 x 0.197)")  # Supp meta TH21 = -1.62 -> exp(-1.62) = 0.198; Table 1 'Healthy volunteers on Vc of KO-739' = x0.197
    e_dis_healthy_vc_ko516           <- -1.87 ; label("Additive shift on log(Vc_KO-516) for healthy-volunteer cohort (unitless; -1.87 -> Vc_KO-516 x 0.154)")  # Supp meta TH22 = -1.87 -> exp(-1.87) = 0.154; Table 1 'Healthy volunteers on Vc of KO-516' = x0.154

    # ---------------- Metabolite IIV ----------------
    # OMEGA final estimates from the metabolite NONMEM output. All diagonals
    # (V3 for both metabolites and FM_ko516 are FIX 0). %CV verified as
    # sqrt(exp(omega^2) - 1) * 100.
    #   OMEGA(1,1) = 1.29 (Q_KO-739)  -> 162% CV
    #   OMEGA(2,2) = 0.79 (V2_KO-739) -> 110% CV
    #   OMEGA(4,4) = 0.443 (CL_KO-739) -> 74.7% CV
    #   OMEGA(5,5) = 0.874 (Q_KO-516) -> 118% CV
    #   OMEGA(6,6) = 1.54 (V2_KO-516) -> 191% CV
    #   OMEGA(8,8) = 0.0928 (CL_KO-516) -> 31.2% CV
    #   OMEGA(9,9) = 0.280 (FM) -> 56.8% CV
    etalq_ko739  ~ 1.29
    etalvc_ko739 ~ 0.79
    etalcl_ko739 ~ 0.443
    etalq_ko516  ~ 0.874
    etalvc_ko516 ~ 1.54
    etalcl_ko516 ~ 0.0928
    etalogitfm   ~ 0.280

    # ---------------- Metabolite residual error ----------------
    # KO-739 uses a combined additive + proportional structure; KO-516 is
    # proportional-only (additive term fixed at 0 in the final metabolite
    # NONMEM run). SD magnitudes copied from Supp Final Parameter Estimate
    # (THETAs 13, 14, 15).
    propSd_ko739 <- 0.452; label("KO-739 proportional residual SD (fraction)")     # Table 1 Proportional KO-739 = 45.2%; Supp meta TH13 = 0.452
    addSd_ko739  <- 0.128; label("KO-739 additive residual SD (ng/mL)")            # Table 1 Additive KO-739 = 0.128 ng/mL; Supp meta TH14 = 0.128
    propSd_ko516 <- 0.364; label("KO-516 proportional residual SD (fraction)")     # Table 1 Proportional KO-516 = 36.4%; Supp meta TH15 = 0.364
  })

  model({
    # ==================== Individual PK parameters ====================
    # Parent (ziftomenib).
    # Log-scale covariate effects enter as additive shifts inside exp().
    ka   <- exp(lka + e_conmed_ppi_ka * CONMED_PPI)
    cl   <- exp(lcl + etalcl +
                e_dis_healthy_cl * DIS_HEALTHY +
                e_conmed_cyp3a4_inh_strong_cl * CONMED_CYP3A4_INH_STRONG)
    vc   <- exp(lvc + etalvc)
    q    <- exp(lq)
    vp   <- exp(lvp)
    tlag <- exp(ltlag + e_fed_ltlag * FED)

    # Fraction of oral dose that is absorbed (F1). Logit link keeps F1 in
    # (0, 1) regardless of covariate combinations; FED and PPI enter as
    # additive shifts on the logit scale (paper NONMEM PK block).
    logitfdepot_i <- logitfdepot +
               e_fed_logitfdepot        * FED +
               e_conmed_ppi_logitfdepot * CONMED_PPI
    fd      <- 1 / (1 + exp(-logitfdepot_i))

    # Metabolites (KO-739 and KO-516).
    q_ko739  <- exp(lq_ko739  + etalq_ko739)
    vc_ko739 <- exp(lvc_ko739 + etalvc_ko739 + e_dis_healthy_vc_ko739 * DIS_HEALTHY)
    vp_ko739 <- exp(lvp_ko739)
    cl_ko739 <- exp(lcl_ko739 + etalcl_ko739 +
                    e_conmed_cyp3a4_inh_strong_cl_ko739 * CONMED_CYP3A4_INH_STRONG)

    q_ko516  <- exp(lq_ko516  + etalq_ko516)
    vc_ko516 <- exp(lvc_ko516 + etalvc_ko516 + e_dis_healthy_vc_ko516 * DIS_HEALTHY)
    vp_ko516 <- exp(lvp_ko516)
    cl_ko516 <- exp(lcl_ko516 + etalcl_ko516 +
                    e_conmed_cyp3a4_inh_strong_cl_ko516 * CONMED_CYP3A4_INH_STRONG)

    # Fraction of parent CL routed to metabolites (FM). Logit scale keeps FM
    # in (0, 1) regardless of covariate + eta combinations. Healthy-volunteer
    # covariate enters as an additive shift on the logit scale (paper NONMEM
    # metabolite PK block). The parent-to-metabolite split FM_KO516 = 0.5
    # (1:1 KO-739 : KO-516) is fixed based on in-vitro human-liver-microsome
    # relative abundance (Kura Oncology Inc., data on file per Table 1
    # footnote d); no covariate and no IIV.
    logitfm_ind   <- logitfm + e_dis_healthy_logitfm * DIS_HEALTHY + etalogitfm
    fm            <- 1 / (1 + exp(-logitfm_ind))
    fm_ko516_frac <- 0.5

    # ==================== Micro-rate constants ====================
    # Follows the metabolite NONMEM $PK block verbatim (K25/K52 parent
    # distribution; K20 parent direct elimination; K23/K24 parent flux to
    # metabolite central compartments; K36/K63 and K47/K74 metabolite
    # distributions; K30 and K40 metabolite elimination).
    k12 <- q / vc                     # parent central -> peripheral1
    k21 <- q / vp                     # peripheral1 -> parent central
    k20 <- (cl / vc) * (1 - fm)                          # parent direct elimination
    k23 <- (cl / vc) * fm * (1 - fm_ko516_frac)          # parent -> KO-739 central
    k24 <- (cl / vc) * fm * fm_ko516_frac                # parent -> KO-516 central

    k36 <- q_ko739 / vc_ko739         # KO-739 central -> peripheral1_ko739
    k63 <- q_ko739 / vp_ko739         # peripheral1_ko739 -> KO-739 central
    k30 <- cl_ko739 / vc_ko739        # KO-739 elimination

    k47 <- q_ko516 / vc_ko516         # KO-516 central -> peripheral1_ko516
    k74 <- q_ko516 / vp_ko516         # peripheral1_ko516 -> KO-516 central
    k40 <- cl_ko516 / vc_ko516        # KO-516 elimination

    # ==================== ODE system ====================
    # Absorption chain: depot -> central (parent), with lag time tlag and
    # first-order rate ka.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot -
                      k20 * central -
                      k23 * central -
                      k24 * central -
                      k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # KO-739 central + one peripheral. Mass flux from parent central (k23 *
    # central) enters KO-739 central directly (no molecular-weight scaling in
    # the paper's model; the fixed FM 0.535 and 50/50 split absorb any
    # implicit MW conversion into the identifiability constraints).
    d/dt(central_ko739)     <-  k23 * central -
                                k30 * central_ko739 -
                                k36 * central_ko739 + k63 * peripheral1_ko739
    d/dt(peripheral1_ko739) <-  k36 * central_ko739 - k63 * peripheral1_ko739

    # KO-516 central + one peripheral.
    d/dt(central_ko516)     <-  k24 * central -
                                k40 * central_ko516 -
                                k47 * central_ko516 + k74 * peripheral1_ko516
    d/dt(peripheral1_ko516) <-  k47 * central_ko516 - k74 * peripheral1_ko516

    # ==================== Bioavailability ====================
    # F1 (parent oral bioavailability) applies to the depot compartment. The
    # paper fixes F1 = 0.129 at the reference condition; FED and PPI shift
    # F1 via the logit-scale additive effects above.
    f(depot) <- fd

    # ==================== Observations ====================
    # Dosing is in mg; volumes are in L. amount(mg) / volume(L) = mg/L =
    # 1000 * ng/mL. Multiply by 1000 to report ng/mL consistent with the
    # bioanalytical assay reporting units (LLOQ 0.2 ng/mL).
    Cc       <- 1000 * central       / vc
    Cc_ko739 <- 1000 * central_ko739 / vc_ko739
    Cc_ko516 <- 1000 * central_ko516 / vc_ko516

    Cc       ~ prop(propSd)
    Cc_ko739 ~ add(addSd_ko739) + prop(propSd_ko739)
    Cc_ko516 ~ prop(propSd_ko516)
  })
}
