Lee_2023_patritumab <- function() {
  description <- "Integrated population PK model for the conjugated (anti-HER3-ac-DXd) and unconjugated payload (DXd) of patritumab deruxtecan (HER3-DXd, U3-1402, anti-HER3 antibody-drug conjugate) in adult cancer patients (Lee 2023 ACoP poster). anti-HER3-ac-DXd disposition is a two-compartment model with three parallel elimination pathways from the central compartment: a transient time-decaying linear clearance CL_t(time) = CL_T * exp(-Kdes * time), a non-specific time-dependent linear clearance CL_ns(time) that declines sigmoidally from CL_ss * (1 + Emax) at time = 0 to CL_ss at infinity via a Hill function CL_ns(time) = CL_ss * (1 + Emax * T50^hill / (T50^hill + time^hill)), and a Michaelis-Menten saturable clearance CL_mm = Vmax / (Km + Cc). DXd is a one-compartment model with parallel linear and Michaelis-Menten elimination; its formation rate equals the sum of the three anti-HER3-ac-DXd elimination rates each scaled by a dimensionless fractional-conversion factor Frac_ns / Frac_t / Frac_mm (Frac_ns fixed at 1 as the identifiability anchor)."
  reference <- "Lee M, Wang S, Byrne R, Joshi R, Abutarif M, Garimella T, Li L. Integrated Population Pharmacokinetic Analysis of Conjugated and Unconjugated Payload of Patritumab Deruxtecan in Cancer Patients. American Conference on Pharmacometrics (ACoP) 14 Poster, Oct 2023. https://metrumrg.com/wp-content/uploads/2023/11/34023850_ACoP-2023_Pop_PK_Poster_V01-20Oct2023-1.pdf"
  vignette <- "Lee_2023_patritumab"
  units <- list(
    time          = "hour",
    dosing        = "nmol",
    concentration = "nmol/L"
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Allometric power exponents applied to all clearance terms (e_wt_cl on CL_t, CL_ss, Q, CL_DXd: exponent 0.343, Lee 2023 Table 6 'CL wt') and to all volume terms (e_wt_vc on V1, V2, V1_DXd: exponent 0.475, Lee 2023 Table 6 'V wt'). A second power effect on the DXd fractional-conversion anchor Frac_ns (e_wt_fracns: exponent 0.139, Lee 2023 Table 6 'Frac ns wt'). Reference 60 kg per Lee 2023 Figure 6 caption reference subject. Population median in PPK analysis 57.4 kg; range 32.4 - 113 kg (Lee 2023 Table 2 'All data').",
      source_name        = "Weight"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value (SI units). Power effects on the non-specific (steady-state) clearance CL_ss (e_alb_cl_ss: exponent -0.490, Lee 2023 Table 6 'CLinf alb') and on the DXd fractional-conversion anchor Frac_ns (e_alb_fracns: exponent -0.271, Lee 2023 Table 6 'Frac ns alb'). Reference 40 g/L per Lee 2023 Figure 6 caption reference subject. Population median 39 g/L; range 17 - 49 g/L (Lee 2023 Table 2).",
      source_name        = "Albumin"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate (CKD-EPI)",
      units              = "mL/min/1.73m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "BSA-normalised CKD-EPI estimated glomerular filtration rate (Lee 2023 Table 2 'Estimated glomular filtration (eGFR) (mL/min/1.73m^2)'). Stored under the canonical CRCL register name. Power effects on CL_ss (e_egfr_cl_ss: exponent -0.302, Lee 2023 Table 6 'CLinf egfr') and on Frac_ns (e_egfr_fracns: exponent 0.0350, Lee 2023 Table 6 'Frac ns egfr'). Reference 90 mL/min/1.73m^2 per Lee 2023 Figure 6 caption reference subject. Population median 91.6 mL/min/1.73m^2; range 31.9 - 134 mL/min/1.73m^2 (Lee 2023 Table 2).",
      source_name        = "eGFR"
    ),
    TUMSZ = list(
      description        = "Baseline sum of the longest diameters of target lesions (RECIST)",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value; paper reports as SLD in cm (Lee 2023 Table 2 'Baseline SLD (cm)'). Canonical unit is mm; stored as `TUMSZ_mm = SLD_cm * 10`. Power effects on CL_t (e_sld_cl_time: exponent -0.0758, Lee 2023 Table 6 'CLt sld') and on Kdes (e_sld_kdes: exponent -0.139, Lee 2023 Table 6 'Kdes sld'). Reference 60 mm (= 6 cm) per Lee 2023 Figure 6 caption reference subject. Population median 6.10 cm; range 1.00 - 31.0 cm (Lee 2023 Table 2).",
      source_name        = "Baseline SLD"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male) in the canonical column. The paper's own reference category is female (see notes).",
      notes              = "Time-fixed. Lee 2023 Table 6 reports the male:female ratio as exp(theta) for V1 (1.18; theta26), CL_ss (1.30; theta28), and Frac_ns (0.848; theta38). To store under the canonical SEXF (1 = female, 0 = male), the male-indicator is applied as (1 - SEXF) so SEXF = 1 yields multiplier 1 (the paper's female reference) and SEXF = 0 yields the paper's male multipliers. Reference category in paper covariate forest plots (Lee 2023 Figure 6) is female.",
      source_name        = "Sex"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator, 1 = Asian, 0 = non-Asian",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. Multiplicative effects applied directly via exp(theta * RACE_ASIAN) on V1 (Asian:non-Asian = 0.927; e_asian_vc = log(0.927) = -0.0758), CL_ss (Asian:non-Asian = 1.02), and Frac_ns (Asian:non-Asian = 1.06). Reference category is non-Asian (Lee 2023 Table 6 'Race effect on ... (Asian:non-Asian)' and Figure 6 reference subject 'non-Asian female ...').",
      source_name        = "Race"
    ),
    ECOG_GE1 = list(
      description        = "ECOG performance status >= 1 indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. Multiplicative effects applied via exp(theta * ECOG_GE1) on CL_ss (>0:0 = 1.04; e_ecog_cl_ss = log(1.04)) and Frac_ns (>0:0 = 1.05; e_ecog_fracns = log(1.05)). The paper distinguishes only ECOG > 0 versus ECOG = 0 (Lee 2023 Table 6 'ECOG effect ... (>0:0)'). No subjects had ECOG >= 2 in the analysis dataset (Lee 2023 Table 3).",
      source_name        = "ECOG status"
    ),
    TUMTP_BREAST = list(
      description        = "Breast-cancer tumor-type indicator, 1 = breast cancer, 0 = NSCLC (reference)",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. Multiplicative effects via exp(theta * TUMTP_BREAST) on CL_t (BC:NSCLC = 0.896; e_bc_cl_time = log(0.896)), CL_ss (BC:NSCLC = 0.937), and Frac_ns (BC:NSCLC = 1.04) per Lee 2023 Table 6. The analysis cohort is exactly two tumor types (BC = U31402-J101, NSCLC = U31402-U102; Lee 2023 Table 3) so the indicator captures the full tumor-type contrast. NSCLC is the reference per Figure 6 reference subject 'NSCLC patient ...'.",
      source_name        = "Tumor type"
    ),
    HEPIMP = list(
      description        = "Hepatic-impairment indicator (NCI ODWG mild or moderate vs normal)",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. Multiplicative effects via exp(theta * HEPIMP) on CL_ss (impaired:normal = 0.906; e_hep_cl_ss = log(0.906)) and on Frac_ns (impaired:normal = 1.09; e_hep_fracns = log(1.09)). The paper pools mild and moderate impairment into a single impaired indicator (Lee 2023 Methods: 'hepatic impairment (based on NCIODWG criteria: mild and moderate vs normal)'). No subjects had severe impairment in the analysis dataset (Lee 2023 Table 3). Reference category is normal hepatic function (Lee 2023 Figure 6 reference subject).",
      source_name        = "Hepatic function category (NCI-ODWG)"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 401L,
    n_studies       = 2L,
    n_observations_ac_dxd = 6869L,
    n_observations_dxd    = 6821L,
    age_range       = "29 - 83 years (median 60)",
    age_median      = "60 years",
    weight_range    = "32.4 - 113 kg (median 57.4)",
    weight_median   = "57.4 kg (population median; Lee 2023 Figure 6 reference subject is 60 kg)",
    sex_female_pct  = 80.0,
    race_ethnicity  = c(Asian = 61.1, White = 32.2, Multiple = 3.5, Black = 2.2, AI_AN = 0.5, NH_PI = 0.5),
    disease_state   = "Advanced or metastatic solid tumors. 182 / 401 (45.4%) breast cancer (HER3-expressing, Study U31402-A-J101); 219 / 401 (54.6%) non-small-cell lung cancer (Study U31402-A-U102).",
    dose_range      = "1.6 - 8.0 mg/kg patritumab deruxtecan IV every 3 weeks (Q3W); the Phase 1 dose-escalation NSCLC arm also enrolled an up-titration regimen (Cycle 1: 3.2 mg/kg, Cycle 2: 4.8 mg/kg, Cycle 3+: 6.4 mg/kg).",
    regions         = "Multi-regional: U31402-A-J101 enrolled in Japan with global expansion; U31402-A-U102 enrolled globally including US/EU/Asia (Lee 2023 Methods).",
    ecog_pct        = c(`0` = 52.1, `1` = 47.9, `2` = 0.0),
    hepatic_function = c(Normal = 71.8, Mild = 26.7, Moderate = 0.7, Severe = 0.0, Missing = 0.7),
    renal_function   = c(Normal = 55.1, Mild = 34.7, Moderate = 10.2, Severe = 0.0, Missing = 0.0),
    reference_subject = "Non-Asian female NSCLC patient, age 60, weight 60 kg, ECOG 0, eGFR 90 mL/min/1.73m^2, albumin 40 g/L, baseline SLD 60 mm (= 6 cm), normal hepatic function (Lee 2023 Figure 6 caption).",
    notes           = "ACoP 14 (Oct 2023) conference poster. NONMEM v7.5 sequential fit: anti-HER3-ac-DXd base model, then integrated anti-HER3-ac-DXd + DXd base model, then anti-HER3-ac-DXd covariates, then full integrated covariate model. Full covariate model: covariate point estimates and 95% confidence intervals reported without stepwise hypothesis-testing selection. 500-replicate Monte Carlo VPCs stratified by tumor type. PK observations: 6,869 anti-HER3-ac-DXd and 6,821 DXd samples from 401 cancer subjects across studies U31402-A-J101 and U31402-A-U102. Concentrations reported in molar units (nmol/L) in the structural-parameter tables (Vmax and Km units nmol/L/hr and nmol/L respectively) - dosing therefore needs a mg-to-nmol conversion through the patritumab MW (approx 150,000 g/mol) at the simulation step. Distinct from the peer-reviewed Lu 2022 (J Clin Pharmacol 2023) patritumab deruxtecan popPK model in `Lu_2022_patritumab.R`, which used a different structural form (parallel linear+M-M intact clearance with a within-cycle DAR-modulated DXd release rate) on a partly different cohort (425 subjects from studies J101 / U102 / U202)."
  )

  ini({
    # ============================================================
    # anti-HER3-ac-DXd (intact ADC) structural parameters
    # ============================================================
    # Time unit kept as HOURS so the Hill-on-time CL_ns parameterization
    # (T50 = 1380 hr) and the CL_t decay (Kdes = 0.217 1/hr) preserve the
    # paper's tabulated values exactly. Concentrations are in nmol/L
    # because Vmax (nmol/L/hr) and Km (nmol/L) are reported in molar
    # units (Lee 2023 Tables 4 and 5).
    #
    # Reference subject for typical-value parameters: non-Asian female
    # NSCLC patient with WT 60 kg, ALB 40 g/L, eGFR 90 mL/min/1.73m^2,
    # baseline SLD 60 mm, ECOG 0, normal hepatic function (Lee 2023
    # Figure 6 caption).
    lcl_time     <- log(0.0858);  label("Transient linear clearance CL_T of anti-HER3-ac-DXd at time = 0 (L/hr)")  # Lee 2023 Table 4 exp(theta1) = 0.0858
    lvc          <- log(2.91);    label("Central volume of distribution V1 of anti-HER3-ac-DXd at reference (L)") # Lee 2023 Table 4 exp(theta2) = 2.91
    lq           <- log(0.0221);  label("Intercompartmental clearance Q of anti-HER3-ac-DXd (L/hr)")               # Lee 2023 Table 4 exp(theta3) = 0.0221
    lvp          <- log(3.17);    label("Peripheral volume of distribution V2 of anti-HER3-ac-DXd at reference (L)") # Lee 2023 Table 4 exp(theta4) = 3.17
    lkdes        <- log(0.217);   label("Rate constant of CL_T exponential decline Kdes (1/hr)")                   # Lee 2023 Table 4 exp(theta5) = 0.217
    lcl_ss       <- log(0.0136);  label("Non-specific linear clearance at infinity CL_inf of anti-HER3-ac-DXd (L/hr)") # Lee 2023 Table 4 exp(theta6) = 0.0136
    lemaxclns    <- log(0.603);   label("Max relative increase of CL_ns above CL_ss at time = 0 (unitless)")        # Lee 2023 Table 4 exp(theta7) = 0.603 (CLinf,EMAX)
    lt50clns     <- log(1380);    label("Time of half-maximal CL_ns decline T50 (hr)")                              # Lee 2023 Table 4 exp(theta8) = 1.38e3
    lhillclns    <- log(3.75);    label("Hill exponent of CL_ns sigmoidal time-decay gamma (unitless)")             # Lee 2023 Table 4 exp(theta9) = 3.75
    lvmax        <- log(2.15);    label("Vmax of anti-HER3-ac-DXd Michaelis-Menten elimination (nmol/L/hr)")        # Lee 2023 Table 4 exp(theta10) = 2.15
    lkm          <- log(45.8);    label("Km of anti-HER3-ac-DXd Michaelis-Menten elimination (nmol/L)")             # Lee 2023 Table 4 exp(theta11) = 45.8

    # ============================================================
    # DXd (unconjugated payload) structural parameters
    # ============================================================
    # DXd is a 1-compartment model (Lee 2023 Model Schematic). The paper
    # omits thetas 14 and 15 (no Q_DXd / V2_DXd in the final model).
    lcl_dxd      <- log(4.42);    label("Linear clearance of DXd CL_DXd (L/hr)")                                    # Lee 2023 Table 5 exp(theta12) = 4.42
    lvc_dxd      <- log(5.96);    label("Central volume of distribution of DXd V1_DXd (L)")                          # Lee 2023 Table 5 exp(theta13) = 5.96
    lvmax_dxd    <- log(6.18);    label("Vmax of DXd Michaelis-Menten elimination (nmol/L/hr)")                       # Lee 2023 Table 5 exp(theta16) = 6.18
    lkm_dxd      <- log(0.483);   label("Km of DXd Michaelis-Menten elimination (nmol/L)")                            # Lee 2023 Table 5 exp(theta17) = 0.483
    # Frac_ns is the identifiability anchor for the DXd-formation fractions
    # and is fixed at exp(theta18) = 1 per Lee 2023 Table 5 ("FIXED"). IIV
    # on Frac_ns is estimated (Table 7 Omega(18,18) = 0.0380).
    lfracns      <- fixed(log(1)); label("Scaling factor for fractional conversion from CL_ns to DXd formation (unitless, fixed identifiability anchor)") # Lee 2023 Table 5 exp(theta18) = 1 FIXED
    lfract       <- log(0.272);   label("Scaling factor for fractional conversion from CL_t to DXd formation (unitless)") # Lee 2023 Table 5 exp(theta21) = 0.272
    lfracmm      <- log(0.272);   label("Scaling factor for fractional conversion from CL_mm to DXd formation (unitless)") # Lee 2023 Table 5 exp(theta22) = 0.272 (identical point estimate and 95% CI 0.0670-1.11 to theta21; likely jointly identified)

    # ============================================================
    # Covariate effects on anti-HER3-ac-DXd parameters
    # ============================================================
    # All covariate-effect coefficients on log-clearance/log-volume scale.
    # Categorical-covariate effects encoded as exp(theta * indicator);
    # the canonical SEXF = 1 = female matches the paper's female-reference
    # category, so SEXF = 0 (male) reproduces the paper's male:female ratio
    # via the (1 - SEXF) male-indicator. RACE_ASIAN, ECOG_GE1, TUMTP_BREAST,
    # and HEPIMP each multiply by exp(theta) when the indicator is 1.
    # Continuous-covariate effects are power-of-ratio with reference values
    # 60 kg / 60 mm / 40 g/L / 90 mL/min/1.73 m^2 (Lee 2023 Figure 6 caption).
    e_bc_cl_time      <- log(0.896); label("log ratio of CL_t for BC vs NSCLC (TUMTP_BREAST)")                       # Lee 2023 Table 6 exp(theta24) = 0.896
    e_sld_cl_time     <- -0.0758;    label("Power exponent of TUMSZ on CL_t (unitless)")                              # Lee 2023 Table 6 theta25 = -0.0758
    e_male_vc         <- log(1.18);  label("log ratio of V1 for male vs female (applied via (1 - SEXF))")             # Lee 2023 Table 6 exp(theta26) = 1.18
    e_asian_vc        <- log(0.927); label("log ratio of V1 for Asian vs non-Asian (RACE_ASIAN)")                     # Lee 2023 Table 6 exp(theta27) = 0.927
    e_male_cl_ss      <- log(1.30);  label("log ratio of CL_ss for male vs female (applied via (1 - SEXF))")          # Lee 2023 Table 6 exp(theta28) = 1.30
    e_asian_cl_ss     <- log(1.02);  label("log ratio of CL_ss for Asian vs non-Asian (RACE_ASIAN)")                  # Lee 2023 Table 6 exp(theta29) = 1.02
    e_ecog_cl_ss      <- log(1.04);  label("log ratio of CL_ss for ECOG >= 1 vs ECOG = 0 (ECOG_GE1)")                 # Lee 2023 Table 6 exp(theta30) = 1.04
    e_bc_cl_ss        <- log(0.937); label("log ratio of CL_ss for BC vs NSCLC (TUMTP_BREAST)")                       # Lee 2023 Table 6 exp(theta31) = 0.937
    e_hep_cl_ss       <- log(0.906); label("log ratio of CL_ss for hepatic impairment (mild or moderate) vs normal (HEPIMP)") # Lee 2023 Table 6 exp(theta33) = 0.906
    e_egfr_cl_ss      <- -0.302;     label("Power exponent of CRCL on CL_ss (unitless)")                              # Lee 2023 Table 6 theta34 = -0.302
    e_alb_cl_ss       <- -0.490;     label("Power exponent of ALB on CL_ss (unitless)")                               # Lee 2023 Table 6 theta35 = -0.490
    e_sld_kdes        <- -0.139;     label("Power exponent of TUMSZ on Kdes (unitless)")                              # Lee 2023 Table 6 theta36 = -0.139

    # Covariate effects on Frac_ns (DXd identifiability anchor)
    e_wt_fracns       <-  0.139;     label("Power exponent of WT on Frac_ns (unitless)")                              # Lee 2023 Table 6 theta37 = 0.139
    e_male_fracns     <- log(0.848); label("log ratio of Frac_ns for male vs female (applied via (1 - SEXF))")         # Lee 2023 Table 6 exp(theta38) = 0.848
    e_asian_fracns    <- log(1.06);  label("log ratio of Frac_ns for Asian vs non-Asian (RACE_ASIAN)")                # Lee 2023 Table 6 exp(theta39) = 1.06
    e_ecog_fracns     <- log(1.05);  label("log ratio of Frac_ns for ECOG >= 1 vs ECOG = 0 (ECOG_GE1)")               # Lee 2023 Table 6 exp(theta40) = 1.05
    e_bc_fracns       <- log(1.04);  label("log ratio of Frac_ns for BC vs NSCLC (TUMTP_BREAST)")                     # Lee 2023 Table 6 exp(theta41) = 1.04
    e_hep_fracns      <- log(1.09);  label("log ratio of Frac_ns for hepatic impairment vs normal (HEPIMP)")           # Lee 2023 Table 6 exp(theta43) = 1.09
    e_egfr_fracns     <-  0.0350;    label("Power exponent of CRCL on Frac_ns (unitless)")                            # Lee 2023 Table 6 theta44 = 0.0350
    e_alb_fracns      <- -0.271;     label("Power exponent of ALB on Frac_ns (unitless)")                              # Lee 2023 Table 6 theta45 = -0.271

    # Allometric weight effects shared across clearance and volume terms
    # ("clearance parameters" and "volume parameters" in Lee 2023 Table 6).
    # Applied to: CL_t, CL_ss, Q, CL_DXd (CL); V1, V2, V1_DXd (V).
    # Vmax / Vmax_DXd / Km / Km_DXd do not carry an explicit weight effect
    # in the paper and are left as size-independent (see Errata).
    e_wt_cl           <-  0.343;     label("Allometric power exponent of WT on linear clearance parameters (unitless)") # Lee 2023 Table 6 theta46 = 0.343
    e_wt_vc           <-  0.475;     label("Allometric power exponent of WT on volume parameters (unitless)")           # Lee 2023 Table 6 theta47 = 0.475

    # ============================================================
    # Inter-individual variability (Lee 2023 Table 7)
    # ============================================================
    # Reported as Omega(i,i) values on the log-scale (variances of the
    # NONMEM ETAs). CV% = sqrt(exp(Omega) - 1) * 100. No correlation
    # blocks are reported in the poster, so all etas are uncorrelated.
    # No IIV reported on Kdes (Omega(5,5)), Emax (Omega(7,7)), gamma
    # (Omega(9,9)), Vmax / Km (Omega(10,10) / Omega(11,11)), Vmax_DXd /
    # Km_DXd (Omega(16,16) / Omega(17,17)), Frac_t (Omega(21,21)), or
    # Frac_mm (Omega(22,22)).
    etalcl_time   ~ 0.299    # Lee 2023 Table 7 Omega(1,1) = 0.299 (CV% 59.0)
    etalvc        ~ 0.0205   # Lee 2023 Table 7 Omega(2,2) = 0.0205 (CV% 14.4)
    etalq         ~ 0.431    # Lee 2023 Table 7 Omega(3,3) = 0.431 (CV% 73.4)
    etalvp        ~ 0.112    # Lee 2023 Table 7 Omega(4,4) = 0.112 (CV% 34.5)
    etalcl_ss     ~ 0.127    # Lee 2023 Table 7 Omega(6,6) = 0.127 (CV% 36.9)
    etalt50clns   ~ 0.702    # Lee 2023 Table 7 Omega(8,8) = 0.702 (CV% 101)
    etalcl_dxd    ~ 0.108    # Lee 2023 Table 7 Omega(12,12) = 0.108 (CV% 33.7)
    etalvc_dxd    ~ 0.00515  # Lee 2023 Table 7 Omega(13,13) = 0.00515 (CV% 7.18)
    etalfracns    ~ 0.0380   # Lee 2023 Table 7 Omega(18,18) = 0.0380 (CV% 19.7)

    # ============================================================
    # Residual error (Lee 2023 Table 7)
    # ============================================================
    # NONMEM Y = IPRED * (1 + EPS) with Sigma reported on the linear-scale
    # proportional error variance (CV% = sqrt(Sigma) * 100). Maps to a
    # proportional residual model in nlmixr2 with propSd = sqrt(Sigma).
    propSd     <- sqrt(0.0342);  label("Proportional residual error on anti-HER3-ac-DXd Cc (fraction)")  # Lee 2023 Table 7 Sigma(1,1) = 0.0342 (CV% 18.5)
    propSd_dxd <- sqrt(0.0814);  label("Proportional residual error on DXd Cc_dxd (fraction)")           # Lee 2023 Table 7 Sigma(2,2) = 0.0814 (CV% 28.5)
  })

  model({
    # ============================================================
    # Derived covariate terms
    # ============================================================
    # Lee 2023 encodes sex as male-indicator (female = 0 reference) and
    # tumor type / race / hepatic / ECOG as 1 = test category. Canonical
    # SEXF = 1 = female, so (1 - SEXF) reproduces the paper's male
    # multiplier while keeping SEXF in its canonical 1=female form.
    sex_male <- 1 - SEXF

    # ============================================================
    # Individual anti-HER3-ac-DXd parameters
    # ============================================================
    # Reference covariates: WT 60 kg, ALB 40 g/L, CRCL 90 mL/min/1.73m^2,
    # TUMSZ 60 mm, SEXF = 1 (female), RACE_ASIAN = 0, ECOG_GE1 = 0,
    # TUMTP_BREAST = 0 (NSCLC), HEPIMP = 0 (normal). All categorical
    # multipliers reduce to 1 at the reference; all continuous power
    # terms reduce to (cov / cov_ref)^0 = 1.
    cl_time_typ <- exp(lcl_time + etalcl_time) *
      (WT    / 60)^e_wt_cl *
      (TUMSZ / 60)^e_sld_cl_time *
      exp(e_bc_cl_time * TUMTP_BREAST)

    vc <- exp(lvc + etalvc) *
      (WT / 60)^e_wt_vc *
      exp(e_male_vc  * sex_male) *
      exp(e_asian_vc * RACE_ASIAN)

    q <- exp(lq + etalq) *
      (WT / 60)^e_wt_cl

    vp <- exp(lvp + etalvp) *
      (WT / 60)^e_wt_vc

    kdes <- exp(lkdes) *
      (TUMSZ / 60)^e_sld_kdes

    cl_ss <- exp(lcl_ss + etalcl_ss) *
      (WT   / 60)^e_wt_cl *
      (CRCL / 90)^e_egfr_cl_ss *
      (ALB  / 40)^e_alb_cl_ss *
      exp(e_male_cl_ss  * sex_male) *
      exp(e_asian_cl_ss * RACE_ASIAN) *
      exp(e_ecog_cl_ss  * ECOG_GE1) *
      exp(e_bc_cl_ss    * TUMTP_BREAST) *
      exp(e_hep_cl_ss   * HEPIMP)

    emax_clns <- exp(lemaxclns)
    t50_clns  <- exp(lt50clns + etalt50clns)
    hill_clns <- exp(lhillclns)

    vmax <- exp(lvmax)
    km   <- exp(lkm)

    # ============================================================
    # Individual DXd parameters
    # ============================================================
    cl_dxd <- exp(lcl_dxd + etalcl_dxd) *
      (WT / 60)^e_wt_cl

    vc_dxd <- exp(lvc_dxd + etalvc_dxd) *
      (WT / 60)^e_wt_vc

    vmax_dxd <- exp(lvmax_dxd)
    km_dxd   <- exp(lkm_dxd)

    fracns <- exp(lfracns + etalfracns) *
      (WT   / 60)^e_wt_fracns *
      (CRCL / 90)^e_egfr_fracns *
      (ALB  / 40)^e_alb_fracns *
      exp(e_male_fracns  * sex_male) *
      exp(e_asian_fracns * RACE_ASIAN) *
      exp(e_ecog_fracns  * ECOG_GE1) *
      exp(e_bc_fracns    * TUMTP_BREAST) *
      exp(e_hep_fracns   * HEPIMP)

    fract  <- exp(lfract)
    fracmm <- exp(lfracmm)

    # ============================================================
    # Time-varying clearance terms (Lee 2023 Methods)
    # ============================================================
    # CL_t exponential decay: CL_t(time) = CL_T * exp(-Kdes * time).
    # CL_ns Hill-on-time decline: CL_ns(time) = CL_ss * (1 + Emax * T50^hill / (T50^hill + time^hill)).
    # At time = 0 the Hill factor is 1 (T50^hill / T50^hill) so
    # CL_ns(0) = CL_ss * (1 + Emax). As time -> infinity the Hill factor
    # vanishes and CL_ns -> CL_ss. With the typical-value estimates this
    # reproduces the paper's text 0.0136 * (1 + 0.603) = 0.0218 ~ 0.0217
    # L/hr at time = 0 declining toward CL_ss = 0.0136 L/hr.
    # `time` is the rxode2 simulation clock (hours since first dose).
    cl_t_now  <- cl_time_typ * exp(-kdes * time)
    cl_ns_now <- cl_ss * (1 + emax_clns * t50_clns^hill_clns /
                              (t50_clns^hill_clns + time^hill_clns))

    # ============================================================
    # ODE system: 2-cmt anti-HER3-ac-DXd + 1-cmt DXd
    # ============================================================
    # central in nmol, vc in L, so Cc = central / vc is in nmol/L which
    # matches the Vmax / Km units in the paper's structural tables.
    Cc <- central / vc

    # Pathway elimination rates (nmol/hr). The Michaelis-Menten rate is
    # written in amount form vmax * central / (km + Cc) so the resulting
    # mass-flux units (nmol/hr) are consistent with the linear pathways
    # (cl * Cc, also nmol/hr).
    rate_t  <- cl_t_now  * Cc
    rate_ns <- cl_ns_now * Cc
    rate_mm <- vmax * central / (km + Cc)

    # DXd formation rate (nmol/hr). Each pathway's elimination rate is
    # scaled by its dimensionless conversion fraction (Frac_ns fixed at 1
    # as the identifiability anchor; Frac_t and Frac_mm estimated).
    formation_dxd <- fracns * rate_ns + fract * rate_t + fracmm * rate_mm

    # DXd elimination (nmol/hr): linear CL_DXd plus Michaelis-Menten.
    Cc_dxd     <- central_dxd / vc_dxd
    rate_d_lin <- cl_dxd * Cc_dxd
    rate_d_mm  <- vmax_dxd * central_dxd / (km_dxd + Cc_dxd)

    d/dt(central)     <- -rate_t - rate_ns - rate_mm -
                          q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <-  q / vc * central - q / vp * peripheral1
    d/dt(central_dxd) <-  formation_dxd - rate_d_lin - rate_d_mm

    # ============================================================
    # Observations and error model
    # ============================================================
    Cc     ~ prop(propSd)
    Cc_dxd ~ prop(propSd_dxd)
  })
}
