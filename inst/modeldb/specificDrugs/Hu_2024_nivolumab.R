Hu_2024_nivolumab <- function() {
  description <- "Two-compartment population PK model with time-varying clearance for intravenous nivolumab (anti-PD-1 IgG4) in a pooled adult and pediatric (1-17 years) oncology population, with tumor-type-dependent pediatric effects on baseline clearance beyond body size (Hu 2024)"
  reference <- "Hu Z, Liu S, Zhao Y, Du S, Hamuro L, Shen J, Roy A, Zhu L. Nivolumab and ipilimumab population pharmacokinetics in support of pediatric dose recommendations-Going beyond the body-size effect. CPT Pharmacometrics Syst Pharmacol. 2024;13(3):476-493. doi:10.1002/psp4.13098"
  vignette <- "Hu_2024_nivolumab_ipilimumab_pediatric"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central     = list(analyte = "nivolumab", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "nivolumab", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on baseline CL with reference 75 kg, and on Q with the same exponent (Hu 2024 Results equations for CL0TV and Q; the Q exponent was constrained equal to the CL exponent per Figure 1 caption). Reference 75 kg is the approximate dataset median (Table 4 Note). Missing values were imputed to the reference in the source NONMEM code (Appendix File S2).",
      source_name        = "WTB"
    ),
    LBM = list(
      description        = "Baseline lean body mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on VC with reference 55 kg, and on VP with the same exponent (Hu 2024 Results equations for VC and VP; the VP exponent was constrained equal to the VC exponent per Figure 1 caption). Table 1 Note: LBM was estimated using the Boer equation (adults) and the Peter equation (children). Missing values were imputed to the reference in the source NONMEM code (Appendix File S2).",
      source_name        = "ELBMB"
    ),
    CRCL = list(
      description        = "Baseline estimated glomerular filtration rate, BSA-normalized",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on baseline CL with reference 90 mL/min/1.73 m^2 (Hu 2024 Results, CL0TV equation and the reference-value list following it). Stored under the canonical CRCL; the source column name is eGFR / EGFRB. The eGFR estimating equation is not restated in Hu 2024; it is inherited from the Bajaj 2017 adult base model that this analysis re-estimated.",
      source_name        = "EGFRB"
    ),
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Exponential effect on baseline CL (exp(-0.0998) = 0.905) and on VC (exp(0.0195) = 1.020). The source NONMEM column SEXN is coded 1 = male (reference), 2 = female; SEXF = 1 reproduces the paper's female-vs-male contrast (Appendix File S2: 'IF (SEX_I .EQ. 2) TVCL = TVCL * EXP(CL_FEMALE); reference is SEX=1 (Male)').",
      source_name        = "SEXN"
    ),
    ECOG_GE1 = list(
      description        = "Baseline ECOG performance-status indicator (1 if performance status > 0, else 0)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ECOG performance status 0; fully active)",
      notes              = "Exponential effect on baseline CL (exp(0.166) = 1.181) and additive effect on the time-varying-CL Emax parameter (-0.157). Table 1 Note: performance status was based on the Eastern Cooperative Oncology Group (ECOG) Performance Status Scale. Hu 2024 collapses the four-level scale to PS = 0 versus PS > 0; missing values were imputed to PS = 1 in the source NONMEM code (Appendix File S2).",
      source_name        = "PS"
    ),
    RACE_BLACK = list(
      description        = "Race indicator: Black / African American",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White or Other; Hu 2024 pools White and Other into a single reference group)",
      notes              = "Exponential effect on baseline CL (CL_RAAA = 0.0693; 95% CI -0.0237 to 0.159, i.e. not statistically significant). Source column RACEN coded 1 = White/Other (reference), 2 = African American, 3 = Asian. Renamed to the canonical RACE_BLACK per covariate-columns.md.",
      source_name        = "RACEN"
    ),
    RACE_ASIAN = list(
      description        = "Race indicator: Asian",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White or Other; same composite reference as RACE_BLACK)",
      notes              = "Exponential effect on baseline CL (CL_RAAS = 0.00333; RSE 1090%, 95% CI -0.0693 to 0.0786, i.e. effectively null). Source column RACEN level 3. Renamed to the canonical RACE_ASIAN per covariate-columns.md.",
      source_name        = "RACEN"
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used only through the three age bands that define the paper's patient-population and volume covariates: young pediatric (< 12 years), adolescent (12-17 years), adult (>= 18 years). Appendix File S2 applies the splits as AGE <= 11 and 12 <= AGE <= 17. The numeric power effect of age on CL (CL_AGE, THETA(8)) was fixed to zero in the selected Full1c model, so age enters this model only categorically.",
      source_name        = "AGE"
    ),
    TUMTP_LYMPH = list(
      description        = "Tumor-type indicator: lymphoma (pooled classical Hodgkin, Hodgkin, and non-Hodgkin lymphoma)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (melanoma, CNS tumor, or other solid tumor)",
      notes              = "Exponential effects on baseline CL: adults exp(-0.382) = 0.682 (32% lower than adult melanoma); pediatric subjects (< 18 years) exp(-0.411) = 0.663 (34% lower). Additive effect on Emax in adults only (+0.132); pediatric lymphoma Emax was assumed equal to adult melanoma and fixed at zero (Methods; Appendix File S2 THETA(39) = 0 FIX). Hu 2024 Methods: 'cHL, HL, and non-HL were combined into a single HL category because of the small sample size' (adult n = 274 from Table 1; pediatric n = 46 comprising classical Hodgkin lymphoma n = 31, non-Hodgkin lymphoma n = 9, Hodgkin lymphoma n = 6).",
      source_name        = "POPN"
    ),
    TUMTP_CNS_PRIM = list(
      description        = "Tumor-type indicator: primary central nervous system tumor",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (melanoma, lymphoma, or other solid tumor)",
      notes              = "Exponential effects on baseline CL: adults exp(-0.578) = 0.561 (44% lower than adult melanoma); pediatric subjects (< 18 years) exp(-0.801) = 0.449 (55% lower). The adult CNS-tumor cohort is essentially all glioblastoma (Table 1: adult GBM n = 556, of which GBM n = 542), which is why Table 4 labels the coefficient CL_GBM. Pediatric CNS tumors additionally carry an Emax offset of +0.696, and adult CNS tumors have Emax set to exactly zero (stationary CL) in the source NONMEM code (Appendix File S2: 'IF (POP_I .EQ. 4) EMAX = 0'). Pediatric composition (Table 1 footnote f): diffuse intrinsic pontine glioma n = 34, medulloblastoma n = 28, high-grade glioma n = 20, ependymoma n = 19, diffuse midline glioma n = 9, atypical teratoid rhabdoid tumor n = 7, and others.",
      source_name        = "POPN"
    ),
    TUMTP_OTHER = list(
      description        = "Tumor-type indicator: other solid tumors (the residual non-melanoma, non-lymphoma, non-CNS pool)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (melanoma, lymphoma, or CNS tumor)",
      notes              = "Exponential effects on baseline CL: adults exp(0.00699) = 1.007 (no meaningful difference from adult melanoma); adolescents 12-17 years exp(-0.223) = 0.800 (20% lower); young pediatric subjects < 12 years exp(-0.580) = 0.560 (44% lower). Additive effect on Emax in adults only (+0.118); pediatric solid-tumor Emax was assumed equal to adult melanoma and fixed at zero (Methods; Appendix File S2 THETA(38) = 0 FIX). Per-paper composition -- adult other (n = 227, Table 1 footnote c): NSCLC n = 139, renal cell carcinoma n = 35, colorectal cancer n = 18, prostate cancer n = 11, others; pediatric other solid tumors (n = 79, Table 1 footnote d): osteosarcoma n = 19, neuroblastoma n = 18, rhabdomyosarcoma n = 17, Ewing sarcoma n = 10, melanoma n = 1, other n = 14.",
      source_name        = "POPN"
    ),
    CONMED_IPI_1Q3W = list(
      description        = "Coadministration regimen: nivolumab + ipilimumab 1 mg/kg every 3 weeks",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (nivolumab monotherapy or any other combination regimen)",
      notes              = "Exponential effect on baseline CL (exp(0.0973) = 1.102; 95% CI -0.00235 to 0.208). Source column COMBO level 1 (Appendix File S2). Table 1: n = 122 of 2325 received nivolumab + ipilimumab 1 mg/kg q3w.",
      source_name        = "COMBO"
    ),
    CONMED_IPI_3Q3W = list(
      description        = "Coadministration regimen: nivolumab + ipilimumab 3 mg/kg every 3 weeks",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (nivolumab monotherapy or any other combination regimen)",
      notes              = "Exponential effect on baseline CL (exp(0.349) = 1.418; the largest single covariate effect on nivolumab CL in this model). Source column COMBO level 2 (Appendix File S2). Table 1: n = 385 of 2325.",
      source_name        = "COMBO"
    ),
    CONMED_BRENTUXIMAB = list(
      description        = "Coadministration regimen: nivolumab + brentuximab vedotin",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (nivolumab monotherapy or any other combination regimen)",
      notes              = "Exponential effect on baseline CL (exp(0.132) = 1.141; 95% CI -0.0572 to 0.282). Source column COMBO level 4 (Appendix File S2). The regimen is brentuximab vedotin 1.8 mg/kg q3w for 4 doses with nivolumab 3 mg/kg q3w in the pediatric / young-adult classical Hodgkin lymphoma study CA209744 (Table S1). Table 1: n = 44 of 2325.",
      source_name        = "COMBO"
    ),
    CONMED_IPI_ANY = list(
      description        = "Any-ipilimumab-coadministration indicator (regimen-agnostic)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no ipilimumab coadministration)",
      notes              = "Additive effect on the time-varying-CL Emax parameter (-0.124; 95% CI -0.361 to 0.0274). Source column IPICO / COMBO levels 1-3 (Appendix File S2: 'IF (COMBO_I.LE.3.AND.COMBO_I.GE.1) EMAX = EMAX + EMAX_IPICO'). Coexists with CONMED_IPI_1Q3W and CONMED_IPI_3Q3W because those act on baseline CL rather than on Emax.",
      source_name        = "COMBO"
    )
  )

  # Covariate columns present in the Hu 2024 nivolumab analysis dataset that
  # carry no coefficient in the selected Full1c model -- either because they
  # define the reference level of a retained factor, or because their
  # coefficient was declared 0 FIX in the source control stream. Documented
  # for provenance only; not referenced in model().
  covariatesDataExcluded <- list(
    TUMTP_MEL = list(
      description = "Tumor-type indicator: melanoma -- the reference level of the patient-population factor",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Fourth level of the patient-population factor (melanoma / lymphoma / CNS tumor /",
        "other solid tumor) that Hu 2024 crosses with the age bands. Adult melanoma is the",
        "reference patient population, so TUMTP_MEL carries no CL or Emax coefficient at any",
        "age in the nivolumab model and therefore does not appear in model(); a subject with",
        "all of TUMTP_LYMPH, TUMTP_CNS_PRIM, and TUMTP_OTHER equal to 0 is a melanoma subject.",
        "Table 1 tumor-type counts: MEL 994, HL 320, CNST 706, other ST 305 of N = 2325.",
        "The single pediatric melanoma subject in the nivolumab analysis dataset was assigned",
        "to the pediatric solid-tumor group by the analysis plan (Table 1 footnote d), so a",
        "pediatric melanoma record does not occur in the source data. The companion model",
        "Hu_2024_ipilimumab.R does estimate a pediatric-melanoma CL effect and keeps",
        "TUMTP_MEL in covariateData."
      )
    ),
    CONMED_IPI_1Q6W = list(
      description = "Coadministration regimen: nivolumab + ipilimumab 1 mg/kg every 6 weeks",
      units       = "(binary)",
      type        = "binary",
      notes       = "COMBO level 3 in the source NONMEM code (Appendix File S2); THETA(30) CL_I1Q6 was declared 0 FIX and does not appear in Table 4. Retained as a data column but with no estimated effect in this analysis."
    ),
    CONMED_RADIOTHERAPY = list(
      description = "Coadministration regimen: nivolumab + radiotherapy (adult glioblastoma studies)",
      units       = "(binary)",
      type        = "binary",
      notes       = "COMBO level 5 in the source NONMEM code (Appendix File S2); THETA(32) CL_RT was declared 0 FIX and does not appear in Table 4. Table 1 records n = 314 receiving nivolumab + radiotherapy."
    ),
    CONMED_RADIOTHERAPY_TEMOZOLOMIDE = list(
      description = "Coadministration regimen: nivolumab + radiotherapy + temozolomide (adult glioblastoma)",
      units       = "(binary)",
      type        = "binary",
      notes       = "COMBO level 6 in the source NONMEM code (Appendix File S2); THETA(33) CL_RTTMZ was declared 0 FIX and does not appear in Table 4. Table 1 records n = 52 receiving nivolumab + radiotherapy + temozolomide."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 2325L,
    n_studies      = 13L,
    n_observations = 13104L,
    n_pediatric    = 275L,
    age_range      = "1-90 years (pediatric 1-17 years; adult >= 18 years)",
    age_median     = "55 years overall; pediatric solid tumors 12 y, pediatric lymphoma 15 y, pediatric CNS tumors 10 y",
    weight_range   = "9.3-168 kg",
    weight_median  = "76.2 kg overall; pediatric solid tumors 43.2 kg, pediatric lymphoma 58.2 kg, pediatric CNS tumors 33 kg",
    sex_female_pct = 36.5,
    race_ethnicity = c(White = 90.3, `Black/African American` = 2.8, Asian = 3.7, Other = 3.2),
    disease_state  = "Advanced melanoma (n = 994), Hodgkin / non-Hodgkin lymphoma (n = 320), central nervous system tumors including adult glioblastoma (n = 706), and other solid tumors (n = 305). Analysis populations: adult MEL n = 993, adult HL n = 274, adult GBM n = 556, adult other ST n = 227, pediatric ST n = 79, pediatric HL n = 46, pediatric CNST n = 150.",
    dose_range     = "Nivolumab 0.1-20 mg/kg or 240/480 mg flat intravenous infusion q2w, q3w, or q4w, alone or combined with ipilimumab 1 or 3 mg/kg q3w, brentuximab vedotin, radiotherapy, or radiotherapy plus temozolomide",
    regions        = "Pooled global phase I, I/II, II, and III studies (13 trials; Table S1)",
    performance_status = "ECOG PS 0 50.5%, PS 1 45.1%, PS 2 4.3%, PS 3 0.04%",
    renal_function = "Baseline eGFR mean 93.8 (SD 22.8) mL/min/1.73 m^2, median 92.9",
    body_composition = "Baseline lean body mass mean 54.5 (SD 13.8) kg, median 55.8 kg; estimated by the Boer equation (adults) and the Peter equation (children)",
    notes          = "Baseline demographics per Hu 2024 Table 1 (N = 2325 across 13 nivolumab studies, of whom 275 were pediatric). Three phase I/II studies contributed the pediatric data (ADVL1412 / CA209070, CA209908, CA209744; Table S1)."
  )

  ini({
    # Structural parameters at the reference covariates: baseline body weight
    # 75 kg, lean body mass 55 kg, performance status 0, baseline eGFR
    # 90 mL/min/1.73 m^2, male, White/Other race, adult melanoma on nivolumab
    # monotherapy (Hu 2024 Results, reference-value list; Table 4 Note).
    # CL and Q are reported in mL/h in Table 4 and are converted to L/day
    # (x 24 / 1000) because this model keeps time in days, matching the
    # sibling nivolumab models in nlmixr2lib.
    lcl <- log(9.66 * 24 / 1000); label("Baseline clearance CL0_REF at the reference covariates (L/day)")     # Hu 2024 Table 4: CL0_REF = 9.66 mL/h
    lvc <- log(4.01);             label("Central volume VC_REF at the reference covariates (L)")               # Hu 2024 Table 4: VC_REF = 4.01 L
    lq  <- log(35.9 * 24 / 1000); label("Intercompartmental clearance Q_REF at the reference covariates (L/day)") # Hu 2024 Table 4: Q_REF = 35.9 mL/h
    lvp <- log(2.77);             label("Peripheral volume VP_REF at the reference covariates (L)")            # Hu 2024 Table 4: VP_REF = 2.77 L

    # Time-varying clearance (Hu 2024 Results, CL_i(t) equation):
    #   CL(t) = CL0TV * exp(Emax * t^HILL / (T50^HILL + t^HILL))
    # Emax is on the linear scale (negative = CL decreases with time); the
    # steady-state-to-baseline CL ratio is exp(Emax) (Figure 1 caption).
    # T50 is reported in hours in Table 4 and converted to days.
    cl_time_max   <-      -0.298;      label("Reference Emax of the time-varying-CL function (unitless; negative = CL decreases with time)") # Hu 2024 Table 4: EMAX_REF = -0.298
    lcl_t50       <- log(2670 / 24);   label("log T50 - time at which the change in CL is 50% of Emax (log days)")                           # Hu 2024 Table 4: T50 = 2670 h
    lcl_time_hill <- log(2.32);        label("log HILL - sigmoidicity of the time-on-CL function (log unitless)")                            # Hu 2024 Table 4: HILL = 2.32

    # Continuous covariate effects (power form (cov / reference)^exponent).
    e_wt_cl   <- 0.630;  label("Power exponent of baseline body weight on CL, also applied to Q (unitless)") # Hu 2024 Table 4: CL_WTB = 0.630
    e_egfr_cl <- 0.0982; label("Power exponent of baseline eGFR on CL (unitless)")                           # Hu 2024 Table 4: CL_eGFR = 0.0982
    e_lbm_vc  <- 0.932;  label("Power exponent of lean body mass on VC, also applied to VP (unitless)")      # Hu 2024 Table 4: VC_LBM = 0.932

    # Categorical covariate effects on baseline CL (exponential form
    # exp(theta * indicator); all indicators are 0 / 1).
    e_sexf_cl     <- -0.0998; label("Exponential coefficient of female sex on CL (unitless)")                              # Hu 2024 Table 4: CL_SEX = -0.0998
    e_ecog_cl     <-  0.166;  label("Exponential coefficient of ECOG performance status > 0 on CL (unitless)")             # Hu 2024 Table 4: CL_PS = 0.166
    e_black_cl    <-  0.0693; label("Exponential coefficient of Black / African American race on CL (unitless)")           # Hu 2024 Table 4: CL_RAAA = 0.0693 (95% CI spans zero)
    e_asian_cl    <-  0.00333; label("Exponential coefficient of Asian race on CL (unitless)")                             # Hu 2024 Table 4: CL_RAAS = 0.00333 (RSE 1090%)
    e_ipi1q3w_cl  <-  0.0973; label("Exponential coefficient of ipilimumab 1 mg/kg q3w coadministration on CL (unitless)") # Hu 2024 Table 4: CL_I1Q3 = 0.0973
    e_ipi3q3w_cl  <-  0.349;  label("Exponential coefficient of ipilimumab 3 mg/kg q3w coadministration on CL (unitless)") # Hu 2024 Table 4: CL_I3Q3 = 0.349
    e_bv_cl       <-  0.132;  label("Exponential coefficient of brentuximab vedotin coadministration on CL (unitless)")    # Hu 2024 Table 4: CL_BVCO = 0.132

    # Patient-population (tumor type x age band) effects on baseline CL,
    # relative to the adult melanoma reference (Hu 2024 Results, CL0TV
    # equation; Appendix File S2 assigns each term by POPN and the AGE
    # bands <= 11 and 12-17 years).
    e_hl_cl      <- -0.382; label("Exponential coefficient of adult lymphoma on CL (unitless)")                          # Hu 2024 Table 4: CL_HL = -0.382
    e_gbm_cl     <- -0.578; label("Exponential coefficient of adult CNS tumor (glioblastoma) on CL (unitless)")          # Hu 2024 Table 4: CL_GBM = -0.578
    e_oth_cl     <-  0.00699; label("Exponential coefficient of adult other solid tumors on CL (unitless)")              # Hu 2024 Table 4: CL_OTH = 0.00699 (RSE 529%)
    e_pedst_cl   <- -0.580; label("Exponential coefficient of young pediatric (< 12 y) solid tumors on CL (unitless)")   # Hu 2024 Table 4: CL_PEDST = -0.580
    e_adost_cl   <- -0.223; label("Exponential coefficient of adolescent (12-17 y) solid tumors on CL (unitless)")       # Hu 2024 Table 4: CL_ADOST = -0.223
    e_pedhl_cl   <- -0.411; label("Exponential coefficient of pediatric (< 18 y) lymphoma on CL (unitless)")             # Hu 2024 Table 4: CL_PEDHL = -0.411
    e_pedcnst_cl <- -0.801; label("Exponential coefficient of pediatric (< 18 y) CNS tumors on CL (unitless)")           # Hu 2024 Table 4: CL_PEDCNST = -0.801

    # Categorical covariate effects on VC (exponential form). The age-band
    # effects on VC apply to every tumor type (Appendix File S2 gates them on
    # AGE alone).
    e_sexf_vc <-  0.0195; label("Exponential coefficient of female sex on VC (unitless)")                       # Hu 2024 Table 4: VC_SEX = 0.0195 (RSE 95.5%)
    e_ped_vc  <- -0.277;  label("Exponential coefficient of young pediatric age (< 12 y) on VC (unitless)")     # Hu 2024 Table 4: VC_PED = -0.277
    e_ado_vc  <- -0.273;  label("Exponential coefficient of adolescent age (12-17 y) on VC (unitless)")         # Hu 2024 Table 4: VC_ADO = -0.273

    # Additive covariate effects on Emax (Hu 2024 Results, EMAX_i equation).
    # Pediatric solid-tumor and pediatric lymphoma Emax offsets were assumed
    # equal to adult melanoma and fixed at zero (Methods; Appendix File S2
    # THETA(38) and THETA(39) = 0 FIX), so no term appears for them here.
    e_ecog_cl_time_max    <- -0.157; label("Additive effect of ECOG performance status > 0 on Emax (unitless)")      # Hu 2024 Table 4: EMAX_PS = -0.157
    e_ipico_cl_time_max   <- -0.124; label("Additive effect of any ipilimumab coadministration on Emax (unitless)")  # Hu 2024 Table 4: EMAX_IPICO = -0.124
    e_hl_cl_time_max      <-  0.132; label("Additive effect of adult lymphoma on Emax (unitless)")                   # Hu 2024 Table 4: EMAX_HL = 0.132
    e_oth_cl_time_max     <-  0.118; label("Additive effect of adult other solid tumors on Emax (unitless)")         # Hu 2024 Table 4: EMAX_OTH = 0.118
    e_pedcnst_cl_time_max <-  0.696; label("Additive effect of pediatric (< 18 y) CNS tumors on Emax (unitless)")    # Hu 2024 Table 4: EMAX_PEDCNST = 0.696

    # Inter-individual variability. CL and VC form a log-normal 2x2 block.
    # Q and VP form a second 2x2 block declared "$OMEGA BLOCK(2) SAME" in the
    # source control stream (Appendix File S2), i.e. constrained to the same
    # variances and covariance as the CL / VC block; that is why Table 4 does
    # not list separate omega^2 values for Q and VP.
    # Emax carries an independent ADDITIVE (normally distributed) eta on the
    # linear scale (Hu 2024 Results, EMAX_i equation).
    etalcl + etalvc ~ c(0.108,
                        0.0220, 0.0751)   # Hu 2024 Table 4: omega^2_CL = 0.108, cov(CL,VC) = 0.0220, omega^2_VC = 0.0751
    etalq + etalvp  ~ c(0.108,
                        0.0220, 0.0751)   # Appendix File S2: $OMEGA BLOCK(2) SAME - Q / VP block constrained equal to the CL / VC block
    etacl_time_max  ~ 0.160               # Hu 2024 Table 4: omega^2_EMAX = 0.160

    # Residual error. Table 4 reports the proportional term as 0.199. The
    # source control stream fixes $SIGMA to 1 and forms the weight as
    # sqrt(F^2 * PERR^2 + AERR^2) with AERR = 0 FIX (Appendix File S2), so
    # 0.199 is the proportional standard deviation (a fraction), not a variance.
    propSd <- 0.199; label("Proportional residual error (fraction)")  # Hu 2024 Table 4: Proportional error term = 0.199
  })
  model({
    # 1. Age-band indicators. Appendix File S2 applies the pediatric splits as
    # AGE <= 11 and 12 <= AGE <= 17; adults are the complement.
    ped   <- 1 - (AGE >= 12)                  # young pediatric, < 12 years
    ado   <- (AGE >= 12) * (1 - (AGE >= 18))  # adolescent, 12-17 years
    adult <- (AGE >= 18)                      # adult, >= 18 years
    child <- 1 - adult                        # pediatric, < 18 years

    # Patient-population effect on baseline CL, relative to adult melanoma.
    # Adult tumor types are contrasted against adult melanoma; the pediatric
    # solid-tumor effect is split at 12 years while the pediatric lymphoma and
    # pediatric CNS-tumor effects pool the whole < 18-year range.
    cl_pop <-
      e_hl_cl      * TUMTP_LYMPH    * adult +
      e_gbm_cl     * TUMTP_CNS_PRIM * adult +
      e_oth_cl     * TUMTP_OTHER    * adult +
      e_pedst_cl   * TUMTP_OTHER    * ped +
      e_adost_cl   * TUMTP_OTHER    * ado +
      e_pedhl_cl   * TUMTP_LYMPH    * child +
      e_pedcnst_cl * TUMTP_CNS_PRIM * child

    # 2. Individual parameters.
    cl0 <- exp(lcl + etalcl) *
      (WT   / 75)^e_wt_cl *
      (CRCL / 90)^e_egfr_cl *
      exp(e_sexf_cl    * SEXF) *
      exp(e_ecog_cl    * ECOG_GE1) *
      exp(e_black_cl   * RACE_BLACK) *
      exp(e_asian_cl   * RACE_ASIAN) *
      exp(e_ipi1q3w_cl * CONMED_IPI_1Q3W) *
      exp(e_ipi3q3w_cl * CONMED_IPI_3Q3W) *
      exp(e_bv_cl      * CONMED_BRENTUXIMAB) *
      exp(cl_pop)

    vc <- exp(lvc + etalvc) *
      (LBM / 55)^e_lbm_vc *
      exp(e_sexf_vc * SEXF) *
      exp(e_ped_vc  * ped) *
      exp(e_ado_vc  * ado)

    # Q inherits the body-weight exponent of CL and VP the lean-body-mass
    # exponent of VC (Hu 2024 Figure 1 caption: "estimates of which were fixed
    # to be similar to those of CL and VC, respectively").
    q  <- exp(lq  + etalq)  * (WT  / 75)^e_wt_cl
    vp <- exp(lvp + etalvp) * (LBM / 55)^e_lbm_vc

    # 3. Time-varying clearance. Adult CNS-tumor (glioblastoma) subjects have
    # stationary CL: the source control stream zeroes the whole Emax term,
    # including its eta and covariate offsets (Appendix File S2:
    # "IF (POP_I .EQ. 4) EMAX = 0").
    emax_stationary <- TUMTP_CNS_PRIM * adult

    cl_time_max_i <- (1 - emax_stationary) *
      (cl_time_max +
         e_ecog_cl_time_max    * ECOG_GE1 +
         e_ipico_cl_time_max   * CONMED_IPI_ANY +
         e_hl_cl_time_max      * TUMTP_LYMPH    * adult +
         e_oth_cl_time_max     * TUMTP_OTHER    * adult +
         e_pedcnst_cl_time_max * TUMTP_CNS_PRIM * child +
         etacl_time_max)

    cl_t50       <- exp(lcl_t50)
    cl_time_hill <- exp(lcl_time_hill)
    cl <- cl0 * exp(cl_time_max_i * t^cl_time_hill / (cl_t50^cl_time_hill + t^cl_time_hill))

    # 4. Two-compartment micro-constants and ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # 5. Observation. Dose in mg and volumes in L give mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
