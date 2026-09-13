Wang_2026_serplulimab <- function() {
  description <- "Two-compartment population PK model with sigmoidal time-varying clearance for intravenous serplulimab (anti-PD-1 IgG4) in adults with advanced solid tumours, with tumour-type-specific baseline clearance and central volume estimated across seven histologies pooled from eleven Phase I-III trials (Wang 2026)"
  reference <- "Wang K, Shen Y, Hu C, Xu F, Kwok Z, Wang Q, Lin Y, Gao Y, Zhou L. Population Pharmacokinetics of Serplulimab and Quantitative Assessment of Transitioning From Weight-Based to Flat-Dosing Strategy. CPT Pharmacometrics Syst Pharmacol. 2026;15:e70204. doi:10.1002/psp4.70204. Updates and supersedes the eight-trial analysis of Wang et al. (2025) Clin Transl Sci 18(9):e70322; see modellib('Wang_2025_serplulimab')."
  vignette <- "Wang_2026_serplulimab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Wang 2026 Section 2.2 (a two-compartment
  # model with time-varying clearance) and Section 3.1 ("14,687 serplulimab
  # serum concentration measurements from 2110 subjects").
  compartmentData <- list(
    central     = list(analyte = "serplulimab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "serplulimab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on baseline CL (exponent 0.514) and on Vc (exponent 0.470); no weight effect on Vp was retained. Reference 62 kg, taken verbatim from the printed final-model equations in Wang 2026 Section 3.1 ('WT/62'); it equals the PK-dataset median of 62.0 kg in Table 1. Weight is the single most influential covariate: the Figure 1 sensitivity forest plot reports Cavg,ss changes of -14.66% at 85 kg and +16.03% at 46 kg relative to the 62 kg reference under a flat 300 mg Q3W regimen, i.e. exposure DECREASES with increasing weight under flat dosing.",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on baseline CL (exponent -0.714), on Vc (exponent -0.320) and on Vp (exponent -1.05). Reference 41.4 g/L, taken verbatim from the printed final-model equations in Wang 2026 Section 3.1 ('ALB/41.4'), which equals the PK-dataset median in Table 1 exactly. Source paper reports albumin in g/L (SI convention), matching the canonical unit. Higher albumin lowers CL and therefore raises exposure (Figure 1: Cavg,ss -15.42% at 32.8 g/L and +10.68% at 47.7 g/L). Wang 2026 Section 4 attributes the weight and albumin effects to FcRn-mediated recycling.",
      source_name        = "ALB"
    ),
    ALP = list(
      description        = "Baseline alkaline phosphatase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on baseline CL only (exponent 0.0553). Reference 94 U/L, taken verbatim from the printed final-model equation in Wang 2026 Section 3.1 ('ALP/94'), equal to the PK-dataset median in Table 1. Higher ALP raises CL and lowers exposure, but the effect is the weakest retained covariate: Figure 1 reports Cavg,ss -4.75% at 236.06 U/L and +2.66% at 57 U/L, and Wang 2026 Section 3.2 calls the maximum change 'not exceeding 3%' for the interquartile span. New to the 2026 analysis; ALP was not retained in Wang 2025.",
      source_name        = "ALP"
    ),
    TUM_SLD = list(
      description        = "Baseline tumour burden, expressed as the sum of longest diameters of target lesions (RECIST)",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on baseline CL (exponent 0.0548) and on Vp (exponent 0.107). Reference 73 mm, taken verbatim from the printed final-model equations in Wang 2026 Section 3.1 ('TUMBUR/73'); the PK-dataset median in Table 1 is 72.9 mm, so 73 is the rounded population median. Higher tumour burden raises CL and lowers exposure (Figure 1: Cavg,ss -4.8% at 177.31 mm and +7.96% at 18.21 mm). Wang 2026 Section 4 reads this as disease-related proteolysis and TMDD-like target-mediated loss. The source column is named TUMBUR and is reported in mm, i.e. a RECIST sum-of-longest-diameters length, which is the TUM_SLD canonical rather than the volumetric TUM_VOL. New to the 2026 analysis; tumour burden was screened but not retained in Wang 2025.",
      source_name        = "TUMBUR"
    ),
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Exponential effect on baseline CL (coefficient -0.145) and on Vc (coefficient -0.14). Wang 2026 Section 3.1 defines 'SEX represents gender (sex = 0 for male, sex = 1 for female)', which matches canonical SEXF directly with no value transformation. Female subjects therefore have exp(-0.145) = 0.865 of the male CL and exp(-0.14) = 0.869 of the male Vc, giving higher exposure in females (Figure 1: Cavg,ss +15.09%, Cmax,ss +15.86%, Cmin,ss +17.76%). The PK dataset was 20.62% female (435/2110; Table 1). Wang 2025 retained a sex effect on Vc only; the 2026 analysis adds the CL effect.",
      source_name        = "SEX"
    ),
    TUMTP_NSCLC_SQUAM = list(
      description        = "Squamous non-small-cell-lung-cancer tumour-type indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any other tumour type in this model; when all six TUMTP_* indicators are 0 the subject is non-squamous NSCLC, which is the model's reference histology)",
      notes              = "Selects the Sq. NSCLC baseline clearance (0.204 L/day) and central volume (3.38 L) from Wang 2026 Table 2. 441/2110 (20.90%) of the PK dataset. Wang 2026 Section 2.2 states that 'separate parameter estimates were obtained for each relevant tumor type rather than pooling them', so tumour type enters as a per-stratum CL0 and Vc rather than as a covariate coefficient.",
      source_name        = "TUMTP"
    ),
    TUMTP_HCC = list(
      description        = "Hepatocellular carcinoma tumour-type indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any other tumour type in this model; all six indicators 0 = non-squamous NSCLC reference)",
      notes              = "Selects the HCC baseline clearance (0.204 L/day) and central volume (3.20 L) from Wang 2026 Table 2. 125/2110 (5.92%) of the PK dataset. HCC shares the Sq. NSCLC CL0 to three significant figures but has its own Vc, so the two histologies are not interchangeable.",
      source_name        = "TUMTP"
    ),
    TUMTP_CRC = list(
      description        = "Colorectal cancer tumour-type indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any other tumour type in this model; all six indicators 0 = non-squamous NSCLC reference)",
      notes              = "Selects the CRC baseline clearance (0.182 L/day) and central volume (3.19 L) from Wang 2026 Table 2. 151/2110 (7.16%) of the PK dataset, contributed largely by the HLX10-015-mCRC301 Phase II/III trial newly added in the 2026 analysis.",
      source_name        = "TUMTP"
    ),
    TUMTP_SCLC = list(
      description        = "Small cell lung cancer tumour-type indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any other tumour type in this model; all six indicators 0 = non-squamous NSCLC reference)",
      notes              = "Selects the SCLC baseline clearance (0.171 L/day) and central volume (3.45 L) from Wang 2026 Table 2. 390/2110 (18.48%) of the PK dataset, from the Phase III ASTRUM-005 trial (HLX10-005-SCLC301 / NCT04063163). SCLC has the LOWEST clearance of the seven histologies and therefore the highest exposure (Figure 1: Cavg,ss +7.17%, Cmin,ss +11.77% vs the non-squamous NSCLC reference).",
      source_name        = "TUMTP"
    ),
    TUMTP_ESCC = list(
      description        = "Oesophageal squamous cell carcinoma tumour-type indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any other tumour type in this model; all six indicators 0 = non-squamous NSCLC reference)",
      notes              = "Selects the ESCC baseline clearance (0.178 L/day) and central volume (3.48 L) from Wang 2026 Table 2. 389/2110 (18.44%) of the PK dataset, from the HLX10-007-EC301 Phase III trial newly added in the 2026 analysis. ESCC has the largest central volume of the seven histologies. Wang 2026 Figure 1 labels this row 'EC'.",
      source_name        = "TUMTP"
    ),
    TUMTP_OTHER = list(
      description        = "Residual 'other tumour types' indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any named histology in this model; all six indicators 0 = non-squamous NSCLC reference)",
      notes              = "Selects the 'other' baseline clearance (0.211 L/day) and central volume (3.19 L) from Wang 2026 Table 2. 116/2110 (5.50%) of the PK dataset: the residual bucket complementing Sq. NSCLC, non-sq. NSCLC, SCLC, ESCC, HCC and CRC, populated mainly by the MSI-H/dMMR solid-tumour, cervical-cancer and head-and-neck cohorts of the Phase I/II trials (Table S1). This group has the HIGHEST clearance and therefore the lowest exposure (Figure 1: Cavg,ss -12.17% vs the non-squamous NSCLC reference).",
      source_name        = "TUMTP"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 2110L,
    n_studies      = 11L,
    n_observations = 14687L,
    age_range      = "23.0-83.0 years",
    age_median     = "61.0 years",
    weight_range   = "32.9-131 kg",
    weight_median  = "62.0 kg",
    height_median  = "167 cm (128-191 cm)",
    bmi_median     = "22.6 kg/m^2 (13.3-42.3 kg/m^2)",
    bsa_median     = "1.70 m^2 (1.19-2.55 m^2)",
    sex_female_pct = 20.62,
    race_ethnicity = c(Asian = 87.44, `Non-Asian` = 12.56),
    disease_state  = "Adults with advanced solid tumours. Tumour-type mix in the PK dataset (n = 2110; Table 1): non-squamous NSCLC 498 (23.60%), squamous NSCLC 441 (20.90%), SCLC 390 (18.48%), ESCC 389 (18.44%), colorectal cancer 151 (7.16%), hepatocellular carcinoma 125 (5.92%), other tumour types 116 (5.50%).",
    dose_range     = "Serplulimab 0.3-10 mg/kg IV plus flat doses of 200 mg Q2W, 300 mg Q3W and 400 mg Q4W across the pooled trials (Table S1). The weight-based regimens taken forward were 3 mg/kg Q2W and 4.5 mg/kg Q3W; the flat regimens evaluated for the dosing transition were 200 mg Q2W and 300 mg Q3W, each as a 1-h infusion.",
    regions        = "Predominantly China (Asian 87.44%); the ASTRUM-004 / ASTRUM-005 Phase III programmes were multinational and contribute the 12.56% non-Asian subjects.",
    ada_status     = "ADA-negative 1977 (93.70%); ADA-positive 133 (6.30%). Wang 2026 Section 3.2 reports no more than a 3% reduction in steady-state exposure among ADA-positive patients; ADA was not retained in the final model.",
    ecog_status    = "ECOG performance status 0 in 567 (26.87%), 1 in 1539 (72.94%), 2 in 3, missing in 1.",
    albumin        = "41.4 g/L median (23.8-67.9 g/L range).",
    tumour_burden  = "72.9 mm median (0.00-350 mm range), as the RECIST sum of longest target-lesion diameters.",
    renal_function = "Creatinine clearance 87.3 mL/min median (26.4-291 mL/min, Cockcroft-Gault, not BSA-normalised); serum creatinine 68.7 umol/L median (23.0-162 umol/L). Neither was retained. Wang 2026 Section 3.2 reports 8.92%-21.5% HIGHER exposure in mild or moderate renal impairment than in normal renal function, judged not clinically meaningful.",
    hepatic_function = "Wang 2026 Section 3.2 reports 1.06%-6.64% lower exposure with mild hepatic impairment than with normal hepatic function; the moderate-impairment subgroup (n = 4) was too small to interpret.",
    concomitant_therapy = "Concomitant chemotherapy in 1640 (77.73%); combination antibody-based anti-tumour therapy in 482 (22.84%).",
    notes          = "Baseline demographics per Wang 2026 Table 1 (PK dataset column). Pooled dataset spans eleven serplulimab (HLX10) trials listed in Table S1 with their per-trial N: HLX10-001 / NCT03468751 (Phase I, n = 57), HLX10HLX04-001 / NCT03757936 (Phase I, n = 26), HLX10-008-HCC201 / NCT03973112 (Phase II, n = 123), HLX10-010-MSI201 / NCT03941574 (Phase II, n = 108), HLX10-011-CC201 / NCT04150575 (Phase II, n = 21), HLX10HLX07-001 / NCT04297995 (Phase II, n = 13), HLX10-002-NSCLC301 / NCT03952403 (Phase III, n = 491), HLX10-004-NSCLC303 / NCT04033354 (Phase III, n = 439), HLX10-005-SCLC301 / NCT04063163 (Phase III, n = 389), HLX10-007-EC301 / NCT03958890 (Phase III, n = 389), HLX10-015-mCRC301 / NCT04547166 (Phase II/III, n = 64). The last three are the trials newly added relative to the eight-trial Wang 2025 analysis; the other eight carry updated cut-off dates. 14,687 serum concentrations from 2110 subjects entered the final PopPK dataset. Fit in NONMEM with stepwise forward inclusion (p < 0.01) and backward elimination (p < 0.001); no covariate was removed during backward elimination. Model performance assessed by bootstrap resampling and prediction-corrected VPC stratified by tumour type."
  )

  # Covariates collected and screened by Wang 2026 but NOT retained in the final
  # model. Documented here to preserve the provenance of the covariate screen
  # without triggering a "declared but not referenced" convention warning.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age",
      units       = "years",
      type        = "continuous",
      notes       = "Collected (Table 1: median 61.0 y, range 23.0-83.0 y) and not retained. Wang 2026 Section 3.2 reports an exploratory EBE-based subgroup contrast of only -0.266% to +1.14% between patients <= 65 y and > 65 y, so no coefficient exists to encode."
    ),
    HT = list(
      description = "Baseline height",
      units       = "cm",
      type        = "continuous",
      notes       = "Collected (Table 1: median 167 cm) and not retained; body weight was the retained body-size descriptor."
    ),
    BMI = list(
      description = "Baseline body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Collected (Table 1: median 22.6 kg/m^2) and not retained; correlated with the retained WT."
    ),
    BSA = list(
      description = "Baseline body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Collected (Table 1: median 1.70 m^2) and not retained; correlated with the retained WT."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Not involved in covariate screening. Wang 2026 Section 3.2 reports an exploratory EBE-based contrast of 10.1%-11.8% LOWER exposure in White than in Asian patients (Figure 2), judged not clinically meaningful. No model coefficient is published, so the effect cannot be encoded."
    ),
    ADA_POS = list(
      description = "Anti-drug-antibody positivity status",
      units       = "(binary)",
      type        = "binary",
      notes       = "Not retained. 133/2110 (6.30%) ADA-positive; Wang 2026 Section 3.2 reports no more than a 3% exposure reduction in ADA-positive patients (Figure 2). No model coefficient is published."
    ),
    LDH = list(
      description = "Baseline lactate dehydrogenase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Collected (Table 1: median 208 U/L) and not retained in the final PK model."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Collected (Table 1: median 21.0 U/L) and not retained. Note that the related ALP WAS retained, on CL."
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Collected (Table 1: median 18.0 U/L) and not retained."
    ),
    TBILI = list(
      description = "Baseline total bilirubin (source column BILI)",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Collected (Table 1: median 10.3 umol/L) and not retained. Reported in SI units (umol/L), matching the canonical."
    ),
    CREAT = list(
      description = "Baseline serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Collected (Table 1: median 68.7 umol/L) and not retained."
    ),
    CRCL = list(
      description = "Baseline creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Collected (Table 1: median 87.3 mL/min) and not retained. Reported as an absolute Cockcroft-Gault clearance in mL/min, NOT the BSA-normalised mL/min/1.73 m^2 of the canonical CRCL entry. Consistent with a 148 kDa IgG4 mAb not being renally cleared."
    ),
    ECOG_GE1 = list(
      description = "Eastern Cooperative Oncology Group performance status >= 1 (source column ECOG, levels 0 / 1 / 2)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Collected (Table 1: ECOG 0 in 26.87%, 1 in 72.94%, 2 in 3 subjects) and not retained in the final PK model."
    ),
    CONMED_CHEMO = list(
      description = "Concomitant chemotherapy indicator (source column COMB)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Collected (Table 1: 1640/2110, 77.73%) and not retained in the final PK model."
    )
  )

  ini({
    # ---- Tumour-type-specific baseline clearance (Wang 2026 Table 2) --------
    # Wang 2026 Section 2.2: "separate parameter estimates were obtained for
    # each relevant tumor type rather than pooling them for the convenience of
    # application." Seven parallel CL0 estimates and seven parallel Vc
    # estimates come out of ONE joint fit, so every stratum carries an explicit
    # suffix and none keeps the bare canonical name (parameter-names.md
    # 'Stratum-suffixed parameters').
    #
    # UNITS. The Section 3.1 equation prints CL0 in L/day with a trailing
    # 'x 24', i.e. the equation's theta is an L/h value and 24 converts it to
    # L/day; the same 'x 24' appears on Q, where 0.0169 x 24 = 0.4056 reproduces
    # the 0.405 L/day of Table 2 exactly. The values below are therefore the
    # already-converted L/day numbers of Table 2, NOT Table 2's numbers times
    # 24 again. Three independent checks confirm this reading:
    #   (a) terminal half-life. With CL0 0.184 L/day and Vss 6.23 L the
    #       terminal t1/2 is 26 days, the textbook value for an IgG4 mAb;
    #       CL0 x 24 would give 5.6 days.
    #   (b) the Figure 1 forest plot's reference Cavg,ss. Back-solving the two
    #       printed 90% bounds, 51 ug/mL (-35.92%) and 118 ug/mL (+49.65%),
    #       puts the reference typical subject at 79.2 ug/mL; the model gives
    #       300 mg / (0.184 x 0.969 L/day x 21 day) = 80.1 ug/mL.
    #   (c) single-dose trough. The model gives Cmin1 = 20.5 ug/mL for
    #       4.5 mg/kg Q3W and 16.9 ug/mL for 3 mg/kg Q2W at 62 kg, against
    #       medians of about 21 and 17.5 ug/mL in Figure 3A,B.
    #
    # The reference subject for every covariate term below is WT 62 kg,
    # ALB 41.4 g/L, TUMBUR 73 mm, ALP 94 U/L, male, non-squamous NSCLC --
    # the rounded PK-dataset medians of Table 1 and the histology Figure 1
    # labels '+0%'.
    lcl_sqnsclc    <- log(0.204); label("Baseline clearance CL0 at reference covariates, squamous NSCLC (L/day)")     # Wang 2026 Table 2: 0.204 (RSE 2.04%, 95% CI 0.196-0.212); bootstrap median 0.204 (0.195-0.215)
    lcl_hcc        <- log(0.204); label("Baseline clearance CL0 at reference covariates, hepatocellular carcinoma (L/day)") # Wang 2026 Table 2: 0.204 (RSE 2.37%, 95% CI 0.195-0.214); bootstrap median 0.205 (0.196-0.216)
    lcl_crc        <- log(0.182); label("Baseline clearance CL0 at reference covariates, colorectal cancer (L/day)")  # Wang 2026 Table 2: 0.182 (RSE 2.84%, 95% CI 0.172-0.192); bootstrap median 0.182 (0.172-0.195)
    lcl_nonsqnsclc <- log(0.184); label("Baseline clearance CL0 at reference covariates, non-squamous NSCLC (L/day)") # Wang 2026 Table 2: 0.184 (RSE 1.54%, 95% CI 0.179-0.190); bootstrap median 0.185 (0.178-0.192)
    lcl_sclc       <- log(0.171); label("Baseline clearance CL0 at reference covariates, small cell lung cancer (L/day)") # Wang 2026 Table 2: 0.171 (RSE 1.8%, 95% CI 0.165-0.177); bootstrap median 0.171 (0.165-0.179)
    lcl_escc       <- log(0.178); label("Baseline clearance CL0 at reference covariates, oesophageal squamous cell carcinoma (L/day)") # Wang 2026 Table 2: 0.178 (RSE 1.98%, 95% CI 0.171-0.185); bootstrap median 0.179 (0.172-0.187)
    lcl_other      <- log(0.211); label("Baseline clearance CL0 at reference covariates, other tumour types (L/day)") # Wang 2026 Table 2: 0.211 (RSE 3.18%, 95% CI 0.198-0.224); bootstrap median 0.211 (0.197-0.226)

    # ---- Tumour-type-specific central volume (Wang 2026 Table 2) -----------
    lvc_sqnsclc    <- log(3.38); label("Central volume Vc at reference covariates, squamous NSCLC (L)")     # Wang 2026 Table 2: 3.38 (RSE 1.1%, 95% CI 3.31-3.45); bootstrap median 3.38 (3.31-3.46)
    lvc_hcc        <- log(3.20); label("Central volume Vc at reference covariates, hepatocellular carcinoma (L)") # Wang 2026 Table 2: 3.20 (RSE 2.02%, 95% CI 3.08-3.33); bootstrap median 3.21 (3.08-3.33)
    lvc_crc        <- log(3.19); label("Central volume Vc at reference covariates, colorectal cancer (L)")  # Wang 2026 Table 2: 3.19 (RSE 1.6%, 95% CI 3.09-3.29); bootstrap median 3.18 (3.08-3.28)
    lvc_nonsqnsclc <- log(3.25); label("Central volume Vc at reference covariates, non-squamous NSCLC (L)") # Wang 2026 Table 2: 3.25 (RSE 0.961%, 95% CI 3.18-3.31); bootstrap median 3.24 (3.18-3.31)
    lvc_sclc       <- log(3.45); label("Central volume Vc at reference covariates, small cell lung cancer (L)") # Wang 2026 Table 2: 3.45 (RSE 1.22%, 95% CI 3.37-3.54); bootstrap median 3.45 (3.37-3.54)
    lvc_escc       <- log(3.48); label("Central volume Vc at reference covariates, oesophageal squamous cell carcinoma (L)") # Wang 2026 Table 2: 3.48 (RSE 1.51%, 95% CI 3.38-3.59); bootstrap median 3.48 (3.37-3.59)
    lvc_other      <- log(3.19); label("Central volume Vc at reference covariates, other tumour types (L)") # Wang 2026 Table 2: 3.19 (RSE 1.68%, 95% CI 3.09-3.30); bootstrap median 3.19 (3.08-3.31)

    # ---- Tumour-type-independent disposition (Wang 2026 Table 2) -----------
    lq  <- log(0.405); label("Intercompartmental clearance Q (L/day)")                # Wang 2026 Table 2: Q = 0.405 L/day (RSE 6.1%, 95% CI 0.359-0.456). The Section 3.1 equation prints the same value per hour, 'Q_i = 0.0169 x exp(eta_Q) x 24', and 0.0169 x 24 = 0.4056 reproduces it.
    lvp <- log(2.98);  label("Peripheral volume Vp at reference covariates (L)")      # Wang 2026 Table 2: Vp = 2.98 L (RSE 2.77%, 95% CI 2.82-3.14); the same 2.98 is the printed leading constant of the Section 3.1 Vp equation

    # ---- Time-varying clearance (Wang 2026 Section 3.1) --------------------
    # Printed form:
    #   CL_i = CL0_i * exp[ Emax_i * (time/T50)^lambda / (1 + (time/T50)^lambda) ]
    # Multiplying numerator and denominator by T50^lambda gives the equivalent
    #   CL_i = CL0_i * exp( Emax_i * t^lambda / (T50^lambda + t^lambda) )
    # used in model() below. The bracket rises from 0 at t = 0 to 1 as
    # t >> T50, so CL falls from CL0 to CL0 * exp(Emax). With Emax = -0.0926
    # the asymptote is exp(-0.0926) = 0.9115 of baseline, i.e. an 8.8%
    # reduction at full saturation -- which is exactly the
    # "Maximum change ratio in clearance, exp(Emax) = 0.912" row of Table 2,
    # an independent confirmation of the sign and scale. Wang 2026 Section 2.2
    # is explicit that this is an EMPIRICAL descriptor of disease-related
    # change, adopted because mechanistic TMDD forms "did not improve model
    # performance" given rapid PD-1 receptor saturation.
    cl_time_max   <- -0.0926;  label("Maximum log-scale change in CL from baseline (Emax; unitless, negative = CL decreases over time)") # Wang 2026 Section 3.1: Emax_i = -0.0926 + eta_Emax,i; Table 2 reports the back-transform exp(Emax) = 0.912 (RSE 2.11%, 95% CI 0.875-0.95)
    lcl_t50       <- log(221); label("log T50 - time at which half of the maximum CL change is reached (log days)")                     # Wang 2026 Table 2: T50 = 221 day (RSE 12.0%, 95% CI 169-273); the same 221 is printed in the Section 3.1 equation block
    lcl_time_hill <- log(2.43); label("log lambda - sigmoidicity (Hill coefficient) of the time-on-CL function (log unitless)")          # Wang 2026 Table 2: lambda = 2.43 (RSE 6.21%, 95% CI 2.13-2.73); the same 2.43 is printed in the Section 3.1 equation block

    # ---- Covariate effects (Wang 2026 Table 2 + Section 3.1 equations) -----
    # Continuous covariates are printed as power terms on the covariate/reference
    # ratio, e.g. (WT/62)^0.514, and are written in that same power form below.
    e_wt_cl      <-  0.514;  label("Power exponent of body weight on baseline CL (unitless)")     # Wang 2026 Table 2 CLwt = 0.514 (RSE 6.46%, 95% CI 0.449-0.579); Section 3.1 equation: (WT/62)^0.514
    e_alb_cl     <- -0.714;  label("Power exponent of albumin on baseline CL (unitless)")         # Wang 2026 Table 2 CLalb = -0.714 (RSE 9.39%, 95% CI -0.845 to -0.582); Section 3.1 equation: (ALB/41.4)^-0.714
    e_tum_sld_cl <-  0.0548; label("Power exponent of baseline tumour burden on baseline CL (unitless)") # Wang 2026 Table 2 CLtmb = 0.0548 (RSE 20.8%, 95% CI 0.0325-0.0771); Section 3.1 equation: (TUMBUR/73)^0.0548
    e_alp_cl     <-  0.0553; label("Power exponent of alkaline phosphatase on baseline CL (unitless)")   # Wang 2026 Table 2 CLalp = 0.0553 (RSE 29.9%, 95% CI 0.0229-0.0877); Section 3.1 equation: (ALP/94)^0.0553
    e_wt_vc      <-  0.470;  label("Power exponent of body weight on Vc (unitless)")              # Wang 2026 Table 2 Vwt = 0.47 (RSE 5.81%, 95% CI 0.416-0.523); Section 3.1 equation: (WT/62)^0.470
    e_alb_vc     <- -0.320;  label("Power exponent of albumin on Vc (unitless)")                  # Wang 2026 Table 2 Vcalb = -0.32 (RSE 14.2%, 95% CI -0.409 to -0.231); Section 3.1 equation: (ALB/41.4)^-0.320
    e_alb_vp     <- -1.05;   label("Power exponent of albumin on Vp (unitless)")                  # Wang 2026 Table 2 Vpalb = -1.05 (RSE 15.5%, 95% CI -1.37 to -0.732); Section 3.1 equation: (ALB/41.4)^-1.05
    e_tum_sld_vp <-  0.107;  label("Power exponent of baseline tumour burden on Vp (unitless)")   # Wang 2026 Table 2 VPtmb = 0.107 (RSE 24.6%, 95% CI 0.0554-0.158); Section 3.1 equation: (TUMBUR/73)^0.107

    # Categorical covariates enter as exp(coefficient * indicator).
    e_sexf_cl <- -0.145; label("Exponential coefficient of female sex on baseline CL (unitless)") # Wang 2026 Table 2 CLsex = -0.145 (RSE 10.8%, 95% CI -0.175 to -0.114); Section 3.1 equation: exp[-0.145 * (Female)] with Female = 1 for female
    e_sexf_vc <- -0.14;  label("Exponential coefficient of female sex on Vc (unitless)")          # Wang 2026 Table 2 Vcsex = -0.14 (RSE 8.6%, 95% CI -0.163 to -0.116); Section 3.1 equation: exp[-0.14 * (Female)] with Female = 1 for female

    # ---- Inter-individual variability (Wang 2026 Table 2) ------------------
    # Table 2 note: "IIV for CL, Vc, Q, Vp, Emax, and residual are reported as
    # approximate CV%". Each reported percentage is omega itself (the SD on the
    # estimation scale) x 100, NOT the variance. Two independent arguments
    # settle this.
    #
    # (1) Every published 95% CI is symmetric on the VARIANCE scale, which is
    #     what a NONMEM covariance step delivers for an OMEGA element:
    #       CL    24.0% -> 0.240^2 = 0.05760; CI (22.6, 25.3) -> (0.05108, 0.06401), midpoint 0.05754
    #       Vc    16.3% -> 0.163^2 = 0.02657; CI (15.0, 17.5) -> (0.02250, 0.03063), midpoint 0.02656
    #       Q     54.3% -> 0.543^2 = 0.29485; CI (43.7, 63.1) -> (0.19097, 0.39816), midpoint 0.29457
    #       Vp    45.9% -> 0.459^2 = 0.21068; CI (41.5, 49.8) -> (0.17223, 0.24800), midpoint 0.21011
    #       Emax  34.1% -> 0.341^2 = 0.11628; CI (27.3, 39.7) -> (0.07453, 0.15761), midpoint 0.11607
    #       sigma 17.6% -> 0.176^2 = 0.03098; CI (17.0, 18.2) -> (0.02890, 0.03312), midpoint 0.03101
    #     A variance reading would put every midpoint far from the point
    #     estimate; the SD reading is the only one consistent with the CIs.
    #
    # (2) The CL-Vc covariance is reported on the raw OMEGA-block scale as
    #     0.014 (its own 95% CI 0.011-0.017 is symmetric about it). Pairing
    #     0.014 with variances of 24.0 and 16.3 would imply a correlation of
    #     0.0007; pairing it with SDs of 0.240 and 0.163 gives 0.358, which is
    #     the only physically sensible reading.
    #
    # CL, Vc, Q and Vp carry log-normal etas. Emax carries an ADDITIVE eta on
    # the linear scale, per the printed Section 3.1 relation
    # Emax_i = -0.0926 + eta_Emax,i. Its omega is still 0.341 and NOT
    # 0.341 x |Emax|: the reporting routine emitted sqrt(omega^2) x 100
    # uniformly for every row of the table, and the resulting spread
    # (exp(Emax_i) roughly 0.52-1.60 over a 90% interval, i.e. some subjects'
    # clearance rising and others' falling over treatment) is the documented
    # behaviour for checkpoint inhibitors that this empirical term exists to
    # capture. ETA shrinkage was 14.9% (CL), 27.0% (Vc), 67.3% (Q), 35.8% (Vp)
    # and 46.7% (Emax).
    etalcl + etalvc ~ c(0.05760,
                        0.014, 0.02657)  # Wang 2026 Table 2: IIV CL 24.0%, IIV Vc 16.3%, Covariance (CL_Vc) 0.014 (RSE 11%, 95% CI 0.011-0.017)
    etalq          ~ 0.29485             # Wang 2026 Table 2: IIV Q 54.3% (RSE 8.98%, 95% CI 43.7-63.1)
    etalvp         ~ 0.21068             # Wang 2026 Table 2: IIV Vp 45.9% (RSE 4.6%, 95% CI 41.5-49.8)
    etacl_time_max ~ 0.11628             # Wang 2026 Table 2: IIV Emax 34.1% (RSE 9.11%, 95% CI 27.3-39.7); additive eta on the linear-scale Emax

    # ---- Residual error (Wang 2026 Table 2) --------------------------------
    # A single residual term of 17.6% is reported, described in the Table 2
    # note as an approximate CV%. A constant-CV residual on a concentration is
    # proportional error in nlmixr2's linear space.
    propSd <- 0.176; label("Proportional residual error (fraction)")  # Wang 2026 Table 2: Residual error sigma = 17.6% (RSE 1.75%, 95% CI 17-18.2); bootstrap median 17.6 (17-18.2)
  })
  model({
    # 1. Tumour-type selection. Exactly one histology applies to any subject.
    #    Non-squamous NSCLC is the complement of the other six indicators --
    #    it is the largest group (498/2110) and the histology Wang 2026
    #    Figure 1 labels "+0%" for every tumour-type contrast -- so a subject
    #    with all six indicators at 0 is a non-squamous-NSCLC subject and gets
    #    that stratum's CL0 and Vc.
    tumtpNonsq <- 1 - TUMTP_NSCLC_SQUAM - TUMTP_HCC - TUMTP_CRC -
      TUMTP_SCLC - TUMTP_ESCC - TUMTP_OTHER

    lcl0 <- lcl_nonsqnsclc * tumtpNonsq +
      lcl_sqnsclc * TUMTP_NSCLC_SQUAM +
      lcl_hcc     * TUMTP_HCC +
      lcl_crc     * TUMTP_CRC +
      lcl_sclc    * TUMTP_SCLC +
      lcl_escc    * TUMTP_ESCC +
      lcl_other   * TUMTP_OTHER

    lvc0 <- lvc_nonsqnsclc * tumtpNonsq +
      lvc_sqnsclc * TUMTP_NSCLC_SQUAM +
      lvc_hcc     * TUMTP_HCC +
      lvc_crc     * TUMTP_CRC +
      lvc_sclc    * TUMTP_SCLC +
      lvc_escc    * TUMTP_ESCC +
      lvc_other   * TUMTP_OTHER

    # 2. Individual baseline CL (CL at t = 0), Vc, Q and Vp, in the power form
    #    printed in Wang 2026 Section 3.1.
    cl0 <- exp(lcl0 + etalcl) *
      (WT      / 62)^e_wt_cl *
      (ALB     / 41.4)^e_alb_cl *
      (TUM_SLD / 73)^e_tum_sld_cl *
      (ALP     / 94)^e_alp_cl *
      exp(e_sexf_cl * SEXF)

    vc <- exp(lvc0 + etalvc) *
      (WT  / 62)^e_wt_vc *
      (ALB / 41.4)^e_alb_vc *
      exp(e_sexf_vc * SEXF)

    q <- exp(lq + etalq)

    vp <- exp(lvp + etalvp) *
      (ALB     / 41.4)^e_alb_vp *
      (TUM_SLD / 73)^e_tum_sld_vp

    # 3. Time-varying CL: sigmoidal in time since the first dose.
    #    t^lambda / (T50^lambda + t^lambda) rises from 0 at t = 0 to 1 as
    #    t >> T50, so CL falls from cl0 to cl0 * exp(Emax_i).
    cl_t50        <- exp(lcl_t50)
    cl_time_hill  <- exp(lcl_time_hill)
    cl_time_max_i <- cl_time_max + etacl_time_max
    cl <- cl0 * exp(cl_time_max_i * t^cl_time_hill / (cl_t50^cl_time_hill + t^cl_time_hill))

    # 4. Two-compartment micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # 5. Dose in mg, volumes in L => central / vc has units mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
