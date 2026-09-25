Li_2019_ceftazidime <- function() {
  description <- paste(
    "Two-compartment IV population PK model for the ceftazidime component of",
    "ceftazidime-avibactam, fitted to 9,155 plasma concentrations from 1,975",
    "adults pooled across 11 Phase 1 studies, 2 Phase 2 studies and 5 Phase 3",
    "trials in complicated intra-abdominal infection (cIAI), complicated",
    "urinary tract infection (cUTI) and nosocomial pneumonia (NP), including",
    "ventilator-associated pneumonia (Li 2019). Creatinine clearance is the",
    "key covariate on clearance, entering as a hinged LINEAR relationship",
    "that is steep below 100 mL/min and shallow above it. Central volume",
    "carries body weight, infection type, acute pyelonephritis, Asian race",
    "and mechanical-ventilation effects. Inter-individual variability is a",
    "full 4x4 OMEGA block across CL, Vc, Vp and Q. Residual variability is",
    "switched between a Phase 1 and a pooled Phase 2/3 stratum. Companion to",
    "Li_2019_avibactam; the two analytes were fitted as separate models to",
    "separate data sets and are combined only in the joint PK/PD target",
    "attainment analysis."
  )
  reference <- paste(
    "Li J, Lovern M, Green ML, Chiu J, Zhou D, Comisar C, Xiong Y, Hing J,",
    "MacPherson M, Wright JG, Riccobene T, Carrothers TJ, Das S.",
    "Ceftazidime-Avibactam Population Pharmacokinetic Modeling and",
    "Pharmacodynamic Target Attainment Across Adult Indications and Patient",
    "Subgroups. Clin Transl Sci. 2019;12(2):151-163. doi:10.1111/cts.12585.",
    "Fixed-effect, random-effect and residual-error estimates are Table 1.",
    "The functional FORMS of every covariate relationship, and the",
    "hard-coded creatinine-clearance slopes, are the final ceftazidime",
    "NONMEM control stream printed in supplementary Data S1",
    "(CTS-12-151-s001, 'NONMEM control file for ceftazidime final PopPK",
    "model').",
    sep = " "
  )
  vignette <- "Li_2019_ceftazidime_avibactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central = list(analyte = "ceftazidime", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ceftazidime", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated with the Cockcroft-Gault equation from chronological serum creatinine records (Li 2019 Methods, 'Analysis data and model construction')",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Raw Cockcroft-Gault mL/min, NOT BSA-normalized. Stored under the canonical CRCL column, which accepts raw mL/min when the source paper applies no BSA normalization (same convention as the sibling combination-product pair Chen_2025_ceftazidime.R / Chen_2025_avibactam.R and as Chandorkar_2015_ceftolozane.R). Observed range in the ceftazidime data set 8-488 mL/min (Li 2019 Methods). The effect on CL is a two-segment LINEAR spline hinged at 100 mL/min, NOT the power-then-linear hinge used by the avibactam companion model: below the hinge the factor is slope1 * CRCL and above it slope1 * 100 + slope2 * (CRCL - 100). Because the lower segment passes through the ORIGIN rather than through 1 at the hinge, the factor equals 1 at CRCL = 1 / 0.0103036 = 97.1 mL/min, which is therefore the reference renal function to which the ini() clearance refers. Ceftazidime carries no separate ESRD, dialysis or augmented-renal-clearance term; the spline alone spans the observed range.",
      source_name = "CLCR"
    ),
    WT = list(
      description = "Body weight at baseline (Li 2019 Table 1, theta13 'WT effect on Vc')",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Centred at 70 kg, the cohort median quoted in Li 2019 Results, 'Ceftazidime'; the control stream form is V1WT = (WT/70)^theta13. Weight acts on Vc ONLY -- the final ceftazidime model carries no allometric term on CL, Q or Vp. CAUTION, the Results prose does NOT reproduce from the tabulated exponent: it reports '24% lower and 26% higher Vc' at the 10th (50 kg) and 90th (94 kg) weight percentiles, whereas (50/70)^1.01 = 0.712 and (94/70)^1.01 = 1.347, i.e. 28.8% lower and 34.7% higher. Recovering -24%/+26% from an exponent of 1.01 would need weights of about 53 kg and 88 kg, and recovering them at 50/94 kg would need an exponent near 0.80, which is not the tabulated value. The equation is taken as authoritative over the prose (the two other Vc effects quoted in the same sentence, 27% for Asian race and 29.7% for NPv, are verbatim theta12 and theta14, so the sentence is a mix of exact and approximate figures). The COMPANION avibactam model has no such problem: its (51/70)^1.08 = 0.710 and (95/70)^1.08 = 1.391 reproduce its quoted 29% lower / 39% higher exactly, which is what identifies the ceftazidime sentence as the unreliable one. Recorded in the validation vignette's Assumptions and deviations section.",
      source_name = "WT"
    ),
    DIS_CIAI = list(
      description = "Complicated intra-abdominal infection infection-type indicator (Li 2019 control stream POP = 3)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not cIAI; the shared all-zero reference across the DIS_* set is the Phase 1 healthy volunteer",
      notes = "781 of 1,975 subjects (39.5%), Li 2019 Results, 'Analysis populations'. Acts on BOTH CL (multiplicative factor 1.16, Table 1 theta5) and Vc (multiplicative factor 1.14 shared with NP, Table 1 theta10). Mutually exclusive with DIS_CUTI, DIS_HABP and DIS_VABP.",
      source_name = "POP = 3"
    ),
    DIS_CUTI = list(
      description = "Complicated urinary tract infection infection-type indicator (Li 2019 control stream POP = 2)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not cUTI; shared all-zero reference is the Phase 1 healthy volunteer",
      notes = "696 of 1,975 subjects (35.2%), Li 2019 Results, 'Analysis populations'. Acts on Vc only (multiplicative factor 1.03, Table 1 theta9); the final model carries NO cUTI effect on ceftazidime CL, so cUTI subjects share the healthy-volunteer clearance. Subjects enrolled under the combined 'cUTI including acute pyelonephritis' protocol additionally carry DIS_AP, which stacks a further Vc shift on top of this one.",
      source_name = "POP = 2"
    ),
    DIS_HABP = list(
      description = "Hospital-acquired bacterial pneumonia infection-type indicator (component of Li 2019 control stream POP = 4, 'NP')",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not HABP; shared all-zero reference is the Phase 1 healthy volunteer",
      notes = "Li 2019 models nosocomial pneumonia as a SINGLE level (POP = 4) covering both hospital-acquired and ventilator-associated pneumonia, so DIS_HABP and DIS_VABP always share one estimated coefficient here and enter model() as their sum. The columns are kept distinct per the DIS_VABP register entry ('when a source paper merges the pneumonia arms into one coefficient, keep the covariate COLUMNS distinct and apply the shared coefficient to their sum'), because the successor models in this lineage separate them. 412 of 1,975 subjects (20.9%) had NP of either kind. Whether a given NP subject was ventilated ON THE PK SAMPLING DAY is carried separately by MECH_VENT, which is NOT the same contrast as HABP vs VABP.",
      source_name = "POP = 4"
    ),
    DIS_VABP = list(
      description = "Ventilator-associated bacterial pneumonia infection-type indicator (component of Li 2019 control stream POP = 4, 'NP')",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not VABP; shared all-zero reference is the Phase 1 healthy volunteer",
      notes = "See DIS_HABP: Li 2019 carries one pooled NP coefficient, applied to DIS_HABP + DIS_VABP. Of the 413 NP subjects in the avibactam data set, 138 had VAP and 275 did not (Li 2019 Table 3).",
      source_name = "POP = 4"
    ),
    DIS_AP = list(
      description = "Acute pyelonephritis indicator within the cUTI cohort (Li 2019 control stream POP1 = 1)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = cUTI without acute pyelonephritis, or any non-cUTI subject",
      notes = "Li 2019 Table 1 theta11, 'Population effect on Vc (cUTI/acute pyelonephritis)' = -0.185, i.e. 18.5% LOWER Vc. Applied as a proportional shift that STACKS multiplicatively on top of the DIS_CUTI factor (control stream: V1COV = V1POP * V1POP1 * V1RCE * V1POP5), so an AP subject carries 1.03 * (1 - 0.185) = 0.839 relative to a healthy volunteer. This is one of only two ceftazidime fixed effects with RSE above 27% (41.2%; Li 2019 Results, 'Ceftazidime') and was retained as a covariate of particular clinical interest rather than on the >=20% clinical-relevance criterion.",
      source_name = "POP1 = 1"
    ),
    RACE_ASIAN_OTH = list(
      description = "Non-Chinese, non-Japanese Asian race indicator (Li 2019 control stream RCE = 3, abbreviated ASN)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = White / Black / Other non-Asian (the RCE levels that take no race factor)",
      notes = "The dominant reference grouping is the pooled non-Asian population (Li 2019 Table 3 reports 1,209 White/other subjects against 248 non-Chinese Asian, 262 Chinese/Taiwanese and 45 Japanese). Li 2019 partitions Asian heritage into three mutually exclusive levels, so RACE_ASIAN_OTH is used in preference to the undivided RACE_ASIAN per that entry's Notes. Carries a clearance shift of its own (-0.161, Table 1 theta7 = 16% lower CL) AND a share of the pooled Asian volume shift (-0.27, theta12), which applies to all three Asian levels.",
      source_name = "RCE = 3"
    ),
    RACE_CHINESE = list(
      description = "Chinese-heritage race indicator (Li 2019 control stream RCE = 14, abbreviated CHN)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = non-Chinese",
      notes = "262 Chinese and Taiwanese subjects (Li 2019 Table 3). Carries a clearance shift of -0.0855 (Table 1 theta8, 9% lower CL) AND a share of the pooled Asian volume shift (-0.27, theta12). Note the clearance shift is SMALLER in magnitude than the non-Chinese Asian one, so the two Asian clearance levels must not be collapsed.",
      source_name = "RCE = 14"
    ),
    RACE_JAPANESE = list(
      description = "Japanese-heritage race indicator (Li 2019 control stream RCE = 13, abbreviated JPN)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = non-Japanese",
      notes = "45 subjects (Li 2019 Table 3). Japanese subjects take a share of the pooled Asian VOLUME shift (-0.27, theta12) but carry NO clearance shift -- the control stream's CLRCE block assigns a factor only to RCE = 3 and RCE = 14, leaving RCE = 13 at the non-Asian reference. Do not borrow the Chinese clearance coefficient for Japanese subjects; the resulting higher exposure is what Li 2019 Results describes as 'Japanese patients had higher ceftazidime and avibactam exposure than the white/other reference population'.",
      source_name = "RCE = 13"
    ),
    MECH_VENT = list(
      description = "Presence of a ventilator in the hospital room on the day of PK sampling (Li 2019 control stream POP5 = 1, abbreviated NPv)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = no ventilator present on the PK sampling day",
      notes = "Li 2019 defines NPv as 'the presence of a ventilator in the hospital room on the day of PK sampling, which includes patients with VAP or HAP who were ventilated on the day of sampling' (Methods, 'Selection of covariates'). Time-fixed at the PK sampling day rather than re-evaluated across the ICU stay, and NOT the same contrast as VABP vs HABP -- a HABP patient ventilated on the sampling day has MECH_VENT = 1 and DIS_VABP = 0. Raises Vc by 29.7% (Table 1 theta14), the second of the two ceftazidime fixed effects with RSE above 27% (45.4%).",
      source_name = "POP5 = 1"
    ),
    STUDY_CAZAVI_PHASE2 = list(
      description = "Phase 2 stratum indicator of the Li 2019 pooled ceftazidime-avibactam analysis (residual-error switch)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 with STUDY_CAZAVI_PHASE3 also 0 selects the Phase 1 stratum",
      notes = "Two Phase 2 studies (cIAI and cUTI) contributed to the pooled data set (Li 2019 Methods, 'Analysis data and model construction'). The ceftazidime residual-error model pools Phase 2 with Phase 3 into ONE stratum (control stream $ERROR: 'IF (PHASE.EQ.2.OR.PHASE.EQ.3)'), so this model uses only the SUM of the two indicators; they are kept as separate columns so the same simulated covariate set drives the avibactam companion, which does distinguish them.",
      source_name = "PHASE = 2"
    ),
    STUDY_CAZAVI_PHASE3 = list(
      description = "Phase 3 stratum indicator of the Li 2019 pooled ceftazidime-avibactam analysis (residual-error switch)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 with STUDY_CAZAVI_PHASE2 also 0 selects the Phase 1 stratum",
      notes = "Five Phase 3 trials: RECLAIM 1/2, RECLAIM 3, RECAPTURE 1/2, REPRISE and REPROVE (Li 2019 Methods). Pooled with Phase 2 for the ceftazidime residual error; see STUDY_CAZAVI_PHASE2.",
      source_name = "PHASE = 3"
    )
  )

  # Screened during stepwise covariate selection but NOT retained in the
  # final ceftazidime model (Li 2019 Methods, 'Selection of covariates';
  # Results, 'Ceftazidime'). Documented here so the paper's covariate screen
  # is preserved without declaring covariates model() never references.
  covariatesDataExcluded <- list(
    APACHE_II_SEV = list(
      description = "Elevated-APACHE-II (score > 10) severity stratum indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = APACHE II <= 10",
      notes = "Screened on ceftazidime CL and Vc; not retained. The companion avibactam model DOES retain it (-0.197 on CL), so the asymmetry is a real feature of the pair rather than an incomplete table.",
      source_name = "APACHE"
    ),
    DIS_BACTEREMIA = list(
      description = "Baseline bacteremia indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = no bacteremia at baseline",
      notes = "88 of 1,553 phase III subjects with exposure data (Li 2019 Table 3). Screened; not retained. Exposures differed by <= 25% between bacteremic and non-bacteremic patients.",
      source_name = "BBACTERM"
    ),
    RENALIMP_ESRD = list(
      description = "End-stage renal disease indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = renal function above ESRD",
      notes = "The ceftazidime data set 'lacked data for subjects with severe renal impairment', so individual CL estimates from patients with renal insufficiency reported in the literature were incorporated into the base model (Li 2019 Methods, Data S2) rather than an ESRD indicator being fitted. The continuous CRCL spline alone spans the anuric range. The companion avibactam model, whose data set did include ESRD subjects, carries an explicit ESRD term.",
      source_name = "ESRD"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = male",
      notes = "Screened; not retained (Li 2019 Methods, 'Selection of covariates').",
      source_name = "SEX"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained. Li 2019 Results: 'Age-related or obesity-related changes in exposure seemed to be adequately captured by changes in CrCL.'",
      source_name = "AGE"
    ),
    BMI = list(
      description = "Body mass index (obesity status)",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as obesity status; not retained. See AGE note.",
      source_name = "BMI"
    ),
    WBC = list(
      description = "White blood cell count",
      units = "cells/uL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened at a 12,000/uL cutoff as a marker of systemic disturbance; not retained.",
      source_name = "WBC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1975,
    n_observations = 9155,
    n_studies = 18,
    disease_state = "Adults with complicated intra-abdominal infection (781, 39.5%), complicated urinary tract infection including acute pyelonephritis (696, 35.2%) or nosocomial pneumonia including ventilator-associated pneumonia (412, 20.9%), plus 86 healthy volunteers (4.4%)",
    renal_function = "Estimated Cockcroft-Gault creatinine clearance 8-488 mL/min. The data set lacked subjects with severe renal impairment, so individual clearance estimates from patients with renal insufficiency reported in the literature were incorporated into the base model (Li 2019 Methods, Data S2).",
    dose_range = "Ceftazidime 2000 mg every 8 hours as a 2-hour intravenous infusion in subjects with CrCL > 50 mL/min, with label-recommended reductions to 1000 mg q8h (CrCL 31-50), 750 mg q12h (CrCL 16-30), 750 mg q24h (CrCL 6-15) and 750 mg q48h (ESRD). Given in a fixed 4:1 ratio with avibactam.",
    regions = "Global, including dedicated Chinese, Japanese, Korean, Taiwanese and Vietnamese subgroups",
    studies = "11 Phase 1 studies, 2 Phase 2 studies (cIAI and cUTI) and 5 Phase 3 trials: RECLAIM 1/2, RECLAIM 3, RECAPTURE 1/2, REPRISE and REPROVE",
    unbound_fraction = "0.85 for ceftazidime (Li 2019 Methods, 'Exposure-response analysis': free plasma concentrations were taken to be 85% of total). The model predicts TOTAL plasma concentration.",
    notes = "Estimation used FOCE-INTER in NONMEM 7.2 for model building; the final models were re-estimated with SAEM followed by importance sampling to improve prediction at the 10th percentile. Outliers with conditional weighted residual > 4, and ceftazidime concentrations above 750 mg/L, were excluded from the final model."
  )

  ini({
    # =====================================================================
    # Structural parameters (Li 2019 Table 1). Values refer to the
    # reference subject: CRCL = 97.1 mL/min (where the clearance spline
    # equals 1), WT = 70 kg, healthy volunteer (all DIS_* indicators 0),
    # non-Asian, not mechanically ventilated.
    # =====================================================================
    lcl <- log(6.95); label("Ceftazidime clearance at CRCL = 97.1 mL/min (L/h)")     # Table 1: theta1 CL = 6.95 L/h (RSE 1.7%)
    lvc <- log(10.5); label("Ceftazidime central volume of distribution at WT = 70 kg (L)") # Table 1: theta2 Vc = 10.5 L (RSE 13.1%)
    lq <- log(31.5); label("Ceftazidime inter-compartmental clearance (L/h)")        # Table 1: theta3 Q = 31.5 L/h (RSE 18.8%)
    lvp <- log(7.57); label("Ceftazidime peripheral volume of distribution (L)")     # Table 1: theta4 Vp = 7.57 L (RSE 9%)

    # =====================================================================
    # Creatinine-clearance effect on clearance: a two-segment LINEAR spline
    # hinged at 100 mL/min.
    #
    #   CRCL <  100:  slope1 * CRCL
    #   CRCL >= 100:  slope1 * 100 + slope2 * (CRCL - 100)
    #
    # The lower segment passes through the ORIGIN, not through 1 at the
    # hinge, so the factor is 1 at CRCL = 1 / slope1 = 97.1 mL/min and
    # 1.030 at the hinge itself. Both slopes are hard-coded constants in
    # the final control stream (supplementary Data S1, '$PK' block:
    # 'SLOPE1=0.01030360' / 'SLOPE2=0.00125182' immediately under the
    # comment ';SCALED SO THAT TVCL = 7.07 (HS)'), so they are encoded
    # fixed() here even though Li 2019 Table 1 quotes an RSE for each --
    # the RSEs belong to the earlier estimation step that produced them.
    # See the validation vignette's Assumptions and deviations section.
    # =====================================================================
    e_crcl_cl_lt100 <- fixed(0.0103036); label("Ceftazidime CRCL linear slope on CL below 100 mL/min (per mL/min)")
    # Table 1 row 'Slope 1: CrCL < 100 mL/min, slope1*CrCL' = 0.0103036 (RSE 0.409%);
    # identical to the control-stream constant SLOPE1.
    e_crcl_cl_ge100 <- fixed(0.00125182); label("Ceftazidime CRCL linear slope on CL at or above 100 mL/min (per mL/min)")
    # Control-stream constant SLOPE2 = 0.00125182; Li 2019 Table 1 prints the
    # same number rounded to 0.001252 (RSE 8.84%). Results: '12.5% increase in
    # CL per 100 mL/min increase in CrCL above 100 mL/min'. That figure is
    # 100 * 0.00125182 = 0.125 read as a percentage, i.e. the ABSOLUTE
    # increment in the renal factor. Because this model's factor is 1.030 at
    # the hinge rather than 1 (the lower segment passes through the origin),
    # the increment expressed as a RELATIVE rise in clearance is
    # 0.125182 / 1.03036 = 12.1%, not 12.5%. Both readings describe the same
    # slope; the vignette gates on the slope itself. Contrast the avibactam
    # companion, whose factor IS exactly 1 at its hinge, so its analogous
    # '27.9% per 100 mL/min' is exactly 100 * 0.00279 under either reading.

    # =====================================================================
    # Infection-type and population effects.
    #
    # NOTE ON FORM: unlike the avibactam companion (and unlike the
    # successor Das 2024 / Xie 2025 models, which write every population
    # effect as X*(1 + theta)), the ceftazidime control stream applies the
    # cIAI / NP / cUTI effects as BARE MULTIPLIERS -- 'IF(POP.EQ.3)
    # CLPOP=THETA(5)', not '(1 + THETA(5))'. The values below are
    # therefore multiplicative factors, where 1 means no effect. The
    # pyelonephritis, race and ventilation effects in the SAME control
    # stream do use the (1 + theta) form; each is commented accordingly.
    # =====================================================================
    e_ciai_cl <- 1.16; label("Multiplicative factor on ceftazidime CL for cIAI (unitless factor, not a 1+theta shift)")
    # Table 1: theta5 'Population effect on CL (cIAI)' = 1.16 (RSE 2.2%).
    # Results: '16% higher CL for patients with cIAI vs. healthy subjects and
    # patients with cUTI'.
    e_habp_vabp_cl <- 0.999; label("Multiplicative factor on ceftazidime CL for nosocomial pneumonia, HAP or VAP (unitless factor)")
    # Table 1: theta6 'Population effect on CL (NP)' = 0.999 (RSE 2.4%), i.e.
    # NP clearance is indistinguishable from the healthy-volunteer reference.
    # One coefficient covers both pneumonia arms; applied to DIS_HABP + DIS_VABP.

    e_cuti_vc <- 1.03; label("Multiplicative factor on ceftazidime Vc for cUTI (unitless factor)")   # Table 1: theta9 'Population effect on Vc (cUTI)' = 1.03 (RSE 11.1%)
    e_ciai_habp_vabp_vc <- 1.14; label("Multiplicative factor on ceftazidime Vc for cIAI or nosocomial pneumonia (unitless factor)")
    # Table 1: theta10 'Population effect on Vc (cIAI or NP)' = 1.14 (RSE 9.9%).
    # A single coefficient shared by cIAI, HAP and VAP; applied to
    # DIS_CIAI + DIS_HABP + DIS_VABP.

    e_dis_ap_vc <- -0.185; label("Proportional shift in ceftazidime Vc for acute pyelonephritis within the cUTI cohort (fraction)")
    # Table 1: theta11 'Population effect on Vc (cUTI/acute pyelonephritis)' =
    # -0.185 (RSE 41.2%). Control stream applies it as (1 + theta11), stacking
    # multiplicatively on the cUTI factor.

    # =====================================================================
    # Race effects. Li 2019 splits Asian heritage into three mutually
    # exclusive levels. Clearance distinguishes non-Chinese Asian from
    # Chinese and leaves Japanese at the non-Asian reference; central
    # volume applies ONE pooled shift to all three levels.
    # =====================================================================
    e_race_asian_oth_cl <- -0.161; label("Proportional shift in ceftazidime CL for non-Chinese, non-Japanese Asian race (fraction)")
    # Table 1: theta7 'Race effect on CL (ASN)' = -0.161 (RSE 11.8%).
    # Results: 'non-Chinese, non-Japanese Asians had 16% lower CL'.
    e_race_chinese_cl <- -0.0855; label("Proportional shift in ceftazidime CL for Chinese race (fraction)")
    # Table 1: theta8 'Race effect on CL' = -0.0855 (RSE 27%).
    # Results: 'Chinese patients had 9% lower CL'.
    e_race_asian_vc <- -0.27; label("Proportional shift in ceftazidime Vc for any Asian race, pooled across the three levels (fraction)")
    # Table 1: theta12 'Race effect on Vc (ASN, CHN, JPN)' = -0.27 (RSE 18.6%).
    # Results: '27% lower Vc for Asian compared with non-Asian patients'.
    # Control stream applies it when RCE is 3, 13 or 14, i.e. to the sum
    # RACE_ASIAN_OTH + RACE_CHINESE + RACE_JAPANESE.

    e_wt_vc <- 1.01; label("Body-weight power exponent on ceftazidime Vc, WT/70 (unitless)")
    # Table 1: theta13 'WT effect on Vc' = 1.01 (RSE 12.6%). Results: weights at
    # the 10th (50 kg) and 90th (94 kg) percentiles give 24% lower and 26%
    # higher Vc than the 70 kg median, which reproduce from (WT/70)^1.01.

    e_mech_vent_vc <- 0.297; label("Proportional shift in ceftazidime Vc when a ventilator is present on the PK sampling day (fraction)")
    # Table 1: theta14 'Population effect on Vc (NPv)' = 0.297 (RSE 45.4%).
    # Results: '29.7% higher Vc for patients with NPv than for non-NPv patients'.

    # =====================================================================
    # Inter-individual variability: a full 4x4 OMEGA block over CL, Vc, Vp
    # and Q, in the control stream's ETA order (ETA1 CL, ETA2 V1, ETA3 V2,
    # ETA4 Q).
    #
    # SCALE. Li 2019 Table 1 footnote b marks every one of these rows
    # 'Reported as variance', so the diagonal entries are log-scale
    # variances and the off-diagonals are covariances -- they are used
    # here verbatim, with no CV%-to-variance conversion. The table's
    # separate 'BSV (CV%)' column carries etashrinkage for these rows, not
    # a second copy of the variance.
    #
    # VERIFICATION. All six correlations printed in Table 1 reproduce from
    # these ten numbers: r(Vc,CL) = -0.189/sqrt(0.179*1.10) = -0.43
    # (table -0.42); r(Vp,CL) = 0.383/sqrt(0.179*1.21) = 0.82 (0.82);
    # r(Vp,Vc) = -0.972/sqrt(1.10*1.21) = -0.84 (-0.84);
    # r(Q,CL) = 0.883/sqrt(0.179*6.70) = 0.81 (0.81);
    # r(Q,Vc) = -0.643/sqrt(1.10*6.70) = -0.24 (-0.24);
    # r(Q,Vp) = 1.73/sqrt(1.21*6.70) = 0.61 (0.61).
    # =====================================================================
    etalcl + etalvc + etalvp + etalq ~ c(
      0.179,
      -0.189, 1.10,
      0.383, -0.972, 1.21,
      0.883, -0.643, 1.73, 6.70
    )

    # =====================================================================
    # Residual unexplained variability, switched by study phase.
    #
    # SCALE. Li 2019 Table 1 footnote b marks all four residual rows
    # 'Reported as variance', and the control stream's $SIGMA block holds
    # EPS variances with Y = F + F*EPS(1) + EPS(2). The nlmixr2 propSd /
    # addSd parameters are STANDARD DEVIATIONS, so each table value is
    # square-rooted below.
    #
    # UNITS. The control stream sets S1 = V1/1000 with doses in mg, so its
    # predictions -- and therefore its ADDITIVE error terms -- are in
    # ng/mL. This model works in mg/L (see units, above), so the additive
    # standard deviations are divided by 1000. The proportional terms are
    # dimensionless and carry over unchanged.
    # =====================================================================
    propSdPhase1 <- 0.2; label("Proportional residual standard deviation, Phase 1 studies (fraction)")
    # Table 1: 'Proportional error, phase I' = 0.04 as a variance; sqrt(0.04) = 0.2.
    addSdPhase1 <- 0.16276; label("Additive residual standard deviation, Phase 1 studies (mg/L)")
    # Table 1: 'Additive error, phase I' = 26489 as a variance in (ng/mL)^2;
    # sqrt(26489) = 162.76 ng/mL = 0.16276 mg/L.
    propSdPhase23 <- 0.337639; label("Proportional residual standard deviation, Phase 2 and Phase 3 studies (fraction)")
    # Table 1: 'Proportional error, phase II and phase III' = 0.114 as a
    # variance; sqrt(0.114) = 0.337639.
    addSdPhase23 <- 0.0042895; label("Additive residual standard deviation, Phase 2 and Phase 3 studies (mg/L)")
    # Table 1: 'Additive error, phase II and phase III' = 18.4 as a variance in
    # (ng/mL)^2; sqrt(18.4) = 4.2895 ng/mL = 0.0042895 mg/L. Li 2019 reports
    # RSE 447% for this term, so it is effectively unidentified and the Phase
    # 2/3 residual is proportional in all but name.
  })

  model({
    # ------------------------------------------------------------------
    # 1. Renal-function factor on clearance. Two linear segments hinged at
    #    100 mL/min, written with 0/1 indicator arithmetic rather than a
    #    branch so both arms stay finite for every positive CRCL.
    #
    #    Unlike the avibactam companion's power-then-linear hinge, the
    #    lower arm here passes through the origin, so this factor is NOT 1
    #    at the hinge -- it is 1.0030 there and 1 at CRCL = 97.1.
    # ------------------------------------------------------------------
    renal_cl <- e_crcl_cl_lt100 * CRCL * (CRCL < 100) +
      (e_crcl_cl_lt100 * 100 + e_crcl_cl_ge100 * (CRCL - 100)) * (CRCL >= 100)

    # ------------------------------------------------------------------
    # 2. Infection-type factors. The four indicators are mutually
    #    exclusive, so a subject picks up at most one term from each
    #    bracket and an all-zero subject keeps the healthy-volunteer
    #    reference of 1. Written in the paper's bare-multiplier form:
    #    (1 - indicator) + factor * indicator.
    # ------------------------------------------------------------------
    infect_cl <- (1 - DIS_CIAI - DIS_HABP - DIS_VABP) +
      e_ciai_cl * DIS_CIAI +
      e_habp_vabp_cl * (DIS_HABP + DIS_VABP)

    infect_vc <- (1 - DIS_CUTI - DIS_CIAI - DIS_HABP - DIS_VABP) +
      e_cuti_vc * DIS_CUTI +
      e_ciai_habp_vabp_vc * (DIS_CIAI + DIS_HABP + DIS_VABP)

    # ------------------------------------------------------------------
    # 3. Individual PK parameters. Control stream:
    #      CL = CLCRCL * exp(log(THETA(1) * CLPOP * CLRCE) + ETA(1))
    #      V1 = V1WT   * exp(log(THETA(2) * V1POP * V1POP1 * V1RCE * V1POP5) + ETA(2))
    #      V2 = THETA(4) * exp(ETA(3))
    #      Q  = THETA(3) * exp(ETA(4))
    #    Neither Q nor Vp carries any covariate.
    # ------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * renal_cl * infect_cl *
      (1 + e_race_asian_oth_cl * RACE_ASIAN_OTH + e_race_chinese_cl * RACE_CHINESE)

    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc * infect_vc *
      (1 + e_dis_ap_vc * DIS_AP) *
      (1 + e_race_asian_vc * (RACE_ASIAN_OTH + RACE_CHINESE + RACE_JAPANESE)) *
      (1 + e_mech_vent_vc * MECH_VENT)

    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)

    # ------------------------------------------------------------------
    # 4. Micro-constants and the two-compartment IV disposition. Dosing is
    #    a zero-order intravenous infusion into the central compartment.
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ------------------------------------------------------------------
    # 5. Observation. Dose in mg, volumes in L -> concentration in mg/L.
    #    This is the TOTAL plasma concentration; multiply by 0.85 to
    #    obtain the free concentration the paper's 50% fT > 8 mg/L target
    #    is defined on (see population$unbound_fraction).
    #
    #    Residual variability is switched by study phase; both phase
    #    indicators 0 selects the Phase 1 stratum.
    # ------------------------------------------------------------------
    Cc <- central / vc

    phase23 <- STUDY_CAZAVI_PHASE2 + STUDY_CAZAVI_PHASE3
    propSd <- propSdPhase1 * (1 - phase23) + propSdPhase23 * phase23
    addSd <- addSdPhase1 * (1 - phase23) + addSdPhase23 * phase23

    Cc ~ add(addSd) + prop(propSd)
  })
}
