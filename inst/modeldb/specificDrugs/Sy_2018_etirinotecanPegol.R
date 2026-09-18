Sy_2018_etirinotecanPegol <- function() {
  description <- paste(
    "Integrated five-analyte population PK model for etirinotecan pegol (EP,",
    "NKTR-102) and its metabolic cascade in adults with advanced solid tumors.",
    "EP is a four-arm polyethylene-glycol conjugate of irinotecan; the model",
    "follows EP -> irinotecan -> (SN-38 | APC) and SN-38 -> SN-38 glucuronide,",
    "with a two-compartment disposition model for each of the five analytes.",
    "Metabolite conversion fractions and volumes are not separately",
    "identifiable, so every metabolite volume is an aggregate ratio and only",
    "EP's central volume is a true volume."
  )
  reference <- paste(
    "Sy SKB, Chia YL, Gordi T, Hoch U, Eldon MA (2018).",
    "Integrated population pharmacokinetics of etirinotecan pegol and its four",
    "metabolites in cancer patients with solid tumors.",
    "Cancer Chemother Pharmacol 81(5):897-909.",
    "doi:10.1007/s00280-018-3562-3."
  )
  vignette <- "Sy_2018_etirinotecanPegol"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list(
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on EP clearance, normalised to 60 years (Sy 2018 Eq. 5).",
        "Median 60 years, range 25-81 (Table 2). Retained in the final model",
        "but judged to have no clinical impact: Figure 4 shows cumulative EP",
        "AUC over six cycles of 6.41 mg.h/mL at 45 years and 7.22 mg.h/mL at",
        "75 years, against 6.87 mg.h/mL for the 60-year reference."
      ),
      source_name = "Age"
    ),
    BSA = list(
      description = "Baseline body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on both EP clearance and EP central volume, normalised",
        "to 1.86 m^2 (Sy 2018 Eqs. 5 and 6). Median 1.86 m^2, range 1.36-2.74",
        "(Table 2). EP is dosed per m^2, so this covariate is already absorbed",
        "by the BSA-based dosing scheme (Sy 2018 Discussion, Conclusions)."
      ),
      source_name = "BSA"
    ),
    CRCL = list(
      description = paste(
        "Baseline estimated glomerular filtration rate, BSA-normalised.",
        "Sy 2018 computed it from serum creatinine with the IDMS-traceable",
        "4-variable MDRD equation printed as their Eq. 1:",
        "GFR = 175 * SCr^-1.154 * Age^-0.203 * 0.742^I(female) *",
        "1.212^I(Black), SCr in mg/dL."
      ),
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on EP clearance, normalised to 84.1 mL/min/1.73 m^2",
        "(Sy 2018 Eq. 5, corroborated by the simulation reference population",
        "in Results 'Clinical impact of significant covariates', which states",
        "'renal function (84.1 mL/min)'). NOTE the Table 4 caption instead",
        "says the typical patient has eGFR = 79.1 mL/min, which is the Table 2",
        "median; the printed equation wins per the standing text-vs-equation",
        "rule and the difference is 1.2% on CL. NOTE also that Table 1 labels",
        "the covariate 'eGFR CG' and its Comments column says 'Calculated from",
        "serum creatinine using Corrected CKD-EPI', but the equation actually",
        "printed as Eq. 1 is MDRD, not Cockcroft-Gault and not CKD-EPI.",
        "Cohort median 79.1, range 34.9-216.6 (Table 2); 31% normal, 48% mild",
        "and 19% moderate renal impairment."
      ),
      source_name = "eGFR"
    ),
    SNP_UGT1A1_RS8175347_HOM = list(
      description = "UGT1A1*28 (rs8175347) homozygous-variant TA(7)/TA(7) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = paste(
        "0 = wild-type TA(6)/TA(6) or heterozygous TA(6)/TA(7). Sy 2018 fit a",
        "single indicator for the homozygous genotype only, so wild-type and",
        "heterozygous subjects share the reference k3e (Eq. 7)."
      ),
      notes = paste(
        "Exponential effect on the SN-38 elimination rate constant k3e",
        "(Sy 2018 Eq. 7). Genotype distribution (Table 2): 38.6% none, 48.2%",
        "one copy, 10.8% two copies, 2.4% indeterminate; the two indeterminate",
        "subjects were pooled into the wild-type category (Methods, 'Covariate",
        "screening'). Only nine homozygous subjects were available, which the",
        "authors flag as limiting generalisation (Results and Discussion)."
      ),
      source_name = "UGT1A1*28 TA(7)/TA(7)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened as a body-size descriptor; BSA was retained instead (Table 1). Median 72.3 kg, range 43.9-153.6 (Table 2)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened on k3e through the WAM algorithm (Table 3 ranks 3, 6, 8, 9) but not retained in the final model. 45.8% female (Table 2)."
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Baseline value collected (Methods, 'Covariate screening') but not carried into the formal WAM screen. Median 36 g/L, range 14-45 (Table 2)."
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Table 1 disposition: 'Removed; bilirubin used as indicator of hepatic function' because ALT correlates with bilirubin. Median 21 U/L, range 7-153 (Table 2)."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Table 1 disposition: 'Removed; bilirubin used as indicator of hepatic function' because AST correlates with ALT. Median 24 U/L, range 11-130 (Table 2)."
    ),
    DBIL = list(
      description = "Baseline total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Table 1 disposition: 'Retained as indicator of hepatic function' for screening, but it does not appear in any of the ten WAM candidate models of Table 3 and is absent from the final model. Median 10.3 umol/L, range 3.4-27.4 (Table 2). Table 1 quotes it in mg/dL; Table 2 reports SI umol/L."
    ),
    CREAT = list(
      description = "Baseline serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Enters the model only through the MDRD equation that produces CRCL (Eq. 1); it is not itself a model covariate. Median 79.6 umol/L, range 35.4-132.6 (Table 2)."
    ),
    SMOKE = list(
      description = "Current smoking indicator",
      units = "(binary)",
      type = "binary",
      notes = "Table 1 disposition: 'Retained as potential clinical covariate' for screening; it does not appear in any Table 3 candidate model and is absent from the final model."
    ),
    RACE_BLACK = list(
      description = "Black race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Table 1 disposition: 'Enrollment of very few non-Caucasian patients precluded meaningful analysis' (94% White, Table 2). Race enters only through the MDRD equation of Eq. 1, not as a model covariate."
    )
  )

  compartmentData <- list(
    central = list(analyte = "etirinotecan pegol", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "etirinotecan pegol", units = "nmol", specimen = "plasma", verified = TRUE),
    central_irinotecan = list(analyte = "irinotecan", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1_irinotecan = list(analyte = "irinotecan", units = "nmol", specimen = "plasma", verified = TRUE),
    central_sn38 = list(analyte = "SN-38", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1_sn38 = list(analyte = "SN-38", units = "nmol", specimen = "plasma", verified = TRUE),
    central_sn38g = list(analyte = "SN-38 glucuronide", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1_sn38g = list(analyte = "SN-38 glucuronide", units = "nmol", specimen = "plasma", verified = TRUE),
    central_apc = list(analyte = "APC", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1_apc = list(analyte = "APC", units = "nmol", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 83L,
    n_studies = 2L,
    age_range = "25-81 years",
    age_median = "60 years",
    weight_range = "43.9-153.6 kg",
    weight_median = "72.3 kg",
    bsa_range = "1.36-2.74 m^2",
    bsa_median = "1.86 m^2",
    sex_female_pct = 45.8,
    race_ethnicity = c(White = 94, Other = 6),
    disease_state = "advanced solid tumors (colorectal 22.9%, lung 18.1%, pancreas 10.8%, ovarian 6%, breast 3.6%, other 38.6%)",
    dose_range = "etirinotecan pegol by 90-min intravenous infusion; three weekly doses every 4 weeks, once every 2 weeks, or once every 3 weeks (06-IN-IR001 dose escalation); 100 or 125 mg/m^2 once every 3 weeks with cetuximab (07-PIR-02)",
    renal_function = "eGFR median 79.1 mL/min/1.73 m^2 (range 34.9-216.6); 31% normal, 48% mild impairment, 19% moderate impairment",
    genotype = "UGT1A1*28: 38.6% none, 48.2% one copy, 10.8% two copies, 2.4% indeterminate",
    co_medication = "weak/moderate CYP3A4 inhibitors 14.5%, weak/moderate CYP3A4 inducers 44.6%; no strong inhibitors or inducers",
    notes = paste(
      "Pooled phase 1 data from 06-IN-IR001 (67 of 76 enrolled patients) and",
      "07-PIR-02 (16 of 18 patients, etirinotecan pegol plus cetuximab);",
      "baseline characteristics from Sy 2018 Table 2. The analysis data set",
      "contained 1414, 1777, 1769, 1731 and 1167 quantifiable concentrations",
      "for EP, irinotecan, SN-38, SN-38G and APC respectively; samples below",
      "the limit of quantification (22.6, 3.8, 4.1, 6.1 and 36.4% of records)",
      "were excluded as missing. Estimation used SAEM followed by importance",
      "sampling in Monolix 2016."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Stage 1 of a two-stage sequential fit (Sy 2018 '3-analyte
    # pharmacokinetic model'). EP, irinotecan and SN-38 parameters were
    # estimated first; their individual estimates were then HELD FIXED while
    # the SN-38G and APC data were added to build the 5-analyte model, with
    # the exception of Corr(CL, V1) and Corr(k2e, V2*), which the software
    # could not fix and which were re-estimated. Every value below carries a
    # standard error in Table 4, so none is encoded with fixed().
    # ---------------------------------------------------------------------

    # ---- Etirinotecan pegol (parent) ----
    lcl <- log(0.237); label("EP clearance (L/h)") # Sy 2018 Table 4, 3-analyte model, EP row 'CL (L/h)': 0.237 +/- 0.008
    lvc <- log(5.05); label("EP central volume, a true volume (L)") # Sy 2018 Table 4, EP row 'V1 (L)': 5.05 +/- 0.14. The Abstract quotes 5.5 L; Table 4 is the parameter table and wins.
    lk12 <- log(6.78e-3); label("EP central-to-peripheral rate constant k1p (1/h)") # Sy 2018 Table 4, EP row 'k1p (h-1)': 6.78e-3 +/- 5.5e-4
    lk21 <- log(5.8e-4); label("EP peripheral-to-central rate constant kp1 (1/h)") # Sy 2018 Table 4, EP row 'kp1 (h-1)': 5.8e-4 +/- 4.1e-5

    e_age_cl <- -0.271; label("Power exponent of AGE/60 on EP clearance (unitless)") # Sy 2018 Table 4, EP row 'theta CL,AGE': -0.271 +/- 0.084; applied per Eq. 5
    e_bsa_cl <- 1.32; label("Power exponent of BSA/1.86 on EP clearance (unitless)") # Sy 2018 Table 4, EP row 'theta CL,BSA': 1.32 +/- 0.23; applied per Eq. 5
    e_crcl_cl <- 0.2; label("Power exponent of CRCL/84.1 on EP clearance (unitless)") # Sy 2018 Table 4, EP row 'theta CL,eGFR': 0.2 +/- 0.068; applied per Eq. 5
    e_bsa_vc <- 1.1; label("Power exponent of BSA/1.86 on EP central volume (unitless)") # Sy 2018 Table 4, EP row 'theta V1,BSA': 1.1 +/- 0.19; applied per Eq. 6

    # ---- Irinotecan ----
    lvc_irinotecan <- log(1.8); label("Irinotecan aggregate central volume V2* = V2/F12 (L)") # Sy 2018 Table 4, Irinotecan row 'V2* = V2/F12 (L)': 1.8 +/- 0.13
    lkel_irinotecan <- log(27.6); label("Irinotecan total elimination rate constant k2e (1/h)") # Sy 2018 Table 4, Irinotecan row 'k2e (h-1)': 27.6 +/- 1.7
    lk12_irinotecan <- log(18.8); label("Irinotecan central-to-peripheral rate constant k2p (1/h)") # Sy 2018 Table 4, Irinotecan row 'k2p (h-1)': 18.8 +/- 1.2
    lk21_irinotecan <- log(3.2e-3); label("Irinotecan peripheral-to-central rate constant kp2 (1/h)") # Sy 2018 Table 4, Irinotecan row 'kp2 (h-1)': 3.2e-3 +/- 1.5e-4

    # ---- SN-38 ----
    lvc_sn38 <- log(80); label("SN-38 aggregate central volume V3** = V3/(F12*F23) (L)") # Sy 2018 Table 4, SN38 row 'V3** = V3/(F12 F23) (L)': 80 +/- 190. See the vignette 'Assumptions and deviations': this printed value is the one parameter in Table 4 whose standard error exceeds its estimate, and it is inconsistent by a factor of 4.6 with the absolute SN-38 exposures the same paper prints in Figure 5 and with the SN-38 concentration range of the Figure 3 goodness-of-fit panels.
    lkel_sn38 <- log(0.0602); label("SN-38 total elimination rate constant k3e (1/h)") # Sy 2018 Table 4, SN38 row 'k3e (h-1)': 0.0602 +/- 0.0042
    lk12_sn38 <- log(0.23); label("SN-38 central-to-peripheral rate constant k3p (1/h)") # Sy 2018 Table 4, SN38 row 'k3p (h-1)': 0.23 +/- 0.023
    lk21_sn38 <- log(8.75e-3); label("SN-38 peripheral-to-central rate constant kp3 (1/h)") # Sy 2018 Table 4, SN38 row 'kp3 (h-1)': 8.75e-3 +/- 5.4e-4
    e_snp_ugt1a1_rs8175347_hom_kel_sn38 <- -0.67; label("Exponential effect of UGT1A1*28 TA(7)/TA(7) on SN-38 k3e (unitless)") # Sy 2018 Table 4, SN38 row 'theta k3e,UGT1A1': -0.67 +/- 0.21; applied per Eq. 7

    # ---- SN-38 glucuronide ----
    lvc_sn38g <- log(11.6); label("SN-38G aggregate central volume V4* = V4/(F12*(F23+F25)*F34) (L)") # Sy 2018 Table 4, 5-analyte model, SN38G row 'V4* (L)': 11.6 +/- 0.95
    lkel_sn38g <- log(1.41); label("SN-38G total elimination rate constant k4e (1/h)") # Sy 2018 Table 4, SN38G row 'k4e (h-1)': 1.41 +/- 0.093
    lk12_sn38g <- log(0.548); label("SN-38G central-to-peripheral rate constant k4p (1/h)") # Sy 2018 Table 4, SN38G row 'k4p (h-1)': 0.548 +/- 0.056
    lk21_sn38g <- log(0.104); label("SN-38G peripheral-to-central rate constant kp4 (1/h)") # Sy 2018 Table 4, SN38G row 'kp4 (h-1)': 0.104 +/- 0.014

    # ---- APC ----
    lfm_sn38 <- log(0.631); label("Fraction of converted irinotecan that becomes SN-38 rather than APC (unitless)") # Sy 2018 Table 4, 5-analyte model, APC row '(F irinotecan->SN38)': 0.631 +/- 0.017
    lkel_apc <- log(0.0235); label("APC total elimination rate constant k5e (1/h)") # Sy 2018 Table 4, APC row 'k5e (h-1)': 0.0235 +/- 0.0027
    lk12_apc <- log(0.0236); label("APC central-to-peripheral rate constant k5p (1/h)") # Sy 2018 Table 4, APC row 'k5p (h-1)': 0.0236 +/- 0.003
    lk21_apc <- log(1.39e-3); label("APC peripheral-to-central rate constant kp5 (1/h)") # Sy 2018 Table 4, APC row 'kp5 (h-1)': 1.39e-3 +/- 1.8e-4

    # ---- Interindividual variability ----
    # Log-normal IIV throughout (Methods, 'Pharmacokinetic model structure').
    # Table 4's right-hand column reports the VARIANCE with its SE, and the
    # parenthesised CV% is sqrt(variance) for every row, so the tabulated
    # numbers are used directly as omega variances on the log scale.
    # Off-diagonals use the 5-analyte 'Correlations' block at the foot of
    # Table 4, which supersedes the 3-analyte values of 0.759 and -0.751.
    # cov = corr * sqrt(var_a * var_b): 0.713 * sqrt(0.07654 * 0.051) and
    # -0.755 * sqrt(0.204 * 0.354).
    etalcl + etalvc ~ c(0.07654, 0.04454, 0.051) # Sy 2018 Table 4 EP rows 'CL' 0.07654 +/- 0.014 (27%) and 'V1' 0.051 +/- 0.01 (22%); Corr(CL, V1) 0.713 +/- 0.056 from the 5-analyte 'Correlations' block
    etalk12 ~ 0.295 # Sy 2018 Table 4, EP row 'k1p': 0.295 +/- 0.07 (54%)
    etalkel_irinotecan + etalvc_irinotecan ~ c(0.204, -0.20289, 0.354) # Sy 2018 Table 4 Irinotecan rows 'k2e' 0.204 +/- 0.045 (45%) and 'V2*' 0.354 +/- 0.062 (59%); Corr(k2e, V2) -0.755 +/- 0.046 from the 5-analyte 'Correlations' block
    etalk12_irinotecan ~ 0.214 # Sy 2018 Table 4, Irinotecan row 'k2p': 0.214 +/- 0.04 (46%)
    etalvc_sn38 ~ 0.14 # Sy 2018 Table 4, SN38 row 'V3**': 0.14 +/- 0.028 (37%)
    etalkel_sn38 ~ 0.224 # Sy 2018 Table 4, SN38 row 'k3e': 0.224 +/- 0.048 (47%)
    etalk12_sn38 ~ 0.561 # Sy 2018 Table 4, SN38 row 'k3p': 0.561 +/- 0.12 (75%)
    etalvc_sn38g ~ 0.415 # Sy 2018 Table 4, SN38G row 'V4*': 0.415 +/- 0.081 (64%)
    etalkel_sn38g ~ 0.232 # Sy 2018 Table 4, SN38G row 'k4e': 0.232 +/- 0.043 (48%)
    etalfm_sn38 ~ 0.0528 # Sy 2018 Table 4, APC row 'F irinotecan->SN38': 0.0528 +/- 0.012 (23%)
    etalkel_apc ~ 0.479 # Sy 2018 Table 4, APC row 'k5e': 0.479 +/- 0.13 (69%)
    etalk12_apc ~ 0.664 # Sy 2018 Table 4, APC row 'k5p': 0.664 +/- 0.2 (81%)

    # ---- Residual error ----
    # Methods describe 'a mixture of additive and proportional error models';
    # in the event Table 4 reports a multiplicative term for every analyte and
    # an additive term for SN-38 only.
    propSd <- 0.289; label("EP proportional residual error (fraction)") # Sy 2018 Table 4, EP row 'Multiplicative error': 0.289 +/- 0.0072
    propSd_irinotecan <- 0.376; label("Irinotecan proportional residual error (fraction)") # Sy 2018 Table 4, Irinotecan row 'Multiplicative error': 0.376 +/- 0.0081
    propSd_sn38 <- 0.34; label("SN-38 proportional residual error (fraction)") # Sy 2018 Table 4, SN38 row 'Multiplicative error': 0.34 +/- 0.0089
    addSd_sn38 <- 0.383; label("SN-38 additive residual error (nmol/L)") # Sy 2018 Table 4, SN38 row 'Additive error': 0.383 +/- 0.033
    propSd_sn38g <- 0.392; label("SN-38G proportional residual error (fraction)") # Sy 2018 Table 4, SN38G row 'Multiplicative error': 0.392 +/- 0.0083
    propSd_apc <- 0.356; label("APC proportional residual error (fraction)") # Sy 2018 Table 4, APC row 'Multiplicative error': 0.356 +/- 0.0092
  })

  model({
    # ---- Individual parameters -----------------------------------------
    # EP: continuous covariates enter as a power model normalised to the
    # population median, Sy 2018 Eqs. 5 and 6.
    cl <- exp(lcl + etalcl) * (AGE / 60)^e_age_cl * (BSA / 1.86)^e_bsa_cl *
      (CRCL / 84.1)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (BSA / 1.86)^e_bsa_vc
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21)
    kel <- cl / vc
    # q and vp are algebraic restatements of k12 / k21 and carry no new
    # information. They are present so that rxode2 does not match the
    # cl / vc pair against its one-compartment analytic kernel and silently
    # discard the explicit d/dt() system below.
    q <- k12 * vc
    vp <- vc * k12 / k21

    # Irinotecan: no covariates were retained.
    vc_irinotecan <- exp(lvc_irinotecan + etalvc_irinotecan)
    kel_irinotecan <- exp(lkel_irinotecan + etalkel_irinotecan)
    k12_irinotecan <- exp(lk12_irinotecan + etalk12_irinotecan)
    k21_irinotecan <- exp(lk21_irinotecan)

    # SN-38: UGT1A1*28 homozygotes carry an exponential effect on k3e,
    # Sy 2018 Eq. 7. exp(-0.67) = 0.512, i.e. homozygotes retain about half
    # the SN-38 elimination capacity of the wild-type / heterozygous
    # reference.
    vc_sn38 <- exp(lvc_sn38 + etalvc_sn38)
    kel_sn38 <- exp(
      lkel_sn38 + etalkel_sn38 +
        e_snp_ugt1a1_rs8175347_hom_kel_sn38 * SNP_UGT1A1_RS8175347_HOM
    )
    k12_sn38 <- exp(lk12_sn38 + etalk12_sn38)
    k21_sn38 <- exp(lk21_sn38)

    # SN-38 glucuronide
    vc_sn38g <- exp(lvc_sn38g + etalvc_sn38g)
    kel_sn38g <- exp(lkel_sn38g + etalkel_sn38g)
    k12_sn38g <- exp(lk12_sn38g)
    k21_sn38g <- exp(lk21_sn38g)

    # APC shares SN-38's aggregate central volume, because Sy 2018 assumed
    # V5 = V3 and therefore V5* = V3* (Results, '5-Analyte pharmacokinetic
    # model'). With that shared volume the SN-38 : APC formation split is
    # fm_sn38 : (1 - fm_sn38), so the APC arm carries the flux ratio
    # (1 - fm_sn38) / fm_sn38 relative to the SN-38 arm. This is the factor
    # that makes F(irinotecan->SN38) identifiable; without it Table 4's
    # 0.631 +/- 0.017 would not enter the model at all.
    fm_sn38 <- exp(lfm_sn38 + etalfm_sn38)
    vc_apc <- vc_sn38
    kel_apc <- exp(lkel_apc + etalkel_apc)
    k12_apc <- exp(lk12_apc + etalk12_apc)
    k21_apc <- exp(lk21_apc)
    r_apc <- (1 - fm_sn38) / fm_sn38

    # ---- ODE system ------------------------------------------------------
    # States hold aggregate amounts: state / volume reproduces the measured
    # plasma concentration of that analyte, but only EP's state and volume
    # are on the true molar scale. Sy 2018 Eqs. 2-4 (3-analyte) and 8-9
    # (5-analyte) are written in concentration; multiplying each by its own
    # aggregate central volume gives the amount form below, in which every
    # metabolic step is a plain first-order flux out of the precursor.
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    d/dt(central_irinotecan) <- kel * central -
      (kel_irinotecan + k12_irinotecan) * central_irinotecan +
      k21_irinotecan * peripheral1_irinotecan
    d/dt(peripheral1_irinotecan) <- k12_irinotecan * central_irinotecan -
      k21_irinotecan * peripheral1_irinotecan

    d/dt(central_sn38) <- kel_irinotecan * central_irinotecan -
      (kel_sn38 + k12_sn38) * central_sn38 + k21_sn38 * peripheral1_sn38
    d/dt(peripheral1_sn38) <- k12_sn38 * central_sn38 -
      k21_sn38 * peripheral1_sn38

    # Sy 2018 Eq. 8 prints k4e in the SN-38G formation term; that is a
    # typographical error for k3e. Both the surrounding text ('Glucuronidation
    # of SN38 to SN38G is governed by the rate constant F34 k3e') and Figure 2b
    # ('k34 = k3e') give k3e, and mass balance forbids a metabolite being
    # formed at its own elimination rate constant.
    d/dt(central_sn38g) <- kel_sn38 * central_sn38 -
      (kel_sn38g + k12_sn38g) * central_sn38g + k21_sn38g * peripheral1_sn38g
    d/dt(peripheral1_sn38g) <- k12_sn38g * central_sn38g -
      k21_sn38g * peripheral1_sn38g

    d/dt(central_apc) <- r_apc * kel_irinotecan * central_irinotecan -
      (kel_apc + k12_apc) * central_apc + k21_apc * peripheral1_apc
    d/dt(peripheral1_apc) <- k12_apc * central_apc - k21_apc * peripheral1_apc

    # ---- Observations ----------------------------------------------------
    Cc <- central / vc
    Cc_irinotecan <- central_irinotecan / vc_irinotecan
    Cc_sn38 <- central_sn38 / vc_sn38
    Cc_sn38g <- central_sn38g / vc_sn38g
    Cc_apc <- central_apc / vc_apc

    Cc ~ prop(propSd)
    Cc_irinotecan ~ prop(propSd_irinotecan)
    Cc_sn38 ~ add(addSd_sn38) + prop(propSd_sn38)
    Cc_sn38g ~ prop(propSd_sn38g)
    Cc_apc ~ prop(propSd_apc)
  })
}
