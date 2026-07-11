Overgaard_2016_liraglutide <- function() {
  description <- "Liraglutide 3.0 mg population PK model in overweight and obese adults with and without type 2 diabetes (Overgaard 2016 SCALE Obesity/Prediabetes + SCALE Diabetes pooled analysis)"
  reference <- "Overgaard RV, Petri KC, Jacobsen LV, Jensen CB. Liraglutide 3.0 mg for Weight Management: A Population Pharmacokinetic Analysis. Clin Pharmacokinet. 2016;55(11):1413-1422. doi:10.1007/s40262-016-0410-7"
  vignette <- "Overgaard_2016_liraglutide"
  units <- list(time = "hr", dosing = "mg", concentration = "nmol/L")
  # Time in hours; SC doses in mg; observed concentrations in nmol/L (nM).
  # Liraglutide MW = 3751.2 g/mol; the vignette uses the MW to convert
  # simulated mg/L concentrations to nM before reporting AUC24 in nM*h to
  # match Overgaard 2016 Fig 3 / Table S1 units.

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL/F only; reference weight 100 kg per Overgaard 2016 Methods (Sect. 2.2, 'rounded value close to mean body weight'). Empirical range 60-234 kg per Table 1.",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female) is the paper's reference category",
      notes              = "Overgaard 2016 Table S1 reports the log-scale coefficient Cov.male = +0.27 (E_male = exp(0.27) = 1.310). Implemented as exp(e_male_cl * (1 - SEXF)) so SEXF = 1 (female) evaluates to 1.00 (reference) and SEXF = 0 (male) evaluates to 1.310, matching the paper's female reference. Cov.male sign inverts when the source column is SEXM: keep the canonical SEXF encoding and let the model do the inversion.",
      source_name        = "SEXM"
    ),
    AGE_GE70 = list(
      description        = "Indicator for baseline age >= 70 years",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (age < 70 years)",
      notes              = "Overgaard 2016 Methods Sect. 2.2 defines the paper's age band as < or >= 70 years (final analysis; the pre-specified 75-year cut-off was lowered to 70 y for adequate cell sizes per Online Resource Sect. S2). Only 2.5 % of subjects (73 of 2923) fell in the >= 70 y stratum.",
      source_name        = "AGE_GE70"
    ),
    RACE_BLACK = list(
      description        = "1 = Black / African American, 0 = other",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White is the paper's reference; the sibling indicators RACE_ASIAN and RACE_OTHER are also 0 for White subjects)",
      notes              = "Overgaard 2016 Table 1 reports 9.7 % Black / African American across the pooled pharmacokinetic dataset (283 of 2923).",
      source_name        = "RACE_BLACK"
    ),
    RACE_ASIAN = list(
      description        = "1 = Asian, 0 = other",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White reference; small cohort, 3.3 % of pooled pharmacokinetic dataset)",
      notes              = "Overgaard 2016 Table 1 reports 3.3 % Asian across the pooled pharmacokinetic dataset (96 of 2923). Coefficient reported at very high RSE (913 %; Table S1) because of small cell size; the paper's own analysis concluded the Asian effect is not pharmacokinetically relevant.",
      source_name        = "RACE_ASIAN"
    ),
    RACE_OTHER = list(
      description        = "1 = race category 'Other', 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White reference)",
      notes              = "Overgaard 2016 Table 1 footnote a: 'Other' pools American Indian / Alaskan Native, Native Hawaiian / other Pacific Islander, and any 'Other'. 2.1 % of pooled pharmacokinetic dataset (62 of 2923).",
      source_name        = "RACE_OTHER"
    ),
    RACE_HISPANIC = list(
      description        = "1 = Hispanic / Latino ethnicity, 0 = non-Hispanic",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Hispanic / -Latino reference)",
      notes              = "Overgaard 2016 reports Hispanic ethnicity as a covariate dimension distinct from race (paper's Methods list: 'race (White, Black or African American, Asian, Other), ethnicity (Hispanic or Latino, non-Hispanic or -Latino)'). The numerical column encoding (1 = Hispanic vs 0 = non-Hispanic) is identical to the canonical RACE_HISPANIC form. 10.5 % of pooled pharmacokinetic dataset (306 of 2923). See RACE_HISPANIC register entry for the ethnicity-vs-race semantic and the OMB-classification note.",
      source_name        = "RACE_HISPANIC"
    ),
    DIS_PREDIAB = list(
      description        = "Indicator for baseline prediabetes glycaemic status",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not prediabetic; either normoglycaemic reference or T2DM which is carried by DIS_DIAB)",
      notes              = "Overgaard 2016 stratifies baseline glycaemic status into three levels: normoglycaemic reference, prediabetic (DIS_PREDIAB = 1), and T2DM (DIS_DIAB = 1). Prediabetic subjects are 49.3 % of the pooled pharmacokinetic dataset (1442 of 2923), all from Trial 1 (SCALE Obesity and Prediabetes).",
      source_name        = "DIS_PREDIAB"
    ),
    DIS_DIAB = list(
      description        = "Indicator for baseline type 2 diabetes mellitus",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normoglycaemic reference; DIS_PREDIAB carries the prediabetic level)",
      notes              = "Overgaard 2016 T2DM subjects are 20.0 % of the pooled pharmacokinetic dataset (584 of 2923), all from Trial 2 (SCALE Diabetes). The paper notes the T2DM effect is confounded with Trial 2 (Sect. 3.3.5); the 90 % CI narrowly falls within the bioequivalence limits so the paper concludes the effect is not pharmacokinetically relevant.",
      source_name        = "DIS_DIAB"
    ),
    DOSE_1P8MG = list(
      description        = "Indicator for the 1.8 mg liraglutide dose arm",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (3.0 mg liraglutide reference)",
      notes              = "Overgaard 2016 Trial 2 (SCALE Diabetes) enrolled 191 subjects on 1.8 mg liraglutide alongside 393 subjects on 3.0 mg. The 1.8 mg arm is 6.5 % of the pooled pharmacokinetic dataset (191 of 2923). The paper uses this indicator to formally test dose linearity across the SCALE-Diabetes 1.8-3.0 mg range; the reported log-scale coefficient (+0.02, 95 % CI -0.03 to 0.08) confirms dose proportionality.",
      source_name        = "DOSE_1P8MG"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 2923L, # Overgaard 2016 Table 1 pooled Trials 1 + 2 (2339 + 584)
    n_studies      = 2L,    # Trials 1 and 2 (Trial 3 is separately used for dose-proportionality only)
    age_range      = "18-84 years (Trial 1) and 18-79 years (Trial 2); mean 47.1 y, SD 12.3 y; 2.5 % (73 / 2923) aged >= 70 y",
    weight_range   = "60-234 kg overall; mean 106 kg, SD 21 kg (Table 1)",
    sex_female_pct = 72.3, # 2112 / 2923 female per Overgaard 2016 Table 1
    race_ethnicity = c(White = 84.9, `Black or African American` = 9.7, Asian = 3.3, Other = 2.1, Hispanic = 10.5), # Overgaard 2016 Table 1
    disease_state  = "Overweight and obese adults with or without T2DM. 30.7 % normoglycaemic, 49.3 % prediabetic, 20.0 % T2DM. Mean BMI 38 kg/m^2, SD 7 kg/m^2 (Table 1).",
    dose_range     = "1.8 mg once daily SC (191 subjects, Trial 2 only) or 3.0 mg once daily SC (2732 subjects, Trials 1 + 2), each preceded by a 4-week 0.6 mg/week dose escalation.",
    regions        = "Multinational SCALE Obesity and Prediabetes and SCALE Diabetes phase IIIa trials (Trial 1: NN8022-1839 [Pi-Sunyer 2015]; Trial 2: NN8022-1922 [Davies 2015]).",
    trials         = c("NN8022-1839", "NN8022-1922"),
    renal_function = "Normal (51 %), mild impairment (44 %), moderate (5 %), severe (< 0.1 %) per CKD-EPI eGFR bands (Overgaard 2016 Table 1).",
    notes          = "Cohort is the SCALE Obesity + Prediabetes + Diabetes pooled analysis at N = 2923 subjects (2339 + 584 after exclusion of records with incomplete dosing history). Data from Trial 3 (Astrup 2009 dose-finding, N = 331) were used only for the dose-proportionality post-hoc using post hoc Bayes estimates with parameters fixed; the population PK model itself was fit to Trials 1 + 2. Overgaard 2016 Sect. 2.1.2."
  )

  ini({
    lka <- fixed(log(0.0806)); label("Absorption rate constant (1/h)") # Overgaard 2016 Methods Sect. 2.2 (fixed to a value estimated from a prior obese-subject multiple-dose popPK, data on file per ref [8]). Reported as 0.0806 h^-1 in the paper text (0.09 h^-1 rounded in Table S1).
    lcl <- log(0.86); label("Apparent clearance (L/h) for the reference subject") # Overgaard 2016 Online Resource Table S1 (RSE 2 %, 95 % CI 0.83-0.90). Reference: female, 100 kg, White, non-Hispanic, non-diabetic, < 70 y, on 3.0 mg once daily.
    lvc <- log(24.60); label("Apparent central volume of distribution (L)") # Overgaard 2016 Online Resource Table S1 (RSE 9 %, 95 % CI 20.3-28.8).

    e_wt_cl <- 0.68; label("Body weight power exponent on CL/F") # Overgaard 2016 Online Resource Table S1 "Cov. weight" (RSE 5 %, 95 % CI 0.61-0.75); applied as (WT/100)^e_wt_cl.
    e_male_cl <- 0.27; label("Log-scale male sex effect on CL/F") # Overgaard 2016 Online Resource Table S1 "Cov. male" (RSE 6 %, 95 % CI 0.24-0.30); applied as exp(e_male_cl * (1 - SEXF)); E_male = exp(0.27) = 1.310.
    e_age_ge70_cl <- -0.10; label("Log-scale age >=70 y effect on CL/F") # Overgaard 2016 Online Resource Table S1 "Cov. age >=70 years" (RSE 45 %, 95 % CI -0.18 to -0.01); E_age_ge70 = exp(-0.10) = 0.905.
    e_race_black_cl <- -0.09; label("Log-scale race-Black effect on CL/F") # Overgaard 2016 Online Resource Table S1 "Cov. Black" (RSE 26 %, 95 % CI -0.13 to -0.04); E_black = exp(-0.09) = 0.914.
    e_race_asian_cl <- -0.001; label("Log-scale race-Asian effect on CL/F") # Overgaard 2016 Online Resource Table S1 "Cov. Asian" (RSE 913 %, 95 % CI -0.09 to 0.08); E_asian = exp(-0.001) approximately 0.999 (paper's own analysis finds this effect not pharmacokinetically relevant).
    e_race_other_cl <- -0.08; label("Log-scale race-Other effect on CL/F") # Overgaard 2016 Online Resource Table S1 "Cov. Other" (RSE 57 %, 95 % CI -0.17 to -0.01); E_other = exp(-0.08) = 0.923.
    e_race_hispanic_cl <- 0.08; label("Log-scale ethnicity-Hispanic effect on CL/F") # Overgaard 2016 Online Resource Table S1 "Cov. Hispanic" (RSE 26 %, 95 % CI 0.04-0.12); E_hispanic = exp(0.08) = 1.083.
    e_dis_prediab_cl <- 0.00; label("Log-scale prediabetes effect on CL/F") # Overgaard 2016 Online Resource Table S1 "Cov. prediabetes" (RSE 904 %, 95 % CI -0.05 to 0.06); effectively no effect.
    e_dis_diab_cl <- 0.18; label("Log-scale T2DM effect on CL/F") # Overgaard 2016 Online Resource Table S1 "Cov. diabetes" (RSE 14 %, 95 % CI 0.13-0.23); E_diab = exp(0.18) = 1.197. Confounded with Trial 2 per Sect. 3.3.5.
    e_dose_1p8mg_cl <- 0.02; label("Log-scale 1.8 mg dose effect on CL/F") # Overgaard 2016 Online Resource Table S1 "Cov. 1.8 mg" (RSE 119 %, 95 % CI -0.03 to 0.08); confirms dose proportionality across the 1.8-3.0 mg range.

    # IIV reported as %CV on the log-normal random effect; convert with omega^2 = log(1 + CV^2).
    etalcl ~ log(1 + 0.247^2) # Overgaard 2016 Online Resource Table S1: IIV-CL/F = 24.70 %CV (shrinkage 23.9 %).
    etalvc ~ log(1 + 0.347^2) # Overgaard 2016 Online Resource Table S1: IIV-V/F = 34.70 %CV. High shrinkage (83.2 %) noted in Section S1: individual V/F estimates are relatively inaccurate under sparse sampling; retained per faithful implementation.
    propSd <- 0.154; label("Proportional residual error (fraction)") # Overgaard 2016 Online Resource Table S1: Sigma = 15.40 %CV (shrinkage 9.37 %).
  })
  model({
    ka <- exp(lka)
    cl_wt <- (WT / 100)^e_wt_cl # Overgaard 2016 Table S1: power effect on CL/F, reference 100 kg per Methods Sect. 2.2.
    cl_covs <- exp(
      e_male_cl       * (1 - SEXF) +   # SEXF = 1 (female) is the reference; (1 - SEXF) is the male indicator
      e_age_ge70_cl   * AGE_GE70 +
      e_race_black_cl * RACE_BLACK +
      e_race_asian_cl * RACE_ASIAN +
      e_race_other_cl * RACE_OTHER +
      e_race_hispanic_cl * RACE_HISPANIC +
      e_dis_prediab_cl * DIS_PREDIAB +
      e_dis_diab_cl    * DIS_DIAB +
      e_dose_1p8mg_cl  * DOSE_1P8MG
    ) # Overgaard 2016 Methods Sect. 2.2 CL/F equation: CL/F_i = TVCL * E_weight * E_dose * E_sex * E_age * E_ethnicity * E_disease_status * exp(eta_i). Race (E_race) is treated in the printed equation as part of the covariate list but not spelled out as a factor -- restored here from Table S1 which lists Cov. Black / Cov. Asian / Cov. Other as separate log-scale coefficients.
    cl <- exp(lcl + etalcl) * cl_wt * cl_covs
    vc <- exp(lvc + etalvc)

    Cc <- linCmt()
    Cc ~ prop(propSd)
  })
}
