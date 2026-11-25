Ueno_2018_daclatasvir_asunaprevir_svr12 <- function() {
  description <- paste0(
    "Logistic-regression exposure-efficacy model for sustained virologic ",
    "response 12 weeks after the end of treatment (SVR12) in Japanese ",
    "adults with chronic hepatitis C virus genotype 1b infection treated ",
    "for 24 weeks with daclatasvir 60 mg once daily plus asunaprevir ",
    "(100 mg softgel capsule or 200 / 600 mg tablet twice daily) ",
    "(Ueno 2018, n = 265, phase 2 AI447017 and phase 3 AI447026). The ",
    "logit of SVR12 is linear and uncentred in the average steady-state ",
    "plasma concentrations of both drugs with a product (interaction) ",
    "term and NO intercept, plus additive shifts for the baseline NS5A ",
    "Y93H resistance-associated variant and for a missing Y93H result. ",
    "There is no PK layer and no ODE: both Cav,ss values are supplied as ",
    "data, and the source took them from the companion population PK ",
    "models Osawa_2018_daclatasvir and Osawa_2018_asunaprevir."
  )
  reference <- paste(
    "Ueno T, Osawa M, Imai Y, Ishikawa H, Garimella T.",
    "Exposure-Response (Efficacy) Analysis of Daclatasvir and Asunaprevir",
    "in Japanese Patients With Hepatitis C Virus Infection.",
    "J Clin Pharmacol. 2018;58(11):1479-1488. doi:10.1002/jcph.1262.",
    "PMID 30063245, PMCID PMC6175176.",
    sep = " "
  )
  vignette <- "Ueno_2018_daclatasvir_asunaprevir"

  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV_DACLATASVIR and CAV_ASUNAPREVIR covariate columns)",
    concentration = "prob_svr12 (probability of sustained virologic response at post-treatment week 12, 0-1)"
  )

  covariateData <- list(
    CAV_DACLATASVIR = list(
      description = "Individual average steady-state plasma concentration of daclatasvir, Cav,ss = AUC(tau,ss) / tau, predicted from the companion population PK model. Supplied as data: this model has no PK layer.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Ueno 2018 Methods, 'Exposure-response end point and model development': exposures were 'predicted from previously described popPK models' and merged with the clinical data. Reproduce the column with modellib('Osawa_2018_daclatasvir'): Cav,ss = 60 mg / (individual CL/F x 24 h), x 1000 for ng/mL. Enters the logit LINEARLY and UNCENTRED (no intercept), and also through the product term with CAV_ASUNAPREVIR. Analysis-set median 557.4 ng/mL (range 148.8-1486.0) on 60 mg once daily (Table 1).",
      source_name = "Cavgss,DCV"
    ),
    CAV_ASUNAPREVIR = list(
      description = "Individual average steady-state plasma concentration of asunaprevir, Cav,ss = AUC(tau,ss) / tau, predicted from the companion population PK model. Supplied as data: this model has no PK layer.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Ueno 2018 Methods: predicted from the companion popPK model, reproducible with modellib('Osawa_2018_asunaprevir'): Cav,ss = daily dose / (individual CL/F x 24 h), x 1000 for ng/mL, with the formulation, AST and cirrhosis effects of that model. Enters the logit LINEARLY and UNCENTRED, and through the product term with CAV_DACLATASVIR. Analysis-set medians (Table 1): 163.7 ng/mL (range 59.3-947.1) on the 100 mg twice-daily softgel capsule, 268.4 (91.0-440.3) on the 200 mg twice-daily tablet and 590.2 (379.4-2111.2) on the 600 mg twice-daily tablet.",
      source_name = "Cavgss,ASV"
    ),
    HCV_NS5A_Y93H = list(
      description = "Baseline HCV NS5A Y93H resistance-associated variant indicator: 1 = Y93H detected at baseline, 0 = not detected (or result missing; see HCV_NS5A_Y93H_MISSING).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no baseline Y93H)",
      notes = "The only covariate retained in the final model (Ueno 2018 Table 2, 'Y93H'). 40 of 265 patients (15.1%) carried the variant; observed SVR12 47.5% with versus 90.3% without. Set to 0 for a patient whose Y93H result is missing and carry that patient in HCV_NS5A_Y93H_MISSING.",
      source_name = "Y93H1"
    ),
    HCV_NS5A_Y93H_MISSING = list(
      description = "Baseline HCV NS5A Y93H result missing: 1 = no baseline Y93H result, 0 = result available.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Y93H result available)",
      notes = "Ueno 2018 coded the 8 of 265 patients (3.0%) without a baseline Y93H result as a separate category whose estimate is 'the effect of missing group relative to patients without Y93H mutation, which is the reference category' (Table 2, 'Y93H Missing', RSE 217%). Describes who was left ungenotyped rather than viral biology; do not extrapolate.",
      source_name = "Y93H2"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested in the full model as (AGE - 62)/62 (Supplementary Table 3, CON2 = -0.0667, RSE 2410%) and removed by backward elimination at the 1% LRT level."
    ),
    WT_BASE = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested in the full model as (WT - 55)/55 (Supplementary Table 3, CON4 = 1.69, RSE 95.3%); not retained."
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Tested in the full model as (ALT - 55)/55 (Supplementary Table 3, CON3 = 0.620, RSE 55.5%); not retained."
    ),
    CRCL_BASE = list(
      description = "Baseline creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Tested in the full model as (CRCL - 84.6)/84.6 (Supplementary Table 3, CON5 = -1.86, RSE 57.5%); not retained."
    ),
    HCV_VLOAD = list(
      description = "Baseline HCV RNA",
      units = "IU/mL",
      type = "continuous",
      notes = "Tested in the full model on the log10 scale as (log10 HCV_VLOAD - 6.8)/6.8 (Supplementary Table 3, CON1 = -3.25, RSE 96.3%), and also dichotomised at 800,000 IU/mL; not retained in either form."
    ),
    SEXF = list(
      description = "Female sex",
      units = "(binary)",
      type = "binary",
      notes = "Tested in the full model (Supplementary Table 3, Sex = 0.0141, RSE 3910%); not retained. The supplement does not state which sex was the reference."
    ),
    DIS_CIRRHOSIS = list(
      description = "Compensated cirrhosis at baseline",
      units = "(binary)",
      type = "binary",
      notes = "Tested in the full model (Supplementary Table 3, Cirrhosis = 0.239, RSE 396%); not retained. The remaining screened covariates -- IL28B rs12979860 genotype (-0.901, RSE 54.3%), patient type (nonresponder versus IFN-ineligible naive/intolerant, -0.143, RSE 355%), study (AI447017 versus AI447026, 0.138, RSE 364%) and OATP1B1 haplotype (tested univariately, not significant) -- were likewise not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 265L,
    n_studies = 2L,
    n_observations = "265 binary SVR12 records (one per patient; landmark analysis, no repeated measures)",
    age_range = "24-75 years (median 62)",
    weight_range = "36.0-93.4 kg (median 55.0)",
    sex_female_pct = 65.7,
    race_ethnicity = "Japanese (100%)",
    disease_state = "Chronic hepatitis C virus genotype 1b infection; nonresponders (null or partial) to peginterferon alfa / ribavirin or interferon beta / ribavirin (40.8%) and interferon-ineligible naive or intolerant patients (59.2%); compensated cirrhosis 8.3%; baseline NS5A Y93H 15.1% (missing 3.0%)",
    dose_range = "Daclatasvir 60 mg tablet once daily plus asunaprevir 100 mg softgel capsule twice daily (phase 3, n = 222) or 600 / 200 mg tablet twice daily (phase 2, n = 43), for 24 weeks",
    regions = "Japan",
    notes = paste0(
      "Baseline characteristics are in Ueno 2018 Table 1. SVR12 = HCV RNA ",
      "below the limit of quantification (target detected or not detected) ",
      "at post-treatment week 12 on the Roche COBAS TaqMan assay (LLOQ ",
      "15 IU/mL), modified intent-to-treat; missing week-12 values counted as ",
      "non-response unless both the week-8 and week-24 values were below the ",
      "limit. Observed SVR12 was 90.3% without and 47.5% with baseline Y93H."
    )
  )

  ini({
    # ==================================================================
    # Ueno 2018 Table 2, 'Final Model Parameter Estimates for SVR12'.
    # The fitted logit is (Results, 'Model Development'):
    #   mu = b1 * Cav,ASV + b2 * Cav,DCV + b3 * Cav,ASV * Cav,DCV
    #        + CAT(Y93H) + CAT(Y93H missing)
    # with no intercept: the intercept was dropped because its RSE was
    # 63% and it was not significant, and the no-intercept linear model
    # with the interaction had the lowest BIC.
    #
    # b3 is printed in a row headed '(b x 10^-3)': the printed -0.00912
    # is -0.00912e-3 = -9.12e-6 per (ng/mL)^2. The paper's own
    # model-based SVR12 rates (Discussion: 91.2% / 45.5% at the softgel
    # medians, 92.7% / 50.2% at the 200 mg tablet medians, without /
    # with Y93H) reproduce only with the 10^-3 scaling; without it the
    # predicted SVR12 at the median exposures is essentially 0.
    # ==================================================================
    e_cav_asunaprevir_logit <- 0.0069; label("Log-odds of SVR12 per ng/mL of asunaprevir Cav,ss (unitless logit per ng/mL)") # Ueno 2018 Table 2, 'Slope of ASV (b1)': 0.0069 (SE 0.00236, RSE 34.2%)
    e_cav_daclatasvir_logit <- 0.00368; label("Log-odds of SVR12 per ng/mL of daclatasvir Cav,ss (unitless logit per ng/mL)") # Ueno 2018 Table 2, 'Slope of DCV (b2)': 0.00368 (SE 0.000587, RSE 16.0%)
    e_cav_asunaprevir_daclatasvir_logit <- -0.00000912; label("Log-odds of SVR12 per (ng/mL)^2 of the asunaprevir x daclatasvir Cav,ss product (unitless logit per (ng/mL)^2)") # Ueno 2018 Table 2, 'Interaction between ASV and DCV (b3) (b x 10^-3)': -0.00912 x 10^-3 (RSE 27.7%)
    e_hcv_ns5a_y93h_logit <- -2.53; label("Additive logit shift for baseline NS5A Y93H versus no Y93H (unitless logit)") # Ueno 2018 Table 2, 'Y93H (CAT1)': -2.53 (SE 0.415, RSE 16.4%)
    e_hcv_ns5a_y93h_missing_logit <- -0.508; label("Additive logit shift for a missing baseline Y93H result versus no Y93H (unitless logit)") # Ueno 2018 Table 2, 'Y93H Missing (CAT2)': -0.508 (SE 1.10, RSE 217%)

    # ==================================================================
    # The source fits a binary logistic regression in NONMEM with an
    # exact Bernoulli likelihood: no between-subject random effect and
    # no residual error. rxode2 requires an observation declaration, so
    # the deterministic probability carries a tiny placeholder additive
    # residual, as in Sangro_2023_nivolumab_ipilimumab_orr.R.
    # ==================================================================
    addSd_prob_svr12 <- fixed(0.001); label("Placeholder additive residual SD on the SVR12 probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Linear, uncentred, no-intercept logit: with both exposures at 0
    # the model returns logit 0 (probability 0.5), an extrapolation far
    # below the observed exposure range that should not be read as an
    # untreated response rate.
    logit_svr12 <- e_cav_asunaprevir_logit * CAV_ASUNAPREVIR +
      e_cav_daclatasvir_logit * CAV_DACLATASVIR +
      e_cav_asunaprevir_daclatasvir_logit * CAV_ASUNAPREVIR * CAV_DACLATASVIR +
      e_hcv_ns5a_y93h_logit * HCV_NS5A_Y93H +
      e_hcv_ns5a_y93h_missing_logit * HCV_NS5A_Y93H_MISSING
    prob_svr12 <- expit(logit_svr12)

    # Deterministic probability of SVR12. Sample binary outcomes with
    # rbinom(n, 1, prob_svr12) on the rxSolve output.
    prob_svr12 ~ add(addSd_prob_svr12)
  })
}
