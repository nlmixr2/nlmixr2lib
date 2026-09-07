Bai_2024_tacrolimus <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption and linear elimination for oral tacrolimus in Chinese adult liver transplant recipients (Bai 2024); direct bilirubin is a power-form covariate on apparent clearance and body weight a power-form covariate on apparent central volume."
  reference <- "Bai H, Yun J, Wang Z, Ma Y, Liu W. Population pharmacokinetics study of tacrolimus in liver transplant recipients: a comparison between patients with or without liver cancer before surgery. Front Pharmacol. 2024;15:1449535. doi:10.3389/fphar.2024.1449535"
  vignette <- "Bai_2024_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot   = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    DBIL = list(
      description        = "Direct (conjugated) serum bilirubin concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on apparent oral clearance: CL_typical = 37.6 * (DBIL / 22.1)^-0.188. The reference 22.1 umol/L is the cohort median direct bilirubin (Bai 2024 Table 1: median 22.1, mean 36.2 (SD 49.1), range 3.40-345 umol/L) and is the value printed inside the final-model equation in Section 3.2. The negative exponent means that rising direct bilirubin (impaired biliary excretion / cholestasis after transplantation) lowers apparent clearance.",
      source_name        = "DBIL"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on apparent central volume: Vc_typical = 1710 * (WT / 70)^1.4. The reference 70 kg is the cohort median weight (Bai 2024 Table 1: median 70.0, mean 70.7 (SD 13.9), range 41-129 kg) and is the value printed inside the final-model equation in Section 3.2. The exponent 1.4 was estimated, not fixed at an allometric 1.0.",
      source_name        = "WT"
    )
  )

  # Bai 2024 Section 2.4 screened 19 candidate covariates and carried
  # only two into the final model (WT on Vc/F, DBIL on CL/F). The
  # remaining screened covariates are recorded here for provenance; the
  # paper reports no coefficient for any of them, so none can be
  # implemented. Two of the 19 have no canonical register column and are
  # therefore described in population$notes instead of listed here: the
  # preoperative disease diagnosis (liver cancer vs non-liver-cancer,
  # the paper's primary comparison) and prothrombin time in seconds.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units       = "unitless",
      type        = "binary",
      notes       = "Screened but not retained by the stepwise forward-inclusion / backward-elimination procedure (Bai 2024 Section 2.4). Cohort 148 male / 48 female (Table 1). No coefficient is reported."
    ),
    AGE = list(
      description = "Subject age",
      units       = "year",
      type        = "continuous",
      notes       = "Screened but not retained (Bai 2024 Section 2.4). Cohort median 55.0, range 16.0-76.0 years (Table 1). No coefficient is reported."
    ),
    POD = list(
      description = "Post-operative day (days since liver transplantation)",
      units       = "days",
      type        = "continuous",
      notes       = "Screened but not retained (Bai 2024 Section 2.4). Cohort median 5.50, range 0.500-30.9 days (Table 1). The Discussion attributes the null result to the short early-postoperative follow-up window and calls for longer follow-up. No coefficient is reported."
    ),
    WBC = list(
      description = "White blood cell count",
      units       = "10^9/L",
      type        = "continuous",
      notes       = "Screened but not retained (Bai 2024 Section 2.4). Cohort median 6.43, range 1.44-22.2 10^9/L (Table 1). No coefficient is reported."
    ),
    LYMPH_ABS = list(
      description = "Absolute peripheral-blood lymphocyte count",
      units       = "10^9/L",
      type        = "continuous",
      notes       = "Screened but not retained (Bai 2024 Section 2.4). Cohort median 0.545, range 0.0500-6.45 10^9/L (Table 1). No coefficient is reported."
    ),
    HCT = list(
      description = "Haematocrit",
      units       = "%",
      type        = "continuous",
      notes       = "Screened but not retained (Bai 2024 Section 2.4). Cohort median 26.3%, range 13.5-40.6% (Table 1). Tacrolimus partitions extensively into erythrocytes (85-95% of whole-blood drug) and haematocrit is a commonly retained covariate in other tacrolimus models, so the null result here is notable. No coefficient is reported."
    ),
    HGB = list(
      description = "Haemoglobin concentration",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained (Bai 2024 Section 2.4). Cohort median 87.5, range 48.0-141 g/L (Table 1). No coefficient is reported."
    ),
    ALB = list(
      description = "Serum albumin concentration",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained (Bai 2024 Section 2.4). Cohort median 36.0, range 28.4-57.0 g/L (Table 1). No coefficient is reported."
    ),
    TPRO = list(
      description = "Total serum protein concentration",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained (Bai 2024 Section 2.4). Cohort median 59.3, range 40.4-139 g/L (Table 1). No coefficient is reported."
    ),
    AST = list(
      description = "Aspartate aminotransferase activity",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained (Bai 2024 Section 2.4). Cohort median 39.0, range 11.0-5970 U/L (Table 1). No coefficient is reported."
    ),
    ALT = list(
      description = "Alanine aminotransferase activity",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained (Bai 2024 Section 2.4). Cohort median 117, range 25.0-1740 U/L (Table 1). No coefficient is reported."
    ),
    TBILI = list(
      description = "Total serum bilirubin concentration",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained (Bai 2024 Section 2.4); the direct (conjugated) fraction DBIL was retained instead. Cohort median 38.3, range 0.290-457 umol/L (Table 1). No coefficient is reported for total bilirubin."
    ),
    ALP = list(
      description = "Alkaline phosphatase activity",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained (Bai 2024 Section 2.4). Cohort median 133, range 51.0-570 U/L (Table 1). No coefficient is reported."
    ),
    BUN = list(
      description = "Serum urea concentration",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened but not retained (Bai 2024 Section 2.4). Bai 2024 Table 1 reports serum urea in SI mmol/L (median 11.9, range 4.07-48.1), not urea nitrogen in mg/dL; the two differ by the factor BUN_mgdL = urea_mmolL * 2.8. No coefficient is reported."
    ),
    CREAT = list(
      description = "Serum creatinine concentration",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained (Bai 2024 Section 2.4). Cohort median 60.5, range 28.0-291 umol/L (Table 1). No coefficient is reported."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 196L,
    n_studies      = 1L,
    n_observations = 802L,
    age_range      = "16-76 years",
    age_median     = "55.0 years",
    age_mean       = "54.4 +/- 10.1 years",
    weight_range   = "41-129 kg",
    weight_median  = "70.0 kg",
    weight_mean    = "70.7 +/- 13.9 kg",
    sex_female_pct = 24.5,
    race_ethnicity = c(Chinese = 100),
    disease_state  = "Adults in the early period after a first orthotopic liver transplantation, sampled from the first postoperative day to hospital discharge (postoperative day median 5.50, range 0.500-30.9 days). The cohort is split by the preoperative diagnosis into a liver-cancer group (N = 118) and a non-liver-cancer group (N = 78); the final model does not distinguish the two. Post-transplant hepatic recovery spans a wide range (direct bilirubin median 22.1, range 3.40-345 umol/L; ALT median 117, range 25.0-1740 U/L) and the cohort is anaemic (haematocrit median 26.3%, haemoglobin median 87.5 g/L).",
    dose_range     = "Oral tacrolimus (Prograf) twice daily on an empty stomach at 08:00 and 20:00, started 6-48 h after transplantation and titrated to trough concentration; daily dose mean 3.90 +/- 0.84 mg, median 3.92 mg, range 1.40-6.00 mg. Co-administered with mycophenolate mofetil and methylprednisolone.",
    regions        = "Single center: Beijing YouAn Hospital of Capital Medical University, Beijing, China.",
    notes          = "Retrospective therapeutic-drug-monitoring analysis, November 2021 - December 2023. Baseline demographics per Bai 2024 Table 1; final population PK parameter estimates and bootstrap validation per Bai 2024 Table 3. Whole-blood tacrolimus by validated LC-MS over a 1.0-50.0 ng/mL calibrated range; observed concentrations mean 4.53 +/- 3.30 ng/mL, median 3.80, range BQL-28.8 ng/mL. Sampling was overwhelmingly pre-morning-dose trough (0.5-1 h before the 08:00 dose), which the authors note makes Ka poorly identified. CYP3A5 genotype was not available and was not tested as a covariate (stated limitation). Section 2.4 lists the full 19-covariate screen: preoperative disease diagnosis, gender, age, weight, postoperative days, white blood cell count, lymphocytes, haematocrit, haemoglobin, albumin, total protein, AST, ALT, total bilirubin, direct bilirubin, alkaline phosphatase, urea, creatinine, and prothrombin time. Only weight (on Vc/F) and direct bilirubin (on CL/F) survived stepwise selection; the other 17 carry no reported coefficient and are recorded in covariatesDataExcluded, except the preoperative diagnosis and prothrombin time (seconds; cohort median 11.5, range 7.90-27.1 s), which have no canonical register column. Bai 2024 Table 2 reports a rejected base model in which CL/F, Vc/F and Ka were estimated separately in the liver-cancer and non-liver-cancer groups (CL/F 38.9 vs 37.7 L/h, Vc/F 1603.6 vs 1772.2 L, Ka 0.354 vs 0.221 1/h); despite a 7.50-point OFV drop the authors rejected it because the parameters were similar and IIV_Ka was not estimable in the non-liver-cancer group (RSE 4959%), so only the pooled final model is implemented here."
  )

  ini({
    # Structural parameters. The typical-value reference subject has
    # direct bilirubin 22.1 umol/L and body weight 70 kg (the cohort
    # medians printed inside the final-model equation, Bai 2024
    # Section 3.2):
    #   CL_i  (L/h) = 37.6 * exp(eta_CL_i)  * (DBIL / 22.1)^-0.188
    #   Vc_i  (L)   = 1710 * exp(eta_VC_i)  * (WT / 70)^1.4
    #   Ka    (1/h) = 0.352 * exp(eta_KA_i)
    # Both CL and Vc are apparent (CL/F, Vc/F): tacrolimus was given
    # orally only, so bioavailability is not separately identifiable and
    # no F parameter is estimated.
    lka <- log(0.352); label("Absorption rate constant (1/h)")                                # Bai 2024 Table 3 final model: Ka = 0.352 1/h (95% CI 0.15-0.826)
    lcl <- log(37.6);  label("Apparent oral clearance at DBIL = 22.1 umol/L (L/h)")           # Bai 2024 Table 3 final model: CL/F = 37.6 L/h (95% CI 34-41.7)
    lvc <- log(1710);  label("Apparent central volume at WT = 70 kg (L)")                     # Bai 2024 Table 3 final model: Vc/F = 1710 L (95% CI 1400-2100)

    # Covariate effects, both power-form on the median-normalised
    # covariate. Bai 2024 Table 3 labels these two rows as -Effect of WT
    # on CL/F- (-0.188) and -Effect of DBIL on Vc/F- (1.4), which
    # transposes the two covariates relative to every other statement in
    # the paper. The printed final-model equation in Section 3.2, the
    # Results text -the final model included weight as the covariate of
    # Vc and DBIL as the covariate of CL-, and the Discussion and
    # Conclusion all agree that -0.188 belongs to DBIL on CL and 1.4 to
    # WT on Vc. The reference constants settle it independently: 22.1 is
    # the Table 1 DBIL median and 70 the Table 1 weight median, and each
    # appears in the equation beside its own covariate. The equation is
    # followed here and the Table 3 row labels are treated as a
    # typesetting transposition; see the vignette Errata.
    e_dbil_cl <- -0.188; label("Power exponent for direct bilirubin on CL/F (unitless)")      # Bai 2024 Section 3.2 final-model equation (95% CI -0.315 to -0.0616 per Table 3)
    e_wt_vc   <-  1.4;   label("Power exponent for body weight on Vc/F (unitless)")           # Bai 2024 Section 3.2 final-model equation (95% CI 0.512-2.28 per Table 3)

    # IIV. Exponential on all three structural parameters, with CL and Vc
    # correlated (NONMEM OMEGA BLOCK(2)); Ka is diagonal. Bai 2024
    # Table 3 reports the diagonals as percentages (66.6, 126 and 236)
    # and the off-diagonal as a raw covariance (0.226).
    #
    # The percentages are log-normal CV, i.e. omega^2 = log(1 + CV^2),
    # not the NONMEM approximation omega = CV. Table 2 pins this: it
    # reports the same three IIVs for the preceding base model as raw
    # variances (0.314 and 0.528 for CL, 0.993 and 0.934 for Vc, 3.95
    # and 0.025 for Ka in the liver-cancer and non-liver-cancer groups
    # respectively). Under the log-normal reading the final variances
    # are 0.367 for CL, 0.951 for Vc and 1.882 for Ka -- each at or just
    # below the corresponding base-model value, which is what adding a
    # covariate should do. Under omega = CV they would be 0.444, 1.588
    # and 5.570, i.e. Vc and Ka would have gained 60% and 130% more
    # unexplained variance on adding a covariate, which is not possible.
    #
    # omega^2 CL = log(1 + 0.666^2) = 0.367110
    # omega^2 Vc = log(1 + 1.26^2)  = 0.950731
    # omega^2 Ka = log(1 + 2.36^2)  = 1.882453
    # implied CL-Vc correlation = 0.226 / sqrt(0.367110 * 0.950731) = 0.383
    etalcl + etalvc ~ c(0.367110,
                        0.226, 0.950731)                                                     # Bai 2024 Table 3 final model: IIV_CL/F 66.6% CV, IIV_Vc/F 126% CV, covariance 0.226
    etalka ~ 1.882453                                                                        # Bai 2024 Table 3 final model: IIV_Ka 236% CV (95% CI upper bound 1070%; lower bound not estimable)

    # Combined proportional-plus-additive residual error on the
    # whole-blood concentration scale. The additive term is reported in
    # ug/L, which is identical to the ng/mL used everywhere else in the
    # paper.
    propSd <- 0.356; label("Proportional residual error (fraction)")                          # Bai 2024 Table 3 final model: proportional residual variability 35.6% (95% CI 32-38.9)
    addSd  <- 0.355; label("Additive residual error (ng/mL)")                                 # Bai 2024 Table 3 final model: additive residual variability 0.355 ug/L (95% CI 0.166-0.473)
  })

  model({
    # Reference (median) covariate values printed inside the final-model
    # equation, Bai 2024 Section 3.2.
    ref_dbil <- 22.1
    ref_wt   <- 70

    # Individual parameters. Power-form covariate on each of CL and Vc.
    dbil_cl <- (DBIL / ref_dbil)^e_dbil_cl
    wt_vc   <- (WT / ref_wt)^e_wt_vc

    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * dbil_cl
    vc <- exp(lvc + etalvc) * wt_vc

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Amounts are in mg and vc is in L, so central / vc is mg/L. The
    # paper reports whole-blood tacrolimus in ng/mL (= ug/L), so scale by
    # 1000 mg/L per ug/L.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
