Zhang_2022_ferricPyrophosphateCitrate_iv <- function() {
  description <- "One-compartment population PK model for ferric pyrophosphate citrate (FPC, Triferic) given as a 4-hour intravenous infusion to healthy Asian and non-Asian adults (Zhang 2022, model M1). Output Cc is the drug-derived (baseline-corrected) serum total iron concentration. Volume of distribution carries a power effect of lean body mass, a linear effect of the 6-hour pre-baseline average serum total iron, and a female-sex multiplier of 2.8."
  reference <- "Zhang L, Gan L, Li K, Xie P, Tan Y, Wei G, Yuan X, Pratt R, Zhou Y, Hui AM, Fang Y, Zuo L, Zheng Q. Ethnicity evaluation of ferric pyrophosphate citrate among Asian and Non-Asian populations: a population pharmacokinetics analysis. Eur J Clin Pharmacol. 2022;78(9):1421-1434. doi:10.1007/s00228-022-03328-9"
  vignette <- "Zhang_2022_ferricPyrophosphateCitrate"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    LBM = list(
      description = "Lean body mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline value. Power effect on Vd normalized to the M1 median 55.24 kg (Table 2 footnote a). The paper does not name the LBM formula; M1 cohort mean (SD) by study 51.8 (5.0), 59.1 (8.8) and 54.0 (6.8) kg, overall range 37.6-70.4 kg (Supplementary Table 4).",
      source_name = "LBM"
    ),
    IRON_BL = list(
      description = "Baseline serum total iron concentration, here the average over the 6-hour baseline window before the dosing day",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed. The paper's Fe.av: the average serum total iron over the 6 hours of the drug-free baseline-period sampling that precedes the administration day. Healthy subjects show a large 24-h fluctuation of serum iron, so the authors used this average rather than the single pre-dose value used for the CKD-5HD models. Linear effect on Vd centered at the M1 median 1110 ng/mL (Table 2 footnote a). M1 cohort mean (SD) by study 946.4 (221.4), 1170.8 (352.8) and 1425.6 (554.7) ng/mL, range 490.7-2828.6 ng/mL (Supplementary Table 4).",
      source_name = "Fe.av"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Table 2 footnote a: COVsex = 1 for males and 1 + theta_sex = 2.8 for females, multiplying Vd. Only 5 of 40 M1 subjects were female; the authors caution against drawing conclusions about sex from this estimate (Discussion).",
      source_name = "Sex"
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "iron (from ferric pyrophosphate citrate)",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 40L,
    n_studies = 3L,
    age_range = "19-62 years",
    age_median = "not reported (study means 30.8, 38.4 and 39.1 years)",
    weight_range = "53.3-104.9 kg",
    weight_median = "not reported (study means 68.1, 84.5 and 73.7 kg)",
    sex_female_pct = 12.5,
    race_ethnicity = c(Asian = 35, `Non-Asian` = 65),
    disease_state = "Healthy adult volunteers",
    dose_range = "Single 4-hour IV infusion of FPC: 6.5 mg (CHN-FPC-14), 6 mg (USA-FPC-12) or 6.6 mg (USA-FPC-18) iron",
    regions = "China; United States",
    notes = "Model M1 of Zhang 2022. 14 healthy Chinese subjects from CHN-FPC-14 (CTR20181113) and 26 healthy non-Asian subjects (African American and Caucasian) from USA-FPC-12 (NCT02636049; 12) and USA-FPC-18 (14). Demographics from Table 1 and Supplementary Table 4; dosing from Supplementary Tables 1-2."
  )

  ini({
    lcl <- log(0.477)
    label("Clearance CL (L/h)") # Table 2, M1 'CL, L/h' = 0.477 (RSE 8.2%)
    lvc <- log(3.62)
    label("Volume of distribution Vd (L) for a male with LBM 55.24 kg and Fe.av 1110 ng/mL") # Table 2, M1 'Vd, L' = 3.62 (RSE 7.0%)

    e_lbm_vc <- 3.26
    label("Power exponent of LBM/55.24 on Vd (unitless)") # Table 2, M1 'LBM on Vd' = 3.26 (RSE 18.7%); footnote a
    e_iron_bl_vc <- 0.000771
    label("Linear slope of (Fe.av - 1110) on Vd (per ng/mL)") # Table 2, M1 'Fe.av on Vd (x10^-4)' = 7.71 (RSE 17.5%); footnote a prints 0.000771
    e_sexf_vc <- 1.80
    label("Fractional increase in Vd for females (unitless; female Vd = 2.8 x male)") # Table 2, M1 'theta sex on Vd' = 1.80 (RSE 30.1%); footnote a COVfemale = 1 + theta = 2.8

    # IIV: Table 2 'omega (CL, %)' and 'omega (Vd, %)' read as 100 x the SD of
    # eta (variance = (P/100)^2); the printed 95% CI for the M1 CL row is on
    # the same SD-fraction scale (0.235, 0.621).
    etalcl ~ 0.183184 # Table 2, M1 omega1 (CL) = 42.8% (RSE 23%); 0.428^2
    etalvc ~ 0.114244 # Table 2, M1 omega2 (Vd) = 33.8% (RSE 16.7%); 0.338^2

    addSd <- 174
    label("Additive residual error (ng/mL)") # Table 2, M1 'sigma (add), ng/ml' = 174 (RSE 8.9%)
  })

  model({
    # Table 2 footnote a:
    # Vind = VTV x (1 + 0.000771 x [Fe.av - 1110]) x (LBM/55.24)^3.26 x COVsex
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) *
      (1 + e_iron_bl_vc * (IRON_BL - 1110)) *
      (LBM / 55.24)^e_lbm_vc *
      (1 + e_sexf_vc * SEXF)

    kel <- cl / vc
    d / dt(central) <- -kel * central

    # Dose in mg and volume in L give mg/L; x 1000 converts to ng/mL.
    Cc <- central / vc * 1000
    Cc ~ add(addSd)
  })
}
