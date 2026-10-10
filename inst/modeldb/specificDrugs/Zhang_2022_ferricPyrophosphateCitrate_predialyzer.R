Zhang_2022_ferricPyrophosphateCitrate_predialyzer <- function() {
  description <- "One-compartment population PK model for ferric pyrophosphate citrate (FPC, Triferic) infused over 3 hours into the pre-dialyzer blood circuit of Asian and non-Asian adults with haemodialysis-dependent stage 5 chronic kidney disease (Zhang 2022, model M3). Output Cc is the drug-derived (baseline-corrected) serum total iron concentration. Apparent clearance CL/F decreases linearly with pre-dose serum total iron; apparent volume Vd/F carries a power effect of lean body mass."
  reference <- "Zhang L, Gan L, Li K, Xie P, Tan Y, Wei G, Yuan X, Pratt R, Zhou Y, Hui AM, Fang Y, Zuo L, Zheng Q. Ethnicity evaluation of ferric pyrophosphate citrate among Asian and Non-Asian populations: a population pharmacokinetics analysis. Eur J Clin Pharmacol. 2022;78(9):1421-1434. doi:10.1007/s00228-022-03328-9"
  vignette <- "Zhang_2022_ferricPyrophosphateCitrate"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    LBM = list(
      description = "Lean body mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline value. Power effect on Vd/F normalized to the M3 median 55.85 kg (Table 2 footnote c). The paper does not name the LBM formula; M3 cohort mean (SD) by study 50.6 (6.6), 62.4 (8.5) and 57.1 (5.7) kg, overall range 37.8-76.8 kg (Supplementary Table 4).",
      source_name = "LBM"
    ),
    IRON_BL = list(
      description = "Baseline serum total iron concentration, here the single value measured at 0 h before administration",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed per administration. The paper's Fe baseline: serum total iron at 0 h before the FPC dose. Linear effect on CL/F centered at the M3 median 660.80 ng/mL (Table 2 footnote c). M3 cohort mean (SD) by study 794.7 (223.0), 594.6 (206.2) and 656.2 (246.5) ng/mL, range 230-1320 ng/mL (Supplementary Table 4).",
      source_name = "Fe baseline"
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
    n_subjects = 51L,
    n_studies = 3L,
    age_range = "25-77 years",
    age_median = "not reported (study means 54.3, 49.2 and 53.7 years)",
    weight_range = "53-128.4 kg",
    weight_median = "not reported (study means 69.3, 98.0 and 84.3 kg)",
    sex_female_pct = 19.6,
    race_ethnicity = c(Asian = 23.5, `Non-Asian` = 76.5),
    disease_state = "Haemodialysis-dependent stage 5 chronic kidney disease (CKD-5HD)",
    dose_range = "6.5 mg (CHN-FPC-21, USA-FPC-20) or 6.6 mg (USA-FPC-16) iron as FPC infused into the pre-dialyzer blood circuit over the first 3 hours of a dialysis session",
    regions = "China; United States",
    notes = "Model M3 of Zhang 2022. 12 Chinese patients from CHN-FPC-21 (CTR20181119) and 39 non-Asian patients from USA-FPC-16 (NCT02739100; 13) and USA-FPC-20 (NCT02767128; 26). Demographics from Table 1 and Supplementary Table 4; dosing from Supplementary Tables 1-2."
  )

  ini({
    lcl <- log(1.02)
    label("Apparent clearance CL/F (L/h) at Fe baseline 660.8 ng/mL") # Table 2, M3 'CL/F, L/h' = 1.02 (RSE 5.9%)
    lvc <- log(3.57)
    label("Apparent volume of distribution Vd/F (L) at LBM 55.85 kg") # Table 2, M3 'Vd/F, L' = 3.57 (RSE 3.8%)

    e_iron_bl_cl <- -0.000702
    label("Linear slope of (Fe baseline - 660.8) on CL/F (per ng/mL)") # Table 2, M3 'Fe baseline on CL/F (x10^-4)' = -7.02 (RSE 23.6%); footnote c prints 1 - 0.000702
    e_lbm_vc <- 1.14
    label("Power exponent of LBM/55.85 on Vd/F (unitless)") # Table 2, M3 'LBM on Vd/F' = 1.14 (RSE 28.2%); footnote c

    # IIV: Table 2 'omega (CL, %)' and 'omega (Vd, %)' read as 100 x the SD of
    # eta (variance = (P/100)^2), consistent with the M1 rows whose printed CI
    # is on the SD-fraction scale.
    etalcl ~ 0.133956 # Table 2, M3 omega1 (CL) = 36.6% (RSE 8.9%); 0.366^2
    etalvc ~ 0.0441 # Table 2, M3 omega2 (Vd) = 21.0% (RSE 23.2%); 0.210^2

    propSd <- 0.271
    label("Proportional residual error (fraction)") # Table 2, M3 'sigma1 (prop), %' = 27.1 (RSE 7.8%)
    addSd <- fixed(0.0710)
    label("Additive residual error (ng/mL)") # Table 2, M3 'sigma2 (add), ng/ml' = 0.0710 with RSE 0.0, no CI, bootstrap median 0.0710 with no PI: held constant
  })

  model({
    # Table 2 footnote c:
    # CLind = CLTV x (1 - 0.000702 x [Fe baseline - 660.80])
    # Vind  = VTV x (LBM/55.85)^1.14
    cl <- exp(lcl + etalcl) * (1 + e_iron_bl_cl * (IRON_BL - 660.8))
    vc <- exp(lvc + etalvc) * (LBM / 55.85)^e_lbm_vc

    kel <- cl / vc
    d / dt(central) <- -kel * central

    # Dose in mg and volume in L give mg/L; x 1000 converts to ng/mL.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
