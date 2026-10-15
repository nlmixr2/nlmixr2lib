Zhang_2022_ferricPyrophosphateCitrate_dialysate <- function() {
  description <- "One-compartment population PK model for ferric pyrophosphate citrate (FPC, Triferic) delivered via the dialysate over a 4-hour haemodialysis session to Asian and non-Asian adults with haemodialysis-dependent stage 5 chronic kidney disease (Zhang 2022, model M2). Output Cc is the drug-derived (baseline-corrected) serum total iron concentration. Apparent clearance CL/F decreases linearly with pre-dose serum total iron; apparent volume Vd/F carries a power effect of lean body mass."
  reference <- "Zhang L, Gan L, Li K, Xie P, Tan Y, Wei G, Yuan X, Pratt R, Zhou Y, Hui AM, Fang Y, Zuo L, Zheng Q. Ethnicity evaluation of ferric pyrophosphate citrate among Asian and Non-Asian populations: a population pharmacokinetics analysis. Eur J Clin Pharmacol. 2022;78(9):1421-1434. doi:10.1007/s00228-022-03328-9"
  vignette <- "Zhang_2022_ferricPyrophosphateCitrate"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    LBM = list(
      description = "Lean body mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline value. Power effect on Vd/F normalized to the M2 median 56.28 kg (Table 2 footnote b). The paper does not name the LBM formula; M2 cohort mean (SD) by study 50.6 (6.6), 62.3 (8.6) and 57.3 (5.9) kg, overall range 37.7-77.0 kg (Supplementary Table 4).",
      source_name = "LBM"
    ),
    IRON_BL = list(
      description = "Baseline serum total iron concentration, here the single value measured at 0 h before administration",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed per administration. The paper's Fe baseline: serum total iron at 0 h before the FPC dose. Serum iron in CKD-5HD patients does not fluctuate materially over 24 h, so the authors used the pre-dose value directly. Linear effect on CL/F centered at the M2 median 642.00 ng/mL (Table 2 footnote b). M2 cohort mean (SD) by study 798.9 (226.9), 563.9 (194.2) and 613.2 (246.2) ng/mL, range 210-1204 ng/mL (Supplementary Table 4).",
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
    n_subjects = 50L,
    n_studies = 3L,
    age_range = "25-77 years",
    age_median = "not reported (study means 54.3, 49.2 and 53.9 years)",
    weight_range = "51.3-129 kg",
    weight_median = "not reported (study means 69.4, 97.8 and 84.9 kg)",
    sex_female_pct = 20,
    race_ethnicity = c(Asian = 24, `Non-Asian` = 76),
    disease_state = "Haemodialysis-dependent stage 5 chronic kidney disease (CKD-5HD)",
    dose_range = "FPC added to the bicarbonate dialysate at 95 ug/L iron (CHN-FPC-21) or 2 uM (about 110 ug/L) iron (USA-FPC-16, USA-FPC-20) for a 4-hour dialysis session; modelled as a 6.5 mg zero-order input over 4 h (Figure 1)",
    regions = "China; United States",
    notes = "Model M2 of Zhang 2022. 12 Chinese patients from CHN-FPC-21 (CTR20181119) and 38 non-Asian patients from USA-FPC-16 (NCT02739100; 13) and USA-FPC-20 (NCT02767128; 25; one USA-FPC-20 patient lacked dialysate PK data). Demographics from Table 1 and Supplementary Table 4; dosing from Supplementary Tables 1-2."
  )

  ini({
    lcl <- log(0.982)
    label("Apparent clearance CL/F (L/h) at Fe baseline 642 ng/mL") # Table 2, M2 'CL/F, L/h' = 0.982 (RSE 6.3%)
    lvc <- log(3.32)
    label("Apparent volume of distribution Vd/F (L) at LBM 56.28 kg") # Table 2, M2 'Vd/F, L' = 3.32 (RSE 3.5%)

    e_iron_bl_cl <- -0.000728
    label("Linear slope of (Fe baseline - 642) on CL/F (per ng/mL)") # Table 2, M2 'Fe baseline on CL/F (x10^-4)' = -7.28 (RSE 6.4%); footnote b prints 1 - 0.000728
    e_lbm_vc <- 0.726
    label("Power exponent of LBM/56.28 on Vd/F (unitless)") # Table 2, M2 'LBM on Vd/F' = 0.726 (RSE 35.4%); footnote b

    # IIV: Table 2 'omega (CL, %)' and 'omega (Vd, %)' read as 100 x the SD of
    # eta (variance = (P/100)^2), consistent with the M1 rows whose printed CI
    # is on the SD-fraction scale.
    etalcl ~ 0.172225 # Table 2, M2 omega1 (CL) = 41.5% (RSE 11.9%); 0.415^2
    etalvc ~ 0.032761 # Table 2, M2 omega2 (Vd) = 18.1% (RSE 21.0%); 0.181^2

    propSd <- 0.236
    label("Proportional residual error (fraction)") # Table 2, M2 'sigma1 (prop), %' = 23.6 (RSE 8.9%)
    addSd <- fixed(0.0220)
    label("Additive residual error (ng/mL)") # Table 2, M2 'sigma2 (add), ng/ml' = 0.0220 with RSE 0.0, no CI, bootstrap median 0.0220 with no PI: held constant
  })

  model({
    # Table 2 footnote b:
    # CLind = CLTV x (1 - 0.000728 x [Fe baseline - 642.00])
    # Vind  = VTV x (LBM/56.28)^0.726
    cl <- exp(lcl + etalcl) * (1 + e_iron_bl_cl * (IRON_BL - 642))
    vc <- exp(lvc + etalvc) * (LBM / 56.28)^e_lbm_vc

    kel <- cl / vc
    d / dt(central) <- -kel * central

    # Dose in mg and volume in L give mg/L; x 1000 converts to ng/mL.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
