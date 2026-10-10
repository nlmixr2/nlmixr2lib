Hosey_2022_t2dm_growth <- function() {
  description <- paste(
    "PBPK system-parameter (anthropometry) model: sex-specific height and",
    "body-weight versus age growth curves for youth who develop type 2",
    "diabetes mellitus (Hosey 2022 Clin Transl Sci). This is a",
    "physiology layer for pediatric PBPK scaling, not a drug model: it",
    "has no drug, no dosing, no compartments and no ODEs, and both",
    "outputs are algebraic functions of the rxode2 time variable, which",
    "the model interprets as chronological age in YEARS (intended range",
    "2 to 18 years). Height is the CDC three-logistic stature function",
    "refitted to electronic-medical-record data from 356 children with",
    "type 2 diabetes; weight is the CDC tenth-degree (male) / ninth-degree",
    "(female) polynomial refitted to the same cohort, using the",
    "six-significant-figure coefficients of Supplemental Information S1.",
    "Residual variability is proportional with the sex- and",
    "output-specific coefficients of variation of the validation set",
    "(Table 2), the same bands drawn in Figure 2. SEXF is required.",
    sep = " "
  )
  reference <- paste(
    "Hosey CM, Halpin K, Shakhnovich V, Bi C, Sweeney B, Yan Y, Leeder JS.",
    "Pediatric growth patterns in youth-onset type 2 diabetes mellitus:",
    "Implications for physiologically-based pharmacokinetic models.",
    "Clin Transl Sci. 2022;15(4):912-922.",
    "doi:10.1111/cts.13207.",
    sep = " "
  )
  vignette <- "Hosey_2022_t2dm_growth"
  units <- list(
    time = "year (chronological age; intended range 2 to 18)",
    dosing = "n/a (no exogenous dosing; PBPK system-physiology model)",
    concentration = "n/a (no drug concentration). Outputs are height (cm) and body weight (kg)."
  )

  covariateData <- list(
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Time-fixed and REQUIRED: Hosey 2022 fits a separate height",
        "equation and a separate weight polynomial for each sex (Results,",
        "'Truncated versions of the final equations'; Supplemental",
        "Information S1), and the residual coefficients of variation are",
        "sex-specific (Table 2). The model selects the female branch when",
        "SEXF is 1 and the male branch when SEXF is 0; any other value",
        "linearly blends the two sexes and is not meaningful."
      ),
      source_name = "sex"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 356L,
    n_studies = 1L,
    age_range = "2 to 18 years (target population); data from 19 to 25 years were also used in model development to stabilise the curves at older ages (Methods, Dataset)",
    sex_female_pct = 59.3,
    race_ethnicity = "Self-identified: female 41.7% White, 32.2% Black, 15.7% Hispanic, 2.8% Asian, 5.2% multiple, 1.9% other/unknown; male 53.1% White, 20.7% Black, 14.5% Hispanic, 4.1% Asian, 3.5% multiple, 4.1% other/unknown (Supplemental Information S4)",
    disease_state = paste(
      "Youth with a diagnosis of type 2 diabetes mellitus (ICD-9 / ICD-10",
      "codes), most with overweight or obesity; children with type 1",
      "diabetes were excluded. Visits both before and after the type 2",
      "diabetes diagnosis were included (Methods, Dataset)."
    ),
    dose_range = "n/a (no drug)",
    regions = "United States (Children's Mercy Kansas City, Kansas and Missouri)",
    n_female = 211L,
    n_male = 145L,
    n_height_points = "female 3973 training + 458 validation; male 2610 training + 306 validation (Table 1)",
    n_weight_points = "female 5189 training + 576 validation; male 3350 training + 380 validation (Table 1)",
    notes = paste(
      "Deidentified electronic-medical-record height, weight, age and sex",
      "from visits recorded 2010-01-01 to 2016-08-24 (with some older",
      "records converted back to 1997). Data were binned in 6-month",
      "intervals and pruned with the modified Thompson-Tau algorithm;",
      "90% of the data were used for development (100 x 5-fold cross",
      "validation with an iteratively reweighted least squares robust",
      "M-estimator, robustbase::nlrob) and 10% were held out for",
      "validation. Weight coefficients are the average over the 500",
      "cross-validation fits; the height curve was refitted through the",
      "averaged monthly predictions (Methods, Model development)."
    )
  )

  ini({
    # Height: Ht(age) = a*q/(1+exp(-b1*(age-c1))) + a*p/(1+exp(-b2*(age-c2)))
    #                   + (f-a)/(1+exp(-b3*(age-c3)))
    # (Methods, CDC stature function; Results, 'Truncated versions of the
    # final equations'). The paper reports these to three significant
    # figures only; no higher-precision height coefficients were published.
    ht_a_male <- fixed(166)
    label("Height curve coefficient a, male (cm)") # Results, 'Male height (cm, age in years)' equation: 166 multiplies 0.283 and 0.717 and is subtracted from 178
    ht_q_male <- fixed(0.283)
    label("Height curve fraction q of a in the first logistic, male (fraction)") # Results, 'Male height' equation: 166 * 0.283
    ht_p_male <- fixed(0.717)
    label("Height curve fraction p of a in the second logistic, male (fraction)") # Results, 'Male height' equation: 166 * 0.717
    ht_b1_male <- fixed(5.18)
    label("Height curve first-logistic rate b1, male (1/year)") # Results, 'Male height' equation: exp(-5.18(age-1.44))
    ht_c1_male <- fixed(1.44)
    label("Height curve first-logistic midpoint c1, male (year)") # Results, 'Male height' equation: exp(-5.18(age-1.44))
    ht_b2_male <- fixed(0.245)
    label("Height curve second-logistic rate b2, male (1/year)") # Results, 'Male height' equation: exp(-0.245(age-5.11))
    ht_c2_male <- fixed(5.11)
    label("Height curve second-logistic midpoint c2, male (year)") # Results, 'Male height' equation: exp(-0.245(age-5.11))
    ht_f_male <- fixed(178)
    label("Height curve coefficient f, male (cm)") # Results, 'Male height' equation: numerator (178 - 166)
    ht_b3_male <- fixed(1.8)
    label("Height curve third-logistic (pubertal) rate b3, male (1/year)") # Results, 'Male height' equation: exp(-1.8(age-12.4))
    ht_c3_male <- fixed(12.4)
    label("Height curve third-logistic (pubertal) midpoint c3, male (year)") # Results, 'Male height' equation: exp(-1.8(age-12.4))

    ht_a_female <- fixed(160)
    label("Height curve coefficient a, female (cm)") # Results, 'Female height (cm, age in years)' equation: 160 multiplies 0.396 and 0.604 and is subtracted from 209
    ht_q_female <- fixed(0.396)
    label("Height curve fraction q of a in the first logistic, female (fraction)") # Results, 'Female height' equation: 160 * 0.396
    ht_p_female <- fixed(0.604)
    label("Height curve fraction p of a in the second logistic, female (fraction)") # Results, 'Female height' equation: 160 * 0.604
    ht_b1_female <- fixed(0.537)
    label("Height curve first-logistic rate b1, female (1/year)") # Results, 'Female height' equation: exp(-0.537(age-2.23))
    ht_c1_female <- fixed(2.23)
    label("Height curve first-logistic midpoint c1, female (year)") # Results, 'Female height' equation: exp(-0.537(age-2.23))
    ht_b2_female <- fixed(0.0039)
    label("Height curve second-logistic rate b2, female (1/year)") # Results, 'Female height' equation: exp(-0.0039(age-1.77)); small as printed, and reproduces Figure 5 (see vignette)
    ht_c2_female <- fixed(1.77)
    label("Height curve second-logistic midpoint c2, female (year)") # Results, 'Female height' equation: exp(-0.0039(age-1.77))
    ht_f_female <- fixed(209)
    label("Height curve coefficient f, female (cm)") # Results, 'Female height' equation: numerator (209 - 160)
    ht_b3_female <- fixed(0.473)
    label("Height curve third-logistic (pubertal) rate b3, female (1/year)") # Results, 'Female height' equation: exp(-0.473(age-9.54))
    ht_c3_female <- fixed(9.54)
    label("Height curve third-logistic (pubertal) midpoint c3, female (year)") # Results, 'Female height' equation: exp(-0.473(age-9.54))

    # Weight: Wt(age) = a*age^10 (males only) + b*age^9 + ... + j*age + k
    # (Methods). Coefficients are the six-significant-figure values of
    # Supplemental Information S1 (Results: 'at least six significant
    # figures appear to be needed'); the three-figure main-text values are
    # not usable (they give about 420 kg for a 10-year-old boy).
    bw_a_male <- fixed(3.44781e-08)
    label("Weight polynomial coefficient a (age^10), male (kg/year^10)") # Supplemental Information S1, row 'a', column 'Male' = 3.44781E-08 (main text 3.45e-8)
    bw_b_male <- fixed(-4.57950e-06)
    label("Weight polynomial coefficient b (age^9), male (kg/year^9)") # Supplemental Information S1, row 'b', column 'Male' = -4.57950E-06 (main text -4.58e-6)
    bw_c_male <- fixed(2.60885e-04)
    label("Weight polynomial coefficient c (age^8), male (kg/year^8)") # Supplemental Information S1, row 'c', column 'Male' = 2.60885E-04 (main text 2.61e-4)
    bw_d_male <- fixed(-8.35633e-03)
    label("Weight polynomial coefficient d (age^7), male (kg/year^7)") # Supplemental Information S1, row 'd', column 'Male' = -8.35633E-03 (main text -8.36e-3)
    bw_e_male <- fixed(0.165867)
    label("Weight polynomial coefficient e (age^6), male (kg/year^6)") # Supplemental Information S1, row 'e', column 'Male' = 1.65867E-01 (main text 0.166)
    bw_f_male <- fixed(-2.12088)
    label("Weight polynomial coefficient f (age^5), male (kg/year^5)") # Supplemental Information S1, row 'f', column 'Male' = -2.12088E+00 (main text -2.12)
    bw_g_male <- fixed(17.5837)
    label("Weight polynomial coefficient g (age^4), male (kg/year^4)") # Supplemental Information S1, row 'g', column 'Male' = 1.75837E+01 (main text 17.6)
    bw_h_male <- fixed(-92.5981)
    label("Weight polynomial coefficient h (age^3), male (kg/year^3)") # Supplemental Information S1, row 'h', column 'Male' = -9.25981E+01 (main text -92.6)
    bw_i_male <- fixed(293.823)
    label("Weight polynomial coefficient i (age^2), male (kg/year^2)") # Supplemental Information S1, row 'i', column 'Male' = 293.823 (main text 294)
    bw_j_male <- fixed(-499.980)
    label("Weight polynomial coefficient j (age), male (kg/year)") # Supplemental Information S1, row 'j', column 'Male' = -499.980 (main text -500)
    bw_k_male <- fixed(358.440)
    label("Weight polynomial intercept k, male (kg)") # Supplemental Information S1, row 'k', column 'Male' = 358.440 (main text 358)

    bw_b_female <- fixed(-4.81428e-07)
    label("Weight polynomial coefficient b (age^9), female (kg/year^9)") # Supplemental Information S1, row 'b', column 'Female' = -4.81428E-07 (main text -4.81e-7); the female polynomial has no age^10 term (row 'a' blank)
    bw_c_female <- fixed(5.00192e-05)
    label("Weight polynomial coefficient c (age^8), female (kg/year^8)") # Supplemental Information S1, row 'c', column 'Female' = 5.00192E-05 (main text 5.00e-5)
    bw_d_female <- fixed(-2.19327e-03)
    label("Weight polynomial coefficient d (age^7), female (kg/year^7)") # Supplemental Information S1, row 'd', column 'Female' = -2.19327E-03 (main text -2.19e-3)
    bw_e_female <- fixed(0.0529444)
    label("Weight polynomial coefficient e (age^6), female (kg/year^6)") # Supplemental Information S1, row 'e', column 'Female' = 5.29444E-02 (main text 0.0529)
    bw_f_female <- fixed(-0.770305)
    label("Weight polynomial coefficient f (age^5), female (kg/year^5)") # Supplemental Information S1, row 'f', column 'Female' = -7.70305E-01 (main text -0.770)
    bw_g_female <- fixed(6.95625)
    label("Weight polynomial coefficient g (age^4), female (kg/year^4)") # Supplemental Information S1, row 'g', column 'Female' = 6.95625E+00 (main text 6.96)
    bw_h_female <- fixed(-38.7195)
    label("Weight polynomial coefficient h (age^3), female (kg/year^3)") # Supplemental Information S1, row 'h', column 'Female' = -3.87195E+01 (main text -38.7)
    bw_i_female <- fixed(127.524)
    label("Weight polynomial coefficient i (age^2), female (kg/year^2)") # Supplemental Information S1, row 'i', column 'Female' = 1.27524E+02 (main text 128)
    bw_j_female <- fixed(-221.140)
    label("Weight polynomial coefficient j (age), female (kg/year)") # Supplemental Information S1, row 'j', column 'Female' = -2.21140E+02 (main text -221)
    bw_k_female <- fixed(164.303)
    label("Weight polynomial intercept k, female (kg)") # Supplemental Information S1, row 'k', column 'Female' = 1.64303E+02 (main text 164)

    # Residual variability: coefficient of variation of the validation set
    # fitted to the T2DM model (Table 2, column 'Validation set fit to
    # growth model: T2DM'). Figure 2's SD bands are the model curve times
    # (1 +/- CV), i.e. a proportional error.
    propSd_ht_male <- fixed(0.0587)
    label("Proportional residual SD of height, male (fraction)") # Table 2, 'Coefficients of variation (%)', Height, Male, T2DM column = 5.87
    propSd_ht_female <- fixed(0.0675)
    label("Proportional residual SD of height, female (fraction)") # Table 2, 'Coefficients of variation (%)', Height, Female, T2DM column = 6.75
    propSd_bw_male <- fixed(0.378)
    label("Proportional residual SD of body weight, male (fraction)") # Table 2, 'Coefficients of variation (%)', Weight, Male, T2DM column = 37.8
    propSd_bw_female <- fixed(0.459)
    label("Proportional residual SD of body weight, female (fraction)") # Table 2, 'Coefficients of variation (%)', Weight, Female, T2DM column = 45.9
  })

  model({
    # rxode2 time is chronological age in years
    age <- time

    # CDC stature function, refitted per sex (Methods; Results)
    ht_male <- ht_a_male * ht_q_male / (1 + exp(-ht_b1_male * (age - ht_c1_male))) +
      ht_a_male * ht_p_male / (1 + exp(-ht_b2_male * (age - ht_c2_male))) +
      (ht_f_male - ht_a_male) / (1 + exp(-ht_b3_male * (age - ht_c3_male)))
    ht_female <- ht_a_female * ht_q_female / (1 + exp(-ht_b1_female * (age - ht_c1_female))) +
      ht_a_female * ht_p_female / (1 + exp(-ht_b2_female * (age - ht_c2_female))) +
      (ht_f_female - ht_a_female) / (1 + exp(-ht_b3_female * (age - ht_c3_female)))

    # CDC weight polynomial, refitted per sex (Methods; Supplemental Information S1)
    bw_male <- bw_a_male * age^10 + bw_b_male * age^9 + bw_c_male * age^8 +
      bw_d_male * age^7 + bw_e_male * age^6 + bw_f_male * age^5 +
      bw_g_male * age^4 + bw_h_male * age^3 + bw_i_male * age^2 +
      bw_j_male * age + bw_k_male
    bw_female <- bw_b_female * age^9 + bw_c_female * age^8 +
      bw_d_female * age^7 + bw_e_female * age^6 + bw_f_female * age^5 +
      bw_g_female * age^4 + bw_h_female * age^3 + bw_i_female * age^2 +
      bw_j_female * age + bw_k_female

    ht <- SEXF * ht_female + (1 - SEXF) * ht_male
    bw <- SEXF * bw_female + (1 - SEXF) * bw_male

    # Sex-specific proportional residual SD (Table 2)
    propSd_ht <- SEXF * propSd_ht_female + (1 - SEXF) * propSd_ht_male
    propSd_bw <- SEXF * propSd_bw_female + (1 - SEXF) * propSd_bw_male

    ht ~ prop(propSd_ht)
    bw ~ prop(propSd_bw)
  })
}
