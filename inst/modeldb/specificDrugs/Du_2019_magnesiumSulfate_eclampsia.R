Du_2019_magnesiumSulfate_eclampsia <- function() {
  description <- paste0(
    "Landmark logistic-regression exposure-response model for the ",
    "probability of eclampsia in women with preeclampsia treated with ",
    "magnesium sulfate (Du 2019, n = 10280 women pooled from the Magpie ",
    "Trial and a Thai study). The exposure metric is the area under the ",
    "CHANGE-FROM-BASELINE serum magnesium concentration-time curve ",
    "(AUC, mg*h/L of elemental magnesium), which the source computes in ",
    "closed form as total dose divided by clearance rather than from an ",
    "ODE solve: AUC = (DOSE_IV + F * DOSE_IM) / CL, with CL carrying the ",
    "serum-creatinine and body-weight covariate effects of the companion ",
    "population PK model. The model therefore has NO ODE and no time ",
    "dimension - it is evaluated once per woman from her total ",
    "administered dose by route, body weight, serum creatinine and age. ",
    "The age effect is PIECE-WISE LINEAR with a knot at 22 years; the ",
    "two published branches are NOT continuous at the knot (see the ini ",
    "comments and the vignette). Doses are supplied in grams of the ",
    "administered salt MgSO4-7H2O and converted to elemental magnesium ",
    "inside the model. Blood-pressure/urinary-protein severity level and ",
    "previous anticonvulsant use were screened and not retained."
  )
  reference <- paste(
    "Du L, Wenning LA, Carvalho B, Duley L, Brookfield KF, Witjes H,",
    "de Greef R, Lumbiganon P, Titapant V, Kongwattanakul K, Long Q,",
    "Sangkomkamhang US, Gulmezoglu AMG, Oladapo OT.",
    "Alternative magnesium sulfate dosing regimens for women with",
    "preeclampsia: a population pharmacokinetic exposure-response",
    "modeling and simulation study.",
    "J Clin Pharmacol. 2019;59(11):1519-1526. doi:10.1002/jcph.1448.",
    "PMCID PMC6790709.",
    "The clearance function that generates the exposure metric is",
    "restated in this paper (Methods, Pharmacokinetic Exposure and",
    "Supplemental Table S1) from the companion population PK analysis:",
    "Du L, Wenning L, Migoya E, et al. Population pharmacokinetic",
    "modeling to evaluate standard magnesium sulfate treatments and",
    "alternative dosing regimens for women with preeclampsia.",
    "J Clin Pharmacol. 2019;59(3):374-385. doi:10.1002/jcph.1328.",
    sep = " "
  )
  vignette <- "Du_2019_magnesiumSulfate_eclampsia"

  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "g MgSO4-7H2O (TOTAL dose per route over the treatment course, supplied as the DOSE_MGSO4_IV_G and DOSE_MGSO4_IM_G covariate columns; there are no dose events)",
    concentration = "prob_eclampsia (probability of eclampsia occurrence, 0-1); the internal exposure metric auc_mg is in mg*h/L of ELEMENTAL magnesium"
  )

  covariateData <- list(
    AGE = list(
      description = "Maternal age at trial entry.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Du 2019 Table S3: Magpie median 26 (P5-P95 17-39, n = 9890), Thailand median 28 (P5-P95 17-40, n = 389). Enters UNCENTRED on a piece-wise linear scale with a knot at 22 years, the knot value selected during model building (Results, Exposure-Response Model Results: 'The final logistic E-R model included a piece-wise linear age effect with a knot point at 22 years', AIC 24 points better than a linear age effect). Du 2019 Discussion explicitly disclaims a mechanism: 'We are not aware of any clinical significance to the knot point of 22 years in relation to the risk of eclampsia.'",
      source_name = "Age"
    ),
    WT = list(
      description = "Maternal body weight at trial entry.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters only through the clearance function, normalised to 85 kg (the Stanford population PK cohort median). Du 2019 Table S3: Thailand mean 75.3 (SD 14.9), median 73.0 (P5-P95 54-103), available in 383 of 389 women; NOT collected in the Magpie Trial, where Du 2019 Methods (Estimation of AUC) assigned the Stanford median of 85 kg to every woman. The paper's own sensitivity check reports that imputing 80 kg and 0.7 mg/dL instead moved AUC by only about 5%.",
      source_name = "WT"
    ),
    CREAT = list(
      description = "Maternal serum creatinine at trial entry.",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "mg/dL, NOT umol/L: the clearance function in Du 2019 Methods is written explicitly as (0.8 mg/dL / Cr_i)^0.731. Reference value 0.8 mg/dL is the Stanford population PK cohort median. Du 2019 Table S3: Thailand mean 0.66 (SD 0.19), median 0.61 (P5-P95 0.43-1.07), available in 168 of 389 women; NOT collected in the Magpie Trial, where the Stanford median of 0.8 mg/dL was assigned to every woman. Missing Thai values were imputed at the Thai median of 0.6 mg/dL (n = 221).",
      source_name = "Cr"
    ),
    DOSE_MGSO4_IV_G = list(
      description = "TOTAL intravenous magnesium sulfate dose administered over the treatment course, expressed as grams of the administered salt MgSO4-7H2O (loading plus maintenance combined). Zero for intramuscular-only and placebo regimens.",
      units = "g MgSO4-7H2O",
      type = "continuous",
      reference_category = NULL,
      notes = "Du 2019 Methods (Pharmacokinetic Exposure) DOSE_IV,i. A TOTAL over the course, not a per-administration amount: 'the total MgSO4 dose administered up to the day of first eclampsia was used in the AUC calculation'. Table S4 enumerates the regimens; e.g. the standard Zuspan regimen (4 g over 20 min then 1 g/h for 24 h) is 28 g. Converted to elemental magnesium inside model() using the Table 1 footnote 'Dosage form in simulated dosing regimen was MgSO4-7H2O, which contains ~10% of magnesium'.",
      source_name = "DOSE_IV"
    ),
    DOSE_MGSO4_IM_G = list(
      description = "TOTAL intramuscular magnesium sulfate dose administered over the treatment course, expressed as grams of the administered salt MgSO4-7H2O (loading plus maintenance combined). Zero for intravenous-only and placebo regimens.",
      units = "g MgSO4-7H2O",
      type = "continuous",
      reference_category = NULL,
      notes = "Du 2019 Methods (Pharmacokinetic Exposure) DOSE_IM,i. Scaled by the intramuscular absolute bioavailability before entering the AUC, which is why it is a separate column from DOSE_MGSO4_IV_G. Table S4 enumerates the regimens; e.g. the standard Pritchard regimen (4 g IV plus 10 g IM loading, then 5 g IM every 4 h for 5 doses) is 4 g IV and 35 g IM, total 39 g, matching the Table 1 total-dose column.",
      source_name = "DOSE_IM"
    )
  )

  covariatesDataExcluded <- list(
    PRIOR_ANTICONV = list(
      description = "Previous anticonvulsant use before magnesium sulfate treatment (yes/no).",
      units = "(binary)",
      type = "binary",
      notes = "Screened as a candidate covariate on the eclampsia logit and NOT retained: Du 2019 Results, Exposure-Response Model Results -- 'Level of blood pressure/urinary protein and anticonvulsant drug use before MgSO4 treatment did not show a significant relationship with eclampsia'. No point estimate exists on disk. Table S3 reports the Magpie distribution as No 8969 (90.7%), Yes 861 (8.7%), Unknown 61 (0.6%); not collected in the Thai study."
    ),
    BP_PROTEINURIA_LEVEL = list(
      description = "Composite preeclampsia severity category built from blood pressure and urinary protein: level 2 = diastolic BP at least 110 mmHg on two occasions or systolic BP at least 170 mmHg on two occasions plus at least 3+ proteinuria; level 1 = diastolic BP at least 100 mmHg on two occasions or systolic BP at least 150 mmHg on two occasions plus at least 2+ proteinuria; level 0 = neither.",
      units = "(0, 1, 2)",
      type = "categorical",
      notes = "Screened and NOT retained (same Results sentence as PRIOR_ANTICONV). No point estimate exists on disk. Table S3 Magpie distribution: level 2 2459 (24.9%), level 1 4252 (43.0%), level 0 3176 (32.1%), unknown 4 (0.04%); not collected in the Thai study. Du 2019 notes the three levels were derived for this analysis and are only similar to, not identical with, the Magpie Trial's own dichotomous severe/not-severe classification, because signs and symptoms of imminent eclampsia were unavailable."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 10280L,
    n_studies = 2L,
    n_observations = "10280 binary eclampsia records (one per woman; landmark analysis, no repeated measures). 127 women developed eclampsia: 37 on magnesium sulfate (36 Magpie, 1 Thai) and 90 on placebo (all Magpie).",
    age_range = "Magpie median 26 years (P5-P95 17-39); Thailand median 28 years (P5-P95 17-40)",
    weight_range = "Thailand median 73 kg (P5-P95 54-103); not collected in Magpie, imputed at the Stanford median of 85 kg",
    disease_state = "preeclampsia; Magpie women had not given birth or were within 24 h postpartum (median gestational age 37 weeks, P5-P95 27-41), Thai women had median gestational age 38 weeks (P5-P95 31-40)",
    dose_range = "Magpie: 4 g IV loading then either 1 g/h IV for 24 h, or 10 g IM loading plus 5 g IM every 4 h for 24 h, or matching placebo. Thailand: 4 g IV loading then 1 g/h IV infusion. Simulated regimens span 6-54 g total per 24 h (Table S4).",
    regions = "Magpie Trial: multicentre randomised trial across 33 countries (1998-2001). Thai study: Srinagarind Hospital (Khon Kaen University) and Siriraj Hospital (Mahidol University), Thailand (2010-2013).",
    notes = paste0(
      "This is the exposure-response analysis set (Du 2019 Table S2). ",
      "The clearance function that generates the exposure metric was ",
      "estimated on a DIFFERENT and much smaller cohort: 92 women with ",
      "preeclampsia in the Stanford study (2012-2014, 19-44 years, ",
      "median gestational age 36 weeks) with rich serial magnesium ",
      "sampling, reported in the companion paper doi:10.1002/jcph.1328 ",
      "and restated in Table S1 here. No serum magnesium samples were ",
      "collected in the Magpie Trial and only sparse post-dose samples ",
      "in the Thai study, so no individual magnesium concentration in ",
      "the exposure-response analysis set informs its own AUC."
    )
  )

  ini({
    # ==================================================================
    # EXPOSURE LAYER
    #
    # Du 2019 Methods, Pharmacokinetic Exposure, prints the clearance
    # function and the AUC definition verbatim:
    #
    #   CL_i = 3.72 * (0.8 mg/dL / Cr_i)^0.731 * (WT_i / 85 kg)^0.75
    #   AUC_i = DOSE_IV,i / CL_i   and   AUC_i = F * DOSE_IM,i / CL_i
    #
    # Supplemental Table S1 gives the same function in the equivalent
    # form (Cr_i / 0.8)^theta with theta = -0.731 (RSE 14.2%), which is
    # how the exponent is signed below. These values are FIXED here:
    # they are the companion population PK model's estimates
    # (doi:10.1002/jcph.1328, restated in Table S1), carried into this
    # analysis as constants rather than re-estimated.
    # ==================================================================
    lcl <- fixed(log(3.72))
    label("Typical magnesium clearance at 85 kg and 0.8 mg/dL serum creatinine (L/h)")
    e_creat_cl <- fixed(-0.731)
    label("Power exponent on (CREAT / 0.8 mg/dL) for clearance (unitless)")
    e_wt_cl <- fixed(0.75)
    label("Power exponent on (WT / 85 kg) for clearance (unitless)")
    lfdepot <- fixed(log(0.862))
    label("Absolute bioavailability of intramuscular magnesium sulfate (fraction)")

    # ==================================================================
    # EXPOSURE-RESPONSE LAYER
    #
    # Du 2019 Results, Exposure-Response Model Results, prints the final
    # model as two branches (the piece-wise linear age effect, knot at
    # 22 years):
    #
    #   Age <= 22 years: Logit(P) = -6.29 - 0.00164 * AUC + 0.154 * age
    #   Age >  22 years: Logit(P) = -4.72 - 0.00164 * AUC + 0.010 * age
    #
    # NOTE ON DISCONTINUITY. The Methods equation is written in the
    # continuous knot form b0 + b1*AUC + b2*Age + delta*(Age-k)*I(Age>k),
    # which is continuous at the knot by construction, but the two
    # printed branches are NOT: at age 22 and AUC 0 the lower branch
    # gives a logit of -2.902 and the upper branch -4.500, a jump of
    # 1.598 logit units (about a 5-fold drop in predicted risk). The
    # two branches as printed are nevertheless what generated the
    # paper's own Table 1 -- they reproduce all 18 regimen rows of both
    # predicted-eclampsia-rate columns to the printed precision (see the
    # vignette) -- so they are encoded exactly as printed rather than
    # reconciled to the continuous form. No erratum exists.
    #
    # Standard errors are not reported for any of these coefficients.
    # ==================================================================
    logit_ref_agele22 <- -6.29
    label("Intercept of the eclampsia logit for women aged 22 years or younger (unitless logit)")
    logit_ref_agegt22 <- -4.72
    label("Intercept of the eclampsia logit for women older than 22 years (unitless logit)")
    e_auc_eclampsia <- -0.00164
    label("Change in the eclampsia logit per 1 mg*h/L of change-from-baseline magnesium AUC (unitless logit per mg*h/L)")
    e_age_eclampsia_agele22 <- 0.154
    label("Change in the eclampsia logit per year of age for women aged 22 years or younger (unitless logit per year)")
    e_age_eclampsia_agegt22 <- 0.010
    label("Change in the eclampsia logit per year of age for women older than 22 years (unitless logit per year)")
    age_knot_eclampsia <- fixed(22)
    label("Knot point of the piece-wise linear age effect on the eclampsia logit (years)")

    # ==================================================================
    # The source fits a logistic regression with an exact Bernoulli
    # likelihood: no between-subject random effect and no residual error
    # are estimated. rxode2 requires an observation declaration, so the
    # deterministic probability carries a tiny placeholder additive
    # residual, following Babel_2026_telisotuzumab_orr.R and the
    # Fukae_2024_valemetostat_* family.
    # ==================================================================
    addSd_prob_eclampsia <- fixed(0.001)
    label("Placeholder additive residual SD on the typical-value eclampsia probability; the source likelihood is Bernoulli (no source residual)")
  })

  model({
    # ---- Exposure ---------------------------------------------------
    # Individual magnesium clearance. Written as (CREAT / 0.8)^-0.731 to
    # match Supplemental Table S1's signed exponent; algebraically
    # identical to the Methods form (0.8 / CREAT)^0.731.
    cl <- exp(lcl) * (CREAT / 0.8)^e_creat_cl * (WT / 85)^e_wt_cl

    # Elemental-magnesium content of the administered salt. Table 1
    # footnote: 'Dosage form in simulated dosing regimen was MgSO4-7H2O,
    # which contains ~10% of magnesium.' The exact mass fraction is
    # 24.305 / 246.47 = 0.098612 g Mg per g MgSO4-7H2O (atomic mass of
    # magnesium over the formula mass of the heptahydrate); the 1000
    # converts the covariate columns from grams of salt to milligrams of
    # elemental magnesium so that auc_mg lands in the paper's mg*h/L.
    doseMg <- (DOSE_MGSO4_IV_G + exp(lfdepot) * DOSE_MGSO4_IM_G) * 1000 * 0.098612

    # Total (0 to infinity) AUC of the CHANGE FROM BASELINE in serum
    # magnesium. Du 2019 computes this in closed form as dose over
    # clearance, not by integrating a solved profile; for a linear
    # disposition the two are identical, which is why no volume,
    # intercompartmental clearance or absorption rate constant appears
    # in this model.
    auc_mg <- doseMg / cl

    # ---- Exposure-response ------------------------------------------
    # Piece-wise linear in age about the 22-year knot; see the ini
    # comments for why the two branches are encoded as printed.
    if (AGE <= age_knot_eclampsia) {
      logit_eclampsia <- logit_ref_agele22 + e_age_eclampsia_agele22 * AGE + e_auc_eclampsia * auc_mg
    } else {
      logit_eclampsia <- logit_ref_agegt22 + e_age_eclampsia_agegt22 * AGE + e_auc_eclampsia * auc_mg
    }
    prob_eclampsia <- expit(logit_eclampsia)

    # Deterministic probability of eclampsia occurrence. Downstream
    # callers can sample binary outcomes with
    # rbinom(n, 1, prob_eclampsia) on the rxSolve output.
    prob_eclampsia ~ add(addSd_prob_eclampsia)
  })
}
