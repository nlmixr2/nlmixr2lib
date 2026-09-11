Johnson_2025_osimertinib <- function() {
  description <- "Joint one-compartment population PK model for osimertinib and its active metabolite AZ5104 in patients with advanced EGFR-mutation-positive non-small cell lung cancer (NSCLC) (Johnson 2025). First-order oral absorption of osimertinib into a parent central compartment is followed by a metabolite compartment in series; the fraction of parent converted to AZ5104 is fixed at 0.25 per the publication. This is the AURA3 + FLAURA update of the earlier AURA / AURA2 analysis (see modellib('Brown_2017_osimertinib')), externally validated against the adjuvant ADAURA study. Baseline body weight (power form on parent CL/F and Vc/F and on AZ5104 CL/F), baseline serum albumin (power form on parent CL/F and Vc/F and on AZ5104 CL/F and Vc/F), and grouped race (Chinese, Japanese, Asian-other, and non-Asian non-White exponential factors on AZ5104 CL/F) were retained in the final covariate model. Random effects on the two apparent clearances are correlated."
  reference <- paste(
    "Johnson M, Lin YW, Schmidt H, Sunnaker M, Van Maanen E, Huang X,",
    "Rukazenkov Y, Tomkinson H, Vishwanathan K. Population",
    "pharmacokinetics of osimertinib in patients with non-small cell",
    "lung cancer. Pharmacol Res Perspect. 2025;13(3):e70098.",
    "doi:10.1002/prp2.70098. PMC12035414.",
    "Covariate model equations 2.1 and 2.2 are from the Supplementary",
    "Appendix (Supplementary Methods, 'Modelling Techniques').",
    sep = " "
  )
  vignette <- "Johnson_2025_osimertinib"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Johnson 2025 Methods section 2.2 states that plasma
  # samples were assayed for both osimertinib and AZ5104, so the analyte /
  # specimen assignment of both central compartments is verified against
  # the source.
  compartmentData <- list(
    depot          = list(analyte = "osimertinib", units = "mg", specimen = "administration site", verified = TRUE),
    central        = list(analyte = "osimertinib", units = "mg", specimen = "plasma", verified = TRUE),
    central_az5104 = list(analyte = "AZ5104", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline total body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form (median-normalized) effects on parent CL/F (exponent 0.421), parent Vc/F (exponent 0.814), and AZ5104 CL/F (exponent 0.822). Reference weight 61 kg, stated in Johnson 2025 Results section 3.3 ('the median body weight (61 kg)') and in the Figure 2 caption ('Typical patient: White, 61 kg body weight, 39 g/L baseline serum albumin'). 61 kg is the median of the three model-building cohorts (Table 1 medians 60.5 kg for AURA + AURA2, 60.0 kg for AURA3, 62.0 kg for FLAURA); the Table 1 'Overall' median of 62.0 kg additionally includes the ADAURA external-validation set, which was excluded from model building. Baseline (not time-varying): Johnson 2025 Table 1 reports weight as a baseline characteristic.",
      source_name        = "Weight"
    ),
    ALB = list(
      description        = "Baseline serum albumin concentration.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form (median-normalized) effects on parent CL/F (exponent 0.825), parent Vc/F (exponent 2.27), AZ5104 CL/F (exponent 0.928), and AZ5104 Vc/F (exponent -0.831). Reference albumin 39 g/L, stated in Johnson 2025 Results section 3.3 ('the median albumin levels (39 g/L)') and in the Figure 2 caption. Values are already in SI g/L in the source (Table 1 'Baseline albumin (g/L)'), so no g/dL conversion is applied. 39 g/L is the AURA + AURA2 median; AURA3 and FLAURA medians are 40.0 g/L, and the Table 1 'Overall' median of 40.0 g/L includes the excluded ADAURA set.",
      source_name        = "Baseline albumin"
    ),
    RACE_CHINESE = list(
      description        = "Chinese-heritage grouped-race indicator (1 = Chinese, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Chinese; the paper-defined reference category of the grouped-race covariate is White).",
      notes              = "Exponential effect exp(0.076 * RACE_CHINESE) on AZ5104 CL/F, i.e. ~7.9 percent higher AZ5104 apparent clearance and ~7.3 percent lower AZ5104 AUCss than White patients. Johnson 2025 'grouped race' (Table 1 footnote a) splits the Asian population into Japanese, Chinese, and non-Japanese/non-Chinese Asian, and carries White (reference), Chinese, Japanese, Asian-other, and Other as five mutually exclusive categories; a White patient is encoded as 0 for every RACE_* indicator. Race had no effect on osimertinib (parent) apparent clearance (Results section 3.3).",
      source_name        = "Grouped race: Chinese (Johnson 2025 Table 2)"
    ),
    RACE_JAPANESE = list(
      description        = "Japanese-heritage grouped-race indicator (1 = Japanese, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Japanese; the paper-defined reference category of the grouped-race covariate is White).",
      notes              = "Exponential effect exp(0.184 * RACE_JAPANESE) on AZ5104 CL/F, i.e. ~20.2 percent higher AZ5104 apparent clearance and ~16.8 percent lower AZ5104 AUCss than White patients.",
      source_name        = "Grouped race: Japanese (Johnson 2025 Table 2)"
    ),
    RACE_ASIAN_OTH = list(
      description        = "Asian-other grouped-race indicator (1 = Asian heritage other than Chinese or Japanese, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not Asian-other; the paper-defined reference category of the grouped-race covariate is White).",
      notes              = "Exponential effect exp(0.182 * RACE_ASIAN_OTH) on AZ5104 CL/F, i.e. ~20.0 percent higher AZ5104 apparent clearance and ~16.6 percent lower AZ5104 AUCss than White patients. Maps onto Johnson 2025's 'Asian (non-Chinese, non-Japanese)' grouped-race category. The dominant reference grouping here is White (not Chinese), matching the sibling model Brown_2017_osimertinib.R. This category is one arm of the paper's 'extreme case' covariate combination (43 kg, 29.3 g/L albumin, Asian non-Chinese non-Japanese).",
      source_name        = "Grouped race: Asian (non-Chinese, non-Japanese) (Johnson 2025 Table 2)"
    ),
    RACE_OTHER = list(
      description        = "Grouped-race category 'Other' indicator (1 = race grouped as Other, i.e. non-Asian and non-White, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the paper-defined reference category of the grouped-race covariate is White).",
      notes              = "Exponential effect exp(0.090 * RACE_OTHER) on AZ5104 CL/F, i.e. ~9.4 percent higher AZ5104 apparent clearance and ~8.6 percent lower AZ5104 AUCss than White patients. Johnson 2025 Table 2 labels this row 'Non-Asian, non-White population effect on CLmetabolite/F'; it is the 'Other' arm of the five-level grouped-race covariate in Table 1 (which pools Black/African American, Native Hawaiian/Pacific Islander, American Indian/Alaska Native, and Other, n = 145, 8.6 percent overall). Same mapping as the sibling model Brown_2017_osimertinib.R, whose corresponding category is 'non-Asian non-Caucasian'.",
      source_name        = "Grouped race: Other / 'Non-Asian, non-White' (Johnson 2025 Table 2)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the full covariate model but not retained: Johnson 2025 Discussion states age was 'tested and not found to have an impact on the PK that would require [inclusion] in the final model', and the Conclusion states no dose adjustment is required for age. No point estimate is reported, so the effect cannot be reconstructed."
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained (Johnson 2025 Discussion and Conclusion); no point estimate reported."
    ),
    SMOKER = list(
      description = "Current-smoker indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained (Johnson 2025 Discussion and Conclusion); no point estimate reported. Table 1 reports never / current / former smoking status."
    ),
    CRCL = list(
      description = "Baseline creatinine clearance (renal function).",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Renal impairment status was screened but not retained (Johnson 2025 Discussion); no point estimate reported."
    ),
    WHO_PS = list(
      description = "World Health Organization performance status (0 or 1).",
      units       = "(score)",
      type        = "categorical",
      notes       = "Screened but not retained (Johnson 2025 Discussion); no point estimate reported."
    ),
    LINE_OF_THERAPY = list(
      description = "Line of therapy (first-line, second-line, third-line onwards, adjuvant).",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened on both CLparent/F and CLmetabolite/F and explicitly rejected: Johnson 2025 Discussion states line of therapy 'was not statistically significant and that the median variability was clearly below the 20 percent threshold for clinical relevance'. No point estimate reported."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1364,
    n_studies      = 4,
    age_range      = "25-91 years",
    age_median     = "62 years",
    weight_range   = "29-122 kg",
    weight_median  = "61 kg (model-building set; 62 kg including the ADAURA validation set)",
    sex_female_pct = 65.3,
    race_ethnicity = c(
      White                            = 28.2,
      Asian_non_Chinese_non_Japanese   = 23.7,
      Chinese                          = 21.7,
      Japanese                         = 17.8,
      Other                            = 8.6
    ),
    disease_state  = "Advanced EGFR-mutation-positive non-small cell lung cancer (NSCLC). Model building used 599 patients from AURA / AURA extension (phase 1/2, NCT01802632), 210 from AURA2 (phase 2, NCT02094261), 277 from AURA3 (phase 3, NCT02151981), and 278 from FLAURA (phase 3, NCT02296125). A further 325 patients with resected stage IB-IIIA EGFR-mutation-positive NSCLC from the adjuvant ADAURA study (phase 3, NCT02511106) were held out for external validation and are NOT part of the 1364 model-building subjects.",
    dose_range     = "Osimertinib 20-240 mg once daily orally in the AURA phase 1 dose escalation (capsule at 20, 40, 80, 160, 240 mg and tablet at 80 mg); 80 mg once daily (tablet) in AURA extension, AURA2, AURA3, FLAURA, and ADAURA. Dose reductions were permitted in AURA2, AURA3, and FLAURA. 80 mg once daily is the recommended dose.",
    regions        = "Multiregional (AURA, AURA2, AURA3, FLAURA, and ADAURA enrolled across North America, Europe, and Asia).",
    n_observations = "41 461 plasma concentration samples from 1364 patients. Three patients (two in AURA extension, one in FLAURA) whose only two observations were both below the lower limit of quantification were removed before analysis. LLOQ = 0.05 nmol/L (osimertinib) and 0.0515 nmol/L (AZ5104).",
    hepatic_function = "1513 of 1689 (89.6 percent) normal, 166 (9.8 percent) at least mild impairment (Johnson 2025 Table 1).",
    renal_function   = "605 of 1689 (35.8 percent) normal, 742 (43.9 percent) at least mild impairment, 328 (19.4 percent) at least moderate impairment (Johnson 2025 Table 1).",
    notes          = "Demographic counts and percentages are the 'Overall (N = 1689)' column of Johnson 2025 Table 1, which pools the 1364 model-building patients with the 325 ADAURA external-validation patients; the paper does not print a model-building-only demographic column. Race percentages are the 'Grouped race' rows (Table 1 footnote a: the Asian population is split into Japanese, Chinese, and non-Japanese / non-Chinese Asian) because that is the categorisation the covariate model uses; the ungrouped 'Race' rows report 63.1 percent Asian and 34.1 percent White. Age and weight ranges are the Table 1 'Overall' ranges."
  )

  ini({
    # Structural parameters. Typical values are for the reference patient:
    # White, 61 kg body weight, 39 g/L baseline serum albumin
    # (Johnson 2025 Figure 2 caption).

    lka        <- log(0.196); label("First-order oral absorption rate constant for osimertinib (1/h)")                          # Johnson 2025 Table 2: ka = 0.196 1/h (RSE 5.09%)
    lcl        <- log(14.3);  label("Apparent osimertinib clearance, CLparent/F, in the reference patient (L/h)")               # Johnson 2025 Table 2: CLparent/F = 14.3 L/h (RSE 1.39%)
    lvc        <- log(918);   label("Apparent osimertinib central volume, Vparent/F, in the reference patient (L)")             # Johnson 2025 Table 2: Vparent/F = 918 L (RSE 3.30%)
    lcl_az5104 <- log(31.3);  label("Apparent AZ5104 clearance, CLmetabolite/F, in the reference patient (L/h)")                # Johnson 2025 Table 2: CLmetabolite/F = 31.3 L/h (RSE 2.00%)
    lvc_az5104 <- log(143);   label("Apparent AZ5104 central volume, Vmetabolite/F, in the reference patient (L)")              # Johnson 2025 Table 2: Vmetabolite/F = 143 L (RSE 4.22%)

    # Fraction of osimertinib converted to AZ5104. Held fixed, not
    # estimated: Johnson 2025 Discussion states "Similar to the previously
    # published model, the conversion of osimertinib to metabolite was
    # fixed at 25%." The value is arbitrary in the sense that it is
    # confounded with Vmetabolite/F and CLmetabolite/F in a joint
    # parent/metabolite model, so simulations are insensitive to it
    # provided the same value is used at simulation as at estimation.
    fm_az5104 <- fixed(0.25); label("Fraction of osimertinib converted to AZ5104 (unitless)")                                   # Johnson 2025 Discussion: conversion fixed at 25%

    # Continuous-covariate power exponents. Supplementary Appendix
    # equation 2.1: P_i = TP * (Cov_i / Ref)^beta * exp(eta_i), where Ref
    # is the median of the covariate in the analysis dataset. All
    # estimated (each has a nonzero RSE in Table 2).
    e_wt_cl          <- 0.421;  label("Power exponent for body weight on CLparent/F (unitless)")                                # Johnson 2025 Table 2: effect of baseline body weight on CLparent/F = 0.421 (RSE 13.9%)
    e_alb_cl         <- 0.825;  label("Power exponent for baseline albumin on CLparent/F (unitless)")                           # Johnson 2025 Table 2: effect of baseline albumin on CLparent/F = 0.825 (RSE 11.7%)
    e_wt_vc          <- 0.814;  label("Power exponent for body weight on Vparent/F (unitless)")                                 # Johnson 2025 Table 2: effect of baseline body weight on Vparent/F = 0.814 (RSE 16.9%)
    e_alb_vc         <- 2.27;   label("Power exponent for baseline albumin on Vparent/F (unitless)")                            # Johnson 2025 Table 2: effect of baseline albumin on Vparent/F = 2.27 (RSE 9.82%)
    e_wt_cl_az5104   <- 0.822;  label("Power exponent for body weight on CLmetabolite/F (unitless)")                            # Johnson 2025 Table 2: effect of baseline body weight on CLmetabolite/F = 0.822 (RSE 8.01%)
    e_alb_cl_az5104  <- 0.928;  label("Power exponent for baseline albumin on CLmetabolite/F (unitless)")                       # Johnson 2025 Table 2: effect of baseline albumin on CLmetabolite/F = 0.928 (RSE 11.4%)
    e_alb_vc_az5104  <- -0.831; label("Power exponent for baseline albumin on Vmetabolite/F (unitless)")                        # Johnson 2025 Table 2: effect of baseline albumin on Vmetabolite/F = -0.831 (RSE 28.6%)

    # Categorical-covariate exponential coefficients. Supplementary
    # Appendix equation 2.2: P_i = TP * exp(beta * 1[CAT_i = x]) *
    # exp(eta_i). White is the reference grouped-race category, so all
    # four indicators are 0 for a White patient. Race enters CLmetabolite/F
    # only; Johnson 2025 Results section 3.3 states "There was no effect of
    # race on the apparent clearance of osimertinib."
    e_race_asian_oth_cl_az5104 <- 0.182; label("Exponential coefficient for Asian-other vs White on CLmetabolite/F (unitless)")     # Johnson 2025 Table 2: Asian (non-Chinese, non-Japanese) population effect on CLmetabolite/F = 0.182 (RSE 11.2%)
    e_race_chinese_cl_az5104   <- 0.076; label("Exponential coefficient for Chinese vs White on CLmetabolite/F (unitless)")         # Johnson 2025 Table 2: Chinese population effect on CLmetabolite/F = 0.076 (RSE 27.2%)
    e_race_japanese_cl_az5104  <- 0.184; label("Exponential coefficient for Japanese vs White on CLmetabolite/F (unitless)")        # Johnson 2025 Table 2: Japanese population effect on CLmetabolite/F = 0.184 (RSE 13.1%)
    e_race_other_cl_az5104     <- 0.090; label("Exponential coefficient for non-Asian non-White vs White on CLmetabolite/F (unitless)")  # Johnson 2025 Table 2: Non-Asian, non-White population effect on CLmetabolite/F = 0.090 (RSE 26.4%)

    # Between-patient variability. Johnson 2025 Table 2 reports these under
    # the heading "Between-patient variability (% coefficient of
    # variation)", so each tabulated percentage is a coefficient of
    # variation, not a log-scale SD. ini() takes log-scale variances, so
    # each entry below is omega^2 = log(1 + CV^2):
    #   CLparent/F      44.9% -> log(1 + 0.449^2) = 0.183656  (omega 0.4286)
    #   CLmetabolite/F  49.7% -> log(1 + 0.497^2) = 0.220741  (omega 0.4698)
    #   Vparent/F       90.9% -> log(1 + 0.909^2) = 0.602425  (omega 0.7762)
    #   Vmetabolite/F   78.2% -> log(1 + 0.782^2) = 0.477281  (omega 0.6909)
    #   ka             109%   -> log(1 + 1.09^2)  = 0.783073  (omega 0.8849)
    # The CV reading (rather than reading the tabulated percentages
    # directly as log-scale SDs) is confirmed by the 90% prediction
    # intervals plotted in Johnson 2025 Figure 2: the typical-patient AUCss
    # interval spans 0.50-2.04 fold, implying omega = ln(2.04/0.50)/(2 *
    # 1.645) = 0.425, which matches 0.4286 (CV reading) and not 0.449
    # (log-SD reading). See the validation vignette for the full check.
    #
    # The off-diagonal covariance is rho * omega_CLparent *
    # omega_CLmetabolite = 0.885 * 0.428551 * 0.469831 = 0.178195.
    etalcl + etalcl_az5104 ~ c(0.183656, 0.178195, 0.220741)  # Johnson 2025 Table 2: Omega(CLparent/F) 44.9% CV, Omega(CLmetabolite/F) 49.7% CV, Correlation(CLparent/F, CLmetabolite/F) 0.885
    etalvc        ~ 0.602425                                  # Johnson 2025 Table 2: Omega(Vparent/F) = 90.9% CV
    etalvc_az5104 ~ 0.477281                                  # Johnson 2025 Table 2: Omega(Vmetabolite/F) = 78.2% CV
    etalka        ~ 0.783073                                  # Johnson 2025 Table 2: Omega(ka) = 109% CV

    # Residual error. Johnson 2025 Table 2 reports a separate combined
    # proportional-plus-additive error model for each analyte, with the
    # additive terms in nmol/L (the assay scale; LLOQ 0.05 nmol/L for
    # osimertinib and 0.0515 nmol/L for AZ5104, Methods section 2.2).
    # This model works in mass units, so each additive term is converted
    # to mg/L as nmol/L * molecular weight (g/mol) / 1e6. Molecular
    # weights are not reported in Johnson 2025 and are computed from the
    # published molecular formulae (see the mw_parent / mw_az5104
    # constants in model()).
    propSd        <- 0.205;      label("Proportional residual error on osimertinib (fraction)")                                          # Johnson 2025 Table 2: proportional error osimertinib = 0.205 (RSE 0.424%)
    addSd         <- 0.0150386;  label("Additive residual error on osimertinib (mg/L) = 30.1 nmol/L * 499.62 g/mol / 1e6")               # Johnson 2025 Table 2: additive error osimertinib = 30.1 nmol/L (RSE 1.00%)
    propSd_az5104 <- 0.215;      label("Proportional residual error on AZ5104 (fraction)")                                               # Johnson 2025 Table 2: proportional error AZ5104 = 0.215 (RSE 0.357%)
    addSd_az5104  <- 2.5056e-4;  label("Additive residual error on AZ5104 (mg/L) = 0.516 nmol/L * 485.59 g/mol / 1e6")                   # Johnson 2025 Table 2: additive error AZ5104 = 0.516 nmol/L (RSE 2.48%)
  })

  model({
    # Molecular weights, used only to (a) render the 1:1 molar
    # parent -> AZ5104 conversion in mass units and (b) convert the
    # nmol/L additive residual errors of Table 2 to mg/L. NOT reported in
    # Johnson 2025; computed from the published molecular formulae using
    # IUPAC standard atomic weights, and identical to the values used by
    # the sibling model Brown_2017_osimertinib.R.
    #   osimertinib C28H33N7O2 = 28*12.011 + 33*1.008 + 7*14.007 + 2*15.999 = 499.62 g/mol
    #   AZ5104 (N-desmethyl-osimertinib) C27H31N7O2                          = 485.59 g/mol
    mw_parent <- 499.62
    mw_az5104 <- 485.59

    # Reference covariate values: the medians of the model-building
    # analysis dataset, per Johnson 2025 Results section 3.3 and the
    # Figure 2 caption ("Typical patient: White, 61 kg body weight,
    # 39 g/L baseline serum albumin").
    ref_wt  <- 61
    ref_alb <- 39

    # Individual parameters. Continuous covariates enter in the
    # median-normalized power form of Supplementary Appendix equation 2.1;
    # the grouped-race indicators enter in the exponential form of
    # equation 2.2 with White as the reference category.
    ka <- exp(lka + etalka)

    cl <- exp(lcl + etalcl) *
          (WT / ref_wt)^e_wt_cl *
          (ALB / ref_alb)^e_alb_cl

    vc <- exp(lvc + etalvc) *
          (WT / ref_wt)^e_wt_vc *
          (ALB / ref_alb)^e_alb_vc

    cl_az5104 <- exp(lcl_az5104 + etalcl_az5104) *
                 (WT / ref_wt)^e_wt_cl_az5104 *
                 (ALB / ref_alb)^e_alb_cl_az5104 *
                 exp(e_race_asian_oth_cl_az5104 * RACE_ASIAN_OTH +
                     e_race_chinese_cl_az5104   * RACE_CHINESE +
                     e_race_japanese_cl_az5104  * RACE_JAPANESE +
                     e_race_other_cl_az5104     * RACE_OTHER)

    vc_az5104 <- exp(lvc_az5104 + etalvc_az5104) *
                 (ALB / ref_alb)^e_alb_vc_az5104

    kel        <- cl / vc
    kel_az5104 <- cl_az5104 / vc_az5104

    # ODE system. Amounts are in mg; the 1:1 molar parent -> AZ5104
    # stoichiometry is rendered in mass units via mw_az5104 / mw_parent so
    # that central_az5104 holds true AZ5104 mass and Cc_az5104 is a true
    # AZ5104 mass concentration.
    d/dt(depot)          <- -ka * depot
    d/dt(central)        <-  ka * depot - kel * central
    d/dt(central_az5104) <-  fm_az5104 * kel * central * (mw_az5104 / mw_parent) -
                             kel_az5104 * central_az5104

    # Observations (mg/L). Johnson 2025 reports concentrations in nmol/L;
    # multiply by 1e6 / molecular weight to convert.
    Cc        <- central / vc
    Cc_az5104 <- central_az5104 / vc_az5104

    Cc        ~ prop(propSd) + add(addSd)
    Cc_az5104 ~ prop(propSd_az5104) + add(addSd_az5104)
  })
}
