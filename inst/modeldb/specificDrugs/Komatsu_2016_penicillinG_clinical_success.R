Komatsu_2016_penicillinG_clinical_success <- function() {
  description <- paste0(
    "Static landmark logistic exposure-response model for clinical success ",
    "of penicillin G therapy in infective endocarditis caused by viridans ",
    "group streptococci (Komatsu 2016, n = 21 of the 25 patients in the ",
    "companion population PK analysis -- the 21 with an isolate and a ",
    "measured MIC). The probability of a positive clinical outcome is ",
    "expit(logit_ref + e_cminmic_success * CTROUGH / mic), i.e. the ",
    "published form 1 / {1 + exp(1.609 - 0.0524 x <penicillin G ",
    "Cmin/MIC>)}, where the exposure driver is the RATIO of the trough ",
    "serum penicillin G concentration to the MIC of the infecting isolate. ",
    "There is no PK layer and no ODE: the trough is supplied as the ",
    "CTROUGH data column, which a user can generate from the companion ",
    "population PK model Komatsu_2016_penicillinG. The MIC is carried as a ",
    "fixed model parameter rather than as a data column so that the model ",
    "can be re-targeted to an isolate of different susceptibility; it ",
    "defaults to 0.06 ug/mL, the lower of the two values Komatsu 2016 fixed ",
    "in its dosing simulations ('because these MICs are seen frequently in ",
    "our hospital'), with 0.12 ug/mL the other. The model was fitted by ",
    "logistic regression in JMP 6.03, separately from the NONMEM ",
    "population PK fit, so it is packaged as its own file. Note that the ",
    "paper's operative decision rule is NOT this curve but the ROC cut-off ",
    "derived from it: a Cmin/MIC ratio of 60, at which this logistic ",
    "returns 0.82 but at which every observed patient in the cohort ",
    "responded (sensitivity 68 %, specificity 100 %)."
  )
  reference <- paste(
    "Komatsu T, Inomata T, Watanabe I, Kobayashi M, Kokubun H, Ako J, Atsuda K.",
    "Population pharmacokinetic analysis and dosing regimen optimization of",
    "penicillin G in patients with infective endocarditis.",
    "J Pharm Health Care Sci. 2016 Apr 5;2:9.",
    "doi:10.1186/s40780-016-0043-x. PMCID PMC4820900.",
    sep = " "
  )
  vignette <- "Komatsu_2016_penicillinG"

  units <- list(
    time          = "n/a (static landmark exposure-response model; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the CTROUGH covariate column)",
    concentration = "CTROUGH and mic in ug/mL; prob_clinical_success is a probability (0-1)"
  )

  covariateData <- list(
    CTROUGH = list(
      description        = "Minimum (pre-dose trough) serum penicillin G concentration over the dosing interval. TOTAL, not unbound.",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Komatsu 2016 writes this quantity 'penicillin G Cmin' and always ",
        "uses it divided by the isolate MIC; 'Cmin' is a registered source ",
        "alias of CTROUGH. Landmark: the measured trough during penicillin ",
        "G therapy -- Methods, Data source: 'Blood samples were obtained ",
        "from a brachial vein immediately before and 2 or 3 h after ",
        "administration of penicillin G'. The fitted values are therefore ",
        "MEASURED troughs, not empirical-Bayes predictions from the ",
        "companion population PK model; Table 2 reports the observed trough ",
        "as 7.8 +/- 9.0 ug/mL overall, 10.1 +/- 9.8 in the 15 patients who ",
        "responded and 2.2 +/- 1.2 in the 6 who failed. To simulate this ",
        "column, take the end-of-interval concentration of ",
        "modellib('Komatsu_2016_penicillinG'). TOTAL concentration: the ",
        "HPLC assay of Methods, 'Assay of penicillin G concentrations', ",
        "measures drug in serum after a methanol protein precipitation, so ",
        "it is not a free-fraction measurement even though only unbound ",
        "penicillin G is microbiologically active."
      ),
      source_name        = "Cmin"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 21L,
    n_studies      = 1L,
    n_observations = "21 binary outcome records, one per patient (landmark analysis, no repeated measures)",
    age_range      = "21-83 years (the parent cohort of 25; Komatsu 2016 Table 1)",
    sex_female_pct = 36,
    race_ethnicity = "Not reported; single-centre Japanese cohort (Kitasato University Hospital, Sagamihara)",
    disease_state  = "Infective endocarditis with viridans group streptococci isolated. Species distribution (Komatsu 2016 Table 2): Streptococcus sanguis 5, S. gordonii 4, species not determined 4, S. agalactiae 2, S. intermedius 2, S. mutans 1, S. mitis 1, S. oralis 1, S. constellatus 1. Observed penicillin G MIC 0.11 +/- 0.20 ug/mL overall, 0.07 +/- 0.02 in responders and 0.21 +/- 0.38 in failures. Vegetation size 11.7 +/- 6.0 mm.",
    dose_range     = "Penicillin G potassium intravenously; individual clinical regimens are not tabulated",
    regions        = "Japan (Kitasato University Hospital, Sagamihara, Kanagawa); patients treated between January 1997 and April 2013",
    notes          = paste0(
      "Outcome definition, Komatsu 2016 Methods, Data source: 'Failure of ",
      "penicillin G treatment was defined as persistence of fever and/or ",
      "bacteremia by the causative pathogen requiring a change in ",
      "antibiotic therapy, and infection-related mortality within 30 ",
      "days.' Clinical success = 1, failure = 0. 15 of 21 succeeded, a ",
      "29 % failure rate that the Discussion compares with the 30 % ",
      "mortality rate reported by Garcia-Cabrera et al. The analysis set ",
      "is 21, not 25, because 4 of the 25 patients in the PK analysis had ",
      "no viridans-streptococcal isolate and therefore no MIC."
    )
  )

  ini({
    # ==================================================================
    # Komatsu 2016 Results: 'Penicillin G Cmin/MIC was a significant
    # predictor of the following clinical outcome equation: probability
    # of a positive clinical outcome = 1/{1 + exp(1.609 - 0.0524 x
    # <penicillin G Cmin/MIC>)}'. The same equation is restated in the
    # Determination-of-the-dosing-regimen column of page 5.
    #
    # 1/(1 + exp(a - b*x)) == expit(-a + b*x), so the printed intercept
    # +1.609 enters the linear predictor with its sign REVERSED and the
    # printed slope 0.0524 keeps its sign. Sign check against the data:
    # Table 2 reports a mean trough of 10.1 ug/mL in responders against
    # 2.2 ug/mL in failures, so the probability must RISE with exposure,
    # which requires a positive slope on the ratio -- as encoded.
    #
    # The source reports no standard errors, no confidence intervals and
    # no between-subject random effect for this regression; JMP fits an
    # exact Bernoulli likelihood. Neither coefficient is FIXED in the
    # source sense -- both were estimated -- so neither is wrapped in
    # fixed().
    # ==================================================================
    logit_ref          <- -1.609;  label("Logit of clinical success at a penicillin G Cmin/MIC ratio of zero (unitless logit)")               # Komatsu 2016 Results, clinical outcome equation: the printed intercept is +1.609 inside exp(1.609 - 0.0524 x ratio), which is -1.609 on the logit scale
    e_cminmic_success  <-  0.0524; label("Log-odds of clinical success per one-unit increase in the penicillin G Cmin/MIC ratio (unitless logit)")  # Komatsu 2016 Results, clinical outcome equation: the printed slope 0.0524 inside exp(1.609 - 0.0524 x ratio)

    # ==================================================================
    # MIC of the infecting isolate. Carried as a model parameter, not a
    # data column, following Chen_2023_tilmicosin.R, Beredaki_2023_*
    # and Olivo_2025_vancomycin_invitro.R: the AUCMIC_TYLO register
    # entry directs that a paper which reports the isolate MIC should
    # split the PK/PD index into an exposure column divided by a model
    # `mic` parameter, so that the model can be re-targeted to an
    # isolate of different susceptibility. Komatsu 2016 does report
    # MICs (Table 2) and fixes the value in its simulations.
    #
    # Default 0.06 ug/mL is the first of the two values used for the
    # Figure 4 target-attainment panels and the Figure 5 nomogram; set
    # mic to 0.12 to reproduce panel (b) and the right-hand branch of
    # the nomogram, or to a measured isolate MIC for clinical use.
    # ==================================================================
    mic <- fixed(0.06); label("Minimum inhibitory concentration of the infecting isolate (ug/mL)")  # Komatsu 2016 Determination of the dosing regimen: 'The MIC was fixed at 0.06 and 0.12 mg/L because these MICs are seen frequently in our hospital'

    # ==================================================================
    # The source likelihood is Bernoulli, so there is no residual error
    # to transcribe. rxode2 requires an observation declaration, so the
    # deterministic probability is emitted with a tiny placeholder
    # additive residual, mirroring Babel_2026_telisotuzumab_orr.R and
    # the Fukae_2024_valemetostat_* family.
    # ==================================================================
    addSd_prob_clinical_success <- fixed(0.001); label("Placeholder additive residual SD on the typical-value clinical-success probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # The PK/PD index the regression was fitted on. CTROUGH is the
    # measured trough serum penicillin G concentration in ug/mL and mic
    # is in the same units, so the ratio is dimensionless. The paper's
    # efficacy target is cminmic > 60.
    cminmic <- CTROUGH / mic

    # Linear predictor. The ratio is UNCENTRED, so logit_ref is the
    # logit at a ratio of zero (probability 0.167) rather than at a
    # typical patient's exposure.
    logit_clinical_success <- logit_ref + e_cminmic_success * cminmic
    prob_clinical_success  <- expit(logit_clinical_success)

    # Deterministic probability of a positive clinical outcome.
    # Downstream callers can sample binary outcomes with
    # rbinom(n, 1, prob_clinical_success) on the rxSolve output.
    prob_clinical_success ~ add(addSd_prob_clinical_success)
  })
}
