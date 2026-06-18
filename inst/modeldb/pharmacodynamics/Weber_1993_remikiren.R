Weber_1993_remikiren <- function() {
  description <- paste(
    "Population pharmacodynamic Hill / Emax model for the orally active renin inhibitor remikiren (Ro 42-5892) in patients with mild-to-moderate essential hypertension.",
    "Inhibition of the angiotensin I production rate (APR; the net inhibition of plasma renin activity corrected for changes in immunoreactive renin) is described as a saturable function of observed plasma remikiren concentration via the Hill equation with the Hill coefficient fixed at 1.",
    "PD-only model: plasma remikiren concentration is supplied as a time-varying covariate CP_REM_NGML (ng/mL).",
    "The source publication characterised remikiren PK with model-independent NCA only (Cmax / tmax / AUC0-t in Table 1) and did not develop a structural population PK model, so the PD model has no coupled PK component.",
    "Population: 144 patients with mild-to-moderate essential hypertension across three multi-dose clinical pharmacology studies (oral solution or 100 mg capsules; 100-800 mg po qd for 8 days)."
  )
  reference <- paste(
    "Weber C, Birnbock H, Leube J, Kobrin I, Kleinbloesem CH, van Brummelen P (1993).",
    "Multiple dose pharmacokinetics and concentration effect relationship of the orally active renin inhibitor remikiren (Ro 42-5892) in hypertensive patients.",
    "Br J Clin Pharmacol 36(6):547-555.",
    "doi:10.1111/j.1365-2125.1993.tb00413.x.",
    sep = " "
  )
  vignette <- "Weber_1993_remikiren"
  units <- list(
    time = "h",
    dosing = "(none; PD-only model fed by an external remikiren plasma-concentration covariate)",
    concentration = "(observation APR is fractional inhibition of the angiotensin I production rate, unitless 0-1; the source paper reports % which is APR * 100; driving covariate CP_REM_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_REM_NGML = list(
      description        = "Instantaneous remikiren (Ro 42-5892) plasma concentration at the time of each PD observation, supplied as a time-varying covariate from observed plasma samples or an upstream PK source.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying per event row. Drives the Hill / Emax expression APR = Imax * CP_REM_NGML^hill / (IC50^hill + CP_REM_NGML^hill).",
        "In Weber 1993 this was the observed plasma remikiren concentration measured by HPLC with fluorescence or coulometric electrochemical detection or by a protein-binding assay (Methods, 'Remikiren concentrations' section); the source NONMEM column name is Cp.",
        "Reference values observed: mean Cmax was 4-6 ng/mL after 200 mg oral, 23-27 ng/mL after 300 mg oral, 47-83 ng/mL after 600-800 mg oral, with substantial intersubject variability (Weber 1993 Table 1).",
        "Set to 0 outside the drug-exposure window (the inhibition term then collapses to 0)."
      ),
      source_name        = "Cp"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 144L,
    n_studies      = 3L,
    age_range      = "24-69 years across the three studies (Study A: 30-64 y; Study B: 24-69 y; Study C: 38-64 y)",
    age_median     = NA_character_,
    weight_range   = NA_character_,
    weight_median  = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Mild-to-moderate essential hypertension (entry sitting diastolic BP in the range 95-115 mm Hg depending on study)",
    dose_range     = paste(
      "Study A: 100, 200, 400, or 800 mg po qd for 8 days as oral solution in orange juice, always on an empty stomach (n = 47 evaluable on active treatment across four dose levels).",
      "Study B: 300 or 600 mg po qd for 8 days as 100 mg capsules on an empty stomach (n = 17 per active dose level), with a 100 mg i.v. dose 4 h after the last oral dose (n = 16 placebo).",
      "Study C: 600 mg po qd for 8 days as 100 mg capsules 1-2 h after breakfast (n = 12 active, n = 12 placebo), with a 100 mg i.v. dose or i.v. placebo 4 h after the last oral dose."
    ),
    regions        = NA_character_,
    notes          = paste(
      "Three double-blind, randomized, placebo-controlled, parallel-group clinical pharmacology studies conducted by F. Hoffmann-La Roche (Basel, Switzerland).",
      "Study A: 70 male volunteers enrolled; Study B: 51 volunteers (26 male, 25 female); Study C: 24 male volunteers.",
      "Plasma renin activity (PRA) and immunoreactive renin (IRR) were measured up to 24 h post-dose on the first and last days of oral treatment, alongside plasma remikiren concentrations.",
      "PD analysis was performed using NONMEM III level 1.2 (Weber 1993 Methods, 'Concentration-effect modelling').",
      "Random effects were 'linked to the fixed effects as constant coefficients of variance' (CCV parameterisation in NONMEM); reported IIV values are CV%. This extraction translates the reported CV% to log-normal omega^2 via omega^2 = log(1 + CV^2)."
    )
  )

  ini({
    # Hill / Emax PD model fit by Weber 1993 to inhibition of the angiotensin I
    # production rate (APR) versus observed plasma remikiren concentration.
    # The paper reports that the Hill coefficient was set to 1, so the
    # working form is the simple Emax: APR = Imax * Cp / (IC50 + Cp).
    # Typical-value parameters from Weber 1993 Results, paragraph 5 (p. 552).

    limax <- log(0.903)
    label("Maximum fractional inhibition of angiotensin I production rate (Imax, unitless, 0-1)")
    # Weber 1993, Results p. 552: 'Imax was estimated as 90.3% (CV = 1.4%)'.
    # The 1.4% CV is the parameter precision (RSE), not IIV.

    lic50 <- log(0.5)
    label("Plasma remikiren concentration producing 50% of maximum APR inhibition (IC50, ng/mL)")
    # Weber 1993, Results p. 552: 'IC50 was 0.5 ng ml-1 (CV = 32.3%)'.
    # The 32.3% CV is parameter precision; the value also matches the
    # in-vitro IC50 of 0.5 ng/mL (0.8 nM) reported in the Discussion p. 553.

    lhill <- fixed(log(1))
    label("Hill coefficient governing sigmoidicity (FIXED to 1 per Weber 1993)")
    # Weber 1993, Results p. 552 ('the Emax-model with the Hill coefficient
    # set to 1') and Methods, 'Concentration-effect modelling'.

    # Interindividual variability (IIV). Weber 1993 fit the random effects
    # 'as constant coefficients of variance' in NONMEM III (Methods,
    # 'Concentration-effect modelling') and reports IIV as a percentage CV.
    # Translated to nlmixr2 log-normal omega^2 via the standard small-CV
    # approximation omega^2 = log(1 + CV^2):
    #   Imax IIV  7.4% -> omega^2 = log(1 + 0.074^2) = 0.005461
    #   IC50 IIV 111%  -> omega^2 = log(1 + 1.110^2) = 0.802923
    etalimax ~ 0.005461
    # Weber 1993, Results p. 552: 'Its interindividual variability of 7.4%
    # was very small' (Imax).

    etalic50 ~ 0.802923
    # Weber 1993, Results p. 552: 'Its interindividual variability was 111%'
    # (IC50).

    # Residual error. The source publication does not tabulate the
    # residual SD / SIGMA -- only Imax, IC50, and their IIV CV% are
    # reported. The APR observation is a fractional inhibition (0-1)
    # so a proportional error structure on APR is the natural choice.
    # The placeholder magnitude of 0.30 (30% CV) below is documented in
    # the validation vignette Errata as non-paper-derived; typical-value
    # predictions are independent of this choice, only the VPC envelope
    # is affected.
    propSd <- 0.30
    label("Proportional residual error on APR-Inhibition (placeholder; Weber 1993 does not tabulate the residual SD)")
  })

  model({
    # Per-subject Imax and IC50 via log-normal IIV.
    imax <- exp(limax + etalimax)
    ic50 <- exp(lic50 + etalic50)
    hill <- exp(lhill)

    # Hill / Emax model for inhibition of the angiotensin I production rate.
    # CP_REM_NGML is supplied as a time-varying covariate (ng/mL). APR is
    # the fractional inhibition (0-1 scale; the source paper reports % which
    # is APR * 100).
    APR <- imax * CP_REM_NGML^hill / (ic50^hill + CP_REM_NGML^hill)
    APR ~ prop(propSd)
  })
}
