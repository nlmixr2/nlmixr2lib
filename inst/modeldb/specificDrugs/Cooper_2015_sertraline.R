Cooper_2015_sertraline <- function() {
  description <- "One-compartment first-order absorption population PK model for sertraline in overdose (Cooper 2015). Apparent clearance is increased 1.92-fold in subjects who received single-dose activated charcoal; the model holds relative bioavailability F at 1 and a shifted lag time at 1 h, with between-subject variability on F, ts_lag, ka, Vc, and CL absorbing the overdose-specific dose-amount and dose-time uncertainty."
  reference   <- "Cooper JM, Duffull SB, Saiao AS, Isbister GK. The pharmacokinetics of sertraline in overdose and the effect of activated charcoal. Br J Clin Pharmacol. 2015 May;79(5):307-15. doi:10.1111/bcp.12500"
  vignette    <- "Cooper_2015_sertraline"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CONMED_CHARCOAL = list(
      description        = "Single-dose activated charcoal administered for gut decontamination",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no activated charcoal)",
      notes              = "Cooper 2015 Table 1 / Methods: 7/28 patients (25%) received single-dose activated charcoal (Carbosorb, 50 g suspension) between 1.5 and 4 h post-overdose (median 3 h). Treated as a per-subject time-fixed binary indicator in the source data; the paper does not time-resolve the SDAC bolus within subject, so the fractional CL effect is applied uniformly across the observation record. Enters the structural model via the paper's final formula CL = theta_CL * f_CL-char, where f_CL-char = 1.92 when CONMED_CHARCOAL = 1 (Table 2, Model 2 final column) and 1 when CONMED_CHARCOAL = 0. The paper also evaluated an effect on relative bioavailability F (Model 1 f_F-char = 0.731) but the f_CL-char form was selected as the final model on objective-function and biological-plausibility grounds (significant effect on CL with P < 0.05).",
      source_name        = "SDAC"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 28L,
    n_studies        = 1L,
    age_range        = "15-55 years",
    age_median       = "32 years",
    sex_female_pct   = 75,
    race_ethnicity   = "Not reported (single regional toxicology unit, Hunter Region, New South Wales, Australia)",
    disease_state    = "Acute sertraline overdose presenting to a regional toxicology unit between February 2001 and February 2010; 21/28 (75%) co-ingested other substances (alcohol most commonly, then analgesics, antihistamines, antipsychotics, benzodiazepines); none of the co-ingestants known to inhibit or induce sertraline metabolism. 7/28 patients developed serotonin toxicity; 4/28 had a Glasgow Coma Score < 15. No deaths or major complications.",
    dose_range       = "250-5000 mg single self-administered oral overdose (median 1550 mg)",
    regions          = "Australia (Hunter Region, New South Wales)",
    co_medication    = "10/28 alcohol; 6/28 analgesics; 3/28 antihistamines; 3/28 antipsychotics; 4/28 benzodiazepines; 2/28 antihypertensives; 2/28 decongestants; 1/28 each: warfarin, sodium valproate, valerian, unknown substance; 9/28 no co-ingestants",
    n_concentrations = 77L,
    notes            = "Baseline demographics from Cooper 2015 Table 1. Reported overdose was corroborated where possible by friends / relatives / ambulance officer record / tablet counts and assigned a 5-point veracity score; 14/28 score 1 (best), 11/28 score 2, 3/28 score 3, 0/28 score 4. The veracity score was tested as a covariate constraint on the between-subject variance of relative bioavailability but did not improve the model (paper Methods Uncertainty in overdose history) and was not retained in the final model. Median 2 concentrations per patient (range 1-6); sampling times 1.17-68 h post-overdose."
  )

  ini({
    # Structural parameters (Cooper 2015 Table 2, Model 2 'Final SDAC f_CL-char'
    # column). Cooper et al. fitted in Monolix v4.2 using SAEM; the typical-value
    # parameters below are taken from the final-model column (col 3 of Table 2).
    # The paper applied a +1 h time shift to every dose / observation record to
    # permit a 'negative' lag for subjects who ingested the overdose before
    # the reported time; ts_lag = 1 h fixed then recovers the typical reported-
    # time alignment, and BSV on ts_lag absorbs the per-subject dose-time
    # uncertainty (paper Methods 'Uncertainty in overdose history'). The
    # tlag canonical is used here for the lag-time semantic; the in-file
    # comments below document the dose-time-uncertainty interpretation.
    ltlag   <- fixed(log(1));    label("Shifted lag time -- dose-time-uncertainty anchor (h)")          # Cooper 2015 Table 2 final: t_s,lag = 1 (fixed)
    lka     <- log(0.895);       label("First-order absorption rate constant ka (1/h)")                  # Cooper 2015 Table 2 final: K_a = 0.895 h^-1 (RSE 42%)
    lvc     <- log(5340);        label("Apparent volume of distribution Vc/F (L)")                       # Cooper 2015 Table 2 final: V = 5340 L (RSE 18%)
    lcl     <- log(130);         label("Apparent oral clearance CL/F at CONMED_CHARCOAL = 0 (L/h)")      # Cooper 2015 Table 2 final: theta_CL = 130 L/h (RSE 34%)
    lfdepot <- fixed(log(1));    label("Relative bioavailability F (fraction) -- dose-amount-uncertainty anchor") # Cooper 2015 Table 2 final: F = 1 (fixed)

    # Covariate effect of single-dose activated charcoal on CL (final-model
    # parameterisation CL = theta_CL * f_CL-char from Cooper 2015 Eq.).
    e_charcoal_cl <- 1.92;       label("Fractional effect on CL when CONMED_CHARCOAL = 1 (f_CL-char)")   # Cooper 2015 Table 2 final: f_CL-char = 1.92 (RSE 70%)

    # Between-subject variability (Cooper 2015 Table 2 final-model column,
    # 'Between subject variance (omega)' rows). The paper Methods state Monolix
    # computes 'between subject variances' for the PK parameters and the table
    # header explicitly labels these as variances (omega^2), so the values are
    # used directly as omega^2 here.
    etaltlag   ~ 0.922   # BSV on ts_lag (variance 0.922; RSE 52%)
    etalka     ~ 1.02    # BSV on ka     (variance 1.02;  RSE 77%)
    etalvc     ~ 0.085   # BSV on Vc     (variance 0.085; RSE 219%)
    etalcl     ~ 0.126   # BSV on CL     (variance 0.126; RSE 103%)
    etalfdepot ~ 0.303   # BSV on F      (variance 0.303; RSE 48%)

    # Residual error (Cooper 2015 Table 2 final-model column: CV% = 11.7 for the
    # proportional residual error).
    propSd <- 0.117;             label("Proportional residual error (fraction)")                         # Cooper 2015 Table 2 final: residual CV = 11.7%
  })
  model({
    # Individual parameters. The lag-time canonical name tlag carries the
    # paper's 'shifted lag time' interpretation; BSV on it absorbs per-subject
    # overdose-time uncertainty.
    tlag   <- exp(ltlag + etaltlag)
    ka     <- exp(lka   + etalka)
    vc     <- exp(lvc   + etalvc)
    fdepot <- exp(lfdepot + etalfdepot)
    cl     <- exp(lcl   + etalcl) * (1 + (e_charcoal_cl - 1) * CONMED_CHARCOAL)

    kel <- cl / vc

    # One-compartment first-order absorption with linear elimination.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Bioavailability F (fixed-typical = 1, BSV-only) and shifted lag-time.
    f(depot)    <- fdepot
    alag(depot) <- tlag

    # Sertraline assay reports concentrations in mg/L. Dose in mg, vc in L ->
    # central/vc has units mg/L directly.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
