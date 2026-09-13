Sunnaker_2026_mitiperstat <- function() {
  description <- "Two-compartment population PK model with first-order absorption and linear elimination for the myeloperoxidase inhibitor mitiperstat (AZD4831) in healthy volunteers, patients with heart failure with preserved or mildly reduced ejection fraction, and patients with severe renal impairment (Sunnaker 2026). Apparent clearance is scaled by baseline BSA-normalized eGFR and baseline body weight (power models) and shifted by Asian race and by heart-failure disease status (linear models); apparent central volume is scaled by age (power model). Bioavailability could not be estimated because no intravenous data were available, so all clearances and volumes are apparent (CL/F, Vc/F, Q/F, Vp/F)."
  reference <- "Sunnaker M, Leander J, Ericsson H. Population Pharmacokinetics of the Novel Myeloperoxidase Inhibitor Mitiperstat. Pharmacol Res Perspect. 2026;14(3):e70259. doi:10.1002/prp2.70259"
  vignette <- "Sunnaker_2026_mitiperstat"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Sunnaker 2026 reports every mitiperstat concentration in nmol/L and every
  # dose in mg, but never states the molar mass, so no exact mg <-> nmol
  # conversion is available from any on-disk source. The model is linear, so it
  # is encoded in self-consistent mass units: a dose in mg gives compartment
  # amounts in mg and Cc in mg/L (= ug/mL). Dosing the same model in nmol
  # returns Cc directly in nmol/L, which is the scale the paper prints.
  compartmentData <- list(
    depot       = list(analyte = "mitiperstat", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "mitiperstat", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "mitiperstat", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Baseline estimated glomerular filtration rate, normalized to a body surface area of 1.73 m^2",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Estimated with the CKD-EPI equation (Sunnaker 2026 Methods, Structural Base Model; Levey 2009 reference 17). Power effect on CL/F centered on 99 mL/min/1.73 m^2, the population median (Table 3 footnote). This was the only covariate carried in the base model, because renal excretion accounts for roughly 32-45 percent of mitiperstat elimination. Cohort means by study (Table 2): SAD 106, MAD 104, JCMAD 112, SATELLITE 69, renal-impairment cohort 23 and its group-matched controls 97. The paper re-estimated the final model with non-BSA-normalized eGFR (correlation 0.95 with the normalized form) and obtained similar parameter estimates.",
      source_name        = "baseline eGFR"
    ),
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline (time-fixed), not time-varying. Power effect on CL/F centered on 77.95 kg, the population median (Sunnaker 2026 Table 3 footnote). The exponent 0.78 was ESTIMATED, not held at an allometric 0.75; the paper explicitly tried fixed allometric scaling instead and reports that it slightly worsened the fit (Results, Final Model; Table S2). Baseline BMI correlates strongly with body weight (Pearson 0.8) and was therefore excluded from the covariate search; see covariatesDataExcluded.",
      source_name        = "baseline body weight"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on Vc/F centered on 41 years, the population median (Sunnaker 2026 Table 3 footnote). Age and eGFR are negatively correlated in this pooled data set (Pearson -0.7), and age was high only in the SATELLITE cohort, which was also the only cohort with heart failure; the paper states the Vc/F-age association should therefore be interpreted with caution (Discussion).",
      source_name        = "age"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian)",
      notes              = "Linear fractional increase of CL/F in Asian relative to non-Asian participants; the canonical 1 = Asian orientation matches the paper's coding, so no value flip is needed. Of the 26 Asian participants, 24 came from the JCMAD study (Japanese and Chinese volunteers, defined as having both parents and four grandparents of that ethnicity) and only 2 from the MAD study, so the paper cautions that the race effect cannot be cleanly separated from other between-study differences (Discussion, limitation 2).",
      source_name        = "race (Asian or non-Asian)"
    ),
    DIS_HFPEF = list(
      description        = "Heart failure with preserved or mildly reduced ejection fraction indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer, or a patient enrolled for renal impairment rather than heart failure)",
      notes              = "The paper's disease-status covariate, contrasting the 25 SATELLITE patients (symptomatic heart failure, left ventricular ejection fraction at or above 40 percent, elevated B-type natriuretic peptides) against the 103 participants without heart failure. Both the healthy volunteers of the SAD, MAD and JCMAD studies AND the severe-renal-impairment cohort take the value 0, because Table 2 classifies the renal-impairment participants as 'No HFpEF/HFmrEF'. The renal impairment of that cohort enters separately through CRCL, so the two covariates are not redundant.",
      source_name        = "disease status (healthy volunteers or patients with HFpEF/HFmrEF)"
    )
  )

  # Screened by the authors but NOT retained in the final model. Documented so
  # the provenance of the covariate search survives, without raising a
  # declared-but-unreferenced convention warning.
  covariatesDataExcluded <- list(
    BMI = list(
      description = "Baseline body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Excluded from the stepwise covariate search a priori because of its strong correlation with baseline body weight (Pearson 0.8); Sunnaker 2026 Results, Covariate Model. Cohort means by study (Table 2): SAD 24.3, MAD 25.2, JCMAD 23.3, SATELLITE 27.3, renal impairment 29.3."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Excluded from the stepwise covariate search a priori because only 22 of 128 participants were female and sex was confounded with both body weight and formulation - the only two studies that enrolled women (SATELLITE and renal impairment) were also the only two that used the tablet; Sunnaker 2026 Results, Covariate Model."
    ),
    FORM_TABLET = list(
      description = "Film-coated tablet versus oral suspension formulation indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested as a covariate on the absorption rate constant and found not significant (Sunnaker 2026 Results, Covariate Model and Discussion limitation 1). The reference oral liquid is the oral suspension used in the SAD, MAD and JCMAD studies; SATELLITE and the renal-impairment study used a film-coated tablet. The paper notes that absorption-phase data for the tablet came only from the renal-impairment study, so the power to detect a formulation effect on ka was low."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 128,
    n_studies      = 5,
    n_observations = 2856,
    age_range      = "18-85 years",
    age_median     = "not reported; study-level means 33.9-35.5 years in the healthy-volunteer studies, 57.1 years in the renal-impairment study and 75.2 years in SATELLITE (Table 2). The covariate model centers age at a population median of 41 years (Table 3 footnote).",
    weight_range   = "50-100 kg in the healthy-volunteer studies; at least 50 kg in the renal-impairment study; 54-113 kg observed in SATELLITE",
    weight_median  = "77.95 kg (the covariate-model centering value, Table 3 footnote)",
    sex_female_pct = 17.2,
    race_ethnicity = c(Asian = 20.3, `Non-Asian` = 79.7),
    disease_state  = "Healthy volunteers (83 participants), patients with heart failure with preserved or mildly reduced ejection fraction (25 participants, SATELLITE), and patients with severe renal impairment plus their group-matched normal-renal-function controls (20 participants)",
    renal_function = "Baseline eGFR study means 104-112 mL/min/1.73 m^2 in healthy volunteers, 69 mL/min/1.73 m^2 in SATELLITE, 23 mL/min/1.73 m^2 in the severe-renal-impairment cohort (eGFR at least 15 and below 30, not on dialysis) and 97 mL/min/1.73 m^2 in its group-matched controls",
    dose_range     = "Single oral doses of 2.5-405 mg; once-daily oral doses of 2.5-45 mg for 10-14 days, and 2.5 mg for 10 days uptitrated to 5 mg for a further 80 days in SATELLITE",
    formulation    = "Oral suspension in the SAD, MAD and JCMAD studies; film-coated tablet in SATELLITE and the renal-impairment study",
    regions        = "Not reported by region; the JCMAD study enrolled Japanese and Chinese volunteers, the remaining studies enrolled a predominantly non-Asian population",
    notes          = "Pooled from five trials: SAD NCT02712372, MAD NCT03136991, JCMAD NCT04232345, phase 2a SATELLITE NCT03756285 and the severe-renal-impairment study NCT04949438. Participant counts, demographics and baseline characteristics: Sunnaker 2026 Tables 1 and 2. Placebo recipients were excluded, as were 139 samples below the 2 nmol/L (0.2 nmol/L in the renal-impairment study) limit of quantification, 4.9 percent of the total."
  )

  ini({
    # Structural parameters - typical values for the reference participant:
    # eGFR 99 mL/min/1.73 m^2, body weight 77.95 kg, age 41 years, non-Asian,
    # no heart failure. All are apparent (divided by the unestimable
    # bioavailability F), because no intravenous data were available
    # (Sunnaker 2026 Results, Base Model).
    lcl <- log(22.1); label("Apparent clearance (L/h)")                          # Sunnaker 2026 Table 3, final model: CL/F = 22.1 L/h (RSE 3.3)
    lvc <- log(742);  label("Apparent central volume of distribution (L)")       # Sunnaker 2026 Table 3, final model: Vc/F = 742 L (RSE 6.4)
    lq  <- log(67.1); label("Apparent intercompartmental clearance (L/h)")       # Sunnaker 2026 Table 3, final model: Q/F = 67.1 L/h (RSE 4.8)
    lvp <- log(834);  label("Apparent peripheral volume of distribution (L)")    # Sunnaker 2026 Table 3, final model: Vp/F = 834 L (RSE 3.6)
    lka <- log(1.94); label("First-order absorption rate constant (1/h)")        # Sunnaker 2026 Table 3, final model: Ka = 1.94 1/h (RSE 10)

    # Covariate effects. The two functional forms are given in Sunnaker 2026
    # Methods, Covariate Model:
    #   continuous  theta_x,i = theta_x * (C_i / C_median)^beta_x
    #   categorical theta_x,i = theta_x * (1 + theta_xcov * COV_i)
    # Centering values are in the Table 3 footnote.
    e_crcl_cl       <-  0.45; label("Power exponent on baseline eGFR relative to its population median for apparent clearance (unitless)")                                                                  # Sunnaker 2026 Table 3, final model: effect of eGFR on CL/F = 0.45 (RSE 12); median 99 mL/min/1.73 m^2
    e_wt_cl         <-  0.78; label("Power exponent on baseline body weight relative to its population median for apparent clearance (unitless)")                                                           # Sunnaker 2026 Table 3, final model: effect of baseline body weight on CL/F = 0.78 (RSE 19); median 77.95 kg
    e_age_vc        <-  0.54; label("Power exponent on age relative to its population median for apparent central volume of distribution (unitless)")                                                       # Sunnaker 2026 Table 3, final model: effect of age on Vc/F = 0.54 (RSE 29); median 41 years
    e_race_asian_cl <-  0.27; label("Fractional change in apparent clearance for Asian relative to non-Asian participants (unitless)")                                                                      # Sunnaker 2026 Table 3, final model: effect of race on CL/F = 0.27 (RSE 31)
    e_hfpef_cl      <- -0.23; label("Fractional change in apparent clearance for patients with heart failure with preserved or mildly reduced ejection fraction relative to healthy volunteers (unitless)")  # Sunnaker 2026 Table 3, final model: effect of disease status on CL/F = -0.23 (RSE 19)

    # Interindividual variability. Sunnaker 2026 Methods, Structural Base Model,
    # defines theta_i = theta * exp(eta_i) with eta_i ~ N(0, omega^2), and the
    # Table 3 note states that the tabulated IIV CV percentages were computed as
    # the square root of the variance, so omega = CV/100 and omega^2 = (CV/100)^2.
    # No interindividual variability was estimated for Q/F.
    etalcl ~ 0.0484  # Sunnaker 2026 Table 3, final model: IIV CV for CL/F = 22% (RSE 8.0); 0.22^2
    etalvc ~ 0.25    # Sunnaker 2026 Table 3, final model: IIV CV for Vc/F = 50% (RSE 8.6); 0.50^2
    etalvp ~ 0.0441  # Sunnaker 2026 Table 3, final model: IIV CV for Vp/F = 21% (RSE 13); 0.21^2
    etalka ~ 0.8836  # Sunnaker 2026 Table 3, final model: IIV CV for Ka = 94% (RSE 9.1); 0.94^2

    # Residual error. Sunnaker 2026 Methods, Structural Base Model:
    # log(y_ij) = log(yhat_ij) + e with e ~ N(0, sigma^2), i.e. an additive
    # error on the log scale, which is a log-normal residual error on the
    # linear scale.
    expSd <- 0.22; label("Log-scale additive residual error (unitless)")  # Sunnaker 2026 Table 3, final model: residual error sigma in log-space = 0.22 (RSE 0.5)
  })

  model({
    # Sunnaker 2026 Table 3 footnote and Methods, Covariate Model:
    #   CL/F = 22.1 * (eGFR/99)^0.45 * (WT/77.95)^0.78
    #               * (1 + 0.27 * Asian) * (1 - 0.23 * HFpEF)
    #   Vc/F = 742  * (AGE/41)^0.54
    #   Q/F  = 67.1                      (no covariates, no IIV)
    #   Vp/F = 834                       (no covariates)
    crcl_cl <- (CRCL / 99)^e_crcl_cl
    wt_cl   <- (WT / 77.95)^e_wt_cl
    race_cl <- 1 + e_race_asian_cl * RACE_ASIAN
    dis_cl  <- 1 + e_hfpef_cl * DIS_HFPEF
    age_vc  <- (AGE / 41)^e_age_vc

    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * crcl_cl * wt_cl * race_cl * dis_cl
    vc <- exp(lvc + etalvc) * age_vc
    q  <- exp(lq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
