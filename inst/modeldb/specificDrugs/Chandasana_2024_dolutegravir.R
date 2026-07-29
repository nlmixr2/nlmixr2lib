Chandasana_2024_dolutegravir <- function() {
  description <- "One-compartment oral population PK model with first-order absorption and elimination for dolutegravir 50 mg once daily as the dolutegravir/lamivudine fixed-dose combination in virologically suppressed adults living with HIV-1 (TANGO; Chandasana 2024)"
  reference <- "Chandasana H, Singh R, Adkison K, Ait-Khaled M, Pene Dumitrescu T. Population pharmacokinetic modeling of dolutegravir/lamivudine to support a once-daily fixed-dose combination regimen in virologically suppressed adults living with HIV-1. Antimicrob Agents Chemother. 2024;68(5):e01504-23. doi:10.1128/aac.01504-23"
  vignette <- "Chandasana_2024_dolutegravir_lamivudine"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline weight; power effect on CL/F and V/F, reference 79 kg (population median per Chandasana 2024 Table 1). Range 50.2-153.0 kg in the analysis population.",
      source_name        = "WT"
    ),
    TBILI = list(
      description        = "Baseline total bilirubin",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL/F, reference 8 umol/L (population median per Chandasana 2024 Table 1). Range 2.00-34.00 umol/L. Direction (negative exponent) is consistent with UGT1A1-mediated competition (higher endogenous bilirubin lowers dolutegravir CL).",
      source_name        = "BILI"
    ),
    RACE_HISPANIC = list(
      description        = "Hispanic / Latino ethnicity indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Hispanic / Latino)",
      notes              = "Exponential effect on CL/F encoded as (theta)^ETHN in Chandasana 2024 Table 2 footnote (0 = non-Hispanic/Latino, 1 = Hispanic/Latino). The paper reports ~15.6% lower CL/F in Hispanic/Latino subjects; the fitted multiplier theta = 0.844 reproduces this reduction. In the U.S. OMB scheme Hispanic is an ethnicity rather than a race, but the paper treats it as a race-like binary covariate on CL/F; encoded here with the canonical RACE_HISPANIC column.",
      source_name        = "ETHN"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline",
      units       = "years",
      type        = "continuous",
      notes       = "Screened but not retained in the final dolutegravir model (Chandasana 2024 Results and Discussion)."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not significant for dolutegravir (Chandasana 2024 Results); differs from Zhang 2015 which reported a female-sex effect on bioavailability."
    ),
    SMOKE = list(
      description = "Smoking status (current vs former/never)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not significant for dolutegravir in TANGO (Chandasana 2024 Discussion); Zhang 2015 previously reported a smoking effect."
    ),
    ALB = list(
      description = "Serum albumin at baseline",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not significant (Chandasana 2024 Discussion)."
    ),
    SCR = list(
      description = "Serum creatinine at baseline",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened but not significant (Chandasana 2024 Discussion)."
    ),
    CRCL = list(
      description = "eGFR (CKD-EPI) at baseline",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Screened but not significant for dolutegravir (Chandasana 2024 Discussion)."
    ),
    CDC_HIV = list(
      description = "CDC HIV classification (1 vs 2/3)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened but not significant (Chandasana 2024 Discussion)."
    ),
    RACE = list(
      description = "Race (White / Black or African American / Asian / American Indian or Alaskan / Native Hawaiian or Other Pacific / Other)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened but not significant for dolutegravir (Chandasana 2024 Discussion)."
    ),
    CONMED_METAL = list(
      description = "Metal-cation-containing product co-administration",
      units       = "(binary)",
      type        = "binary",
      notes       = "Concomitant medication class evaluated (used by 11-13% of subjects); not retained as a significant covariate (Chandasana 2024 Discussion)."
    ),
    CONMED_CYP3A4_INHIBITOR = list(
      description = "CYP3A inhibitor co-administration (weak / moderate / strong pooled)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Concomitant medication class evaluated (used by 17-19% of subjects); not retained as a significant covariate (Chandasana 2024 Discussion)."
    ),
    CONMED_UGT_INHIBITOR = list(
      description = "UGT inhibitor co-administration",
      units       = "(binary)",
      type        = "binary",
      notes       = "Concomitant medication class evaluated (Table S1); not retained as a significant covariate."
    ),
    HCV = list(
      description = "Hepatitis C co-infection at baseline",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (n = 16, 4.4%); not retained as a significant covariate on dolutegravir (Chandasana 2024 Discussion)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 362L,
    n_studies      = 1L,
    n_observations = 2629L,
    age_range      = "20-74 years (median 40)",
    weight_range   = "50.2-153.0 kg (median 78.8)",
    sex_female_pct = 6.9,
    race_ethnicity = c(White = 80.4, `Black or African American` = 13.5, Asian = 3.59,
                       `American Indian or Alaskan` = 1.93,
                       `Native Hawaiian or other Pacific` = 0.28, Other = 0.28),
    ethnicity      = c(`Non-Hispanic or Latino` = 80.7, `Hispanic or Latino` = 19.3),
    disease_state  = "Virologically-suppressed adults living with HIV-1 switching to dolutegravir 50 mg / lamivudine 300 mg fixed-dose combination once daily (TANGO phase 3 PK substudy).",
    dose_range     = "Dolutegravir 50 mg once daily as the dolutegravir/lamivudine 50/300 mg fixed-dose combination tablet.",
    regions        = "Multinational phase 3 study (NCT03446573).",
    smoking_pct    = "Never 49.7%, current 31.5%, former 18.8% (Chandasana 2024 Table 1)",
    hcv_pct        = "4.4% Hepatitis C co-infection at baseline (Chandasana 2024 Table 1)",
    notes          = "TANGO PK substudy: sparse sampling in all subjects at weeks 4, 8, 12, 24, 36, 48 plus intensive sampling in 30 subjects at week 4 (pre-dose, 0.5, 1, 1.5, 2, 3, 4, 6, 10, 24 h post-dose). 362 subjects contributed 2,629 dolutegravir samples for the population PK analysis (one additional subject was excluded from the final analysis owing to an incorrect race code; Chandasana 2024 Table 1 footnote b)."
  )

  ini({
    # Structural parameters; reference covariates:
    #   WT = 79 kg, TBILI = 8 umol/L, RACE_HISPANIC = 0 (non-Hispanic/Latino).
    # From Chandasana 2024 Table 2 final dolutegravir model.
    lka  <- log(2.15);  label("Absorption rate constant Ka (1/h)")                             # Chandasana 2024 Table 2: Ka = 2.15 h^-1 (95% CI 1.72, 2.59)
    lcl  <- log(0.858); label("Apparent oral clearance CL/F (L/h)")                             # Chandasana 2024 Table 2: CL/F = 0.858 L/h (95% CI 0.826, 0.890)
    lvc  <- log(16.7);  label("Apparent volume of distribution V/F (L)")                        # Chandasana 2024 Table 2: V/F = 16.7 L (95% CI 15.8, 17.6)

    # Continuous covariate power exponents (Chandasana 2024 Table 2 footnote).
    #   CL/F: (WT/79)^e_wt_cl * (TBILI/8)^e_tbili_cl
    #   V/F:  (WT/79)^e_wt_vc
    e_wt_cl    <-  0.427; label("Power exponent of body weight on CL/F (unitless)")            # Chandasana 2024 Table 2: CL/F~WT = 0.427 (95% CI 0.293, 0.562)
    e_wt_vc    <-  0.917; label("Power exponent of body weight on V/F (unitless)")             # Chandasana 2024 Table 2: V/F~WT = 0.917 (95% CI 0.711, 1.12)
    e_tbili_cl <- -0.153; label("Power exponent of total bilirubin on CL/F (unitless)")        # Chandasana 2024 Table 2: CL/F~bilirubin = -0.153 (95% CI -0.223, -0.0827)

    # Ethnicity effect on CL/F -- encoded in Chandasana 2024 Table 2 footnote as
    # (CL/F~ETHN)^ETHN with theta = 0.844 (95% CI 0.777, 0.911), i.e. a
    # multiplicative fold-change of 0.844 for Hispanic/Latino vs the
    # non-Hispanic/Latino reference. Encoded here as log(0.844) so the
    # in-model expression uses exp(e_race_hispanic_cl * RACE_HISPANIC),
    # matching the canonical (theta^indicator) parameterisation.
    e_race_hispanic_cl <- log(0.844); label("log-scale multiplier of Hispanic/Latino ethnicity on CL/F (unitless; exp(e_race_hispanic_cl) = 0.844)") # Chandasana 2024 Table 2: CL/F~ethnicity = 0.844 (95% CI 0.777, 0.911)

    # IIV -- omega^2 on the log scale per Chandasana 2024 Table 2.
    # The paper's residual model carries an IIV term on the proportional
    # error magnitude itself: Y = IPRED * [1 + theta_prop * exp(eta4) * eps1],
    # with sigma^2(eps1) fixed to 1. To reproduce this structure in nlmixr2
    # the per-subject residual SD is built inside model() as
    # propSd_i = propSd * exp(etapropSd), and Cc ~ prop(propSd_i) is used
    # for the residual, so etapropSd corresponds to the paper's eta4.
    etalcl   ~ 0.0682   # Chandasana 2024 Table 2: omega^2 CL/F      = 0.0682 (CV 26.6%)
    etapropSd ~ 0.0567  # Chandasana 2024 Table 2: omega^2 prop RUV  = 0.0567 (CV 24.2%)

    # Residual error -- proportional only, but with IIV on its magnitude
    # (see etapropSd above). Table 2 reports sigma^2_prop = 0.341;
    # propSd = sqrt(0.341) = 0.5840, i.e. the typical residual CV before
    # the per-subject exp(etapropSd) modifier.
    propSd <- sqrt(0.341); label("Proportional residual error scale (fraction; CV ~58% before IIV on RUV)")  # Chandasana 2024 Table 2: sigma^2 prop = 0.341
  })

  model({
    # Individual PK parameters. Reference covariates: WT = 79 kg,
    # TBILI = 8 umol/L, RACE_HISPANIC = 0 (non-Hispanic/Latino).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (WT    / 79)^e_wt_cl *
      (TBILI /  8)^e_tbili_cl *
      exp(e_race_hispanic_cl * RACE_HISPANIC)
    vc <- exp(lvc) * (WT / 79)^e_wt_vc

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Concentration: dose in mg, V in L => central / vc is in mg/L = ug/mL
    # (matches the paper's plasma concentration units in Table 3).
    Cc <- central / vc

    # Per-subject proportional residual error carries the paper's IIV on
    # the proportional error term (eta4 in the Chandasana 2024 Table 2
    # footnote; encoded here as etapropSd).
    propSd_i <- propSd * exp(etapropSd)
    Cc ~ prop(propSd_i)
  })
}
