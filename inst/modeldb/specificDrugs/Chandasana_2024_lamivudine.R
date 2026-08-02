Chandasana_2024_lamivudine <- function() {
  description <- "Two-compartment oral population PK model with first-order absorption and elimination for lamivudine 300 mg once daily as the dolutegravir/lamivudine fixed-dose combination in virologically suppressed adults living with HIV-1 (TANGO; Chandasana 2024)"
  reference <- "Chandasana H, Singh R, Adkison K, Ait-Khaled M, Pene Dumitrescu T. Population pharmacokinetic modeling of dolutegravir/lamivudine to support a once-daily fixed-dose combination regimen in virologically suppressed adults living with HIV-1. Antimicrob Agents Chemother. 2024;68(5):e01504-23. doi:10.1128/aac.01504-23"
  vignette <- "Chandasana_2024_dolutegravir_lamivudine"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline weight; fixed allometric scaling with exponent 0.75 on CL/F and Q/F, exponent 1.0 on V2/F (=V3/F), reference 70 kg. Range 50.2-153.0 kg in the analysis population (Chandasana 2024 Table 1).",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate (CKD-EPI)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "eGFR calculated with the CKD-EPI equation and BSA-normalised to mL/min/1.73 m^2 (Chandasana 2024 Table 1 footnote). Power effect on CL/F, reference 99 mL/min/1.73 m^2 (population median per Chandasana 2024 Table 1). Range 44.0-147.0. Canonical column name CRCL covers the BSA-normalised CKD-EPI eGFR variant used here.",
      source_name        = "eGFR"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Black / African American; pools White, Asian, American Indian or Alaskan, Native Hawaiian or other Pacific, and Other races)",
      notes              = "Exponential effect on CL/F encoded as (theta)^RACE in Chandasana 2024 Table 2 footnote (0 = non-Black/African American, 1 = Black/African American). The paper reports 21% lower CL/F in Black/African American subjects; the fitted multiplier theta = 0.789 reproduces this reduction.",
      source_name        = "RACE"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline",
      units       = "years",
      type        = "continuous",
      notes       = "Screened but not retained in the final lamivudine model (Chandasana 2024 Discussion)."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not significant for lamivudine (Chandasana 2024 Discussion)."
    ),
    RACE_HISPANIC = list(
      description = "Hispanic / Latino ethnicity indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Ethnicity was screened as a lamivudine covariate but not retained; only race (Black/African American vs other, encoded via RACE_BLACK) reached significance (Chandasana 2024 Discussion)."
    ),
    SCR = list(
      description = "Serum creatinine at baseline",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened but not retained separately once eGFR was included (Chandasana 2024 Discussion)."
    ),
    CDC_HIV = list(
      description = "CDC HIV classification (1 vs 2/3)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened but not significant (Chandasana 2024 Discussion)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 362L,
    n_studies      = 1L,
    n_observations = 2611L,
    age_range      = "20-74 years (median 40)",
    weight_range   = "50.2-153.0 kg (median 78.8)",
    sex_female_pct = 6.9,
    race_ethnicity = c(White = 80.4, `Black or African American` = 13.5, Asian = 3.59,
                       `American Indian or Alaskan` = 1.93,
                       `Native Hawaiian or other Pacific` = 0.28, Other = 0.28),
    ethnicity      = c(`Non-Hispanic or Latino` = 80.7, `Hispanic or Latino` = 19.3),
    disease_state  = "Virologically-suppressed adults living with HIV-1 switching to dolutegravir 50 mg / lamivudine 300 mg fixed-dose combination once daily (TANGO phase 3 PK substudy).",
    dose_range     = "Lamivudine 300 mg once daily as the dolutegravir/lamivudine 50/300 mg fixed-dose combination tablet.",
    regions        = "Multinational phase 3 study (NCT03446573).",
    renal_function = "eGFR (mL/min/1.73 m^2, CKD-EPI) 99 median, range 44.0-147.0 (Chandasana 2024 Table 1). Distribution: normal (>90) 82.3%, mild (60-89) 17.4%, moderate (50-59) 0.28%.",
    notes          = "TANGO PK substudy: sparse sampling in all subjects at weeks 4, 8, 12, 24, 36, 48 plus intensive sampling in 30 subjects at week 4 (pre-dose, 0.5, 1, 1.5, 2, 3, 4, 6, 10, 24 h post-dose). 362 subjects contributed 2,611 lamivudine samples for the population PK analysis. A base one-compartment lamivudine model (Moore 1999-style twice-daily backbone) was extended to two compartments here because the once-daily sampling captured absorption, distribution, and elimination phases (Chandasana 2024 Discussion); V3/F was set equal to V2/F because V3/F alone was unreliably estimated on the full data set."
  )

  ini({
    # Structural parameters; reference covariates:
    #   WT = 70 kg, CRCL = 99 mL/min/1.73 m^2, RACE_BLACK = 0.
    # From Chandasana 2024 Table 2 final lamivudine model.
    lka  <- log(2.30); label("Absorption rate constant Ka (1/h)")                              # Chandasana 2024 Table 2: Ka = 2.30 h^-1 (95% CI 1.6, 2.99)
    lcl  <- log(19.6); label("Apparent oral clearance CL/F (L/h)")                             # Chandasana 2024 Table 2: CL/F = 19.6 L/h (95% CI 18.8, 20.5)
    lvc  <- log(105);  label("Apparent central volume V2/F (L; V3/F = V2/F)")                  # Chandasana 2024 Table 2: V2/F = V3/F = 105 L (95% CI 96.9, 113)
    lq   <- log(2.97); label("Apparent intercompartmental clearance Q/F (L/h)")                # Chandasana 2024 Table 2: Q/F = 2.97 L/h (95% CI 2.59, 3.34)

    # Allometric weight scaling (Chandasana 2024 Table 2 footnote; fixed
    # exponents per Methods and Results). CL/F and Q/F share the WT^0.75
    # exponent; V2/F (and V3/F, which is set equal to V2/F) uses WT^1.0.
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of WT shared across CL/F and Q/F (unitless)")  # Chandasana 2024 Results / Table 2 footnote
    e_wt_vc   <- fixed(1.0);  label("Allometric exponent of WT on V2/F = V3/F (unitless)")              # Chandasana 2024 Results / Table 2 footnote

    # Estimated covariate effects on CL/F.
    #   CL/F: (CRCL/99)^e_crcl_cl * (theta_race)^RACE_BLACK
    # In-file encoding: (CRCL/99)^e_crcl_cl * exp(e_race_black_cl * RACE_BLACK)
    # with e_race_black_cl <- log(0.789) reproducing the paper's theta.
    e_crcl_cl       <-  0.533;      label("Power exponent of eGFR on CL/F (unitless)")                                                # Chandasana 2024 Table 2: CL/F~eGFR = 0.533 (95% CI 0.336, 0.731)
    e_race_black_cl <- log(0.789);  label("log-scale multiplier of Black/African American race on CL/F (unitless; exp(e_race_black_cl) = 0.789)") # Chandasana 2024 Table 2: CL/F~race = 0.789 (95% CI 0.691, 0.887)

    # IIV -- omega^2 on the log scale per Chandasana 2024 Table 2.
    # A between-subject covariance (BSC) was reported between eta_CL and eta_V2:
    #   omega^2 CL   = 0.0883 (CV 30.4%)
    #   omega^2 V2   = 0.158  (CV 41.4%)
    #   BSC(CL, V2)  = -0.0531
    # The residual model additionally carries an IIV term on the proportional
    # error magnitude itself: Y = IPRED * [1 + theta_prop * exp(eta5) * eps1],
    # with sigma^2(eps1) fixed to 1. Encoded in model() as
    # propSd_i = propSd * exp(etapropSd) so etapropSd corresponds to eta5.
    etalcl + etalvc ~ c(0.0883,
                        -0.0531, 0.158)  # Chandasana 2024 Table 2 (omega^2 CL, BSC CL~V2, omega^2 V2)
    etapropSd ~ 0.247                    # Chandasana 2024 Table 2: omega^2 prop RUV = 0.247 (CV 52.9%)

    # Residual error -- proportional only, with IIV on its magnitude
    # (see etapropSd above). Table 2 reports sigma^2_prop = 0.359;
    # propSd = sqrt(0.359) = 0.5992.
    propSd <- sqrt(0.359); label("Proportional residual error scale (fraction; typical CV ~60% before IIV on RUV)")  # Chandasana 2024 Table 2: sigma^2 prop = 0.359
  })

  model({
    # Individual PK parameters. Reference covariates: WT = 70 kg,
    # CRCL = 99 mL/min/1.73 m^2, RACE_BLACK = 0. Fixed allometric
    # exponents (0.75 for CL and Q, 1.0 for V) are applied inline.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (WT   / 70)^e_wt_cl_q *
      (CRCL / 99)^e_crcl_cl *
      exp(e_race_black_cl * RACE_BLACK)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- vc
    q  <- exp(lq) * (WT / 70)^e_wt_cl_q

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Concentration: dose in mg, V in L => central / vc is in mg/L = ug/mL
    # (matches the paper's plasma concentration units in Table 3).
    Cc <- central / vc

    # Per-subject proportional residual error carries the paper's IIV on
    # the proportional error term (eta5 in the Chandasana 2024 Table 2
    # footnote; encoded here as etapropSd).
    propSd_i <- propSd * exp(etapropSd)
    Cc ~ prop(propSd_i)
  })
}
