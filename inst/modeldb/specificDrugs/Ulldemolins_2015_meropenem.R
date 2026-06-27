Ulldemolins_2015_meropenem <- function() {
  description <- "One-compartment IV population PK model for meropenem in 30 critically ill adults with septic shock and continuous renal replacement therapy (Ulldemolins 2015). Clearance is the sum of a constant CRRT-mediated baseline (3.68 L/h at zero residual diuresis) and an additive linear contribution from 24-hour residual diuresis (0.22 L/h per 100 mL/24h); central volume scales with body weight by power exponent 2.07 around the population-median 73 kg. CRRT intensity, blood flow, filter type, and serum albumin were tested but not retained."
  reference <- "Ulldemolins M, Soy D, Llaurado-Serra M, Vaquer S, Castro P, Rodriguez AH, Pontes C, Calvo G, Torres A, Martin-Loeches I. Meropenem population pharmacokinetics in critically ill patients with septic shock and continuous renal replacement therapy: influence of residual diuresis on dose requirements. Antimicrob Agents Chemother. 2015;59(9):5520-5528. doi:10.1128/AAC.00712-15"
  vignette <- "Ulldemolins_2015_meropenem"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight on admission",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Ulldemolins 2015 Table 1: median 72.8 kg (range 49-126) across the 30-subject cohort. Reference value 73 kg (population median rounded) used in the power-scaling term V = theta_V * (WT/73)^theta_WT (Results, p. 5523).",
      source_name        = "WT"
    ),
    URINE_VOL_24H = list(
      description        = "Residual diuresis: total urine volume collected over the 24 hours of the natural day of the PK study",
      units              = "mL/24h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Ulldemolins 2015 Methods (p. 5521): 'the level of residual diuresis (defined as the volume of urine collected over the 24 h of the natural day of the study)'. Time-fixed per study day; collected by nursing staff. Median 137.5 mL/24h (range <10-2050) overall. Clinical cutoffs used in the paper: anuria <100, oliguria 100-500, preserved diuresis >500 mL/24h. Source-paper alias 'residual diuresis' / theta_DIUR. Used in the additive linear effect CL = theta_CL + theta_DIUR * (URINE_VOL_24H / 100) centered at the 100 mL/24h anuria cutoff.",
      source_name        = "residual diuresis"
    )
  )

  covariatesDataExcluded <- list(
    CRRT_INTENSITY = list(
      description = "CRRT intensity (filtrate + dialysate flow rate) / ideal body weight",
      units       = "mL/kg/h",
      type        = "continuous",
      notes       = "Screened on CL; not retained. Ulldemolins 2015 Results (p. 5523): 'exploratory and regression analyses on the effects of covariates on individual CL did not show any visual or statistical trend between intensity and the estimates of individual CL'. Median 34.7 mL/kg/h (range 18.7-60.1) in the cohort. Consistent with Roberts et al. who also failed to identify intensity as a meropenem CL modifier (Discussion, p. 5526)."
    ),
    BFR = list(
      description = "Blood flow rate through the CRRT circuit",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened on CL; not retained. Ulldemolins 2015 Results (p. 5523). Median 200 mL/min (range 130-250) in the cohort."
    ),
    ALB = list(
      description = "Serum albumin concentration",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened on both CL and V; not retained. Ulldemolins 2015 Results (p. 5523). Median 23.4 g/L (range 12.4-38) in the cohort."
    ),
    AGE = list(
      description = "Age on admission",
      units       = "year",
      type        = "continuous",
      notes       = "Screened on V; not retained. Ulldemolins 2015 Results (p. 5523). Median 66.5 years (range 34-85)."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 30L,
    n_studies        = 1L,
    age_range        = "34-85 years",
    age_median       = "66.5 years",
    weight_range     = "49-126 kg",
    weight_median    = "72.8 kg",
    sex_female_pct   = 46.7,
    race_ethnicity   = "Not reported (multicenter Spanish ICU cohort)",
    disease_state    = "Critically ill adults with septic shock (Surviving Sepsis Campaign criteria) and continuous renal replacement therapy. CVVHDF in 26/30 patients, CVVHF in 4/30. Median APACHE II 24 (range 5-44), median SOFA 12 (range 4-19) on the day of study. Hepatic impairment (LFTs > 2 ULN) in 6/30. Vasopressors in 28/30. Sources of infection: intra-abdominal 13, respiratory 7, bloodstream 4, urinary tract 2, CNS 2, unknown 2.",
    dose_range       = "Meropenem IV 500 mg q6h-q12h or 1000 mg q8h-q12h or 2000 mg q8h. Administered as a 30-min bolus or a 3-h / 4-h extended infusion at the discretion of the treating physician. Median duration of therapy 10 days (range 4-28).",
    regions          = "Spain (three centers: Hospital Parc Tauli, Sabadell; Hospital Clinic, Barcelona; Hospital Joan XXIII, Tarragona). Enrolled January 2012-May 2014.",
    renal_function   = "All patients on CRRT. Median residual diuresis 137.5 mL/24h (range <10-2050); 14 anuric (<100), 11 oliguric (100-500), 5 with preserved diuresis (>500). Median serum creatinine 1.4 mg/dL (range 0.4-2.6). Median urea 61.7 mg/dL (range 22-168).",
    n_concentrations = 153L,
    notes            = "Demographics from Ulldemolins 2015 Table 1 (overall n=30 column). Model development used 24 subjects with 124 samples; external validation used 6 subjects with the remaining 29 samples. Total meropenem concentration measured by LC-MS/MS, linear over 0.4-300 mg/L. NONMEM 7.3 FOCE-INTER. Bootstrap 200 replicates (Table 3). Single-vs-multicenter design: multicenter ICU cohort, single protocol."
  )

  ini({
    # Structural parameters (Ulldemolins 2015 Table 3 final-model "Estimate"
    # column; the structural equations are stated verbatim in Results p. 5523:
    # CL = 3.68 + 0.22 * (residual diuresis / 100), L/h
    # V  = 33.00 * (weight / 73)^2.07,             L
    # Reference subject: URINE_VOL_24H = 0 mL/24h (anuric); WT = 73 kg.
    lcl <- log(3.68);  label("Clearance intercept at URINE_VOL_24H = 0 (L/h)")  # Ulldemolins 2015 Table 3: theta_CL = 3.68 (RSE 11%)
    lvc <- log(33.00); label("Central volume at WT = 73 kg (L)")                # Ulldemolins 2015 Table 3: theta_V  = 33.00 (RSE 10%)

    # Covariate effects (Ulldemolins 2015 Table 3 + Results p. 5523 equations):
    #   CL (L/h) = exp(lcl) + e_urine_vol_24h_cl * (URINE_VOL_24H / 100)  -- additive linear
    #   V  (L)   = exp(lvc) * (WT / 73)^e_wt_vc                            -- power
    e_urine_vol_24h_cl <- 0.22; label("Additive slope of CL on URINE_VOL_24H/100 (L/h per 100 mL/24h)")  # Ulldemolins 2015 Table 3: theta_DIUR = 0.22 (RSE 47%)
    e_wt_vc            <- 2.07; label("Power exponent on (WT/73) for V")                                 # Ulldemolins 2015 Table 3: theta_WT  = 2.07 (RSE 24%)

    # Inter-individual variability (Ulldemolins 2015 Table 3 "IIV_CL (% CV)"
    # and "IIV_V (% CV)" rows, log-normal model theta_i = theta_TV * exp(eta_i);
    # for log-normal etas omega^2 = log(CV^2 + 1)).
    etalcl ~ 0.12834  # log(0.37^2 + 1); 37% CV on CL (RSE 27%)
    etalvc ~ 0.18441  # log(0.45^2 + 1); 45% CV on V  (RSE 61%)

    # Combined additive + proportional residual error (Ulldemolins 2015 Table 3).
    # The published "Proportional residual error" entry is "-0.258 (10% RSE)";
    # the negative sign is a NONMEM display convention for the residual SD
    # (the variance estimated is sigma^2 > 0). Encoded here as the magnitude.
    addSd  <- 0.0002; label("Additive residual error (mg/L)")                                            # Ulldemolins 2015 Table 3: additive    = 0.0002 mg/L (RSE 42.76%)
    propSd <- 0.258;  label("Proportional residual error (fraction)")                                    # Ulldemolins 2015 Table 3: proportional = -0.258 magnitude (RSE 10%)
  })
  model({
    # Individual PK parameters.
    #
    # CL has an ADDITIVE LINEAR effect of residual diuresis (not multiplicative):
    # the typical value is exp(lcl) + e_urine_vol_24h_cl * (URINE_VOL_24H / 100)
    # in L/h, and IIV is log-normal multiplicative on that typical value
    # (Ulldemolins 2015 Methods p. 5522: "Continuous covariates were assessed
    # as a proportional or a power function"; Results p. 5523: "CL = 3.68 +
    # 0.22 * (residual diuresis / 100)").
    #
    # V scales by (WT / 73)^2.07 (Results p. 5523).
    cl_typical <- exp(lcl) + e_urine_vol_24h_cl * (URINE_VOL_24H / 100)
    cl <- cl_typical * exp(etalcl)
    vc <- exp(lvc + etalvc) * (WT / 73)^e_wt_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, vc in L -> central/vc has units mg/L (matches the paper's
    # observed concentrations in mg/L; LC-MS/MS method linear 0.4-300 mg/L).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
