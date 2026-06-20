Valade_2014_emtricitabine <- function() {
  description <- "Two-compartment oral population PK model for emtricitabine (FTC) in HIV-infected pregnant and non-pregnant women, with first-order absorption and elimination. Creatinine clearance (Cockcroft-Gault, raw mL/min) on apparent oral clearance via the power model CL/F = 22.3 * (CRCL/135)^0.33 captures the 18% CL/F increase observed during pregnancy as a manifestation of the pregnancy-associated 50% rise in estimated glomerular filtration rate; pregnancy itself, gestational age, age, weight, serum creatinine and co-medication were screened but not retained after CLcr inclusion (Valade 2014, BJCP)."
  reference <- paste(
    "Valade E, Treluyer JM, Dabis F, Arrive E, Pannier E, Benaboud S,",
    "Fauchet F, Bouazza N, Foissac F, Urien S, Hirt D. (2014).",
    "Modified renal function in pregnancy: impact on emtricitabine",
    "pharmacokinetics. Br J Clin Pharmacol 78(6):1378-1386.",
    "doi:10.1111/bcp.12457"
  )
  vignette <- "Valade_2014_emtricitabine"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CLcr. Computed by the Cockcroft-Gault equation in raw mL/min (NOT BSA-normalized to mL/min/1.73 m^2). Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization, with the per-model description recording the assay form). Reference value 135 mL/min (population median across all 179 women, Valade 2014 Table 1). Applied as a power-scaling covariate on CL/F: CL/F = 22.3 * (CRCL / 135)^0.33. The CLcr increase during pregnancy (median 160 vs 106 mL/min, non-pregnant) is the load-bearing driver of the 18% CL/F rise; CLcr therefore subsumed the effects of pregnancy, gestational age, age, weight, serum creatinine and BW/Scr on CL/F (Results paragraph after the final covariate equation).",
      source_name        = "CLcr"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Tested as a continuous covariate on CL/F per the equation CL/F = TVCL * (AGE / median_AGE)^beta (Valade 2014 Methods, equation (i) for continuous covariates). Not retained: subsumed by CRCL after CLcr inclusion in the final model (Valade 2014 Results, paragraph following Table 2)."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Tested as a continuous covariate on CL/F (Valade 2014 Methods). Not retained: subsumed by CRCL after CLcr inclusion (Results)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Tested as a continuous covariate on CL/F (Valade 2014 Methods). Not retained: subsumed by CRCL after CLcr inclusion (Results)."
    ),
    PREG = list(
      description = "Pregnancy indicator (1 = pregnant or parturient, 0 = non-pregnant)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested as a binary covariate on CL/F per the equation CL/F = TVCL * beta^PREG (Valade 2014 Methods, equation (ii) for binary covariates). Not retained: subsumed by CRCL after CLcr inclusion. The pregnancy-associated 50% rise in glomerular filtration rate (median CLcr 160 mL/min pregnant vs 106 mL/min non-pregnant) is the load-bearing physiologic driver; encoding the effect mechanistically via CRCL on CL/F is preferred to a pregnancy indicator (Results)."
    ),
    GA = list(
      description = "Gestational age (pregnant or parturient women only; 0 otherwise)",
      units       = "weeks",
      type        = "continuous",
      notes       = "Tested as continuous, categorical (by trimester), and Hill-saturation covariates on CL/F (Valade 2014 Methods, equation (iii) for influence of pregnancy). Not retained: subsumed by CRCL after CLcr inclusion (Results). The increase in glomerular filtration rate begins early after conception and remains constant throughout pregnancy (Discussion); GA-keyed effects are therefore well approximated by the binary pregnancy state, which is in turn well approximated by CRCL."
    ),
    CONMED_ARV = list(
      description = "Associated antiretroviral co-medication (binary indicators for TDF, PI, NNRTI, NRTI)",
      units       = "(binary set)",
      type        = "categorical",
      notes       = "Tested per the binary-covariate equation CL/F = TVCL * beta^CONMED (Valade 2014 Methods, equation (ii)). Not retained: none of the four co-medication classes improved the model after CRCL inclusion (Results)."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 179L,
    n_observations   = 457L,
    n_studies        = 2L,
    age_range        = "16-72 years",
    age_median       = "35 years",
    weight_range     = "37-130 kg",
    weight_median    = "69 kg",
    sex_female_pct   = 100,
    race_ethnicity   = "Not reported (combined Paris-region TDM cohort + multi-country TEmAA ANRS 12109 cohort: Ivory Coast, Cambodia, South Africa)",
    disease_state    = "HIV-infected women receiving FTC for HIV treatment or for prevention of mother-to-child transmission (PMTCT)",
    dose_range       = "Oral FTC: 200 mg chronic daily (TDM cohort, n = 148 women) or 400 mg single dose at the start of labour (TEmAA ANRS 12109 cohort, n = 31 women)",
    pregnancy_status = "103 non-pregnant + 83 pregnant women (1 first trimester, 26 second trimester, 67 third trimester including 48 on the day of delivery; 7 women appear in both pregnant and non-pregnant strata on different occasions)",
    creatinine_range = "Serum creatinine median 60 umol/L (range 27-183); Cockcroft-Gault CrCl median 135 mL/min (range 35-335); CrCl median 160 mL/min in pregnant women, 106 mL/min in non-pregnant women",
    gestational_age  = "Pregnant cohort: median 37 weeks (range 5-41); the first-trimester stratum has only 1 woman",
    regions          = "France (Cochin Hospital TDM cohort) + Ivory Coast, Cambodia, South Africa (TEmAA ANRS 12109 labour study)",
    co_medication    = "Non-pregnant: 94% TDF, 61% PI, 22% NNRTI. Pregnant: 97% TDF, 57% PI, 40% NNRTI, 51% NRTI.",
    notes            = "Baseline demographics from Valade 2014 Table 1. 457 total FTC plasma concentrations; 14 below the LOQ (0.01 mg/L) were set to half-LOQ (3% of records). Pooled across the Cochin TDM study and the TEmAA ANRS 12109 study at labour."
  )

  ini({
    # Structural PK parameters (Valade 2014 Table 2 final-model estimates).
    # Two-compartment oral model with first-order absorption and elimination.
    # All PK parameters are apparent (.../F): bioavailability F is folded into
    # CL/F, Vc/F, Q/F and Vp/F since FTC was administered only orally and F
    # is not separately identifiable. Times in hours, doses in mg, central
    # volume in L, so the implied concentration unit is mg/L (matching the
    # paper's units throughout).
    lka <- log(0.616); label("Absorption rate constant ka (1/h)")                   # Valade 2014 Table 2: ka = 0.616 1/h (RSE 32%)
    lcl <- log(22.3);  label("Apparent oral clearance CL/F (L/h) at median CRCL = 135 mL/min") # Valade 2014 Table 2: CL/F = 22.3 L/h (RSE 3%)
    lvc <- log(100);   label("Apparent central volume of distribution Vc/F (L)")     # Valade 2014 Table 2: Vc/F = 100 L (RSE 16%)
    lq  <- log(5.89);  label("Apparent inter-compartmental clearance Q/F (L/h)")     # Valade 2014 Table 2: Q/F = 5.89 L/h (RSE 22%)
    lvp <- log(76.1);  label("Apparent peripheral volume of distribution Vp/F (L)")  # Valade 2014 Table 2: Vp/F = 76.1 L (RSE 28%)

    # Covariate effect: power-scaling of CL/F on CRCL with reference 135 mL/min
    # (population median across all 179 women, Valade 2014 Table 1). The final
    # covariate equation in the paper is CL/F = 22.3 * (CLcr / 135)^0.33
    # (Valade 2014 Results, formula immediately following Table 2 description).
    e_crcl_cl <- 0.33; label("Power-scaling exponent of CRCL on CL/F (unitless)")    # Valade 2014 Table 2: beta_CLcr on CL/F = 0.33 (RSE 16%)

    # Inter-individual variability. The paper uses Monolix 4.1.3 with an
    # exponential IIV model (P_i = TVP * exp(eta_i), eta_i ~ N(0, omega^2))
    # and reports the omega = SD of eta on the log scale (Valade 2014 Methods,
    # "Inter-individual variabilities (IIV or eta) were assumed to be
    # exponential"). Convert to nlmixr2 variance via variance = omega^2.
    # IIV was significant only for ka and CL/F (Valade 2014 Results).
    etalka ~ 0.253009  # Valade 2014 Table 2: omega_ka  = 0.503 -> variance = 0.503^2 = 0.253009; RSE 36%
    etalcl ~ 0.022801  # Valade 2014 Table 2: omega_CL/F = 0.151 -> variance = 0.151^2 = 0.022801; RSE 14%

    # Residual error. The paper reports a proportional error model with a
    # separate sigma for each of the two contributing studies (Valade 2014
    # Results: "Estimating a residual variability for each study resulted in
    # a 12.4 point decrease in the OFV and improved the goodness of fit").
    # The Monolix proportional model y_obs = y_pred * (1 + eps), eps ~ N(0,
    # sigma^2), maps directly to nlmixr2 `prop(propSd)` with propSd = sigma.
    # Default propSd is the TDM-study (200 mg chronic) sigma_1 = 0.505,
    # which characterises 148/179 (83%) of the cohort and the standard
    # chronic dosing scenario. The TEmAA labour-study (400 mg single)
    # sigma_2 = 0.373 is documented in the vignette's Assumptions and
    # deviations and can be substituted with `mod$propSd <- 0.373` if
    # simulating labour-day exposures.
    propSd <- 0.505; label("Proportional residual error on Cc (fraction; TDM study, 200 mg chronic dose)") # Valade 2014 Table 2: sigma_1 = 0.505 (RSE 6%) for the TDM study; sigma_2 = 0.373 (RSE 6%) for the TEmAA labour study (Assumptions and deviations in the vignette).
  })
  model({
    # Individual PK parameters
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (CRCL / 135)^e_crcl_cl
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # Micro-rate constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment oral PK with first-order absorption from a depot.
    d/dt(depot)       <- -ka  * depot
    d/dt(central)     <-  ka  * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central                                - k21 * peripheral1

    # Plasma FTC concentration (dose mg / Vc L gives mg/L, matching the
    # paper's units throughout Tables 2, 3 and 4 and Figures 1-3).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
