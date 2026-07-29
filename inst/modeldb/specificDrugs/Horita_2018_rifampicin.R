Horita_2018_rifampicin <- function() {
  description <- "One-compartment population pharmacokinetic model with sequential zero-order then first-order absorption and first-order elimination for oral rifampin (rifampicin) in Ghanaian children with active tuberculosis (Horita 2018); allometric weight scaling on CL/F (fixed 0.75) and V/F (fixed 1.0) normalised to the cohort median 14.3 kg."
  reference <- "Horita Y, Alsultan A, Kwara A, Antwi S, Enimil A, Ortsin A, Dompreh A, Yang H, Wiesner L, Peloquin CA. Evaluation of the Adequacy of WHO Revised Dosages of the First-Line Antituberculosis Drugs in Children with Tuberculosis Using Population Pharmacokinetic Modeling and Simulations. Antimicrob Agents Chemother. 2018;62(9):e00008-18. doi:10.1128/AAC.00008-18"
  vignette <- "Horita_2018_rifampicin"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL/F and V/F with fixed exponents 0.75 and 1.0 respectively (Horita 2018 Table 2 / Results 'RIF' paragraph 1: 'The fixed exponents were 0.75 for CL/F and 1.0 for V/F'). Reference weight is the cohort median 14.3 kg (Table 1: median weight 14.3 kg, IQR 9.70 to 20.1 kg). The paper does not state the reference weight explicitly; the cohort median was inferred by back-computation of the typical-value CL/F against the published adult value of pyrazinamide (Alsultan 2017 5.06 L/h at 70 kg), which matches Horita's PZA typical value when normalised at 14.3 kg.",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 113L,
    n_studies      = 1L,
    age_range      = "3 months to 14 years (median 5.00 years, IQR 2.17 to 8.25)",
    age_median     = "5.00 years",
    weight_range   = "5-30 kg (median 14.3, IQR 9.70 to 20.1)",
    weight_median  = "14.3 kg",
    sex_female_pct = 44.2,
    hiv_positive_pct = 52.2,
    nat2_slow_pct  = 45.1,
    disease_state  = "Ghanaian children with active tuberculosis (HIV-positive and HIV-negative). 21.2% under 2 years of age.",
    dose_range     = "Rifampin 10-20 mg/kg orally daily (median 15.8 mg/kg, IQR 13.6-18.8). Administered as part of standard four-drug anti-TB regimen (RIF, INH, PZA, EMB).",
    regions        = "Ghana (Komfo Anokye Teaching Hospital, Kumasi).",
    notes          = "Patients enrolled October 2012-August 2015. PK sampling after at least 4 weeks of anti-TB treatment (steady state). Blood samples at 0, 1, 2, 4, 8 h postdose. RIF concentrations 0.117-30 ug/mL by LC-MS/MS. ClinicalTrials.gov NCT01687504. Demographics from Horita 2018 Table 1; structural model and parameters from Table 2."
  )

  ini({
    # Structural PK parameters -- Horita 2018 Table 2 final population pharmacokinetic
    # model for rifampin. Typical values are at the cohort median weight 14.3 kg (see
    # the WT covariate notes for the reference-weight derivation).
    logitfr_zo <- logit(0.0878); label("Fraction of dose absorbed via zero-order (proportion)")   # Table 2: Fr (proportion) = 0.0878 (RSE 34%) -- 'About 9% was explained by zero-order absorption'
    ltk0       <- log(0.342);    label("Duration of zero-order absorption phase, Tk0 (h)")        # Table 2: Tk0 = 0.342 h (RSE 19%)
    lka        <- log(0.645);    label("First-order absorption rate constant after zero-order phase, ka (1/h)")  # Table 2: ka = 0.645 1/h (RSE 14%)
    lvc        <- log(13.8);     label("Apparent central volume of distribution V/F at WT = 14.3 kg (L)")  # Table 2: V/F = 13.8 L (RSE 10%)
    lcl        <- log(7.53);     label("Apparent oral clearance CL/F at WT = 14.3 kg (L/h)")             # Table 2: CL/F = 7.53 L/h (RSE 5%)

    # Allometric exponents on body weight -- fixed at canonical theoretical values.
    # Horita 2018 Results 'RIF' paragraph 1: 'The fixed exponents were 0.75 for CL/F and 1.0 for V/F.'
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F (fixed, unitless)")   # Horita 2018 Results 'RIF' paragraph 1: fixed exponent
    e_wt_vc <- fixed(1.0);  label("Allometric exponent on V/F  (fixed, unitless)")   # Horita 2018 Results 'RIF' paragraph 1: fixed exponent

    # Inter-individual variability. Table 2 IIV column reports 'omega (CV%)' on the
    # log scale (back-transformed via CV% = sqrt(exp(omega^2) - 1) * 100). Encode
    # variance = omega^2.
    etalogitfr_zo ~ 1.166   # Table 2: 1.08 (148.7% CV)  -- 1.08^2 = 1.166 (logit-scale variance for the zero-order fraction)
    etaltk0       ~ 0.835   # Table 2: 0.914 (114.3% CV) -- 0.914^2 = 0.835
    etalka        ~ 0.215   # Table 2: 0.464 (49% CV)    -- 0.464^2 = 0.215
    etalvc        ~ 0.0471  # Table 2: 0.217 (22% CV)    -- 0.217^2 = 0.0471
    etalcl        ~ 0.299   # Table 2: 0.547 (59.1% CV)  -- 0.547^2 = 0.299

    # Combined residual error. Table 2: 'Constant a' = 0.0476 (RSE 14%), 'Slope b'
    # = 0.202 (RSE 7%). Table 4 footnote b: combined-1 form y = f + (a + bf)*eps;
    # the nlmixr2 convention `add(addSd) + prop(propSd)` uses the Pythagorean
    # combined-2 form. The difference is small in practice; see the vignette
    # Errata for the deviation.
    addSd  <- 0.0476; label("Additive residual SD (ug/mL)")                            # Table 2: constant a = 0.0476
    propSd <- 0.202;  label("Proportional residual SD (fraction)")                     # Table 2: slope b    = 0.202
  })

  model({
    # Individual PK parameters with allometric weight scaling (reference 14.3 kg).
    fr_zo <- expit(logitfr_zo + etalogitfr_zo)
    tk0   <- exp(ltk0 + etaltk0)
    ka    <- exp(lka  + etalka)
    cl    <- exp(lcl  + etalcl) * (WT / 14.3)^e_wt_cl
    vc    <- exp(lvc  + etalvc) * (WT / 14.3)^e_wt_vc

    kel <- cl / vc

    # Sequential zero-order then first-order absorption (Horita 2018 Results 'RIF'
    # paragraph 1: 'sequential zero- and first-order absorption'). The zero-order
    # arm absorbs the fraction `fr_zo` of the dose at constant rate
    # `fr_zo * podo(depot) / tk0` over the interval [0, tk0]; the first-order arm
    # absorbs the remaining (1 - fr_zo) fraction at rate `ka * depot` for t > tk0.
    # This follows the rxode2 idiom used in Cirincione_2017_exenatide.R.
    mtime(tk_switch) <- tk0

    kzero <- fr_zo * podo(depot) / tk0
    if (tad(depot) > tk0) kzero <- 0.0

    ka_eff <- ka
    if (tad(depot) <= tk0) ka_eff <- 0.0

    d/dt(depot)   <- -ka_eff * depot - kzero
    d/dt(central) <-  ka_eff * depot + kzero - kel * central

    # Concentration: dose mg / Vc L -> mg/L = ug/mL (matches the source paper).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
