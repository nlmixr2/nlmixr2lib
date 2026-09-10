He_2025_YPEGrhGH <- function() {
  description <- paste(
    "Population PK/PD model for YPEG-rhGH (Y-shape branched 40 kDa",
    "PEGylated recombinant human growth hormone, Xiamen Amoytop) in",
    "healthy elderly subjects and healthy young adults. The PK is a",
    "two-compartment model with first-order subcutaneous absorption",
    "and Michaelis-Menten elimination from the central compartment,",
    "with all disposition parameters apparent (V1/F, V2/F, CL2/F)",
    "because bioavailability was not identifiable. Age (centred at",
    "35.5 years) is a power covariate on Ka, V1/F and Vmax, and body",
    "weight (centred at 70 kg) is a power covariate on Vmax. The",
    "pharmacodynamic layer is an indirect-response model in which",
    "YPEG-rhGH stimulates production of insulin-like growth factor 1",
    "(IGF-1) with a Hill coefficient fixed to 1; the PD state is the",
    "IGF-1 concentration expressed as a ratio to each subject's own",
    "baseline, so its typical baseline is Kin/Kout = 1. The PD",
    "parameters are the elderly-subject estimates."
  )
  reference <- paste(
    "He Y, Hu J, Zeng X, Yang Q, You Q, Huang J, Zhang Y, Si L,",
    "Zhai X. Population pharmacokinetics/pharmacodynamics and safety",
    "of YPEG-rhGH in elderly subjects.",
    "Front Pharmacol. 2025 Nov 25;16:1651323.",
    "doi:10.3389/fphar.2025.1651323.",
    "Final-model equations (including the covariate centring values)",
    "are given in the Supplementary Material, section",
    "'PopPK model equations' / 'PopPK/PD model equations'."
  )
  vignette <- "He_2025_YPEGrhGH"
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "ng/mL"
  )

  covariateData <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power covariate centred at 35.5 years (Supplementary",
        "Material, PopPK model equations). Enters three parameters:",
        "Ka with exponent -0.6957, V1/F with exponent 1.470 and Vmax",
        "with exponent 0.7989. The centring value 35.5 years is the",
        "median age of the pooled 52-subject analysis dataset (16",
        "elderly subjects aged 65-74 y plus 36 healthy young adult",
        "males aged 21-44 y, Table 2); it is not tabulated in the",
        "main text and is read from the supplement's final-model",
        "equations."
      ),
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline body weight, power covariate on Vmax with exponent",
        "2.202 centred at 70 kg (Supplementary Material, PopPK model",
        "equations). The pooled cohort medians were 59.0 kg (elderly)",
        "and 66.4 kg (healthy adults) per Table 2, so 70 kg is a",
        "rounded reference rather than the dataset median."
      ),
      source_name        = "WEIGHT"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "YPEG-rhGH", units = "ug",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "YPEG-rhGH", units = "ug",
      specimen = "serum", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "YPEG-rhGH", units = "ug",
      specimen = "serum", verified = TRUE
    ),
    igf1 = list(
      analyte = "insulin-like growth factor 1",
      units = "fraction of baseline",
      specimen = "serum", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 52L,
    n_studies      = 2L,
    age_range      = "65-74 years (elderly, n = 16); 21-44 years (healthy adults, n = 36)",
    age_median     = "66.0 years (elderly); 29.0 years (healthy adults)",
    weight_range   = "51.1-76.1 kg (elderly); 52.4-79.6 kg (healthy adults)",
    weight_median  = "59.0 kg (elderly); 66.4 kg (healthy adults)",
    sex_female_pct = 23,
    disease_state  = paste(
      "Healthy volunteers. The elderly cohort (study TB2208GH) was",
      "screened to a baseline IGF-1 standard-deviation score below 0",
      "and a BMI of 18-30 kg/m2, so it models the growth-hormone",
      "-deficient elderly target population without carrying a",
      "diagnosis of adult growth hormone deficiency."
    ),
    dose_range     = paste(
      "Elderly (TB2208GH): 30 ug/kg SC once every 2 weeks for 23",
      "weeks (12 injections). Healthy adults (TB1010GH): single SC",
      "doses of 10, 30, 60, 120 or 200 ug/kg."
    ),
    regions        = "China",
    notes          = paste(
      "Baseline demographics are Table 2 of the source paper; the",
      "study designs and sampling schedules are Table 1 and",
      "Supplementary Table 3. 813 of 884 PK samples (71 BQL samples",
      "excluded) and 980 IGF-1 measurements entered the analysis. The",
      "PK model was fitted to the pooled elderly + healthy-adult data;",
      "the indirect-response PD parameters in this file are the",
      "elderly-subject estimates (Table 4), fitted sequentially with",
      "individual PK parameters taken from the final PK model.",
      "Sex was NOT retained as a covariate: the two studies are",
      "confounded with sex (all 36 adults were male, 12 of 16 elderly",
      "subjects were female), so the authors kept AGE and WEIGHT",
      "instead (Discussion, limitations)."
    )
  )

  ini({
    # PK structural parameters - Table 3, reproduced with their
    # covariate terms in the Supplementary Material 'PopPK model
    # equations'. V1/F, V2/F and CL2/F are apparent parameters:
    # bioavailability was not estimated, so no f(depot) is applied.
    lka   <- log(0.01086); label("First-order SC absorption rate Ka (1/h) at AGE = 35.5 years")   # Table 3 Ka = 0.01086 (RSE 8.49%)
    lvc   <- log(2.387);   label("Apparent central volume V1/F (L) at AGE = 35.5 years")          # Table 3 V1/F = 2.387 (RSE 15.54%)
    lvp   <- log(21.42);   label("Apparent peripheral volume V2/F (L)")                           # Table 3 V2/F = 21.42 (RSE 20.48%)
    lq    <- log(0.05575); label("Apparent intercompartmental clearance CL2/F (L/h)")             # Table 3 CL2/F = 0.05575 (RSE 10.15%)
    lvmax <- log(80.13);   label("Michaelis-Menten maximum elimination rate Vmax (ug/h) at AGE = 35.5 years and WT = 70 kg") # Table 3 Vmax = 80.13 (RSE 7.82%)
    lkm   <- log(70.08);   label("Michaelis-Menten constant Km (ug/L = ng/mL)")                   # Table 3 Km = 70.08 (RSE 13.78%)

    # Covariate effects - power model theta_i = theta_TV *
    # (Cov_i / Cov_median)^theta_x (Equation 1), with the centring
    # values taken from the supplement's final-model equations
    # (AGE / 35.5 and WEIGHT / 70).
    e_age_ka   <- -0.6957; label("Power exponent of AGE/35.5 on Ka (unitless)")     # Table 3 'AGE on Ka' = -0.6957 (RSE 28.49%)
    e_age_vc   <- 1.470;   label("Power exponent of AGE/35.5 on V1/F (unitless)")   # Table 3 'AGE on V1/F' = 1.470 (RSE 24.69%)
    e_age_vmax <- 0.7989;  label("Power exponent of AGE/35.5 on Vmax (unitless)")   # Table 3 'AGE on Vmax' = 0.7989 (RSE 18.37%)
    e_wt_vmax  <- 2.202;   label("Power exponent of WT/70 on Vmax (unitless)")      # Table 3 'WEIGHT on Vmax' = 2.202 (RSE 18.05%)

    # PD structural parameters - indirect response on the IGF-1
    # ratio to baseline, elderly subjects (Table 4). The final PD
    # ODE is printed in the supplement:
    # dR/dt = 0.023 * (1 + 2.245 * C / (50.74 + C)) - 0.023 * R.
    lkin  <- log(0.023);  label("IGF-1 zero-order production rate Kin (baseline-ratio units per h)") # Table 4 Kin = 0.023 (RSE 14.79%)
    lkout <- log(0.023);  label("IGF-1 first-order degradation rate Kout (1/h)")                     # Table 4 Kout = 0.023 (RSE 15.50%)
    lemax <- log(2.245);  label("Maximum fractional stimulation of IGF-1 production Emax (unitless)") # Table 4 Emax = 2.245 (RSE 21.80%)
    lec50 <- log(50.74);  label("YPEG-rhGH concentration giving half-maximal stimulation EC50 (ng/mL)") # Table 4 EC50 = 50.74 (RSE 37.68%)
    hill  <- fixed(1);    label("Hill coefficient gamma on the stimulation term (unitless)")          # Results 3.3.1: 'when gamma was fixed at 1'

    # IIV - reported as a percentage in Table 3 / Table 4 and read
    # here as the standard deviation of the log-normal random effect
    # (omega x 100), so variance = (percent / 100)^2. Supplementary
    # Table 1 reports omega(Ka) and omega(Km) as 54.77% with no RSE
    # and no confidence interval in BOTH the final-model and the
    # 1,000-run bootstrap columns, which identifies them as held
    # constant rather than estimated; 0.5477^2 = 0.300.
    etalka   ~ fixed(0.3)       # Supplementary Table 1 'omega(Ka), %' = 54.77 in both final and bootstrap columns, no RSE / no CI
    etalvc   ~ 0.94245          # Table 3 'IIV V1/F' = 97.08% (RSE 24.30%): 0.9708^2
    etalkm   ~ fixed(0.3)       # Supplementary Table 1 'omega(Km), %' = 54.77 in both final and bootstrap columns, no RSE / no CI
    etalkout ~ 0.010547         # Table 4 'IIV Kout' = 10.27% (RSE 3.901%): 0.1027^2
    etalec50 ~ 0.620314         # Table 4 'IIV EC50' = 78.76% (RSE 30.66%): 0.7876^2

    # Residual error. Both tables print the point estimate on a
    # percentage-style scale while the main-text 95% CI column prints
    # the same quantity as a fraction (PK 28.94 vs 0.2730-0.3057;
    # PD 13.20 vs 0.1161-0.1479), so both are entered as fractions.
    propSd            <- 0.2894; label("Proportional residual error on YPEG-rhGH serum concentration (fraction)") # Table 3 'Residual error sigma (prop), %' = 28.94 (RSE 2.88%)
    addSd_IGF1ratio   <- 0.1320; label("Additive residual SD on the IGF-1 ratio to baseline (fraction of baseline)") # Table 4 'Residual error, sigma (add)' = 13.20 (RSE 6.13%)
  })

  model({
    # 1. Covariate terms - power model centred at the values printed
    #    in the supplement's final-model equations.
    age_norm <- AGE / 35.5
    wt_norm  <- WT / 70

    # 2. Individual parameters
    ka   <- exp(lka   + etalka) * age_norm^e_age_ka
    vc   <- exp(lvc   + etalvc) * age_norm^e_age_vc
    vp   <- exp(lvp)
    q    <- exp(lq)
    vmax <- exp(lvmax) * age_norm^e_age_vmax * wt_norm^e_wt_vmax
    km   <- exp(lkm   + etalkm)

    kin  <- exp(lkin)
    kout <- exp(lkout + etalkout)
    emax <- exp(lemax)
    ec50 <- exp(lec50 + etalec50)

    # 3. Micro-constants for the two-compartment disposition
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system. Serum concentration is in ug/L (= ng/mL), so the
    #    Michaelis-Menten term Vmax * Cc / (Km + Cc) is in ug/h and
    #    matches the amount units of the central compartment.
    Cc <- central / vc

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - k12 * central + k21 * peripheral1 -
                          vmax * Cc / (km + Cc)
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. IGF-1 indirect response on the ratio-to-baseline scale.
    #    Production is stimulated by serum YPEG-rhGH; the typical
    #    baseline ratio is Kin / Kout = 1.
    stim <- emax * Cc^hill / (ec50^hill + Cc^hill)

    d/dt(igf1) <- kin * (1 + stim) - kout * igf1
    igf1(0)    <- kin / kout

    # 6. Observations and residual error
    IGF1ratio <- igf1

    Cc        ~ prop(propSd)
    IGF1ratio ~ add(addSd_IGF1ratio)
  })
}
