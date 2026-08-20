Shen_2025_voriconazole <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption and first-order elimination for oral voriconazole in immunocompromised children under 2 years of age (Shen 2025); body weight enters apparent clearance and apparent volume through fixed allometric exponents, and the absorption rate constant is fixed to a published literature value"
  reference <- "Shen L, Hu M, Xu X, Zhou Y, Wu W, Ge X, Wang G, Wang Y, Li Z. Precision dosing of voriconazole in immunocompromised children under 2 years: integrated machine learning and population pharmacokinetic modeling. Front Pharmacol. 2025;16:1671652. doi:10.3389/fphar.2025.1671652"
  vignette <- "Shen_2025_voriconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Both states verified against Shen 2025 section 2.1
  # (oral voriconazole) and section 2.2 (plasma HPLC assay).
  compartmentData <- list(
    depot   = list(analyte = "voriconazole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "voriconazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Sole retained covariate. Enters CL/F as (WT/70)^0.75 and V/F as (WT/70)^1 with a 70 kg reference weight (Shen 2025 Equations 1 and 2). Both exponents were fixed at the canonical allometric values rather than estimated: Supplementary Table S3 reports no RSE, no bootstrap median, and no confidence interval for either exponent, and the Methods state only that 'Using allometric scaling, we examined the effects of body weight on PK parameters.' The 70 kg reference and the two exponents are confirmed exactly by the paper's own individual empirical-Bayes estimates: an 8.00 kg subject gives 788 * (8/70) = 90.06 L against the reported median V of 90.07 L/70kg (Supplementary Table S1), and an 11 kg subject gives 17.9 * (11/70)^0.75 = 4.47 L/h and 788 * (11/70) = 123.8 L against the reported CL/F 4.46 L/h and V/F 123.83 L in the section 3.5 worked example.",
      source_name        = "WT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Postnatal age",
      units       = "months",
      type        = "continuous",
      notes       = "Screened in the stepwise forward-addition / backward-elimination covariate search (Shen 2025 section 2.3) but not retained in the final PopPK model; age was retained only as an input feature of the downstream XGBoost machine-learning layer, which is not part of this pharmacokinetic model."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Collected as part of the laboratory panel (Shen 2025 Table 1) and screened during covariate selection, but not retained in the final PopPK model."
    ),
    EGFR = list(
      description = "Estimated glomerular filtration rate (Schwartz formula, k = 0.45 for age < 1 year and k = 0.413 for age 1-2 years)",
      units       = "mL/min/1.73m2",
      type        = "continuous",
      notes       = "Derived and reported for the study population (Shen 2025 Table 1 and section 2.1) and screened during covariate selection, but not retained in the final PopPK model. Voriconazole is hepatically cleared, so no renal-function effect on CL/F was expected or found."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 76L,
    n_studies      = 1L,
    n_observations = 110L,
    age_range      = "under 24 months (newborns and preterm infants excluded)",
    age_median     = "11.0 months (IQR 7.38-17.00; mean 12.08, SD 5.94)",
    weight_range   = "IQR 6.95-9.00 kg (full range not reported)",
    weight_median  = "8.05 kg (mean 8.00, SD 1.91)",
    sex_female_pct = 23.7,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Immunocompromised hospitalized infants and toddlers (haematologic malignancy or post-haematopoietic stem cell transplantation) receiving oral voriconazole for suspected or documented invasive fungal infection.",
    dose_range     = "Oral voriconazole 4-10 mg/kg every 12 h, individualized by the treating physician against the 9 mg/kg q12h regimen recommended for children 2 to <12 years. Median total daily dose 100.00 mg (IQR 100.00-133.25); median therapy duration at TDM 9.50 days (IQR 5.00-17.75).",
    regions        = "Single center: Children's Hospital of Fudan University, National Children's Medical Center, Shanghai, China.",
    co_medication  = c(glucocorticoids_pct = 61.8, proton_pump_inhibitor_pct = 44.7, tacrolimus_pct = 32.9, cyclosporine_A_pct = 13.3, sirolimus_pct = 2.6),
    notes          = "Retrospective observational single-center study, January 2020 - June 2025. 110 steady-state trough therapeutic drug monitoring (TDM) samples drawn within the 30 min immediately preceding the next scheduled dose. HPLC-UV assay, linear range 0.25-10 mg/L; concentrations below the 0.25 mg/L LLOQ (<4% of observations) were imputed at LLOQ/2 = 0.125 mg/L. Baseline demographics per Shen 2025 Table 1; final PopPK parameter estimates and bootstrap validation per Shen 2025 Supplementary Table S3. A separate cohort of 10 patients sampled July-August 2025 was used for external validation of the downstream machine-learning layer only. The paper's XGBoost concentration-prediction model, which consumes the empirical-Bayes CL and V produced by this PopPK model as input features, is not a pharmacokinetic structural model and is not represented here."
  )

  ini({
    # Structural parameters. Reference values are for a 70 kg individual;
    # the study population is 6.95-9.00 kg (IQR), so the reference is an
    # extrapolated allometric anchor rather than an observed subject.
    lka <- fixed(log(1.19)); label("Absorption rate constant (1/h)")                        # Shen 2025 Supplementary Table S3 (Ka 1.19, marked fixed); fixed because the trough-only sampling design could not identify it, value taken by Shen 2025 from Gastine 2018
    lcl <- log(17.9);        label("Apparent clearance CL/F at 70 kg (L/h)")           # Shen 2025 Supplementary Table S3 (17.9, RSE 10.8%) and Equation 1
    lvc <- log(788);         label("Apparent volume of distribution V/F at 70 kg (L)") # Shen 2025 Supplementary Table S3 (788, RSE 15.4%) and Equation 2

    # Allometric exponents. Fixed at the canonical values: Supplementary
    # Table S3 reports no RSE, bootstrap median, or confidence interval for
    # either exponent, and both appear as literal constants in the printed
    # equations.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F (unitless)") # Shen 2025 Equation 1: CL/F = 17.9 * (WT/70)^0.75
    e_wt_vc <- fixed(1);    label("Allometric exponent on V/F (unitless)")  # Shen 2025 Equation 2: V/F = 788 * (WT/70)

    # IIV. Shen 2025 estimated exponential interindividual variability on
    # CL/F only; section 3.2 states IIV was NOT estimated for V/F because
    # shrinkage exceeded 90% and the RSE exceeded 30%, so V/F carries no
    # eta here. The reported quantity is already a variance on the log
    # scale (omega^2 = 0.674, 'indicating a standard deviation of
    # approximately 0.821'), so no CV-to-variance conversion is needed.
    etalcl ~ 0.674 # Shen 2025 Supplementary Table S3 (omega^2 CL/F = 0.674, RSE 10.7%, bootstrap median 0.672, 95% CI 0.175-0.976) and section 3.2

    # Residual error. Supplementary Table S3 reports a single residual
    # variance, sigma^2 = 0.16 (RSE 18.4%; bootstrap median 0.152, 95% CI
    # 0.064-0.282). Exactly one sigma is tabulated, which rules out the
    # COMBINED model arithmetically and leaves proportional vs additive;
    # either way the magnitude is sqrt(0.16) = 0.4. Shen 2025 section 2.3
    # says proportional, additive and combined were all evaluated but never
    # names the one retained, and no $ERROR block is published.
    #
    # Read as PROPORTIONAL on a falsifier in Supplementary Figure S2 panel
    # A (observed vs INDIVIDUAL-predicted), which discriminates the two
    # forms directly because the residual is the only source of scatter
    # left in that panel. The ~20 observations with IPRED < 0.5 mg/L lie on
    # the identity line within roughly 0.15 mg/L. An additive SD of
    # 0.4 mg/L would scatter those same points by about +/-0.8 mg/L and
    # would visibly truncate them at zero; every one of them landing inside
    # 0.4 SD is not a plausible draw. A 40% CV predicts an SD of 0.12 mg/L
    # at IPRED = 0.3, which is what the panel shows. The same panel's
    # scatter then grows with concentration (about +/-0.8 mg/L around
    # IPRED = 3.5, i.e. ~23%), i.e. roughly constant on a RELATIVE scale,
    # which is the proportional signature and not the additive one.
    # Supporting, non-decisive: ~4% of observations were BQL and imputed at
    # 0.125 mg/L, and the observations span 0.125 to ~5.3 mg/L (>40-fold),
    # a range no constant additive SD covers. See the vignette Errata.
    propSd <- 0.4; label("Proportional residual error (fraction)") # Shen 2025 Supplementary Table S3 (sigma^2 = 0.16); propSd = sqrt(0.16)
  })

  model({
    # Individual parameters. Allometric scaling on body weight with a 70 kg
    # reference (Shen 2025 Equations 1 and 2). Exponential IIV on CL/F only.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    # One-compartment disposition with first-order absorption; NONMEM
    # ADVAN2 TRANS2 (Shen 2025 section 2.3). CL and V are apparent
    # (CL/F, V/F) because only oral dosing was studied, so bioavailability
    # is not separately identifiable and is absorbed into both parameters.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose in mg divided by volume in L gives mg/L, which is the unit of
    # the reported trough concentrations and equals ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
