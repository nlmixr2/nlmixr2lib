Ye_2025_mycophenolic_acid <- function() {
  description <- "Population PK model for total mycophenolic acid (MPA, the active moiety of mycophenolate mofetil MMF) in paediatric patients with lupus nephritis receiving oral MMF twice daily (Ye 2025). Two-compartment disposition with first-order absorption, an absorption lag time and linear elimination, fitted in Phoenix NLME. Body weight is the only retained covariate and enters the apparent peripheral volume as a power term Vp/F = 1287.12 * (WT/41.13)^2.05; no covariate was retained on clearance, so steady-state exposure in this model is independent of body weight. Inter-individual variability is log-normal on Vc/F, Vp/F and Q/F only -- the random effects on ka, CL/F and Tlag were dropped for high eta-shrinkage. Residual error is combined proportional plus additive. Doses are MMF mass (mg) with no MMF-to-MPA molecular-weight conversion: CL/F, Vc/F, Vp/F and Q/F are apparent parameters absorbing both the molecular-weight ratio and oral bioavailability."
  reference <- paste(
    "Ye C, Liu B, Chen L, Zhang L, Zheng Y, Tang K, Jiang X, Chen P.",
    "Impact of body weight on mycophenolic acid population pharmacokinetics",
    "in paediatric lupus nephritis: a pharmacogenomic integration study.",
    "Lupus Science & Medicine. 2025;12(1):e001535.",
    "doi:10.1136/lupus-2025-001535.",
    sep = " "
  )
  vignette <- "Ye_2025_mycophenolic_acid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Ye 2025 Methods (plasma MPA by
  # LC-MS/MS) and the two-compartment-plus-lag structure in Results.
  compartmentData <- list(
    depot       = list(analyte = "mycophenolic acid", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "mycophenolic acid", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "mycophenolic acid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The only covariate retained in the final model. Power effect on the apparent peripheral volume of distribution: Vp/F = 1287.12 * (WT/41.13)^2.05. The normalising weight 41.13 kg is printed inside the covariate equation in Ye 2025 online supplemental table 5 and equals the model-group mean weight reported in Table 1 (41.13 +/- 11.31 kg). Ye 2025 online supplemental table 5 prints the relationship as V2,i = V2,pop * (body weight / 41.13) with the estimated coefficient omitted; the same table omits the estimated coefficient from every other screened relationship as well (e.g. V_i = V_pop * e^(dsDNA) with no theta), so the omission is that table's notation and not a claim that the exponent is 1. Table 3 supplies the coefficient as theta_v2,weight = 2.05 (RSE 11.99%, bootstrap 95% CI 0.39-2.60), which is read here as the power exponent -- see vignette Errata. Body weight was NOT retained on clearance, so steady-state AUC in this model does not depend on weight. Ye 2025 also screened gender, age, white blood cell count, blood neutrophil count, platelet count, lymphocyte count, haematocrit, haemoglobin, ALT, AST, alkaline phosphatase, total bilirubin, albumin, glucose, blood urea nitrogen, SLEDAI-2K and cystatin C; none was retained. Creatinine clearance was listed as a candidate but could not be computed because 24-hour urine collection was impractical in this largely outpatient cohort (Discussion).",
      source_name        = "body weight"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 51L,
    n_studies      = 1L,
    n_observations = 1170L,
    n_profiles     = 146L,
    age_mean       = "12.22 +/- 2.33 years (model group); 13.90 +/- 2.18 years (external validation group)",
    weight_mean    = "41.13 +/- 11.31 kg (model group); 39.54 +/- 7.14 kg (external validation group)",
    sex_female_pct = 90.2,
    race_ethnicity = "Not tabulated. Single-centre Chinese study (The First Affiliated Hospital, Sun Yat-sen University, Guangzhou); candidate SNPs were selected on a minor-allele-frequency > 5% threshold in Han Chinese reference data, so the cohort is Han Chinese.",
    disease_state  = "Paediatric lupus nephritis. All patients met the 2019 American College of Rheumatology SLE classification criteria and the diagnostic criteria for lupus nephritis (persistent proteinuria >= 0.5 g/day, active cellular casts, or biopsy evidence). SLEDAI-2K < 5 (inactive) in 78/146 profiles (53.4%), >= 5 (active) in 44/146 (30.1%), not graded in 24/146 (16.4%). Mean plasma albumin 41.76 g/L, with only 5/146 profiles (3.4%) below the 31 g/L threshold at which free rather than total MPA monitoring is recommended.",
    dose_range     = "Oral mycophenolate mofetil 125-750 mg every 12 hours as capsule (CellCept, 41 profiles) or dispersible tablet (Saikeping, 105 profiles); mean total daily dose 739.81 +/- 247.37 mg/day. Doses in this model file are MMF mass in mg with no MMF-to-MPA molecular-weight conversion applied.",
    regions        = "Single centre: The First Affiliated Hospital, Sun Yat-sen University, Guangzhou, China.",
    co_medication  = "Corticosteroids in 51/51 patients (100%), mean 12.86 mg/kg/day; tacrolimus 9/51 (17.6%); hydroxychloroquine 18/51 (35.3%); biological agents 7/51 (13.7%).",
    notes          = "Prospective study, September 2021 - January 2023, with a separate nine-patient external validation cohort recruited under the same criteria (MDPE 9.09%, MAPE 24.82%, F20 52.5%, F30 72.5%). Sampling was at steady state (MMF for at least 7 days; mean treatment duration 48.66 days) pre-dose and at 0.5, 1.5, 2.5, 4, 6, 9 and 12 hours post-dose -- eight samples per profile, consistent with the reported 1170 samples across 146 profiles (1170/146 = 8.01). Plasma total MPA was assayed by validated LC-MS/MS, calibration range 0.1-50 ug/mL, LLOQ 0.1 ug/mL. Estimation was first-order conditional estimation-extended least squares in Phoenix NLME 8.3. Twenty-nine SNPs across 13 candidate genes were genotyped and eight variants in UGT1A9 (rs6717546, rs13418420, rs7586110, rs2070959, rs6759892), UGT2B7 (rs7438135), ABCC2 (rs7910642) and CES1 (rs12149373) were associated with MPA exposure in the univariate analysis, but no genetic covariate survived the stepwise popPK covariate search and none appears in the final model. Enterohepatic recirculation of MPA could not be modelled (no samples around the second peak, few patients with a pronounced EHC process). Baseline demographics per Ye 2025 Table 1; final-model parameter estimates per Ye 2025 Table 3; covariate functional form per Ye 2025 online supplemental table 5."
  )

  ini({
    # Structural parameters. Ye 2025 Table 3, "Final model" column.
    #
    # Table 3 heads this row "theta_Ka (hours)", but ka is a first-order
    # absorption rate constant and the Discussion compares it against
    # literature values quoted "from 0.39 to 5.21 hours^-1". The header
    # unit is a typo for 1/h; the value is used unchanged.
    lka   <- log(1.54);    label("Absorption rate constant (1/h)")                        # Ye 2025 Table 3, theta_Ka = 1.54 (RSE 10.25%); bootstrap median 1.69 (95% CI 1.42-2.06)

    lvc   <- log(7.99);    label("Apparent central volume of distribution (L)")           # Ye 2025 Table 3, theta_V/F = 7.99 L (RSE 22.76%); bootstrap median 8.31 (95% CI 1.67-12.68)

    # Typical value at the covariate reference weight of 41.13 kg, at
    # which the weight factor equals 1.
    lvp   <- log(1287.12); label("Apparent peripheral volume of distribution at WT = 41.13 kg (L)")  # Ye 2025 Table 3, theta_V2/F = 1287.12 L (RSE 8.65%); bootstrap median 1116.28 (95% CI 841.72-1682.86)

    # No covariate was retained on clearance (Ye 2025 Results, "'body
    # weight' emerged as the sole significant covariate in the final
    # model, demonstrating a positive relationship with V2").
    lcl   <- log(15.23);   label("Apparent clearance (L/h)")                              # Ye 2025 Table 3, theta_CL/F = 15.23 L/h (RSE 1.85%); bootstrap median 16.27 (95% CI 14.39-18.40)

    lq    <- log(37.65);   label("Apparent intercompartmental clearance (L/h)")           # Ye 2025 Table 3, theta_CL2/F = 37.65 L/h (RSE 3.32%); bootstrap median 38.86 (95% CI 30.26-46.83)

    ltlag <- log(0.40);    label("Absorption lag time (h)")                               # Ye 2025 Table 3, Tlag = 0.40 h (RSE 3.29%); bootstrap median 0.42 (95% CI 0.17-0.49)

    # Covariate effect: power of body weight on the apparent peripheral
    # volume. Ye 2025 online supplemental table 5 prints the form as
    #   V2,i = V2,pop * (body weight / 41.13)
    # with the estimated coefficient omitted, exactly as it omits the
    # coefficient from every other screened relationship in that table
    # (the genotype rows print V_i = V_pop * e^(dsDNA) with no theta
    # either). Table 3 supplies the coefficient, theta_v2,weight = 2.05,
    # which is read here as the exponent of the printed ratio. See the
    # vignette Errata for the alternative readings considered and why a
    # centred-linear form was rejected (it turns Vp negative below
    # 41.13 * (1 - 1/2.05) = 21.1 kg, inside this cohort's weight range).
    e_wt_vp <- 2.05;       label("Power exponent on (WT/41.13) for the apparent peripheral volume (unitless)")  # Ye 2025 Table 3, theta_v2,weight = 2.05 (RSE 11.99%); bootstrap median 2.05 (95% CI 0.39-2.60); form per online supplemental table 5

    # Inter-individual variability. Ye 2025 Table 3 labels these rows
    # omega^2, and Phoenix NLME reports the Omega matrix of an
    # exponential IIV model as variances, so the tabulated values are
    # used directly as eta variances.
    #
    # WHICH PARAMETERS CARRY AN ETA. The text states that "the typical
    # value of absorption rate constant (Ka), typical value of clearance
    # (CL) and Tlag were fixed due to their high eta-shrinkage values",
    # explaining that "high shrinkage ... suggests that the random effect
    # may be negligible". That reasoning removes a RANDOM EFFECT, not a
    # THETA, and all three thetas carry non-zero RSEs and bootstrap
    # medians that differ from the final estimates (1.69 vs 1.54, 16.27
    # vs 15.23, 0.42 vs 0.40) -- which a genuinely fixed THETA cannot.
    # So the etas on ka, CL/F and Tlag were fixed to zero (dropped), and
    # the three thetas were estimated. That leaves exactly the three
    # omega^2 rows Table 3 prints: V/F, V2/F and CL2/F.
    #
    # Table 3's shrinkage column is internally inconsistent on the third
    # eta: the "Final model" column prints a shrinkage beside theta_CL/F
    # (18.30) and none beside theta_CL2/F, while the "Bootstrap" column
    # does the opposite (13.60 beside CL2/F, none beside CL/F). The IIV
    # rows and the table footnote both say CL2/F, and the text says
    # CL's random effect was dropped, so the eta is placed on Q/F here
    # and the stray "(18.30)" is treated as a mis-set cell. Recorded in
    # the vignette Errata.
    etalvc ~ 1.02   # Ye 2025 Table 3, omega^2 V/F = 1.02 (RSE 14.89%); bootstrap median 1.10 (95% CI 0.98-1.22)
    etalvp ~ 1.24   # Ye 2025 Table 3, omega^2 V2/F = 1.24 (RSE 13.80%); bootstrap median 1.25 (95% CI 1.01-1.48)
    etalq  ~ 0.20   # Ye 2025 Table 3, omega^2 CL2/F = 0.20 (RSE 17.51%); bootstrap median 0.23 (95% CI 0.16-0.30)

    # Residual error: the mixed (combined additive plus multiplicative)
    # structure selected in Ye 2025 Results, "we selected a mixed error
    # structure to characterise intraindividual variability (AIC:
    # 4659.74 for mixed error, compared with 4738.81 and 4742.22 for
    # additive and multiplicative errors, respectively)". Table 3's
    # footnote glosses the two rows as "MultStdev, multiplicative SE"
    # and "stdev, additive SE", i.e. both are standard deviations.
    propSd <- 0.50;        label("Proportional residual error (fraction)")                # Ye 2025 Table 3, MultStdev = 0.50 (RSE 4.29%); bootstrap median 0.52 (95% CI 0.48-0.56)
    addSd  <- 0.38;        label("Additive residual error (mg/L)")                       # Ye 2025 Table 3, stdev = 0.38 (RSE 9.33%); bootstrap median 0.42 (95% CI 0.20-0.56)
  })

  model({
    # No IIV on ka, CL/F or Tlag -- see the ini() comment on which
    # parameters carry an eta.
    ka   <- exp(lka)
    cl   <- exp(lcl)
    tlag <- exp(ltlag)

    vc <- exp(lvc + etalvc)

    # Power effect of body weight on the apparent peripheral volume,
    # normalised to the 41.13 kg reference printed in Ye 2025 online
    # supplemental table 5.
    vp <- exp(lvp + etalvp) * (WT / 41.13)^e_wt_vp

    q  <- exp(lq + etalq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Bioavailability is not identifiable from an oral-only dataset, so F
    # is absorbed into the apparent parameters CL/F, Vc/F, Vp/F and Q/F
    # and no f(depot) term is applied. Doses are MMF mass in mg; the
    # MMF-to-MPA molecular-weight ratio is likewise absorbed into the
    # apparent parameters (Ye 2025 applies no molecular-weight
    # conversion -- see vignette Errata).
    #
    # Observation: dose in mg, vc in L -> central/vc in mg/L, numerically
    # identical to the ug/mL used in Ye 2025 Table 2.
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
