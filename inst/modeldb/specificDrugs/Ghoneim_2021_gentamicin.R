Ghoneim_2021_gentamicin <- function() {
  description <- "Two-compartment intravenous population PK model for gentamicin in non-critically ill pediatric inpatients aged 1 month to 6 years, fitted to routine therapeutic-drug-monitoring peak and trough concentrations (Ghoneim 2021). Clearance and central volume scale as power functions of total body weight normalized to 70 kg (estimated exponents 0.71 and 0.93); peripheral volume and intercompartmental clearance carry no covariate. Between-subject variability on clearance and central volume; additive residual error."
  reference <- "Ghoneim RH, Thabit AK, Lashkar MO, Ali AS. Optimizing gentamicin dosing in different pediatric age groups using population pharmacokinetics and Monte Carlo simulation. Ital J Pediatr. 2021;47:167. doi:10.1186/s13052-021-01114-4"
  vignette <- "Ghoneim_2021_gentamicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # What each ODE state holds, in what amount units, in what biological
  # matrix. Gentamicin was given as a 30-minute IV infusion, so the dose
  # enters `central` directly and there is no depot state. The specimen is
  # verified: Ghoneim 2021 Abstract and Methods ('Patients and data') describe
  # the data as 'plasma gentamicin concentration data' / 'Plasma concentration
  # data', assayed by PETINIA immunoassay (Dimension, Dade Behring).
  compartmentData <- list(
    central = list(analyte = "gentamicin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "gentamicin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Ghoneim 2021 Table 1: mean (SD) 10.13 (5.25) kg, range 3.98-17.7 kg (n = 22). Reference weight 70 kg is the normalizing constant printed in both final-model equations in the Results ('6 x (weight in kg/70)^0.71' and '15.87 x (weight in kg/70)^0.93') and in the Table 3 / Results units 'L/hr./70 kg' and 'L/70 kg'. Weight enters CL and Vc only; Vp and Q carry no weight term (Abstract: 'weight incorporated as a significant covariate for both clearance and volume of distribution'). Adding weight on CL and Vc lowered -2LL by 61 points and AIC by 57 points (Results). The cohort weight range is 3.98-17.7 kg, so the 70 kg reference lies far outside the data.",
      source_name = "weight in kg"
    )
  )

  # Screened in the Ghoneim 2021 covariate analysis (Methods, 'Population
  # pharmacokinetic model development and evaluation': 'covariates including
  # age, gender, serum creatinine, and body weight were examined') but NOT
  # retained in the final model. Documentation only; not referenced in
  # model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (stepwise forward addition p < 0.01, backward deletion p < 0.001); not retained. Table 1 reports age in MONTHS: mean (SD) 34.88 (31.9), range 1-72 months (i.e. 0.08-6 years). Age is strongly collinear with weight over this range. The Discussion notes that post-gestational age was not documented for newborns, who were all recorded as 4 weeks old.",
      source_name = "age"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "male",
      notes = "Screened as 'gender'; not retained. Table 1: 13 of 22 (59.1%) male, so 40.9% female. The source reports male counts; canonical SEXF is the complement (SEXF = 1 - SEXM).",
      source_name = "gender"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained. Table 1 prints mean (SD) '0.39 +/- 0.82 (0.27 - 0.51)' mg/dL; the SD of 0.82 cannot be reconciled with a range of 0.27-0.51 and is presumably a typographical error (perhaps 0.082). Renal function was assessed by serum creatinine only; the Discussion notes that urine creatinine was not available.",
      source_name = "serum creatinine"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 22L,
    n_studies = 1L,
    n_sites = 1L,
    age_range = "1 to 72 months (Table 1 range; mean (SD) 34.88 (31.9) months). Inclusion allowed neonates to 12 years, but no enrolled patient was older than 6 years.",
    weight_range = "3.98 to 17.7 kg (Table 1; mean (SD) 10.13 (5.25) kg)",
    sex_female_pct = 40.9,
    race_ethnicity = "Not reported. Single-centre Saudi Arabian cohort.",
    disease_state = "Non-critically ill pediatric inpatients receiving intravenous gentamicin for empiric treatment of Gram-negative infections. Excluded: pediatric intensive care admissions, surgical prophylaxis, and co-administration of other nephrotoxic drugs. Serum creatinine mean 0.39 mg/dL (range 0.27-0.51), i.e. normal renal function.",
    dose_range = "Weight-based intravenous gentamicin given as a 30-minute infusion; Table 1: 2.26 (0.33) mg/kg per dose (range 1.78-2.73; the Results prose says 2.75 (0.33)), median dosing interval 8 h (range 8-12 h).",
    regions = "Saudi Arabia (King Abdulaziz University Hospital, Jeddah)",
    notes = "Retrospective chart review, February to November 2015. Sparse TDM sampling: a peak drawn 30 min after the end of the 30-min infusion of the third dose and a trough drawn just before the fourth dose. Observed peak 5.45 (1.08) mg/L and trough 0.58 (0.28) mg/L (Table 1). Fitted in Phoenix NLME 8.2; one- and two-compartment structures with additive, multiplicative and mixed residual error were tested. Model evaluation by 1000-replicate bootstrap (Table 3) and a 1000-replicate visual predictive check (Figure 3)."
  )

  ini({
    # Structural parameters: Ghoneim 2021 Table 3 'Parameter estimates for
    # the final model and bootstrap', column 'Mean estimate of the final
    # model +/- SE'. CL and Vc are the typical values for a 70 kg subject
    # (Results: 'clearance was estimated to be 4.64 L/hr./70 kg'; 'The
    # average volume of the central compartment was 15.87 L/70 kg').
    #
    # NOTE ON THE CLEARANCE INTERCEPT (4.64 versus 6). The same Results
    # sentence that states 'clearance was estimated to be 4.64 L/hr./70 kg'
    # goes on to print the equation '6 x (weight in kg/70)^0.71'. This
    # extraction uses 4.64 because every other printing agrees on it:
    #   (1) Table 3 point estimate 4.64 +/- 0.56 L/h;
    #   (2) Table 3 BOOTSTRAP MEDIAN 4.64 (95% CI 3.58-7.34), an independent
    #       computation that matches the point estimate to two decimals for
    #       every row (Vc 15.87 vs 15.86, Vp 4.11 vs 4.11, Q 0.62 vs 0.63);
    #   (3) the prose value 4.64 L/hr./70 kg in the same sentence;
    #   (4) the Discussion, 'our final estimates of 4.6 L/hr. and 15 L'.
    # The companion Vc equation prints the Table 3 intercept (15.87)
    # unchanged, so the equations restate the table intercepts and '6' is the
    # lone deviation. See the vignette section 'Adjudicating the clearance
    # intercept' and Assumptions and deviations.
    lcl <- log(4.64); label("Clearance for a 70 kg subject (L/h)")                  # Ghoneim 2021 Table 3: CL 4.64 +/- 0.56 L/h (bootstrap median 4.64); Results: 4.64 L/hr./70 kg
    lvc <- log(15.87); label("Central volume of distribution for a 70 kg subject (L)") # Ghoneim 2021 Table 3: Vc 15.87 +/- 3.99 L (bootstrap median 15.86); Results: 15.87 L/70 kg
    lvp <- log(4.11); label("Peripheral volume of distribution (L)")                # Ghoneim 2021 Table 3: Vp 4.11 +/- 0.99 L (bootstrap median 4.11)
    lq <- log(0.62); label("Intercompartmental clearance (L/h)")                    # Ghoneim 2021 Table 3: Q 0.62 +/- 0.11 L/h (bootstrap median 0.63)

    # Body-weight exponents, printed only inside the Results equations
    # '(weight in kg/70)^0.71' (CL) and '(weight in kg/70)^0.93' (Vc).
    # Neither is a theoretical allometric value (0.75 / 1) and the Results
    # say weight was 'incorporated' as an estimated covariate effect, so they
    # are treated as estimated and NOT wrapped in fixed(). Table 3 does not
    # list them, so no SE is available.
    e_wt_cl <- 0.71; label("Power exponent on (WT/70) for CL (unitless)")  # Ghoneim 2021 Results: CL equation (weight in kg/70)^0.71
    e_wt_vc <- 0.93; label("Power exponent on (WT/70) for Vc (unitless)")  # Ghoneim 2021 Results: Vc equation (weight in kg/70)^0.93

    # Between-subject variability. Table 3 reports 'Between subject
    # variability associated with Vc (%)' = 37.80% and '... with CL (%)' =
    # 27.89%. Read as a coefficient of variation of a log-normal parameter,
    # so omega^2 = log(CV^2 + 1). No covariance between the two etas is
    # reported, so they are independent.
    etalcl ~ 0.074908  # log(0.2789^2 + 1); Ghoneim 2021 Table 3: BSV CL 27.89% (bootstrap 27.00%)
    etalvc ~ 0.133555  # log(0.378^2 + 1); Ghoneim 2021 Table 3: BSV Vc 37.80% (bootstrap 35.53%)

    # Residual variability. Table 3 'Additive error (mg/L)' = 0.011 (bootstrap
    # 0.011); the Results and Abstract state 'additive residual error'.
    # Phoenix NLME reports the standard deviation of the additive epsilon in
    # concentration units, so this is addSd. (The base model's additive error
    # was 0.28 mg/L, Table 2.)
    addSd <- 0.011; label("Additive residual error (mg/L)")  # Ghoneim 2021 Table 3: additive error 0.011 mg/L (bootstrap 0.011)
  })

  model({
    # Individual PK parameters. Ghoneim 2021 Results:
    #   CL = 4.64 * (WT/70)^0.71   (printed with intercept 6; see ini() note)
    #   Vc = 15.87 * (WT/70)^0.93
    # Vp and Q carry no covariate and no BSV (Table 3).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment IV infusion model with first-order elimination
    # (Results: 'A two-compartment IV infusion model with additive residual
    # error and first-order elimination'). The infusion rate or duration is
    # set on the dose record.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, vc in L, so central/vc is mg/L.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
