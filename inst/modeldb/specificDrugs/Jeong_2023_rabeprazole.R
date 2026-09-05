Jeong_2023_rabeprazole <- function() {
  description <- "Co-linked population PK-PD model for a single oral 10 mg rabeprazole enteric-coated tablet in 45 healthy Korean adults (24 men, 21 women) from a two-way crossover bioequivalence study (Jeong 2023). PK is a two-compartment disposition model (central + peripheral1) fed by a chain of three sequential first-order absorption compartments (depot -> transit1 -> transit2 -> central, with distinct rate constants Ka1, Ka2 and Ka3) preceded by an absorption lag time Tlag; volumes and clearances are apparent (Vc/F, CLc/F, Vp/F, CLp/F). Two covariates were retained: female sex increases Tlag by a fraction of 0.73 (a linear, not exponential, effect), and body surface area lowers the third absorption rate constant Ka3 through a power model centred on the cohort median BSA of 1.75 m^2. Residual error is log-additive, i.e. log-normal on the linear concentration scale. PD is a direct-response sigmoid Emax model with baseline for the rise in intragastric pH driven by the model-predicted plasma rabeprazole concentration: pH = E0 + Emax * C^gamma / (EC50^gamma + C^gamma). The PD observations were digitised from Chen 2006 rather than measured in this cohort, and no covariate or inter-individual variability was estimated on any PD parameter."
  reference <- paste(
    "Jeong SH, Jang JH, Lee YB.",
    "Exploring Differences in Pharmacometrics of Rabeprazole between Genders via Population Pharmacokinetic-Pharmacodynamic Modeling.",
    "Biomedicines. 2023 Nov 10;11(11):3021.",
    "doi:10.3390/biomedicines11113021.",
    "The gastric-pH observations used to fit the pharmacodynamic sub-model were digitised from",
    "Chen ZY, Xie HT, Zheng QS, Sun RY, Hu G.",
    "Pharmacokinetic and pharmacodynamic population modeling of orally administered rabeprazole in healthy Chinese volunteers by the NONMEM method.",
    "Eur J Drug Metab Pharmacokinet. 2006;31(1):27-33.",
    sep = " "
  )
  vignette <- "Jeong_2023_rabeprazole"

  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot       = list(analyte = "rabeprazole", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "rabeprazole", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "rabeprazole", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "rabeprazole", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "rabeprazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    SEXF = list(
      description        = "Biological sex indicator; 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Time-fixed. Supplementary Information S7 prints the retained relationship as",
        "Tlag = tvTlag * (1 + dTlagdGender * (if female = 1 and male = 0)) * exp(eta_Tlag),",
        "so the sex term is LINEAR in (1 + 0.73 * SEXF) and NOT an exponential",
        "exp(0.73 * SEXF) factor; the typical female Tlag is 1.56 * 1.73 = 2.70 h against",
        "1.56 h in men. Adding sex to Tlag lowered the objective function by 157.97",
        "(Table 4 model 5), by far the largest single-covariate drop. Sex was also screened",
        "on Ka1, Ka2, Ka3, Vc/F and CLc/F and retained on none of them (Table 4 models 2, 3,",
        "4, 6 and 7, all dOFV between +0.70 and -1.87). The authors attribute the longer",
        "female lag to slower gastric emptying of solids delaying arrival of the",
        "enteric-coated tablet at its small-intestinal absorption window (Discussion).",
        sep = " "
      ),
      source_name        = "Gender (Jeong 2023 Table 5 'dTlagdGender'; Supplementary Information S7)"
    ),
    BSA = list(
      description        = "Body surface area computed by the Mosteller formula, sqrt(height (cm) * body weight (kg) / 3600)",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline. Supplementary Information S7 prints",
        "Ka3 = tvKa3 * (BSA / median BSA)^dKa3dBSA * exp(eta_Ka3), i.e. a power model centred",
        "on the median BSA of the observed population, which Table S1 gives as 1.75 m^2",
        "(mean 1.76, range 1.41-2.16). The exponent is negative (-1.11), so Ka3 rises as BSA",
        "falls; with the male and female median BSA of 1.87 and 1.58 m^2 used in the paper's",
        "own simulations (Results 3.4) the female Ka3 is about 21 percent higher than the",
        "male value. BSA was screened on Ka1, Ka2, Tlag, Vc/F and CLc/F as well; the Vc/F and",
        "CLc/F effects met the forward-selection criterion (dOFV -3.95 and -3.98) but failed",
        "backward elimination at p < 0.01 and were dropped (Table 4 models 8, 9, 11, 12, 13).",
        "Supplementary Information S1 states the Mosteller derivation; the printed formula",
        "lost its radical in the PDF text layer, but the reported mean height 167.83 cm and",
        "weight 66.66 kg give sqrt(167.83 * 66.66 / 3600) = 1.76 m^2, matching Table S1.",
        sep = " "
      ),
      source_name        = "BSA (Jeong 2023 Table 5 'dKa3dBSA'; Table S1; Supplementary Information S1 and S7)"
    )
  )

  covariatesDataExcluded <- list(
    BMI = list(
      description        = "Body mass index at baseline",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened on Vc/F and CLc/F (Table 4 models 17 and 18; dOFV -3.11 and -3.08). Both",
        "met the forward-selection threshold of 3.84 only marginally and neither survived",
        "backward elimination at p < 0.01, so BMI was not retained. Cohort mean (SD)",
        "23.54 (2.93) kg/m^2, range 18.20-29.50 (Table S1).",
        sep = " "
      ),
      source_name        = "BMI (Jeong 2023 Table 4; Table S1)"
    ),
    GGT = list(
      description        = "Serum gamma-glutamyl transpeptidase activity at baseline",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened on Vc/F and CLc/F (Table 4 models 15 and 16; dOFV +0.74 and -0.18) and not",
        "retained. GGT differed significantly between the sexes (24.13 vs 13.71 U/L,",
        "p = 1.73e-3; Table S2), but the Discussion concludes that liver-function indicators",
        "within the normal range of a healthy cohort have a negligible effect on rabeprazole",
        "pharmacokinetics. Cohort mean (SD) 19.27 (11.80) U/L, range 6.00-62.00 (Table S1).",
        sep = " "
      ),
      source_name        = "GTP (Jeong 2023 Table 4; Table S1)"
    ),
    CRCL = list(
      description        = "Creatinine clearance by the Cockcroft-Gault equation, not body-surface-area normalised",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened on Ka1, Ka2 and Ka3 (Table 4 models 19, 20 and 21; dOFV +2.27, +1.46 and",
        "-0.51) and not retained. Supplementary Information S1 gives the derivation as",
        "((140 - age) * body weight (kg) / serum creatinine (mg/dL) * 72), i.e. raw",
        "Cockcroft-Gault in mL/min rather than the BSA-normalised mL/min/1.73 m^2 form that",
        "the CRCL canonical usually carries. Cohort mean (SD) 123.70 (24.25) mL/min",
        "(Table S1). The Discussion notes that renal function was not expected to matter:",
        "rabeprazole pharmacokinetics are unchanged even in severe renal failure.",
        sep = " "
      ),
      source_name        = "CrCL (Jeong 2023 Table 4; Table S1; Supplementary Information S1)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 45L,
    n_studies      = 1L,
    age_range      = "20-49 years (mean (SD) 32.31 (8.40); median 31.00)",
    age_median     = "31 years",
    weight_range   = "45.60-94.10 kg (mean (SD) 66.66 (11.98); median 66.00; men 74.05 (9.90), women 58.21 (7.90))",
    weight_median  = "66 kg",
    height_range   = "150.50-183.90 cm (mean (SD) 167.83 (8.59); median 168.10)",
    bmi_range      = "18.20-29.50 kg/m^2 (mean (SD) 23.54 (2.93))",
    bsa_range      = "1.41-2.16 m^2 (mean (SD) 1.76 (0.19); median 1.75; men 1.89 (0.15), women 1.61 (0.12))",
    sex_female_pct = 46.7,
    race_ethnicity = "Korean",
    disease_state  = "Healthy adult volunteers",
    dose_range     = "Single oral 10 mg rabeprazole enteric-coated tablet",
    regions        = "Republic of Korea (single centre; Ministry of Food and Drug Safety approval MB22-002)",
    renal_function = "Normal; creatinine clearance mean (SD) 123.70 (24.25) mL/min, MDRD eGFR 106.89 (18.25) mL/min/1.73 m^2 (Table S1)",
    notes          = paste(
      "Randomised, single-dose, open-label, two-way crossover bioequivalence study with a",
      "7-day washout (Supplementary Information S3). All 45 subjects fasted for more than",
      "10 h and took the tablet with 150 mL of water. Plasma sampling at 0, 1, 2, 2.5, 3,",
      "3.5, 4, 4.5, 5, 5.5, 6, 7, 8 and 10 h post-dose. Sex was designed into the trial, so",
      "the male (n = 24) and female (n = 21) strata are balanced by construction. Rabeprazole",
      "was not detected in the plasma of ANY woman at the 1 h sample, whereas it was",
      "measurable in men at that time (Results 3.1) -- the observation the Tlag sex effect",
      "encodes. The pharmacodynamic sub-model was NOT fitted to data from this cohort: the",
      "gastric-pH versus plasma-concentration observations were digitised with WebPlotDigitizer",
      "4.6 from Chen 2006 (Methods 2.5).",
      sep = " "
    )
  )

  ini({
    # ========================================================================
    # PK. Two-compartment disposition fed by three sequential first-order
    # absorption compartments with an absorption lag time. Structure from
    # Jeong 2023 Results 3.2 and Table 3 (model 02-03-02-09); estimates from
    # Table 5; equations from Supplementary Information S7. All disposition
    # parameters are apparent (divided by the unknown oral bioavailability F),
    # so no separate bioavailability term is estimated.
    #
    # The chain is dose -> depot -(Ka1)-> transit1 -(Ka2)-> transit2 -(Ka3)->
    # central, matching the Table 3 footnote *** legend "Ka1: dosing
    # depot-depot 1, Ka2: depot 1-depot 2, Ka3: depot 2-central compartment".
    # The three rate constants are estimated separately, so this is a
    # sequential absorption chain rather than a Savic equal-ktr transit chain.
    # ========================================================================
    lka1  <- log(1.91);  label("First sequential absorption rate constant Ka1, depot to transit1 (1/h)")            # Jeong 2023 Table 5: tvKa1 = 1.91 1/h (SE 0.40, CV 21.17%; bootstrap median 1.89, 95% CI 1.37-3.11)
    lka2  <- log(2.43);  label("Second sequential absorption rate constant Ka2, transit1 to transit2 (1/h)")        # Jeong 2023 Table 5: tvKa2 = 2.43 1/h (SE 0.93, CV 38.40%; bootstrap median 2.03, 95% CI 1.66-4.59)
    lka3  <- log(3.34);  label("Third sequential absorption rate constant Ka3, transit2 to central (1/h)")          # Jeong 2023 Table 5: tvKa3 = 3.34 1/h (SE 0.82, CV 24.44%; bootstrap median 3.26, 95% CI 1.97-5.16)
    ltlag <- log(1.56);  label("Absorption lag time Tlag in men (h)")                                               # Jeong 2023 Table 5: tvTlag = 1.56 h (SE 0.21, CV 13.37%; bootstrap median 1.55, 95% CI 1.18-1.98)
    lvc   <- log(10.31); label("Apparent central volume of distribution Vc/F (L)")                                  # Jeong 2023 Table 5: tvVc/F = 10.31 L (SE 1.85, CV 17.90%; bootstrap median 10.06, 95% CI 6.83-14.37)
    lcl   <- log(25.65); label("Apparent central clearance CLc/F (L/h)")                                            # Jeong 2023 Table 5: tvCLc/F = 25.65 L/h (SE 1.36, CV 5.30%; bootstrap median 25.46, 95% CI 23.26-28.03)
    lvp   <- log(11.46); label("Apparent peripheral volume of distribution Vp/F (L)")                               # Jeong 2023 Table 5: tvVp/F = 11.46 L (SE 1.20, CV 10.49%; bootstrap median 11.30, 95% CI 9.43-13.94)
    lq    <- log(5.45);  label("Apparent inter-compartmental clearance CLp/F (L/h)")                                # Jeong 2023 Table 5: tvCLp/F = 5.45 L/h (SE 0.68, CV 12.54%; bootstrap median 5.38, 95% CI 4.02-6.88)

    # ---- Retained covariate effects ----------------------------------------
    # Supplementary Information S7:
    #   Ka3  = tvKa3  * (BSA / median BSA)^dKa3dBSA * exp(eta_Ka3)
    #   Tlag = tvTlag * (1 + dTlagdGender * SEXF)   * exp(eta_Tlag)
    # Note the two forms differ: BSA is a power effect on the median-normalised
    # covariate, while sex is a LINEAR fractional effect. Writing the sex term
    # as exp(0.73 * SEXF) would give a female Tlag of 3.24 h instead of 2.70 h.
    e_bsa_ka3   <- -1.11; label("Power exponent on (BSA / 1.75 m^2) for Ka3 (unitless)")           # Jeong 2023 Table 5: dKa3dBSA = -1.11 (SE 0.48, CV 43.24%; bootstrap median -0.95, 95% CI -1.89 to -0.01)
    e_sexf_tlag <-  0.73; label("Fractional increase in Tlag for women, applied as 1 + e * SEXF (unitless)")  # Jeong 2023 Table 5: dTlagdGender = 0.73 (SE 0.25, CV 34.57%; bootstrap median 0.67, 95% CI 0.32-1.24)

    # ========================================================================
    # PD. Direct-response sigmoid Emax model with baseline for intragastric pH
    # driven by the plasma rabeprazole concentration. Jeong 2023 Table 6 prints
    #   E = E0 + (Emax * C^gamma) / (EC50^gamma + C^gamma)
    # Estimates from Table 6. The paper's sigmoidicity factor gamma is a Hill
    # exponent in a sigmoid Emax function, so it takes the canonical name hill.
    # No IIV and no residual error were reported for any PD parameter
    # (Discussion: the PD model "does not reflect the separate covariates in
    # the drug efficacy variability between individuals").
    # ========================================================================
    le0   <- log(2.50);  label("Baseline intragastric pH without rabeprazole E0 (pH units)")                    # Jeong 2023 Table 6: E0 = 2.50 (SE 0.29, RSE 11.60%, 95% CI 1.94-3.07)
    lemax <- log(4.72);  label("Maximal rabeprazole-induced increase in intragastric pH Emax (pH units)")       # Jeong 2023 Table 6: Emax = 4.72 (SE 0.88, RSE 18.64%, 95% CI 2.98-6.46)
    lec50 <- log(51.58); label("Plasma rabeprazole concentration giving half of Emax, EC50 (ng/mL)")            # Jeong 2023 Table 6: EC50 = 51.58 (SE 5.20, RSE 10.08%, 95% CI 41.29-61.87)
    lhill <- log(5.04);  label("Sigmoidicity (Hill) exponent gamma of the pH response (unitless)")              # Jeong 2023 Table 6: gamma = 5.04 (SE 2.14, RSE 42.46%, 95% CI 0.81-9.27)

    # ========================================================================
    # Inter-individual variability. Supplementary Information S7 specifies an
    # exponential model P_i = P_tv * exp(eta_i) with eta_i ~ N(0, omega^2), and
    # Table 5 reports omega^2 directly (the rows are labelled "omega 2 <par>"),
    # so the variances below are transcribed without conversion.
    #
    # Table 3 model 02-03-02-09 is the selected IIV structure: IIV on Vc/F,
    # CLc/F, Ka1, Ka2, Ka3 and Tlag, and NO IIV on Vp/F or CLp/F (removing the
    # latter two improved the model, d-2LL = -2.51).
    # ========================================================================
    etalvc   ~ 2.04         # Jeong 2023 Table 5: omega^2 Vc/F  = 2.04 (SE 0.61, CV 30.01%; bootstrap median 1.55, 95% CI 0.35-2.75)
    etalcl   ~ 0.14         # Jeong 2023 Table 5: omega^2 CLc/F = 0.14 (SE 0.04, CV 25.16%; bootstrap median 0.14, 95% CI 0.07-0.20)
    etalka1  ~ fixed(0)     # Jeong 2023 Table 5: omega^2 Ka1 reported as 0.00 (SE 0.00; bootstrap median 0.00, 95% CI 0.00-0.00). Retained in the model (Table 3 model 02-03-02-05: removing it worsened -2LL by 24.18) but estimated below the two-decimal reporting precision, so it is carried as an explicitly zero variance rather than invented.
    etalka2  ~ fixed(0)     # Jeong 2023 Table 5: omega^2 Ka2 reported as 0.00 (SE 0.00; bootstrap median 0.00, 95% CI 0.00-0.00). Retained in the model (Table 3 model 02-03-02-06: removing it worsened -2LL by 24.97) but estimated below the two-decimal reporting precision.
    etalka3  ~ 1.74         # Jeong 2023 Table 5: omega^2 Ka3  = 1.74 (SE 0.78, CV 44.67%; bootstrap median 1.82, 95% CI 0.30-3.34)
    etaltlag ~ 0.23         # Jeong 2023 Table 5: omega^2 Tlag = 0.23 (SE 0.07, CV 32.12%; bootstrap median 0.23, 95% CI 0.08-0.38)

    # ========================================================================
    # Residual error. Table 3 selects the "log additive" residual model
    # (02-03-02), which in Phoenix NLME is an additive error on the natural-log
    # scale, log(Cobs) = log(Cpred) + eps -- equivalently a log-normal residual
    # on the linear concentration scale, so it maps to `~ lnorm(...)`.
    #
    # The gastric-pH endpoint has NO reported residual error (Table 6 lists
    # only E0, Emax, gamma and EC50), so it is carried as fixed(0). See the
    # vignette Assumptions and deviations section.
    # ========================================================================
    expSd            <- 0.37;     label("Log-scale residual SD on rabeprazole Cc (log-additive error)")  # Jeong 2023 Table 5: epsilon = 0.37 (SE 0.10, CV 26.38%; bootstrap median 0.35, 95% CI 0.24-0.58)
    addSd_gastric_ph <- fixed(0); label("Additive residual SD on intragastric pH (pH units)")            # Jeong 2023 Table 6 reports no residual-error term for the pharmacodynamic model
  })

  model({
    # ---- 1. Individual PK parameters ---------------------------------------
    # Supplementary Information S7 covariate forms; see the ini() note on why
    # the BSA effect is a power model and the sex effect is linear.
    ka1  <- exp(lka1 + etalka1)
    ka2  <- exp(lka2 + etalka2)
    ka3  <- exp(lka3 + etalka3) * (BSA / 1.75)^e_bsa_ka3
    tlag <- exp(ltlag + etaltlag) * (1 + e_sexf_tlag * SEXF)
    vc   <- exp(lvc + etalvc)
    cl   <- exp(lcl + etalcl)
    vp   <- exp(lvp)
    q    <- exp(lq)

    # ---- 2. Micro-constants ------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 3. ODE system -----------------------------------------------------
    # The enteric-coated tablet does not dissolve in the stomach, so nothing
    # reaches the absorption window until the lag time has elapsed; the three
    # sequential compartments then describe the drawn-out 3-5 h absorption
    # phase seen in Figure 1.
    d/dt(depot)       <- -ka1 * depot
    d/dt(transit1)    <-  ka1 * depot    - ka2 * transit1
    d/dt(transit2)    <-  ka2 * transit1 - ka3 * transit2
    d/dt(central)     <-  ka3 * transit2 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                   k12 * central - k21 * peripheral1

    # ---- 4. Lag time -------------------------------------------------------
    alag(depot) <- tlag

    # ---- 5. Observations ---------------------------------------------------
    # Central amount is in mg and vc in L, so central/vc is mg/L; multiply by
    # 1000 to report ng/mL, the unit in which EC50 was estimated.
    Cc <- (central / vc) * 1000

    # ---- 6. Individual PD parameters and the direct-effect pH response -----
    e0   <- exp(le0)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Jeong 2023 Table 6: E = E0 + (Emax * C^gamma) / (EC50^gamma + C^gamma),
    # with C the plasma rabeprazole concentration. With Cc = 0 this returns the
    # unstimulated baseline pH of E0 = 2.50.
    gastric_ph <- e0 + emax * Cc^hill / (ec50^hill + Cc^hill)

    # ---- 7. Residual error -------------------------------------------------
    Cc         ~ lnorm(expSd)
    gastric_ph ~ add(addSd_gastric_ph)
  })
}
