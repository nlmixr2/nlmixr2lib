Randell_2024_metronidazole <- function() {
  description <- "One-compartment intravenous population PK model for metronidazole in critically ill preterm and term infants (Randell 2024, 'optimized model'). Clearance scales linearly with body weight and matures with postmenstrual age through a sigmoidal Emax (Hill) function (TM50 25.6 weeks, Hill 15.7); central volume scales with body weight through an estimated power exponent (0.763). The model was developed by externally validating and then optimizing the PTN_METRO model of Cohen-Wolkowiez 2013 using estimated plasma concentrations derived from opportunistically collected dried blood spots."
  reference <- paste(
    "Randell RL, Balevic SJ, Greenberg RG, et al. Opportunistic dried blood spot sampling validates and optimizes a pediatric population pharmacokinetic model of metronidazole.",
    "Antimicrob Agents Chemother. 2024;68(4):e01533-23. doi:10.1128/aac.01533-23.",
    "Erratum: Antimicrob Agents Chemother. 2025;69(9):e00972-25. doi:10.1128/aac.00972-25 (Table 1 Hill parameter: '7' should read '15.7')."
  )
  vignette <- "Randell_2024_metronidazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central = list(analyte = "metronidazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Combined analysis cohort mean (SD) 1.8 (1.0) kg (Randell 2024 Table S1). Enters CL linearly (fixed exponent 1, the PTN_METRO base structural relationship) and V through an estimated power exponent theta_WT-V = 0.763. No reference weight is used: both relationships are written on absolute body weight in kg, so theta_CL is a per-kg clearance and theta_V is the volume at WT = 1 kg.",
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age plus chronologic age)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. The paper reports PMA in weeks (combined cohort mean (SD) 33.2 (5.3) weeks; supported range 23-48 weeks). Drives the sigmoidal Emax (Hill) maturation function on CL with TM50 = 25.6 weeks. Inside model() the canonical PAGE (months) is converted back to weeks via PMA_wk = PAGE * 4.35 so the published TM50 in weeks applies directly.",
      source_name        = "PMA"
    )
  )

  covariatesDataExcluded <- list(
    BGA = list(
      description = "Gestational age at birth",
      units       = "weeks",
      type        = "continuous",
      notes       = "Screened during the covariate re-evaluation for the optimized model (Randell 2024 Methods 'Optimization'; selection process in Table S3) but not retained: the postmenstrual-age maturation function on CL captured the developmental effect. Combined cohort mean (SD) 30.6 (5.2) weeks, supported range 22.7-41.0 weeks."
    ),
    PNA = list(
      description = "Postnatal age",
      units       = "days",
      type        = "continuous",
      notes       = "Screened but not retained in the optimized model (Randell 2024 Methods 'Optimization', Table S3). Combined cohort mean (SD) 18.6 (17.3) days, supported range 0-80 days."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained (Randell 2024 Methods 'Optimization', Table S3). Combined cohort 81/146 (55%) male, i.e., 45% female (Table S1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened but not retained (Randell 2024 Methods 'Optimization', Table S3). Combined cohort mean (SD) 0.6 (0.4) mg/dL (Table S1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Screened but not retained (Randell 2024 Methods 'Optimization', Table S3). Available in 101/146 infants; mean (SD) 2.4 (0.6) g/dL (Table S1)."
    ),
    CONMED_CYP3A_INDUCER = list(
      description = "Concomitant administration of a CYP3A inducer",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained (Randell 2024 Methods 'Optimization', Table S3). Metronidazole is metabolized in part by CYP3A, so induction and inhibition were biologically plausible covariates on CL; neither met the forward-inclusion criterion."
    ),
    CONMED_CYP3A_INHIBITOR = list(
      description = "Concomitant administration of a CYP3A inhibitor",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained (Randell 2024 Methods 'Optimization', Table S3)."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 146L,
    n_studies        = 2L,
    age_range        = "22.7-41.0 weeks gestational age at birth; 0-80 days postnatal age; 23-48 weeks postmenstrual age",
    age_median       = "Means (SD) in the combined cohort: gestational age 30.6 (5.2) weeks; postnatal age 18.6 (17.3) days; postmenstrual age 33.2 (5.3) weeks",
    weight_range     = "Combined cohort mean (SD) 1.8 (1.0) kg",
    weight_median    = "1.8 kg (mean)",
    sex_female_pct   = 45,
    race_ethnicity   = c(White = 60, Other = 40),
    disease_state    = "Critically ill preterm and term infants with suspected or confirmed complicated intra-abdominal infection (including necrotizing enterocolitis) receiving intravenous metronidazole in the neonatal intensive care unit.",
    dose_range       = "SCAMP trial: metronidazole given as a 30-minute intravenous infusion, 15 mg/kg loading dose followed by a 7.5 mg/kg maintenance dose at 24 hours, then 7.5 mg/kg at postmenstrual-age-defined intervals (every 12 h for PMA 23-<34 weeks, every 8 h for PMA 34-40 weeks, every 6 h for PMA >40 weeks). PTN_METRO doses were those of the source study (Cohen-Wolkowiez 2013).",
    regions          = "United States and Canada (multicenter Pediatric Trials Network sites)",
    renal_function   = "Serum creatinine mean (SD) 0.6 (0.4) mg/dL in the combined cohort",
    hepatic_function = "ALT mean (SD) 21.0 (20.5) U/L (n = 75), AST 51.0 (47.6) U/L (n = 73), total bilirubin 1.0 (1.0) mg/dL (n = 85); measured only in SCAMP",
    ethnicity        = "30/146 (21%) Hispanic ethnicity",
    notes            = "The optimized model was fit to a combined data set of plasma metronidazole concentrations from the Metronidazole Pharmacokinetics in Premature Infants study (PTN_METRO, NCT01222585) and plasma concentrations estimated from dried blood spots (ePlasma) collected in the Antibiotic Safety in Infants with Complicated Intra-Abdominal Infections trial (SCAMP, NCT01994993). Demographics per Randell 2024 Table S1 (SCAMP n = 122, PTN_METRO n = 24, combined n = 146). DBS were converted to ePlasma with the comparability regression Cplasma = 1.11 * CDBS + 253 ng/mL (R-squared 0.86, 42 paired samples from 21 PTN_METRO infants). NONMEM 7.4 FOCE-I, ADVAN1 TRANS2. Eta shrinkage 8% (CL) and 29% (V); epsilon shrinkage 17%. Bootstrap 1000 replicates, 99.7% converged."
  )

  ini({
    # Structural parameters (Randell 2024 Table 1, 'optimized model').
    # CL = theta_CL * WT * [PMA^Hill / (TM50^Hill + PMA^Hill)]
    # V  = theta_V  * WT^theta_WT-V
    lcl <- log(0.036);  label("Weight-normalized clearance at full maturation (theta_CL, L/h/kg)")  # Randell 2024 Table 1: theta_CL = 0.036 L/kg/h (RSE 4%; bootstrap 2.5th-97.5th 0.033-0.039)
    lvc <- log(0.853);  label("Volume of distribution at WT = 1 kg (theta_V, L)")                   # Randell 2024 Table 1: theta_V = 0.853 L/kg (RSE 3%; bootstrap 0.802-0.921)

    # Body-weight exponents. The exponent on CL is the fixed linear relationship
    # carried over from the PTN_METRO base structural model (CLi = CL * WT^1);
    # the exponent on V was estimated during the optimization.
    e_wt_cl <- fixed(1);  label("Exponent of body weight on CL (unitless; linear)")  # Randell 2024 Methods 'Optimization': CLi = CL * WT^1 (fixed exponent) and Table 1 structural model 'CL = theta_CL * WT'
    e_wt_vc <- 0.763;     label("Exponent of body weight on V (unitless)")                 # Randell 2024 Table 1: theta_WT-V = 0.763 (RSE 8%; bootstrap 0.645-0.876)

    # Sigmoidal Emax (Hill) maturation of CL on postmenstrual age, in weeks.
    pma50_cl <- 25.6; label("Postmenstrual age at 50% of mature CL (TM50, weeks)")           # Randell 2024 Table 1: TM50 = 25.6 weeks (RSE 2%; bootstrap 24.4-26.5)
    hill_cl  <- 15.7; label("Hill coefficient of the PMA maturation function on CL (unitless)")  # Randell 2024 Table 1 as corrected by the 2025 erratum (doi:10.1128/aac.00972-25): Hill = 15.7, not the 7 printed in Table 1 row 7 (RSE 37%; bootstrap 8.4-41.0, which excludes 7)

    # Interindividual variability (Randell 2024 Table 1, reported as %CV).
    # omega^2 = log(CV^2 + 1) for a log-normal eta on a log-transformed parameter.
    etalcl ~ 0.11933  # 35.6 %CV -> log(1 + 0.356^2) = 0.11933; Randell 2024 Table 1 (RSE 21%; bootstrap 28.4-43.4 %CV; eta shrinkage 8%)
    etalvc ~ 0.06252  # 25.4 %CV -> log(1 + 0.254^2) = 0.06252; Randell 2024 Table 1 (RSE 31%; bootstrap 16.5-31.6 %CV; eta shrinkage 29%)

    # Residual error.
    propSd <- 0.191; label("Proportional residual error (fraction)")  # Randell 2024 Table 1: proportional error = 19.1% (RSE 13%; bootstrap 16.5-21.1%); epsilon shrinkage 17%
  })
  model({
    # Convert canonical PAGE (months) back to the source-paper PMA (weeks) so
    # the published TM50 in weeks applies directly.
    pma_wk <- PAGE * 4.35

    # Sigmoidal Emax maturation of clearance on postmenstrual age.
    fmat_cl <- pma_wk^hill_cl / (pma50_cl^hill_cl + pma_wk^hill_cl)

    # Individual PK parameters.
    cl <- exp(lcl + etalcl) * WT^e_wt_cl * fmat_cl
    vc <- exp(lvc + etalvc) * WT^e_wt_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, volume in L -> central/vc has units mg/L (1 mg/L = 1000 ng/mL,
    # the units in which the source paper reports concentrations).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
