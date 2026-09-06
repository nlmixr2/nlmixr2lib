Stott_2023_flucytosine <- function() {
  description <- paste(
    "Four-compartment population PK model for oral flucytosine in adults with",
    "HIV-associated cryptococcal meningoencephalitis (Stott 2023): a gut",
    "absorption compartment feeding a central compartment that exchanges with a",
    "CSF/CNS compartment and with a peripheral compartment through asymmetric",
    "first-order transfer rate constants, with first-order elimination from",
    "central. Non-parametric (Pmetrics NPAG) population estimates; no covariate",
    "was retained in the final model.",
    sep = " "
  )
  reference <- paste(
    "Stott KE, Ahmadu A, Kajanga C, Moyo M, Gondwe E, Chimang'anga W,",
    "Chasweka M, Unsworth J, Jimenez-Valverde A, Jagota B, Shah RV,",
    "Lawrence DS, Lalloo DG, Harrison T, Jarvis JN, Hope W, Mwandumba HC.",
    "Population pharmacokinetics and CSF penetration of flucytosine in adults",
    "with HIV-associated cryptococcal meningoencephalitis.",
    "J Antimicrob Chemother. 2023;78(4):1015-1022. doi:10.1093/jac/dkad038.",
    sep = " "
  )
  vignette <- "Stott_2023_flucytosine"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Stott 2023 "Population PK covariate screening" and Results: age, weight,
  # sex, baseline ALT, baseline serum creatinine and baseline Cockcroft-Gault
  # CrCl were screened by bidirectional stepwise multivariate linear regression
  # against the Bayesian posterior parameter estimates. None correlated
  # significantly; a variant model with CL scaled to CrCl gave no improvement in
  # fit, so the base model was retained as the final model. These covariates are
  # therefore documented but never referenced in model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened; not retained. Stott 2023 Results: 'patient weight did not",
        "significantly correlate with volume of distribution'. Median 50 kg",
        "(IQR 47-56 kg). Weight still sets the mg/kg dose but does not scale",
        "any model parameter.",
        sep = " "
      ),
      source_name = "weight"
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened; not retained. Median 36 years (IQR 33-41 years).",
      source_name = "age"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = male",
      notes              = paste(
        "Screened; not retained. The source reports sex as the count of female",
        "patients (24 of 64, 37%), so SEXF is the natural canonical encoding.",
        sep = " "
      ),
      source_name = "sex"
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Screened; not retained. No baseline summary statistic is reported;",
        "the source states only that no patient had a raised ALT at enrolment.",
        sep = " "
      ),
      source_name = "baseline ALT"
    ),
    CREAT = list(
      description = "Baseline serum creatinine",
      units       = "not reported",
      type        = "continuous",
      notes       = paste(
        "Screened; not retained. Stott 2023 Results: 'neither plasma creatinine",
        "level nor CrCl significantly correlated with flucytosine clearance'.",
        "No baseline summary statistic or unit is reported.",
        sep = " "
      ),
      source_name = "baseline serum creatinine"
    ),
    CRCL = list(
      description = "Baseline creatinine clearance (Cockcroft-Gault)",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Screened; not retained. Raw Cockcroft-Gault CrCl in mL/min, NOT",
        "BSA-normalized to mL/min/1.73 m^2. Median 105.15 mL/min (IQR",
        "80.63-123.05 mL/min) at enrolment. The source refit the model with CL",
        "scaled to CrCl because of consistent prior literature, but no",
        "improvement in fit could be demonstrated and the term was dropped; the",
        "Discussion attributes the null result to the narrow CrCl range in this",
        "cohort. CrCl did govern the prescribed dosing interval (q12h for CrCl",
        "20-40 mL/min, q24h for CrCl 10-20 mL/min), which is a dosing-record",
        "property rather than a model parameter.",
        sep = " "
      ),
      source_name = "CrCl"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "flucytosine", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "flucytosine", units = "mg", specimen = "plasma", verified = TRUE),
    csf         = list(analyte = "flucytosine", units = "mg", specimen = "CSF", verified = TRUE),
    peripheral1 = list(analyte = "flucytosine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 64,
    n_studies      = 1,
    age_median     = "36 years (IQR 33-41 years)",
    weight_median  = "50 kg (IQR 47-56 kg)",
    sex_female_pct = 37,
    disease_state  = "HIV-associated cryptococcal meningoencephalitis with advanced immunosuppression (median CD4 39 cells/mm^3, IQR 21-83, n = 57)",
    dose_range     = "flucytosine 100 mg/kg/day orally or by nasogastric tube, given as 25 mg/kg every 6 h; the interval was extended to q12h (CrCl 20-40 mL/min) or q24h (CrCl 10-20 mL/min) on deterioration of renal function",
    renal_function = "Cockcroft-Gault CrCl at enrolment median 105.15 mL/min (IQR 80.63-123.05 mL/min)",
    regions        = "Malawi (Queen Elizabeth Central Hospital, Blantyre)",
    notes          = paste(
      "PK substudy of the Phase III AMBIsome Therapy Induction OptimisatioN",
      "(AMBITION-cm) trial; 64 patients recruited November 2018 to October 2019,",
      "31 in the amphotericin B deoxycholate control arm and 33 in the",
      "single-high-dose liposomal amphotericin B intervention arm. 595 plasma",
      "and 209 CSF flucytosine observations (means of 9.2 and 3.5 samples per",
      "patient); 34 plasma and 82 CSF observations were below the 0.1 mg/L limit",
      "of quantitation and were set to 0.05 mg/L for estimation. Fitted with the",
      "non-parametric adaptive grid algorithm of Pmetrics 2.0.0. Baseline",
      "demographics are in Stott 2023 Results, 'Patients'.",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- Stott 2023 Table 1, "Mean" column. The source
    # fits a non-parametric (NPAG) support-point distribution and reports its
    # Mean, SD, Median and CV% per parameter; the text states that mean
    # estimates fitted the data better than median estimates and were the ones
    # carried into the Bayesian posterior step. The mean set is what reproduces
    # the paper's own reported terminal half-life (14.4 h computed here vs
    # "approximately 14.5 h" reported) and CSF:plasma AUC ratio (0.72 computed
    # here vs a posterior median of 0.69); the median set reproduces neither.
    # All volumes and clearances are APPARENT (per the paper's CL/F, V/F,
    # Vcns/F notation) because flucytosine was given orally / by nasogastric
    # tube and no bioavailability term was estimated.
    lka            <- log(1.77);  label("Absorption rate constant Ka, gut to central (1/h)")                    # Stott 2023 Table 1: Ka mean 1.77 /h (SD 1.81, median 1.01, CV 102.06%)
    lcl            <- log(5.88);  label("Apparent clearance CL/F from central (L/h)")                           # Stott 2023 Table 1: CL/F mean 5.88 L/h (SD 3.35, median 5.31, CV 56.93%)
    lvc            <- log(17.50); label("Apparent central volume V/F (L)")                                      # Stott 2023 Table 1: V/F mean 17.50 L (SD 9.99, median 14.76, CV 57.06%)
    lk_central_csf <- log(15.55); label("Transfer rate constant central to CSF, paper's K23 (1/h)")             # Stott 2023 Table 1: K23 mean 15.55 /h (SD 4.97, median 18.57, CV 31.97%)
    lk_csf_central <- log(9.02);  label("Transfer rate constant CSF to central, paper's K32 (1/h)")             # Stott 2023 Table 1: K32 mean 9.02 /h (SD 5.43, median 6.69, CV 60.21%)
    lk12           <- log(5.68);  label("Transfer rate constant central to peripheral1, paper's K24 (1/h)")     # Stott 2023 Table 1: K24 mean 5.68 /h (SD 7.60, median 1.09, CV 133.77%)
    lk21           <- log(1.38);  label("Transfer rate constant peripheral1 to central, paper's K42 (1/h)")     # Stott 2023 Table 1: K42 mean 1.38 /h (SD 2.33, median 0.24, CV 168.39%)
    lvcsf          <- log(41.73); label("Apparent CSF/CNS compartment volume Vcns/F (L)")                       # Stott 2023 Table 1: Vcns/F mean 41.73 L (SD 13.66, median 44.55, CV 32.74%)

    # IIV. The source is a non-parametric fit, so it publishes a support-point
    # distribution rather than an OMEGA matrix; Table 1 reports only a per-
    # parameter SD and CV%. The authors themselves approximated that
    # distribution by a log-normal for their confirmatory ADAPT 5 simulations
    # ("values for the system parameters were randomly selected from a log
    # normal distribution, with mean and the diagonal of the covariance matrix
    # derived from the original model fit in Pmetrics"), and this encoding
    # follows them: each variance is omega^2 = log(1 + CV^2) using the
    # four-significant-figure CV% column, and the matrix is DIAGONAL because the
    # source used only the diagonal of the covariance matrix. Off-diagonal
    # covariances between support-point parameters are not reported.
    etalka            ~ 0.71375   # Stott 2023 Table 1: CV(Ka)     = 102.06% -> log(1 + 1.0206^2)
    etalcl            ~ 0.28073   # Stott 2023 Table 1: CV(CL/F)   =  56.93% -> log(1 + 0.5693^2)
    etalvc            ~ 0.28185   # Stott 2023 Table 1: CV(V/F)    =  57.06% -> log(1 + 0.5706^2)
    etalk_central_csf ~ 0.09732   # Stott 2023 Table 1: CV(K23)    =  31.97% -> log(1 + 0.3197^2)
    etalk_csf_central ~ 0.30934   # Stott 2023 Table 1: CV(K32)    =  60.21% -> log(1 + 0.6021^2)
    etalk12           ~ 1.02584   # Stott 2023 Table 1: CV(K24)    = 133.77% -> log(1 + 1.3377^2)
    etalk21           ~ 1.34430   # Stott 2023 Table 1: CV(K42)    = 168.39% -> log(1 + 1.6839^2)
    etalvcsf          ~ 0.10183   # Stott 2023 Table 1: CV(Vcns/F) =  32.74% -> log(1 + 0.3274^2)

    # Residual error. Stott 2023 Methods states that "for each model, additive
    # and multiplicative error models were assessed and optimized", but neither
    # the winning form nor any coefficient of the Pmetrics error polynomial is
    # reported anywhere in the paper, and there is no supplement. Both forms are
    # therefore carried for both outputs and fixed to zero rather than invented;
    # the only precision figures the paper gives are LC-MS-MS assay quality-
    # control CVs (interday 5.82-8.58%, intraday 6.01-7.01%), which measure the
    # assay and not the model's residual variability.
    propSd      <- fixed(0); label("Proportional residual SD on plasma Cc (fraction; 0 -- not reported in the source)")    # Stott 2023 Methods, "Population PK modelling": error-model form assessed but no estimate published
    addSd       <- fixed(0); label("Additive residual SD on plasma Cc (mg/L; 0 -- not reported in the source)")            # Stott 2023 Methods, "Population PK modelling": error-model form assessed but no estimate published
    propSd_Ccsf <- fixed(0); label("Proportional residual SD on CSF Ccsf (fraction; 0 -- not reported in the source)")     # Stott 2023 Methods, "Population PK modelling": error-model form assessed but no estimate published
    addSd_Ccsf  <- fixed(0); label("Additive residual SD on CSF Ccsf (mg/L; 0 -- not reported in the source)")             # Stott 2023 Methods, "Population PK modelling": error-model form assessed but no estimate published
  })

  model({
    # Individual parameters. No covariate enters any of them: the final model
    # is the base model (Stott 2023 Results, "Population PK analysis").
    ka            <- exp(lka            + etalka)
    cl            <- exp(lcl            + etalcl)
    vc            <- exp(lvc            + etalvc)
    vcsf          <- exp(lvcsf          + etalvcsf)
    k_central_csf <- exp(lk_central_csf + etalk_central_csf)
    k_csf_central <- exp(lk_csf_central + etalk_csf_central)
    k12           <- exp(lk12           + etalk12)
    k21           <- exp(lk21           + etalk21)
    kel           <- cl / vc

    # Stott 2023 Equations 1-4. The source numbers its states 1 = gut, 2 =
    # circulation, 3 = CNS, 4 = peripheral, so its K23 / K32 are the central
    # <-> CSF pair and its K24 / K42 are the central <-> peripheral pair. The
    # CSF and peripheral legs are each parameterised by an asymmetric pair of
    # first-order rate constants rather than by a single inter-compartmental
    # clearance; K23 * V (272 L/h) does not equal K32 * Vcns (376 L/h), so the
    # system genuinely cannot be re-expressed with a symmetric q.
    d/dt(depot)       <- -ka * depot                                                                       # Equation 1: dX(1)/dt = -Ka * X(1)
    d/dt(central)     <-  ka * depot - (k_central_csf + k12 + kel) * central +
                          k_csf_central * csf + k21 * peripheral1                                          # Equation 2: dX(2)/dt = Ka*X(1) - (K23 + K24 + CL/V)*X(2) + K32*X(3) + K42*X(4)
    d/dt(csf)         <-  k_central_csf * central - k_csf_central * csf                                    # Equation 3: dX(3)/dt = K23*X(2) - K32*X(3)
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1                                                # Equation 4: dX(4)/dt = K24*X(2) - K42*X(4)

    Cc   <- central / vc                                                                                   # Equation 5: Y(1) = X(2)/V
    Ccsf <- csf / vcsf                                                                                     # Equation 6: Y(2) = X(3)/Vcns

    Cc   ~ add(addSd)      + prop(propSd)
    Ccsf ~ add(addSd_Ccsf) + prop(propSd_Ccsf)
  })
}
