Chung_2026_vancomycin <- function() {
  description <- "One-compartment intravenous population PK model for vancomycin in preterm and term neonates (Chung 2026 Equations 1-2; the model was developed and internally evaluated in Chung 2023 and is reproduced in full by Chung 2026, which externally validated it). Clearance scales linearly with body weight normalized to 70 kg, matures with postmenstrual age through a sigmoidal Emax (Hill) function (TM50 47.7 weeks, Hill 0.739), and falls as a power function of serum creatinine (exponent -0.653, reference 34 umol/L); volume of distribution scales linearly with body weight normalized to 70 kg and carries no other covariate. In a head-to-head external validation against 32 other published neonatal vancomycin models in 366 neonates it was the only one of the 33 to meet all five predefined predictive-performance criteria. Chung 2026 publishes only the typical-value (a priori) clearance and volume equations, so between-subject and residual variability are encoded as zero."
  reference <- paste(
    "Chung E, Tabbara N, Seto W, Shah V. Vancomycin therapeutic drug monitoring, clinical outcomes and population pharmacokinetic model evaluation in neonates.",
    "Children (Basel). 2026;13(5):649. doi:10.3390/children13050649.",
    "Model originally developed in: Chung E, Seto W. Using population pharmacokinetics to optimize initial vancomycin dosing guidelines for neonates to treat sepsis caused by coagulase-negative staphylococcus.",
    "Pharmacotherapy. 2023;43(12):1262-1276. doi:10.1002/phar.2865."
  )
  vignette <- "Chung_2026_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Vancomycin is given as a 1 h intravenous infusion
  # (Chung 2026 Methods, "External predictive performance evaluation"), so the
  # dose enters `central` directly and there is no depot state. Chung 2026
  # Methods, "Vancomycin dosing and monitoring practice" states that
  # "vancomycin serum concentrations were measured using an immunoturbidimetric
  # assay on the Roche cobas platform", so the specimen is serum.
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at the start of vancomycin therapy",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Chung 2026 Equations (1) and (2) write both CL and V on (WT/70)^1, i.e. linear scaling normalized to a 70 kg reference, so theta_CL = 13.9 L/h and theta_V = 65.5 L are the values a 70 kg reference subject would take. The exponent of exactly 1 is printed in both equations and carries no uncertainty, so it is encoded with fixed(). Model-development cohort (Chung 2026 Table 1, 'Previous cohort [19] used for popPK model development') mean (SD) 1.88 (0.99) kg; external-validation cohort 1.06 (0.60) kg. The extraction is calibrated in the neonatal weight range only: nothing in either cohort approaches 70 kg, so the 70 kg normalization is a scaling convention rather than a supported extrapolation.",
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age plus postnatal age)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Drives the sigmoidal Emax (Hill) maturation function on CL in Chung 2026 Equation (1), PMA^0.739 / (PMA^0.739 + 47.7^0.739). The source paper writes PMA in WEEKS and its TM50 of 47.7 weeks is the Rhodin 2009 renal-maturation half-time; the canonical PAGE column is in months, so model() converts back with pma_wk = PAGE * 4.35 and the published TM50 in weeks applies unchanged. Note that Chung 2023 kept the Rhodin TM50 of 47.7 weeks but estimated a much shallower Hill of 0.739 (Rhodin uses 3.4), which flattens the maturation profile across the neonatal range. Model-development cohort mean (SD) 34.3 (5.03) weeks; external-validation cohort 28.9 (3.81) weeks. Both cohorts excluded PMA of 44 weeks and above (Chung 2026 Methods, 'Study design and population').",
      source_name        = "PMA"
    ),
    CREAT = list(
      description        = "Serum creatinine, measured by the enzymatic creatinine method on the Roche cobas platform",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Chung 2026 Methods, 'PopPK model identification' states explicitly that SCr in Equation (1) is in umol/L, and the model enters it as (SCr/34)^-0.653. The reference of 34 umol/L (about 0.38 mg/dL) is close to the Chung 2023 model-development cohort median of 29 umol/L (IQR 21-43, Chung 2026 Table 1); the external-validation cohort is markedly more renally impaired at a median of 56 umol/L (IQR 39.0-79.5). The negative exponent means a higher serum creatinine (poorer glomerular filtration) lowers clearance. About 30% of the validation cohort lacked a documented baseline SCr and those neonates could not be used in the popPK evaluation (Chung 2026 Limitations).",
      source_name        = "SCr"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 648L,
    n_studies        = 1L,
    n_sites          = 1L,
    age_range        = "Postmenstrual age below 44 weeks and postnatal age below 28 days (eligibility); model-development cohort mean (SD) postmenstrual age 34.3 (5.03) weeks and gestational age 29.9 (5.52) weeks",
    age_median       = "Postnatal age median 21.3 days (IQR 9.33-44.3) in the model-development cohort",
    weight_range     = "Model-development cohort mean (SD) 1.88 (0.99) kg; birth weight 1.49 (1.04) kg",
    weight_median    = "1.88 kg (mean)",
    sex_female_pct   = 41.2,
    disease_state    = "Neonates in a tertiary (Level IV) neonatal intensive care unit receiving intravenous vancomycin, predominantly for suspected or confirmed late-onset sepsis caused by coagulase-negative staphylococci. Comorbidity in the model-development cohort: intra-abdominal infection 51.2%, patent ductus arteriosus 24.7%, other congenital heart disease 12.8%.",
    dose_range       = "Weight- and postmenstrual-age-based institutional dosing with therapeutic drug monitoring. The external-validation site's standard initial regimen was 10 mg/kg per dose as a 1 h intravenous infusion, every 12 h for postmenstrual age below 30 weeks and every 8 h for postmenstrual age of 30 weeks and above (Chung 2026 Methods, 'Vancomycin dosing and monitoring practice').",
    regions          = "Canada (Toronto, Ontario)",
    renal_function   = "Serum creatinine median 29.0 umol/L (IQR 21.0-43.0) in the model-development cohort; last-24-h urine output mean (SD) 3.83 (1.68) mL/kg/h",
    co_medication    = "Concurrent nephrotoxic drugs in the model-development cohort: aminoglycoside (gentamicin or tobramycin) 52.2%, furosemide 30.6%, acyclovir 4.3%, amphotericin B 3.1%, hydrochlorothiazide-spironolactone 1.5%",
    notes            = "The parameter values reproduced here are those of the Chung 2023 model, developed in neonates admitted to the SickKids Level IV NICU in Toronto and reproduced verbatim as Equations (1) and (2) of Chung 2026. Chung 2026 Table 1 gives the development cohort as N = 648; the Chung 2026 Introduction instead says the model was developed 'based on 442 neonates', a discrepancy that cannot be resolved from Chung 2026 alone (see the vignette Errata). The externally-validating cohort of Chung 2026 is a separate, more preterm and more renally-impaired population: 366 neonates with 661 vancomycin concentrations (median 1 per neonate) from the Mount Sinai Hospital Level III NICU in Toronto, treated between 1 October 2016 and 31 December 2021, with mean (SD) postmenstrual age 28.9 (3.81) weeks, weight 1.06 (0.60) kg and median serum creatinine 56.0 umol/L. Against that cohort the model met all five predefined a priori criteria (mean error -0.46 mg/L, relative mean error 12.2%, relative median error 4.2%, RMSE 8.00, 49.3% of predictions within 30% of observed; Chung 2026 Table 5) and was the best of the 33 models compared. Concentrations were assayed by immunoturbidimetry on the Roche cobas platform with a lower limit of quantification of 4 mg/L. Model fitting for Chung 2023 was performed in Phoenix WinNonlin/NLME."
  )

  ini({
    # ===== Structural PK: Chung 2026 Equation (1) and Equation (2) =====
    #   CL_i (L/h) = 13.9 * (WT/70)^1 * PMA^0.739 / (PMA^0.739 + 47.7^0.739) * (SCr/34)^-0.653
    #   V_i  (L)   = 65.5 * (WT/70)^1
    # with WT in kg, PMA in weeks and SCr in umol/L. The reference subject is
    # therefore 70 kg, fully mature, with SCr = 34 umol/L. No standard errors,
    # RSEs or confidence intervals are reported for any of these values in
    # Chung 2026; Chung 2023 holds the estimation detail.
    #
    # UNITS OF 13.9 AND 65.5 ARE INFERRED, NOT PRINTED. Chung 2026 states the
    # units of the three covariates but not of CL and V. L/h and L are the only
    # reading consistent with the paper's own numbers, and two independent
    # cross-checks confirm it (both reproduced in the vignette): (a) the
    # model-development cohort typical subject gives a half-life of 6.7 h
    # against the median of 7.3 h that Chung 2026 Methods quotes from
    # Chung 2023, and (b) the validation cohort's standard 10 mg/kg regimen
    # gives a typical pre-4th-dose trough near the observed cohort median of
    # 8.9 mg/L. mL/min would put a 70 kg subject at 0.83 L/h, roughly an order
    # of magnitude below any published adult vancomycin clearance.
    lcl <- log(13.9); label("Clearance at WT = 70 kg, full PMA maturation and SCr = 34 umol/L (CL, L/h)")  # Chung 2026 Methods 'PopPK model identification', Equation (1): coefficient 13.9
    lvc <- log(65.5); label("Volume of distribution at WT = 70 kg (V, L)")                                  # Chung 2026 Methods 'PopPK model identification', Equation (2): coefficient 65.5

    # Body-weight exponents. Both equations print an exponent of exactly 1 with
    # no uncertainty, i.e. linear (not 0.75-allometric) weight scaling that was
    # imposed rather than estimated, so both are wrapped in fixed().
    e_wt_cl <- fixed(1); label("Exponent of (WT / 70 kg) on CL (unitless; linear)")  # Chung 2026 Equation (1): (WT/70)^1
    e_wt_vc <- fixed(1); label("Exponent of (WT / 70 kg) on V (unitless; linear)")   # Chung 2026 Equation (2): (WT/70)^1

    # Sigmoidal Emax (Hill) maturation of CL on postmenstrual age, in weeks.
    # TM50 = 47.7 weeks is the Rhodin 2009 renal-maturation half-time, but the
    # Hill coefficient of 0.739 is far shallower than Rhodin's 3.4.
    pma50_cl <- 47.7;  label("Postmenstrual age at 50% of mature CL (TM50, weeks)")                   # Chung 2026 Equation (1): 47.7 in the PMA^0.739 / (PMA^0.739 + 47.7^0.739) term
    hill_cl  <- 0.739; label("Hill coefficient of the PMA maturation function on CL (unitless)")      # Chung 2026 Equation (1): exponent 0.739 on both PMA and 47.7

    # Renal-function effect on CL. Negative, so a higher serum creatinine
    # (poorer glomerular filtration) lowers clearance.
    e_creat_cl <- -0.653; label("Exponent of (SCr / 34 umol/L) on CL (unitless)")  # Chung 2026 Equation (1): (SCr/34)^-0.653

    # ===== Variability =====
    # Chung 2026 reproduces only the typical-value (a priori) equations, which
    # is all its external-validation workflow needed: population predictions
    # were computed per observed sampling time from each neonate's covariates
    # (Methods, 'External predictive performance evaluation'). Chung 2023 must
    # carry between-subject and residual variance terms -- Chung 2026 reports a
    # 1000-simulation visual predictive check (Supplementary Figure S2) and
    # individual predictions (Figure 2b), neither of which is computable
    # without them -- but no omega and no sigma is printed anywhere in
    # Chung 2026 or its supplementary file list, and Chung 2023 is paywalled
    # and not on disk. The magnitudes are therefore encoded as zero rather
    # than invented; simulations from this model are typical-value profiles.
    # See the vignette Errata.
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- not reported in the source)")  # Chung 2026: no residual-error estimate published
    addSd  <- fixed(0); label("Additive residual SD (mg/L; 0 -- not reported in the source)")          # Chung 2026: no residual-error estimate published
  })
  model({
    # ----- Derived covariate terms -----
    # Convert canonical PAGE (months) back to the source-paper PMA (weeks) so
    # the published TM50 of 47.7 weeks applies directly.
    pma_wk <- PAGE * 4.35

    # Sigmoidal Emax maturation of clearance on postmenstrual age
    # (Chung 2026 Equation 1, middle term):
    #   F_PMA = PMA^Hill / (PMA^Hill + TM50^Hill)
    fmat_cl <- pma_wk^hill_cl / (pma50_cl^hill_cl + pma_wk^hill_cl)

    # Power effect of serum creatinine on clearance (Chung 2026 Equation 1,
    # right-hand term), with the published reference of 34 umol/L.
    creat_factor <- (CREAT / 34)^e_creat_cl

    # ----- Individual PK parameters -----
    # The reference body weight of 70 kg is the normalization printed in both
    # Chung 2026 Equation (1) and Equation (2).
    cl <- exp(lcl) * (WT / 70)^e_wt_cl * fmat_cl * creat_factor
    vc <- exp(lvc) * (WT / 70)^e_wt_vc

    # ----- Micro-constants -----
    kel <- cl / vc

    # ----- ODE system -----
    # Vancomycin is given as a 1 h intravenous infusion into the central
    # compartment; there is no depot state.
    d/dt(central) <- -kel * central

    # ----- Output -----
    # Dose in mg, vc in L, so central/vc is mg/L -- the units Chung 2026
    # reports vancomycin serum concentrations in throughout.
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
