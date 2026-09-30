Wang_2020_polymyxinB <- function() {
  description <- "Two-compartment intravenous population PK model for polymyxin B (sum of polymyxin B1 and B2) in Chinese adults with multidrug-resistant Gram-negative bacterial infections, sampled at steady state on day 4 of therapy (Wang 2020). Cockcroft-Gault creatinine clearance is the sole retained covariate, entering clearance as a power term normalized to 105.9 mL/min with exponent 0.362. Correlated inter-individual variability on V, CL and V2 plus independent variability on Q. Proportional residual error."
  reference <- paste(
    "Wang P, Zhang Q, Zhu Z, Feng M, Sun T, Yang J, Zhang X.",
    "Population Pharmacokinetics and Limited Sampling Strategy for",
    "Therapeutic Drug Monitoring of Polymyxin B in Chinese Patients With",
    "Multidrug-Resistant Gram-Negative Bacterial Infections.",
    "Front Pharmacol. 2020;11:829.",
    "doi:10.3389/fphar.2020.00829. PMCID PMC7289991.",
    sep = " "
  )
  vignette <- "Wang_2020_polymyxinB"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # What each ODE state holds, in what amount units, in what biological
  # matrix. Verified against Wang 2020: plasma polymyxin B1 and B2 were
  # quantified by LC-MS/MS and summed to a total polymyxin B concentration
  # (Methods, 'Quantification of Polymyxin B Concentrations'), and the final
  # model is two-compartment (Results, Equ. 4-7: V central, V2 peripheral).
  compartmentData <- list(
    central = list(analyte = "polymyxinB", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "polymyxinB", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance calculated with the Cockcroft-Gault equation using body weight. NOT body-surface-area normalised: the source reports raw Cockcroft-Gault mL/min.",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Wang 2020 Methods 'Patients Demographics': 'CrCL was calculated using the Cockcroft-Gault equation with body weight (Janmahasatian et al., 2005).' Applied as a power covariate on CL, Cl = 1.786 * (CrCl / 105.9)^0.362 * exp(eta_Cl) (Results, Equ. 6). The Results text states that '105.9 ml/min was the median of CrCL', but Table 1 reports the cohort CrCL as median 89.3 (range 15.6-315.2) mL/min; the normalizing constant of the printed equation (105.9) is used here, and the discrepancy is recorded in the vignette. Stored under the canonical CRCL column with raw mL/min units, following the raw Cockcroft-Gault precedents listed in the CRCL register entry (e.g. Shu_2024_posaconazole.R, Delattre_2010_amikacin.R). Time-fixed per subject in this analysis (all samples collected on day 4).",
      source_name = "CrCL"
    )
  )

  covariatesDataExcluded <- list(
    SCREENED_NOT_RETAINED = list(
      description = "Candidate covariates screened by stepwise forward selection / backward elimination and not retained",
      units = "(various)",
      type = "continuous",
      notes = "Wang 2020 Methods 'Population Pharmacokinetics Analysis': age, sex, body weight, ALT, AST, GGT, ALP, urea nitrogen, serum creatinine, serum uric acid, serum proteins, serum albumin, total bilirubin, direct bilirubin and CrCL were evaluated. Only CrCL on CL was retained (dOFV = 8.26, P < 0.01); Results: 'Age, gender, and other laboratory data had no significant effect on population PK parameters.' No per-covariate dOFV is reported for the rejected candidates. The Discussion attributes the absence of a body-weight effect to the narrow weight range (45-98 kg)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 46L,
    n_studies = 1L,
    age_range = "18-94 years (median 46)",
    weight_range = "45-98 kg (median 70)",
    sex_female_pct = 15.22,
    race_ethnicity = "Chinese (single-centre cohort in Zhengzhou, China)",
    disease_state = "Adults (>= 18 years) with documented multidrug-resistant Gram-negative bacterial infections receiving intravenous polymyxin B sulfate for >= 72 hours; patients on renal replacement therapy were excluded. Sputum was the most common primary infection site. Pathogens: Klebsiella pneumoniae (20), Acinetobacter baumannii (19), Pseudomonas aeruginosa (7), Escherichia coli (4), others (2); 6 patients carried two MDR Gram-negative organisms.",
    dose_range = "Intravenous polymyxin B sulfate; maintenance 50-100 mg twice daily with clinical loading doses of 100-150 mg. Daily dose 100 mg (n = 25), 150 mg (n = 15) or 200 mg (n = 8); daily dose per body weight median 1.91 (range 1.18-3.33) mg/kg. Infusion duration 0.5 h (n = 2), 1 h (n = 37) or 2 h (n = 7).",
    regions = "China (First Affiliated Hospital of Zhengzhou University), April 2018 to November 2019",
    renal_function = "Cockcroft-Gault CrCL median 89.3 (range 15.6-315.2) mL/min (Table 1); serum creatinine median 73.0 (21.0-387.0) umol/L. Renal replacement therapy was an exclusion criterion.",
    n_observations = "331 plasma polymyxin B concentrations (sum of B1 and B2) from 46 patients. On day 4 of therapy one pre-dose sample (C0h) plus five to seven post-dose samples (mainly 0.5, 1, 1.5, 2, 4, 6 and 8 h) were collected within one dosing interval.",
    notes = "Single-centre prospective study. Phoenix NLME 7.0 with first-order conditional estimation. Model evaluation by goodness-of-fit plots, prediction-corrected VPC (200 replicates) and a 1000-sample bootstrap (Table 2). Plasma quantified by LC-MS/MS (polymyxin B1 0.2-10 ug/mL, B2 0.05-2.5 ug/mL calibration ranges). The paper also develops Bayesian and multiple-linear-regression limited sampling strategies for AUC0-12h (Tables 4-5); those are TDM estimation tools, not part of the structural model."
  )

  ini({
    # Structural parameters -- Wang 2020 Table 2 (final model) and the
    # final-model equations Equ. 4-7 in the Results. The typical clearance is
    # the value at CrCL = 105.9 mL/min.
    lvc <- log(6.218)  ; label("Central volume of distribution V (L)")                     # Table 2: tvV = 6.218 L (SE 0.83, CV 13.33%; bootstrap median 5.960); Equ. 4
    lvp <- log(11.922) ; label("Peripheral volume of distribution V2 (L)")                 # Table 2: tvV2 = 11.922 L (SE 1.74, CV 14.62%; bootstrap median 12.073); Equ. 5
    lcl <- log(1.786)  ; label("Clearance CL (L/h) at CRCL = 105.9 mL/min")                # Table 2: tvCl = 1.786 L/h (SE 0.12, CV 6.75%; bootstrap median 1.771); Equ. 6
    lq  <- log(13.518) ; label("Intercompartmental clearance Q (L/h)")                     # Table 2: tvQ = 13.518 L/h (SE 3.35, CV 24.82%; bootstrap median 14.427); Equ. 7

    # Covariate effect on CL -- Wang 2020 Equ. 6:
    #   Cl (L/h) = 1.786 * (CrCl / 105.9)^0.362 * exp(eta_Cl)
    # Estimated (Table 2 gives an SE and a bootstrap 95% CI), so not held
    # constant in ini().
    e_crcl_cl <- 0.362 ; label("Power exponent on (CRCL / 105.9 mL/min) for CL (unitless)") # Table 2: dCldCrCL = 0.362 (SE 0.09; bootstrap median 0.357, 95% CI 0.196-0.513); Equ. 6

    # Inter-individual variability -- Wang 2020 Methods Equ. 1,
    # P_i = theta * exp(eta_i), eta ~ N(0, omega^2). Table 2 reports the
    # 'w2' rows as variances (footnote: 'variance of inter-individual
    # variability'). V, Cl and V2 share a block (Results: 'a strong
    # correlation between Cl, Cl2, and V was observed, and then incorporated
    # it into non-diagonal random effects', dOFV = 30.60) with correlations
    # V-Cl 0.713, V-V2 0.667, Cl-V2 0.571; covariance = r * sqrt(var_i * var_j).
    # Q is not in the block. Covariances: V-Cl 0.713*sqrt(0.318*0.208),
    # V-V2 0.667*sqrt(0.318*0.690), Cl-V2 0.571*sqrt(0.208*0.690).
    etalvc + etalcl + etalvp ~ c(
      0.318,
      0.18337273, 0.208,
      0.31243806, 0.21631783, 0.690
    ) # Table 2: w2 V 0.318, w2 Cl 0.208, w2 V2 0.690; CorrV-Cl 0.713, CorrV-V2 0.667, CorrCl-V2 0.571
    etalq ~ 1.508 # Table 2: w2 Q = 1.508 (SE 0.46)

    # Residual variability -- Wang 2020 Results: 'a two-compartment model
    # with a proportional option was chosen'; Methods defines the proportional
    # model as Cobs = Cpred * (1 + eps). Table 2 row 'Residual variability
    # (s) stdev0' = 0.110, a standard deviation.
    propSd <- 0.110 ; label("Proportional residual error (fraction)") # Table 2: stdev0 = 0.110 (SE 0.01; bootstrap 95% CI 0.093-0.127)
  })

  model({
    # Individual parameters -- Wang 2020 Equ. 4-7.
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    cl <- exp(lcl + etalcl) * (CRCL / 105.9)^e_crcl_cl
    q <- exp(lq + etalq)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment disposition with first-order elimination from the
    # central compartment. Polymyxin B is given as an intravenous infusion
    # (0.5-2 h in the cohort; 50 mg/h in the Table 3 simulations), so drug
    # enters the central compartment directly.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
