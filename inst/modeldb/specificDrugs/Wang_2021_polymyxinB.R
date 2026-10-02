Wang_2021_polymyxinB <- function() {
  description <- "Two-compartment intravenous population PK model for polymyxin B (sum of polymyxin B1 and B2) in obese Chinese adults (BMI >= 30) with resistant Gram-negative infections, sampled at steady state (Wang 2021). No covariate was retained: age, total, ideal and adjusted body weight, BMI, sex, SOFA score, three creatinine clearance estimates, serum creatinine and GFR were screened and rejected. Log-normal inter-individual variability on CL, V2 and Q (none on V); proportional residual error."
  reference <- paste(
    "Wang P, Zhang Q, Feng M, Sun T, Yang J, Zhang X.",
    "Population Pharmacokinetics of Polymyxin B in Obese Patients for",
    "Resistant Gram-Negative Infections.",
    "Front Pharmacol. 2021;12:754844.",
    "doi:10.3389/fphar.2021.754844. PMCID PMC8645997.",
    sep = " "
  )
  vignette <- "Wang_2021_polymyxinB"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # What each ODE state holds, in what amount units, in what biological
  # matrix. Verified against Wang 2021: plasma polymyxin B1 and B2 were
  # quantified by HPLC-MS/MS and summed to a total polymyxin B concentration
  # (Methods, 'Polymyxin B Administration and Assay'), and the final model
  # is two-compartment with central volume V and peripheral volume V2
  # (Results, 'Population PK Model'; Table 2).
  compartmentData <- list(
    central = list(analyte = "polymyxinB", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "polymyxinB", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    SCREENED_NOT_RETAINED = list(
      description = "Candidate covariates screened by stepwise forward addition / backward elimination and not retained",
      units = "(various)",
      type = "continuous",
      notes = "Wang 2021 Methods 'Population PK Modeling': age, TBW, BMI, IBW, ABW, sex, SOFA score, CrCL, adjusted CrCL, ideal CrCL, serum creatinine and GFR were tested (forward dOFV > 6.63, backward dOFV > 10.83). Results: 'age, TBW, BMI, IBW, ABW, sex, SOFA score, CrCLs, serum creatinine, and GFR had no systematic relationship with PK parameters.' IBW = [height (cm)/2.54 - 60] x 2.3 kg + 50 kg (male) or 45.5 kg (female); ABW = IBW + 0.4 x (TBW - IBW); the three CrCL variants are Cockcroft-Gault with TBW, ABW or IBW. The Discussion attributes the absence of a body-weight effect to the narrow weight range (75-125 kg). The body-weight measures enter only the mg/kg dosing of the paper's Monte Carlo regimens, not the model."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 26L,
    n_studies = 1L,
    age_range = "18-83 years",
    age_median = "52 years",
    weight_range = "75-125 kg (total body weight)",
    weight_median = "90 kg",
    sex_female_pct = 34.62,
    race_ethnicity = "Chinese (single-centre cohort in Zhengzhou, China)",
    disease_state = "Obese adults (BMI >= 30; range 30.04-40.35, median 32.65) receiving intravenous polymyxin B sulfate for multidrug-resistant Gram-negative infections (lung 24, bloodstream 9, abdomen 4, intracranial 1; Klebsiella pneumoniae 14, Acinetobacter baumannii 13, Pseudomonas aeruginosa 3, Escherichia coli 2). Patients on CRRT or ECMO were excluded. SOFA score median 10 (5-17).",
    dose_range = "Intravenous polymyxin B sulfate, 100-200 mg loading dose then 50-100 mg twice daily, infused over at least 1 h; daily dose per body weight median 1.65 (range 0.92-2.70) mg/kg.",
    regions = "China (First Affiliated Hospital of Zhengzhou University), April 2018 to March 2021",
    renal_function = "Cockcroft-Gault CrCL (TBW) median 84.04 (range 21.35-239.99) mL/min; adjusted CrCL 71.83 (19.57-201.67); ideal CrCL 65.94 (13.50-125.71); serum creatinine 76 (34-368) umol/L; GFR 80.90 (11.26-149.57) mL/min/1.73 m^2 (Table 1).",
    ibw_abw = "IBW median 65.94 (48.76-74.99) kg; ABW median 75.15 (59.26-92.82) kg (Table 1).",
    n_observations = "142 plasma polymyxin B concentrations (sum of B1 and B2). 10 patients (67 samples) had 4-7 samples within one steady-state dosing interval on day 4 (mainly 0, 1, 1.5, 2, 4, 6 and 8 h); 16 patients (75 samples) had routine two-point TDM samples (pre-dose and 2 h after the start of a 1-h infusion).",
    notes = "Retrospective single-centre study. Phoenix NLME 7.0, FOCE-ELS. Evaluation by goodness-of-fit plots, a prediction-corrected VPC (1,000 replicates) and a 1,000-sample bootstrap (Table 2). Assay: HPLC-MS/MS, polymyxin B1 0.2-10.0 ug/mL and B2 0.05-2.5 ug/mL. Demographics from Table 1."
  )

  ini({
    # Structural parameters -- Wang 2021 Table 2 (final model). No covariate
    # was retained, so these are the population values for every patient.
    lvc <- log(11.24) ; label("Central volume of distribution V (L)")          # Table 2: tvV = 11.24 L (SE 1.56, RSE 13.87%; bootstrap median 11.46)
    lvp <- log(39.70) ; label("Peripheral volume of distribution V2 (L)")      # Table 2: tvV2 = 39.70 L (SE 12.09, RSE 30.46%; bootstrap median 40.46)
    lcl <- log(2.86)  ; label("Clearance CL (L/h)")                            # Table 2: tvCL = 2.86 L/h (SE 0.25, RSE 8.62%; bootstrap median 2.73)
    lq  <- log(7.36)  ; label("Intercompartmental clearance Q (L/h)")          # Table 2: tvQ = 7.36 L/h (SE 1.57, RSE 21.38%; bootstrap median 7.71)

    # Inter-individual variability -- Table 2 'w2' rows, footnoted as the
    # 'variance of inter-individual variability'. Encoded as log-normal
    # (P_i = theta * exp(eta_i)), the Phoenix NLME default, which reproduces
    # the paper's Monte Carlo Table 3 (see the vignette). No IIV on V
    # (Results: 'Because of shrinkage factor > 0.5, the random effect of V was
    # not taken into the model'); no correlations ('No correlation between
    # random effects was found during modeling').
    etalcl ~ 0.17 # Table 2: w2 CL = 0.17 (SE 0.05; bootstrap median 0.17)
    etalvp ~ 1.00 # Table 2: w2 V2 = 1.00 (SE 0.47; bootstrap median 0.87)
    etalq ~ 0.43  # Table 2: w2 Q = 0.43 (SE 0.14; bootstrap median 0.59)

    # Residual variability -- Table 2 row 'Residual variability (s) stdev0'
    # = 0.24, a standard deviation. The Results call the error model
    # 'proportional'; the Figure 2C observed-vs-IPRED scatter widens with
    # concentration, consistent with a proportional form.
    propSd <- 0.24 ; label("Proportional residual error (fraction)") # Table 2: stdev0 = 0.24 (SE 0.02; bootstrap 95% CI 0.18-0.30)
  })

  model({
    # Individual parameters -- Wang 2021 Table 2
    vc <- exp(lvc)
    vp <- exp(lvp + etalvp)
    cl <- exp(lcl + etalcl)
    q <- exp(lq + etalq)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment disposition with first-order elimination from the
    # central compartment. Polymyxin B is given as an intravenous infusion
    # (at least 1 h in the cohort; 50 mg/h in the paper's simulations), so
    # drug enters the central compartment directly.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
