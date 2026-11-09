Tang_2023_polymyxinB <- function() {
  description <- "Two-compartment intravenous population PK model for polymyxin B in critically ill adults with carbapenem-resistant organism nosocomial pneumonia, developed from a prospective Chinese ICU cohort (Tang 2023). No covariate was retained. Inter-individual variability on CL and Vc only, with a 100% CL-Vc correlation encoded as a single shared random effect scaled onto Vc; variability on CLd and Vp was held at zero. Proportional residual error."
  reference <- paste(
    "Tang T, Li Y, Xu P, Zhong Y, Yang M, Ma W, Xiang D, Zhang B, Zhou Y.",
    "Optimization of polymyxin B regimens for the treatment of",
    "carbapenem-resistant organism nosocomial pneumonia: a real-world",
    "prospective study. Crit Care. 2023;27:164.",
    "doi:10.1186/s13054-023-04448-z. PMCID PMC10142183.",
    sep = " "
  )
  vignette <- "Tang_2023_polymyxinB"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # What each ODE state holds. Plasma total polymyxin B (B1 + B2, combined on a
  # molar basis) was quantified by HPLC-MS/MS (Methods, 'Drug administration and
  # concentration determination'); the disposition model is the two-compartment
  # structure of the Results 'Population PK Model' final-model equations.
  compartmentData <- list(
    central = list(analyte = "polymyxinB", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "polymyxinB", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      notes = "Tang 2023 Methods: screened by stepwise covariate modeling (forward p < 0.05, backward p < 0.001). Results: 'no covariate was statistically significant to PK parameters'. Cohort median 65 years (IQR 55-76), Table 1."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "categorical",
      notes = "Tang 2023 Methods: sex screened by stepwise covariate modeling and not retained. 29 of 105 patients (27.6%) were female (Table 1)."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Tang 2023 Methods: screened by stepwise covariate modeling and not retained. Cohort median 55.0 kg (IQR 50.0-65.0), Table 1. The Discussion attributes the null weight effect to the small sample and the narrow weight distribution."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Tang 2023 Methods: screened by stepwise covariate modeling and not retained. Cohort median 30.7 g/L (IQR 27.9-34.2), Table 1."
    ),
    CRCL = list(
      description = "Cockcroft-Gault creatinine clearance (not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      notes = "Tang 2023 Methods: 'Creatinine clearance (CrCL) was calculated using the Cockcroft-Gault equation', and CrCL was screened by stepwise covariate modeling and not retained. Baseline median 72.2 mL/min (IQR 50.5-111.5), Table 1. Patients on renal replacement therapy were excluded."
    ),
    OTHER_LIVER_LABS = list(
      description = "Liver-function laboratory values screened as candidate covariates",
      units = "(various)",
      type = "continuous",
      notes = "Tang 2023 Methods lists alanine aminotransferase (ALT), aspartate aminotransferase (AST), total bilirubin (TBIL) and direct bilirubin (DBIL) among the screened covariates; none was retained. Grouped into one entry because the paper reports no per-analyte result and no cohort values for these laboratories."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 105L,
    n_studies = 1L,
    age_median = "65 years (IQR 55-76)",
    weight_median = "55.0 kg (IQR 50.0-65.0)",
    sex_female_pct = 27.6,
    race_ethnicity = "Not reported (single-center Chinese ICU cohort)",
    disease_state = "Critically ill adults (>= 18 years) with hospital-acquired (nosocomial) pneumonia caused by carbapenem-resistant organisms: carbapenem-resistant Acinetobacter baumannii 81.9%, Klebsiella pneumoniae 40.0% and Pseudomonas aeruginosa 21.9% (patients could carry more than one pathogen). Median APACHE II 18 (IQR 12-25); sepsis 39.0%; mechanical ventilation 77.1%. All patients received combination therapy with another antimicrobial. Exclusions: lung cancer with obstructive pneumonitis or cystic fibrosis, solid-organ transplant, hematologic malignancy or hematopoietic cell transplant, and renal replacement therapy.",
    dose_range = "Intravenous polymyxin B: loading dose 100-200 mg (median 100 mg, IQR 100-150; given in 87 of 105 patients) followed by 40-100 mg every 12 h for at least 3 days, infused over at least 1 h; median daily dose 150 mg (IQR 100-150), 2.3 mg/kg/day; median treatment duration 12 days. 34 patients also received nebulized polymyxin B (25 or 50 mg twice daily), which the PK model does not describe. Monte Carlo simulations used 100-200 mg loading doses followed by 50-150 mg q12h as 2-hour infusions.",
    regions = "China (two intensive care units of the Second Xiangya Hospital of Central South University, Changsha), January 2020 to December 2021",
    renal_function = "Baseline Cockcroft-Gault creatinine clearance median 72.2 mL/min (IQR 50.5-111.5); chronic renal dysfunction 23.8%. Renal replacement therapy was an exclusion criterion.",
    n_observations = "295 plasma concentrations from 105 patients (about three per patient), sampled randomly (two to six samples) immediately before the seventh dose and at 0, 1, 2, 4, 6, 8 and 10 h after the end of infusion. All concentrations were above the 0.03 mg/L quantification limit.",
    notes = "Prospective study (ChiCTR1900022231). Phoenix NLME 8.1, FOCE-ELS. Model evaluation by goodness-of-fit plots, prediction-corrected VPC (Figure 2) and bootstrap (Table 2). Total polymyxin B computed from measured polymyxin B1 and B2 on a molar basis. The paper's PK/PD analysis (multivariable logistic regression of clinical efficacy on AUCss,24h/MIC, daily dose and inhaled polymyxin B) reports adjusted odds ratios without an intercept and is not encoded here."
  )

  ini({
    # Structural parameters -- Tang 2023 Table 2 (final model) and the
    # final-model equations in Results 'Population PK Model and polymyxin B
    # exposure': CL = 1.56 * exp(eta_CL), Vc = 12.5 * exp(eta_Vc), CLd = 2.41,
    # Vp = 29.9. No covariates were retained.
    lcl <- log(1.56); label("Clearance CL (L/h)")                             # Table 2: CL = 1.56 L/h (SE 4.50%; bootstrap 1.56)
    lvc <- log(12.5); label("Central volume of distribution Vc (L)")          # Table 2: Vc = 12.5 L (SE 5.67%; bootstrap 12.5)
    lq  <- log(2.41); label("Intercompartmental clearance CLd (L/h)")         # Table 2: CLd = 2.41 L/h (SE 11.8%; bootstrap 2.41)
    lvp <- log(29.9); label("Peripheral volume of distribution Vp (L)")       # Table 2: Vp = 29.9 L (SE 17.3%; bootstrap 29.9)

    # Inter-individual variability. Table 2 reports 'omega CL' = 39.4% and
    # 'omega Vc' = 28.6% under the heading 'Random effect (%)'. The percentages
    # are read on the omega (SD) scale, omega = pct / 100, the reading the
    # sibling polymyxin B models Wang_2026_polymyxinB.R and
    # Yang_2025_polymyxinB.R use for an 'omega (%)' row; the alternative
    # sqrt(log(1 + CV^2)) reading gives 0.381 and 0.281 (3-2% lower).
    etalcl ~ 0.155236 # Table 2: omega CL = 39.4% (SE 6.4%) -> omega = 0.394 -> variance 0.394^2

    # Table 2 'omega CL*Vc' = 100 (the CL-Vc correlation; Results: 'the
    # estimated covariance between CL and Vc is 100%'). A 2x2 omega block with
    # correlation 1 is singular and cannot be sampled, so the perfect
    # correlation is encoded as one shared random effect scaled onto Vc:
    # eta_Vc = vc_eta_scale * etalcl, with vc_eta_scale = omega_Vc / omega_CL =
    # 0.286 / 0.394 = 0.7259. This reproduces SD(eta_Vc) = 0.286 and
    # corr(eta_CL, eta_Vc) = 1 exactly. The ratio is derived, not estimated.
    vc_eta_scale <- fixed(0.7259); label("Scaling of the shared CL random effect onto Vc, reproducing their 100% correlation (unitless)") # derived: Table 2 omega Vc 28.6% / omega CL 39.4%; Table 2 'omega CL*Vc' = 100

    # Results: 'the shrinkage of ... CLd and Vp were more than 50%, so the
    # inter-individual variability (IIV) of CLd and Vp were fixed 0'. Kept as
    # explicit zero-variance random effects so that status stays visible.
    etalq  ~ fixed(0) # Results 'Population PK Model': IIV on CLd held at 0 (shrinkage > 50%)
    etalvp ~ fixed(0) # Results 'Population PK Model': IIV on Vp held at 0 (shrinkage > 50%)

    # Residual variability: proportional (Results: 'A proportional error model
    # was used to evaluate the residual variability'). Phoenix NLME reports the
    # proportional epsilon as its standard deviation.
    propSd <- 0.2786; label("Proportional residual error (fraction)") # Table 2: proportional residual = 27.86% (SE 9.11%)
  })

  model({
    # Individual parameters -- Tang 2023 Results final-model equations, with
    # the 100% CL-Vc correlation carried by the shared etalcl.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + vc_eta_scale * etalcl)
    q  <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment disposition with first-order elimination from the
    # central compartment; polymyxin B is infused intravenously into central
    # (Methods: infusion time at least 1 h; simulations used 2-h infusions).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
