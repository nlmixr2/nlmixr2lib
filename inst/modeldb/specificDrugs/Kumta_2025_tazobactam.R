Kumta_2025_tazobactam <- function() {
  description <- "Two-compartment population PK model for tazobactam in the plasma of critically ill neurosurgical adults with an external ventricular drain (Kumta 2025): linear elimination from the central compartment, no covariates retained. CSF tazobactam was below the limit of quantification in 47% of samples and showed no change across the dosing interval, so the paper deliberately did NOT model a CSF compartment for tazobactam."
  reference <- "Kumta N, Heffernan AJ, Cotta MO, Liu X, Parker S, Wallis S, Livermore A, Starr T, Wong WT, Joynt GM, Lipman J, Roberts JA. Population pharmacokinetics of piperacillin-tazobactam in the plasma and cerebrospinal fluid of critically ill patients. Antimicrob Agents Chemother. 2025;69(2):e00601-24. doi:10.1128/aac.00601-24"
  vignette <- "Kumta_2025_piperacillin_tazobactam"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L"
  )

  # Issue #482: what molecule each compartment holds, in what units, in what
  # biological matrix. Verified against Kumta 2025 Results ("A two-compartment
  # model with first order elimination best described the PK of tazobactam in
  # plasma") and the Table 3 parameter definitions.
  compartmentData <- list(
    central     = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  # Kumta 2025 Results, "Pharmacokinetic model": "Likewise, no covariate effect
  # was identified." The final tazobactam model is the base model and carries
  # no covariateData. The screened-but-not-retained set is documented here so
  # the paper's covariate search is preserved without triggering an
  # unused-covariate warning. The screen was the same one applied to
  # piperacillin -- see modellib("Kumta_2025_piperacillin").
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened by stepwise forward inclusion / backward elimination (Methods, 'Population pharmacokinetic model development') and not retained. Cohort median 59 years, range 42-75 (Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened and not retained. 5 of 8 patients (62.5%) were female (Table 1)."
    ),
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened and not retained; the final model carries no allometric term. Cohort median 70 kg, range 47-110 (Table 1). The Discussion nevertheless attributes the low tazobactam V1 of 7.6 L (against 19 L in healthy volunteers) to the cohort's lower mean weight of 71.5 kg, so body size is invoked narratively even though no weight term survived the covariate screen."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened and not retained. Recorded prospectively (Methods, 'Patient population') but no summary value is tabulated."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened and not retained. Recorded prospectively but no summary value is tabulated."
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation score at ICU admission",
      units       = "points",
      type        = "continuous",
      notes       = "Screened and not retained. Table 1 reports a median of 18 (range 12-27) under the label 'APACHE scores'; the paper does not state the APACHE version, and the score is mapped here to the APACHE_II canonical because that is the version in routine adult ICU use and the observed range is consistent with it. Because the covariate is documentation-only, the mapping carries no modelling consequence."
    ),
    SOFA = list(
      description = "Sequential Organ Failure Assessment score",
      units       = "points",
      type        = "continuous",
      notes       = "Screened and not retained. Cohort median 7.5, range 2-12 (Table 1). No SOFA canonical exists in inst/references/covariate-columns.md; because this entry is documentation-only and the column is never referenced in model(), no register entry was created."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened and not retained. Cohort median 25.5 g/L, range 22-32 (Table 1)."
    ),
    CRCL = list(
      description = "Creatinine clearance, Cockcroft-Gault on total body weight, BSA-normalized",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Screened and not retained. Methods: 'creatinine clearance (calculated using the Cockcroft-Gault equation using total body weight and expressed in mL/min/1.73 m^2)'. Cohort median 84, range 52-163 (Table 1 / Table 2). Discussion: 'clearance of piperacillin or tazobactam was not influenced by creatinine clearance or any other covariates'. The cohort excluded renal replacement therapy and plasma creatinine above 200 umol/L."
    ),
    ALT = list(
      description = "Serum alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened and not retained. No summary value is tabulated; the cohort excluded pre-existing hepatic dysfunction defined as gamma-glutamyl transferase above 200 IU/L."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 8L,
    n_studies      = 1L,
    age_range      = "42-75 years",
    age_median     = "59 years",
    weight_range   = "47-110 kg",
    weight_median  = "70 kg",
    sex_female_pct = 62.5,
    disease_state  = "Critically ill neurosurgical ICU adults with an external ventricular drain in situ and either a ventriculostomy-associated infection (n = 1, 12.5%) or an extracranial infection (pneumonia, n = 7, 87.5%)",
    renal_function = "Creatinine clearance median 84 mL/min/1.73 m^2 (range 52-163); 1 patient (12.5%) with augmented renal clearance (>= 130 mL/min/1.73 m^2). Renal replacement therapy and plasma creatinine > 200 umol/L were exclusion criteria.",
    dose_range     = "Piperacillin-tazobactam 4.5 g every 6 h by intermittent intravenous infusion (7 of 8 patients; 0.5 g tazobactam per dose); one patient received a continuous infusion totalling 13.1 g of the combination over the sampling period",
    regions        = "Two university-associated tertiary ICUs: Royal Brisbane and Women's Hospital (Australia) and the Chinese University of Hong Kong (Hong Kong, China)",
    notes          = "45 plasma and 30 CSF samples from 8 patients (Results, 'Study population'); the plasma samples supported this model. CSF tazobactam was below the 0.625 mg/L limit of quantification in 14 of 30 samples (47%) and showed no change across the dosing interval, so Methods state that 'these data were not used for population PK analysis' and Results that it was 'not possible to model CSF exposures'. Total (not unbound) tazobactam was assayed by UHPLC-MS/MS. Estimation used SAEM in Monolix 2023R1; accuracy was assessed by a 1,000-run bootstrap (Rsmlx 2023.1.1)."
  )

  ini({
    # ---- Plasma disposition (Kumta 2025 Table 3, tazobactam "Estimate (%RSE)" column) ----
    lcl <- log(11.7); label("Clearance from the central compartment (Cl, L/h)")                       # Table 3: Cl 11.7 L/h (%RSE 12.7); bootstrap median 11.4 (95% CI 9.00-15.5)
    lvc <- log(7.64); label("Central compartment volume (V1, L)")                                    # Table 3: V1 7.64 L (%RSE 43.2); bootstrap median 7.17 (95% CI 0.18-18.2)
    lq  <- log(46.5); label("Inter-compartmental clearance central <-> peripheral (Q, L/h)")         # Table 3: Q 46.5 L/h (%RSE 64.5); bootstrap median 48.8 (95% CI 6.18-216)
    lvp <- log(12.0); label("Peripheral compartment volume (V2, L)")                                 # Table 3: V2 12.0 L (%RSE 23.6); bootstrap median 15.7 (95% CI 5.08-20.4)

    # ---- Between-subject variability ----
    # Methods: BSV used the exponential model theta_j = theta_p * exp(eta_j)
    # with eta_j ~ N(0, omega^2). Table 3 reports the "Random effect" rows as
    # percentages ("BSV_Cl (%)" etc.), which is the Monolix coefficient-of-
    # variation convention CV% = sqrt(exp(omega^2) - 1) * 100. They are
    # converted here with omega^2 = log(CV^2 + 1). See the vignette
    # "Assumptions and deviations". BSV was estimated on Cl and V1 only --
    # Q and V2 carry none.
    etalcl ~ 0.106966  # log(0.336^2 + 1); Table 3: BSV_Cl 33.6% (%RSE 28.7), bootstrap median 33.4% (95% CI 9.85-40.4)
    etalvc ~ 0.079152  # log(0.287^2 + 1); Table 3: BSV_V1 28.7% (%RSE 67.7), bootstrap median 25.4% (95% CI 4.41-49.2)

    # ---- Residual error ----
    # Results: "Residual variability was best described by a combined additive
    # and proportional error model."
    addSd  <- 0.47; label("Additive residual SD on plasma tazobactam Cc (mg/L)")                     # Table 3: Additive residual_plasma 0.47 mg/L (%RSE 28.6); bootstrap median 0.51 (95% CI 0.22-0.88)
    propSd <- 0.15; label("Proportional residual SD on plasma tazobactam Cc (fraction)")             # Table 3: Proportional_plasma 0.15 (%RSE 20.0); bootstrap median 0.14 (95% CI 0.08-0.18)
  })

  model({
    # Individual parameters. No covariate enters any of them: the final model
    # is the base model (Results, "Pharmacokinetic model"). BSV was estimated
    # on Cl and V1 only, so Q and V2 have no eta.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # Concentrations (mg/L). States are amounts in mg.
    Cc <- central     / vc
    Cp <- peripheral1 / vp

    # Two-compartment disposition with first-order elimination from the
    # central compartment (Results, "Pharmacokinetic model").
    d/dt(central)     <- -cl * Cc - q * Cc + q * Cp
    d/dt(peripheral1) <-  q * Cc - q * Cp

    Cc ~ add(addSd) + prop(propSd)
  })
}
