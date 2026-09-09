Kumta_2025_piperacillin <- function() {
  description <- "Three-compartment population PK model for piperacillin in the plasma and cerebrospinal fluid of critically ill neurosurgical adults with an external ventricular drain (Kumta 2025): two-compartment plasma disposition with linear elimination from the central compartment, plus a small CSF compartment exchanging with the central compartment through a very low inter-compartmental clearance and with NO elimination of its own."
  reference <- "Kumta N, Heffernan AJ, Cotta MO, Liu X, Parker S, Wallis S, Livermore A, Starr T, Wong WT, Joynt GM, Lipman J, Roberts JA. Population pharmacokinetics of piperacillin-tazobactam in the plasma and cerebrospinal fluid of critically ill patients. Antimicrob Agents Chemother. 2025;69(2):e00601-24. doi:10.1128/aac.00601-24"
  vignette <- "Kumta_2025_piperacillin_tazobactam"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L"
  )

  # Issue #482: what molecule each compartment holds, in what units, in what
  # biological matrix. Verified against Kumta 2025 Figure 1 (schematic of the
  # final piperacillin PK model) and the Table 3 parameter definitions.
  compartmentData <- list(
    central     = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE),
    csf         = list(analyte = "piperacillin", units = "mg", specimen = "CSF",    verified = TRUE)
  )

  # Kumta 2025 Results, "Pharmacokinetic model": "No covariates improved model
  # diagnostics significantly and thus could be retained in the final model."
  # The final model is therefore the base model and carries NO covariateData.
  # The screened-but-not-retained set is documented below so the paper's
  # covariate search is preserved without triggering an unused-covariate
  # warning.
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
      notes       = "Screened and not retained; the final model carries no allometric term. Cohort median 70 kg, range 47-110 (Table 1); the Discussion quotes a mean of 71.5 kg."
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
      notes       = "Screened and not retained. Cohort median 25.5 g/L, range 22-32 (Table 1) -- uniformly hypoalbuminaemic, as is typical of critical illness."
    ),
    CRCL = list(
      description = "Creatinine clearance, Cockcroft-Gault on total body weight, BSA-normalized",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Screened and not retained -- the headline negative result of the paper. Methods: 'creatinine clearance (calculated using the Cockcroft-Gault equation using total body weight and expressed in mL/min/1.73 m^2)'. Cohort median 84, range 52-163 (Table 1 / Table 2); one patient (12.5%) met the augmented-renal-clearance threshold of 130. Discussion: 'clearance of piperacillin or tazobactam was not influenced by creatinine clearance or any other covariates'. The cohort excluded renal replacement therapy and plasma creatinine above 200 umol/L, so the model carries no information about renal impairment."
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
    dose_range     = "Piperacillin-tazobactam 4.5 g every 6 h by intermittent intravenous infusion (7 of 8 patients; 4 g piperacillin + 0.5 g tazobactam per dose); one patient received a continuous infusion totalling 13.1 g over the sampling period",
    regions        = "Two university-associated tertiary ICUs: Royal Brisbane and Women's Hospital (Australia) and the Chinese University of Hong Kong (Hong Kong, China)",
    notes          = "45 plasma and 30 CSF samples from 8 patients (Results, 'Study population'). Plasma sampled 0.5, 1, 1.5, 2, 4 and 6 h and CSF 0.5, 2, 4 and 6 h after the start of infusion for intermittent-bolus patients; the single continuous-infusion patient was sampled 17-26 h after commencement. Total (not unbound) piperacillin was assayed by UHPLC-MS/MS with a lower limit of quantification of 0.5 mg/L. Estimation used SAEM in Monolix 2023R1; accuracy was assessed by a 1,000-run bootstrap (Rsmlx 2023.1.1). Other baseline characteristics (Table 1): albumin median 25.5 g/L (22-32), CSF protein median 0.73 g/L (0.19-1.8), CSF protein:serum albumin median 0.016 (0.006-0.056), CSF volume drained from the EVD median 37.5 mL (0-67), SOFA median 7.5 (2-12), APACHE median 18 (12-27)."
  )

  ini({
    # ---- Plasma disposition (Kumta 2025 Table 3, piperacillin "Estimate (%RSE)" column) ----
    lcl   <- log(12.7);    label("Clearance from the central compartment (Cl, L/h)")                                  # Table 3: Cl 12.7 L/h (%RSE 11.1); bootstrap median 12.6 (95% CI 10.3-15.8)
    lvc   <- log(13.4);    label("Central compartment volume (V1, L)")                                                # Table 3: V1 13.4 L (%RSE 15.0); bootstrap median 13.5 (95% CI 9.35-17.6)
    lq    <- log(7.25);    label("Inter-compartmental clearance central <-> peripheral (Q, L/h)")                     # Table 3: Q 7.25 L/h (%RSE 3.41); bootstrap median 6.85 (95% CI 1.85-70.7)
    lvp   <- log(4.99);    label("Peripheral compartment volume (V2, L)")                                             # Table 3: V2 4.99 L (%RSE 28.5); bootstrap median 4.98 (95% CI 2.94-12.6)

    # ---- CSF compartment (Kumta 2025 Fig. 1 and Table 3) ----
    # Fig. 1 draws the CSF compartment exchanging with the CENTRAL compartment
    # only, with rate constants Q3/V1 (out of central) and Q3/V3 (out of CSF),
    # and draws NO elimination arrow from the CSF box. Results: "a three-
    # compartment model without clearance from the CSF compartment".
    lqcsf <- log(0.00024); label("Inter-compartmental clearance central <-> CSF (Q3, L/h)")                           # Table 3: Q3 0.00024 L/h (%RSE 85.3); bootstrap median 0.00025 (95% CI 0.00013-0.00047)
    lvcsf <- log(0.16);    label("CSF compartment volume (V3, L)")                                                    # Table 3: V3 0.16 L (%RSE 80.2); bootstrap median 0.16 (95% CI 0.12-0.25)

    # ---- Between-subject variability ----
    # Methods: BSV used the exponential model theta_j = theta_p * exp(eta_j)
    # with eta_j ~ N(0, omega^2). Table 3 reports the "Random effect" rows as
    # percentages ("BSV_Cl (%)" etc.), which is the Monolix coefficient-of-
    # variation convention CV% = sqrt(exp(omega^2) - 1) * 100. They are
    # converted here with omega^2 = log(CV^2 + 1). See the vignette
    # "Assumptions and deviations" for the alternative reading (that the
    # printed percentages are omega itself x 100) and why it was not adopted.
    # BSV was estimated on Cl, V1 and Q3 only -- V2, Q and V3 carry none.
    etalcl   ~ 0.083445  # log(0.295^2 + 1); Table 3: BSV_Cl 29.5% (%RSE 27.8), bootstrap median 27.1% (95% CI 10.6-36.1)
    etalvc   ~ 0.070365  # log(0.270^2 + 1); Table 3: BSV_V1 27.0% (%RSE 36.0), bootstrap median 24.8% (95% CI 4.25-35.8)
    etalqcsf ~ 0.556622  # log(0.863^2 + 1); Table 3: BSV_Q3 86.3% (%RSE 28.4), bootstrap median 75.6% (95% CI 22.7-124)

    # ---- Residual error ----
    # Results: "Residual variability was best described by a combined (additive
    # plus proportional) error model for plasma concentrations, while a
    # proportional error model was selected for CSF concentrations."
    addSd       <- 2.10; label("Additive residual SD on plasma piperacillin Cc (mg/L)")                               # Table 3: Additive residual_plasma 2.10 mg/L (%RSE 31.7); bootstrap median 2.26 (95% CI 0.45-3.75)
    propSd      <- 0.08; label("Proportional residual SD on plasma piperacillin Cc (fraction)")                       # Table 3: Proportional_plasma 0.08 (%RSE 31.1); bootstrap median 0.07 (95% CI 0.01-0.12)
    propSd_Ccsf <- 0.30; label("Proportional residual SD on CSF piperacillin Ccsf (fraction)")                        # Table 3: Proportional_CSF 0.30 (%RSE 16.3); bootstrap median 0.30 (95% CI 0.18-0.42)
  })

  model({
    # Individual parameters. No covariate enters any of them: the final model
    # is the base model (Results, "Pharmacokinetic model"). BSV was estimated
    # on Cl, V1 and Q3 only, so Q, V2 and V3 have no eta.
    cl   <- exp(lcl   + etalcl)
    vc   <- exp(lvc   + etalvc)
    q    <- exp(lq)
    vp   <- exp(lvp)
    qcsf <- exp(lqcsf + etalqcsf)
    vcsf <- exp(lvcsf)

    # Concentrations (mg/L). States are amounts in mg.
    Cc   <- central     / vc    # Fig. 1: Cc,   plasma concentration in the central compartment
    Cp   <- peripheral1 / vp    # Fig. 1: Cp,   concentration in the peripheral compartment
    Ccsf <- csf         / vcsf  # Fig. 1: CCSF, concentration in CSF

    # Kumta 2025 Fig. 1. Elimination (Cl/V1) leaves the central compartment
    # only; the central <-> peripheral leg is the symmetric pair Q/V1 and
    # Q/V2; the central <-> CSF leg is the symmetric pair Q3/V1 and Q3/V3.
    # There is deliberately no elimination term on the CSF compartment.
    d/dt(central)     <- -cl * Cc - q * Cc + q * Cp - qcsf * Cc + qcsf * Ccsf
    d/dt(peripheral1) <-  q * Cc - q * Cp
    d/dt(csf)         <-  qcsf * Cc - qcsf * Ccsf

    Cc   ~ add(addSd) + prop(propSd)
    Ccsf ~ prop(propSd_Ccsf)
  })
}
