Shiau_2026_vancomycin <- function() {
  description <- "Two-compartment IV population PK model for vancomycin in critically ill children with multiple organ dysfunction syndrome (MODS), 1 month to 17 years (Shiau 2026). Clearance and intercompartmental clearance scale allometrically with body weight (exponent 0.75, reference 28.4 kg), and clearance additionally scales as a power function of CKiD Under-25 (U25) estimated GFR (exponent 0.85, reference 96 mL/min/1.73 m^2); central and peripheral volumes scale linearly with body weight (exponent 1, reference 28.4 kg). Between-subject variability is exponential on all four structural parameters and residual variability is proportional. The source is a conference poster abstract, but its Table 1 reports the complete Monolix parameter set together with the individual-parameter equations, so no value in this file is inferred."
  reference <- "Shiau J, Amajor V, Marianski S, Rhodes NJ, Bwint A, Sharova A, Hall M, Pai MP, Wen B, Downes KJ, Scheetz MH. P-1237. Vancomycin Population Pharmacokinetics and Toxicity-Exposure Relationships in Children with Multiple Organ Dysfunction Syndrome. Open Forum Infect Dis. 2026;13(Suppl 1):S810. doi:10.1093/ofid/ofaf695.1429. IDWeek 2025 poster abstract (Session 148, PK/PD Studies); PMCID PMC12791793."
  vignette <- "Shiau_2026_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Vancomycin is given intravenously, so the dose enters
  # `central` directly and there is no depot state.
  #
  # Both entries are left unverified on `specimen`. The abstract's Methods
  # state only that "up to 15 PK samples were collected via volumetric
  # absorptive microsampling over 3 days" -- volumetric absorptive
  # microsampling (VAMS) collects capillary WHOLE BLOOD, but the abstract
  # never states whether the assayed and modelled concentrations are
  # whole-blood values or converted plasma-equivalent values, and Figure 1
  # labels its axis only "Observations (mg/L)". "plasma" therefore follows the
  # repository default for a PK central compartment and is NOT a
  # paper-sourced claim; a reader matching this model to whole-blood VAMS
  # data should confirm the matrix against the full publication when it
  # appears. `peripheral1` is a mathematical distribution compartment whose
  # matrix the abstract never discusses.
  compartmentData <- list(
    central     = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference weight 28.4 kg is the normalizing constant printed in all four individual-parameter equations of Shiau 2026 Table 1, and is named in the Table 1 caption as 'weight adjusted (individual weight/28.4kg)'. The abstract's Results report a cohort median weight of 30 kg (range 3-214 kg), so 28.4 kg is close to but not identical with the reported median; the abstract does not say how it was chosen. The abstract does not state whether weight was treated as time-varying over the 3-day sampling window.",
      source_name        = "Wt"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate calculated with the CKiD Under-25 (U25) equation, BSA-normalized",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Shiau 2026 Results name the estimating equation ('glomerular function (calculated using CKiD Under 25 [U25] equation)') and Table 1's clearance equation prints the reference value 96 as (U25/96)^Theta_GFR. The abstract does NOT reproduce the U25 equation itself, so this file does not restate it; the U25 equation is a paediatric creatinine- and/or cystatin-C-based eGFR estimator published by Pierce et al. and is applied when preparing the covariate column, not inside the model. The abstract reports no cohort eGFR distribution. Stored under canonical CRCL, which covers BSA-normalized creatinine-based GFR estimates; document the estimator per model because size normalizations are not interchangeable.",
      source_name        = "U25"
    )
  )

  population <- list(
    species          = "human",
    n_subjects        = 66L,
    n_studies         = 1L,
    age_range         = "1 month to 17 years",
    age_median        = "10 years (range 1 month - 17 years)",
    weight_range      = "3 - 214 kg",
    weight_median     = "30 kg (range 3 - 214)",
    disease_state     = "Critically ill children with multiple organ dysfunction syndrome (MODS) receiving intravenous vancomycin for suspected or documented severe gram-positive infection. MODS severity was classified with the Proulx score; the abstract reports no baseline renal-function strata.",
    dose_range        = "Not reported. The abstract states only that 'most initial VAN dosing is guided via weight-based approaches' in this population; no dose amounts, infusion durations or dosing intervals are given.",
    regions           = "United States (multicenter; PALISI network sites, including Children's Hospital of Philadelphia and Nationwide Children's Hospital)",
    renal_function    = "Reference CKiD U25 eGFR 96 mL/min/1.73 m^2 is the normalizing constant in the Shiau 2026 Table 1 clearance equation. The abstract reports no cohort eGFR median or range.",
    exposure_observed = "Empirical-Bayes-estimate AUCs computed in Simulx 2024R1: median AUC0-24 454 mg*h/L (range 194-1569) and median AUC24-48 505 mg*h/L (range 8-1994).",
    toxicity          = "7 of 66 subjects met criteria for ICU-emergent acute kidney injury (a 0.3 mg/dL or 50% rise in serum creatinine from ICU baseline). An AUC24-48 of 465.8 mg*h/L retained 100% sensitivity for AKI; sensitivity fell to 20% or below at AUC24-48 of 545.5 mg*h/L or above. In a stepwise multivariable logistic regression only the Proulx MODS score was significant (odds ratio 6.675, 95% CI 1.09-40.81, p = 0.04); AUC0-24 was not (odds ratio 0.998, 95% CI 0.99-1.00, p = 0.343).",
    notes             = "Multicenter prospective observational PK study (AMPLE) embedded in a larger study of critically ill children with MODS (PARADIGM), conducted within the Pediatric Acute Lung Injury and Sepsis Investigators (PALISI) network. Up to 15 PK samples per child were collected by volumetric absorptive microsampling over 3 days. Parametric population PK modeling was performed in Monolix 2024R1 with covariate inclusion based on objective function and physiologic relevance; individual AUCs were computed from empirical Bayes estimates in Simulx 2024R1. Figure 1 shows censored (below-limit-of-quantification) observations, so the fit used Monolix's censored-data likelihood; the censoring limit is not reported. Demographics are from the abstract Results paragraph -- the abstract has no baseline-demographics table."
  )

  ini({
    # Structural parameters (Shiau 2026 Table 1, 'Fixed Effects' block). The
    # reference subject weighs 28.4 kg and has a CKiD U25 eGFR of
    # 96 mL/min/1.73 m^2, per the individual-parameter equations printed at the
    # foot of the same table.
    lcl <- log(2.03); label("Clearance at WT=28.4 kg and CRCL=96 mL/min/1.73 m^2 (L/h)")  # Shiau 2026 Table 1 'Cl (L/hr)': 2.03 (S.E. 0.15, R.S.E. 7.18%, P2.5-P97.5 1.76-2.33)
    lvc <- log(7.97); label("Central volume at WT=28.4 kg (L)")                           # Shiau 2026 Table 1 'V1 (L)': 7.97 (S.E. 1.25, R.S.E. 15.71%, P2.5-P97.5 5.89-10.79)
    lq  <- log(2.54); label("Intercompartmental clearance at WT=28.4 kg (L/h)")            # Shiau 2026 Table 1 'Q (L/hr)': 2.54 (S.E. 0.46, R.S.E. 18.08%, P2.5-P97.5 1.8-3.59)
    lvp <- log(9.59); label("Peripheral volume at WT=28.4 kg (L)")                         # Shiau 2026 Table 1 'V2 (L)': 9.59 (S.E. 1.45, R.S.E. 15.09%, P2.5-P97.5 7.17-12.83)

    # Allometric exponents. Shiau 2026 Table 1 prints them inside the
    # individual-parameter equations (0.75 on Cl and Q; an explicit exponent of
    # 1 on V1 and V2) and gives them no S.E. / R.S.E. row, and the Results text
    # names the structure "clearance adjusted for allometric scaling
    # (weight^0.75)". Both are therefore canonical fixed allometric exponents,
    # not estimated parameters.
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent on (WT/28.4) for CL and Q (unitless)")   # Shiau 2026 Table 1 equations: (Wt/28.4)^0.75 on Cl and on Q
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent on (WT/28.4) for Vc and Vp (unitless)")  # Shiau 2026 Table 1 equations: (Wt/28.4)^1 on V1 and on V2

    # Renal-function effect on clearance. Estimated (it carries its own S.E. /
    # R.S.E. row), so it is not wrapped in fixed().
    e_crcl_cl <- 0.85; label("Power exponent on (CRCL/96) for CL (unitless)")  # Shiau 2026 Table 1 'Theta_U25 GFR': 0.85 (S.E. 0.064, R.S.E. 7.6%, P2.5-P97.5 0.72-0.97)

    # Between-subject variability. Shiau 2026 Table 1 heads this block
    # 'Standard Deviation of the Random Effects' and prints omega directly as a
    # log-scale SD alongside a 'C.V. (%)' column; nlmixr2 wants the VARIANCE,
    # so each entry below is omega^2. The two columns identify the convention
    # unambiguously as lognormal CV = sqrt(exp(omega^2) - 1): back-solving each
    # printed CV% recovers the printed omega to within its 2-decimal rounding
    # (57.37% -> 0.5334, 75.05% -> 0.6684, 76.89% -> 0.6815, 84.70% -> 0.7354),
    # which rules out reading the omega column as a variance or the CV column
    # as omega x 100. Table 1 reports no correlations between random effects,
    # so the block is diagonal.
    etalcl ~ 0.2809  # 0.53^2; Shiau 2026 Table 1 'omega_Cl (L/hr)' = 0.53 (C.V. 57.37%, R.S.E. 9.79%)
    etalvc ~ 0.4489  # 0.67^2; Shiau 2026 Table 1 'omega_V1 (L)'    = 0.67 (C.V. 75.05%, R.S.E. 17.38%)
    etalq  ~ 0.4624  # 0.68^2; Shiau 2026 Table 1 'omega_Q (L/hr)'  = 0.68 (C.V. 76.89%, R.S.E. 20.08%)
    etalvp ~ 0.5476  # 0.74^2; Shiau 2026 Table 1 'omega_V2 (L)'    = 0.74 (C.V. 84.70%, R.S.E. 16.81%)

    # Residual variability. Shiau 2026 Table 1's 'Error Model Parameters'
    # block reports a single parameter, b = 0.25. In Monolix's error-model
    # vocabulary b alone (with no a) is the PROPORTIONAL model
    # y = f + b * f * e, which is exactly nlmixr2's prop(). A constant or
    # combined model would have reported an 'a' as well, and none is present.
    propSd <- 0.25; label("Proportional residual error (fraction)")  # Shiau 2026 Table 1 'b': 0.25 (S.E. 0.0092, R.S.E. 3.66%, P2.5-P97.5 0.23-0.27)
  })
  model({
    # Individual parameters, transcribed from the equations printed at the foot
    # of Shiau 2026 Table 1:
    #   Cl_i = Cl * (Wt/28.4)^0.75 * (U25/96)^Theta_GFR * exp(eta_CL)
    #   V1_i = V1 * (Wt/28.4)^1                         * exp(eta_V1)
    #   Q_i  = Q  * (Wt/28.4)^0.75                      * exp(eta_Q)
    #   V2_i = V2 * (Wt/28.4)^1                         * exp(eta_V2)
    cl <- exp(lcl + etalcl) * (WT / 28.4)^e_wt_cl_q * (CRCL / 96)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 28.4)^e_wt_vc_vp
    q  <- exp(lq  + etalq)  * (WT / 28.4)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 28.4)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L, so central/vc is mg/L == ug/mL, matching the
    # "Observations (mg/L)" axis of Shiau 2026 Figure 1.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
