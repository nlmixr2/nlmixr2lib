Petric_2023_vinpocetine <- function() {
  description <- "Two-compartment population PK model for apovincaminic acid (AVA), the active de-esterified metabolite of vinpocetine, in healthy adult male volunteers dosed with vinpocetine (Petric 2023). Only AVA (not the parent vinpocetine) is modelled; the reported CL/F, V1/F, Q/F, V2/F apparent parameters fold the fraction of the vinpocetine dose that appears in plasma as AVA and the vinpocetine oral bioavailability into the /F term. Absorption is described as a zero-order input of duration Tk0 into the central compartment preceded by an absorption lag Tlag. Formulation is the only significant covariate: the sustained-release beta-cyclodextrin complex Ultra Vinca is the reference; the Cavinton immediate-release tablet and the extemporaneous 10 mg / 5 mL oral solution enter as log-additive shifts on Tk0 and V1/F. Between-subject variability is placed on Tlag, Tk0, CL/F, V1/F, Q/F, V2/F with a correlation of 0.72 between CL/F and V1/F. Residual error is proportional."
  reference   <- "Petric Z, Paixao P, Filipe A, Guimaraes Morais J. Clinical Pharmacology of Vinpocetine: Properties Revisited and Introduction of a Population Pharmacokinetic Model for Its Metabolite, Apovincaminic Acid (AVA). Pharmaceutics. 2023;15(10):2502. doi:10.3390/pharmaceutics15102502."
  vignette    <- "Petric_2023_vinpocetine"
  units       <- list(time = "h", dosing = "mg", concentration = "ug/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "vinpocetine", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "vinpocetine", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    FORM_VINP_IR = list(
      description        = "Vinpocetine immediate-release tablet indicator (1 = Cavinton immediate-release 5 mg tablet, Organon / Gedeon Richter; 0 = sustained-release beta-cyclodextrin complex Ultra Vinca 10 mg tablet, Tecnimede -- the Petric 2023 reference formulation). Per-dose-occasion indicator (crossover design: each subject received all three formulations across occasions with a 7-day washout).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Ultra Vinca sustained-release beta-cyclodextrin complex; the typical-value Tk0 and V1/F reference in Petric 2023 Table 1).",
      notes              = "Paired with FORM_SOLUTION to encode the three-level Petric 2023 formulation stratification {Ultra Vinca SR beta-cyclodextrin (reference), Cavinton immediate-release, oral solution}. Both indicators = 0 selects the Ultra Vinca SR reference. Enters multiplicatively on the log-scale of Tk0 (exp(beta = -0.4)) and V1/F (exp(beta = -1.26)) per Petric 2023 Table 1.",
      source_name        = "Formulation#2 (categorical level of the three-level Formulation covariate)"
    ),
    FORM_SOLUTION = list(
      description        = "Oral-solution formulation indicator (1 = 10 mg / 5 mL extemporaneous oral solution prepared at the hospital pharmacy from the pure API supplied by Tecnimede; 0 = Ultra Vinca sustained-release beta-cyclodextrin 10 mg tablet, the Petric 2023 reference formulation). Per-dose-occasion indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Ultra Vinca sustained-release beta-cyclodextrin complex; paired with FORM_VINP_IR = 0 in the reference state).",
      notes              = "Paired with FORM_VINP_IR to encode the three-level Petric 2023 formulation stratification. Enters multiplicatively on the log-scale of Tk0 (exp(beta = -0.68)) and V1/F (exp(beta = -1.24)) per Petric 2023 Table 1. The Petric 2023 SR-tablet reference is a solid oral form (beta-cyclodextrin complex tablet), which matches the FORM_SOLUTION canonical requirement that the 0 level is a non-solution comparator; document the specific SR reference per this per-model note.",
      source_name        = "Formulation#3 (categorical level of the three-level Formulation covariate)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age at study entry (years).",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a potential covariate on Tlag, Tk0, CL/F, V1/F, Q/F, and V2/F but not retained in the final model (Petric 2023 Section 3.3: 'none of the continuous covariates showed a significant impact on parameter variability'). Distribution in the study cohort: min 20, median 23, max 35 years (Petric 2023 Supplementary Table S1)."
    ),
    HEIGHT = list(
      description = "Subject height (m).",
      units       = "m",
      type        = "continuous",
      notes       = "Screened as a potential covariate but not retained (Petric 2023 Section 3.3). Distribution: min 1.64, median 1.74, max 1.83 m (Supplementary Table S1)."
    ),
    BMI = list(
      description = "Body-mass index (kg / m^2). Petric 2023 screened weight as BMI rather than as raw kg.",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened as a potential covariate but not retained (Petric 2023 Section 3.3). Distribution: min 20.9, median 24.9, max 32.56 kg / m^2 (Supplementary Table S1)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 12L,
    n_studies      = 1L,
    age_range      = "20-35 years (median 23)",
    weight_range   = "66-93 kg (median 73.5)",
    sex_female_pct = 0,
    disease_state  = "healthy adult male volunteers",
    dose_range     = "20 mg single oral dose of vinpocetine (either 2 x 10 mg Ultra Vinca SR beta-cyclodextrin tablets, 4 x 5 mg Cavinton IR tablets, or 10 mL of a 10 mg / 5 mL extemporaneous oral solution -- one formulation per occasion, three occasions per subject)",
    regions        = "Portugal (Hospital Pulido Valente, Lisbon)",
    notes          = "Open crossover relative-bioavailability study with a 7-day washout between formulations (Petric 2023 Section 2.2). Plasma AVA sampled between 0.25 h and 10 h post-dose. LLOQ 5 ng/mL (equivalent to 5 ug/L). BLQ observations were retained (not censored) per Section 3.3. Analytical method (HPLC-UV at 254 nm) validated in a prior publication (Petric 2023 ref [26])."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Petric 2023 Table 1 (Stochastic
    # Approximation final estimates). Monolix categorical-covariate
    # coefficients (beta) are log-additive on the individual parameter
    # for lognormal distributions, so log(Tk0_i) = log(Tk0_pop) +
    # beta_Tk0_c * I(cov = c) + eta_Tk0_i and analogously for V1/F.
    # Formulation#1 (Ultra Vinca SR beta-cyclodextrin) is the reference
    # (beta = 0). The FORM_VINP_IR indicator selects Formulation#2
    # (Cavinton IR); FORM_SOLUTION selects Formulation#3 (oral solution).
    # ------------------------------------------------------------------
    ltlag <- log(0.17);   label("Absorption lag time before zero-order input begins (h)")               # Petric 2023 Table 1 Tlag_pop = 0.17 h, RSE 7.66 %
    ld1   <- log(1.35);   label("Duration of zero-order input into central at Formulation#1 reference (h)")  # Petric 2023 Table 1 Tk0_pop = 1.35 h, RSE 8.87 %
    lcl   <- log(56.15);  label("Apparent clearance CL/F of AVA at reference covariates (L/h)")        # Petric 2023 Table 1 Cl/F_pop = 56.15 L/h, RSE 3.83 %
    lvc   <- log(208.34); label("Apparent central volume V1/F of AVA at Formulation#1 reference (L)")  # Petric 2023 Table 1 V1/F_pop = 208.34 L, RSE 5.31 %
    lq    <- log(14.63);  label("Apparent inter-compartmental clearance Q/F of AVA (L/h)")             # Petric 2023 Table 1 Q/F_pop = 14.63 L/h, RSE 11.6 % (unit typo in the printed table cell "L" corrected to "L/h" -- see vignette Errata)
    lvp   <- log(76.15);  label("Apparent peripheral volume V2/F of AVA (L)")                          # Petric 2023 Table 1 V2/F_pop = 76.15 L, RSE 14.9 %

    # Formulation covariate effects on Tk0 (log-scale, Monolix beta
    # convention). Ultra Vinca SR = reference (beta = 0).
    e_form_ir_ld1  <- -0.40; label("Log-additive effect of Cavinton immediate-release tablet on Tk0 (unitless)")  # Petric 2023 Table 1 beta_Tk0_Formulation#2 = -0.4, RSE 29.0 %
    e_form_sol_ld1 <- -0.68; label("Log-additive effect of the 10 mg / 5 mL oral solution on Tk0 (unitless)")    # Petric 2023 Table 1 beta_Tk0_Formulation#3 = -0.68, RSE 18.5 %

    # Formulation covariate effects on V1/F (log-scale, Monolix beta
    # convention). Because CL/F does not vary with formulation but V1/F
    # does, these coefficients absorb the between-formulation differences
    # in vinpocetine bioavailability and vinpocetine-to-AVA conversion.
    e_form_ir_lvc  <- -1.26; label("Log-additive effect of Cavinton immediate-release tablet on V1/F (unitless)")  # Petric 2023 Table 1 beta_V1/F_Formulation#2 = -1.26, RSE 5.44 %
    e_form_sol_lvc <- -1.24; label("Log-additive effect of the 10 mg / 5 mL oral solution on V1/F (unitless)")    # Petric 2023 Table 1 beta_V1/F_Formulation#3 = -1.24, RSE 5.71 %

    # ------------------------------------------------------------------
    # IIV -- Petric 2023 Table 1 "Standard Deviation of the Random Effects
    # (omega)". Monolix omega values ARE the SD of the underlying
    # zero-mean normal random effect on the log-scale for lognormal
    # parameters, so variance = omega^2 directly (no log(1 + CV^2)
    # transformation needed).
    #
    # Block correlation between etalcl and etalvc: Petric 2023 Table 1
    # reports corr_V1/F_Cl/F = 0.72. Log-scale covariance is therefore
    # 0.72 * 0.21 * 0.15 = 0.0227.
    # ------------------------------------------------------------------
    etaltlag         ~ 0.1444        # Petric 2023 Table 1 omega_Tlag = 0.38 -> var = 0.38^2 = 0.1444, RSE 21.4 %
    etald1           ~ 0.0576        # Petric 2023 Table 1 omega_Tk0  = 0.24 -> var = 0.24^2 = 0.0576, RSE 14.8 %
    etalcl + etalvc ~ c(0.0441,
                        0.0227,
                        0.0225)
    # Table 1: omega_Cl/F = 0.21 -> var 0.0441, RSE 13.0 %
    #          omega_V1/F = 0.15 -> var 0.0225, RSE 22.4 %
    #          corr_V1/F_Cl/F = 0.72 -> cov = 0.72 * 0.21 * 0.15 = 0.0227, RSE 19.0 %
    etalq            ~ 0.2209        # Petric 2023 Table 1 omega_Q/F  = 0.47 -> var = 0.47^2 = 0.2209, RSE 21.1 %
    etalvp           ~ 0.0576        # Petric 2023 Table 1 omega_V2/F = 0.24 -> var = 0.24^2 = 0.0576, RSE 39.1 %

    # ------------------------------------------------------------------
    # Residual error -- Petric 2023 Table 1 "Error Model Parameters". The
    # single reported parameter is "b" = 0.14 for a proportional error
    # model (Monolix "b" = proportional coefficient acting on the
    # individual prediction f), i.e., y = f * (1 + b * eps).
    # ------------------------------------------------------------------
    propSd <- 0.14;  label("Proportional residual error on plasma AVA Cc (fraction)")  # Petric 2023 Table 1 b = 0.14, RSE 5.28 %
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual PK parameters. Formulation effects enter log-
    # additively on Tk0 and V1/F only (per Petric 2023 Section 3.3, no
    # continuous covariate was retained; Formulation was retained only
    # on Tk0 and V1/F).
    # ------------------------------------------------------------------
    tlag <- exp(ltlag + etaltlag)
    d1   <- exp(ld1   + etald1  +
                e_form_ir_ld1  * FORM_VINP_IR +
                e_form_sol_ld1 * FORM_SOLUTION)
    cl   <- exp(lcl   + etalcl)
    vc   <- exp(lvc   + etalvc  +
                e_form_ir_lvc  * FORM_VINP_IR +
                e_form_sol_lvc * FORM_SOLUTION)
    q    <- exp(lq    + etalq)
    vp   <- exp(lvp   + etalvp)

    # 2. Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ------------------------------------------------------------------
    # 3. ODE system -- two-compartment with zero-order input into
    # central. The dose event carries the vinpocetine amount in mg
    # (target compartment = central, rate = -2 with dur() supplied by
    # the model, or an explicit rate = amount / d1 in the event table).
    # An absorption lag Tlag delays the onset of the zero-order input.
    # ------------------------------------------------------------------
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    dur(central)  <- d1
    alag(central) <- tlag

    # ------------------------------------------------------------------
    # 4. Observation. Dose is in mg (vinpocetine); central is in mg;
    # V1/F is in L. central / vc gives mg/L; multiply by 1000 to display
    # in ug/L (= ng/mL), matching Petric 2023 Figures 2, 4, 5 and the
    # 5 ng/mL LLOQ (5 ug/L).
    # ------------------------------------------------------------------
    Cc <- central / vc * 1000

    Cc ~ prop(propSd)
  })
}
