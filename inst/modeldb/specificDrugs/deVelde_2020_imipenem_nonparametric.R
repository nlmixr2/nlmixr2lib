deVelde_2020_imipenem_nonparametric <- function() {
  description <- paste(
    "Two-compartment IV population PK model for imipenem in 26 critically",
    "ill adults treated with imipenem-cilastatin in a Geneva intensive care",
    "unit (de Velde 2020, non-parametric Pmetrics NPAG arm), parameterised",
    "in micro-constants: an elimination rate constant, two distribution rate",
    "constants and a central volume. The elimination rate constant scales as",
    "a power of absolute (BSA-unadjusted) CKD-EPI eGFR in mL/min. NPAG places",
    "every parameter, including the covariate exponent, in a discrete joint",
    "density; that density is approximated here by independent log-normal",
    "marginals matched to the published mean and CV%. Residual error is the",
    "Pmetrics gamma-scaled assay polynomial, a linear sum of additive and",
    "proportional terms. The parametric NONMEM arm fitted to the same data",
    "is the sibling model deVelde_2020_imipenem.",
    sep = " "
  )
  reference <- paste(
    "de Velde F, de Winter BCM, Neely MN, Yamada WM, Koch BCP,",
    "Harbarth S, von Dach E, van Gelder T, Huttner A, Mouton JW, on behalf",
    "of COMBACTE-NET consortium.",
    "Population pharmacokinetics of imipenem in critically ill patients: a",
    "parametric and nonparametric model converge on CKD-EPI estimated",
    "glomerular filtration rate as an impactful covariate.",
    "Clin Pharmacokinet. 2020;59(7):885-898. doi:10.1007/s40262-020-00859-1",
    sep = " "
  )
  vignette <- "deVelde_2020_imipenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Glomerular filtration rate estimated by the CKD-EPI equation and",
        "then DE-NORMALISED to an absolute per-patient rate in mL/min by",
        "multiplying the BSA-normalised value by the patient's body surface",
        "area (de Velde 2020 Sect. 2.5 and Table 1 footnote: 'CKD-EPI-abs",
        "absolute CKD-EPI (i.e. CKD-EPI multiplied by BSA)')."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "UNIT HAZARD. This column is an ABSOLUTE rate in mL/min, NOT the",
        "mL/min/1.73 m^2 that a CKD-EPI calculator returns by default and",
        "that the CRCL canonical otherwise carries. The reference of 119",
        "mL/min is the cohort median CKD-EPI-abs at inclusion (Table 1).",
        "Enters the elimination rate constant as (CRCL/119)^Ke(cov) (Eq.",
        "10). TIME-VARYING in the source analysis; Pmetrics applies",
        "covariates at each dose event with LOCF and linear interpolation",
        "(Sect. 2.5). MDRD, MDRD-abs, Jelliffe-abs and CKD-EPI-abs were",
        "statistically indistinguishable on -2LL; CKD-EPI-abs was chosen",
        "for the lowest bias and imprecision (Sect. 3.4)."
      ),
      source_name = "CKD-EPI-abs"
    )
  )

  # Univariately significant but removed at backward elimination (Sect. 3.4:
  # 'TBW, IBW and LBW on Ke ... and CG, Jelliffe and Jelliffe-abs on V ...
  # After backward elimination at p < 0.001, none of these six covariates
  # remained in the final model').
  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Significant on Ke univariately (d-2LL 4.0-8.8 across TBW/IBW/LBW), removed at backward elimination (Sect. 3.4). Cohort median 75 kg (Table 1)."
    ),
    IBW = list(
      description = "Ideal body weight",
      units = "kg",
      type = "continuous",
      notes = "Significant on Ke univariately, removed at backward elimination (Sect. 3.4)."
    ),
    LBW = list(
      description = "Lean body weight",
      units = "kg",
      type = "continuous",
      notes = "Significant on Ke univariately, removed at backward elimination (Sect. 3.4). Not a registered canonical in inst/references/covariate-columns.md; recorded here as documentation only, since it is never referenced in model()."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 26L,
    n_studies = 1L,
    age_median = "51 years (IQR 39-54); inclusion 18-60 years",
    weight_median = "75 kg (IQR 66-85)",
    sex_female_pct = 30.8,
    race_ethnicity = NULL,
    disease_state = paste(
      "Critically ill adults with suspected or documented severe bacterial",
      "infection (lower respiratory tract 62%, intra-abdominal 15%,",
      "bloodstream 12%, other 12%); APACHE II median 22 (IQR 17-27); no",
      "continuous renal replacement therapy. Exclusion: Cockcroft-Gault",
      "eGFR < 60 mL/min, BMI < 18 or > 30 kg/m^2, pregnancy."
    ),
    dose_range = paste(
      "Imipenem/cilastatin 500 mg/500 mg four times daily as a 30-min",
      "intermittent IV infusion (Sect. 2.1)."
    ),
    regions = "Switzerland (Geneva University Hospitals ICU)",
    renal_function = "CKD-EPI 116 mL/min/1.73 m^2 (IQR 104-124); CKD-EPI-abs 119 mL/min (IQR 110-139); Cockcroft-Gault 146 mL/min (IQR 123-170) at inclusion (Table 1).",
    n_concentrations = 125L,
    notes = paste(
      "Same 26 patients and 125 above-LOQ concentrations as the parametric",
      "sibling (Sect. 3.2). Pmetrics 1.5.2 NPAG on untransformed",
      "concentrations; IT2B supplied the parameter search ranges (Ke 0-1.5",
      "1/h, V 1-70 L, Kcp and Kpc 0-1 1/h; Sect. 3.4). The final",
      "population is a discrete distribution of 16 support points (Fig. 1).",
      "The unbounded log-normal marginals used here ignore those search",
      "ranges, so simulated tails are wider than the NPAG distribution's:",
      "the paper reports a roughly twofold lower 2.5th percentile than the",
      "parametric model, and this approximation gives roughly three- to",
      "fivefold. Medians are unaffected."
    )
  )

  ini({
    # ===== Structural PK -- de Velde 2020 Table 2, Pmetrics 'Final model',
    # 'Mean parameter estimate' column: the probability-weighted mean of the
    # NPAG support points (Sect. 2.4). Encoded as the typical value (the
    # median of each log-normal marginal), following the repository's
    # precedent for Pmetrics means (Tsai_2023_ceftriaxone.R). A log-normal
    # whose median is the printed mean has a mean exp(omega^2/2) higher
    # (6% for Ke, 29% for Kcp); see the vignette for this approximation. =====
    lkel <- log(0.681); label("Elimination rate constant at CRCL = 119 mL/min (1/h)")  # Table 2 Pmetrics: Ke mean = 0.681 1/h (bootstrap median 0.586, 95% CI 0.533-0.905)
    lk12 <- log(0.374); label("Central-to-peripheral rate constant Kcp (1/h)")         # Table 2 Pmetrics: Kcp mean = 0.374 1/h (bootstrap median 0.347, 95% CI 0.122-0.563)
    lk21 <- log(0.495); label("Peripheral-to-central rate constant Kpc (1/h)")         # Table 2 Pmetrics: Kpc mean = 0.495 1/h (bootstrap median 0.387, 95% CI 0.278-0.846)
    lvc  <- log(31.1);  label("Central volume of distribution (L)")                    # Table 2 Pmetrics: Vc mean = 31.1 L (bootstrap median 35.1, 95% CI 20.1-38.3)

    # ===== Covariate effect -- Eq. 10:
    #   Ke_i = Ke_i,med x (CKD-EPI-abs_i / 119)^Ke(cov)_i,med
    # Ke(cov) is itself a random parameter in NPAG, so it carries its own
    # eta below (log-normal on the coefficient, as Downes_2023_vancomycin_full
    # does for its Pmetrics covariate slope). =====
    e_crcl_kel <- 0.658; label("Power exponent on (CRCL/119) for kel (unitless)")  # Table 2 Pmetrics: Ke(cov) mean = 0.658 (bootstrap median 0.791, 95% CI 0.516-1.000)

    # ===== Between-subject variability. The Pmetrics CV% is SD/mean of the
    # discrete support-point distribution (Sect. 2.4), a genuine descriptive
    # CV, so a log-normal marginal matches it exactly through
    # omega^2 = log(CV^2 + 1). Correlations between parameters are not
    # reported (only 'No large correlation (> 0.95)', Sect. 3.4), so the
    # marginals are independent. =====
    etalkel      ~ 0.109392  # Table 2 Pmetrics: Ke CV 34.0%; log(0.340^2 + 1)
    etalk12      ~ 0.506422  # Table 2 Pmetrics: Kcp CV 81.2%; log(0.812^2 + 1)
    etalk21      ~ 0.417657  # Table 2 Pmetrics: Kpc CV 72.0%; log(0.720^2 + 1)
    etalvc       ~ 0.166765  # Table 2 Pmetrics: Vc CV 42.6%; log(0.426^2 + 1)
    etae_crcl_kel ~ 0.265976 # Table 2 Pmetrics: Ke(cov) CV 55.2%; log(0.552^2 + 1)

    # ===== Residual error -- gamma model, Eqs. 3 and 5:
    #   error = gamma x (C0 + C1 x OBS), with C2 = C3 = 0
    # Final C0 = C1 = 0.05 (Sect. 3.4) and gamma = 3.40 (Table 2), so
    #   SD = 3.40 x 0.05 + 3.40 x 0.05 x C = 0.17 mg/L + 0.17 x C.
    # The two terms add linearly, which is nlmixr2's combined1() form. Eq. 5
    # evaluates the polynomial at the observation; here it is evaluated at
    # the prediction, as nlmixr2 requires. =====
    addSd  <- 0.17; label("Additive residual SD (mg/L)")         # gamma x C0 = 3.40 x 0.05 (Table 2; Sect. 3.4)
    propSd <- 0.17; label("Proportional residual SD (fraction)") # gamma x C1 = 3.40 x 0.05 (Table 2; Sect. 3.4)
  })

  model({
    # ----- Individual PK parameters (Eq. 10) -----
    e_crcl_kel_i <- e_crcl_kel * exp(etae_crcl_kel)
    kel <- exp(lkel + etalkel) * (CRCL / 119)^e_crcl_kel_i
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    vc  <- exp(lvc + etalvc)

    # ----- ODE system -----
    # Imipenem-cilastatin given as a 30-min IV infusion into the central
    # compartment; the infusion duration comes from the event table's
    # rate / dur column.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # ----- Output -----
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
