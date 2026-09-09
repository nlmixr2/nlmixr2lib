DiDeo_2025_tideglusib_dm1 <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption and",
    "dose-dependent bioavailability for tideglusib in adolescent and adult",
    "patients with congenital or juvenile-onset myotonic dystrophy type 1",
    "(Phase II study AMO-02-MD-2-001), estimated using the elderly healthy",
    "subject model as a Bayesian prior"
  )
  reference <- paste(
    "Di Deo A, Oosterholt S, Horrigan J, Evans S, McMorn A, Della Pasqua O.",
    "Population Pharmacokinetics of Tideglusib in Congenital and Childhood",
    "Myotonic Dystrophy Type 1: Influence of Demographic and Clinical Factors",
    "on Systemic Exposure. Pharmaceutics. 2025;17(8):1065.",
    "doi:10.3390/pharmaceutics17081065.",
    "Estimated with NONMEM $PRIOR NWPRI using the healthy-elderly-subject",
    "model of the same paper as the prior; see",
    "modellib('DiDeo_2025_tideglusib_healthy')."
  )
  vignette <- "DiDeo_2025_tideglusib"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")
  # NOTE ON CONCENTRATION UNITS: the paper reports concentrations and NCA
  # metrics in ng/mL. This model is parameterised in the paper's own parameter
  # units (dose mg, CL L/h, V L), so central/vc is natively mg/L == ug/mL.
  # No conversion factor is hardcoded in model(); multiply Cc by 1000 to
  # compare against the paper's ng/mL values (see the vignette).

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling with a reference body weight of 70 kg and FIXED",
        "exponents, per Di Deo 2025 Section 2.4.1: 'This effect was",
        "parameterised using a reference body weight of 70 kg and fixed",
        "allometric exponents on clearance (0.75) and volumes (1).' Applied to",
        "CL, Q, V2 and V3 -- Figure 1 caption names 'the effect of body weight",
        "as a covariate on clearance (CL), intercompartmental clearance (Q),",
        "and volumes of distribution (V2 and V3)', and Section 3.1 states for",
        "the patient model that 'Weight was the only covariate to have a",
        "significant effect on the disposition of tideglusib.' Q takes the",
        "clearance exponent (0.75) because it is a clearance term. Body weight",
        "in this cohort ranged 36.8-122.6 kg (mean 63.6, s.d. 19.6; Section 2.2",
        "and Table 1). Section 2.7 states that NO maturation function and no",
        "age-dependent allometric exponent were included -- the Discussion",
        "limitation reads 'as there is incomplete data on the primary route of",
        "metabolism in vivo, neither a maturation function nor an",
        "age-dependent allometric exponent were included in the simulated",
        "scenarios' -- so weight is the sole size/age descriptor even though",
        "the model is used to extrapolate down to 5 kg."
      ),
      source_name        = "WT"
    ),
    DOSE_HIGH = list(
      description        = paste(
        "High-dose indicator. 1 = the subject's administered tideglusib dose",
        "was greater than 400 mg; 0 = 400 mg or less."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Di Deo 2025 Table 2 splits BOTH the absorption rate constant and the",
        "relative bioavailability at a 400 mg threshold: rows 'Absorption rate",
        "constant, kA (1/h)' vs 'Absorption rate constant, kA (dose > 400 mg)",
        "(1/h)', and 'Bioavailability (doses <= 400 mg)' (FIXED to 1) vs",
        "'Bioavailability (doses > 400 mg)'. The residual error is split on the",
        "SAME threshold, so this one indicator gates three quantities.",
        "Section 3.1 for the patient model: 'dose (strength) was found to have",
        "a statistically significant effect on residual variability in the",
        "patient population, with different parameter estimates for doses",
        "higher than 400 mg.' Discussion: 'bioavailability and absorption rate",
        "constant were dichotomized into two levels, with exposure data",
        "following a 400 mg dose as reference'. In THIS study only two dose",
        "levels were used (Section 2.1): the 400 mg q.d. arm is DOSE_HIGH = 0",
        "and the 1000 mg q.d. arm is DOSE_HIGH = 1, 8 subjects each.",
        "Time-fixed per subject: subjects were assigned one dose for the full",
        "12-week active phase. When simulating the Section 2.7 paediatric",
        "weight-banded regimen, set DOSE_HIGH from the BAND's dose (only the",
        "60+ kg band's 1000 mg and the 45-60 kg band's 800 mg and the",
        "35-45 kg band's 600 mg exceed 400 mg)."
      ),
      source_name        = "derived from the assigned dose level"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened but NOT retained. Di Deo 2025 Section 2.4.1 lists age among",
        "the factors considered and notes that where collinearity exists",
        "between covariates, such as age and weight, the choice was guided by",
        "which factor is most likely the primary cause of variability.",
        "Section 3.1: 'Weight was the only covariate to have a significant",
        "effect on the disposition of tideglusib.' Range in this cohort",
        "13.8-34.9 years (mean 20.9, s.d. 5.9)."
      )
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Screened but NOT retained; named in the Section 2.4.1 covariate list.",
        "Body weight was the only retained covariate."
      )
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Reported in Table 1 (6 of 16 patients female) but not retained as a",
        "covariate in the final model."
      )
    ),
    FED = list(
      description = "Fed-vs-fasted indicator at the dose record",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened but NOT retained, and NOT identifiable as a fed-vs-fasted",
        "contrast in this study: Section 2.1 states every dose was given after",
        "an overnight fast with food restricted for at least one hour",
        "afterwards, so all records are fasted. Sections 2.6 and 3.3 instead",
        "explore the TIME of the first post-dose meal and the meal TYPE",
        "(light / standard / other) against derived AUC and Cmax in 30",
        "records, and find no effect: 'There was no evidence of an effect of",
        "food or meal type on the predicted exposure to tideglusib.' That",
        "analysis is post hoc on NCA metrics, not a covariate in the model,",
        "which is why no meal term appears in model(). Section 2.7 states the",
        "dose-optimisation simulations assume all individuals are fasted."
      )
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "tideglusib", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "tideglusib", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "tideglusib", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 16,
    n_studies      = 1,
    age_range      = "13.8-34.9 years",
    age_mean       = "20.9 years (s.d. 5.9)",
    weight_range   = "36.8-122.6 kg",
    weight_mean    = "63.6 kg (s.d. 19.6)",
    sex_female_pct = 37.5,
    disease_state  = "congenital or juvenile-onset myotonic dystrophy type 1 (DM-1)",
    dose_range     = "400 mg or 1000 mg q.d. oral, 12 weeks (1:1 allocation, 8 subjects per arm)",
    regions        = "United Kingdom (single centre, Newcastle upon Tyne Hospitals NHS Trust)",
    notes          = paste(
      "Study AMO-02-MD-2-001, a Phase II single-blind, placebo-controlled,",
      "fixed-dose study in adolescent and adult patients aged 13-34 years.",
      "All subjects received 2 weeks of single-blind placebo before the",
      "12-week active phase. 51 evaluable plasma samples were available; 9",
      "were below the LLOQ of 1 ng/mL (2 with no peak) and all 9 were pre-dose",
      "trough measurements. Sparse sampling: pre-dose and between 2 and 4 h",
      "post-dose at weeks 2 and 12. Because of the sparse design, parameters",
      "were estimated with NONMEM $PRIOR NWPRI using the healthy-elderly model",
      "as an informative prior (Section 2.4.2). Baseline demographics are",
      "Di Deo 2025 Table 1 (last two rows) and Section 2.2. Section 3.1 notes",
      "the disposition parameters did not differ meaningfully from the healthy",
      "elderly subjects."
    )
  )

  ini({
    # ----- Structural parameters (Di Deo 2025 Table 2, AMO-02-MD-2-001 column) -----
    # NOTE ON TABLE 2 COLUMN ORDER: Table 2's title reads "in healthy elderly
    # subjects and DM-1 patients" but its COLUMNS are ordered AMO-02-MD-2-001
    # (DM-1 patients) FIRST, then NP031112-07A03 (healthy elderly). The column
    # headers are correct as printed; the title order is not the column order.
    # Confirmed against the paper's own precision claims: Section 3.1 reports
    # the DM-1 fit as "%RSE < 23%" for fixed effects and "%RSE < 32%" for IIV,
    # which bound the AMO-02-MD-2-001 RSE column maxima exactly (22.1 and
    # 31.6); the healthy-subject claims ("< 27.4%", "< 41.1%") are the
    # NP031112-07A03 column maxima. This file uses the AMO-02-MD-2-001 column
    # throughout. Independent corroboration: the Discussion's headline numbers
    # are from this column -- "apparent systemic clearance of 341 L/h" and an
    # apparent steady-state volume of "1140 L", which is V2 + V3 = 154 + 986.
    #
    # All disposition parameters are APPARENT (CL/F, V/F): bioavailability is
    # FIXED to 1 at doses <= 400 mg, which anchors the scale.
    lka <- log(0.767); label("Absorption rate constant, doses <= 400 mg (1/h)")  # Table 2, AMO-02-MD-2-001, row "Absorption rate constant, kA (1/h)" = 0.767
    lka_highdose <- log(0.609); label("Absorption rate constant, doses > 400 mg (1/h)")  # Table 2, AMO-02-MD-2-001, row "Absorption rate constant, kA (dose > 400 mg) (1/h)" = 0.609
    lcl <- log(341); label("Apparent clearance CL/F at 70 kg (L/h)")  # Table 2, AMO-02-MD-2-001, row "Clearance, CL (L/h)" = 341
    lvc <- log(154); label("Apparent central volume V2/F at 70 kg (L)")  # Table 2, AMO-02-MD-2-001, row "Volume of distribution central compartment, V2 (L)" = 154
    lq <- log(66.1); label("Apparent intercompartmental clearance Q/F at 70 kg (L/h)")  # Table 2, AMO-02-MD-2-001, row "Intercompartmental clearance, Q (L/h)" = 66.1
    lvp <- log(986); label("Apparent peripheral volume V3/F at 70 kg (L)")  # Table 2, AMO-02-MD-2-001, row "Volume of distribution peripheral compartment, V3 (L)" = 986

    # ----- Dose-dependent relative bioavailability -----
    lfdepot <- fixed(log(1)); label("Relative bioavailability, doses <= 400 mg (fraction)")  # Table 2, row "Bioavailability (doses <= 400 mg)" = 1, reported as FIXED
    lfdepot_highdose <- log(0.88); label("Relative bioavailability, doses > 400 mg (fraction)")  # Table 2, AMO-02-MD-2-001, row "Bioavailability (doses > 400 mg)" = 0.88

    # ----- Allometric exponents (fixed, not estimated) -----
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)")  # Section 2.4.1: "fixed allometric exponents on clearance (0.75) and volumes (1)"
    e_wt_vc_vp <- fixed(1); label("Allometric exponent on V2 and V3 (unitless)")  # Section 2.4.1: "fixed allometric exponents on clearance (0.75) and volumes (1)"

    # ----- Interindividual variability -----
    # Table 2 tabulates VARIANCES, not standard deviations. Proof: this
    # column's variances are printed to 3 significant figures and satisfy
    # CV% = sqrt(variance) * 100 to the digit for all five terms --
    # sqrt(0.186) = 43.1, sqrt(1.15) = 107.2, sqrt(1.01) = 100.5,
    # sqrt(0.321) = 56.7, sqrt(0.069) = 26.3, against the printed 43.1, 107.2,
    # 100.5, 56.7 and 26.4. The log-normal alternative
    # sqrt(exp(variance) - 1) * 100 gives 45.2, 146.9, 132.1, 61.5 and 26.7 and
    # is excluded. (The table's footnote b prints the formula as
    # "CV = sqrt(exp(Omega^2)) x 100", which does not reproduce its own column;
    # the arithmetic that does is sqrt(variance).) nlmixr2 takes variances on
    # the diagonal, so the printed values are used directly.
    etalcl ~ 0.186  # Table 2, AMO-02-MD-2-001, row "eta CL variance" = 0.186 (43.1% CV)
    etalvc ~ 0.321  # Table 2, AMO-02-MD-2-001, row "eta V2 variance" = 0.321 (56.7% CV)
    etalq ~ 1.15  # Table 2, AMO-02-MD-2-001, row "eta Q variance" = 1.15 (107.2% CV)
    etalvp ~ 1.01  # Table 2, AMO-02-MD-2-001, row "eta V3 variance" = 1.01 (100.5% CV)
    etalka ~ 0.069  # Table 2, AMO-02-MD-2-001, row "eta KA variance" = 0.069 (26.4% CV)

    # ----- Residual error (proportional, split on the same 400 mg threshold) -----
    # Table 2's abbreviation footnote defines "sigma = residual variance", and
    # the IIV rows of the same table are demonstrably variances (above), so the
    # printed residual values are variances too and are converted to the
    # standard deviations nlmixr2 expects. sqrt() is kept inline so the printed
    # value stays visible in the source trace.
    propSd <- sqrt(0.54); label("Proportional residual error, doses <= 400 mg (fraction)")  # Table 2, AMO-02-MD-2-001, row "Proportional error (doses <= 400 mg)" = 0.54 (variance)
    propSd_highdose <- sqrt(0.46); label("Proportional residual error, doses > 400 mg (fraction)")  # Table 2, AMO-02-MD-2-001, row "Proportional error (doses > 400 mg)" = 0.46 (variance)
  })

  model({
    # 1. Individual parameters. The dose switch uses the two-typical-value form
    #    lx + (lx_highdose - lx) * DOSE_HIGH, which returns exactly exp(lx) at
    #    DOSE_HIGH = 0 and exp(lx_highdose) at DOSE_HIGH = 1, so both printed
    #    Table 2 values appear literally in ini() (the Kastrissios_2006_apricoxib.R
    #    precedent). A single IIV term applies to ka across both dose levels,
    #    matching the single "eta KA variance" row in Table 2.
    ka <- exp(lka + (lka_highdose - lka) * DOSE_HIGH + etalka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q <- exp(lq + etalq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODE system
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Dose-dependent relative bioavailability, same switch form as ka.
    f(depot) <- exp(lfdepot + (lfdepot_highdose - lfdepot) * DOSE_HIGH)

    # 5. Observation and error. The proportional residual error is itself
    #    dose-dependent (Table 2 residual-error rows), so the error magnitude is
    #    assembled as a model variable and passed to prop().
    Cc <- central / vc
    propSdDose <- propSd + (propSd_highdose - propSd) * DOSE_HIGH
    Cc ~ prop(propSdDose)
  })
}
