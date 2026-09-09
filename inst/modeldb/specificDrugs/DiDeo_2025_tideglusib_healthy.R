DiDeo_2025_tideglusib_healthy <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption and",
    "dose-dependent bioavailability for tideglusib in elderly healthy",
    "subjects (Phase I study NP031112-07A03)"
  )
  reference <- paste(
    "Di Deo A, Oosterholt S, Horrigan J, Evans S, McMorn A, Della Pasqua O.",
    "Population Pharmacokinetics of Tideglusib in Congenital and Childhood",
    "Myotonic Dystrophy Type 1: Influence of Demographic and Clinical Factors",
    "on Systemic Exposure. Pharmaceutics. 2025;17(8):1065.",
    "doi:10.3390/pharmaceutics17081065"
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
        "and volumes of distribution (V2 and V3)', and Section 3.1 states",
        "'Weight was the only covariate to have a significant effect on CL,",
        "V2, V3 and Q.' Q takes the clearance exponent (0.75) because it is a",
        "clearance term. Body weight in this cohort ranged 50.7-98.1 kg",
        "(mean 74.5, s.d. 9.9; Section 2.2 and Table 1)."
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
        "SAME threshold ('Proportional error (doses <= 400 mg)' vs '(doses >",
        "400 mg)'), so this one indicator gates three quantities. Section 3.1:",
        "'a separate parameter was required for doses lower and higher than",
        "400 mg'. Discussion: 'bioavailability and absorption rate constant",
        "were dichotomized into two levels, with exposure data following a",
        "400 mg dose as reference'. The threshold is on the PER-ADMINISTRATION",
        "dose, not the daily dose: the 400 mg b.i.d. arm (800 mg/day) is",
        "DOSE_HIGH = 0 because each administration was 400 mg. In this study",
        "the 300 mg b.i.d. and 400 mg b.i.d. arms are DOSE_HIGH = 0 and the",
        "600, 800, 1000 and 1200 mg q.d. arms are DOSE_HIGH = 1 (Section 2.1",
        "dose groups). Time-fixed per subject: each subject remained on a",
        "single dose level for the whole study."
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
        "the factors considered ('In addition to body weight, several factors",
        "were considered (e.g., age, BMI, dose, and food intake...)') and",
        "states that where collinearity exists between covariates, such as age",
        "and weight, the choice was guided by which factor is most likely the",
        "primary cause of variability. Section 3.1: 'Weight was the only",
        "covariate to have a significant effect on CL, V2, V3 and Q.' Range in",
        "this cohort 60-74 years (mean 64.3, s.d. 3.5)."
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
        "Reported in Table 1 (24 of 54 subjects female) but not retained as a",
        "covariate in the final model."
      )
    ),
    FED = list(
      description = "Fed-vs-fasted indicator at the dose record",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened but NOT retained. Section 2.4.1 lists 'food intake, i.e.,",
        "fasting or fed state' among the candidate covariates. The formal food",
        "effect analysis in this paper was performed only in the DM-1 patients",
        "(Section 2.6/3.3) and found no effect; the Discussion notes a marked",
        "food effect had been seen with the EARLIER F05-052 formulation in",
        "studies CL031112-05A01 and CL031112-06A02, which are NOT part of this",
        "model's dataset (Section 2.4.1 restricts the analysis to study",
        "NP031112-07A03, formulation F06-037F)."
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
    n_subjects     = 54,
    n_studies      = 1,
    age_range      = "60-74 years",
    age_mean       = "64.3 years (s.d. 3.5)",
    weight_range   = "50.7-98.1 kg",
    weight_mean    = "74.5 kg (s.d. 9.9)",
    sex_female_pct = 44.4,
    disease_state  = "elderly healthy volunteers",
    dose_range     = paste(
      "300 mg b.i.d., 400 mg b.i.d., 600 mg q.d., 800 mg q.d., 1000 mg q.d.",
      "and 1200 mg q.d. oral suspension (formulation F06-037F); single dose,",
      "48 h washout, then 14 days of repeat dosing"
    ),
    regions        = "Germany (single centre, Parexel International GmbH, Berlin)",
    notes          = paste(
      "Study NP031112-07A03, a double-blind, randomised, parallel-group,",
      "placebo-controlled multiple-ascending-dose study. 72 subjects were",
      "randomised (12 per dose group, 9 active : 3 placebo); the 54 who",
      "received active tideglusib contributed 1832 plasma samples to the",
      "pharmacokinetic analysis. Baseline demographics are Di Deo 2025 Table 1",
      "(first six rows) and Section 2.2. Sampling on days 1 and 16 at pre-dose,",
      "20 and 40 min, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 12, 16, 24, 30, 36 and 48 h",
      "post-dose, plus troughs on days 5, 7, 9, 11, 13 and 15."
    )
  )

  ini({
    # ----- Structural parameters (Di Deo 2025 Table 2, NP031112-07A03 column) -----
    # NOTE ON TABLE 2 COLUMN ORDER: Table 2's title reads "in healthy elderly
    # subjects and DM-1 patients" but its COLUMNS are ordered AMO-02-MD-2-001
    # (DM-1 patients) FIRST, then NP031112-07A03 (healthy elderly). The column
    # headers are correct as printed; the title order is not the column order.
    # Confirmed against the paper's own precision claims: Section 3.1 reports
    # the healthy-subject fit as "%RSE < 27.4%" for fixed effects and "%RSE <
    # 41.1%" for IIV, which are exactly the maxima of the NP031112-07A03 RSE
    # column (27.4 and 41.1); the DM-1 claims ("<23%", "<32%") match the
    # AMO-02-MD-2-001 column maxima (22.1 and 31.6). This file uses the
    # NP031112-07A03 column throughout.
    #
    # All disposition parameters are APPARENT (CL/F, V/F): bioavailability is
    # FIXED to 1 at doses <= 400 mg, which anchors the scale.
    lka <- log(0.78); label("Absorption rate constant, doses <= 400 mg (1/h)")  # Table 2, NP031112-07A03, row "Absorption rate constant, kA (1/h)" = 0.78
    lka_highdose <- log(0.61); label("Absorption rate constant, doses > 400 mg (1/h)")  # Table 2, NP031112-07A03, row "Absorption rate constant, kA (dose > 400 mg) (1/h)" = 0.61
    lcl <- log(327); label("Apparent clearance CL/F at 70 kg (L/h)")  # Table 2, NP031112-07A03, row "Clearance, CL (L/h)" = 327
    lvc <- log(152); label("Apparent central volume V2/F at 70 kg (L)")  # Table 2, NP031112-07A03, row "Volume of distribution central compartment, V2 (L)" = 152
    lq <- log(66.1); label("Apparent intercompartmental clearance Q/F at 70 kg (L/h)")  # Table 2, NP031112-07A03, row "Intercompartmental clearance, Q (L/h)" = 66.1
    lvp <- log(1010); label("Apparent peripheral volume V3/F at 70 kg (L)")  # Table 2, NP031112-07A03, row "Volume of distribution peripheral compartment, V3 (L)" = 1010

    # ----- Dose-dependent relative bioavailability -----
    lfdepot <- fixed(log(1)); label("Relative bioavailability, doses <= 400 mg (fraction)")  # Table 2, row "Bioavailability (doses <= 400 mg)" = 1, reported as FIXED
    lfdepot_highdose <- log(0.85); label("Relative bioavailability, doses > 400 mg (fraction)")  # Table 2, NP031112-07A03, row "Bioavailability (doses > 400 mg)" = 0.85

    # ----- Allometric exponents (fixed, not estimated) -----
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)")  # Section 2.4.1: "fixed allometric exponents on clearance (0.75) and volumes (1)"
    e_wt_vc_vp <- fixed(1); label("Allometric exponent on V2 and V3 (unitless)")  # Section 2.4.1: "fixed allometric exponents on clearance (0.75) and volumes (1)"

    # ----- Interindividual variability -----
    # Table 2 tabulates VARIANCES, not standard deviations. The proof is in the
    # ADJACENT AMO-02-MD-2-001 column, whose variances are printed to 3
    # significant figures: CV% = sqrt(variance) * 100 reproduces every one of
    # its five printed CV% values to the digit (sqrt(0.186) = 43.1,
    # sqrt(1.15) = 107.2, sqrt(1.01) = 100.5, sqrt(0.321) = 56.7,
    # sqrt(0.069) = 26.3 vs printed 43.1, 107.2, 100.5, 56.7, 26.4), whereas
    # the log-normal form sqrt(exp(variance) - 1) * 100 gives 45.2, 146.9,
    # 132.1, 61.5, 26.7 and does not. The NP031112-07A03 variances used here
    # are printed to only 2 significant figures, so they reproduce their own
    # CV% column to rounding rather than exactly (sqrt(0.18) = 42.4 vs printed
    # 43.5, i.e. the unrounded variance is nearer 0.189). (Footnote b prints
    # the formula as "CV = sqrt(exp(Omega^2)) x 100", which does not reproduce
    # its own column; the arithmetic that does is sqrt(variance).) nlmixr2
    # takes variances on the diagonal, so the printed values are used directly.
    etalcl ~ 0.18  # Table 2, NP031112-07A03, row "eta CL variance" = 0.18 (43.5% CV)
    etalvc ~ 0.33  # Table 2, NP031112-07A03, row "eta V2 variance" = 0.33 (57.8% CV)
    etalq ~ 0.96  # Table 2, NP031112-07A03, row "eta Q variance" = 0.96 (98.2% CV)
    etalvp ~ 0.90  # Table 2, NP031112-07A03, row "eta V3 variance" = 0.90 (95.1% CV)
    etalka ~ 0.06  # Table 2, NP031112-07A03, row "eta KA variance" = 0.06 (26.2% CV)

    # ----- Residual error (proportional, split on the same 400 mg threshold) -----
    # Table 2's abbreviation footnote defines "sigma = residual variance", and
    # the IIV rows of the same table are demonstrably variances (above), so the
    # printed residual values are variances too and are converted to the
    # standard deviations nlmixr2 expects. sqrt() is kept inline so the printed
    # value stays visible in the source trace.
    propSd <- sqrt(0.56); label("Proportional residual error, doses <= 400 mg (fraction)")  # Table 2, NP031112-07A03, row "Proportional error (doses <= 400 mg)" = 0.56 (variance)
    propSd_highdose <- sqrt(0.45); label("Proportional residual error, doses > 400 mg (fraction)")  # Table 2, NP031112-07A03, row "Proportional error (doses > 400 mg)" = 0.45 (variance)
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
