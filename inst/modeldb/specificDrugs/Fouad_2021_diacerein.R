Fouad_2021_diacerein <- function() {
  description <- paste(
    "PBPK reduction (Simcyp version 17, minimal PBPK with ADAM absorption).",
    "One-compartment oral pharmacokinetic model of diacerein in healthy",
    "adults, reduced from the Simcyp model Fouad 2021 used to compare plain",
    "crystalline diacerein with an optimised PEG 8000 solid dispersion",
    "(drug:polymer 1:4 w/w). The source model has no single adjusting",
    "compartment, a predicted Vss of 0.093 L/kg and an in-vivo oral",
    "clearance input of 1.5 L/h (30% CV), with hepatic and gut extraction",
    "both equal to 1, so the plasma profile is one-compartmental. Absorption",
    "is driven by the formulation's in-vitro dissolution profile (released",
    "drug per hour, linearly interpolated, held flat after 1 h) followed by",
    "first-order absorption of dissolved drug at the Simcyp-predicted",
    "first-order equivalent ka (2.098 1/h, fraction absorbed 0.992). No",
    "parameter is fitted. The reduction reproduces the source's simulated",
    "Cmax and AUC0-24 for both formulations to within about 6%, but peaks",
    "earlier because gastric emptying and small-intestinal transit are not",
    "represented. The model is valid for single doses and for repeat doses",
    "given at least 1 h apart. The geriatric simulation used a different",
    "(mechanistic diffusion-layer) dissolution model and Simcyp geriatric",
    "physiology and is not reproducible from this model.",
    sep = " "
  )
  reference <- paste(
    "Fouad SA, Malaak FA, El-Nabarawi MA, Abu Zeid K, Ghoneim AM. (2021).",
    "Preparation of solid dispersion systems for enhanced dissolution of",
    "poorly water soluble diacerein: In-vitro evaluation, optimization and",
    "physiologically based pharmacokinetic modeling.",
    "PLoS ONE 16(1):e0245482. doi:10.1371/journal.pone.0245482.",
    sep = " "
  )
  vignette <- "Fouad_2021_diacerein"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(
      analyte = "diacerein",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    gut_lumen = list(
      analyte = "diacerein",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "diacerein",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Scales the distribution volume linearly, because the source model",
        "predicts Vss per kilogram (Simcyp 'Vss (L/kg)' = 0.0930, S6 File",
        "'Input Sheet'). The Simcyp virtual population's weights are not",
        "reported, so 70 kg is used as the reference weight. Clearance is",
        "NOT weight-scaled: the source enters it as an absolute oral",
        "clearance in L/h."
      ),
      source_name = "body weight (Simcyp Sim-Healthy Volunteers population)"
    ),
    FORM_DCN_SD = list(
      description = paste(
        "1 = dose given as the optimised diacerein solid dispersion (PEG 8000,",
        "drug:polymer 1:4 w/w, prepared by solvent evaporation); 0 = plain",
        "crystalline diacerein powder."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (plain crystalline diacerein)",
      notes = paste(
        "Selects the in-vitro dissolution profile that drives drug release.",
        "Per-dose record: the profile is read when the dose is released, so a",
        "subject can carry different values on different dose records."
      ),
      source_name = "optimized SD vs plain DCN (Fouad 2021 Table 6)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 10L,
    n_studies = 1L,
    age_range = "20-60 years (Simcyp healthy-volunteer trial design; the article text says 40-55 years)",
    weight_median = "70 kg (reference weight for the L/kg volume; the Simcyp population weights are not reported)",
    sex_female_pct = 50,
    disease_state = "Healthy adult volunteers (virtual Simcyp population).",
    dose_range = "Single 50 mg oral dose with water, fasted.",
    regions = "Simcyp Sim-Healthy Volunteers population file (version 17).",
    notes = paste(
      "This is a PBPK simulation analysis, not a population-PK fit. The",
      "article describes ten trials of ten subjects (n = 100), but the",
      "deposited healthy-volunteer workbooks (S6 and S8 Files, which have",
      "the same inputs and results) are one trial of ten subjects;",
      "n_subjects follows the",
      "deposit. The source model was verified against the single-dose",
      "clinical study of Nguyen et al. (the article's reference 1). A",
      "geriatric (65-75 years) simulation of the solid dispersion is also",
      "reported but is outside the scope of this model."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Every parameter is fixed: nothing was estimated in building this
    # reduction. Values are the Simcyp inputs and Simcyp-computed
    # quantities in the deposited workbook (Fouad 2021 S6 File, sheet
    # 'Input Sheet'), which reproduces the article's Table 6 plain-drug
    # row exactly.
    # ------------------------------------------------------------------

    lka <- fixed(log(2.098))
    label("First-order absorption rate constant of dissolved drug ka (1/h)")
    # S6 File 'Input Sheet', Absorption block, 'ka (1/h)' = 2.0979
    # (input type 'Predicted' from 'Peff,man (10-4 cm/s)' = 4.80).

    lfdepot <- fixed(log(0.992))
    label("Fraction of dissolved drug absorbed (fraction)")
    # S6 File 'Input Sheet', Absorption block, 'fa' = 0.9922
    # (input type 'Predicted'). Gut-wall and hepatic availability are both
    # 1 (S6 File 'Summary', 'Fg (Sub)' = 1 and 'Fh (Sub)' = 1).

    # Distribution volume. The minimal-PBPK distribution model is used with
    # no single adjusting compartment (S6 File 'Input Sheet': 'Volume
    # [Vsac] (L/kg)' = 1e-5, SAC kin = kout = 0), so the predicted Vss is
    # the one-compartment plasma volume:
    #   0.0930 L/kg * 70 kg = 6.512 L
    # Table 3 of the article lists 'Vss (L/Kg) 0.23' from its reference 39,
    # but the deposited run used the Simcyp-predicted value
    # ('Vss input type: Predicted', 'Prediction Method: Method 2 (Rodgers
    # et al)', 'Vss (L/kg)' = 0.09303), and only the predicted value
    # reproduces the Table 6 exposures (see the vignette).
    lvc <- fixed(log(6.512))
    label("Central volume vc at the 70 kg reference weight (L)")
    # S6 File 'Input Sheet', Distribution block, 'Vss (L/kg)' = 0.09303,
    # times 70 kg.

    lcl <- fixed(log(1.5))
    label("Clearance (L/h)")
    # Fouad 2021 Table 3, 'CL (L/h) 1.5'; S6 File 'Input Sheet',
    # Elimination block, 'Clearance Type: In Vivo Clearance',
    # 'CL (po) (L/h)' = 1.5.

    etalcl ~ 0.0862
    # S6 File 'Input Sheet', 'CV CL (po) (%)' = 30; log-normal variance
    # log(1 + 0.30^2) = 0.0862.

    # Fouad 2021 is a PBPK simulation analysis, not a population-PK fit,
    # and reports no residual-error model. The residual error is fixed at
    # zero rather than invented.
    propSd <- fixed(0)
    label("Proportional residual error SD (fraction; zero, no error model reported by the source)")
  })

  model({
    # 1. Individual parameters. Body weight scales the volume only.
    ka <- exp(lka)
    fdepot <- exp(lfdepot)
    vc <- exp(lvc) * WT / 70
    cl <- exp(lcl + etalcl)
    kel <- cl / vc

    # 2. Dissolution (release) rate, as a fraction of the dose per hour.
    # The source enters each formulation's in-vitro dissolution profile
    # (phosphate buffer pH 6.8) as a discrete profile, linearly
    # interpolated and held at its last value after 1 h, so the release
    # rate is the slope of each 15-minute segment and zero afterwards.
    #   Plain diacerein, S6 File 'Input Sheet' 'Dissolution Profile':
    #     0, 13.27, 32.90, 42.30, 48.30 % at 0, 0.25, 0.5, 0.75, 1 h
    #   Optimised solid dispersion, S2 File sheet 'Dissolution Data' row
    #   'OPTIMIZED' (the workbook for this run is not deposited; the same
    #   15-minute grid as the plain-drug run is used):
    #     0, 95.197, 96.148, 98.349, 99.127 % at 0, 0.25, 0.5, 0.75, 1 h
    # The clock restarts at each dose, so repeat doses are handled exactly
    # when they are at least 1 h apart. Before the first dose tad() is
    # missing; nothing is being released then.
    tdis <- tad(depot)
    if (is.na(tdis)) tdis <- 1e6
    rdis_plain <- (13.27 / 0.25) * (tdis < 0.25) +
      ((32.90 - 13.27) / 0.25) * (tdis >= 0.25) * (tdis < 0.5) +
      ((42.30 - 32.90) / 0.25) * (tdis >= 0.5) * (tdis < 0.75) +
      ((48.30 - 42.30) / 0.25) * (tdis >= 0.75) * (tdis < 1)
    rdis_sd <- (95.197 / 0.25) * (tdis < 0.25) +
      ((96.148 - 95.197) / 0.25) * (tdis >= 0.25) * (tdis < 0.5) +
      ((98.349 - 96.148) / 0.25) * (tdis >= 0.5) * (tdis < 0.75) +
      ((99.127 - 98.349) / 0.25) * (tdis >= 0.75) * (tdis < 1)
    rdis <- ((1 - FORM_DCN_SD) * rdis_plain + FORM_DCN_SD * rdis_sd) / 100

    # Released amount per hour (mg/h) from the most recent dose.
    doseAmt <- podo(depot)
    if (is.na(doseAmt)) doseAmt <- 0
    rel <- doseAmt * rdis

    # 3. ODEs. Amounts in mg. The undissolved remainder (51.7% of a plain
    # diacerein dose) stays in the depot and is never absorbed, as in the
    # source, where the dissolution profile is held at its last value.
    d/dt(depot) <- -rel
    d/dt(gut_lumen) <- rel - ka * gut_lumen
    d/dt(central) <- fdepot * ka * gut_lumen - kel * central

    # 4. Observation: mg / L (= ug/mL, the units of Table 6).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
