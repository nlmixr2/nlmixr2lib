Ibrahim_2019_oral2cmt_zeroorder_mm <- function() {
  description <- paste(
    "Methodology reference. Ground-truth (data-generating) two-compartment",
    "population PK model with ZERO-ORDER absorption and Michaelis-Menten",
    "elimination for a HYPOTHETICAL drug, taken from the 'Simple PK example'",
    "in the supplementary material of Ibrahim 2019, the paper that introduces",
    "model-based conditional weighted residuals (CWRES) analysis for structural",
    "model assessment. There is no real molecule and no real patients: the",
    "authors DEFINED these six structural parameters, simulated a data set of",
    "100 subjects from them, and then fitted that data set with both the true",
    "model and a deliberately misspecified variant (first-order instead of",
    "zero-order absorption) so that their CWRES-bias diagnostic could be shown",
    "to recover a known prediction bias. Every value here is therefore an",
    "author-chosen simulation constant, not an estimate, and all are encoded",
    "with fixed(). The source reports no units, no dose, no inter-individual",
    "variability magnitudes and no residual-error magnitude; the model is",
    "consequently typical-value-only and its units are placeholders. See the",
    "vignette Errata for the full list of gaps and the one interpretive call",
    "(ka0 read as a zero-order input RATE, per its printed row label).",
    sep = " "
  )
  reference <- paste(
    "Ibrahim MMA, Ueckert S, Freiberga S, Kjellsson MC, Karlsson MO.",
    "Model-Based Conditional Weighted Residuals Analysis for Structural Model",
    "Assessment. AAPS J. 2019 Feb 27;21(2):34.",
    "doi:10.1208/s12248-019-0305-2. PMCID PMC6394649.",
    "Structure and all six parameter values transcribed from Supplementary",
    "Material 1 (12248_2019_305_MOESM1_ESM.docx), section 'Simple PK example',",
    "Table 1 'Simulation specifications and dOFVBias'.",
    sep = " "
  )
  vignette <- "Ibrahim_2019_oral2cmt_zeroorder_mm"
  units <- list(
    # The source reports no units for any of the six constants, for the dose,
    # or for time. Placeholder tokens are used so that downstream consumers
    # cannot mistake an assumed unit for a published one.
    time = "time_unit",
    dosing = "dose_unit",
    concentration = "dose_unit/volume_unit"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. There is no molecule and no specimen in this source:
  # the states are the abstract central and peripheral compartments of the
  # authors' simulation model. verified = FALSE means NOT checked against a
  # source statement, because the source makes none.
  compartmentData <- list(
    central = list(
      analyte = "hypothetical drug",
      units = NA_character_,
      specimen = "not applicable",
      verified = FALSE
    ),
    peripheral1 = list(
      analyte = "hypothetical drug",
      units = NA_character_,
      specimen = "not applicable",
      verified = FALSE
    )
  )

  covariateData <- list()

  population <- list(
    species = "None (methodology paper; simulation-only toy model with no drug, no patients and no fitted estimates).",
    n_subjects = 100L,
    n_studies = 1L,
    disease_state = "N/A (Monte Carlo simulation study demonstrating a model-diagnostic method; not a fit of any real molecule).",
    dose_range = "Not reported. The supplement states only that a data set of 100 subjects was simulated; it gives no dose amount, no dosing route beyond 'zero order absorption', and no sampling schedule.",
    regions = "N/A",
    scope_note = paste(
      "Filed under inst/modeldb/pharmacokinetics/ (not specificDrugs/) because",
      "there is no drug: the model is an author-invented hypothetical used to",
      "demonstrate a diagnostic method. This mirrors the packaging of",
      "Beal_2001_iv1cmt_bql.R and the Schoning_2026_oral1cmt_* family.",
      "The paper's own product is a CWRES-bias regression (a residual",
      "post-processing diagnostic, now shipped in PsN's qa tool from version",
      "4.8.1), which is not a pharmacokinetic model and is not encoded here.",
      "The two integrated glucose-insulin models the main text uses as",
      "demonstration vehicles (the IGI model of Silber 2007,",
      "doi:10.1177/0091270007304457, and the integrated minimal model of",
      "Largajolli 2013, PAGE 22 Abstract 2762) are cited backbones whose",
      "parameter values appear nowhere in this paper or its supplement; they",
      "are separate primary sources, not layers of this extraction.",
      sep = " "
    ),
    notes = paste(
      "Supplementary Material 1, 'Simple PK example': 'In this example we used",
      "a two-compartment PK model with zero order absorption and",
      "Michaelis-Menten elimination to simulate a dataset of 100 subjects. The",
      "simulated data set was then used to fit two models: a true model (same",
      "as the simulation model) and a misspecified model that is the same as",
      "the simulation model except for using 1st order absorption process",
      "instead of the true zero order absorption process.' Table 1 of that",
      "supplement reports dOFVBias = -165.8 for the misspecified fit and lists",
      "the six 'Simulated parameters' encoded below, with the footnote 'Vc",
      "volume of central compartment, Vp volume of peripheral compartment, ka0",
      "zero order absorption rate, Q intercompartmental clearance, KM and VMAX",
      "Michaelis-Menten elimination parameters.' The misspecified companion fit",
      "is NOT encoded as a sibling model file because the supplement publishes",
      "no parameter values for it; see the vignette Errata.",
      sep = " "
    )
  )

  ini({
    # All six values are author-chosen simulation constants from Supplementary
    # Material 1, Table 1, column 'Simulated parameters'. None is an estimate,
    # none carries an uncertainty, and none is reported with a unit.
    lvc <- fixed(log(4.14))
    label("Central volume of distribution Vc (volume_unit)")
    lvp <- fixed(log(7))
    label("Peripheral volume of distribution Vp (volume_unit)")
    lq <- fixed(log(3.24))
    label("Intercompartmental clearance Q (volume_unit/time_unit)")
    lvmax <- fixed(log(9.21))
    label("Michaelis-Menten maximum elimination rate VMAX (dose_unit/time_unit)")
    lkm <- fixed(log(16.34))
    label("Michaelis-Menten constant KM (dose_unit/volume_unit)")
    lr1 <- fixed(log(10.28))
    label("Zero-order input rate ka0 into the central compartment (dose_unit/time_unit)")

    # No inter-individual variability is encoded: the supplement states that
    # 100 subjects were simulated but reports no omega for any parameter, and
    # the deposited files in Supplementary Material 2 belong to a different,
    # 26-subject worked example of the R implementation (its fitted profile is
    # reproduced by a FIRST-order absorption structure, not by the constants
    # above). Inventing variances is not permitted, so the model is
    # typical-value-only.
    addSd <- fixed(0)
    label("Additive residual SD; magnitude not reported by the source")
  })

  model({
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)
    vmax <- exp(lvmax)
    km <- exp(lkm)
    r1 <- exp(lr1)

    Cc <- central / vc
    Cp <- peripheral1 / vp

    # Zero-order absorption. ka0 is read as an input RATE (amount per unit
    # time) rather than a duration, following its row label in the supplement
    # ('ka0 zero order absorption rate') and the library's canonical lr1 / r1
    # encoding for a zero-order input. rxode2 reads rate() at solve time for
    # dose records flagged rate = -1 and turns the bolus into a constant-rate
    # input of duration amt / r1; the dose goes straight into central because
    # the authors describe a two-compartment model with no depot state.
    rate(central) <- r1

    # Michaelis-Menten elimination is written on the central CONCENTRATION, so
    # VMAX carries amount-per-time and KM carries concentration units. There is
    # no linear clearance term: the supplement's parameter list contains no CL,
    # and no cl / vc pair is created here, so rxode2 cannot silently replace
    # these ODEs with an analytic linear solution.
    d/dt(central) <- -q * (Cc - Cp) - vmax * Cc / (km + Cc)
    d/dt(peripheral1) <- q * (Cc - Cp)

    Cc ~ add(addSd)
  })
}
