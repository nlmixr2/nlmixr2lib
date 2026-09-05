Munn_2024_carprofen_raw <- function() {
  description <- paste(
    "Preclinical (sheep).",
    "Three-compartment population PK model for carprofen after a single",
    "4 mg/kg intravenous bolus in merino wethers, fitted to plasma and",
    "subcutaneous tissue-cage fluid simultaneously. Plasma disposition is",
    "a two-compartment model (central + peripheral1) in Monolix",
    "micro-constant form (Vc, CL, k12, k21). Each animal also carried a",
    "6 cm and a 10 cm subcutaneous tissue cage; each cage is a separate",
    "state whose concentration is driven by the central-compartment",
    "concentration through first-order k13 (into the cage) and k31 (out of",
    "the cage), with NO back-transfer into the central compartment -- the",
    "negligible-mass approach of Sheiner that the source inherits from",
    "Munn 2022. Cage length enters as a categorical effect on k13 and k31",
    "with the 6 cm cage as the reference level. Because every volume and",
    "clearance term in the source is normalised to body weight, the model",
    "is coded on a per-kilogram basis: the dosed amount is mg/kg and the",
    "volume is L/kg, so central/vc lands directly in ug/mL, the assay",
    "units. This is the RAW carprofen fit; see",
    "Munn_2024_carprofen_moxidectinCorrected for the companion fit to the",
    "moxidectin-corrected concentrations.",
    sep = " ")
  reference <- paste(
    "Munn R, Whittem T.",
    "Moxidectin is a candidate for use as an in vivo internal standard in",
    "pharmacokinetic studies, as demonstrated with use in simultaneous",
    "tissue cage and ultrafiltration fluid collection.",
    "Front Vet Sci. 2024;11:1332974.",
    "doi:10.3389/fvets.2024.1332974.",
    "Parameter values from Supplementary Table S1, 'Raw Analysis' panel.",
    "The tissue-cage model structure is inherited from",
    "Munn R, Whittem T, Woodward AP.",
    "The surface area to volume ratio changes the pharmacokinetic and",
    "pharmacodynamic parameters in the subcutaneous tissue cage model:",
    "as illustrated by carprofen in sheep.",
    "Front Vet Sci. 2022;9:905797. doi:10.3389/fvets.2022.905797.",
    sep = " ")
  vignette <- "Munn_2024_carprofen_moxidectin_internal_standard"

  # `cage6` / `cage10` are paper-mechanistic states: they hold the drug
  # CONCENTRATION inside a surgically implanted subcutaneous tissue cage of
  # the stated length. They are deliberately NOT mapped to the canonical
  # `isf` compartment -- the source's own Conclusion states that "tissue
  # cage-derived samples are less likely to represent true physiological
  # spaces than samples obtained by ultrafiltration", so calling the cage
  # interstitial fluid would misstate the paper.
  paper_specific_compartments <- c("cage6", "cage10")

  units <- list(
    time          = "h",
    dosing        = "mg/kg",
    concentration = "ug/mL"
  )

  # Per-kg normalisation: the source reports Vc in L/kg and CL in L/h/kg
  # and the dose as 4 mg/kg, so state amounts are mg/kg and central/vc is
  # mg/L == ug/mL, the units of the carprofen assay (calibration range
  # 0.25-100 ug/mL, Methods 2.1).
  compartmentData <- list(
    central     = list(analyte = "carprofen", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "carprofen", units = "mg/kg", specimen = "plasma", verified = TRUE),
    # The cage states hold a concentration, not an amount: they are driven
    # by Cc and never exchange mass with central (see model() note).
    cage6       = list(analyte = "carprofen", units = "ug/mL", specimen = "tissue", verified = TRUE),
    cage10      = list(analyte = "carprofen", units = "ug/mL", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species        = "sheep (merino wether)",
    n_subjects     = 8L,
    n_studies      = 1L,
    age_range      = "approximately 18 months",
    weight_range   = "42-51.5 kg",
    sex_female_pct = 0,
    disease_state  = "Healthy (veterinary clinical examination plus routine haematology and biochemistry before enrolment)",
    dose_range     = paste(
      "Single 4 mg/kg carprofen intravenous bolus into a cephalic vein at",
      "time zero. Separately, 0.2 mg/kg moxidectin (Cydectin) was injected",
      "subcutaneously into a hind limb 14 days before the experiment to",
      "serve as the in vivo internal standard; moxidectin is not a state",
      "in this model."
    ),
    regions        = "Australia (University of Melbourne, Werribee, Victoria)",
    notes          = paste(
      "Two tissue cages (6 cm and 10 cm length) were implanted",
      "subcutaneously in the neck three weeks before the experiment, so",
      "every animal contributes both cage sizes simultaneously. Plasma and",
      "tissue-cage samples were taken at -0.5, 0.5, 1, 2, 3, 4, 5, 7, 24,",
      "36, 48 and 72 h. Ultrafiltration probes were also inserted but",
      "yielded too little sample volume to analyse, so NO ultrafiltration",
      "compartment was fitted and none appears here (Results 3.1).",
      "Experiments were run October-November 2020; ethics approval",
      "UoM 2015111. See Munn 2024 Materials and methods."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Munn 2024 Supplementary Table S1,
    # "Raw Analysis" panel, Fixed Effects. Monolix 2023R1 SAEM.
    #
    # HALF-LIFE CROSS-CHECK. The Discussion states "an estimated half-life
    # of 19.5 h compared to 27.2 h in the previous dataset". Both numbers
    # are reproduced exactly as ln(2)/(CL/Vc):
    #   this study: 0.693 / (0.0016 / 0.045)   = 19.5 h
    #   Munn 2022 : 0.693 / (0.00235 / 0.0924) = 27.2 h
    # so the printed half-life independently confirms BOTH Vc and CL as
    # transcribed here. (Note this is the central-compartment half-life,
    # not the true terminal half-life of the two-compartment system, which
    # is 27.3 h for this parameter set.)
    # ------------------------------------------------------------------
    lvc  <- log(0.045)   ; label("Central volume of distribution Vc (log L/kg)")            # Suppl Table S1 raw: Central Volume = 0.045 L/kg (SE 0.0021, RSE 4.69%)
    lcl  <- log(0.0016)  ; label("Clearance from central CL (log L/h/kg)")                  # Suppl Table S1 raw: Clearance = 0.0016 L/h.kg (SE 0.000088, RSE 5.54%)
    lk12 <- log(0.11)    ; label("Rate constant central -> peripheral1 k12 (log 1/h)")      # Suppl Table S1 raw: k12 = 0.11 /h (SE 0.044, RSE 40.7%)
    lk21 <- log(0.30)    ; label("Rate constant peripheral1 -> central k21 (log 1/h)")      # Suppl Table S1 raw: k21 = 0.3 /h (SE 0.096, RSE 32.2%)

    # ------------------------------------------------------------------
    # Tissue-cage rate constants. The 6 cm cage is the reference level:
    # Supplementary Table S1 prints "Covariate for k13 for Cage Size 6 cm"
    # and "... k31 for Cage Size 6 cm" as 0* (Monolix's fixed-reference
    # marker), with a single estimated offset for the 10 cm cage.
    # ------------------------------------------------------------------
    lk13 <- log(0.052)   ; label("Rate constant central -> tissue cage k13, 6 cm cage (log 1/h)")   # Suppl Table S1 raw: k13 = 0.052 /h (SE 0.013, RSE 25.3%)
    lk31 <- log(0.22)    ; label("Rate constant out of tissue cage k31, 6 cm cage (log 1/h)")       # Suppl Table S1 raw: k31 = 0.22 /h (SE 0.078, RSE 35.8%)

    # ------------------------------------------------------------------
    # Cage-length effects on the cage rate constants, 10 cm vs the 6 cm
    # reference.
    #
    # SCALE OF THESE COEFFICIENTS. They are applied on the LOG scale,
    # k_10cm = k_6cm * exp(beta), which is Monolix's default for a
    # log-normally distributed parameter. The alternative additive reading
    # (k_10cm = k_6cm + beta) is falsified twice over:
    #
    #  1. Internally. k13 would become 0.052 + (-0.056) = -0.004 /h, a
    #     NEGATIVE rate constant, which is impossible.
    #  2. Against the upstream paper. Munn 2022 Table 3 reports the same
    #     covariate as a continuous per-cm effect (k31 -0.147, k13 -0.0378
    #     on k31 = 0.455, k13 = 0.124) and its Table 5 tabulates the
    #     resulting median rate constants for cage sizes 0-18 cm. The log
    #     reading reproduces that table across the whole range (k13 at
    #     10 cm: 0.124*exp(-0.0378*10) = 0.085 vs 0.0862 tabulated; k31 at
    #     18 cm: 0.455*exp(-0.147*18) = 0.032 vs 0.0400 tabulated), while
    #     the additive reading goes negative from 4 cm onward.
    #
    # The "(h-1)" unit shown on these rows of Supplementary Table S1 is a
    # Monolix table-formatting artefact -- the row inherits the unit of the
    # parent parameter. The coefficients themselves are dimensionless.
    # ------------------------------------------------------------------
    e_cage10_k13 <- -0.056 ; label("Log-scale effect of a 10 cm (vs 6 cm) tissue cage on k13 (unitless)")   # Suppl Table S1 raw: Covariate for k13 for Cage Size 10 cm = -0.056 (SE 0.03, RSE 54.0%)
    e_cage10_k31 <- -0.17  ; label("Log-scale effect of a 10 cm (vs 6 cm) tissue cage on k31 (unitless)")   # Suppl Table S1 raw: Covariate for k31 for Cage Size 10 cm = -0.17 (SE 0.045, RSE 26.5%)

    # ------------------------------------------------------------------
    # IIV -- Supplementary Table S1, "Standard Deviation of the Random
    # Effects" block of the Raw Analysis panel. The tabulated values are
    # the SD of eta on the LOG scale; the table's own C.V.(%) column
    # confirms it, since CV = sqrt(exp(omega^2) - 1):
    #   k12: sqrt(exp(0.42^2) - 1) = 43.9%  vs 43.73% tabulated
    #   CL : sqrt(exp(0.13^2) - 1) = 13.1%  vs 13.48% tabulated
    # nlmixr2 takes variances here, so each entry below is SD^2.
    #
    # Vc and k21 were driven to essentially zero by the estimator and are
    # NOT identifiable: their relative standard errors are 1.82e+7% and
    # 2.72e+7% respectively. They are reproduced as reported rather than
    # dropped, so the model states what the paper states, but a user
    # should treat Vc and k21 as having no meaningful between-animal
    # variability in this dataset.
    # ------------------------------------------------------------------
    etalvc  ~ 8.281e-09  # Suppl Table S1 raw: Volume SD = 0.000091 (CV 0.0091%, RSE 1.82e+7%) -> 0.000091^2
    etalk12 ~ 0.1764     # Suppl Table S1 raw: k12 SD = 0.42 (CV 43.73%, RSE 42.8%) -> 0.42^2
    etalk21 ~ 2.704e-07  # Suppl Table S1 raw: k21 SD = 0.00052 (CV 0.052%, RSE 2.72e+7%) -> 0.00052^2
    etalcl  ~ 0.0169     # Suppl Table S1 raw: Clearance SD = 0.13 (CV 13.48%, RSE 31.7%) -> 0.13^2

    # ------------------------------------------------------------------
    # Residual error -- NOT REPORTED. Supplementary Table S1 lists only
    # the fixed effects and the random-effect SDs; it has no error-model
    # rows, and the main text gives no residual-error parameter. The
    # upstream Munn 2022 Table 3 DOES report a proportional error model
    # (b1 = 0.136 for plasma, b2 = 0.468 for tissue cage), but those are
    # estimates from a different dataset and are deliberately NOT imported
    # here. Each residual SD is therefore fixed to zero, which makes the
    # model a typical-value / IIV-only simulator. See the vignette
    # "Assumptions and deviations" section.
    # ------------------------------------------------------------------
    propSd         <- fixed(0) ; label("Proportional residual error, plasma (fraction; not reported by the source)")        # Not reported; see note above
    propSd_Ccage6  <- fixed(0) ; label("Proportional residual error, 6 cm tissue cage (fraction; not reported)")            # Not reported; see note above
    propSd_Ccage10 <- fixed(0) ; label("Proportional residual error, 10 cm tissue cage (fraction; not reported)")           # Not reported; see note above
  })

  model({
    # 1. Individual plasma-disposition parameters (log-normal IIV).
    vc  <- exp(lvc  + etalvc)
    cl  <- exp(lcl  + etalcl)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)

    # 2. Cage rate constants, one pair per implanted cage length. The 6 cm
    # cage is the reference; the 10 cm cage carries the estimated
    # log-scale offset. No IIV was estimated on k13 or k31 in this study.
    k13_cage6  <- exp(lk13)
    k31_cage6  <- exp(lk31)
    k13_cage10 <- exp(lk13 + e_cage10_k13)
    k31_cage10 <- exp(lk31 + e_cage10_k31)

    kel <- cl / vc

    # 3. Plasma disposition -- a self-contained two-compartment model.
    # Munn 2022 Results: the plasma model was built first, WITHOUT the
    # tissue-cage data, and the cages were then appended without altering
    # it.
    d/dt(central)     <- k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <- k12 * central     - k21 * peripheral1

    # 4. Plasma concentration. Amount in mg/kg over volume in L/kg gives
    # mg/L == ug/mL.
    Cc <- central / vc

    # 5. Tissue cages. Each cage state holds a CONCENTRATION and is driven
    # by the central-compartment concentration, gaining drug at k13 and
    # losing it at k31. There is deliberately no term returning drug to
    # `central`: Munn 2022 states the cage kinetics "are driven by the
    # central compartment concentrations without altering the central
    # compartment concentrations", because only a negligible proportion of
    # the dose enters the cages. Both cages are solved simultaneously
    # because each sheep carried both.
    d/dt(cage6)  <- k13_cage6  * Cc - k31_cage6  * cage6
    d/dt(cage10) <- k13_cage10 * Cc - k31_cage10 * cage10

    Ccage6  <- cage6
    Ccage10 <- cage10

    Cc      ~ prop(propSd)
    Ccage6  ~ prop(propSd_Ccage6)
    Ccage10 ~ prop(propSd_Ccage10)
  })
}
