Alidori_2017_fCNT_kidney_cyno <- function() {
  description <- paste(
    "Preclinical (cynomolgus monkey). Renal tissue kinetics of a tracer-radiolabelled",
    "ammonium-functionalized single-walled carbon nanotube ([86Y]fCNT) measured by dynamic",
    "whole-body PET/CT and fitted with a serial (catenary) two-compartment model driven by a",
    "whole-blood input forcing function. Compartment 1 is the renal tissue space in toto and",
    "compartment 2 is the proximal-tubule-cell population that acts as a hold-up volume delaying",
    "washout; drug filtered from blood is either resorbed into the tubule-cell sink or eliminated",
    "intact into the urine. Four first-order rate constants were estimated (SAAM II) with the",
    "renal blood-volume fraction fixed at 15%. Deterministic typical-value model (two animals, no",
    "between-subject random effects or residual-error model reported).",
    sep = " "
  )
  reference <- paste(
    "Alidori S, Asqiriba K, Londero P, Bergkvist M, Betancur M, Cheung PK, et al. (2017).",
    "Carbon nanotubes exhibit fibrillar pharmacology in primates. PLoS One 12(8):e0183902.",
    "doi:10.1371/journal.pone.0183902.",
    sep = " "
  )
  vignette <- "Alidori_2017_fCNT_primate"

  # C1 is the canonical renal `kidney` compartment; C2 (the proximal-tubule-cell hold-up
  # volume) is a genuinely paper-mechanistic state with no canonical role.
  paper_specific_compartments <- c("kidney_ptc")

  units <- list(
    time = "min",
    dosing = "n/a (no exogenous dosing event; tissue uptake is driven by the blood forcing function SUV_BLOOD)",
    concentration = "SUV (standardized uptake value, dimensionless)"
  )

  compartmentData <- list(
    kidney = list(
      analyte = "[86Y]fCNT",
      units = "SUV",
      specimen = "tissue",
      verified = TRUE
    ),
    kidney_ptc = list(
      analyte = "[86Y]fCNT",
      units = "SUV",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    SUV_BLOOD = list(
      description = "Whole-blood [86Y]fCNT activity time course (forcing function)",
      units = "SUV",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying error-free regressor; declared as linear(SUV_BLOOD) so rxode2 interpolates",
        "linearly between rows. Alidori 2017, Compartmental modeling analysis: 'The input",
        "function was derived from the dynamic PET images by calculating the mean radioactivity",
        "concentration within a VOI placed over the heart'; the whole-blood activity is assumed",
        "entirely available for tissue uptake. The per-animal heart-VOI curve is not tabulated in",
        "the paper, so users must supply the observed or an assumed blood SUV profile through SUV_BLOOD.",
        sep = " "
      ),
      source_name = "CB"
    )
  )

  population <- list(
    species = "cynomolgus monkey (Macaca fascicularis)",
    n_subjects = 2,
    n_studies = 1,
    age_range = "3.5 years",
    weight_range = "3.2-4.5 kg at dosing (Table 3)",
    sex_female_pct = 0,
    disease_state = "healthy",
    dose_range = "single IV injection of 5 mg (5 MBq, 0.135 mCi) [86Y]fCNT; approx 1 mg/kg",
    regions = "not reported",
    notes = paste(
      "Two healthy male cynomolgus monkeys (Charles River), 3.5 years old, imaged on a Siemens",
      "Biograph64 mCT PET/CT for a 3 hour dynamic protocol immediately after IV [86Y]fCNT. The",
      "correlative mouse tissue-harvest study in the same paper was not fitted with a",
      "compartmental model and is not packaged here.",
      sep = " "
    )
  )

  ini({
    # ===================== Renal serial two-compartment rate constants =====================
    # Alidori 2017 Fig 3B and Supporting Information S1 Table (SAAM II fit, two cynomolgus
    # monkeys). All four constants were estimated with the renal blood-volume fraction fixed at
    # 15%; they are held fixed here because this packaged model reproduces the published
    # typical-value fit rather than re-estimating it. Natural scale (1/min) rather than log scale
    # so each value source-traces directly to the printed table.
    # S1 Table reports st.dev / CV% / 95% CI for each (parameter uncertainty of the
    # deterministic fit, NOT between-subject variability): k1 7.0% CV, k2 14.5%, k3 13.2%,
    # k4 17.9%.
    kup <- fixed(0.93922)
    label("Blood -> kidney uptake rate constant k1 (1/min)") # S1 Table 'k1' = 0.93922 (Fig 3B 0.939)
    kelim <- fixed(0.15439)
    label("Kidney -> urine elimination rate constant k2 (1/min)") # S1 Table 'k2' = 0.15439 (Fig 3B 0.154)
    kseq <- fixed(0.02875)
    label("Kidney -> proximal-tubule-cell sink rate constant k3 (1/min)") # S1 Table 'k3' = 0.02875 (Fig 3B 0.0289)
    kret <- fixed(0.00382)
    label("Tubule-cell sink -> kidney return rate constant k4 (1/min)") # S1 Table 'k4' = 0.00382 (Fig 3B 0.00382)

    # Renal blood-volume fraction, fixed at 15% in the fit. Alidori 2017 Compartmental modeling
    # analysis and S1 Table caption: 'Assume a constant 15% blood volume.'
    fvb <- fixed(0.15)
    label("Renal blood-volume fraction (dimensionless)") # S1 Table caption: constant 15% renal blood volume

    # Deterministic model: the paper reports no residual-error model for the compartmental fit
    # (weighted least squares over two animals). A small fixed proportional SD is supplied only
    # so the model carries a valid residual endpoint for simulation; it is not a fitted value.
    propSd <- fixed(0.1)
    label("Proportional residual SD on renal SUV (fraction)") # not reported; nominal (see vignette Errata)
  })

  model({
    # Whole-blood [86Y]fCNT activity is an exogenous, error-free forcing function derived from
    # the heart VOI; the model does not describe blood kinetics. Supplied in SUV at every row and
    # linearly interpolated between rows.
    linear(SUV_BLOOD)

    # Serial (catenary) topology, Alidori 2017 Fig 3A:
    #   blood (SUV_BLOOD) --k1--> kidney (C1)
    #   kidney (C1) --k2--> urine                (elimination)
    #   kidney (C1) --k3--> kidney_ptc (C2)      (resorption into tubule-cell sink)
    #   kidney_ptc (C2) --k4--> kidney (C1)      (slow return)
    # Both states hold SUV of the non-blood tissue signal; the hold-up sink C2 delays washout.
    d/dt(kidney) <- kup * SUV_BLOOD - (kelim + kseq) * kidney + kret * kidney_ptc
    d/dt(kidney_ptc) <- kseq * kidney - kret * kidney_ptc
    kidney(0) <- 0
    kidney_ptc(0) <- 0

    # Observed renal-cortex VOI SUV = fixed blood-volume contribution + the two tissue
    # compartments. The 15% blood fraction was held fixed during the fit. The observed activity
    # (SUV) is the model's single concentration-like output, so it carries the canonical name Cc
    # (units$concentration = "SUV"); there is no plasma central compartment in this tissue model.
    Cc <- fvb * SUV_BLOOD + kidney + kidney_ptc
    Cc ~ prop(propSd)
  })
}
