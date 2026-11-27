Alidori_2017_fCNT_liver_cyno <- function() {
  description <- paste(
    "Preclinical (cynomolgus monkey). Hepatic tissue kinetics of a tracer-radiolabelled",
    "ammonium-functionalized single-walled carbon nanotube ([86Y]fCNT) measured by dynamic",
    "whole-body PET/CT and fitted with a parallel (mammillary) two-compartment model driven by a",
    "whole-blood input forcing function. Blood-borne tracer partitions in parallel into the",
    "hepatic tissue space in toto (compartment 1) and the liver sinusoidal endothelium",
    "(compartment 2, a slow hold-up volume that scavenges fCNT via Stabilin receptors); each",
    "compartment has its own loss route back to blood or into bile. Four first-order rate",
    "constants were estimated (SAAM II) with the hepatic blood-volume fraction fixed at 20%.",
    "Deterministic typical-value model (two animals, no between-subject random effects or",
    "residual-error model reported).",
    sep = " "
  )
  reference <- paste(
    "Alidori S, Asqiriba K, Londero P, Bergkvist M, Betancur M, Cheung PK, et al. (2017).",
    "Carbon nanotubes exhibit fibrillar pharmacology in primates. PLoS One 12(8):e0183902.",
    "doi:10.1371/journal.pone.0183902.",
    sep = " "
  )
  vignette <- "Alidori_2017_fCNT_primate"

  # C1 is the canonical `liver` compartment; C2 (the liver sinusoidal endothelial hold-up
  # volume) is a genuinely paper-mechanistic state with no canonical role.
  paper_specific_compartments <- c("liver_lsec")

  units <- list(
    time = "min",
    dosing = "n/a (no exogenous dosing event; tissue uptake is driven by the blood forcing function SUV_BLOOD)",
    concentration = "SUV (standardized uptake value, dimensionless)"
  )

  compartmentData <- list(
    liver = list(
      analyte = "[86Y]fCNT",
      units = "SUV",
      specimen = "tissue",
      verified = TRUE
    ),
    liver_lsec = list(
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
    # ===================== Hepatic parallel two-compartment rate constants =====================
    # Alidori 2017 Fig 4B and Supporting Information S2 Table (SAAM II fit, two cynomolgus
    # monkeys). All four constants were estimated with the hepatic blood-volume fraction fixed at
    # 20%; they are held fixed here because this packaged model reproduces the published
    # typical-value fit rather than re-estimating it. Natural scale (1/min) rather than log scale
    # so each value source-traces directly to the printed table.
    # S2 Table reports st.dev / CV% / 95% CI for each (parameter uncertainty of the
    # deterministic fit, NOT between-subject variability): k1 14.2% CV, k2 14.6%, k3 3.8%,
    # k4 11.2%.
    kup <- fixed(0.21158)
    label("Blood -> liver C1 uptake rate constant k1 (1/min)") # S2 Table 'k1' = 0.21158 (Fig 4B 0.211)
    kelim <- fixed(0.41722)
    label("Liver C1 -> blood/bile loss rate constant k2 (1/min)") # S2 Table 'k2' = 0.41722 (Fig 4B 0.417)
    kseq <- fixed(0.01003)
    label("Blood -> sinusoidal endothelium uptake rate constant k3 (1/min)") # S2 Table 'k3' = 0.01003 (Fig 4B 0.0100)
    kret <- fixed(0.00319)
    label("Sinusoidal endothelium -> blood/bile loss rate constant k4 (1/min)") # S2 Table 'k4' = 0.00319 (Fig 4B 0.00319)

    # Hepatic blood-volume fraction, fixed at 20% in the fit. Alidori 2017 Compartmental modeling
    # analysis and S2 Table caption: 'Assume a constant 20% blood volume.'
    fvb <- fixed(0.20)
    label("Hepatic blood-volume fraction (dimensionless)") # S2 Table caption: constant 20% hepatic blood volume

    # Deterministic model: the paper reports no residual-error model for the compartmental fit
    # (weighted least squares over two animals). A small fixed proportional SD is supplied only
    # so the model carries a valid residual endpoint for simulation; it is not a fitted value.
    propSd <- fixed(0.1)
    label("Proportional residual SD on hepatic SUV (fraction)") # not reported; nominal (see vignette Errata)
  })

  model({
    # Whole-blood [86Y]fCNT activity is an exogenous, error-free forcing function derived from
    # the heart VOI; the model does not describe blood kinetics. Supplied in SUV at every row and
    # linearly interpolated between rows.
    linear(SUV_BLOOD)

    # Parallel (mammillary) topology, Alidori 2017 Fig 4A: blood feeds BOTH tissue compartments
    # directly, each with its own loss route (no exchange between C1 and C2):
    #   blood (SUV_BLOOD) --k1--> liver (C1) --k2--> blood/bile
    #   blood (SUV_BLOOD) --k3--> liver_lsec (C2) --k4--> blood/bile
    # k4 (sinusoidal-endothelium loss) is an order of magnitude smaller than its entrance rate
    # k3, concentrating tracer in that slow hold-up compartment.
    d/dt(liver) <- kup * SUV_BLOOD - kelim * liver
    d/dt(liver_lsec) <- kseq * SUV_BLOOD - kret * liver_lsec
    liver(0) <- 0
    liver_lsec(0) <- 0

    # Observed liver VOI SUV = fixed blood-volume contribution + the two tissue compartments.
    # The 20% blood fraction was held fixed during the fit. The observed activity (SUV) is the
    # model's single concentration-like output, so it carries the canonical name Cc
    # (units$concentration = "SUV"); there is no plasma central compartment in this tissue model.
    Cc <- fvb * SUV_BLOOD + liver + liver_lsec
    Cc ~ prop(propSd)
  })
}
