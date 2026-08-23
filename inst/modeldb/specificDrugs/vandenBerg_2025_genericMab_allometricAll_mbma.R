vandenBerg_2025_genericMab_allometricAll_mbma <- function() {
  description <- paste0(
    "MBMA. Generic mAb two-compartment intravenous population PK model with body-weight ",
    "allometry on ALL FOUR disposition parameters -- 'covariate combination II', the second ",
    "most common weight-covariate configuration in the literature (15 of 160 collected ",
    "models). The disposition parameters are the across-model medians of a systematic ",
    "meta-analysis of 160 published two-compartment population PK models of 69 marketed IgG ",
    "monoclonal antibodies (CL 0.22 L/day, Vc 3.42 L, Vp 2.68 L, Q 0.54 L/day), and the weight ",
    "exponents are the medians within that subset (0.75 on CL, 0.85 on Vc, 0.85 on Vp, 0.77 on ",
    "Q), normalised to a 70 kg reference. Scaling the peripheral parameters too makes terminal ",
    "half-life DECREASE as weight falls (19 days at 40 kg), the opposite of the sibling model ",
    "vandenBerg_2025_genericMab_allometricClVc_mbma (24.4 days), which is the comparison the ",
    "paper draws; single-dose AUC nonetheless differs by only ~10% between the two and the ",
    "difference vanishes at steady state. Deterministic by construction: the authors ran it ",
    "without inter-individual or residual variability. Values are van den Berg 2025 Figure 4h ",
    "and the Supplementary Figure S7 caption."
  )

  reference <- paste(
    "van den Berg SPH, Adolfsen PEA, Dorlo TPC, Rispens T.",
    "Does one model fit all mAbs? An evaluation of population pharmacokinetic models.",
    "mAbs. 2025;17(1):2512217.",
    "doi:10.1080/19420862.2025.2512217",
    sep = " "
  )

  vignette <- "vandenBerg_2025_genericMab"

  units <- list(
    time          = "day",
    dosing        = "mg (intravenous into central; the paper's simulation used a single or three-weekly 100 mg dose)",
    concentration = "mg/L (equivalently ug/mL; central amount in mg divided by Vc in L)"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric power scaling of CL, Vc, Vp and Q, normalised to a 70 kg reference:",
        "P_WT = P_ref * (WT / 70)^exponent, for which the source gives the worked example",
        "CL_40kg = 0.22 * (40/70)^0.75 (van den Berg 2025 Supplementary Figure S7 caption).",
        "The simulated weight range was 40-100 kg."
      ),
      source_name        = "WT"
    )
  )

  compartmentData <- list(
    central     = list(analyte = "genericMab", units = "mg", specimen = "serum", verified = FALSE),
    peripheral1 = list(analyte = "genericMab", units = "mg", specimen = "serum", verified = FALSE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 143094L,
    n_studies      = 160L,
    n_models       = 160L,
    n_mabs         = 69L,
    weight_range   = "40-100 kg simulated (van den Berg 2025 Figure 4h and Figure S7B-E)",
    disease_state  = paste(
      "Mixed; the pooled evidence base spans every indication in which a marketed canonical IgG",
      "mAb has had a population PK model published, plus healthy donors. Paediatric-only models",
      "were excluded (van den Berg 2025 Figure 1a)."
    ),
    dose_range     = "100 mg intravenous, single dose or every 3 weeks (van den Berg 2025 Figure S7C-E)",
    regions        = "Global; the collected models were not stratified by region.",
    design         = paste(
      "Meta-analysis of published population PK model PARAMETER ESTIMATES. The weight exponents",
      "here are the medians of the 15 collected models that placed body weight on all of CL, Vc,",
      "Vp and Q; the disposition parameters are the medians over all 160 models",
      "(van den Berg 2025 Results 'Covariate model' and Supplementary Figure S7 caption)."
    ),
    reference_subject = "A 70 kg individual, at whom CL, Vc, Vp and Q equal the across-model medians.",
    notes = paste(
      "n_subjects is the SUM of per-model patient counts across the 157 of 160 models reporting",
      "one (Supplementary Table S1 column N), with patients double-counted across models fitted",
      "to overlapping datasets. This configuration was simulated deterministically -- no",
      "inter-individual and no residual variability (van den Berg 2025 Figure S7 caption)."
    )
  )

  ini({
    # =====================================================================
    # Covariate combination II of van den Berg 2025: body weight on CL, Vc,
    # Vp and Q. Disposition parameters are the across-model medians (identical
    # to vandenBerg_2025_genericMab_mbma); the exponents are the medians within
    # the subset of models that used this configuration. Nothing here was
    # estimated by these authors -- every value is a meta-analytic median held
    # constant for simulation, so every value is fixed().
    # =====================================================================

    # ----- Structural disposition parameters at the 70 kg reference -----
    lcl <- fixed(log(0.22)); label("Clearance CL at 70 kg (L/day)")                    # van den Berg 2025 Suppl. Figure S7 caption "CL = 0.22 L/d"; Results: median (IQR) 0.22 (0.17; 0.29) L/d
    lvc <- fixed(log(3.42)); label("Central volume of distribution Vc at 70 kg (L)")   # van den Berg 2025 Suppl. Figure S7 caption "VC = 3.42 L"; Results: median (IQR) 3.42 (2.96; 3.99) L
    lvp <- fixed(log(2.68)); label("Peripheral volume of distribution Vp at 70 kg (L)")# van den Berg 2025 Suppl. Figure S7 caption "VP = 2.68 L"; Results: median (IQR) 2.68 (2.09; 3.54) L
    lq  <- fixed(log(0.54)); label("Intercompartmental clearance Q at 70 kg (L/day)")  # van den Berg 2025 Suppl. Figure S7 caption "Q = 0.54 L/d"; Results: median (IQR) 0.54 (0.36; 0.84) L/d

    # ----- Body-weight allometric exponents, combination II -----
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL (unitless)")  # van den Berg 2025 Suppl. Figure S7 caption: "II) CL_WT = 0.75, Vc_WT = 0.85, Vp_WT = 0.85, Q_WT = 0.77"; the same caption's worked example uses 0.75 for CL
    e_wt_vc <- fixed(0.85); label("Allometric exponent of body weight on Vc (unitless)")  # van den Berg 2025 Suppl. Figure S7 caption, combination II
    e_wt_vp <- fixed(0.85); label("Allometric exponent of body weight on Vp (unitless)")  # van den Berg 2025 Suppl. Figure S7 caption, combination II
    e_wt_q  <- fixed(0.77); label("Allometric exponent of body weight on Q (unitless)")   # van den Berg 2025 Suppl. Figure S7 caption, combination II

    # =====================================================================
    # No inter-individual variability and no residual variability: the
    # Figure 4h / S7 simulations were run for a typical individual only
    # ("This simulation was performed with the generic model ... without
    # inter-individual variability or residual variability"). propSd is a
    # structural zero so the model has a well-formed endpoint; it is not an
    # estimate and must be re-specified before fitting data.
    # =====================================================================
    propSd <- fixed(0); label("Proportional residual error SD (fraction; structural zero, no residual variability was modelled)")  # van den Berg 2025 Suppl. Figure S7 caption: "without inter-individual variability or residual variability"
  })

  model({
    # Allometric scaling on all four disposition parameters, normalised to
    # 70 kg (Suppl. Figure S7 caption: P_40kg = P_ref * (WT/70)^P_WT).
    cl <- exp(lcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp) * (WT / 70)^e_wt_vp
    q  <- exp(lq)  * (WT / 70)^e_wt_q

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
