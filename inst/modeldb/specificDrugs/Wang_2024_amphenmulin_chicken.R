Wang_2024_amphenmulin_chicken <- function() {
  description <- "Preclinical (broiler chicken). One-compartment intravenous PK model for amphenmulin, a novel pleuromutilin derivative, in healthy yellow-feathered broiler chickens given a single 20 mg/kg body-weight dose. Wang 2024 analysed the plasma concentration-time data non-compartmentally in Phoenix 8.1 and reported no structural compartmental model, but the paper states that a one-compartment open model with first-order elimination was the structure used to drive the companion in-vitro dynamic system, and the intravenous NCA parameters in Table 1 determine such a model completely and self-consistently: CL = 1.17 L/h/kg and Vss = 3.64 L/kg give kel = 0.321 1/h (Table 1 T1/2Ke = 2.13 h implies 0.325 1/h) and MRT = Vss/CL = 3.11 h (Table 1 MRT = 3.13 h). ONLY the intravenous route is packaged. The oral and intramuscular routes reported in Table 1 are deliberately excluded because they cannot be represented by any one-compartment model built on these intravenous parameters: the published oral Cmax of 0.73 ug/mL is 2.26 times F*Dose/Vss = 0.323 ug/mL, which is the absolute upper bound on the concentration a one-compartment model can reach after an oral dose of 20 mg/kg at F = 5.88%, and the published oral Cmax/AUC ratio of 0.658 1/h exceeds kel = 0.321 1/h, which a one-compartment model also cannot produce. The intramuscular route is feasible in magnitude but Wang 2024 reports no intramuscular absorption rate constant, and its Tmax of 0.38 h is incompatible with its own terminal rate constant of 0.13 1/h (flip-flop absorption at 0.13 1/h would place Tmax at 4.7 h). See the vignette Assumptions and deviations section. Wang 2024 reports mean and SD of the individual NCA estimates but no population model, so no between-subject variability is encoded and the residual error is held at zero for typical-value simulation."
  reference <- paste(
    "Wang W, Yu J, Ji X, Xia X, Ding H. (2024).",
    "Pharmacokinetic/pharmacodynamic integration of amphenmulin:",
    "a novel pleuromutilin derivative against Mycoplasma gallisepticum.",
    "Microbiology Spectrum 12(2):e03675-23.",
    "doi:10.1128/spectrum.03675-23.",
    sep = " "
  )
  vignette <- "Wang_2024_amphenmulin"
  units <- list(
    time = "h",
    dosing = "mg/kg",
    concentration = "ug/mL"
  )

  compartmentData <- list(
    central = list(analyte = "amphenmulin", units = "mg/kg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "chicken (yellow-feathered broiler, healthy adult)",
    n_subjects     = 10L,
    n_studies      = 1L,
    weight_range   = "1.58-2.04 kg",
    disease_state  = "Healthy; no experimental infection. Birds were caged and fed ad libitum without antibiotics or anticoccidial drugs",
    dose_range     = "Single 20 mg/kg body-weight dose of amphenmulin (98% purity). The dose level was chosen from a previous tiamulin report",
    design         = "Thirty birds randomised to three groups of 10, receiving the same 20 mg/kg dose by intravenous injection, intramuscular injection, or oral gavage. Only the intravenous group informs this model",
    sampling       = "Brachial wing-vein blood at 0, 5, 10, 15, 30 and 45 min and 1, 2, 3, 4, 6, 8, 10, 12 and 24 h after intravenous injection; plasma assayed by HPLC-MS/MS (LLOQ 0.005 ug/mL, LOD 0.001 ug/mL, linear 0.005-0.5 ug/mL)",
    regions        = "China (South China Agricultural University, Guangzhou)",
    notes          = "Ethics approval 2023a015 (Animal Ethics Committee of South China Agricultural University). Table 1 reports the intravenous NCA results as mean +/- SD over the 10 birds: Kel 0.35 +/- 0.11 1/h, T1/2Ke 2.13 +/- 0.63 h, AUC0-t 18.76 +/- 6.09 h*ug/mL, AUC0-inf 18.89 +/- 6.03 h*ug/mL, MRT 3.13 +/- 1.05 h, CL 1.17 +/- 0.39 L/h/kg, Vss 3.64 +/- 1.58 L/kg. The absolute bioavailability derived from the AUC ratios was 52.14% intramuscular and 5.88% oral. Amphenmulin was designed and synthesised in the authors' own laboratory and is not a marketed drug."
  )

  ini({
    # =================================================================
    # Wang 2024 Table 1, intravenous column
    # =================================================================
    # CL and Vss are the two directly reported intravenous NCA
    # parameters and are taken as the model's primaries. They are the
    # mean of the per-bird NCA estimates, so they are not exactly
    # reciprocal with the separately averaged AUC: D/CL = 17.09
    # h*ug/mL against the reported mean AUC0-inf of 18.89 h*ug/mL
    # (-9.5%, well inside the reported SD of 6.03). This is the usual
    # mean-of-ratios versus ratio-of-means gap and is not an
    # inconsistency in the data. See the vignette Errata.
    lcl <- log(1.17)
    label("Clearance (L/h/kg)")  # Wang 2024 Table 1, IV column, Cl = 1.17 +/- 0.39 L/h/kg

    # Vss is used as the single-compartment volume. For a one-
    # compartment model Vss is the volume, and the choice is
    # corroborated by the reported MRT: Vss/CL = 3.11 h against the
    # reported MRT of 3.13 h, and log(2)*Vss/CL = 2.16 h against the
    # reported T1/2Ke of 2.13 h.
    lvc <- log(3.64)
    label("Volume of distribution (L/kg)")  # Wang 2024 Table 1, IV column, Vss = 3.64 +/- 1.58 L/kg

    # =================================================================
    # Residual error
    # =================================================================
    # Wang 2024 fitted no population model and reported no residual
    # error magnitude, only the mean +/- SD of the individual NCA
    # estimates. The residual SD is therefore held at zero for
    # deterministic typical-value simulation rather than inventing a
    # variance. The Table 1 SDs are reproduced in the vignette so a
    # user can add variability deliberately.
    addSd <- fixed(0)
    label("Additive residual SD on plasma concentration (ug/mL; not reported in Wang 2024)")  # Wang 2024 reports no residual error; only NCA mean +/- SD
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)

    kel <- cl / vc

    # One-compartment open model with first-order elimination
    # (Wang 2024, Discussion: the structure used for the in-vitro
    # dynamic system). Dose the `central` compartment directly; the
    # oral and intramuscular routes are not represented, see the
    # description and the vignette Assumptions and deviations.
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
