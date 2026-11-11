Liang_2024b_osimertinib_qsp <- function() {
  description <- paste(
    "QSP. Intracranial EGFR-occupancy (EO) model for osimertinib (OSI,",
    "AZD9291) in NSCLC patients with brain metastases: covalent binding of",
    "osimertinib to four EGFR variants (wild type, L858R+T790M, L858R and",
    "C797S) with re-synthesis of free EGFR, as in Liang 2024 Eqs 3-5. Four",
    "independent free-EGFR / osimertinib-EGFR-complex pairs share one driver.",
    "The driver is the free osimertinib concentration in brain tumour (BRT,",
    "the interstitial sub-compartment of the PK-Sim brain), supplied per",
    "record as the time-varying covariate CEFFECT in umol/L; it is already",
    "unbound, so no fraction-unbound factor is applied. The whole-body PK-Sim",
    "11.2 PBPK layer that generates that concentration is a platform model",
    "and is NOT extracted. Deterministic mechanism model with no IIV and no",
    "residual error. IMPORTANT DEVIATION -- the association rate constants",
    "are encoded in 1/(umol/L)/min (x60 into hours), not the 1/(umol/L)/s",
    "printed in Table 1. As printed, the model predicts about 77% C797S",
    "occupancy at the paper's own intracranial concentrations, contradicting",
    "Figure 2B (about 10%) and the paper's C797S-resistance conclusion; the",
    "per-minute reading (PK-Sim's base unit for a second-order rate",
    "constant) reproduces Figure 2B semi-quantitatively. This is an inferred",
    "unit correction; the numeric values are the paper's. See the vignette",
    "'Assumptions and deviations' section.",
    sep = " "
  )
  reference <- paste(
    "Liang F, Zhang Y, Xue Q, Zhang X. (2024). Integrated PBPK-EO modeling of",
    "osimertinib to predict plasma concentrations and intracranial EGFR",
    "engagement in patients with brain metastases. Sci Rep 14:12736.",
    "doi:10.1038/s41598-024-63743-z. PMC11148161.",
    "EGFR-occupancy equations from Methods Eqs 3-5; parameters from Table 1",
    "('EGFR occupancy' block).",
    sep = " "
  )
  vignette <- "Liang_2024b_osimertinib"

  units <- list(
    time = "h",
    dosing = "(not applicable; the free brain-tumour osimertinib concentration is supplied as the time-varying covariate CEFFECT)",
    concentration = "umol/L"
  )

  # The four EGFR variants are point mutants (or the wild type), not
  # anatomical locations, so the canonical `target` / `complex` TMDD names
  # carry a variant suffix and are declared here.
  paper_specific_compartments <- c(
    "target_wt",
    "complex_wt",
    "target_l858r_t790m",
    "complex_l858r_t790m",
    "target_l858r",
    "complex_l858r",
    "target_c797s",
    "complex_c797s"
  )

  # Liang 2024 works in molar concentrations (EGFR T0 = 0.299 umol/L), so the
  # states are concentrations rather than amounts. The occupancy is simulated
  # in brain tumour (Methods: 'the occupancy of T790M and L858R mutations in
  # BRT'), hence specimen "tissue".
  compartmentData <- list(
    target_wt = list(
      analyte = "wild-type EGFR (free, unbound) in brain tumour",
      units = "umol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_wt = list(
      analyte = "osimertinib-wild-type EGFR complex in brain tumour",
      units = "umol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    target_l858r_t790m = list(
      analyte = "EGFR L858R+T790M double mutant (free, unbound) in brain tumour",
      units = "umol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_l858r_t790m = list(
      analyte = "osimertinib-EGFR L858R+T790M complex in brain tumour",
      units = "umol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    target_l858r = list(
      analyte = "EGFR L858R mutant (free, unbound) in brain tumour",
      units = "umol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_l858r = list(
      analyte = "osimertinib-EGFR L858R complex in brain tumour",
      units = "umol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    target_c797s = list(
      analyte = "EGFR C797S mutant (free, unbound) in brain tumour",
      units = "umol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_c797s = list(
      analyte = "osimertinib-EGFR C797S complex in brain tumour",
      units = "umol/L",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    CEFFECT = list(
      description = "Time-varying FREE osimertinib concentration in brain tumour (umol/L) driving EGFR occupancy",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "For this model CEFFECT is Liang 2024's C_OSI,BRT, 'the free OSI",
        "concentration in BRT' (Eqs 3-4), which the paper takes from the",
        "free concentration in the interstitial sub-compartment of the",
        "PK-Sim brain. It is already unbound: do NOT multiply by a fraction",
        "unbound. Supply 0 for the drug-free state.",
        "In the source study CEFFECT comes from a whole-body PK-Sim 11.2",
        "PBPK model (brain-tumour-to-plasma partition coefficient KBRT,p =",
        "2.89, Table 1) that is not extracted. Published anchors for 80 mg",
        "once daily (Liang 2024 Table 3, predicted intracranial free",
        "concentration across five brain-metastasis studies): Cmax 12.7-20.2",
        "nmol/L and Ctrough 10.3-15.9 nmol/L, i.e. 0.0127-0.0202 and",
        "0.0103-0.0159 umol/L; observed CSF Ctrough 4.1-19.5 nmol/L.",
        "Users may supply CEFFECT from any osimertinib PK model; Liang 2024",
        "did not couple this layer to a published popPK model."
      ),
      source_name = "C_OSI,BRT (computed by the PK-Sim PBPK model; not a dataset column)"
    )
  )

  population <- list(
    species = "human (in silico; PK-Sim virtual populations)",
    n_subjects = NA_integer_,
    n_studies = 10L,
    age_range = "31-84 years across the simulated clinical studies (Liang 2024 Table 2)",
    sex_female_pct = NA_real_,
    race_ethnicity = c(Japanese = NA_real_, Chinese = NA_real_, White = NA_real_),
    disease_state = "non-small cell lung cancer (NSCLC) with brain metastases; EGFR wild type, L858R+T790M, L858R and C797S variants",
    dose_range = "osimertinib 20-240 mg once daily or 20-120 mg twice daily (Liang 2024 Figure 4); 80 mg once daily for 14 days in the Figure 2B occupancy simulation",
    regions = "simulated Japanese, Chinese and White populations",
    notes = paste(
      "Liang 2024 is a PBPK-EO simulation study: the virtual populations'",
      "race, dose, subject numbers, age, proportion female, BMI and albumin",
      "were taken from the clinical studies in Table 2 (Planchard, Zhao,",
      "Harvey, Grande for plasma PK; Goldstein, Yamaguchi, Leeuw, Fukuhara,",
      "Ekman for intracranial PK; 4-49 patients each). Dosing-regimen",
      "simulations used 100 virtual patients with Planchard demographics",
      "(Japanese NSCLC, 75% female, median age 62.5 y); DDI simulations used",
      "10. Disease-state physiology: hematocrit 0.33 and albumin 0.31 g/dL",
      "as printed in Table 1 (the Methods text gives 0.31 g/dL for patients",
      "and 0.45 g/dL for healthy subjects). These enter the PBPK layer only;",
      "the EO parameters extracted here are shared by every scenario, so",
      "population enters this model solely through CEFFECT."
    )
  )

  ini({
    # ---- Association rate constants, one per EGFR variant ----------------
    # Liang 2024 Table 1 'kon EGFR (uM-1 s-1)', calculated as kinact/Ki.
    # Encoded as 1/(umol/L)/min, x60 into the hourly time base -- an inferred
    # unit correction of the printed 's-1' label (see description and the
    # vignette). Fixed: literature-derived, not estimated by Liang 2024.
    lkon_wt <- fixed(log(0.028 * 60))
    label("Log association rate constant of osimertinib to wild-type EGFR (1/(umol/L)/h)")
    # Table 1 kon wild type = 0.028, read per min: 1.68 /uM/h.

    lkon_l858r_t790m <- fixed(log(1.40 * 60))
    label("Log association rate constant of osimertinib to EGFR L858R+T790M (1/(umol/L)/h)")
    # Table 1 kon L858R + T790M = 1.40, read per min: 84 /uM/h.

    lkon_l858r <- fixed(log(0.57 * 60))
    label("Log association rate constant of osimertinib to EGFR L858R (1/(umol/L)/h)")
    # Table 1 kon L858R = 0.57, read per min: 34.2 /uM/h.

    lkon_c797s <- fixed(log(0.0026 * 60))
    label("Log association rate constant of osimertinib to EGFR C797S (1/(umol/L)/h)")
    # Table 1 kon C797S = 0.0026, read per min: 0.156 /uM/h.

    # ---- Dissociation rate constant --------------------------------------
    lkoff <- fixed(log(0.001))
    label("Log dissociation rate constant of the osimertinib-EGFR complex (1/h)")
    # Table 1 koff EGFR = 0.001 h-1; Methods: 'a small finite value at 0.001 h-1 is assigned to koff'.

    # ---- EGFR turnover ----------------------------------------------------
    lkdeg <- fixed(log(0.025))
    label("Log degradation (and re-synthesis) rate constant of free EGFR (1/h)")
    # Table 1 kdeg EGFR = 0.025 h-1 (half-life 27.5 h in cell experiments).

    # ---- Baseline EGFR concentration -------------------------------------
    lrbase <- fixed(log(0.299))
    label("Log baseline (initial) EGFR concentration EGFR0 (umol/L)")
    # Table 1 EGFR T0 = 0.299 uM.
  })

  model({
    # 1. Back-transform the fixed structural parameters.
    kon_wt <- exp(lkon_wt)
    kon_l858r_t790m <- exp(lkon_l858r_t790m)
    kon_l858r <- exp(lkon_l858r)
    kon_c797s <- exp(lkon_c797s)
    koff <- exp(lkoff)
    kdeg <- exp(lkdeg)
    rbase <- exp(lrbase)

    # 2. Liang 2024 Eqs 3-4, transcribed as printed, one pair per variant:
    #      dEO/dt       = kon * C_OSI,BRT * EGFRfree - koff * EO
    #      dEGFRfree/dt = (EGFR0 - EGFRfree) * kdeg
    #                     - kon * C_OSI,BRT * EGFRfree + koff * EO
    #    CEFFECT is C_OSI,BRT, already the free concentration. The complex
    #    has no degradation term in the paper, so total EGFR (free + complex)
    #    rises above EGFR0 under exposure; that is the printed structure.
    d/dt(target_wt) <- (rbase - target_wt) * kdeg -
      kon_wt * CEFFECT * target_wt + koff * complex_wt
    d/dt(complex_wt) <- kon_wt * CEFFECT * target_wt - koff * complex_wt

    d/dt(target_l858r_t790m) <- (rbase - target_l858r_t790m) * kdeg -
      kon_l858r_t790m * CEFFECT * target_l858r_t790m +
      koff * complex_l858r_t790m
    d/dt(complex_l858r_t790m) <- kon_l858r_t790m * CEFFECT * target_l858r_t790m -
      koff * complex_l858r_t790m

    d/dt(target_l858r) <- (rbase - target_l858r) * kdeg -
      kon_l858r * CEFFECT * target_l858r + koff * complex_l858r
    d/dt(complex_l858r) <- kon_l858r * CEFFECT * target_l858r -
      koff * complex_l858r

    d/dt(target_c797s) <- (rbase - target_c797s) * kdeg -
      kon_c797s * CEFFECT * target_c797s + koff * complex_c797s
    d/dt(complex_c797s) <- kon_c797s * CEFFECT * target_c797s -
      koff * complex_c797s

    # 3. Drug-free start: all EGFR free at EGFR0, no complex.
    target_wt(0) <- rbase
    target_l858r_t790m(0) <- rbase
    target_l858r(0) <- rbase
    target_c797s(0) <- rbase

    # 4. EGFR occupancy (%), Liang 2024 Eq 5:
    #      occupancy = (EGFR0 - EGFRfree) / EGFR0 * 100
    #    i.e. depletion of free EGFR relative to baseline, not the bound
    #    fraction of total EGFR.
    ro_wt <- (rbase - target_wt) / rbase * 100
    ro_l858r_t790m <- (rbase - target_l858r_t790m) / rbase * 100
    ro_l858r <- (rbase - target_l858r) / rbase * 100
    ro_c797s <- (rbase - target_c797s) / rbase * 100
  })
}
