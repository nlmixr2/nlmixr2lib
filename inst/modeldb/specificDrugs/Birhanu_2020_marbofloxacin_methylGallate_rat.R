Birhanu_2020_marbofloxacin_methylGallate_rat <- function() {
  description <- paste0(
    "Preclinical (rat, Sprague-Dawley). One-compartment PK model for marbofloxacin (MAR) given as a ",
    "single 5 mg/kg intravenous or oral dose TOGETHER WITH 20 mg/kg methyl gallate (MG) by the same ",
    "route to healthy male rats, with an algebraic ex vivo antibacterial-effect output against ",
    "Salmonella Typhimurium (Birhanu 2020). The sibling file Birhanu_2020_marbofloxacin_rat carries ",
    "the same structure fitted to rats given MAR alone. Methyl gallate itself is not modelled: its ",
    "plasma concentrations were below the assay's reach and only marbofloxacin was measured. ",
    "Disposition (V, CL) is the authors' one-compartment fit to the intravenous data in Phoenix ",
    "WinNonlin (Table 3). The oral data were analysed only non-compartmentally (Table 2), so the oral ",
    "input is reconstructed here: F = 89.1% is the published bioavailability (oral AUC / intravenous ",
    "AUC) and ka = 1.085 1/h is DERIVED as the root of ln(ka/kel)/(ka - kel) = Tmax at the published ",
    "Tmax of 2 h; it is corroborated by the published oral Cmax (0.757 predicted vs 0.78 ug/mL) and ",
    "AUC (6.11 vs 6.11 h*ug/mL). The model does not reproduce the oral terminal half-life of 6.01 h ",
    "(it gives the intravenous 3.93 h). Dose `central` for the intravenous route and `depot` for the ",
    "oral route. The output log_cfu is the predicted log10 CFU/mL of S. Typhimurium after 24 h of ex ",
    "vivo incubation (starting inoculum 10^6 CFU/mL) of a plasma sample drawn at the current time, ",
    "from the sigmoid Imax model E = E0 - Imax * Ce^gamma / (IC50^gamma + Ce^gamma) of Table 4 (the ",
    "Methods equation misprints the minus as a plus), where Ce is the sample concentration times ",
    "24 h divided by an MIC of 0.016 ug/mL. That MIC is not printed: it is DERIVED from the paper's ",
    "own Table 5 indices (AUC24/MIC 381.93 = 6.11 / 0.016) and the Figure 4B x-axis, and it differs ",
    "from the 0.008 ug/mL the paper later uses in its Table 5 dosage calculation (see the vignette ",
    "Errata). No between-animal variability or residual error was reported: there are no etas and ",
    "the residual SDs are FIXED at 0, so the model is for typical-value simulation."
  )
  reference <- paste(
    "Birhanu BT, Lee EB, Park SC (2020).",
    "Evaluation of the pharmacokinetic-pharmacodynamic integration of marbofloxacin",
    "in combination with methyl gallate against Salmonella Typhimurium in rats.",
    "PLoS One 15(6):e0234211. doi:10.1371/journal.pone.0234211. PMCID: PMC7272065.",
    sep = " "
  )
  vignette <- "Birhanu_2020_marbofloxacin"
  units <- list(
    time = "h",
    dosing = "mg/kg",
    concentration = "ug/mL"
  )

  compartmentData <- list(
    depot = list(analyte = "marbofloxacin", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "marbofloxacin", units = "mg/kg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = 10L,
    n_studies = 1L,
    sex_female_pct = 0,
    age_range = "five weeks old at purchase, acclimated for one week before dosing",
    weight_range = "180-215 g (average weight range)",
    disease_state = "Healthy (clinically healthy, uninfected rats); the antibacterial effect was measured ex vivo",
    dose_range = "Single 5 mg/kg marbofloxacin with 20 mg/kg methyl gallate, intravenously or orally (this file), or marbofloxacin alone (sibling file)",
    co_medication = "Methyl gallate 20 mg/kg given with every marbofloxacin dose by the same route; not quantified in plasma",
    design = "Four-way crossover in ten rats (MAR IV, MAR PO, MAR + MG IV, MAR + MG PO) at two-week intervals",
    sampling = "300 uL tail-vein blood into EDTA at 0, 0.25, 0.5, 0.75, 1, 2, 4, 8, 12 and 24 h; plasma assayed by HPLC-UV at 293 nm (LOD 0.048 ug/mL, LOQ 0.16 ug/mL, linear 0.05-50 ug/mL, intra- and inter-assay CV < 10%)",
    organism = "Salmonella enterica serovar Typhimurium; marbofloxacin MIC 0.031 ug/mL alone and 0.008 ug/mL in combination with 125 ug/mL methyl gallate (ATCC 14028 and S-2; FIC index 0.5); methyl gallate MIC 500 ug/mL",
    regions = "Republic of Korea (Kyungpook National University; rats from OrientBio Inc.)",
    notes = paste(
      "Ethics approval 2017-0159-1 (Kyungpook National University Animal Ethics Committee).",
      "PK means +/- SD are over n = 10 rats; the ex vivo antibacterial study (Figures 3 and 4,",
      "Table 4) used plasma from n = 6 rats. Ex vivo assay: 25 uL of plasma from each sampling",
      "time added to 225 uL of cation-adjusted Mueller Hinton II broth containing 10^6 CFU/mL",
      "S. Typhimurium, incubated at 37 C and counted at 0, 1, 2, 4, 8, 12 and 24 h. Plasma from",
      "these rats carries whatever methyl gallate was absorbed, so the Table 4 combination",
      "parameters describe the plasma's combined activity indexed to marbofloxacin alone. The",
      "sigmoid Imax parameters were fitted to the plasma of ORALLY dosed rats only."
    )
  )

  ini({
    # =================================================================
    # Disposition -- Birhanu 2020 Table 3, MGM column (intravenous,
    # one-compartment fit in Phoenix WinNonlin 8)
    # =================================================================
    # Consistency: CL / V = 0.7291 / 4.1343 = 0.1764 1/h, matching the
    # printed K10 = 0.18 1/h, K10 half-life 3.93 h and MRT 5.67 h, and
    # Dose / CL = 6.858 h*ug/mL against the printed AUC of 6.86.
    lcl <- log(0.7291)
    label("Clearance (L/h/kg)") # Birhanu 2020 Table 3, MGM column, CL = 729.1 +/- 192 mL/h/kg
    lvc <- log(4.1343)
    label("Volume of distribution (L/kg)") # Birhanu 2020 Table 3, MGM column, V = 4134.3 +/- 358 mL/kg

    # =================================================================
    # Oral input -- Birhanu 2020 Table 2 / Results text (oral data were
    # analysed non-compartmentally; no ka was published)
    # =================================================================
    # DERIVED: ka is the root of ln(ka/kel)/(ka - kel) = Tmax with
    # kel = CL/V = 0.17635 1/h and the Table 2 Tmax of 2 h, giving
    # ka = 1.0846 1/h. Corroborated by Cmax = F*Dose/V*exp(-kel*Tmax)
    # = 0.757 ug/mL against the published 0.78 +/- 0.2 ug/mL, a value not
    # used in the derivation.
    lka <- log(1.0846)
    label("Oral first-order absorption rate constant (1/h)") # DERIVED from Birhanu 2020 Table 2, MGM column, Tmax = 2 h, with kel from Table 3
    # Equals the AUC ratio 6.11 / 6.86 = 0.891 (Table 2 AUCall / Table 3 AUC).
    lfdepot <- log(0.891)
    label("Oral bioavailability (fraction)") # Birhanu 2020 Results: bioavailability (F) of MAR given with MG = 89.1%

    # =================================================================
    # Ex vivo antibacterial effect -- Birhanu 2020 Table 4, MGM column
    # =================================================================
    # E = E0 - Imax * Ce^gamma / (IC50^gamma + Ce^gamma), E = change in
    # log10 CFU/mL over 24 h of ex vivo incubation (sign as in the MAR
    # sibling; Figure 4B falls from about +3.6 to about -7).
    # The Results text prints IC50 = 504.8; Table 4 prints 584.75, which
    # Figure 4B (half-effect near 580 h) and the Table 5 dosages confirm:
    # the index for E = 0, -3, -4 is 546.7, 613.9, 640.4 h, and
    # index * 0.008 ug/mL * (CL/F = 0.75398 L/h/kg) / (F = 0.891) gives
    # 3.70, 4.16 and 4.34 mg/kg, exactly the published values, whereas
    # IC50 = 504.8 gives 3.19, 3.59 and 3.74 mg/kg.
    e0 <- 3.6
    label("Change in S. Typhimurium count over 24 h in drug-free ex vivo plasma, e0 (log10 CFU/mL)") # Birhanu 2020 Table 4, MGM column, E0 = 3.6 +/- 0.67
    limax <- log(10.66)
    label("Maximum reduction from e0, imax (log10 CFU/mL)") # Birhanu 2020 Table 4, MGM column, Imax = 10.66 +/- 1.14
    lic50 <- log(584.75)
    label("Ex vivo AUC0-24/MIC at half the maximum effect, ic50 (h)") # Birhanu 2020 Table 4, MGM column, IC50 = 584.75 +/- 30.06
    lhill <- log(10)
    label("Hill coefficient of the ex vivo AUC0-24/MIC effect curve (unitless)") # Birhanu 2020 Table 4, MGM column, Gamma = 10 +/- 5.51

    # MIC that defines the combination's PK/PD index. DERIVED, not printed:
    # Table 5 AUC24/MIC = 381.93 h = 6.11 / 0.016 and Cmax/MIC = 49.06
    # = 0.785 / 0.016 for this arm, and the Figure 4B x-axis reaches about
    # 1170 h = 0.78 ug/mL * 24 h / 0.016. Neither the single-agent MIC
    # (0.031) nor the combination checkerboard concentration (0.008)
    # reproduces those three numbers.
    mic <- fixed(0.016)
    label("Marbofloxacin MIC used for the combination's PK/PD index (ug/mL)") # DERIVED from Birhanu 2020 Table 5, MGM column, AUC24/MIC = 381.93 with Table 2 AUCall = 6.11

    # Ex vivo starting inoculum, a design input.
    log10_cfu0 <- fixed(6)
    label("Starting S. Typhimurium density of the ex vivo incubation (log10 CFU/mL)") # Birhanu 2020 Methods 'Ex vivo antibacterial effect of plasma samples': final concentration 10^6 CFU/mL

    # =================================================================
    # Residual error -- not reported
    # =================================================================
    addSd <- fixed(0)
    label("Additive residual SD on plasma marbofloxacin (ug/mL; not reported)") # Birhanu 2020 reports only mean +/- SD of per-rat estimates
    addSd_log_cfu <- fixed(0)
    label("Additive residual SD on ex vivo log10 CFU/mL (not reported)") # Birhanu 2020 Table 4 gives parameter SDs only, no residual error
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    ka <- exp(lka)
    kel <- cl / vc

    imax <- exp(limax)
    ic50 <- exp(lic50)
    hill <- exp(lhill)

    # One-compartment model. Dose `central` for IV, `depot` for oral; the
    # bioavailability applies only to the depot input.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    f(depot) <- exp(lfdepot)

    Cc <- central / vc

    # Ex vivo PK/PD index of a plasma sample drawn now: the sample is
    # incubated for 24 h at its fixed concentration, so its AUC0-24 is
    # Cc * 24 h (Birhanu 2020 Methods, 'Ce is the AUC0-24/MIC of the
    # plasma'; Figure 4 x-axis in h). The floor at zero keeps a
    # sub-atol integrator undershoot of Cc from turning the non-integer
    # power below into NaN.
    aucmic_exvivo <- max(Cc, 0) * 24 / mic
    dlog10cfu <- e0 - imax * aucmic_exvivo^hill / (aucmic_exvivo^hill + ic50^hill)
    # Predicted ex vivo count after 24 h of incubation of that sample.
    log_cfu <- log10_cfu0 + dlog10cfu

    Cc ~ add(addSd)
    log_cfu ~ add(addSd_log_cfu)
  })
}
