Birhanu_2020_marbofloxacin_rat <- function() {
  description <- paste0(
    "Preclinical (rat, Sprague-Dawley). One-compartment PK model for marbofloxacin (MAR) given ALONE ",
    "as a single 5 mg/kg intravenous or oral dose to healthy male rats, with an algebraic ex vivo ",
    "antibacterial-effect output against Salmonella Typhimurium (Birhanu 2020). The sibling file ",
    "Birhanu_2020_marbofloxacin_methylGallate_rat carries the same structure fitted to rats given MAR ",
    "together with 20 mg/kg methyl gallate. Disposition (V, CL) is the authors' one-compartment fit to ",
    "the intravenous data in Phoenix WinNonlin (Table 3). The oral data were analysed only ",
    "non-compartmentally (Table 2), so the oral input is reconstructed here: F = 131% is the ",
    "published bioavailability (oral AUC / intravenous AUC) and ka = 7.64 1/h is DERIVED as the root ",
    "of the one-compartment identity ln(ka/kel)/(ka - kel) = Tmax at the published Tmax of 0.5 h; it ",
    "is corroborated by the published oral Cmax (1.00 predicted vs 1.04 ug/mL) and AUC (5.96 vs ",
    "5.98 h*ug/mL). The model does not reproduce the oral terminal half-life of 7.62 h (it gives the ",
    "intravenous 3.77 h). Dose `central` for the intravenous route and `depot` for the oral route. ",
    "The output log_cfu is the predicted log10 CFU/mL of S. Typhimurium after 24 h of ex vivo ",
    "incubation (starting inoculum 10^6 CFU/mL) of a plasma sample drawn at the current time. ",
    "Birhanu 2020 fitted the sigmoid Imax model E = E0 - Imax * Ce^gamma / (IC50^gamma + Ce^gamma) ",
    "(Table 4; the Methods equation misprints the minus as a plus, refuted by Figure 4 and Table 5), ",
    "where Ce is the ex vivo AUC0-24/MIC of the sample, i.e. the sample concentration times 24 h ",
    "divided by the MIC (0.03125 ug/mL, printed as 0.031). No between-animal variability or ",
    "residual error was reported: there are no etas and the residual SDs are FIXED at 0, so the ",
    "model is for typical-value simulation."
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
    dose_range = "Single 5 mg/kg marbofloxacin intravenously or orally, alone (this file) or with 20 mg/kg methyl gallate (sibling file)",
    design = "Four-way crossover in ten rats (MAR IV, MAR PO, MAR + MG IV, MAR + MG PO) at two-week intervals",
    sampling = "300 uL tail-vein blood into EDTA at 0, 0.25, 0.5, 0.75, 1, 2, 4, 8, 12 and 24 h; plasma assayed by HPLC-UV at 293 nm (LOD 0.048 ug/mL, LOQ 0.16 ug/mL, linear 0.05-50 ug/mL, intra- and inter-assay CV < 10%)",
    organism = "Salmonella enterica serovar Typhimurium; marbofloxacin MIC 0.031 ug/mL (ATCC 14028 and field strain S-2), 0.5 ug/mL (field strain S-15)",
    regions = "Republic of Korea (Kyungpook National University; rats from OrientBio Inc.)",
    notes = paste(
      "Ethics approval 2017-0159-1 (Kyungpook National University Animal Ethics Committee).",
      "PK means +/- SD are over n = 10 rats; the ex vivo antibacterial study (Figures 3 and 4,",
      "Table 4) used plasma from n = 6 rats. Ex vivo assay: 25 uL of plasma from each sampling",
      "time added to 225 uL of cation-adjusted Mueller Hinton II broth containing 10^6 CFU/mL",
      "S. Typhimurium, incubated at 37 C and counted at 0, 1, 2, 4, 8, 12 and 24 h. The PK/PD",
      "index is computed by the authors from the undiluted plasma concentration, not the 1:10",
      "in-well concentration (Figure 4 x-axis). The sigmoid Imax parameters were fitted to the",
      "plasma of ORALLY dosed rats only."
    )
  )

  ini({
    # =================================================================
    # Disposition -- Birhanu 2020 Table 3, MAR column (intravenous,
    # one-compartment fit in Phoenix WinNonlin 8)
    # =================================================================
    # Consistency: CL / V = 1.0994 / 5.986 = 0.1837 1/h, matching the
    # printed K10 = 0.18 1/h, K10 half-life 3.77 h and MRT 5.44 h, and
    # Dose / CL = 4.548 h*ug/mL against the printed AUC of 4.55.
    lcl <- log(1.0994)
    label("Clearance (L/h/kg)") # Birhanu 2020 Table 3, MAR column, CL = 1099.4 +/- 345 mL/h/kg
    lvc <- log(5.986)
    label("Volume of distribution (L/kg)") # Birhanu 2020 Table 3, MAR column, V = 5986 +/- 620 mL/kg

    # =================================================================
    # Oral input -- Birhanu 2020 Table 2 / Results text (oral data were
    # analysed non-compartmentally; no ka was published)
    # =================================================================
    # DERIVED: ka is the root of ln(ka/kel)/(ka - kel) = Tmax with
    # kel = CL/V = 0.18366 1/h and the Table 2 Tmax of 0.5 h, giving
    # ka = 7.6397 1/h. Corroborated by Cmax = F*Dose/V*exp(-kel*Tmax)
    # = 0.998 ug/mL against the published 1.04 +/- 0.4 ug/mL, a value not
    # used in the derivation.
    lka <- log(7.6397)
    label("Oral first-order absorption rate constant (1/h)") # DERIVED from Birhanu 2020 Table 2, MAR column, Tmax = 0.5 h, with kel from Table 3
    # Published absolute bioavailability above 1; equals the AUC ratio
    # 5.98 / 4.55 = 1.314 (Table 2 AUCall / Table 3 AUC).
    lfdepot <- log(1.31)
    label("Oral bioavailability (fraction; published value exceeds 1)") # Birhanu 2020 Results: bioavailability (F) of MAR given alone = 131%

    # =================================================================
    # Ex vivo antibacterial effect -- Birhanu 2020 Table 4, MAR column
    # =================================================================
    # E = E0 - Imax * Ce^gamma / (IC50^gamma + Ce^gamma), E = change in
    # log10 CFU/mL over 24 h of ex vivo incubation. The Methods print a
    # plus sign; Figure 4A (curve falls from +4 to about -6.8) and the
    # Table 5 dosages (reproduced exactly only with the minus sign) refute
    # it. Table 5 check: the index for E = 0, -3, -4 is 368.8, 454.9,
    # 493.2 h, and index * MIC * (CL/F) / F with MIC = 0.03125 ug/mL,
    # CL/F = 0.70542 L/h/kg (Table 2) and F = 1.31 gives 6.21, 7.65 and
    # 8.30 mg/kg against the published 6.20, 7.65 and 8.3 mg/kg.
    e0 <- 4.05
    label("Change in S. Typhimurium count over 24 h in drug-free ex vivo plasma, e0 (log10 CFU/mL)") # Birhanu 2020 Table 4, MAR column, E0 = 4.05 +/- 0.29
    limax <- log(10.86)
    label("Maximum reduction from e0, imax (log10 CFU/mL)") # Birhanu 2020 Table 4, MAR column, Imax = 10.86 +/- 0.58
    lic50 <- log(405.99)
    label("Ex vivo AUC0-24/MIC at half the maximum effect, ic50 (h)") # Birhanu 2020 Table 4, MAR column, IC50 = 405.99 +/- 11.78
    lhill <- log(5.41)
    label("Hill coefficient of the ex vivo AUC0-24/MIC effect curve (unitless)") # Birhanu 2020 Table 4, MAR column, Gamma = 5.41 +/- 1.46

    # MIC that defines the PK/PD index. Table 1 prints 0.031 ug/mL; the
    # Methods dilution series starts at 1 ug/mL, so the plate value is
    # 1/32 = 0.03125 ug/mL, which is the value that reproduces the Table 5
    # dosages (0.031 gives 6.16 instead of the published 6.20 mg/kg) and
    # the Table 5 AUC24/MIC of 191.22 h (5.98 / 0.03125 = 191.4).
    mic <- fixed(0.03125)
    label("Marbofloxacin MIC against S. Typhimurium ATCC 14028 (ug/mL)") # Birhanu 2020 Table 1, MIC MAR, ATCC 14028 = 0.031 ug/mL (twofold series from 1 ug/mL, Methods)

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
