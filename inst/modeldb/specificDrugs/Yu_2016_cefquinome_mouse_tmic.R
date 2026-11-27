Yu_2016_cefquinome_mouse_tmic <- function() {
  description <- paste(
    "Preclinical (mouse, lactating female CD-1; Staphylococcus aureus",
    "mastitis). Cefquinome PK/PD model of Yu 2016 driven by %T>MIC.",
    "A one-compartment model with first-order absorption from the",
    "mammary glands describes plasma cefquinome after intramammary",
    "infusion into the L4 and R4 glands (Yu 2016 Results, 'Plasma PK of",
    "cefquinome'); the time the plasma concentration exceeds the MIC",
    "during the first 24 h of treatment (t_above_mic), as a percentage of",
    "24 h, gives the PK/PD index, and the inhibitory sigmoid Emax model of",
    "Yu 2016 equation 2 (Table 2, %T>MIC column) gives the change in gland",
    "S. aureus count over the",
    "24 h of treatment (dlog10cfu, log10 CFU per gland, positive = net",
    "growth). log_cfu adds that change to the mean pre-treatment burden",
    "of 7.61 log10 CFU per gland. Both PD outputs are the paper's 24 h",
    "endpoint and are meant to be read at t = 24 h after the first dose;",
    "before 24 h they reflect only the exposure accumulated so far. The",
    "paper's own parameter names are inverted relative to the usual",
    "convention (its 'Emax' is the drug-free growth, its 'E0' the",
    "maximal kill); this file maps them onto the canonical e0 (effect at",
    "zero exposure) = +2.94 and emax (maximal drug effect) = -2.11 log10",
    "CFU per gland without changing any value. Dosing amounts are the",
    "TOTAL dose per mouse in ug (two glands, so twice the per-gland dose",
    "the paper quotes). ka and kel are the means of the per-dose-group",
    "rate constants of Yu 2016 Table 1; V/F, which the paper does not",
    "print, is recovered from the AUC/MIC column of S2 Table, which the",
    "authors computed from these same parameters (see the vignette).",
    "MIC 0.5 ug/mL (S. aureus Newbould 305 and the MIC90 of 38 bovine",
    "isolates). Plasma concentrations are total drug (protein binding 17",
    "percent in mouse). No between-animal variability or residual error",
    "was reported, so the model is typical-value only and addSd is held",
    "at 0. The parallel AUC0-24/MIC fit is Yu_2016_cefquinome_mouse_aucmic."
  )
  reference <- paste(
    "Yu Y, Zhou YF, Chen MR, Li X, Qiao GL, Sun J, Liao XP, Liu YH. (2016).",
    "In vivo pharmacokinetics/pharmacodynamics of cefquinome in an",
    "experimental mouse model of Staphylococcus aureus mastitis following",
    "intramammary infusion. PLoS One 11(5):e0156273.",
    "doi:10.1371/journal.pone.0156273.",
    sep = " "
  )
  vignette <- "Yu_2016_cefquinome"
  units <- list(
    time = "h",
    dosing = "ug (total cefquinome per mouse; both L4 and R4 glands)",
    concentration = "ug/mL (plasma cefquinome, Cc); log10 CFU per gland (log_cfu, dlog10cfu)"
  )

  paper_specific_compartments <- c("t_above_mic")

  compartmentData <- list(
    depot = list(analyte = "cefquinome", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cefquinome", units = "ug", specimen = "plasma", verified = TRUE),
    t_above_mic = list(analyte = "cefquinome", units = "h", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "mouse (CD-1, lactating female)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    weight_range = "35-45 g",
    sex_female_pct = 100,
    organism = paste(
      "Staphylococcus aureus Newbould 305 (ATCC 29740), a bovine mastitis",
      "isolate; cefquinome MIC 0.5 ug/mL by CLSI broth microdilution.",
      "MIC50 and MIC90 of 38 clinical bovine mastitis isolates from Inner",
      "Mongolia, China were both 0.5 ug/mL (range 0.25-0.5 ug/mL, S1 Table)"
    ),
    disease_state = paste(
      "Mouse model of S. aureus mastitis: lactating CD-1 mice 10-12 days",
      "post partum, offspring removed 1 h before inoculation of about 4 x",
      "10^3 CFU into each of the L4 and R4 glands through the teat canal;",
      "treatment started 9 h later at a mean burden of 7.61 log10 CFU per",
      "gland. PK was measured in healthy lactating mice."
    ),
    dose_range = paste(
      "PK: single intramammary doses of 25, 50, 100, 200 and 400 ug per",
      "gland into both L4 and R4 (6 mice per dose). PD: 21 regimens of",
      "12.5-800 ug per gland every 8, 12 or 24 h for 24 h (4 mice, 8",
      "glands, per regimen)"
    ),
    regions = "China (South China Agricultural University, Guangzhou)",
    notes = paste(
      "PK fitted per dose group in WinNonlin 5.2.1 (one-compartment,",
      "first-order absorption, chosen by AIC); the PK/PD indices of the",
      "21 regimens were simulated by the authors from the single-dose PK",
      "(Yu 2016 equation 1; S2 Table). The sigmoid Emax fit of %T>MIC",
      "had R^2 = 0.8466. Published %T>MIC targets: 7.59 for",
      "bacteriostasis and 35.98 for a 1.8 log10 CFU per gland reduction."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # Plasma PK: one compartment, first-order absorption from the glands
    # -----------------------------------------------------------------
    # Yu 2016 fitted each dose group separately and reports only the
    # half-lives in Table 1. ka and kel are the means of the per-group
    # rate constants ln(2)/T1/2 (the authors averaged rate constants:
    # the S2 Table %T>MIC column rises by ln(2)/kel = 0.389 h per dose
    # doubling, i.e. kel = 1.781 1/h, not ln(2)/0.40 = 1.733 1/h).
    lka <- log(9.727)
    label("Absorption rate constant from the mammary glands ka (1/h)") # Table 1 'T1/2Ka (h)' 0.07, 0.07, 0.11, 0.08, 0.05: mean of ln(2)/T1/2Ka = 9.727; S2 Table 'Cmax/MIC' column implies 9.747
    lkel <- log(1.7803)
    label("Elimination rate constant kel (1/h)") # Table 1 'T1/2Kel (h)' 0.49, 0.41, 0.34, 0.35, 0.39: mean of ln(2)/T1/2Kel = 1.7803
    lvc <- log(27.66)
    label("Apparent volume of distribution V/F, per total dose to both glands (mL)") # back-solved: S2 Table 'AUC/MIC (h)' 194.97 at 800 ug/gland q8h -> AUC per ug/gland = 0.0406188 h/mL; F/V = 0.0406188 * kel = 0.072314 /mL per ug/gland; V/F = 2 / 0.072314 = 27.66 mL per mouse

    # Isolate MIC, a measured input to the index.
    mic <- fixed(0.5)
    label("Cefquinome MIC of S. aureus Newbould 305 by CLSI broth microdilution (ug/mL)") # Results 'MICs of cefquinome...': 0.5 ug/mL against Newbould 305; MIC90 of 38 isolates also 0.5 ug/mL

    # -----------------------------------------------------------------
    # Inhibitory sigmoid Emax on %T>MIC (equation 2, Table 2)
    # -----------------------------------------------------------------
    # Printed form: E = Emax - (Emax - E0) * Ce^N / (EC50^N + Ce^N), with
    # the paper's Emax the 24 h change in drug-free glands (+2.94) and its
    # E0 the change at maximal drug effect (-2.11). Canonical roles:
    #   e0   <- paper 'Log Emax' (effect at zero exposure)
    #   emax <- paper 'Log E0'   (maximal drug effect)
    # Both are signed log10 changes and stay on the natural scale.
    e0 <- 2.94
    label("Change in gland S. aureus count over 24 h without drug (log10 CFU per gland)") # Table 2, %T>MIC column, 'Log Emax (logCFU/gland)' = 2.94 +/- 0.19
    emax <- -2.11
    label("Change in gland S. aureus count over 24 h at maximal drug effect (log10 CFU per gland)") # Table 2, %T>MIC column, 'Log E0 (logCFU/gland)' = -2.11 +/- 0.22
    lec50 <- log(6.12)
    label("Percent of the first 24 h with plasma cefquinome above the MIC producing half the maximal effect EC50 (unitless)") # Table 2, %T>MIC column, 'EC50 (h)' = 6.12 +/- 0.49 (a percentage of 24 h despite the printed unit)
    lhill <- log(1.54)
    label("Hill coefficient N of the time-above-MIC effect curve (unitless)") # Table 2, %T>MIC column, 'Slope (N)' = 1.54 +/- 0.19

    # Pre-treatment burden, a design input.
    log10_cfu0 <- fixed(7.61)
    label("Gland S. aureus burden at the start of treatment (log10 CFU per gland)") # Results 'PD effectiveness of 21 dosing regimens': 7.61 logCFU/gland (mean value) before administration

    addSd <- fixed(0)
    label("Additive residual SD on log_cfu (log10 CFU per gland)") # not reported (R^2 only, Figure 4)
  })

  model({
    ka <- exp(lka)
    kel <- exp(lkel)
    vc <- exp(lvc)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    # Total plasma cefquinome (ug/mL).
    Cc <- central / vc

    # Time (h) with plasma cefquinome above the MIC during the first 24 h
    # after the first dose, as a percentage of 24 h; constant after 24 h,
    # so the index holds its 24 h value from then on.
    d/dt(t_above_mic) <- (Cc > mic) * (t < 24)
    pct_tmic <- 100 * t_above_mic / 24

    # Signed change in log10 CFU per gland over the 24 h of treatment:
    # e0 at zero exposure, approaching emax at high exposure.
    dlog10cfu <- e0 - (e0 - emax) * pct_tmic^hill / (pct_tmic^hill + ec50^hill)
    log_cfu <- log10_cfu0 + dlog10cfu
    log_cfu ~ add(addSd)
  })
}
