vanderMeijden_2023_polymyxinB_EcoliATCC25922_fractionation <- function() {
  description <- paste(
    "Preclinical (mouse, neutropenic female CD-1; thigh infection with",
    "Escherichia coli ATCC 25922, polymyxin B MIC 1 mg/L). Polymyxin B",
    "PK/PD model from the dose-fractionation study of van der Meijden 2023",
    "(Figure 2). The one-compartment subcutaneous PK model of",
    "vanderMeijden_2023_polymyxinB, at its typical values, gives the",
    "unbound AUC over the first 24 h of treatment (fAUC0-24, unbound",
    "fraction 0.2); an inhibitory sigmoid Emax model of fAUC0-24/MIC then",
    "gives the change in thigh bacterial burden between the start of",
    "treatment and 24 h (dlog10cfu, log10 cfu per thigh; positive = net",
    "growth), and the observed output log_cfu adds that change to the",
    "burden at the start of therapy (7.18 log10 cfu per thigh, the 1.5 x",
    "10^7 cfu average of the dose-fractionation studies). Both are the",
    "paper's 24 h endpoint and are meant to be read at t = 24 h; before 24",
    "h they reflect only the exposure accumulated so far. The paper prints",
    "R^2 (0.9075) but no Emax coefficients; the maintainers digitised the",
    "fitted curve of the fAUC/MIC panel of Figure 2 from the 1400-pixel",
    "raster in the published PDF; the recovered curve overlays the drawn",
    "line to within its width. No residual error was reported, so addSd is",
    "held at 0."
  )
  reference <- paste(
    "van der Meijden A, Aranzana-Climent V, van der Spek H, de Winter BCM,",
    "Couet W, Meletiadis J, Muller AE, van den Berg S. Pharmacokinetic and",
    "pharmacodynamic properties of polymyxin B in Escherichia coli and",
    "Klebsiella pneumoniae murine infection models. J Antimicrob Chemother.",
    "2023;78(3):832-839. doi:10.1093/jac/dkad022"
  )
  vignette <- "vanderMeijden_2023_polymyxinB"
  units <- list(
    time = "h",
    dosing = "mg/kg",
    concentration = "mg/L (Cc); log10 cfu per thigh (log_cfu)"
  )

  paper_specific_compartments <- c("fauc_0_24")

  covariateData <- list(
    DOSE = list(
      description = "Polymyxin B dose given at each subcutaneous administration",
      units = "mg/kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Dose per administration in mg/kg (the model is parameterised for",
        "a virtual 1 kg mouse). Enters CL/F and ka as power functions of",
        "DOSE / 10 mg/kg (Table 1 footnote). For a regimen of equal doses",
        "it is the size of each dose, not the total daily dose. Place the",
        "event columns before DOSE in the event table."
      ),
      source_name = "DOSE (polymyxin B dose in mg/kg; Table 1 footnote)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "polymyxin B", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "polymyxin B", units = "mg/kg", specimen = "plasma", verified = TRUE),
    fauc_0_24 = list(analyte = "polymyxin B", units = "mg*h/L", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "mouse (CD-1, outbred, female; neutropenic)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "7-8 weeks",
    weight_range = "20-25 g",
    sex_female_pct = 100,
    organism = "Escherichia coli ATCC 25922 (no ESBL); polymyxin B MIC 1 mg/L by broth microdilution",
    disease_state = paste(
      "Neutropenic (cyclophosphamide) mice infected intramuscularly in each",
      "thigh with about 2.5 x 10^6 cfu, treated from 2 h after infection"
    ),
    dose_range = "16-512 mg/kg/day subcutaneously as q6h, q12h or q24h regimens",
    regions = "The Netherlands (Erasmus MC, Rotterdam)",
    notes = paste(
      "Thigh-infected neutropenic mice, two data points per regimen;",
      "bacterial burden counted 24 h after the first dose. E. coli ATCC",
      "25922 carries no ESBL (Table S1). The fAUC/MIC values were",
      "computed by the authors in MicLab 2.71 from the population PK",
      "model with an unbound fraction of 0.2; the typical PK of this",
      "file reproduces the fAUC/MIC positions of the Figure 3 data",
      "points to within 0.3%."
    )
  )

  ini({
    # Population PK (Table 1), typical values for a virtual 1 kg mouse at
    # the 10 mg/kg reference dose. The authors computed the PK/PD indices
    # from the typical PK, so no between-subject variability is carried.
    lka <- log(1.63)
    label("Absorption rate constant ka at a 10 mg/kg dose (1/h)") # Table 1 'ka (h-1)' = 1.63
    lcl <- log(0.384)
    label("Apparent clearance CL/F at a 10 mg/kg dose (L/h/kg)") # Table 1 'CL/F (L/h/kg)' = 0.384
    lvc <- log(0.774)
    label("Apparent central volume of distribution Vc/F (L/kg)") # Table 1 'Vc/F (L/kg)' = 0.774
    e_dose_ka <- -0.407
    label("Power exponent of dose per administration (DOSE/10 mg/kg) on ka (unitless)") # Table 1 'Dose on ka' = -0.407
    e_dose_cl <- -0.26
    label("Power exponent of dose per administration (DOSE/10 mg/kg) on CL/F (unitless)") # Table 1 'Dose on CL' = -0.26

    # Unbound fraction and isolate MIC, both measured inputs to the index.
    fu <- fixed(0.2)
    label("Unbound fraction of polymyxin B in mouse plasma (fraction)") # Results, 'Plasma protein binding': single factor of 20% used for fCmax and fAUC
    mic <- fixed(1)
    label("Broth-microdilution MIC of E. coli ATCC 25922 (mg/L)") # Table 2 and Table S1, MIC 1 mg/L

    # Inhibitory sigmoid Emax model of fAUC0-24/MIC (GraphPad Prism 8.0).
    # Figure 2 shows the fit only; the paper reports no targets for it.
    le0 <- log(2.012)
    label("24 h change in thigh burden without drug, upper plateau e0 (log10 cfu)") # Figure 2 panel 'E coli ATCC 25922', fAUC/MIC, upper plateau (digitised)
    lemax <- log(2.187)
    label("Span from upper to lower plateau, emax (log10 cfu)") # Figure 2 panel 'E coli ATCC 25922', fAUC/MIC, e0 - lower plateau -0.175 (digitised)
    lec50 <- log(48)
    label("fAUC0-24/MIC at half the maximum effect, ec50 (unitless)") # Figure 2 panel 'E coli ATCC 25922', fAUC/MIC (digitised)
    lhill <- log(3.898)
    label("Hill coefficient of the fAUC0-24/MIC effect curve (unitless)") # Figure 2 panel 'E coli ATCC 25922', fAUC/MIC (digitised)

    # Measured experimental input, not an estimate.
    log10_cfu0 <- fixed(7.18)
    label("Bacterial burden at the start of therapy (log10 cfu per thigh)") # Results, 'Dose fractionation study': average burden at t = 0 of 1.5 x 10^7 cfu/thigh; log10(1.5e7) = 7.18

    # The paper reports R^2 only, no residual SD.
    addSd <- fixed(0)
    label("Additive residual SD on log_cfu (log10 cfu per thigh)") # not reported
  })

  model({
    ka <- exp(lka) * (DOSE / 10)^e_dose_ka
    cl <- exp(lcl) * (DOSE / 10)^e_dose_cl
    vc <- exp(lvc)
    kel <- cl / vc

    e0 <- exp(le0)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    Cc <- central / vc

    # Unbound AUC over the first 24 h after the first dose; the state is
    # constant after 24 h, so the index holds its 24 h value from then on.
    d/dt(fauc_0_24) <- fu * Cc * (t < 24)
    fauc_mic <- fauc_0_24 / mic

    # Change in log10 cfu per thigh from the start of treatment to 24 h;
    # e0 at zero exposure, approaching e0 - emax at high exposure.
    dlog10cfu <- e0 - emax * fauc_mic^hill / (fauc_mic^hill + ec50^hill)
    log_cfu <- log10_cfu0 + dlog10cfu
    log_cfu ~ add(addSd)
  })
}
