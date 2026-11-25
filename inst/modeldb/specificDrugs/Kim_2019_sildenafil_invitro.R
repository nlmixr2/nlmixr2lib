Kim_2019_sildenafil_invitro <- function() {
  description <- paste(
    "In vitro (USP paddle dissolution, pH 1.2). Michaelis-Menten dissolution",
    "model for sildenafil immediate-release (IR, 20 mg) and three",
    "sustained-release (SR, 60 mg: SRfast, SRmedium, SRslow) HPMC matrix",
    "tablets (Kim 2019, Equation 1, Table 3). Undissolved drug in the tablet",
    "(depot) dissolves into the medium (released) at a dose-normalised",
    "Michaelis-Menten rate with a formulation-specific maximum rate and a",
    "shared half-saturation amount; the percent dissolved is",
    "100 * released / dose. Fitted in S-ADAPT with between-tablet variability.",
    "The fitted maximum rates are the inputs of the in vitro-in vivo",
    "correlation in Kim_2019_sildenafil_dog_ivivc. Release-only: no",
    "disposition; residual error is not reported."
  )
  reference <- paste(
    "Kim TH, Shin S, Jeong SW, Lee JB, Shin BS.",
    "Physiologically Relevant In Vitro-In Vivo Correlation (IVIVC) Approach",
    "for Sildenafil with Site-Dependent Dissolution.",
    "Pharmaceutics. 2019;11(6):251. doi:10.3390/pharmaceutics11060251"
  )
  vignette <- "Kim_2019_sildenafil"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "not applicable (no disposition); dissolved amount in mg"
  )

  covariateData <- list(
    FORM_SILDENAFIL_SRFAST = list(
      description = "1 = SRfast sustained-release tablet (10% HPMC 2208), 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (IR tablet when all three FORM_SILDENAFIL_SR* indicators are 0)",
      notes = "Selects the SRfast in vitro maximum dissolution rate (Table 3) and its between-tablet variability.",
      source_name = "SRfast"
    ),
    FORM_SILDENAFIL_SRMEDIUM = list(
      description = "1 = SRmedium sustained-release tablet (30% HPMC 2208), 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (IR tablet when all three FORM_SILDENAFIL_SR* indicators are 0)",
      notes = "Selects the SRmedium in vitro maximum dissolution rate (Table 3) and its between-tablet variability.",
      source_name = "SRmedium"
    ),
    FORM_SILDENAFIL_SRSLOW = list(
      description = "1 = SRslow sustained-release tablet (60% HPMC 2208), 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (IR tablet when all three FORM_SILDENAFIL_SR* indicators are 0)",
      notes = "Selects the SRslow in vitro maximum dissolution rate (Table 3) and its between-tablet variability.",
      source_name = "SRslow"
    )
  )

  dosing <- "depot"

  # Bookkeeping sink integrating the dissolution flux: the dissolution
  # medium, not a biological compartment.
  paper_specific_compartments <- c("released")

  compartmentData <- list(
    depot = list(
      analyte = "sildenafil (undissolved, in the tablet)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    released = list(
      analyte = "sildenafil (dissolved in the medium)",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (sildenafil citrate tablets, USP paddle apparatus, 0.1 N HCl pH 1.2)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    disease_state = "Not applicable (in vitro dissolution)",
    dose_range = "IR tablet 20 mg; SRfast, SRmedium and SRslow tablets 60 mg sildenafil (Table 1)",
    notes = paste(
      "Methods 2.3: Distek 2500 paddle apparatus, 100 rpm, 37 C; samples at",
      "0.25-24 h. Only the pH 1.2 profiles were modelled (Methods 2.6.2,",
      "Results 3.3.1); dissolution at pH 6.8 was incomplete (Figure 2) and is",
      "not described by this model. The number of tablets per formulation is",
      "not reported. Tablet compositions (w/w): HPMC 2208-100 cps 0%, 10%,",
      "30% and 60% for IR, SRfast, SRmedium and SRslow (Table 1)."
    )
  )

  ini({
    # Equation (1), dose-normalised (Methods 2.6.2 and 2.6.3): AM50/dose is a
    # fraction of the dose and Vmax/dose is in 1/h.
    lam50 <- log(0.395); label("Dose-normalised amount in the tablet at half-maximal in vitro dissolution rate AM50,in vitro/dose (fraction)") # Table 3 AM50,in vitro/dose = 0.395
    lvmax <- log(6.58); label("Dose-normalised in vitro maximum dissolution rate, IR tablet (1/h)") # Table 3 Vmax,in vitro,IR/dose = 6.58 1/h
    lvmax_srfast <- log(1.40); label("Dose-normalised in vitro maximum dissolution rate, SRfast tablet (1/h)") # Table 3 Vmax,in vitro,SRfast/dose = 1.40 1/h
    lvmax_srmedium <- log(0.358); label("Dose-normalised in vitro maximum dissolution rate, SRmedium tablet (1/h)") # Table 3 Vmax,in vitro,SRmedium/dose = 0.358 1/h
    lvmax_srslow <- log(0.124); label("Dose-normalised in vitro maximum dissolution rate, SRslow tablet (1/h)") # Table 3 Vmax,in vitro,SRslow/dose = 0.124 1/h

    # Between-tablet variability, Table 3 'BSV' column: SD of the log-scale
    # random effect (S-ADAPT), so omega^2 = BSV^2.
    etalam50 ~ 0.00099856 # Table 3 AM50,in vitro BSV 0.0316 -> 0.0316^2
    etalvmax ~ 0.00080089 # Table 3 Vmax,in vitro,IR BSV 0.0283 -> 0.0283^2
    etalvmax_srfast ~ 0.021609 # Table 3 Vmax,in vitro,SRfast BSV 0.147 -> 0.147^2
    etalvmax_srmedium ~ 0.00061009 # Table 3 Vmax,in vitro,SRmedium BSV 0.0247 -> 0.0247^2
    etalvmax_srslow ~ 0.00110889 # Table 3 Vmax,in vitro,SRslow BSV 0.0333 -> 0.0333^2
  })

  model({
    am50 <- exp(lam50 + etalam50)
    vmax_ir <- exp(lvmax + etalvmax)
    vmax_srfast <- exp(lvmax_srfast + etalvmax_srfast)
    vmax_srmedium <- exp(lvmax_srmedium + etalvmax_srmedium)
    vmax_srslow <- exp(lvmax_srslow + etalvmax_srslow)
    form_ir <- 1 - FORM_SILDENAFIL_SRFAST - FORM_SILDENAFIL_SRMEDIUM - FORM_SILDENAFIL_SRSLOW
    vmax <- form_ir * vmax_ir +
      FORM_SILDENAFIL_SRFAST * vmax_srfast +
      FORM_SILDENAFIL_SRMEDIUM * vmax_srmedium +
      FORM_SILDENAFIL_SRSLOW * vmax_srslow

    # Equation (1) in amounts: d(X/Dose)/dt = -Vmax/Dose * (X/Dose) / (AM50/Dose + X/Dose)
    dose_tab <- podo(depot)
    rdiss <- 0
    if (dose_tab > 0) rdiss <- vmax * depot / (am50 + depot / dose_tab)

    d/dt(depot) <- -rdiss
    d/dt(released) <- rdiss
  })
}
