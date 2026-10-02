Kosinsky_2022_camostat_invitro <- function() {
  description <- paste(
    "In vitro. Semi-mechanistic PD model of FOY-251 (the active camostat",
    "metabolite) reversible covalent inhibition of the host serine",
    "protease TMPRSS2, and the link to SARS-CoV-2 viral entry",
    "(Kosinsky 2022). Fit to two in-vitro experiments of Hoffmann et al.",
    "2020: cell-free recombinant TMPRSS2 activity after 1 h incubation",
    "with FOY-251, and SARS-2-S-driven pseudovirus entry after 2 h",
    "incubation. Active enzyme (target, relative activity 1 at baseline)",
    "is covalently bound by FOY-251 at rate kcat with half-saturation ki",
    "into an inactive complex that recovers at kdis (fixed to a 14 h",
    "complex half-life); there is no enzyme turnover in this cell-free",
    "assay (contrast the in-vivo Kosinsky_2022_camostat_pkpd, which adds",
    "kdeg synthesis/degradation). Remaining TMPRSS2 activity is linked to",
    "viral entry rate by an empirical Hill-Langmuir relationship. The",
    "applied FOY-251 well concentration is the per-record covariate",
    "CONC_FOY251_NM; ki, ki50 and the Hill coefficient are estimated,",
    "kcat and kdis are fixed. The fit implies ~95% TMPRSS2 inhibition is",
    "needed for 50% inhibition of viral entry rate."
  )
  reference <- paste(
    "Kosinsky Y, Peskov K, Stanski DR, Wetmore D, Vinetz J.",
    "Semi-Mechanistic Pharmacokinetic-Pharmacodynamic Model of Camostat",
    "Mesylate-Predicted Efficacy against SARS-CoV-2 in COVID-19.",
    "Microbiol Spectr. 2022;10(2):e02167-21.",
    "doi:10.1128/spectrum.02167-21.",
    "PD parameters in Table 2; the in-vitro ODEs are Equations 3",
    "(TMPRSS2 activity) and 4 (viral entry) of Materials and Methods;",
    "in-vitro data digitised from Hoffmann et al. 2020 (Cell 181:271-280,",
    "Figures 4 and 7).",
    sep = " "
  )
  vignette <- "Kosinsky_2022_camostat"
  units <- list(
    time = "h",
    dosing = "(no PK dosing; CONC_FOY251_NM is the applied in-vitro well concentration in nM)",
    concentration = "nM"
  )

  covariateData <- list(
    CONC_FOY251_NM = list(
      description = paste(
        "Applied FOY-251 concentration in the in-vitro assay well (nM).",
        "Per-record covariate held constant during the 1 h (TMPRSS2",
        "activity) or 2 h (viral entry) incubation; drives the covalent",
        "inhibition rate of TMPRSS2. Set to 0 for the drug-free control",
        "well."
      ),
      units = "nM",
      type = "continuous",
      reference_category = NULL,
      notes = "Applied well concentration from the Hoffmann et al. 2020 in-vitro experiments digitised in Kosinsky 2022 Figure 2A/B.",
      source_name = "C"
    )
  )

  compartmentData <- list(
    target = list(
      analyte = "TMPRSS2 (active recombinant serine protease)",
      units = "fraction of baseline activity",
      specimen = "not applicable",
      verified = TRUE
    ),
    complex = list(
      analyte = "FOY-251-TMPRSS2 covalent complex (inactive enzyme)",
      units = "fraction of baseline activity",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (recombinant TMPRSS2; SARS-2-S pseudovirus entry assay)",
    n_subjects = NA_integer_,
    disease_state = "n/a (cell-free and cell-based in-vitro assays)",
    notes = paste(
      "Two in-vitro experiments of Hoffmann et al. 2020 were combined for",
      "the fit: recombinant TMPRSS2 enzymatic activity versus FOY-251",
      "concentration after 1 h incubation (Kosinsky 2022 Figure 2A), and",
      "SARS-2-S-driven pseudovirus entry rate versus FOY-251 after 2 h",
      "incubation (Figure 2B). Only the kcat/ki ratio is identifiable from",
      "the data, so kcat was fixed at 400 1/h and ki estimated",
      "(Materials and Methods)."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # TMPRSS2 covalent inhibition (Table 2). ki and ki50 and hill are
    # estimated; kcat and kdis are fixed. In-vitro ODEs (Equation 3):
    #   dSP/dt  = -kcat*SP*C/(C+ki) + kdis*SPi
    #   dSPi/dt =  kcat*SP*C/(C+ki) - kdis*SPi,  SP(0)=1, SPi(0)=0.
    # Viral entry (Equation 4): 100% * SP^h / (SP^h + ki50^h).
    # ---------------------------------------------------------------------
    lki <- log(45638.51)
    label("FOY-251 TMPRSS2 half-saturation / inhibition constant ki (nM)") # Table 2 Ki = 45638.51 nM (= 45.6 uM; RSE 7.46%)
    kcat <- fixed(400)
    label("FOY-251 covalent-binding catalytic rate kcat (1/h)") # Table 2 kcat = 400 (fixed)
    kdis <- fixed(0.049)
    label("TMPRSS2 activity recovery rate kdis from covalent complex (1/h)") # Table 2 kdis = 0.049 (fixed; 14 h complex half-life)
    lki50 <- log(0.047)
    label("TMPRSS2 activity at half-maximal viral entry ki50 (fraction)") # Table 2 Ksp = 0.047 (RSE 42.4%)
    hill <- 0.59
    label("Hill coefficient of viral entry vs TMPRSS2 activity (unitless)") # Table 2 h = 0.59 (RSE 16.6%)

    # ---------------------------------------------------------------------
    # Additive (constant) residual error on each in-vitro readout, on the
    # fraction-of-baseline scale (Table 2 reports a1/a2 as percentages).
    # ---------------------------------------------------------------------
    addSd_tmprss2 <- 0.0346
    label("Additive residual error on relative TMPRSS2 activity (fraction)") # Table 2 a1 = 3.46% (RSE 22.5%)
    addSd_viralentry <- 0.0758
    label("Additive residual error on relative viral entry rate (fraction)") # Table 2 a2 = 7.58% (RSE 29.6%)
  })

  model({
    # PD constants.
    ki <- exp(lki)
    ksp <- exp(lki50)

    # Applied FOY-251 concentration in the assay well (nM), static during
    # incubation; 0 for the drug-free control.
    cfoy <- CONC_FOY251_NM

    # TMPRSS2 covalent inhibition without turnover (Equation 3). target and
    # complex are relative to baseline enzyme activity (target(0) = 1).
    bind <- kcat * target * cfoy / (cfoy + ki)
    d/dt(target) <- -bind + kdis * complex
    d/dt(complex) <- bind - kdis * complex
    target(0) <- 1
    complex(0) <- 0

    # Relative TMPRSS2 activity (read at 1 h) and relative SARS-2-S viral
    # entry rate (read at 2 h) by the empirical Hill-Langmuir link.
    tmprss2 <- target
    viralentry <- target^hill / (target^hill + ksp^hill)

    tmprss2 ~ add(addSd_tmprss2)
    viralentry ~ add(addSd_viralentry)
  })
}
