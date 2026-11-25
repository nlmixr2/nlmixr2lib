Mao_2018_pravastatin_invitro_him <- function() {
  description <- paste(
    "In vitro (plated cryopreserved human hepatocytes in protein-free hepatocyte",
    "induction medium). Michaelis-Menten model of the net active",
    "(rifampicin-sensitive, OATP-mediated) hepatic uptake of pravastatin into",
    "plated human hepatocytes incubated in hepatocyte induction medium (HIM)",
    "instead of human plasma. Net active uptake is the cell-associated pravastatin",
    "accumulated over 3 min alone minus that accumulated with 25 uM rifampicin,",
    "each corrected by its 0-min value, and is fitted against the nominal",
    "pravastatin concentration of the medium; the medium concentration is held",
    "constant over the 3-min initial-rate window, as in the fit. The fraction",
    "unbound in the incubation (hepatocyte binding) is carried so the model can",
    "report the unbound Michaelis constant and the unbound uptake intrinsic",
    "clearance that the paper entered in place of the plasma-medium values to show",
    "that the Simcyp v15 full-PBPK model still captured the intravenous and oral",
    "pravastatin profiles without a scaling factor (Supplementary Figures S2 and",
    "S3). The Simcyp PBPK layer itself is not packaged. Sibling model:",
    "Mao_2018_pravastatin_invitro_plasma, the paper's primary assay in 100% human",
    "plasma (Figure 1).",
    sep = " "
  )
  reference <- paste(
    "Mao J, Doshi U, Wright M, Hop CECA, Li AP, Chen Y. Prediction of the",
    "Pharmacokinetics of Pravastatin as an OATP Substrate Using Plateable Human",
    "Hepatocytes With Human Plasma Data and PBPK Modeling. CPT Pharmacometrics",
    "Syst Pharmacol. 2018;7(4):251-258. doi:10.1002/psp4.12283. PMID: 29388346.",
    "PMCID: PMC5915609. Assay design: Supplementary Methods, 'Active Hepatic",
    "Uptake Transporter Studies in Plateable Human Hepatocytes with Human Plasma'",
    "(which describes the HIM incubation). Jmax and Km: Discussion, paragraph 2;",
    "Supplementary Figure S1. fu,inc: Supplementary Figure S2 and S3 captions",
    "(0.98, the simulation input); the Discussion gives 0.96 citing Menochet K et",
    "al. Drug Metab Dispos 2012;40:1744-1756.",
    sep = " "
  )
  vignette <- "Mao_2018_pravastatin_invitro"
  units <- list(
    time = "min",
    dosing = "uM (nominal medium concentration)",
    concentration = "uM",
    uptake = "pmol/million cells"
  )

  # The incubation medium is dosed directly: there is no depot and no
  # central compartment, so the "dose" is the nominal pravastatin
  # concentration placed into the well at time 0.
  dosing <- c("pravastatin")

  # `pravastatin` holds the medium concentration and `hep_pravastatin` the
  # cell-associated net active uptake per million hepatocytes; neither is a
  # canonical compartment role.
  paper_specific_compartments <- c("pravastatin", "hep_pravastatin")

  compartmentData <- list(
    pravastatin = list(
      analyte = "pravastatin",
      units = "uM",
      specimen = "administration site",
      verified = TRUE
    ),
    hep_pravastatin = list(
      analyte = "pravastatin (net rifampicin-sensitive uptake into hepatocytes)",
      units = "pmol/million cells",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (plated cryopreserved human hepatocytes, single donor)",
    n_subjects = 1L,
    n_studies = 1L,
    donor = "lot HH1057, 33-year-old Caucasian female",
    system = paste(
      "Cryopreserved plateable human hepatocytes seeded at 350,000 cells/well on",
      "collagen-coated 24-well plates, attached for 4 h at 37 C, then incubated",
      "with pravastatin in hepatocyte induction medium (HIM) with or without",
      "25 uM rifampicin"
    ),
    concentration_range = "0.14, 0.4, 1.2, 3.7, 11.1, 33.3, 100, 200, 300 and 400 uM pravastatin (nominal)",
    incubation = "uptake terminated at 0 and 3 min",
    replication = "N = 3 (Supplementary Figure S1)",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "Same donor lot and plate format as the plasma-medium assay; the medium",
      "and the rifampicin concentration (25 vs 125 uM) are the differences",
      "stated by the Supplementary Methods."
    )
  )

  ini({
    # Discussion: "when the in vitro data (Jmax 28.49 +/- 3.11 pmol/min/
    # million cells and Km 11.72 +/- 6.83 uM) corrected for hepatocyte
    # binding of pravastatin (fu,inc 0.96) was incorporated into the same
    # model". The +/- values are regression standard errors; Supplementary
    # Figures S2/S3 captions repeat 28.49 +/- 3.11 and 11.72 +/- 6.83.
    lvmax <- log(28.49); label("Log maximum rate of net active hepatic uptake, Jmax (pmol/min/million cells)") # Discussion para 2; Supplementary Figure S2/S3 captions 'Jmax, OATP as 28.49'
    lkm <- log(11.72); label("Log apparent Michaelis constant of net active uptake on nominal concentration in induction medium, Km (uM)") # Discussion para 2; Supplementary Figure S2/S3 captions 'Km, OATP as 11.72'

    # The Supplementary Figure S2 and S3 captions (the simulation inputs) both
    # give fu,inc = 0.98; the Discussion gives 0.96. The as-run 0.98 is
    # carried; the vignette Assumptions records the conflict.
    fu_inc <- fixed(0.98); label("Fraction of pravastatin unbound in the incubation, hepatocyte binding (unitless)") # Supplementary Figure S2/S3 captions 'fu,inc as 0.98'; Discussion para 2 prints 0.96

    # No residual error model is reported for the GraphPad fit; held at 0
    # rather than invented.
    addSd <- fixed(0); label("Additive residual SD on accumulated net uptake (pmol/million cells); not reported by the source") # Methods: GraphPad Prism 7 Michaelis-Menten fit; no residual model reported
  })

  model({
    vmax <- exp(lvmax)
    km <- exp(lkm)

    # Net active uptake rate (pmol/min/million cells) against the medium
    # concentration, the quantity plotted in Supplementary Figure S1.
    uptake <- vmax * pravastatin / (km + pravastatin)

    # Unbound Michaelis constant (uM) and unbound linear-range uptake
    # intrinsic clearance (uL/min/million cells), the values the Simcyp
    # OATP1B1 term effectively used after the fu,inc correction.
    km_u <- km * fu_inc
    clint_uptake <- vmax / km
    clint_uptake_u <- vmax / km_u

    # The fit used nominal concentrations, so the medium is not depleted.
    d/dt(pravastatin) <- 0
    d/dt(hep_pravastatin) <- uptake

    # Accumulated net uptake; hep_pravastatin(t) = uptake * t exactly.
    hep_pravastatin ~ add(addSd)
  })
}
