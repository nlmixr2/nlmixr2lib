Mao_2018_pravastatin_invitro_plasma <- function() {
  description <- paste(
    "In vitro (plated cryopreserved human hepatocytes in 100% human plasma).",
    "Michaelis-Menten model of the net active (rifampicin-sensitive, OATP-mediated)",
    "hepatic uptake of pravastatin into plated human hepatocytes incubated in",
    "modified human plasma medium. Net active uptake is the cell-associated",
    "pravastatin accumulated over 3 min alone minus that accumulated with 125 uM",
    "rifampicin, each corrected by its 0-min value, and is fitted against the",
    "NOMINAL total pravastatin concentration of the medium, so km is an APPARENT",
    "total-concentration constant that already carries the effect of plasma",
    "protein binding. The medium concentration is held constant over the 3-min",
    "initial-rate window, as in the fit. The fraction unbound in the incubation",
    "(set equal to the plasma fraction unbound, 0.485) is carried so the model can",
    "report the unbound Michaelis constant and the unbound uptake intrinsic",
    "clearance that the paper's Simcyp v15 full-PBPK model consumed as its OATP1B1",
    "input without any scaling factor. The Simcyp PBPK layer itself is not",
    "packaged: its distribution, absorption, renal and system parameters are",
    "platform-internal. Sibling model: Mao_2018_pravastatin_invitro_him, the same",
    "assay run in hepatocyte induction medium without plasma (Supplementary Figure",
    "S1).",
    sep = " "
  )
  reference <- paste(
    "Mao J, Doshi U, Wright M, Hop CECA, Li AP, Chen Y. Prediction of the",
    "Pharmacokinetics of Pravastatin as an OATP Substrate Using Plateable Human",
    "Hepatocytes With Human Plasma Data and PBPK Modeling. CPT Pharmacometrics",
    "Syst Pharmacol. 2018;7(4):251-258. doi:10.1002/psp4.12283. PMID: 29388346.",
    "PMCID: PMC5915609. Assay design: Methods, 'Active hepatic uptake transporter",
    "studies in plateable human hepatocytes with human plasma' and 'Modeling of",
    "the active hepatic uptake transporter kinetic parameters Jmax and Km of",
    "pravastatin'. Jmax and Km: Abstract; Results, 'Active hepatic uptake of",
    "pravastatin by OATP transporter'; Figure 1; Table 1. fu,inc = fu,p = 0.485:",
    "Table 1 and Methods, 'Clearance'.",
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
  # concentration placed into the well at time 0 (convention of
  # Hyland_2008_maraviroc_hlm.R and Morse_2012_ghb_rbc_invitro.R).
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
      "with pravastatin in modified human plasma medium (100% human plasma,",
      "HPZ-A) with or without 125 uM rifampicin"
    ),
    concentration_range = "0.14, 0.4, 1.2, 3.7, 11.1, 33.3, 100, 200, 300 and 400 uM pravastatin (nominal)",
    incubation = "uptake terminated at 0 and 3 min; 3 min chosen as within the linear range of active uptake",
    replication = "N = 3 (Figure 1)",
    transporter_abundance = "OATP1B1 18.3 and OATP1B3 1.7 fmol/million cells in the hepatocyte lot (Discussion)",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "Fit by the Michaelis-Menten model in GraphPad Prism 7 to the net uptake",
      "rate against nominal concentration. The paper assigns the whole active",
      "uptake to OATP1B1, the dominant isoform, because the hepatocyte assay",
      "cannot separate OATP1B1 from OATP1B3."
    )
  )

  ini({
    # Results 'Active hepatic uptake of pravastatin by OATP transporter':
    # "The Jmax and the apparent Km for pravastatin were determined as
    # 134.4 +/- 14.3 pmol/min/million cells and 76.77 +/- 27.4 uM". The
    # +/- values are the standard errors of the regression fit, not
    # between-donor variability (one donor lot was used).
    lvmax <- log(134.4); label("Log maximum rate of net active hepatic uptake, Jmax (pmol/min/million cells)") # Abstract; Results; Table 1 'Jmax, OATP' = 134.4
    lkm <- log(76.77); label("Log apparent Michaelis constant of net active uptake on nominal total concentration in 100% human plasma, Km (uM)") # Abstract; Results; Table 1 'Km, OATP (uM)' = 76.77

    # Table 1 'fu,inc 0.485 Same as fu,p'; Methods 'Clearance': "Because the
    # incubation medium is 100% human plasma, the fu,inc was assumed to be the
    # same as the fu,p 0.485". An assumption of the authors, not an estimate.
    fu_inc <- fixed(0.485); label("Fraction of pravastatin unbound in the incubation, set equal to fu,p (unitless)") # Table 1 'fu,inc' = 0.485 'Same as fu,p'

    # No residual error model is reported for the GraphPad fit; held at 0
    # rather than invented.
    addSd <- fixed(0); label("Additive residual SD on accumulated net uptake (pmol/million cells); not reported by the source") # Methods: GraphPad Prism 7 Michaelis-Menten fit; no residual model reported
  })

  model({
    vmax <- exp(lvmax)
    km <- exp(lkm)

    # Net active uptake rate (pmol/min/million cells) against the medium
    # concentration, the quantity plotted in Figure 1.
    uptake <- vmax * pravastatin / (km + pravastatin)

    # Unbound Michaelis constant (uM) and unbound linear-range uptake
    # intrinsic clearance (pmol/min/million cells per uM = uL/min/million
    # cells): Simcyp applies fu,inc to the entered Km, so these are the
    # values its OATP1B1 term effectively used (Discussion: "the apparent Km
    # was corrected with human plasma protein binding").
    km_u <- km * fu_inc
    clint_uptake <- vmax / km
    clint_uptake_u <- vmax / km_u

    # The fit used nominal concentrations, so the medium is not depleted.
    d/dt(pravastatin) <- 0
    d/dt(hep_pravastatin) <- uptake

    # The assay datum is cell-associated pravastatin accumulated over 3 min,
    # divided by 3 min to give the reported rate; since the medium is
    # constant, hep_pravastatin(t) = uptake * t exactly.
    hep_pravastatin ~ add(addSd)
  })
}
