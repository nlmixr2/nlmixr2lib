Sluka_2016_acetaminophen_hepatocyte_qsp <- function() {
  description <- paste(
    "QSP (subcellular reaction kinetics, SBML). The standalone per-hepatocyte",
    "acetaminophen (APAP) metabolism module of the Sluka 2016 liver-centric",
    "multiscale framework: one well-stirred hepatocyte in which APAP is",
    "glucuronidated, sulfated (both Michaelis-Menten) and oxidised by CYP2E1",
    "(Michaelis-Menten) to NAPQI, which conjugates with glutathione (GSH) by",
    "second-order mass action; GSH is resynthesised towards a ceiling. Time",
    "is in seconds and concentrations in mmol/L, as in the source. In the",
    "multiscale model a copy of this network runs in each simulated",
    "hepatocyte and exchanges APAP and its conjugates with the sinusoid",
    "blood; standalone (paper Fig 6) it starts from an intracellular APAP",
    "concentration and has no inflow or outflow. Deterministic, no",
    "variability or residual error. The whole-body module is",
    "Sluka_2016_acetaminophen_pbpk.",
    sep = " "
  )
  reference <- paste(
    "Sluka JP, Fu X, Swat M, Belmonte JM, Cosmanescu A, Clendenon SG,",
    "Wambaugh JF, Glazier JA. A Liver-Centric Multiscale Modeling Framework",
    "for Xenobiotics. PLoS One. 2016;11(9):e0162428.",
    "doi:10.1371/journal.pone.0162428. PMID 27636091. PMCID PMC5026379.",
    "Rate laws from Table 5 and the deposited SBML file",
    "APAP_GSH_Liver_Metab_Sul_Glu.xml (S2 File, BioModels MODEL1509230011);",
    "parameter values from Table 6 (parameter set REFSIM; S1 File Table 10).",
    sep = " "
  )
  vignette <- "Sluka_2016_acetaminophen"
  units <- list(time = "s", dosing = "mmol", concentration = "mmol/L")

  # The intracellular species of one hepatocyte are not organ-level PK
  # compartments and have no canonical names; the paper's SBML species
  # (APAP, APAPconj_Glu, APAPconj_Sul, NAPQI, GSH, NAPQIGSH) map onto
  # these hepatocyte-prefixed states one-to-one.
  paper_specific_compartments <- c("hep_apap", "hep_apapg", "hep_apaps", "hep_napqi", "hep_gsh", "hep_napqigsh")

  # APAP enters the standalone network as an initial intracellular amount
  # (Fig 6: 0.1 mmol/L in the 1 L SBML compartment).
  dosing <- "hep_apap"

  compartmentData <- list(
    hep_apap = list(analyte = "acetaminophen", units = "mmol", specimen = "tissue", verified = TRUE),
    hep_apapg = list(analyte = "acetaminophen glucuronide", units = "mmol", specimen = "tissue", verified = TRUE),
    hep_apaps = list(analyte = "acetaminophen sulfate", units = "mmol", specimen = "tissue", verified = TRUE),
    hep_napqi = list(
      analyte = "N-acetyl-p-benzoquinone imine (NAPQI)",
      units = "mmol",
      specimen = "tissue",
      verified = TRUE
    ),
    hep_gsh = list(analyte = "glutathione (reduced)", units = "mmol", specimen = "tissue", verified = TRUE),
    hep_napqigsh = list(analyte = "NAPQI-glutathione conjugate", units = "mmol", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human (single hepatocyte; in silico)",
    n_subjects = 0L,
    n_studies = 0L,
    disease_state = "Not applicable: a single simulated hepatocyte, not a fit to individual data.",
    dose_range = "Standalone run (Fig 6): initial intracellular APAP 0.1 mmol/L (15 ug/mL).",
    notes = paste(
      "The Michaelis-Menten constants were taken or adjusted from Critchley et al.",
      "2005 (paper reference 58) so that the complete multiscale model reproduces",
      "the in vivo ratio of APAP-glucuronide to APAP-sulfate after a 1.4 g oral",
      "dose (Results, 'Subcellular Reaction Module'); GSH constants are from",
      "paper reference 59 or estimated (Table 6). No subject-level data were",
      "fitted to this module on its own.",
      sep = " "
    )
  )

  ini({
    v_hep <- fixed(1); label("Volume of the hepatocyte reaction compartment (L)") # deposited SBML compartment 'Hepatocyte' size = 1 litre

    lvmax_gluc <- fixed(log(0.001)); label("Maximum rate of APAP glucuronidation (mmol/L/s)") # Table 6 'sc_Vmax_GLUC' 0.001 mmol/L/s; S1 File Table 10 REFSIM 1.E-3
    lkm_gluc <- fixed(log(1)); label("Michaelis constant of APAP glucuronidation (mmol/L)") # Table 6 'sc_Km_GLUC' 1 mmol/L
    lvmax_sulf <- fixed(log(0.000175)); label("Maximum rate of APAP sulfation (mmol/L/s)") # Table 6 'sc_Vmax_SULF' 0.000175 mmol/L/s; S1 File Table 10 REFSIM 1.75E-4
    lkm_sulf <- fixed(log(0.2)); label("Michaelis constant of APAP sulfation (mmol/L)") # Table 6 'sc_Km_SULF' 0.2 mmol/L
    lvmax_cyp2e1 <- fixed(log(0.00002)); label("Maximum rate of CYP2E1 oxidation of APAP to NAPQI (mmol/L/s)") # Table 6 'sc_Vmax_2E1_APAP' 0.00002 mmol/L/s; S1 File Table 10 REFSIM 2.E-5
    lkm_cyp2e1 <- fixed(log(1.29)); label("Michaelis constant of CYP2E1 oxidation of APAP (mmol/L)") # Table 6 'sc_Km_2E1_APAP' 1.29 mmol/L
    k_napqi_gsh <- fixed(0.1); label("Second-order rate constant of NAPQI-GSH conjugation (L/mmol/s)") # Table 6 'sc_kNAPQIGSH' 0.1 (printed unit mmol/L/s); SBML kNapqiGsh units second_order_rate_constant
    ksyn_gsh <- fixed(0.0001); label("First-order rate constant of GSH resynthesis towards its ceiling (1/s)") # Table 6 'sc_kGsh' 0.0001 (printed unit mmol/L/s); SBML kGsh units first_order_rate_constant
    gsh0 <- fixed(9.9); label("Initial cellular GSH concentration (mmol/L)") # Table 6 'GSH0' 9.9 mmol/L
    gshmax <- fixed(10); label("Maximum cellular GSH concentration (mmol/L)") # Table 6 'GSHmax' 10 mmol/L
  })

  model({
    vmax_gluc <- exp(lvmax_gluc)
    km_gluc <- exp(lkm_gluc)
    vmax_sulf <- exp(lvmax_sulf)
    km_sulf <- exp(lkm_sulf)
    vmax_cyp2e1 <- exp(lvmax_cyp2e1)
    km_cyp2e1 <- exp(lkm_cyp2e1)

    c_apap <- hep_apap / v_hep
    c_napqi <- hep_napqi / v_hep
    c_gsh <- hep_gsh / v_hep

    # Table 5 rate equations (mmol/L/s).
    r_gluc <- vmax_gluc * c_apap / (km_gluc + c_apap) # APAP -> APAPG
    r_sulf <- vmax_sulf * c_apap / (km_sulf + c_apap) # APAP -> APAPS
    r_cyp2e1 <- vmax_cyp2e1 * c_apap / (km_cyp2e1 + c_apap) # APAP -> NAPQI
    r_conj <- k_napqi_gsh * c_napqi * c_gsh # NAPQI + GSH -> NAPQI-GSH
    r_gsh <- ksyn_gsh * (gshmax - c_gsh) # phi -> GSH

    d/dt(hep_apap) <- -v_hep * (r_gluc + r_sulf + r_cyp2e1)
    d/dt(hep_apapg) <- v_hep * r_gluc
    d/dt(hep_apaps) <- v_hep * r_sulf
    d/dt(hep_napqi) <- v_hep * (r_cyp2e1 - r_conj)
    d/dt(hep_gsh) <- v_hep * (r_gsh - r_conj)
    d/dt(hep_napqigsh) <- v_hep * r_conj
    hep_gsh(0) <- gsh0 * v_hep

    # Intracellular concentrations (mmol/L), the quantities plotted in Fig 6.
    Cc <- c_apap
    Cc_apapg <- hep_apapg / v_hep
    Cc_apaps <- hep_apaps / v_hep
    Cc_napqi <- c_napqi
    Cc_gsh <- c_gsh
    Cc_napqigsh <- hep_napqigsh / v_hep
  })
}
