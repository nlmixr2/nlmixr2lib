Park_2026_tigecycline_invitro_media <- function() {
  description <- "In vitro (Mycobacterium abscessus ATCC 19977 in Middlebrook 7H9 broth with OADC or in Mueller-Hinton II broth, 96-well broth-microdilution MIC plates at 37 C). Modified mechanism-based PK/PD model of tigecycline against M. abscessus fitted jointly to the 72 h bacterial time courses of the MIC assay run in the two media. The structure is that of Park_2026_tigecycline_invitro (one-compartment well with first-order chemical degradation of tigecycline; logistic bacterial growth at net rate K0 toward Nmax; Emax killing with onset delayed by (1 - exp(-beta * t))), but the alpha term for the time-dependent fall in KC50 is omitted because the 72 h observation window could not estimate it. The medium selects the fixed degradation rate constant (7H9 or MHB II stability study) and scales the growth rate constant by exp(CovK0) in MHB II; random effect on K0 only."
  reference <- paste(
    "Park H, Maloney Norcross SE, Hickey AJ, Gonzalez-Juarrero M, Meibohm B.",
    "In vitro PK/PD modeling and simulation to accurately assess the",
    "antimicrobial activity of tigecycline against Mycobacterium abscessus.",
    "Antimicrob Agents Chemother. 2026;70(2):e01025-25.",
    "doi:10.1128/aac.01025-25. PMID: 41432422. PMCID: PMC12888882.",
    "Structural model: Equations 3-5 and Materials and Methods 'Model application to MIC assay platform'; parameter estimates: Table 3."
  )
  vignette <- "Park_2026_tigecycline"
  units <- list(time = "h", dosing = "ug", concentration = "ug/mL")

  paper_specific_compartments <- c("bact")

  compartmentData <- list(
    central = list(
      analyte = "tigecycline",
      units = "ug",
      specimen = "not applicable",
      verified = TRUE
    ),
    bact = list(
      analyte = "Mycobacterium abscessus ATCC 19977, total viable count",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    MEDIUM_MHB = list(
      description = "Culture-medium indicator of the in vitro experiment: 1 = Mueller-Hinton II broth (MHB II), 0 = Middlebrook 7H9 broth with OADC supplement.",
      units = "(binary)",
      type = "binary",
      source_name = "media",
      reference_category = "0 (Middlebrook 7H9 with OADC)",
      notes = paste(
        "Selects the fixed tigecycline degradation rate constant (Kdeg.7H9 = 0.044 1/h, Kdeg.MHB II = 0.017 1/h; Table 3) and multiplies the bacterial net growth rate constant by exp(CovK0) in MHB II (Table 3 footnote d).",
        "7H9 is the reference because it is the medium of the main dynamic time-kill experiments and carries the typical K0."
      )
    )
  )

  population <- list(
    species = "in vitro (Mycobacterium abscessus ATCC 19977)",
    n_studies = 1L,
    disease_state = "Mycobacterium abscessus reference strain ATCC 19977 in freshly prepared Middlebrook 7H9 broth with OADC or Mueller-Hinton II broth; tigecycline MIC 3.125 mg/L in 7H9 and 0.782 mg/L in MHB II",
    dose_range = "Tigecycline serially diluted from 100 to 0.0976 mg/L final concentration in 200 uL wells of a 96-well broth-microdilution plate (Figure 6 shows 0, 0.098, 0.391, 1.56 and 6.25 mg/L); single addition at time zero",
    regions = "University of Tennessee Health Science Center, Memphis, TN, USA",
    notes = paste(
      "Second plate of the broth-microdilution MIC test, run in duplicate per medium with an inoculum of 1e6 CFU/mL (mixed 1:1 with drug-containing medium); CFU/mL counted at 0, 24, 48 and 72 h.",
      "Degradation half-lives at 37 C: 15.8 h in 7H9 and 40.8 h in MHB II (stability study, Results).",
      "Fitted in MonolixSuite 2022R1 by SAEM to CFU/mL data without log transformation."
    )
  )

  ini({
    # ---- PK of the well: fixed by experimental design (Table 3) ----------
    lvc <- fixed(log(0.2)); label("Log well volume V (mL)") # Table 3: V = 0.2 mL, footnote a 'Fixed on volume per well in 96-well plate'
    lkdeg_7h9 <- fixed(log(0.044)); label("Log first-order degradation rate constant of tigecycline in 7H9, Kdeg.7H9 (1/h)") # Table 3: Kdeg.7H9 = 0.044 1/h, footnote b (stability study)
    lkdeg_mhb <- fixed(log(0.017)); label("Log first-order degradation rate constant of tigecycline in MHB II, Kdeg.MHB II (1/h)") # Table 3: Kdeg.MHB II = 0.017 1/h, footnote b (stability study)

    # ---- Bacterial growth (Table 3) ---------------------------------------
    linoc <- log(5.11e5); label("Log initial bacterial count N0 (CFU/mL)") # Table 3: N0 = 5.11 x 10^5 (RSE 11.3%)
    lkgrow <- log(0.19); label("Log net bacterial growth rate constant K0 in 7H9 (1/h)") # Table 3: K0 = 0.19 1/h (RSE 11.3%)
    e_mhb_kgrow <- -0.73; label("Effect of MHB II medium on log K0, CovK0 (unitless)") # Table 3: CovK0 = -0.73 (RSE 18.7%); footnote d: K0,MHB = K0 x exp(CovK0)
    lbmax <- log(2.68e9); label("Log maximum bacterial count Nmax (CFU/mL)") # Table 3: Nmax = 2.68 x 10^9 (RSE 20%)

    # ---- Drug effect (Table 3) --------------------------------------------
    lkmax <- log(0.18); label("Log maximum bacterial killing rate constant Kkill,max (1/h)") # Table 3: Kkill.max = 0.18 1/h (RSE 10.5%)
    lec50 <- log(0.11); label("Log tigecycline concentration for half-maximal killing, KC50 (ug/mL)") # Table 3: KC50 = 0.11 ug/mL (RSE 41.8%)
    lbeta <- log(0.16); label("Log rate constant of the delayed onset of bacterial killing, beta (1/h)") # Table 3: beta = 0.16 1/h (RSE 32.3%); footnote e: alpha omitted

    # ---- Between-well variability (Table 3) -------------------------------
    etalkgrow ~ 0.1849 # Table 3: omega [K0] = 0.43, Monolix SD of the log-normal random effect -> 0.43^2

    # ---- Residual error (Table 3) -----------------------------------------
    propSd_cfu <- 0.47; label("Proportional residual error on the bacterial count (fraction)") # Table 3: residual error (proportional) N = 0.47 (RSE 5.43%)
  })

  model({
    # ---- Individual parameters --------------------------------------------
    vc <- exp(lvc)
    kdeg <- exp(lkdeg_7h9 * (1 - MEDIUM_MHB) + lkdeg_mhb * MEDIUM_MHB)
    n0 <- exp(linoc)
    kgrow <- exp(lkgrow + e_mhb_kgrow * MEDIUM_MHB + etalkgrow)
    bmax <- exp(lbmax)
    kmax <- exp(lkmax)
    ec50 <- exp(lec50)
    beta <- exp(lbeta)

    # ---- Tigecycline in the well: first-order chemical degradation (Eq 3) --
    d/dt(central) <- -kdeg * central
    Cc <- central / vc

    # ---- Bacteria: Eq 4 with alpha omitted (Table 3 footnote e), so KC50
    # is constant in time. `t` is time since drug addition.
    ekill <- (1 - exp(-beta * t)) * kmax * Cc / (ec50 + Cc)
    d/dt(bact) <- (kgrow * (1 - bact / bmax) - ekill) * bact
    bact(0) <- n0

    # ---- Observation (Eq 5: proportional error on linear scale) ------------
    cfu <- bact
    cfu ~ prop(propSd_cfu)
  })
}
