Park_2026_tigecycline_invitro <- function() {
  description <- "In vitro (Mycobacterium abscessus ATCC 19977 in Middlebrook 7H9 broth with OADC, 24-well microscale dynamic time-kill assay at 37 C). Mechanism-based PK/PD model of tigecycline against M. abscessus. Tigecycline is chemically unstable in 7H9, and Park 2026 used that first-order degradation as the PK driver: each 2 mL well is a one-compartment system with volume fixed to the well volume and an elimination rate constant fixed to the degradation rate constant from the stability study. Bacteria (linear CFU/mL) grow logistically at the net rate K0 toward Nmax and are killed by an Emax term whose onset is delayed by (1 - exp(-beta * t)) and whose KC50 falls with time as KC50 * exp(-alpha * t), the increasing susceptibility seen in the long assays. Random effects on the inoculum N0, Nmax, alpha and beta describe the between-well and between-batch variation; the TGC concentration and the CFU count each carry a proportional residual error. Companion model Park_2026_tigecycline_invitro_media is the modified fit to the 96-well MIC-assay time courses in 7H9 and Mueller-Hinton II broth."
  reference <- paste(
    "Park H, Maloney Norcross SE, Hickey AJ, Gonzalez-Juarrero M, Meibohm B.",
    "In vitro PK/PD modeling and simulation to accurately assess the",
    "antimicrobial activity of tigecycline against Mycobacterium abscessus.",
    "Antimicrob Agents Chemother. 2026;70(2):e01025-25.",
    "doi:10.1128/aac.01025-25. PMID: 41432422. PMCID: PMC12888882.",
    "Structural model: Equations 3-5 and Figure 1; parameter estimates: Table 2."
  )
  vignette <- "Park_2026_tigecycline"
  units <- list(time = "h", dosing = "ug", concentration = "ug/mL")

  # `bact` is the single bacterial state of the logistic-growth time-kill
  # model; declared paper-specific as in Beredaki_2024_amphotericinB_liposomal_calbicans.R.
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

  population <- list(
    species = "in vitro (Mycobacterium abscessus ATCC 19977)",
    n_studies = 3L,
    disease_state = "Mycobacterium abscessus reference strain ATCC 19977 grown in freshly prepared Middlebrook 7H9 broth with OADC; tigecycline MIC 3.125 mg/L in 7H9 and 0.782 mg/L in Mueller-Hinton II broth by broth microdilution",
    dose_range = "Tigecycline spiked into 2 mL wells to initial concentrations of 0, 0.75, 1.5, 3, 6, 12, 24 and 48 mg/L; single dose, or the same dose repeated every 24 h for 5 days",
    regions = "University of Tennessee Health Science Center, Memphis, TN, USA",
    notes = paste(
      "Microscale dynamic time-kill assay in 24-well plates, triplicate wells, three independent assays with CFU counted in duplicate; inoculum about 1e6 CFU/mL of log-phase bacteria; incubation at 37 C; sampling for up to 336 h (single dose) or 624 h (multiple dose). Limit of counting 1e2 CFU/mL.",
      "Tigecycline concentrations in the multiple-dose regimens were measured by LC-MS/MS in parallel uninoculated wells; LLOQ 2.93 ng/mL.",
      "The TGC degradation half-life in 7H9 at 37 C was 15.8 h (stability study, Results); the PK/PD fit fixed Kdeg to 0.041 1/h (Table 2).",
      "Fitted in MonolixSuite 2022R1 by SAEM to CFU/mL data without log transformation."
    )
  )

  ini({
    # ---- PK of the well: fixed by experimental design (Table 2) ----------
    lvc <- fixed(log(2)); label("Log well volume V (mL)") # Table 2: V = 2.0 mL, footnote a 'Fixed based on volume per well in 24-well plate'
    lkdeg <- fixed(log(0.041)); label("Log first-order degradation rate constant of tigecycline in 7H9, Kdeg (1/h)") # Table 2: Kdeg = 0.041 1/h, footnote b 'Fixed based on the results of TGC stability study in 7H9'

    # ---- Bacterial growth (Table 2) ---------------------------------------
    linoc <- log(6.70e5); label("Log initial bacterial count N0 (CFU/mL)") # Table 2: N0 = 6.70 x 10^5 (RSE 10.8%)
    lkgrow <- log(0.12); label("Log net bacterial growth rate constant K0 (1/h)") # Table 2: K0 = 0.12 1/h (RSE 1.83%)
    lbmax <- log(1.40e9); label("Log maximum bacterial count Nmax (CFU/mL)") # Table 2: Nmax = 1.40 x 10^9 (RSE 6.15%)

    # ---- Drug effect (Table 2) --------------------------------------------
    lkmax <- log(0.20); label("Log maximum bacterial killing rate constant Kkill,max (1/h)") # Table 2: Kkill.max = 0.20 1/h (RSE 1.61%)
    lec50 <- log(0.65); label("Log tigecycline concentration for half-maximal killing at time zero, KC50 (ug/mL)") # Table 2: KC50 = 0.65 ug/mL (RSE 1.51%)
    lalpha <- log(0.017); label("Log rate constant of the time-dependent fall in KC50 (increasing susceptibility), alpha (1/h)") # Table 2: alpha = 0.017 1/h (RSE 6.57%)
    lbeta <- log(0.069); label("Log rate constant of the delayed onset of bacterial killing, beta (1/h)") # Table 2: beta = 0.069 1/h (RSE 13.8%)

    # ---- Between-well variability (Table 2, CV%) --------------------------
    # Monolix log-normal random effects reported as CV%; omega^2 = log(1 + CV^2).
    etalinoc ~ 0.1844 # Table 2: [N0] 45.0 CV% -> log(1 + 0.450^2)
    etalbmax ~ 0.01996 # Table 2: [Nmax] 14.2 CV% -> log(1 + 0.142^2)
    etalalpha ~ 0.04038 # Table 2: [alpha] 20.3 CV% -> log(1 + 0.203^2)
    etalbeta ~ 0.2874 # Table 2: [beta] 57.7 CV% -> log(1 + 0.577^2)

    # ---- Residual error (Table 2) -----------------------------------------
    propSd <- 0.19; label("Proportional residual error on the tigecycline concentration (fraction)") # Table 2: residual error (proportional) C = 0.19 (RSE 2.77%)
    propSd_cfu <- 0.47; label("Proportional residual error on the bacterial count (fraction)") # Table 2: residual error (proportional) N = 0.47 (RSE 5.43%)
  })

  model({
    # ---- Individual parameters --------------------------------------------
    vc <- exp(lvc)
    kdeg <- exp(lkdeg)
    n0 <- exp(linoc + etalinoc)
    kgrow <- exp(lkgrow)
    bmax <- exp(lbmax + etalbmax)
    kmax <- exp(lkmax)
    ec50 <- exp(lec50)
    alpha <- exp(lalpha + etalalpha)
    beta <- exp(lbeta + etalbeta)

    # ---- Tigecycline in the well: first-order chemical degradation (Eq 3) --
    d/dt(central) <- -kdeg * central
    Cc <- central / vc

    # ---- Bacteria: logistic growth minus delayed, time-sensitised Emax
    # killing (Eq 4). `t` is time since the start of the experiment.
    ekill <- (1 - exp(-beta * t)) * kmax * Cc / (ec50 * exp(-alpha * t) + Cc)
    d/dt(bact) <- (kgrow * (1 - bact / bmax) - ekill) * bact
    bact(0) <- n0

    # ---- Observations (Eq 5: proportional error on linear scale) -----------
    cfu <- bact
    Cc ~ prop(propSd)
    cfu ~ prop(propSd_cfu)
  })
}
