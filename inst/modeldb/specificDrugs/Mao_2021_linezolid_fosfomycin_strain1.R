Mao_2021_linezolid_fosfomycin_strain1 <- function() {
  description <- "In vitro (Enterococcus faecalis clinical isolate No. 1; linezolid MIC 4 mg/L, fosfomycin MIC 64 mg/L). Semimechanistic PK/PD model of linezolid and fosfomycin alone and in combination in a one-compartment dynamic in vitro PK/PD system. Logistic growth of the log10 viable count toward a stationary-phase maximum, killed by a sigmoid Emax term per drug whose EC50 is inflated over time by a concentration-driven adaptive-resistance factor; the two drug effects combine through an empirical interaction exponent. Linezolid is held at a constant concentration; fosfomycin follows a human two-compartment profile (0.5-h infusion every 8 h). CAUTION: the monotherapy arms reproduce the published fits, but the combination arms (the paper's headline) do not reproduce the published Figure 5-7 curves with the printed parameters under any reading of Equation 5; see the vignette. One of four strain-specific files from this paper"
  reference <- "Mao J, Li T, Zhang N, Wang S, Li Y, Peng Y, Liu H, Yang G, Yan Y, Jiang L, Liu Y, Li J, Huang X. Dose optimization of combined linezolid and fosfomycin against Enterococcus by using an in vitro pharmacokinetic/pharmacodynamic model. Microbiol Spectr. 2021;9(3):e00871-21. doi:10.1128/spectrum.00871-21. PMCID: PMC8635129. Structural model: Materials and Methods, 'Semimechanical PK/PD model', Equations 1-5. Pharmacodynamic estimates: Table 3, column No. 1. Fosfomycin pharmacokinetics: supplemental Table S1. Interaction (Equation 5) encoded in the form of the cited source, Mohamed AF et al. J Antimicrob Chemother 2016;71(5):1279-1290 (reference 17 of the paper), because the form printed in the paper is inconsistent with its own monotherapy fits"
  vignette <- "Mao_2021_linezolid_fosfomycin"
  units <- list(
    time = "h",
    dosing = "mg (fosfomycin, into central); mg/L (linezolid, into lzd)",
    concentration = "mg/L"
  )

  # Fosfomycin is infused into `central` (mg; the volume in L gives mg/L).
  # Linezolid was held constant in the system (added to every reservoir), so it
  # is a bath-concentration state dosed once at time 0 with amt = the target
  # concentration in mg/L.
  dosing <- c("central", "lzd")

  compartmentData <- list(
    central = list(analyte = "fosfomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "fosfomycin", units = "mg", specimen = "plasma", verified = TRUE),
    lzd = list(analyte = "linezolid", units = "mg/L", specimen = "administration site", verified = TRUE),
    log_cfu = list(
      analyte = "Enterococcus faecalis clinical isolate No. 1, total viable count",
      units = "log10 CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "in vitro (Enterococcus faecalis clinical isolate No. 1)",
    n_subjects = 1L,
    n_studies = 1L,
    organism = "Enterococcus faecalis clinical isolate No. 1, isolated from urine or blood at the First Affiliated Hospital of Anhui Medical University",
    mic_values = c(linezolid = "4 mg/L", fosfomycin = "64 mg/L"),
    system = "Dynamic one-compartment in vitro PK/PD model (200 mL central chamber sealed with a 0.45 um filter, peristaltic-pump dilution) in Mueller-Hinton broth with 25 mg/L glucose-6-phosphate; viable counts at 0, 2, 4, 6, 8, 10, 12 and 24 h",
    temperature = "37 C",
    duration = "24 h",
    starting_inoculum = "~1 x 10^6 CFU/mL (not estimated)",
    limit_of_quantification = "limit of detection 100 CFU/mL (2 log10 CFU/mL)",
    dose_range = "linezolid 4 or 8 mg/L constant; fosfomycin 4 or 6 g as a 0.5-h infusion every 8 h (human two-compartment profile); combinations linezolid 4 mg/L + fosfomycin 4 g and linezolid 8 mg/L + fosfomycin 6 g",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "Strain-specific semimechanistic PK/PD model fit in NONMEM 7.4.1 (FOCE-I, ADVAN6) to dynamic in vitro time-kill data.",
      "No between-replicate variability, residual error or initial inoculum is reported.",
      "The printed Hill coefficient for linezolid (187) and the adaptation rates k LZD = k FOF = 0.0001 look like estimates at a bound; they are encoded as printed. The negative interaction exponent is discussed by the authors as the model judging the combination non-synergistic for this strain, contrary to the checkerboard and time-kill results.",
      "Companion files: Mao_2021_linezolid_fosfomycin_ATCC29212, _strain2, _strain6."
    )
  )

  ini({
    # ---- Fosfomycin pharmacokinetics (supplemental Table S1) ---------------
    # One parameter set for all three regimens (4, 6, 8 g). The printed Ke is
    # used as printed; the paper's stated in-system half-life, AUC0-8 and Cmax
    # back-solve to Ke ~0.42 1/h (see vignette).
    lkel <- log(0.22); label("Fosfomycin elimination rate constant from the central compartment (1/h)") # Table S1: Ke = 0.22 h-1
    lvc <- log(21.51); label("Fosfomycin central volume of distribution (L)") # Table S1: V1 = 21.51 L
    lk12 <- log(0.20); label("Fosfomycin central-to-peripheral rate constant (1/h)") # Table S1: K12 = 0.20 h-1 (V2 = 10.26 L is not needed in the rate-constant parameterisation)
    lk21 <- log(0.41); label("Fosfomycin peripheral-to-central rate constant (1/h)") # Table S1: K21 = 0.41 h-1

    # ---- Bacterial growth (Table 3, column No. 1) ---------------------
    lkgrow <- log(0.452); label("Net bacterial growth rate constant (1/h)") # Table 3: Kg = 0.452 h-1
    log10cfumax <- 8.55; label("Bacterial count in the stationary phase (log10 CFU/mL)") # Table 3: Bmax = 8.55 log10 CFU/mL
    log10cfu0 <- fixed(6); label("Initial inoculum (log10 CFU/mL)") # Methods, In vitro PK/PD model: 'the initial inoculation amount of bacteria was ~1 x 10^6 CFU/mL'; not estimated

    # ---- Linezolid drug effect (Table 3) -----------------------------------
    lemax_lzd <- log(0.191); label("Maximum kill rate constant of linezolid (1/h)") # Table 3: Emax-LZD = 0.191 h-1
    lec50_lzd <- log(3.95); label("Linezolid concentration giving half the maximum kill rate before adaptation (mg/L)") # Table 3: EC50-LZD = 3.95 mg/L
    hill_lzd <- 187; label("Hill coefficient of linezolid killing (unitless)") # Table 3: Hill LZD = 187
    lfadapt_lzd <- log(5.29); label("Maximal adaptive-resistance factor for linezolid (unitless)") # Table 3: f LZD = 5.29
    lkadapt_lzd <- log(0.0001); label("Rate of adaptation for linezolid (L/(mg*h))") # Table 3: k LZD = 0.0001; Equation 4 multiplies k by a concentration and time

    # ---- Fosfomycin drug effect (Table 3) ----------------------------------
    lemax_fof <- log(0.294); label("Maximum kill rate constant of fosfomycin (1/h)") # Table 3: Emax-FOF = 0.294 h-1
    lec50_fof <- log(0.0135); label("Fosfomycin concentration giving half the maximum kill rate before adaptation (mg/L)") # Table 3: EC50-FOF = 0.0135 mg/L
    hill_fof <- 3.63; label("Hill coefficient of fosfomycin killing (unitless)") # Table 3: Hill FOF = 3.63
    lfadapt_fof <- log(68600); label("Maximal adaptive-resistance factor for fosfomycin (unitless)") # Table 3: f FOF = 68,600
    lkadapt_fof <- log(0.0001); label("Rate of adaptation for fosfomycin (L/(mg*h))") # Table 3: k FOF = 0.0001

    # ---- Interaction (Table 3, Equation 5) ---------------------------------
    inter <- -1.96; label("Interaction exponent between linezolid and fosfomycin effects (unitless; >0 synergy)") # Table 3: Int = '2 1.96' in the PDF text layer, i.e. -1.96

    # ---- Residual error ----------------------------------------------------
    # Not reported (the paper shows goodness-of-fit and VPC plots only).
    addSd <- fixed(0); label("Additive residual SD on the log10 bacterial count (log10 CFU/mL)") # not reported in the paper; set to zero
  })

  model({
    # ---- 1. Parameters -----------------------------------------------------
    kel <- exp(lkel)
    vc <- exp(lvc)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kgrow <- exp(lkgrow)
    emax_lzd <- exp(lemax_lzd)
    ec50_lzd <- exp(lec50_lzd)
    fadapt_lzd <- exp(lfadapt_lzd)
    kadapt_lzd <- exp(lkadapt_lzd)
    emax_fof <- exp(lemax_fof)
    ec50_fof <- exp(lec50_fof)
    fadapt_fof <- exp(lfadapt_fof)
    kadapt_fof <- exp(lkadapt_fof)

    # ---- 2. Fosfomycin two-compartment PK (Table S1) -----------------------
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(lzd) <- 0
    Cc <- central / vc

    # ---- 3. Adaptive resistance (Equation 4) -------------------------------
    # alpha = 1 + f * (1 - exp(-C * k * t)), with C the current concentration
    # and t the time since the start of the experiment.
    alpha_lzd <- 1 + fadapt_lzd * (1 - exp(-lzd * kadapt_lzd * t))
    alpha_fof <- 1 + fadapt_fof * (1 - exp(-Cc * kadapt_fof * t))

    # ---- 4. Per-drug sigmoid Emax kill (Equation 3) ------------------------
    # Emax * C^H / ((alpha * EC50)^H + C^H) written as Emax / (1 + (alpha *
    # EC50 / C)^H), which is algebraically identical but does not overflow for
    # the large printed Hill coefficients of the clinical strains.
    e_lzd <- 0
    if (lzd > 0) {
      e_lzd <- emax_lzd / (1 + (alpha_lzd * ec50_lzd / lzd)^hill_lzd)
    }
    e_fof <- 0
    if (Cc > 0) {
      e_fof <- emax_fof / (1 + (alpha_fof * ec50_fof / Cc)^hill_fof)
    }

    # ---- 5. Combined effect (Equation 5, cited-source form) ----------------
    # E = E_LZD * (1 + E_FOF / (E_LZD + E_FOF))^int
    #   + E_FOF * (1 + E_LZD / (E_LZD + E_FOF))^int
    # Each bracket carries the OTHER drug's effect, as in the Mohamed 2016
    # control stream the paper cites for Equation 5; it reduces to the single
    # drug's effect in monotherapy. The paper prints the SAME drug in each
    # bracket, which would multiply a monotherapy effect by 2^int.
    e_tot <- e_lzd + e_fof
    ekill <- 0
    if (e_tot > 0) {
      ekill <- e_lzd * (1 + e_fof / e_tot)^inter + e_fof * (1 + e_lzd / e_tot)^inter
    }

    # ---- 6. Bacterial count (Equations 1-2) --------------------------------
    # B is the log10 viable count: Bmax is reported in log10 CFU/mL.
    log_cfu(0) <- log10cfu0
    d/dt(log_cfu) <- (kgrow * (1 - log_cfu / log10cfumax) - ekill) * log_cfu

    log_cfu ~ add(addSd)
  })
}
