Caballero_2023_isavuconazole_micafungin <- function() {
  description <- "In vitro (six Candida auris blood isolates). Semimechanistic PK/PD model of static-concentration time-kill by isavuconazole plus micafungin, alone and in combination, in RPMI medium over 48 h. The log10 fungal count grows logistically toward a maximum density, with growth delayed by (1 - exp(-klag * t)) at the start of the experiment; each drug kills through a sigmoid Emax term and the two kill rates combine through the empirical interaction function of Mohamed et al. 2016 (exponent `inter` > 0 = synergy). Drug concentrations are static bath states (mg/L). No between-isolate or between-occasion variability was retained. One of three combination-specific files from this paper (companions: Caballero_2023_isavuconazole_anidulafungin, Caballero_2023_isavuconazole_caspofungin)"
  reference <- "Caballero U, Eraso E, Quindos G, Vozmediano V, Schmidt S, Jauregizar N. PK/PD modeling and simulation of the in vitro activity of the combinations of isavuconazole with echinocandins against Candida auris. CPT Pharmacometrics Syst Pharmacol. 2023;12(6):770-782. doi:10.1002/psp4.12949. PMCID: PMC10272309. Structural model: Materials and Methods, Equations 1-3 and the NONMEM control stream in Data S3. Parameter estimates: Table 1, column ISV + MCF"
  vignette <- "Caballero_2023_isavuconazole_echinocandins"
  units <- list(time = "h", dosing = "mg/L", concentration = "mg/L")

  # The "dose" in a static time-kill assay is the bath concentration placed
  # directly into each drug state at time 0 (amt in mg/L).
  dosing <- c("cisv", "cmcf")

  # cisv / cmcf are the paper's own NONMEM data-column names (Data S3, $INPUT
  # CMCF CISV) and follow the `c<drug>` bath-concentration family (`cmem`,
  # `ctob`); they are not registered canonicals, so they are declared here.
  paper_specific_compartments <- c("cisv", "cmcf")

  compartmentData <- list(
    log_cfu = list(
      analyte = "Candida auris (six blood isolates pooled), total viable count",
      units = "log10 CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    cisv = list(analyte = "isavuconazole", units = "mg/L", specimen = "administration site", verified = TRUE),
    cmcf = list(analyte = "micafungin", units = "mg/L", specimen = "administration site", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "in vitro (Candida auris, six clinical blood isolates)",
    n_subjects = 6L,
    n_studies = 1L,
    organism = "Candida auris blood isolates CJ94, CJ97, CJ98, CJ99, CJ100 and CJ102 from an outbreak at the Hospital Universitario y Politecnico La Fe (Valencia, Spain); all from one clade closely related to the South African clade",
    mic_values = c(isavuconazole = "0.06 mg/L", micafungin = "0.125 mg/L"),
    system = "Static time-kill curves in flat-bottom microtiter plates, RPMI medium, 200 uL per well; viable counts at 0, 2, 4, 6, 8, 24 and 48 h plated in triplicate on Sabouraud dextrose agar; independent duplicate experiments",
    temperature = "37 C",
    duration = "48 h",
    starting_inoculum = "1-5 x 10^5 CFU/mL (nominal)",
    limit_of_quantification = "lower limit of detection 5 CFU/mL",
    dose_range = "isavuconazole 0.125, 0.25, 2 and 4 mg/L; micafungin 0.5-4 mg/L; combinations chosen from earlier checkerboard assays, each concentration also run as monotherapy, plus drug-free growth control",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "Fit in NONMEM 7.4.3 (FOCE with interaction, ADVAN13) to log10 CFU/mL from the monotherapy, combination and growth-control wells of all six isolates jointly.",
      "Neither between-isolate (IIV) nor between-experiment (IOV) variability improved the fit, so all OMEGAs are fixed to zero and only the additive residual error remains.",
      "All six isolates share the same MIC for each drug, so the model carries no isolate-level covariate."
    )
  )

  ini({
    # ---- Fungal growth (Table 1, column ISV + MCF) --------------------------
    lkgrow <- fixed(log(0.145)); label("Fungal growth rate constant (1/h)") # Table 1: k growth = 0.145 h-1 (fixed); Data S3 $THETA (0.145) FIX
    lklag <- log(0.158); label("Rate constant of the delay in fungal growth, alpha (1/h)") # Table 1: alpha = 0.158 (RSE 4%, 95% CI 0.145-0.171)
    log10cfumax <- fixed(8); label("Maximum fungal density Nmax (log10 CFU/mL)") # Table 1: N max = 8 log CFU/mL (fixed); Data S3 $THETA (8.00) FIX
    log10cfu0 <- fixed(5.35); label("Starting inoculum (log10 CFU/mL)") # Methods: starting inoculum 1-5 x 10^5 CFU/mL; log10 of the geometric midpoint sqrt(1 x 5) x 10^5 = 5.35, not a model estimate (Data S3 enters it per well through AMT)

    # ---- Isavuconazole kill (Table 1) ----------------------------------------
    lemax_isv <- log(0.0176); label("Maximum kill rate constant of isavuconazole (1/h)") # Table 1: E maxISV = 0.0176 h-1 (RSE 3%, 95% CI 0.0164-0.0186)
    lec50_isv <- log(0.0511); label("Isavuconazole concentration giving half the maximum kill rate (mg/L)") # Table 1: EC50ISV = 0.0511 mg/L (RSE 3%, 95% CI 0.0476-0.0543); the Results text quotes 0.0584, which lies outside the table's own bootstrap CI
    hill_isv <- 1.12; label("Hill coefficient of isavuconazole killing (unitless)") # Table 1: h ISV = 1.12 (RSE 6%)

    # ---- Micafungin kill (Table 1) ----------------------------------------
    lemax_mcf <- log(0.025); label("Maximum kill rate constant of micafungin (1/h)") # Table 1: E maxCANDIN = 0.025 h-1 (RSE 4%)
    lec50_mcf <- log(0.171); label("Micafungin concentration giving half the maximum kill rate (mg/L)") # Table 1: EC50CANDIN = 0.171 mg/L (RSE 9%)
    hill_mcf <- fixed(1); label("Hill coefficient of micafungin killing (unitless)") # Table 1: h CANDIN = 1 (fixed); Data S3 $THETA (1.0)FIX ; HILL1

    # ---- Pharmacodynamic interaction (Table 1, Equation 3) -------------------
    inter <- 0.41; label("Interaction exponent between the isavuconazole and micafungin kill rates (unitless; >0 synergy)") # Table 1: Int = 0.41 (RSE 18%, 95% CI 0.28-0.56)

    # ---- Residual error (Table 1) --------------------------------------------
    addSd <- 0.27; label("Additive residual SD on the log10 fungal count (log10 CFU/mL)") # Table 1: sigma = 0.27 log CFU/mL (RSE 6%)
  })

  model({
    # ---- 1. Back-transforms --------------------------------------------------
    kgrow <- exp(lkgrow)
    klag <- exp(lklag)
    emax_isv <- exp(lemax_isv)
    ec50_isv <- exp(lec50_isv)
    emax_mcf <- exp(lemax_mcf)
    ec50_mcf <- exp(lec50_mcf)

    # ---- 2. Static drug exposures (mg/L) -------------------------------------
    d/dt(cisv) <- 0
    d/dt(cmcf) <- 0

    # ---- 3. Per-drug sigmoid Emax kill rates (Equation 1) --------------------
    eff_isv <- emax_isv * cisv^hill_isv / (ec50_isv^hill_isv + cisv^hill_isv)
    eff_mcf <- emax_mcf * cmcf^hill_mcf / (ec50_mcf^hill_mcf + cmcf^hill_mcf)

    # ---- 4. Combined effect (Equation 3; Data S3 kanf + kisv, kanf = micafungin term)
    # Each bracket carries the OTHER drug's effect, so the expression reduces
    # to the single drug's kill rate in monotherapy. The guard replaces the
    # 0/0 of the drug-free control with zero kill; Data S3 reaches the same
    # result by switching the ODE on the TYPE data column.
    eff_tot <- eff_isv + eff_mcf
    ekill <- 0
    if (eff_tot > 0) {
      ekill <- eff_mcf * (1 + eff_isv / eff_tot)^inter +
        eff_isv * (1 + eff_mcf / eff_tot)^inter
    }

    # ---- 5. Log10 fungal count (Equations 1-2) -------------------------------
    # The state N is the log10 count itself (Nmax is in log CFU/mL and the
    # Data S3 prediction IPRED = A(1) is compared directly with log CFU/mL).
    # t is the time since the start of the experiment.
    log_cfu(0) <- log10cfu0
    d/dt(log_cfu) <- (kgrow * (1 - log_cfu / log10cfumax) * (1 - exp(-klag * t)) - ekill) * log_cfu

    log_cfu ~ add(addSd)
  })
}
