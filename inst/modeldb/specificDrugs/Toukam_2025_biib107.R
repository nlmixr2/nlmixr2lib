Toukam_2025_biib107 <- function() {
  description <- "Two-compartment population PK model with parallel linear and Michaelis-Menten elimination, plus direct sigmoidal Emax PK/PD model of alpha-4 integrin receptor saturation, for BIIB107 (humanized aglycosyl anti-alpha-4 integrin IgG4 monoclonal antibody) in healthy adult volunteers (Toukam 2025)."
  reference <- "Toukam M, Karimian N, Bame E, Xu Y. Dose Optimization of BIIB107, an Anti-Alpha-4 Integrin Monoclonal Antibody, Through Population Pharmacokinetic and Pharmacodynamic Modeling. J Clin Pharmacol. 2026;66(1). doi:10.1002/jcph.70109 (PMID 41014552). Study NCT04593121."
  vignette <- "Toukam_2025_biib107"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL", response = "% alpha-4 integrin receptor saturation")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL, V2, V3, and Q with reference weight 70 kg (Toukam 2025 Table 3 and Methods 'Pki = theta_k * (Xij/M(Xj))^theta_j'). Body weight was the only covariate retained after stepwise selection (forward p<0.01 / backward p<0.001).",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = 62L,
    n_studies      = 1L,
    study_name     = "217HV101 (NCT04593121, first-in-human single + multiple ascending dose)",
    age_range      = "19-55 years",
    age_median     = "33.0 years",
    weight_range   = "57.0-99.7 kg",
    weight_median  = "78.7 kg",
    bmi_range      = "19.8-30.0 kg/m^2",
    sex_female_pct = 19.4,
    race_ethnicity = "White 71.0%; Black/African American 22.6%; Asian 4.8%; Other 29.0% (overlap with prior categories per source); Hispanic/Latino 50.0%; missing race 1.6%",
    disease_state  = "Healthy adult volunteers (BIIB107 is being developed for multiple sclerosis but the modeling dataset is the FIH study in healthy participants).",
    dose_range     = "SAD: 3.6, 18, 72, 180, 360, 600 mg SC and 360 mg IV (single dose); MAD 1B: 180 mg SC Q2W loading on Days 1, 15, 29 followed by 360 mg SC maintenance on Day 85; MAD 2B: 600 mg SC Q8W on Days 1 and 57.",
    regions        = "Not stated",
    pk_records     = "1151 PK samples in the modeling dataset (715 quantifiable + 436 BQL); 715 quantifiable used for the final model after exclusion of placebo and BQL records (M3 censored-likelihood approach was tested but caused instability and was dropped). PD records: 1077 alpha-4 integrin saturation observations.",
    ada_status     = "ADA-positive 58.1% (36/62) post-baseline; ADA was not retained as a PK or PD covariate.",
    notes          = "Baseline demographics from Toukam 2025 Table 2. ADA was assessed but not included in the final covariate model (no significant effect on PK or EC50). Renal function: 91.9% normal (CrCl >=90 mL/min), 8.1% mild impairment (60 <= CrCl < 90 mL/min)."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PK parameters (Toukam 2025 Table 3, final population PK
    # model). Reference subject is 70 kg adult; paper reports CL and Q in
    # mL/day, volumes in mL, Vmax in ug/day. Converted to L/day, L, and
    # mg/day for consistency with units list (concentration ug/mL = mg/L).
    # ------------------------------------------------------------------
    lka     <- log(0.288);          label("Absorption rate constant for SC dosing (Ka, 1/day)")             # Toukam 2025 Table 3 (Ka 0.288 1/day)
    lcl     <- log(0.159);          label("Linear clearance for a 70 kg adult (CL, L/day)")                 # Toukam 2025 Table 3 (CL 159 mL/day)
    lvc     <- log(3.01);           label("Central volume of distribution for a 70 kg adult (V2, L)")       # Toukam 2025 Table 3 (V2 3010 mL)
    lvp     <- log(1.18);           label("Peripheral volume of distribution for a 70 kg adult (V3, L)")    # Toukam 2025 Table 3 (V3 1180 mL)
    lq      <- log(0.301);          label("Inter-compartmental clearance for a 70 kg adult (Q, L/day)")     # Toukam 2025 Table 3 (Q 301 mL/day)
    lvmax   <- log(1.89);           label("Maximum rate of saturable target-mediated elimination (Vmax, mg/day)")  # Toukam 2025 Table 3 (Vmax 1890 ug/day)
    lkm     <- fixed(log(0.00435)); label("Michaelis-Menten constant (Km, mg/L = ug/mL; FIXED at in vitro Kd)")    # Toukam 2025 Table 3 (Km 0.00435 ug/mL FIXED)
    lfdepot <- log(0.738);          label("SC bioavailability (F, fraction)")                               # Toukam 2025 Table 3 (F 73.8%)
    lalag   <- fixed(log(0.0793));  label("Absorption lag time for SC doses (Tlag/ALAG1, day; FIXED)")      # Toukam 2025 Table 3 (Tlag 0.0793 day FIXED)

    # Allometric exponents on body weight (reference 70 kg). Toukam 2025
    # Table 3: exponent on CL was estimated at 1.07 (RSE 39.7%); exponents
    # on V2, V3, and Q were fixed at 1, 1, and 0.75 respectively.
    allo_cl <- 1.07;          label("Allometric exponent on CL (estimated, unitless)")  # Toukam 2025 Table 3 (Exponent on CL 1.07 estimated)
    allo_vc <- fixed(1);      label("Allometric exponent on V2 (FIXED, unitless)")      # Toukam 2025 Table 3 (Exponent on V2 1 FIXED)
    allo_vp <- fixed(1);      label("Allometric exponent on V3 (FIXED, unitless)")      # Toukam 2025 Table 3 (Exponent on V3 1 FIXED)
    allo_q  <- fixed(0.75);   label("Allometric exponent on Q (FIXED, unitless)")       # Toukam 2025 Table 3 (Exponent on Q 0.75 FIXED)

    # Inter-individual variability on log-transformed parameters. Paper
    # reports CV%; converted via omega^2 = log(CV^2 + 1). IIV was retained
    # only on CL, V2, and Ka (Toukam 2025 Table 3); no correlations were
    # reported between PK etas.
    etalcl ~ 0.02527  # Toukam 2025 Table 3 (IIV on CL 16% CV; omega^2 = log(1 + 0.16^2) = 0.02527)
    etalvc ~ 0.11556  # Toukam 2025 Table 3 (IIV on V2 35% CV; omega^2 = log(1 + 0.35^2) = 0.11556)
    etalka ~ 0.15536  # Toukam 2025 Table 3 (IIV on Ka 41% CV; omega^2 = log(1 + 0.41^2) = 0.15536)

    # PK residual error. Toukam 2025 Table 3 reports separate additive
    # error terms for SC (8.67) and IV (7.32) administration. Units are
    # not explicitly stated in the table; given concentrations in ug/mL,
    # the literal SC value is used here as a single additive SD. The
    # large magnitude vs. typical mAb assays (LLOQ 0.2 ug/mL) is flagged
    # in the validation vignette's Errata / Assumptions sections; route-
    # specific error coding is intentionally simplified to a single
    # additive term for library use.
    CcaddSd <- 8.67; label("Additive residual error on BIIB107 serum concentration (ug/mL, SC value)")  # Toukam 2025 Table 3

    # ------------------------------------------------------------------
    # Sigmoidal Emax PD model of alpha-4 integrin receptor saturation
    # (Toukam 2025 Table 4 and Methods PK-PD equation:
    #   E = E0 + Emax * Cc^gamma / (EC50^gamma + Cc^gamma)
    # with Cc = predicted central serum concentration in ug/mL).
    # ------------------------------------------------------------------
    a4satE0     <- 17.7;          label("Baseline alpha-4 integrin saturation (E0, %)")             # Toukam 2025 Table 4 (E0 17.7%)
    a4satEmax   <- 77.5;          label("Maximum drug-induced alpha-4 integrin saturation (Emax, % above E0)")  # Toukam 2025 Table 4 (Emax 77.5%)
    la4satEC50  <- log(0.376);    label("BIIB107 concentration producing 50% of Emax (EC50, ug/mL)")            # Toukam 2025 Table 4 (EC50 0.376 ug/mL)
    a4satGamma  <- fixed(1);      label("Hill coefficient for alpha-4 integrin saturation (gamma, FIXED)")     # Toukam 2025 PK-PD equation (gamma fixed at 1)

    # IIV on EC50 (only PD parameter with IIV). Paper reports 10 (with
    # high shrinkage 62%); interpreted as 10% CV consistent with the
    # other reported IIVs in Tables 3-4 (% CV).
    etala4satEC50 ~ 0.00995  # Toukam 2025 Table 4 (IIV on EC50 10% CV; omega^2 = log(1 + 0.10^2) = 0.00995)

    # PD residual error. Toukam 2025 Table 4 reports separate additive
    # error terms for SC (15.8%) and IV (7.93%) administration. The SC
    # value is used here for the single output.
    addSd_a4sat <- 15.8; label("Additive residual error on alpha-4 integrin saturation (%, SC value)")  # Toukam 2025 Table 4
  })

  model({
    # ------------------------------------------------------------------
    # Individual PK parameters with weight-based allometric scaling
    # relative to a 70 kg reference. Vmax is NOT weight-scaled in the
    # source paper (Toukam 2025 Methods: scaling applied to CL, V2, V3,
    # and Q only).
    # ------------------------------------------------------------------
    ka   <- exp(lka + etalka)
    cl   <- exp(lcl + etalcl) * (WT / 70)^allo_cl
    vc   <- exp(lvc + etalvc) * (WT / 70)^allo_vc
    vp   <- exp(lvp)          * (WT / 70)^allo_vp
    q    <- exp(lq)           * (WT / 70)^allo_q
    vmax <- exp(lvmax)
    km   <- exp(lkm)

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Concentration in central compartment (mg/L = ug/mL)
    Cc <- central / vc

    # ------------------------------------------------------------------
    # Two-compartment PK with first-order SC absorption (depot, with
    # bioavailability F and lag time ALAG1) and parallel linear plus
    # Michaelis-Menten elimination from central. IV doses are entered
    # directly into central by the user via the cmt column.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - vmax * Cc / (km + Cc) - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                                       k12 * central - k21 * peripheral1

    f(depot)    <- exp(lfdepot)
    alag(depot) <- exp(lalag)

    # ------------------------------------------------------------------
    # Direct sigmoidal Emax PD: alpha-4 integrin receptor saturation.
    # The PD model has no time delay (instantaneous binding observed in
    # the source data; no hysteresis), so the saturation depends only
    # on the concurrent central concentration.
    # ------------------------------------------------------------------
    a4satEC50_i <- exp(la4satEC50 + etala4satEC50)
    a4sat <- a4satE0 + a4satEmax * Cc^a4satGamma / (a4satEC50_i^a4satGamma + Cc^a4satGamma)

    # ------------------------------------------------------------------
    # Observation and residual-error models.
    # ------------------------------------------------------------------
    Cc    ~ add(CcaddSd)
    a4sat ~ add(addSd_a4sat)
  })
}
