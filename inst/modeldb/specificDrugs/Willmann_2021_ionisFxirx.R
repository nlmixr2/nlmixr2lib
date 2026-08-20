Willmann_2021_ionisFxirx <- function() {
  description <- paste(
    "Population PK/PD model for the unconjugated FXI antisense oligonucleotide",
    "IONIS-FXIRX (BAY2306001) after subcutaneous administration in healthy",
    "volunteers (study ASO-CS1) and patients with end-stage renal disease",
    "requiring hemodialysis (study ASO-CS4), pooled with FXI-LICA data to fit",
    "the FXI-activity indirect-response arm (Willmann 2021).",
    "Pharmacokinetics: two-compartment model with parallel first-order",
    "(fraction F1, rate ka) and zero-order (fraction 1 - F1, duration D2)",
    "subcutaneous absorption and first-order elimination from the central",
    "compartment. Covariates on PK: end-stage renal disease reduces CL by",
    "53% and peripheral volume V3 by 38% (proportional NONMEM form);",
    "body weight enters V2 as a power (V2 propto (WT/70)^0.967).",
    "Pharmacodynamics: an indirect-response model on FXI activity with",
    "sigmoid-Imax inhibition (Imax fixed to 1) of the zero-order",
    "synthesis rate Kin = baseline * kout, driven by an effect-site",
    "concentration linked to the plasma central compartment via first-",
    "order equilibration ke0. An additional multiplicative factor keoPAT",
    "modulates the effect-site concentration in ESRD patients (Ce_ESRD =",
    "keoPAT * Cc at steady state; keoPAT = 0.329 reduces the effective",
    "driving concentration to ~one third of the healthy-volunteer value).",
    "The PK arm was fitted separately per compound; the FXI-activity arm",
    "was fitted simultaneously to pooled ASO-CS1 + ASO-CS4 + LICA-CS1",
    "data with shared kout, baseline, ke0, and Hill exponent (Table S4)."
  )
  reference <- paste(
    "Willmann S, Marostica E, Snelder N, Solms A, Jensen M, Lobmeyer M,",
    "Lensing AWA, Bethune C, Morgan E, Yu RZ, Wang Y, Jung SW, Geary R,",
    "Bhanot S. PK/PD modeling of FXI antisense oligonucleotides to bridge",
    "the dose-FXI activity relation from healthy volunteers to end-stage",
    "renal disease patients. CPT Pharmacometrics Syst Pharmacol.",
    "2021;10(8):890-901. doi:10.1002/psp4.12663"
  )
  vignette <- "Willmann_2021_fxi_aso"
  paper_specific_compartments <- c("fxi")
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "IONIS-FXIRX (BAY2306001)", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "IONIS-FXIRX (BAY2306001)", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "IONIS-FXIRX (BAY2306001)", units = "mg", specimen = "plasma", verified = FALSE),
    fxi         = list(analyte = "FXI activity", units = "mg", specimen = "plasma", verified = FALSE),
    effect      = list(analyte = "FXI activity", units = "mg", specimen = "not applicable", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Willmann 2021 Table S1: healthy-volunteer mean 72.7 (SD 10.5) kg in ASO-CS1; ESRD mean 87.8 (SD 22.7) kg in ASO-CS4. The paper's exploratory covariate analysis retained body weight only on V2 (central volume) as a power relationship of the form V2 = V2_ref * (WT/WT_ref)^e_wt_vc with e_wt_vc = 0.967 (Table S2). The paper does not state the reference weight; a 70 kg reference is used here per the standard NONMEM allometric-scaling convention (see vignette Errata for the audit trail).",
      source_name        = "WT"
    ),
    RRT_HEMODIAL_STATUS = list(
      description        = "Intermittent-hemodialysis treatment-status indicator (1 = subject with end-stage renal disease on hemodialysis, 0 = healthy volunteer)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer)",
      notes              = "Willmann 2021 Methods (study ASO-CS4): patients with ESRD on hemodialysis received IONIS-FXIRX every 28 days over 12 weeks. A dedicated PK cohort demonstrated that hemodialysis itself did not alter the PK or PD (single 300 mg dose given immediately post-dialysis versus immediately pre-dialysis produced equivalent concentration profiles), so ESRD enters the model as a subject-level treatment-status indicator rather than a per-session activity indicator. Retained in the final PK model on CL (proportional reduction of 53%) and on the peripheral volume V3 (proportional reduction of 38%) and in the final PK/PD model on the effect-site driving concentration through the keoPAT factor (0.329, ~ one-third the HV value; see Table S4 and Discussion). Enters the structural PK model as multiplicative factors (1 - e_esrd_cl) on CL and (1 - e_esrd_vp) on V3, and enters the effect-compartment equation as the multiplicative factor exp(le_esrd_effect * RRT_HEMODIAL_STATUS) on the plasma concentration driving the effect site (equivalently keoPAT^RRT_HEMODIAL_STATUS).",
      source_name        = "ESRD"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 137L,
    n_studies        = 2L,
    age_range        = "ASO-CS1 mean 45.4 (SD 11.9) years; ASO-CS4 mean 59.2 (SD 12.9) years",
    weight_range     = "ASO-CS1 mean 72.7 (SD 10.5) kg; ASO-CS4 mean 87.8 (SD 22.7) kg; ESRD cohort spanned 48.5-164 kg per Methods",
    sex_female_pct   = 33.6,
    race_ethnicity   = "Not reported in the primary publication.",
    disease_state    = "Healthy adult volunteers (ASO-CS1; n = 88, 66 active drug + 22 placebo) and adults with end-stage renal disease requiring hemodialysis (ASO-CS4; n = 49, 36 active drug + 13 placebo).",
    dose_range       = "ASO-CS1: single subcutaneous doses of 50, 100, 200, or 300 mg and multiple subcutaneous doses of 50, 100, 200, or 300 mg (three doses in week 1 then once weekly for weeks 2-6; 8 doses total). ASO-CS4: multiple subcutaneous doses of 200 or 300 mg every 28 days for up to 12 weeks.",
    regions          = "Phase I ASO-CS1 and phase II ASO-CS4 sites; specific regions not summarised in Willmann 2021.",
    renal_function   = "Two-strata cohort: ASO-CS1 healthy volunteers with normal renal function; ASO-CS4 adults on maintenance hemodialysis for ESRD. Hemodialysis itself did not alter PK or PD (dedicated PK-cohort sub-study).",
    n_concentrations = 2229L,
    notes            = "Baseline demographics summarised in Willmann 2021 Supplementary Information Table S1. Total PK/PD observations for IONIS-FXIRX studies (ASO-CS1 + ASO-CS4) = 1126 + 1103 = 2229. Only active-drug subjects contribute to the PK/PD fit (66 + 36 = 102 active); placebo subjects (22 + 13 = 35) contributed FXI-activity data to inform baseline and residual variability. Study ASO-CS1 subjects who received placebo had FXI-activity observations pooled with IONIS-FXIRX ASO-CS1 data (Table S4 residual-error rows). Doses were subcutaneous throughout."
  )

  ini({
    # ---------------------------------------------------------------------
    # PK STRUCTURAL PARAMETERS -- Willmann 2021 Supplementary Table S2
    # (final PK model for IONIS-FXIRX, fitted separately to ASO-CS1 + ASO-CS4).
    # ---------------------------------------------------------------------
    lcl  <- log(3.12)  ; label("Clearance CL in healthy volunteers (L/h)")                                    # Willmann 2021 Table S2: CL = 3.12 L/h (95% CI 2.86-3.38)
    lvc  <- log(9.55)  ; label("Central volume V2 at reference body weight 70 kg (L)")                       # Willmann 2021 Table S2: V2 = 9.55 L (95% CI 8.59-10.5)
    lq   <- log(0.202) ; label("Intercompartmental clearance Q (L/h)")                                       # Willmann 2021 Table S2: Q = 0.202 L/h (95% CI 0.170-0.233)
    lvp  <- log(164)   ; label("Peripheral volume V3 in healthy volunteers (L)")                             # Willmann 2021 Table S2: V3 = 164 L (95% CI 123-204)
    lka  <- log(0.191) ; label("First-order subcutaneous absorption rate constant ka (1/h)")                 # Willmann 2021 Table S2: kA = 0.191 1/h (95% CI 0.169-0.212)
    ld1  <- log(4.88)  ; label("Duration of the zero-order subcutaneous absorption arm D2 (h)")              # Willmann 2021 Table S2: D2 = 4.88 h (95% CI 4.53-5.22)
    logitfdepot <- qlogis(0.688) ; label("Logit of the depot (first-order) absorbed fraction F1 (unitless; F1 = plogis(logitfdepot))") # Willmann 2021 Table S2: Logit(F1) = 0.790 -> F1 = 0.688 (95% CI 0.563-0.790)

    # ---------------------------------------------------------------------
    # COVARIATE EFFECTS ON PK -- Willmann 2021 Table S2
    # ---------------------------------------------------------------------
    # ESRD enters CL and V3 as proportional reductions (NONMEM-style):
    #   CL = CL_HV * (1 - e_esrd_cl * RRT_HEMODIAL_STATUS)
    #   V3 = V3_HV * (1 - e_esrd_vp * RRT_HEMODIAL_STATUS)
    # e_esrd_cl and e_esrd_vp are the fractional reductions themselves
    # (0.530 -> CL reduced by 53% for ESRD; 0.379 -> V3 reduced by 38%).
    e_esrd_cl <- fixed(0.530) ; label("Fractional reduction in CL for ESRD subjects (unitless; CL_ESRD = CL_HV * (1 - e_esrd_cl))")   # Willmann 2021 Table S2: ESRD on CL = 0.530 (95% CI 0.454-0.607)
    e_esrd_vp <- fixed(0.379) ; label("Fractional reduction in V3 for ESRD subjects (unitless; V3_ESRD = V3_HV * (1 - e_esrd_vp))")   # Willmann 2021 Table S2: ESRD on V3 = 0.379 (95% CI 0.239-0.518)
    e_wt_vc   <- fixed(0.967) ; label("Power exponent on body weight for V2 (unitless; V2 = V2_ref * (WT/70)^e_wt_vc)")               # Willmann 2021 Table S2: BWT on V2 = 0.967 (95% CI 0.619-1.32); reference weight 70 kg per NONMEM allometric convention (paper does not state reference; see vignette Errata)

    # ---------------------------------------------------------------------
    # PD STRUCTURAL PARAMETERS -- Willmann 2021 Supplementary Table S4
    # (final PK/PD model, simultaneously fitted to pooled FXI-activity data
    # from ASO-CS1 + ASO-CS4 + LICA-CS1; shared baseline, kout, ke0, hill).
    # ---------------------------------------------------------------------
    imax     <- fixed(1.00)    ; label("Maximum fractional inhibition of Kin (unitless)")                      # Willmann 2021 Table S4: Imax = 1.00 (fixed)
    lic50    <- log(167)       ; label("IC50 of IONIS-FXIRX at the effect site (ng/mL)")                                    # Willmann 2021 Table S4: IC50 IONIS-FXIRX = 167 ng/mL (95% CI 139-195)
    lrbase   <- log(0.994)     ; label("Baseline FXI activity Baseline (U/mL)")                                             # Willmann 2021 Table S4: Baseline FXI activity = 0.994 U/mL (95% CI 0.971-1.02)
    lkout    <- log(0.00435)   ; label("First-order elimination rate constant kout for FXI activity (1/h)")                 # Willmann 2021 Table S4: kout = 0.00435 1/h (95% CI 0.00383-0.00488)
    lhill    <- log(1.50)      ; label("Sigmoid Imax Hill / gamma exponent on the effect-site concentration (unitless)")    # Willmann 2021 Table S4: gamma = 1.50 (95% CI 1.41-1.60)
    lke0     <- log(0.00115)   ; label("First-order equilibration rate constant ke0 between plasma and the effect site (1/h)") # Willmann 2021 Table S4: keo = 0.00115 1/h (95% CI 0.00102-0.00129)
    le_esrd_effect <- fixed(log(0.329)) ; label("Log of the ESRD multiplicative factor keoPAT on the effect-site driving concentration (unitless)") # Willmann 2021 Table S4: factor keoPAT = 0.329 (95% CI 0.242-0.415); enters as (keoPAT)^RRT_HEMODIAL_STATUS

    # ---------------------------------------------------------------------
    # INTER-INDIVIDUAL VARIABILITY -- Willmann 2021 Tables S2 (PK) and S4 (PD).
    # Estimates are reported as omega^2 (variance on the log scale) except for
    # F1 which is estimated on the logit scale (variance omega^2_F1 = 1.84 is
    # normal on the logit); adopt the paper's values verbatim.
    # OMEGA correlations were tested for IONIS-FXIRX PK but not retained in
    # the final model (unstable in the OMEGA BLOCK per Results / Table S2
    # footnote), so IIV terms are diagonal here.
    # ---------------------------------------------------------------------
    etalcl         ~ 0.0864     # Willmann 2021 Table S2: omega^2 CL = 0.0864 (%CV = 30.0)
    etalvc         ~ 0.130      # Willmann 2021 Table S2: omega^2 V2 = 0.130 (%CV = 37.2)
    etalq          ~ 0.281      # Willmann 2021 Table S2: omega^2 Q = 0.281 (%CV = 57.0)
    etalvp         ~ 0.348      # Willmann 2021 Table S2: omega^2 V3 = 0.348 (%CV = 64.5)
    etalka         ~ 0.102      # Willmann 2021 Table S2: omega^2 kA = 0.102 (%CV = 32.8)
    etalogitfdepot ~ 1.84       # Willmann 2021 Table S2: omega^2 F1 (on the logit scale) = 1.84 (%CV = 230 -- large IIV in logit units)

    etalrbase ~ 0.0262          # Willmann 2021 Table S4: omega^2 Baseline = 0.0262 (%CV = 16.3)
    etalic50  ~ 0.281           # Willmann 2021 Table S4: omega^2 IC50 = 0.281 (%CV = 57.0)
    etalkout  ~ 0.115           # Willmann 2021 Table S4: omega^2 kout = 0.115 (%CV = 35.0)
    etalke0   ~ 0.109           # Willmann 2021 Table S4: omega^2 keo = 0.109 (%CV = 34.0)

    # ---------------------------------------------------------------------
    # RESIDUAL ERROR -- Willmann 2021 Tables S2 and S4 (per-study proportional).
    # Default residual-error values are the ASO-CS1 healthy-volunteer study
    # (the larger IONIS-FXIRX dataset). Study ASO-CS4 (ESRD) reports higher
    # proportional error on both outputs (see vignette Errata); users
    # simulating an ESRD-specific dataset should override propSd and
    # propSd_fxi with the ASO-CS4 values (0.290 and 0.172, respectively).
    # ---------------------------------------------------------------------
    propSd        <- 0.229    ; label("Proportional residual SD on Cc (fraction; ASO-CS1 healthy volunteers)")                  # Willmann 2021 Table S2: prop error ASO-CS1 = 0.229
    propSd_FXIact <- 0.0851   ; label("Proportional residual SD on FXI activity (fraction; ASO-CS1 healthy volunteers)")        # Willmann 2021 Table S4: prop error ASO-CS1 = 0.0851
  })

  model({
    # -------------------------------------------------------------------
    # Individual PK parameters (log-normal IIV where reported)
    # -------------------------------------------------------------------
    # ESRD subject-level factor on CL, V3, and the effect-site driving
    # concentration. RRT_HEMODIAL_STATUS = 1 for the ASO-CS4 patients with
    # end-stage renal disease on hemodialysis; = 0 for healthy volunteers.
    cl <- exp(lcl + etalcl) * (1 - e_esrd_cl * RRT_HEMODIAL_STATUS)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp) * (1 - e_esrd_vp * RRT_HEMODIAL_STATUS)
    ka <- exp(lka + etalka)
    d1 <- exp(ld1)
    # Depot (first-order absorbed) fraction: F1 = plogis(logitfdepot + eta),
    # bounded to (0, 1). Central (zero-order absorbed) fraction = 1 - F1.
    fdepot <- expit(logitfdepot + etalogitfdepot)

    # -------------------------------------------------------------------
    # Individual PD parameters (log-normal IIV where reported)
    # -------------------------------------------------------------------
    rbase <- exp(lrbase + etalrbase)
    ic50  <- exp(lic50  + etalic50)
    kout  <- exp(lkout  + etalkout)
    ke0   <- exp(lke0   + etalke0)
    hill  <- exp(lhill)

    # Effect-site driving-concentration scaling for ESRD: at steady state
    # d/dt(effect) = 0 gives Ce = keoPAT * Cc. Encoded here as
    # keo_esrd_factor = keoPAT ^ RRT_HEMODIAL_STATUS = exp(le_esrd_effect
    # * RRT_HEMODIAL_STATUS), which equals 1 for HVs and 0.329 for ESRD.
    keo_esrd_factor <- exp(le_esrd_effect * RRT_HEMODIAL_STATUS)

    # Zero-order production of FXI (holds fxi = rbase at baseline when the
    # sigmoid inhibition term is zero).
    kin <- rbase * kout

    # -------------------------------------------------------------------
    # ODE system -- Willmann 2021 Figure 2 and Supplementary Information
    # Section S2 (NONMEM $DES). The paper's compartments 1-5 map to
    # depot, central, peripheral1, fxi, effect respectively; here the
    # `effect` compartment is placed after `fxi` in the d/dt ordering
    # because rxode2 assigns compartment slots by declaration order and
    # observations reference the paper-mechanistic FXI activity pool.
    # -------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot + q / vp * peripheral1 - (cl + q) / vc * central
    d/dt(peripheral1) <-  q / vc * central - q / vp * peripheral1
    d/dt(fxi)         <-  kin * (1 - imax * effect^hill / (ic50^hill + effect^hill)) - kout * fxi
    d/dt(effect)      <-  ke0 * (keo_esrd_factor * (central / vc) * 1000 - effect)

    # Initial conditions: pre-dose steady state -- FXI activity sits at the
    # baseline value; effect-compartment concentration is zero (no drug on
    # board); PK compartments are zero.
    fxi(0)    <- rbase

    # -------------------------------------------------------------------
    # Bioavailability + absorption controls
    # -------------------------------------------------------------------
    # The dose must be entered as TWO event rows per administration: one
    # to `depot` (first-order arm, fraction F1) and one to `central`
    # (zero-order arm, fraction 1 - F1, duration D2). See vignette for a
    # worked example event table.
    f(depot)     <- fdepot
    f(central)   <- 1 - fdepot
    dur(central) <- d1

    # -------------------------------------------------------------------
    # Observations
    # -------------------------------------------------------------------
    # Plasma IONIS-FXIRX concentration in ng/mL (dose enters in mg, vc in
    # L, so central/vc yields mg/L which is 1000 * ng/mL).
    Cc  <- central / vc * 1000

    # FXI activity in U/mL (readout of the paper-specific compartment).
    FXIact <- fxi

    Cc     ~ prop(propSd)
    FXIact ~ prop(propSd_FXIact)
  })
}
