Willmann_2021_fxiLica <- function() {
  description <- paste(
    "Population PK/PD model for the GalNAc-conjugated FXI antisense",
    "oligonucleotide FXI-LICA (BAY2976217, formerly ION-957943) after",
    "subcutaneous administration in healthy volunteers (study LICA-CS1),",
    "pooled with IONIS-FXIRX data to fit the FXI-activity indirect-",
    "response arm (Willmann 2021).",
    "Pharmacokinetics: two-compartment model with parallel first-order",
    "(fraction F1, rate ka) and zero-order (fraction 1 - F1, duration D2)",
    "subcutaneous absorption and first-order elimination from the central",
    "compartment. The final PK model retained an OMEGA block for the",
    "correlations between IIV on CL, V2, and Q (rho ~ 0.6-0.7; Table S3).",
    "Pharmacodynamics: indirect-response model on FXI activity with",
    "sigmoid-Imax inhibition (Imax fixed to 1) of the zero-order",
    "synthesis rate Kin = baseline * kout, driven by an effect-site",
    "concentration linked to the plasma central compartment via first-",
    "order equilibration ke0. FXI-LICA has an approximately 20-fold higher",
    "hepatocyte potency than IONIS-FXIRX (IC50 = 2.59 vs 167 ng/mL) owing",
    "to GalNAc-mediated ASGPR uptake into hepatocytes. The main-analysis",
    "FXI-LICA PK model was fitted to healthy-volunteer data only; the",
    "authors' ESRD simulations for FXI-LICA (Figures 3a, 4) borrow the",
    "IONIS-FXIRX ESRD covariate effects on CL, V3, and the effect-site",
    "driving concentration (keoPAT) under the assumption that the",
    "receptor-mediated hepatic uptake pathway is not altered in ESRD.",
    "These borrowed ESRD factors are exposed as fixed parameters in this",
    "model so that the FXI-LICA-in-ESRD scenario can be reproduced.",
    "The FXI-activity arm was fitted simultaneously with pooled ASO-CS1 +",
    "ASO-CS4 + LICA-CS1 data (shared kout, baseline, ke0, and Hill)."
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
    depot       = list(analyte = "FXI-LICA", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "FXI-LICA", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "FXI-LICA", units = "mg", specimen = "plasma", verified = FALSE),
    fxi         = list(analyte = "FXI activity", units = "mg", specimen = "blood cell", verified = FALSE),
    effect      = list(analyte = "FXI activity", units = "mg", specimen = "not applicable", verified = FALSE)
  )

  covariateData <- list(
    RRT_HEMODIAL_STATUS = list(
      description        = "Intermittent-hemodialysis treatment-status indicator (1 = subject with end-stage renal disease on hemodialysis, 0 = healthy volunteer)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer)",
      notes              = "Willmann 2021 did NOT fit an ESRD covariate on the FXI-LICA main PK model (Table S3; only healthy-volunteer data from LICA-CS1 contributed to the FXI-LICA PK fit). For the FXI-LICA-in-ESRD simulations reported in Figures 3a and 4, the authors borrowed the IONIS-FXIRX-derived ESRD factors on CL, V3, and the effect-site driving concentration under the explicit assumption that the receptor-mediated hepatocyte uptake pathway is not altered in ESRD (Methods / Discussion). This model file exposes the borrowed factors (fixed at the IONIS-FXIRX values) so that the paper's ESRD-simulation scenario can be reproduced by activating RRT_HEMODIAL_STATUS = 1; leaving the covariate at 0 (default) reproduces the pure healthy-volunteer LICA-CS1 fit.",
      source_name        = "ESRD"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 66L,
    n_studies        = 1L,
    age_range        = "mean 50.8 (SD 11.5) years",
    weight_range     = "mean 77.5 (SD 13.2) kg",
    sex_female_pct   = 34.8,
    race_ethnicity   = "Not reported in the primary publication.",
    disease_state    = "Healthy adult volunteers in the first-in-human LICA-CS1 phase I study (n = 66, 48 active drug + 18 placebo).",
    dose_range       = "Single subcutaneous doses of 40, 60, 80, or 120 mg; multiple subcutaneous doses of 10, 20, or 30 mg once weekly for 6 doses; or 80 mg every 4 weeks for 4 doses (Willmann 2021 Table S1).",
    regions          = "Phase I LICA-CS1 site(s); specific regions not summarised in Willmann 2021.",
    renal_function   = "Healthy volunteers with normal renal function.",
    n_concentrations = 773L,
    notes            = "Baseline demographics from Willmann 2021 Supplementary Information Table S1. Only active-drug subjects contribute to the FXI-LICA PK fit (48 active). Placebo subjects (18) contributed FXI-activity data to inform baseline and residual variability. The paper's Discussion notes limited FXI-LICA phase I data and that a full covariate analysis is deferred pending the phase IIb ESRD study."
  )

  ini({
    # ---------------------------------------------------------------------
    # PK STRUCTURAL PARAMETERS -- Willmann 2021 Supplementary Table S3
    # (final PK model for FXI-LICA, fitted separately to LICA-CS1 HV data).
    # ---------------------------------------------------------------------
    lcl  <- log(12.8)  ; label("Clearance CL in healthy volunteers (L/h)")                                    # Willmann 2021 Table S3: CL = 12.8 L/h (95% CI 11.4-14.2)
    lvc  <- log(48.9)  ; label("Central volume V2 (L)")                                                       # Willmann 2021 Table S3: V2 = 48.9 L (95% CI 41.3-56.5)
    lq   <- log(1.82)  ; label("Intercompartmental clearance Q (L/h)")                                       # Willmann 2021 Table S3: Q = 1.82 L/h (95% CI 1.53-2.10)
    lvp  <- log(924)   ; label("Peripheral volume V3 (L)")                                                    # Willmann 2021 Table S3: V3 = 924 L (95% CI 738-1109)
    lka  <- log(0.864) ; label("First-order subcutaneous absorption rate constant ka (1/h)")                 # Willmann 2021 Table S3: kA = 0.864 1/h (95% CI 0.659-1.07)
    ld1  <- log(9.38)  ; label("Duration of the zero-order subcutaneous absorption arm D2 (h)")              # Willmann 2021 Table S3: D2 = 9.38 h (95% CI 8.95-9.81)
    logitfdepot <- qlogis(0.620) ; label("Logit of the depot (first-order) absorbed fraction F1 (unitless; F1 = plogis(logitfdepot))") # Willmann 2021 Table S3: Logit(F1) = 0.490 -> F1 = 0.620 (95% CI 0.523-0.708)

    # ---------------------------------------------------------------------
    # ESRD COVARIATE EFFECTS -- BORROWED FROM IONIS-FXIRX (Willmann 2021
    # Methods / Discussion): the paper's FXI-LICA-in-ESRD simulations
    # apply the IONIS-FXIRX-derived ESRD factors to FXI-LICA under the
    # assumption that ESRD does not alter the receptor-mediated hepatic
    # uptake. These are fixed here (no IIV, values transcribed from
    # Table S2 for PK effects and Table S4 for the effect-site factor).
    # ---------------------------------------------------------------------
    e_esrd_cl <- fixed(0.530) ; label("Borrowed fractional reduction in CL for ESRD subjects (unitless; CL_ESRD = CL_HV * (1 - e_esrd_cl))")     # Willmann 2021 Table S2: ESRD on CL = 0.530; carried into FXI-LICA ESRD simulations
    e_esrd_vp <- fixed(0.379) ; label("Borrowed fractional reduction in V3 for ESRD subjects (unitless; V3_ESRD = V3_HV * (1 - e_esrd_vp))")     # Willmann 2021 Table S2: ESRD on V3 = 0.379; carried into FXI-LICA ESRD simulations
    le_esrd_effect <- fixed(log(0.329)) ; label("Borrowed log of the ESRD multiplicative factor keoPAT on the effect-site driving concentration (unitless)") # Willmann 2021 Table S4: keoPAT = 0.329; carried into FXI-LICA ESRD simulations

    # ---------------------------------------------------------------------
    # PD STRUCTURAL PARAMETERS -- Willmann 2021 Supplementary Table S4
    # (final PK/PD model, simultaneously fitted to pooled FXI-activity data
    # from ASO-CS1 + ASO-CS4 + LICA-CS1; shared baseline, kout, ke0, hill).
    # ---------------------------------------------------------------------
    imax     <- fixed(1.00)   ; label("Maximum fractional inhibition of Kin (unitless)")                       # Willmann 2021 Table S4: Imax = 1.00 (fixed)
    lic50    <- log(2.59)     ; label("IC50 of FXI-LICA at the effect site (ng/mL)")                                        # Willmann 2021 Table S4: IC50 FXI-LICA = 2.59 ng/mL (95% CI 2.18-3.00)
    lrbase   <- log(0.994)    ; label("Baseline FXI activity Baseline (U/mL)")                                              # Willmann 2021 Table S4: Baseline FXI activity = 0.994 U/mL (95% CI 0.971-1.02)
    lkout    <- log(0.00435)  ; label("First-order elimination rate constant kout for FXI activity (1/h)")                  # Willmann 2021 Table S4: kout = 0.00435 1/h (95% CI 0.00383-0.00488)
    lhill    <- log(1.50)     ; label("Sigmoid Imax Hill / gamma exponent on the effect-site concentration (unitless)")    # Willmann 2021 Table S4: gamma = 1.50 (95% CI 1.41-1.60)
    lke0     <- log(0.00115)  ; label("First-order equilibration rate constant ke0 between plasma and the effect site (1/h)") # Willmann 2021 Table S4: keo = 0.00115 1/h (95% CI 0.00102-0.00129)

    # ---------------------------------------------------------------------
    # INTER-INDIVIDUAL VARIABILITY -- Willmann 2021 Tables S3 (PK) and S4 (PD).
    # PK IIV includes an OMEGA BLOCK on CL, V2, Q (correlations rho ~ 0.6-0.7).
    # F1 IIV was tested but is NOT reported in Table S3 (only 5 IIV rows CL,
    # V2, Q, V3, kA are shown; F1 IIV is absent from Table S3), so F1 has
    # no IIV in this model (typical-value F1 for all subjects).
    # ---------------------------------------------------------------------
    # OMEGA BLOCK covariance parameterisation (Table S3):
    #   Estimate on-diagonal = omega^2 (log scale)
    #   Covariance CLxV2 = 0.172 (rho ~ 0.716)
    #   Covariance CLxQ  = 0.127 (rho ~ 0.622)
    #   Covariance V2xQ  = 0.189 (rho ~ 0.608)
    # nlmixr2 tilde syntax for a lower-triangular block: name + name + ... ~
    # c(var, cov, var, cov, cov, var, ...).
    etalcl + etalvc + etalq ~
      c(0.146,
        0.172, 0.312,
        0.127, 0.189, 0.234)      # Willmann 2021 Table S3: OMEGA BLOCK (CL, V2, Q); omega^2_CL = 0.146, cov(CL,V2) = 0.172, omega^2_V2 = 0.312, cov(CL,Q) = 0.127, cov(V2,Q) = 0.189, omega^2_Q = 0.234
    etalvp ~ 0.362              # Willmann 2021 Table S3: omega^2 V3 = 0.362 (%CV = 66.1)
    etalka ~ 0.194              # Willmann 2021 Table S3: omega^2 kA = 0.194 (%CV = 46.2)

    etalrbase ~ 0.0262          # Willmann 2021 Table S4: omega^2 Baseline = 0.0262 (%CV = 16.3)
    etalic50  ~ 0.281           # Willmann 2021 Table S4: omega^2 IC50 = 0.281 (%CV = 57.0)
    etalkout  ~ 0.115           # Willmann 2021 Table S4: omega^2 kout = 0.115 (%CV = 35.0)
    etalke0   ~ 0.109           # Willmann 2021 Table S4: omega^2 keo = 0.109 (%CV = 34.0)

    # ---------------------------------------------------------------------
    # RESIDUAL ERROR -- Willmann 2021 Tables S3 (PK) and S4 (PD, LICA-CS1
    # row). Both are proportional.
    # ---------------------------------------------------------------------
    propSd        <- 0.259    ; label("Proportional residual SD on Cc (fraction; LICA-CS1 healthy volunteers)")                  # Willmann 2021 Table S3: prop error = 0.259
    propSd_FXIact <- 0.101    ; label("Proportional residual SD on FXI activity (fraction; LICA-CS1 healthy volunteers)")       # Willmann 2021 Table S4: prop error LICA-CS1 = 0.101
  })

  model({
    # -------------------------------------------------------------------
    # Individual PK parameters (log-normal IIV; CL, V2, Q correlated via
    # the OMEGA BLOCK declared in ini()). ESRD subject-level factor on
    # CL and V3 is applied here for the paper's FXI-LICA-in-ESRD
    # simulation scenario; when RRT_HEMODIAL_STATUS = 0 the factors
    # reduce to 1 and this arm matches the pure LICA-CS1 fit.
    # -------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (1 - e_esrd_cl * RRT_HEMODIAL_STATUS)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp) * (1 - e_esrd_vp * RRT_HEMODIAL_STATUS)
    ka <- exp(lka + etalka)
    d1 <- exp(ld1)
    # Depot (first-order absorbed) fraction: F1 = plogis(logitfdepot),
    # bounded to (0, 1). No IIV on F1 in the FXI-LICA final model
    # (Table S3 does not report an omega^2 for F1). Central (zero-order
    # absorbed) fraction = 1 - F1.
    fdepot <- expit(logitfdepot)

    # -------------------------------------------------------------------
    # Individual PD parameters (log-normal IIV where reported)
    # -------------------------------------------------------------------
    rbase <- exp(lrbase + etalrbase)
    ic50  <- exp(lic50  + etalic50)
    kout  <- exp(lkout  + etalkout)
    ke0   <- exp(lke0   + etalke0)
    hill  <- exp(lhill)

    # Effect-site driving-concentration scaling for ESRD (borrowed from
    # IONIS-FXIRX per Willmann 2021 Methods / Discussion): equals 1 for
    # healthy volunteers and 0.329 for ESRD subjects when the borrowed
    # factor is applied.
    keo_esrd_factor <- exp(le_esrd_effect * RRT_HEMODIAL_STATUS)

    kin <- rbase * kout

    # -------------------------------------------------------------------
    # ODE system -- Willmann 2021 Figure 2 and Supplementary Information
    # Section S2 (NONMEM $DES).
    # -------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot + q / vp * peripheral1 - (cl + q) / vc * central
    d/dt(peripheral1) <-  q / vc * central - q / vp * peripheral1
    d/dt(fxi)         <-  kin * (1 - imax * effect^hill / (ic50^hill + effect^hill)) - kout * fxi
    d/dt(effect)      <-  ke0 * (keo_esrd_factor * (central / vc) * 1000 - effect)

    # Initial conditions: pre-dose FXI activity sits at baseline.
    fxi(0) <- rbase

    # -------------------------------------------------------------------
    # Bioavailability + absorption controls
    # -------------------------------------------------------------------
    f(depot)     <- fdepot
    f(central)   <- 1 - fdepot
    dur(central) <- d1

    # -------------------------------------------------------------------
    # Observations
    # -------------------------------------------------------------------
    # Plasma FXI-LICA concentration in ng/mL (dose enters in mg, vc in L,
    # so central/vc yields mg/L which is 1000 * ng/mL).
    Cc  <- central / vc * 1000

    # FXI activity in U/mL.
    FXIact <- fxi

    Cc     ~ prop(propSd)
    FXIact ~ prop(propSd_FXIact)
  })
}
