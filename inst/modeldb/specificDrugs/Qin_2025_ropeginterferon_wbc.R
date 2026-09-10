Qin_2025_ropeginterferon_wbc <- function() {
  description <- paste0(
    "Sequential pharmacokinetic-pharmacodynamic model for the white ",
    "blood cell count response to subcutaneous ropeginterferon alfa-2b ",
    "(ropeg) in 78 Chinese and Japanese patients with polycythaemia ",
    "vera (Qin 2025, phase II studies A19-201 slow titration and ",
    "A20-202 fast titration). The quasi-equilibrium TMDD population PK ",
    "model of Qin_2025_ropeginterferon drives a one-compartment ",
    "indirect-response turnover model in which total serum ropeg ",
    "inhibits white cell production through an Imax function. Because ",
    "the fit was SEQUENTIAL, every PK parameter below is fixed at its ",
    "Qin 2025 Table 2 value and only the white cell parameters were ",
    "estimated. Imax is FIXED AT 1; the model carries separate initial ",
    "(WBC0 = 11.4) and steady-state (WBCss = 6.72) counts. The Hill ",
    "coefficient was tested and NOT retained. Qin 2025's Discussion ",
    "restates three of these values differently from Table 3; Table 3 ",
    "is used here and the conflict is recorded in the vignette Errata. ",
    "Time is in DAYS despite the source table's h^-1 labels. Companion ",
    "models in the Qin_2025_ropeginterferon_* family."
  )
  reference <- paste(
    "Qin A, Shimoda K, Suo S, Fu R, Kirito K, Wu D, Liao J, Chen H, Wu L,",
    "Su X, Gao Y, Sato T, Li Y, Zhang J, Shen W, Wang W, Zhang L, Jin J,",
    "Komatsu N.",
    "Population pharmacokinetics-pharmacodynamics and exposure-response of",
    "ropeginterferon alfa-2b in Chinese and Japanese patients with",
    "polycythemia vera.",
    "Pharmacol Res Perspect. 2025;13(3):e70109.",
    "doi:10.1002/prp2.70109.",
    sep = " "
  )
  vignette <- "Qin_2025_ropeginterferon"
  units <- list(
    time          = "day (NOT hour: Qin 2025 Tables 2 and 3 label the rate constants h^-1, but the values are per day; Methods 2.4.1 gives kout and kdec 'in day-1'. See the vignette Errata)",
    dosing        = "ug (micrograms of ropeginterferon alfa-2b, subcutaneous)",
    concentration = "ug/L total serum ropeg (Cc), numerically identical to the ng/mL Qin 2025 reports; white blood cell count (circ_wbc) in 10^9/L"
  )


  compartmentData <- list(
    depot        = list(analyte = "ropeginterferon alfa-2b", units = "ug", specimen = "administration site", verified = FALSE),
    central      = list(analyte = "ropeginterferon alfa-2b", units = "ug", specimen = "serum", verified = FALSE),
    total_target = list(analyte = "ropeg target (total, free plus drug-bound binding capacity)", units = "ug/L (= ng/mL)", specimen = "serum", verified = FALSE),
    circ_wbc          = list(analyte = "white blood cell count", units = "10^9/L", specimen = "whole blood", verified = FALSE)
  )

  covariateData <- list(
    BMI = list(
      description        = "Baseline body mass index.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Acts on ropeg clearance in the inherited PK layer only; no",
        "covariate was retained on any white cell parameter (Qin 2025",
        "Table 3 carries no covariate row). Power form P_i = P_TV *",
        "(BMI / 23.1)^0.813 from Equation (1), centred on the pooled",
        "median 23.1 kg/m^2 (Table 1, Overall)."
      ),
      source_name        = "BMI (body mass index)"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant indicator; 1 = healthy volunteer, 0 = patient with polycythaemia vera.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (patient with polycythaemia vera)",
      notes              = paste(
        "Inherited from the PK layer, where it gates the chronic decline",
        "of the target binding capacity (Qin 2025 Methods 2.4.1). SET IT",
        "TO 0 FOR EVERY SUBJECT IN THIS MODEL: the white cell analysis",
        "population is the 78 phase II PV patients. The column is",
        "retained only so the inherited PK layer matches",
        "Qin_2025_ropeginterferon exactly."
      ),
      source_name        = "healthy-volunteer study flag; the source control stream name is not published"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 78L,
    n_studies      = 2L,
    n_observations = "white blood cell count measured every 2 weeks in patients with PV (Qin 2025 Methods 2.3); the record count is not reported",
    age_range      = "median 54.0 years, range 26.0-72.0 (A19-201) and median 56.0 years, range 29.0-70.0 (A20-202) (Qin 2025 Table 1)",
    weight_range   = "median 56.0 kg, range 43.6-76.5 (A19-201) and median 67.9 kg, range 44.0-91.0 (A20-202) (Qin 2025 Table 1)",
    sex_female_pct = 43.6,
    race_ethnicity = "Japanese (A19-201, n = 29) and Chinese (A20-202, n = 49)",
    disease_state  = "polycythaemia vera. Baseline white blood cell count median 14.0 x 10^9/L, range 4.79-34.0 (A19-201) and median 8.93 x 10^9/L, range 3.81-67.4 (A20-202) (Qin 2025 Table 1) -- the A19-201 median is above the 10 x 10^9/L complete-hematologic-response threshold",
    dose_range     = "A19-201 (slow titration): 100 ug every 2 weeks, or 50 ug on prior cytoreductive therapy, titrated in 50 ug steps to a 500 ug maximum. A20-202 (fast titration): 250 ug at week 0, 350 ug at week 2, 500 ug from week 4",
    regions        = "Japan (A19-201) and China (A20-202)",
    notes          = paste0(
      "White blood cell count < 10 x 10^9/L is one of the three ",
      "components of complete hematologic response. As for platelets, ",
      "Qin 2025 Results 3.2.3 found the simulated fast- and ",
      "slow-titration white cell profiles nearly indistinguishable ",
      "(Figure 2C, ANOVA p = 0.105), in contrast to the clear ",
      "separation for hematocrit. A decrease in white blood cell count ",
      "is also one of the ten most frequent treatment-related adverse ",
      "events carried into the exposure-safety screen, though it did ",
      "not reach a significant exposure-response relationship."
    )
  )

  ini({
    # ==================================================================
    # PHARMACOKINETIC LAYER -- ENTIRELY FIXED.
    #
    # Qin 2025 Methods 2.4.1: "Sequential modeling was used to construct
    # separate PopPK-PD models". The PK parameters are carried over and
    # held constant while the PD parameters are estimated, so every
    # value here -- including the inter-individual variances -- is
    # wrapped in fixed(). They are the Qin 2025 Table 2 estimates; see
    # Qin_2025_ropeginterferon for the full source-trace and for the
    # three proofs that the table's h^-1 unit labels should read day^-1.
    # ==================================================================
    lka     <- fixed(log(0.18))    ; label("First-order absorption rate constant from the subcutaneous depot (Ka, 1/day)")   # Qin 2025 Table 2: Ka 0.18 (RSE 1.04%); inherited unchanged by the sequential PD fit
    ltlag   <- fixed(log(0.62))    ; label("Absorption lag time on the subcutaneous depot (ALAG, day)")                      # Qin 2025 Table 2: 0.62 (RSE 0.421%)
    lcl     <- fixed(log(0.753))   ; label("Linear clearance at the median BMI of 23.1 kg/m^2 (CL, L/day)")                  # Qin 2025 Table 2: CL 0.753 (RSE 0.846%)
    lvc     <- fixed(log(3.29))    ; label("Volume of the serum compartment (Vc, L)")                                        # Qin 2025 Table 2: Vc 3.29 (RSE 1.03%)
    lrtot0  <- fixed(log(0.317))   ; label("Baseline maximum target binding capacity (Rtot0, ng/mL)")                        # Qin 2025 Table 2: Rtot0 0.317 (RSE 0.904%)
    lrtot_ss <- fixed(log(0.012))  ; label("Steady-state maximum target binding capacity under chronic dosing (Rtot,SS, ng/mL)")  # Qin 2025 Table 2: Rtot,SS 0.012 (RSE 0.997%)
    lkint   <- fixed(log(0.0223))  ; label("First-order elimination rate constant of the drug-target complex (kint, 1/day)") # Qin 2025 Table 2: kint 0.0223 (RSE 0.88%)
    lkdeg   <- fixed(log(0.51))    ; label("First-order degradation rate constant of free target (kdeg, 1/day)")             # Qin 2025 Table 2: kdeg 0.51 (RSE 0.576%)
    lkd     <- fixed(log(0.0662))  ; label("Equilibrium dissociation constant of ropeg for its target (KD, ng/mL)")          # Qin 2025 Table 2: KD 0.0662 (RSE 0.981%)
    lkdecay <- fixed(log(0.0255))  ; label("First-order rate constant of the decline in binding capacity (kdec, 1/day)")     # Qin 2025 Table 2: kdec 0.0255 (RSE 0.892%); Methods 2.4.1 states this one is "in day-1"
    t_start <- fixed(7)            ; label("Time after the first dose at which target mediation begins in patients (TSTART, day)")  # Qin 2025 Table 2: TSTART 7 (FIX)
    e_bmi_cl <- fixed(0.813)       ; label("Power exponent on (BMI / 23.1) for clearance (unitless)")                        # Qin 2025 Table 2: BMI effect on CL 0.813 (RSE 9.39%)

    etalka    ~ fixed(0.397837) ; label("IIV variance on log Ka")     # Qin 2025 Table 2: 69.9% CV; log(0.699^2 + 1)
    etalcl    ~ fixed(0.117435) ; label("IIV variance on log CL")     # Qin 2025 Table 2: 35.3% CV; log(0.353^2 + 1)
    etalvc    ~ fixed(0.639175) ; label("IIV variance on log Vc")     # Qin 2025 Table 2: 94.6% CV; log(0.946^2 + 1)
    etalrtot0 ~ fixed(2.128041) ; label("IIV variance on log Rtot0")  # Qin 2025 Table 2: 272% CV; log(2.72^2 + 1)
    etalkint  ~ fixed(0.940983) ; label("IIV variance on log kint")   # Qin 2025 Table 2: 125% CV; log(1.25^2 + 1)

    # The PK RESIDUAL error terms of Table 2 are deliberately absent:
    # this is a SEQUENTIAL fit, so the objective function of the white
    # cell step runs over white cell observations only and Cc enters as
    # a derived driver, not as a fitted endpoint.

    # ==================================================================
    # WHITE BLOOD CELL LAYER -- Qin 2025 Table 3, "WBC" block. All RSEs
    # below 15% (Results 3.2.2).
    #
    # STRUCTURE (Qin 2025 Figure 1C, which prints the ODE verbatim):
    #   dWBC/dt = kin,W * (1 - Imax,W*Cs/(IC50,W + Cs)) - kout,W * WBC
    # A plain one-compartment indirect response with production
    # inhibition. Methods 2.4.1 states the zero-order production
    # constant is "calculated", which pins kin,W = kout,W * WBCss.
    #
    # DISCUSSION-VERSUS-TABLE CONFLICT. Qin 2025's Discussion restates
    # the white cell parameters as "the estimated typical initial WBC
    # was 12.1 x 10^9 L-1 and the equilibrium value was 7.32 x 10^9 L-1
    # ... the corresponding IC50,H estimated at 124 ng mL-1", against
    # Table 3's 11.4, 6.72 and 152. Table 3 is used here because it is
    # the designated parameter-estimates table and carries RSEs, and
    # because that Discussion paragraph is demonstrably degraded: it
    # labels the platelet and white-cell quantities with the HEMATOCRIT
    # subscript H ("kin,H", "IC50,H") in all three paragraphs. For
    # hematocrit, where the subscript is correct, the Discussion and
    # Table 3 agree to every printed digit. Recorded in the vignette
    # Errata.
    #
    # UNIT-LABEL TYPO. Table 3 labels WBCss and WBC0 "10^9 L-1 day-1",
    # a RATE. They are counts (10^9 L-1); only kin,W, which the table
    # does not report, carries the day-1.
    #
    # HILL COEFFICIENT: tested and REJECTED. Results 3.2.2: "the
    # decreases in OFV after the inclusion of the Hill coefficient were
    # 0.016 and 5.654" for HCT and WBC respectively, "indicating that
    # the data were adequately described without the inclusion of the
    # Hill coefficient". The final model is a plain Imax with an
    # implicit exponent of 1, and no gamma appears in Table 3.
    # ==================================================================
    lrbase    <- log(11.4)  ; label("Initial white blood cell count at the start of ropeg treatment, the initial condition of the white cell state (WBC0, 10^9/L)")  # Qin 2025 Table 3, WBC block: WBC0 11.4 (RSE 6.95%). The Discussion restates 12.1; see the conflict note above
    lrbase_ss <- log(6.72)  ; label("Drug-free steady-state white blood cell count the system relaxes toward (WBCss, 10^9/L)")                                       # Qin 2025 Table 3, WBC block: WBCss 6.72 (RSE 6.42%). The Discussion restates 7.32; see the conflict note above
    lic50     <- log(152)   ; label("Total serum ropeg concentration producing half of Imax on white cell production (IC50,W, ng/mL)")                               # Qin 2025 Table 3, WBC block: IC50,W 152 (RSE 1.18%). The Discussion restates 124; see the conflict note above. This is the HIGHEST IC50 of the three hematologic endpoints
    limax     <- fixed(log(1)) ; label("Maximum fractional inhibition of white cell production by ropeg (Imax,W, unitless fraction of kin)")                         # Qin 2025 Table 3, WBC block: Imax,W 1 (FIX). Complete suppression of production is attainable
    lkout     <- log(0.0475); label("First-order elimination rate constant of circulating white blood cells (kout,W, 1/day)")                                        # Qin 2025 Table 3, WBC block: kout,W 0.0475 (RSE 5.28%), printed as h-1; per day, which gives a 15-day white cell pool half-life (see the vignette Errata)

    # ----- IIV on the white cell parameters -----
    # Table 3 reports IIV as a per-cent CV for log-normal random effects,
    # so omega^2 = log(CV^2 + 1). WBCss and IC50,W are correlated; Table
    # 3 gives their COVARIANCE directly, not a correlation coefficient.
    # The implied correlation is -0.19 / (sqrt(0.202676) *
    # sqrt(0.713146)) = -0.500, a legal correlation, so the printed
    # number is admissible as stated.
    etalrbase_ss + etalic50 ~ c(0.202676,
                               -0.190, 0.713146) ; label("IIV variances of log WBCss and log IC50,W with their covariance")  # Qin 2025 Table 3, WBC block: WBCss 47.4% CV (RSE 8.09%, shrinkage 0%) -> log(0.474^2 + 1) = 0.202676; IC50,W 102% CV (RSE 11.7%, shrinkage 16%) -> log(1.02^2 + 1) = 0.713146; covariance -0.19
    etalkout  ~ 0.999197 ; label("IIV variance on log kout,W")  # Qin 2025 Table 3, WBC block: 131% CV (RSE 9.53%, shrinkage 9.44%); log(1.31^2 + 1) = 0.999197
    etalrbase ~ 0.329753 ; label("IIV variance on log WBC0")    # Qin 2025 Table 3, WBC block: 62.5% CV (RSE 9.88%, shrinkage 1.44%); log(0.625^2 + 1) = 0.329753

    propSd_circ_wbc <- 0.16 ; label("Proportional residual SD on white blood cell count (fraction)")  # Qin 2025 Table 3, WBC block: proportional residual error 16% (RSE 0.633%, shrinkage 5.29%). Table 3 reports no additive residual term for any of the three hematologic endpoints
  })

  model({
    # ================= Inherited PK layer (fixed) =====================
    ka      <- exp(lka + etalka)
    tlag    <- exp(ltlag)
    cl      <- exp(lcl + etalcl) * (BMI / 23.1)^e_bmi_cl
    vc      <- exp(lvc + etalvc)
    rtot0   <- exp(lrtot0 + etalrtot0)
    rtot_ss <- exp(lrtot_ss)
    kint    <- exp(lkint + etalkint)
    kdeg    <- exp(lkdeg)
    kd      <- exp(lkd)
    kdecay  <- exp(lkdecay)

    tdec     <- max(t - t_start, 0)
    rtot_cap <- rtot0 + (1 - DIS_HEALTHY) * (rtot_ss - rtot0) * (1 - exp(-kdecay * tdec))
    ksyn     <- kdeg * rtot_cap

    ctot    <- central / vc
    disc    <- ctot - total_target - kd
    cfree   <- 0.5 * (disc + sqrt(disc * disc + 4 * kd * ctot))
    complex <- total_target * cfree / (kd + cfree)

    total_target(0) <- rtot0

    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot - cl * cfree - kint * complex * vc
    d/dt(total_target) <-  ksyn - kdeg * (total_target - complex) - kint * complex
    alag(depot)        <- tlag

    Cc <- ctot

    # ================= White blood cell layer ========================
    rbase    <- exp(lrbase + etalrbase)
    rbase_ss <- exp(lrbase_ss + etalrbase_ss)
    ic50     <- exp(lic50 + etalic50)
    imax     <- exp(limax)
    kout     <- exp(lkout + etalkout)

    # Methods 2.4.1: the zero-order production constant is calculated,
    # not estimated. kin,W = kout,W * WBCss makes WBCss the drug-free
    # steady state.
    kin_wbc <- kout * rbase_ss

    # Methods 2.4.1 names the driver CS, "total serum ropeg
    # concentrations", so the Imax term is driven by TOTAL drug and not
    # by the free concentration that drives linear elimination.
    inh_wbc <- imax * Cc / (ic50 + Cc)

    # Qin 2025 Figure 1C, printed verbatim beneath the schematic.
    d/dt(circ_wbc) <- kin_wbc * (1 - inh_wbc) - kout * circ_wbc
    circ_wbc(0)    <- rbase

    # ================= Observation ===================================
    # White blood cell count is the ONLY endpoint: the sequential fit
    # estimated the white cell parameters against white cell
    # observations with the PK held fixed, so Cc above is a derived
    # driver rather than a second fitted endpoint.
    circ_wbc ~ prop(propSd_circ_wbc)
  })
}
