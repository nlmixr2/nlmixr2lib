Qin_2025_ropeginterferon_hct <- function() {
  description <- paste0(
    "Sequential pharmacokinetic-pharmacodynamic model for the ",
    "hematocrit response to subcutaneous ropeginterferon alfa-2b ",
    "(ropeg) in 78 Chinese and Japanese patients with polycythaemia ",
    "vera (Qin 2025, phase II studies A19-201 slow titration and ",
    "A20-202 fast titration). The quasi-equilibrium TMDD population ",
    "PK model of Qin_2025_ropeginterferon drives an indirect-response ",
    "turnover model in which total serum ropeg inhibits hematocrit ",
    "production through an Imax function, with a single transit ",
    "compartment interposed between production and the observed ",
    "hematocrit. Because the fit was SEQUENTIAL, every PK parameter ",
    "below is fixed at its Qin 2025 Table 2 value and only the ",
    "hematocrit parameters were estimated. The model carries ",
    "SEPARATE initial (HCT0 = 0.459) and steady-state (HCTss = ",
    "0.489) hematocrit values, so an untreated patient drifts upward ",
    "while ropeg drives the hematocrit down. Hematocrit is a ",
    "FRACTION, not a percentage. The Hill coefficient was tested and ",
    "NOT retained. Time is in DAYS despite the source table's h^-1 ",
    "labels; see the vignette Errata. Companion models in the ",
    "Qin_2025_ropeginterferon_* family."
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
    concentration = "ug/L total serum ropeg (Cc), numerically identical to the ng/mL Qin 2025 reports; hematocrit (hct) is a unitless volume FRACTION, 0-1"
  )


  compartmentData <- list(
    depot        = list(analyte = "ropeginterferon alfa-2b", units = "ug", specimen = "administration site", verified = FALSE),
    central      = list(analyte = "ropeginterferon alfa-2b", units = "ug", specimen = "serum", verified = FALSE),
    total_target = list(analyte = "ropeg target (total, free plus drug-bound binding capacity)", units = "ug/L (= ng/mL)", specimen = "serum", verified = FALSE),
    transit1     = list(analyte = "hematocrit (unobserved transit pool preceding the measured compartment)", units = "unitless volume fraction (0-1)", specimen = "whole blood", verified = FALSE),
    hct          = list(analyte = "hematocrit", units = "unitless volume fraction (0-1)", specimen = "whole blood", verified = FALSE)
  )

  covariateData <- list(
    BMI = list(
      description        = "Baseline body mass index.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Acts on ropeg clearance in the inherited PK layer only; no",
        "covariate was retained on any hematocrit parameter (Qin 2025",
        "Table 3 carries no covariate row). Power form P_i = P_TV *",
        "(BMI / 23.1)^0.813 from Equation (1), centred on the pooled",
        "median 23.1 kg/m^2 (Table 1, Overall). The 78 PV patients of",
        "this analysis had median BMI 21.2 (A19-201) and 23.9 (A20-202)."
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
        "TO 0 FOR EVERY SUBJECT IN THIS MODEL: the hematocrit analysis",
        "population is the 78 phase II PV patients, and no healthy",
        "volunteer contributed a hematocrit observation. The column is",
        "retained only so the inherited PK layer is byte-for-byte the",
        "same as Qin_2025_ropeginterferon."
      ),
      source_name        = "healthy-volunteer study flag; the source control stream name is not published"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 78L,
    n_studies      = 2L,
    n_observations = "hematocrit measured every 2 weeks in patients with PV (Qin 2025 Methods 2.3); the record count is not reported",
    age_range      = "median 54.0 years, range 26.0-72.0 (A19-201) and median 56.0 years, range 29.0-70.0 (A20-202) (Qin 2025 Table 1)",
    weight_range   = "median 56.0 kg, range 43.6-76.5 (A19-201) and median 67.9 kg, range 44.0-91.0 (A20-202) (Qin 2025 Table 1)",
    sex_female_pct = 43.6,
    race_ethnicity = "Japanese (A19-201, n = 29) and Chinese (A20-202, n = 49)",
    disease_state  = "polycythaemia vera. All A20-202 patients and all but two A19-201 patients carried JAK2 V617F; A20-202 enrolled patients resistant to or intolerant of hydroxyurea. Baseline hematocrit median 0.459 (A19-201, range 0.356-0.539)",
    dose_range     = "A19-201 (slow titration): 100 ug every 2 weeks, or 50 ug on prior cytoreductive therapy, titrated in 50 ug steps to a 500 ug maximum. A20-202 (fast titration): 250 ug at week 0, 350 ug at week 2, 500 ug from week 4",
    regions        = "Japan (A19-201) and China (A20-202)",
    notes          = paste0(
      "Complete hematologic response, the primary phase II efficacy ",
      "endpoint, required hematocrit < 45% without phlebotomy in the ",
      "previous 3 months, together with WBC < 10 x 10^9/L and PLT <= ",
      "400 x 10^9/L. Qin 2025 simulated a median time to first ",
      "hematocrit < 45% of 11 weeks under fast titration versus 18.3 ",
      "weeks under slow titration, a median per-patient difference of ",
      "5.43 weeks (95% prediction interval 0.129-9.24). Hematocrit was ",
      "the ONLY one of the three hematologic endpoints for which fast ",
      "titration produced a materially faster response."
    )
  )

  ini({
    # ==================================================================
    # PHARMACOKINETIC LAYER -- ENTIRELY FIXED.
    #
    # Qin 2025 Methods 2.4.1: "Sequential modeling was used to construct
    # separate PopPK-PD models". In a sequential fit the PK parameters
    # are carried over and held constant while the PD parameters are
    # estimated, so every value in this block -- including the
    # inter-individual variances -- is wrapped in fixed(). They are the
    # Qin 2025 Table 2 estimates, identical to the standalone
    # Qin_2025_ropeginterferon model; see that file for the full
    # source-trace and for the three proofs that the table's h^-1 unit
    # labels should read day^-1.
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

    # The PK RESIDUAL error terms of Table 2 are deliberately absent.
    # This is a SEQUENTIAL fit, so the objective function of the
    # hematocrit step runs over hematocrit observations only; the serum
    # concentration Cc enters as a derived driver, not as a fitted
    # endpoint. Use Qin_2025_ropeginterferon when a residual error model
    # for the concentration itself is wanted.

    # ==================================================================
    # HEMATOCRIT LAYER -- Qin 2025 Table 3, "HCT" block. These are the
    # only parameters the sequential fit estimated. All RSEs are below
    # 15% (Results 3.2.2).
    #
    # STRUCTURE (Qin 2025 Figure 1B, which prints both ODEs):
    #   dTransit/dt = kin * (1 - Imax*Cs/(IC50 + Cs)) - ktr * Transit
    #   dHCT/dt     = ktr * Transit - ktr * HCT
    # A single transit compartment sits between the inhibited zero-order
    # production and the observed hematocrit, and the SAME rate constant
    # ktr governs both transfers -- there is no separate kout for
    # hematocrit, and Table 3 accordingly reports no kout row for HCT.
    # Methods 2.4.1 states the production constant is "calculated"
    # rather than estimated, which pins kin = ktr * HCTss (the value
    # that makes HCTss the drug-free steady state of both states).
    #
    # HILL COEFFICIENT: tested and REJECTED. Results 3.2.2: "the
    # decreases in OFV after the inclusion of the Hill coefficient were
    # 0.016 and 5.654" for HCT and WBC respectively, "indicating that
    # the data were adequately described without the inclusion of the
    # Hill coefficient". The final model is therefore a plain Imax with
    # an implicit exponent of 1, and no gamma appears in Table 3.
    # ==================================================================
    lrbase    <- log(0.459) ; label("Initial hematocrit at the start of ropeg treatment, the initial condition of both hematocrit states (HCT0, fraction)")  # Qin 2025 Table 3, HCT block: HCT0 0.459 (RSE 1.41%). Discussion restates it as "the estimated initial HCT was 45.9%"
    lrbase_ss <- log(0.489) ; label("Drug-free steady-state hematocrit the system relaxes toward (HCTss, fraction)")                                         # Qin 2025 Table 3, HCT block: HCTss 0.489 (RSE 0.516%). Discussion restates it as "the equilibrium value was 48.9%". Note HCTss > HCT0, so an untreated patient drifts UPWARD
    lic50     <- log(137)   ; label("Total serum ropeg concentration producing half of Imax on hematocrit production (IC50, ng/mL)")                         # Qin 2025 Table 3, HCT block: IC50 137 (RSE 14.6%). Discussion restates it as "IC50,H estimated at 137 ng mL-1"
    limax     <- log(0.592) ; label("Maximum fractional inhibition of hematocrit production by ropeg (Imax, unitless fraction of kin)")                      # Qin 2025 Table 3, HCT block: Imax 0.592 (RSE 3.42%). Discussion restates it as "maximally reduce kin,H (Imax,H) by 59.2%"
    lktr      <- log(0.023) ; label("First-order transit rate constant governing both the transit-to-hematocrit and hematocrit-loss transfers (ktr, 1/day)")  # Qin 2025 Table 3, HCT block: ktr 0.023 (RSE 10.1%), printed as h-1; per day (see the vignette Errata)

    # ----- IIV on the hematocrit parameters -----
    # Table 3 reports IIV as a per-cent CV for log-normal random effects,
    # so omega^2 = log(CV^2 + 1). HCTss and IC50 are correlated; Table 3
    # gives their COVARIANCE directly (row "Covariance of IIV_HCTss and
    # IIV_IC50"), not a correlation coefficient. The implied correlation
    # is -0.298 / (sqrt(0.039221) * sqrt(2.969902)) = -0.873, which is a
    # legal correlation, so the printed number is admissible as stated.
    etalrbase_ss + etalic50 ~ c(0.039221,
                               -0.298, 2.969902) ; label("IIV variances of log HCTss and log IC50 with their covariance")  # Qin 2025 Table 3, HCT block: HCTss 20% CV (RSE 11.8%, shrinkage 14.8%) -> log(0.20^2 + 1) = 0.039221; IC50 430% CV (RSE 13.2%, shrinkage 15.1%) -> log(4.30^2 + 1) = 2.969902; covariance -0.298
    etalktr   ~ 0.822815 ; label("IIV variance on log ktr")   # Qin 2025 Table 3, HCT block: 113% CV (RSE 13.9%, shrinkage 17.8%); log(1.13^2 + 1) = 0.822815
    etalrbase ~ 0.010758 ; label("IIV variance on log HCT0")  # Qin 2025 Table 3, HCT block: 10.4% CV (RSE 9.27%, shrinkage 1.05%); log(0.104^2 + 1) = 0.010758

    propSd_hct <- 0.0409 ; label("Proportional residual SD on hematocrit (fraction)")  # Qin 2025 Table 3, HCT block: proportional residual error 4.09% (RSE 0.785%, shrinkage 5.64%). Table 3 reports no additive residual term for any of the three hematologic endpoints
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

    # ================= Hematocrit layer ==============================
    rbase    <- exp(lrbase + etalrbase)
    rbase_ss <- exp(lrbase_ss + etalrbase_ss)
    ic50     <- exp(lic50 + etalic50)
    imax     <- exp(limax)
    ktr      <- exp(lktr + etalktr)

    # Methods 2.4.1: the zero-order production constant is calculated,
    # not estimated. Setting kin = ktr * HCTss makes HCTss the drug-free
    # steady state of the transit compartment and hence of hematocrit.
    kin_hct <- ktr * rbase_ss

    # Qin 2025 Methods 2.4.1 names the driver CS, "total serum ropeg
    # concentrations", so the Imax term is driven by TOTAL drug and not
    # by the free concentration that drives linear elimination.
    inh_hct <- imax * Cc / (ic50 + Cc)

    # Qin 2025 Figure 1B, printed verbatim beneath the schematic.
    d/dt(transit1) <- kin_hct * (1 - inh_hct) - ktr * transit1
    d/dt(hct)      <- ktr * transit1 - ktr * hct

    # Both states start at HCT0. The paper reports HCT0 as the initial
    # value of hematocrit but does not print the transit compartment's
    # initial condition; starting the whole chain at HCT0 is the reading
    # that makes an untreated patient relax monotonically from HCT0 to
    # HCTss, and it reproduces the small early rise then fall of the
    # simulated median profile in Figure 2A. See the vignette
    # Assumptions and deviations.
    transit1(0) <- rbase
    hct(0)      <- rbase

    # ================= Observation ===================================
    # Hematocrit is the ONLY endpoint: the sequential fit estimated the
    # hematocrit parameters against hematocrit observations with the PK
    # held fixed, so Cc above is a derived driver rather than a second
    # fitted endpoint.
    hct ~ prop(propSd_hct)
  })
}
