Liesenfeld_2006_dabigatran_ECT <- function() {
  description <- paste(
    "Pharmacodynamic model for the prolongation of ecarin clotting time",
    "(ECT) by dabigatran in orthopaedic surgery patients receiving oral",
    "dabigatran etexilate after total hip replacement (Liesenfeld 2006",
    "BISTRO I PK-PD analysis). The concentration-ECT relationship is a",
    "single linear function whose slope decays exponentially from an",
    "initial value SLO0 to a final value SLO_F with rate constant KM;",
    "the baseline ECT also declines with time-since-surgery via a",
    "proportional inhibitory Emax form. Covariate analysis retained no",
    "demographic, comedication, or laboratory variables. The 2006 paper",
    "does not develop a PK model; the PK component embedded here is the",
    "Liesenfeld 2013 two-compartment dabigatran disposition with all PK",
    "thetas fixed so the model is self-contained for simulation. The",
    "2013 PK was fit in end-stage renal-disease subjects and will",
    "overestimate dabigatran exposure for the 2006 BISTRO I orthopaedic-",
    "surgery population; users targeting BISTRO I-style scenarios",
    "should override the PK thetas or supply observed concentrations",
    "to the PD layer.",
    sep = " "
  )
  reference <- paste(
    "Liesenfeld K-H, Schaefer HG, Troconiz IF, Tillmann C, Eriksson BI, Stangier J.",
    "Effects of the direct thrombin inhibitor dabigatran on ex vivo coagulation",
    "time in orthopaedic surgery patients: a population model analysis.",
    "Br J Clin Pharmacol. 2006 Nov;62(5):527-537.",
    "doi:10.1111/j.1365-2125.2006.02667.x.",
    "PK structure embedded for simulation convenience from Liesenfeld 2013;",
    "see modellib('Liesenfeld_2013_dabigatran').",
    sep = " "
  )
  vignette <- "Liesenfeld_2006_dabigatran"

  paper_specific_etas <- c("etalslope0")

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age (years)",
      units       = "years",
      type        = "continuous",
      notes       = "Screened by GAM + forward inclusion / backward elimination on ECT model parameters and not retained (Liesenfeld 2006 Methods, Results)."
    ),
    WT = list(
      description = "Body weight (kg)",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened and not retained (Liesenfeld 2006 Results)."
    ),
    SEXF = list(
      description = "Biological sex (1 = female)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (paper covariate 'gender') and not retained (Liesenfeld 2006 Results)."
    ),
    HT = list(
      description = "Height (cm)",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened and not retained (Liesenfeld 2006 Results)."
    ),
    BMI = list(
      description = "Body mass index (kg/m^2)",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened and not retained (Liesenfeld 2006 Results)."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened (paper covariate 'serum creatinine clearance') and not retained (Liesenfeld 2006 Results)."
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 287L,
    n_observations  = 5060L,
    n_studies       = 1L,
    age_range       = "35-88 years (mean 67)",
    weight_range    = "49-130 kg (mean 78)",
    sex_female_pct  = 52.6,
    race_ethnicity  = "Not reported in the source paper.",
    disease_state   = "Adults undergoing elective total hip replacement surgery receiving oral dabigatran etexilate for venous thromboembolism prophylaxis 4-8 h after surgery and continuing for 6-10 days.",
    dose_range      = "Oral dabigatran etexilate 12.5, 25, 50, 100, 150, 200, or 300 mg twice daily, or 150 or 300 mg once daily. Per-subgroup N ranged 20-46 (Liesenfeld 2006 Table 1).",
    regions         = "BISTRO I: multicentre, open-label, dose-escalating Phase IIa study (Eriksson et al., reference [8] of the source paper).",
    notes           = "Baseline demographics from Liesenfeld 2006 Table 1. Of 289 patients enrolled in BISTRO I, 287 contributed to the ECT analysis with 5060 paired PK-PD observations. Sampling per protocol: predose + 4 h on the day of surgery; predose + 2 h on each subsequent dosing day; frequent sampling at steady state (typically day 4) with predose + 0.5, 1, 2, 4, 8, 12 h in a 21-patient subset. ECT was measured on a Biomatic B10 coagulometer (Desaga). Covariates screened (gender, age, height, body mass index, serum creatinine clearance, standard clinical laboratory parameters, comedications) -- none retained in the final model."
  )

  ini({
    # =================================================================
    # Embedded Liesenfeld 2013 dabigatran PK (two-compartment with
    # first-order absorption and an absorption lag). ALL PK PARAMETERS
    # FIXED so the model is self-contained for simulation. See
    # modellib('Liesenfeld_2013_dabigatran') for the standalone PK
    # model with estimable parameters. The 2013 ESRD CL/F (12.4 L/h)
    # is low relative to non-ESRD adults and will overestimate exposure
    # for the BISTRO I orthopaedic-surgery population; users targeting
    # that population should override these thetas with population-
    # appropriate values or supply observed concentrations to the PD
    # layer instead.
    # =================================================================
    lcl     <- fixed(log(12.4));  label("Apparent clearance CL/F (L/h); from Liesenfeld 2013 Table 2 (ESRD population)")     # Liesenfeld 2013 Table 2 (fixed)
    lvc     <- fixed(log(531));   label("Apparent central volume V2/F (L); from Liesenfeld 2013 Table 2")                    # Liesenfeld 2013 Table 2 (fixed)
    lq      <- fixed(log(152));   label("Apparent inter-compartmental clearance Q/F (L/h); from Liesenfeld 2013 Table 2")    # Liesenfeld 2013 Table 2 (fixed)
    lvp     <- fixed(log(499));   label("Apparent peripheral volume V3/F (L); from Liesenfeld 2013 Table 2")                 # Liesenfeld 2013 Table 2 (fixed)
    lka     <- fixed(log(0.821)); label("First-order absorption rate constant ka (1/h); from Liesenfeld 2013 Table 2")       # Liesenfeld 2013 Table 2 (fixed)
    ltlag   <- fixed(log(1.67));  label("Absorption lag time (h); from Liesenfeld 2013 Table 2 (fed condition)")             # Liesenfeld 2013 Table 2 (fixed)
    lfdepot <- fixed(log(1.00));  label("Relative bioavailability anchor (fraction); from Liesenfeld 2013 Table 2")          # Liesenfeld 2013 Table 2 (fixed)

    # =================================================================
    # PD ECT model (Liesenfeld 2006 Equations 2, 3, 5; Table 3).
    #   ECT(C, T) = E0(T) + SLOP(T) * C                                (Eq. 2)
    #   E0(T)     = BAS0 * (1 - EM_BA * T / (ET50 + T))                (Eq. 3)
    #   SLOP(T)   = SLO_F + (SLO0 - SLO_F) * exp(-KM * T)              (Eq. 5)
    # T = time since first dose / time after surgery.
    # The paper reports ET50 in days; ET50 is converted to hours
    # (2.86 days * 24 h/day = 68.64 h) in the ini() value below.
    # KM is reported in Table 3 with the unit "h^-1" but the
    # corresponding text states a half-life of 1.1 day = 27 h, which is
    # consistent ONLY with KM in day^-1. Treating the table unit as a
    # transcription error (day^-1 is the value used), KM = 0.617 day^-1
    # = 0.02571 h^-1; the model uses the per-hour value.
    # =================================================================
    lrbase     <- log(28.0);            label("Initial baseline ECT at time 0 on the day of surgery, BAS0 (s)")                                # Liesenfeld 2006 Table 3 (RSE 0.49%)
    lslope0    <- log(0.377);           label("Initial slope of the ECT-concentration line at time 0, SLO0 (s per ng/mL)")                     # Liesenfeld 2006 Table 3 (RSE 2.18%)
    lslope_inf <- log(0.268);           label("Final slope of the ECT-concentration line at time infinity, SLO_F (s per ng/mL)")               # Liesenfeld 2006 Table 3 (RSE 1.49%)
    lkm        <- log(0.617 / 24);      label("Exponential decay rate constant of the slope, KM (1/h)")                                        # Liesenfeld 2006 Table 3 (KM = 0.617 day^-1; table-unit transcription error -- see file comment, RSE 13.55%)
    let50      <- log(2.86 * 24);       label("Time at which baseline ECT decreases by 50%, ET50 (h)")                                         # Liesenfeld 2006 Table 3 (ET50 = 2.86 days, RSE 13.50%)
    emba       <- 0.175;                label("Maximum fractional decrease in baseline ECT over time, EM_BA (unitless)")                       # Liesenfeld 2006 Table 3 (RSE 6.46%)

    # =================================================================
    # Inter-individual variability (Liesenfeld 2006 Table 3). Per the
    # footnote IIV is on E0 and SLOP (the time-varying values), not on
    # the initial-value parameters BAS0 / SLO0; encoding eta on the
    # typical-value lrbase / lslope0 is equivalent because the
    # time-decay factor is shared across subjects.
    #   E0   CV  8.2% -> omega^2 = log(1 + 0.082^2) = 0.006700
    #   SLOP CV 13.7% -> omega^2 = log(1 + 0.137^2) = 0.01859
    # =================================================================
    etalrbase  ~ 0.006700   # Liesenfeld 2006 Table 3 (IIV on E0, CV 8.2%)
    etalslope0 ~ 0.01859    # Liesenfeld 2006 Table 3 (IIV on SLOP, CV 13.7%)

    # =================================================================
    # Residual error. Liesenfeld 2006 Table 3 reports a proportional
    # error on ECT with sigma = 6.63% CV (RSE 6.83%).
    # =================================================================
    propSd <- 0.0663; label("Proportional residual error on ECT (fraction)")   # Liesenfeld 2006 Table 3 (RSE 6.83%)
  })

  model({
    # =================================================================
    # PK: Liesenfeld 2013 two-compartment dabigatran disposition. Doses
    # in mg, vc in L; central/vc has units of mg/L = ug/mL; multiply by
    # 1000 to express Cc in ng/mL to match the PD source values.
    # =================================================================
    ka     <- exp(lka)
    cl     <- exp(lcl)
    vc     <- exp(lvc)
    vp     <- exp(lvp)
    q      <- exp(lq)
    alag_d <- exp(ltlag)
    fdepot <- exp(lfdepot)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)    <- fdepot
    alag(depot) <- alag_d

    Cc <- 1000 * central / vc

    # =================================================================
    # PD ECT. Simulation time `t` (hours) is the time-after-surgery
    # proxy, consistent with the paper's convention of anchoring T to
    # the first dose 4-8 h after surgery. ET50 and KM in ini() are in
    # hour units (their per-day source values are noted in the
    # parameter comments).
    # =================================================================
    rbase     <- exp(lrbase  + etalrbase)
    slope0    <- exp(lslope0 + etalslope0)
    slope_inf <- exp(lslope_inf)
    km        <- exp(lkm)
    et50      <- exp(let50)

    rbase_t <- rbase * (1 - emba * t / (et50 + t))                       # Eq. 3
    slope_t <- slope_inf + (slope0 - slope_inf) * exp(-km * t)           # Eq. 5

    ECT <- rbase_t + slope_t * Cc                                        # Eq. 2
    ECT ~ prop(propSd)
  })
}
