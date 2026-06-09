Liesenfeld_2006_dabigatran_aPTT <- function() {
  description <- paste(
    "Pharmacodynamic model for the prolongation of activated partial",
    "thromboplastin time (aPTT) by dabigatran in orthopaedic surgery",
    "patients receiving oral dabigatran etexilate after total hip",
    "replacement (Liesenfeld 2006 BISTRO I PK-PD analysis). The",
    "concentration-aPTT relationship combines a linear and an Emax model;",
    "the baseline aPTT and the maximum nonlinear effect Emax both decline",
    "with time since surgery via a proportional inhibitory Emax form",
    "sharing a single ET50. Covariate analysis retained no demographic,",
    "comedication, or laboratory variables. The 2006 paper does not",
    "develop a PK model; the PK component embedded here is the Liesenfeld",
    "2013 two-compartment dabigatran disposition with all PK thetas fixed",
    "so the model is self-contained for simulation. The 2013 PK was fit",
    "in end-stage renal-disease subjects and will overestimate dabigatran",
    "exposure for the 2006 BISTRO I orthopaedic-surgery population;",
    "users targeting BISTRO I-style scenarios should override the PK",
    "thetas or supply observed concentrations to the PD layer.",
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

  paper_specific_etas <- c("etalslope")

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age (years)",
      units       = "years",
      type        = "continuous",
      notes       = "Screened by GAM + forward inclusion / backward elimination on aPTT model parameters and not retained (Liesenfeld 2006 Methods, Results)."
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
    n_observations  = 4854L,
    n_studies       = 1L,
    age_range       = "35-88 years (mean 67)",
    weight_range    = "49-130 kg (mean 78)",
    sex_female_pct  = 52.6,
    race_ethnicity  = "Not reported in the source paper.",
    disease_state   = "Adults undergoing elective total hip replacement surgery receiving oral dabigatran etexilate for venous thromboembolism prophylaxis 4-8 h after surgery and continuing for 6-10 days.",
    dose_range      = "Oral dabigatran etexilate 12.5, 25, 50, 100, 150, 200, or 300 mg twice daily, or 150 or 300 mg once daily. Per-subgroup N ranged 20-46 (Liesenfeld 2006 Table 1).",
    regions         = "BISTRO I: multicentre, open-label, dose-escalating Phase IIa study (Eriksson et al., reference [8] of the source paper).",
    notes           = "Baseline demographics from Liesenfeld 2006 Table 1. Of 289 patients enrolled in BISTRO I, 287 contributed to the aPTT analysis (two were excluded for missing plasma concentration data). 4854 paired PK-PD observations were available. Sampling per protocol: predose + 4 h on the day of surgery; predose + 2 h on each subsequent dosing day; frequent sampling at steady state (typically day 4) with predose + 0.5, 1, 2, 4, 8, 12 h in a 21-patient subset. Concentration assay: LC-MS/MS; lower limit of quantification 1.0 ng/mL. Covariates screened (gender, age, height, body mass index, serum creatinine clearance, standard clinical laboratory parameters, comedications including diuretics, opioids, NSAIDs, GI-transit accelerators, acetaminophen, CYP3A4 inhibitors) -- none retained in the final model."
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
    # PD aPTT model (Liesenfeld 2006 Equations 1, 3, 4; Table 2).
    #   aPTT(C, T) = E0(T) + Emax(T) * C / (EC50 + C) + SLOP * C  (Eq. 1)
    #   E0(T)   = BAS0 * (1 - EM_BA * T / (ET50 + T))             (Eq. 3)
    #   Emax(T) = EMA0 * (1 - EM_MX * T / (ET50 + T))             (Eq. 4)
    # T = time since first dose / time after surgery.
    # The paper reports ET50 in days; the model uses simulation time in
    # hours, so ET50 is converted (1.62 days * 24 h/day = 38.88 h) in
    # the ini() value below.
    # =================================================================
    lrbase  <- log(33.4);          label("Initial baseline aPTT at time 0 on the day of surgery, BAS0 (s)")                 # Liesenfeld 2006 Table 2 (RSE 0.63%)
    lemax   <- log(26.9);          label("Initial Emax at time 0 on the day of surgery, EMA0 (s)")                          # Liesenfeld 2006 Table 2 (RSE 12.45%)
    lec50   <- log(94.7);          label("EC50 of the nonlinear Emax part (ng/mL)")                                         # Liesenfeld 2006 Table 2 (RSE 17.11%)
    lslope  <- log(0.0509);        label("Slope of the linear part SLOP (s per ng/mL)")                                     # Liesenfeld 2006 Table 2 (RSE 6.68%)
    let50   <- log(1.62 * 24);     label("Time at which baseline and Emax each decrease by 50%, ET50 (h)")                  # Liesenfeld 2006 Table 2 (ET50 = 1.62 days, RSE 15.99%)
    emba    <- 0.102;              label("Maximum fractional decrease in baseline aPTT over time, EM_BA (unitless)")        # Liesenfeld 2006 Table 2 (RSE 14.41%)
    emmx    <- 0.463;              label("Maximum fractional decrease in Emax over time, EM_MX (unitless)")                 # Liesenfeld 2006 Table 2 (RSE 12.68%)

    # =================================================================
    # Inter-individual variability (Liesenfeld 2006 Table 2). Per the
    # footnote IIV is on E0 and Emax (the time-varying values), not on
    # the initial-value parameters BAS0 / EMA0; encoding eta on the
    # typical-value lrbase / lemax with the (1 - EM_BA*T/(ET50+T))
    # factor applied multiplicatively afterwards is equivalent because
    # the time-decay factor is shared across subjects. Variances on the
    # internal log scale: omega^2 = log(1 + CV^2).
    #   E0   CV  8.7% -> omega^2 = log(1 + 0.087^2) = 0.007541
    #   Emax CV 19.9% -> omega^2 = log(1 + 0.199^2) = 0.03884
    #   EC50 CV 38.5% -> omega^2 = log(1 + 0.385^2) = 0.13821
    #   SLOP CV 15.2% -> omega^2 = log(1 + 0.152^2) = 0.02284
    # =================================================================
    etalrbase ~ 0.007541   # Liesenfeld 2006 Table 2 (IIV on E0, CV 8.7%)
    etalemax  ~ 0.03884    # Liesenfeld 2006 Table 2 (IIV on Emax, CV 19.9%)
    etalec50  ~ 0.13821    # Liesenfeld 2006 Table 2 (IIV on EC50, CV 38.5%)
    etalslope ~ 0.02284    # Liesenfeld 2006 Table 2 (IIV on SLOP, CV 15.2%)

    # =================================================================
    # Residual error. Liesenfeld 2006 Table 2 reports a proportional
    # error on aPTT with sigma = 7.55% CV (RSE 3.53%).
    # =================================================================
    propSd <- 0.0755; label("Proportional residual error on aPTT (fraction)")   # Liesenfeld 2006 Table 2 (RSE 3.53%)
  })

  model({
    # =================================================================
    # PK: Liesenfeld 2013 two-compartment dabigatran disposition. Doses
    # in mg, vc in L, so central/vc has units of mg/L = ug/mL; multiply
    # by 1000 to express Cc in ng/mL to match the PD source values.
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
    # PD aPTT. The model uses simulation time `t` (hours) as the proxy
    # for time-after-surgery T; this is consistent with the source
    # paper's convention of anchoring T to the first dose, which the
    # protocol placed 4-8 h after surgery. ET50 in the ini() block is
    # stored in hours; the day-form value reported in the paper is
    # noted in the parameter comment.
    # =================================================================
    rbase   <- exp(lrbase + etalrbase)
    emax_0  <- exp(lemax  + etalemax)
    ec50    <- exp(lec50  + etalec50)
    slope   <- exp(lslope + etalslope)
    et50    <- exp(let50)

    time_factor <- t / (et50 + t)
    rbase_t <- rbase  * (1 - emba * time_factor)   # Eq. 3
    emax_t  <- emax_0 * (1 - emmx * time_factor)   # Eq. 4

    aPTT <- rbase_t + emax_t * Cc / (ec50 + Cc) + slope * Cc   # Eq. 1
    aPTT ~ prop(propSd)
  })
}
