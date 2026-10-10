Srimani_2022_ixazomib_platelets <- function() {
  description <- "Semi-mechanistic platelet-count model for oral ixazomib plus lenalidomide-dexamethasone (LenDex) in relapsed/refractory multiple myeloma from the phase III TOURMALINE-MM1 trial (Srimani 2022). Ixazomib plasma concentrations come from the three-compartment population PK model of Gupta 2017 (fixed). An adapted Friberg chain (proliferating pool, two maturation transit compartments, circulating platelets) without rebound feedback carries two additive first-order losses from the proliferating pool: an ixazomib effect linear in the plasma concentration plus a scaled cumulative AUC, and a LenDex effect linear in the effect rate of a hypothetical kinetic-pharmacodynamic (K-PD) compartment fed by the lenalidomide doses. The individual baseline is a typical baseline plus a linear effect of the observed baseline platelet count and a random effect whose variance equals the residual variance (Dansirikul 2008)."
  reference <- paste(
    "Srimani JK, Diderichsen PM, Hanley MJ, Venkatakrishnan K, Labotka R, Gupta N.",
    "Population pharmacokinetic/pharmacodynamic joint modeling of ixazomib efficacy",
    "and safety using data from the pivotal phase III TOURMALINE-MM1 study in",
    "multiple myeloma patients. CPT Pharmacometrics Syst Pharmacol.",
    "2022;11(8):1085-1099. doi:10.1002/psp4.12815.",
    "Ixazomib PK layer: Gupta N, Diderichsen PM, Hanley MJ, et al. Clin Pharmacokinet.",
    "2017;56(11):1355-1368. doi:10.1007/s40262-017-0526-4",
    "(also available as modellib('Gupta_2017_ixazomib'))."
  )
  vignette <- "Srimani_2022_ixazomib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL", platelets = "10^9/L")

  compartmentData <- list(
    depot = list(analyte = "ixazomib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ixazomib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ixazomib", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "ixazomib", units = "mg", specimen = "tissue", verified = TRUE),
    auc_central = list(analyte = "ixazomib", units = "ng*h/mL", specimen = "not applicable", verified = TRUE),
    depot_kpd = list(
      analyte = "lenalidomide (hypothetical K-PD amount)",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    prol = list(
      analyte = "platelet precursors (proliferating pool)",
      units = "(fraction of baseline)",
      specimen = "tissue",
      verified = TRUE
    ),
    transit1 = list(
      analyte = "platelet precursors (transit 1)",
      units = "(fraction of baseline)",
      specimen = "tissue",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "platelet precursors (transit 2)",
      units = "(fraction of baseline)",
      specimen = "tissue",
      verified = TRUE
    ),
    circ = list(
      analyte = "circulating platelets",
      units = "(fraction of baseline)",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area at baseline.",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only by the Gupta 2017 ixazomib PK layer: power covariate on the second peripheral volume, reference 1.87 m^2, exponent 2.06.",
      source_name = "BSA"
    ),
    PLT_BASE = list(
      description = "Observed baseline platelet count (BL_i,obs).",
      units = "10^9/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed. Linear covariate on the individual baseline, centred on the safety-dataset median of 197 x 10^9/L (Srimani 2022 Equation 19; Supplementary Table 4). The modelled platelet count PLT is the time-varying output, so the observed baseline needs its own column.",
      source_name = "BL_i,obs"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 720L,
    n_studies = 1L,
    age_range = "30-91 years (median 66)",
    age_median = "66 years",
    sex_female_pct = 43.3,
    race_ethnicity = c(White = 85.0, Black = 1.8, Asian = 8.9, `Not reported` = 2.6, Other = 1.7),
    disease_state = "Relapsed and/or refractory multiple myeloma after 1-3 prior lines of therapy.",
    dose_range = "Ixazomib 4 mg or matching placebo orally on days 1, 8 and 15 of 28-day cycles, both arms with lenalidomide 25 mg (10 mg for reduced creatinine clearance) on days 1-21 and dexamethasone 40 mg on days 1, 8, 15 and 22.",
    regions = "Global (TOURMALINE-MM1, C16010).",
    baseline_labs = "Baseline platelet count 197 x 10^9/L (35.0-666); creatinine clearance 77.9 mL/min (20.2-231); hemoglobin 115 g/L (68.0-170). Medians [range] from Srimani 2022 Supplementary Table 4.",
    notes = "Safety population: all 720 treated patients (361 ixazomib, 359 placebo). Platelets were measured through the first six treatment cycles. Demographics from Srimani 2022 Supplementary Table 4."
  )

  ini({
    # Ixazomib PK -- Gupta 2017 Table 3, fixed (Srimani 2022 Methods).
    lka <- fixed(log(0.34)); label("Ixazomib first-order absorption rate constant (1/h)") # Gupta 2017 Table 3: Ka = 0.34/h
    lcl <- fixed(log(1.86)); label("Ixazomib clearance (L/h)") # Gupta 2017 Table 3: CL = 1.86 L/h
    lvc <- fixed(log(13.7)); label("Ixazomib central volume (L)") # Gupta 2017 Table 3: V2 = 13.7 L
    lfdepot <- fixed(log(0.58)); label("Ixazomib oral bioavailability (fraction)") # Gupta 2017 Table 3: F = 58%
    lq <- fixed(log(5.18)); label("Ixazomib intercompartmental clearance to peripheral1 (L/h)") # Gupta 2017 Table 3: Q3 = 5.18 L/h
    lvp <- fixed(log(309)); label("Ixazomib first peripheral volume (L)") # Gupta 2017 Table 3: V3 = 309 L
    lq2 <- fixed(log(26.1)); label("Ixazomib intercompartmental clearance to peripheral2 (L/h)") # Gupta 2017 Table 3: Q4 = 26.1 L/h
    lvp2 <- fixed(log(205)); label("Ixazomib second peripheral volume (L)") # Gupta 2017 Table 3: V4 = 205 L
    ltlag <- fixed(log(13 / 60)); label("Ixazomib absorption lag time (h)") # Gupta 2017 Table 3: TLAG = 13 min
    e_bsa_vp2 <- fixed(2.06); label("Power exponent of BSA on the second peripheral volume (unitless)") # Gupta 2017 Table 3: V4[BSA] = 2.06
    etalcl + etalfdepot ~ fixed(c(0.17697, 0.22550, 0.42726)) # Gupta 2017 Table 3: CL 44%CV, F 73%CV, rho 0.82
    etalvp2 ~ fixed(0.48492) # Gupta 2017 Table 3: V4 79%CV

    # Platelet PK/PD model -- Srimani 2022 Table 3. NONMEM time unit is the
    # hour: exp(-3.20) * 168 = 6.85/week against the printed 6.88/week.
    lkprp <- -3.20; label("Log platelet maturation rate constant kprp (log 1/h)") # Table 3 kprp = -3.20 (6.88/week)
    lkin <- -3.21; label("Log zero-order production rate kIn of the normalised proliferating pool, also the circulating-platelet loss rate (log 1/h)") # Table 3 kIn = -3.21 (6.78/week)
    bl <- 203; label("Typical baseline platelet count (10^9/L)") # Table 3 BL = 203
    e_plt_base_bl <- 0.837; label("Linear coefficient W of the centred observed baseline on the individual baseline (unitless)") # Table 3 W = 0.837
    slp_ixa <- 0.000859; label("Ixazomib effect slope on the precursor loss rate (1/h per ng/mL)") # Table 3 slpIXA = 0.000859 (0.144/week/(ng/ml))
    k_ixa <- 0.0000818; label("Scaling of the cumulative ixazomib AUC to a concentration in the ixazomib effect (1/h)") # Table 3 kIXA = 0.0000818 (0.0137/week)
    lkel <- log(0.0483); label("Log elimination rate constant of the lenalidomide K-PD compartment (log 1/h)") # Table 3 kLEN = 0.0483 (8.11/week), reported untransformed
    slp_len <- 0.0134; label("LenDex effect slope on the precursor loss rate per unit K-PD effect rate (1/h per mg/h)") # Table 3 slpLEN = 0.0134 (2.24/week/(len conc.))

    # Individual-baseline random effect: eta ~ N(0, 1) scaled by the residual
    # SD at the typical baseline (Equation 19).
    etabl ~ fixed(1) # Equation 19: eta_i has mean 0 and variance 1

    propSd_PLT <- 0.160; label("Proportional residual error on platelet count (fraction)") # Table 3 Prop. Error = -0.160 (sign of an SD theta is immaterial)
    addSd_PLT <- 28.1; label("Additive residual error on platelet count (10^9/L)") # Table 3 Add. Error = 28.1
  })

  model({
    # Ixazomib PK (Gupta 2017)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    q2 <- exp(lq2)
    vp2 <- exp(lvp2 + etalvp2) * (BSA / 1.87)^e_bsa_vp2
    fdepot <- exp(lfdepot + etalfdepot)
    tlag <- exp(ltlag)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (cl / vc) * central -
      (q / vc) * central + (q / vp) * peripheral1 -
      (q2 / vc) * central + (q2 / vp2) * peripheral2
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1
    d/dt(peripheral2) <- (q2 / vc) * central - (q2 / vp2) * peripheral2
    f(depot) <- fdepot
    alag(depot) <- tlag

    # mg / L * 1000 = ng/mL; cumulative AUC in ng*h/mL
    Cc <- 1000 * central / vc
    d/dt(auc_central) <- Cc

    # Lenalidomide K-PD: lenalidomide doses (mg) go into depot_kpd. The effect
    # is driven by the K-PD effect rate kel * depot_kpd (Jacqmin 2007); driving
    # it by the amount itself would push the placebo-arm platelet count to
    # about 13% of baseline within a week, against the ~25% within-cycle dip
    # of Supplementary Figure S14.
    kel <- exp(lkel)
    d/dt(depot_kpd) <- -kel * depot_kpd
    ceff_len <- kel * depot_kpd

    # Drug effects (Equations 20-21), additive first-order losses from prol
    e_ixa <- slp_ixa * (Cc + k_ixa * auc_central)
    e_len <- slp_len * ceff_len

    # Friberg chain without feedback (Equations 15-18), normalised so the
    # drug-free steady state has circ = 1. The paper lists no separate kout;
    # normalisation to baseline with production kIn implies kout = kIn.
    kprp <- exp(lkprp)
    kin <- exp(lkin)
    kout <- kin
    prol(0) <- kin / kprp
    transit1(0) <- kin / kprp
    transit2(0) <- kin / kprp
    circ(0) <- 1
    d/dt(prol) <- kin - (kprp + e_len + e_ixa) * prol
    d/dt(transit1) <- kprp * (prol - transit1)
    d/dt(transit2) <- kprp * (transit1 - transit2)
    d/dt(circ) <- kprp * transit2 - kout * circ

    # Individual baseline (Equation 19); 197 is the safety-dataset median.
    bl_i <- bl + (PLT_BASE - 197) * e_plt_base_bl +
      sqrt((propSd_PLT * bl)^2 + addSd_PLT^2) * etabl
    PLT <- bl_i * circ

    PLT ~ add(addSd_PLT) + prop(propSd_PLT)
  })
}
