`Fiedler-Kelly_2019_fremanezumab` <- function() {
  description <- "Two-compartment population PK model for fremanezumab (anti-CGRP IgG2 delta-a/kappa mAb) with first-order SC absorption, absorption lag time, and route-specific central volume / residual error supporting both IV and SC administration in healthy adults and adults with chronic or episodic migraine (Fiedler-Kelly 2019)."
  reference <- paste(
    "Fiedler-Kelly JB, Cohen-Barak O, Morris DN, Yoon E, Yeo KR, Ludwig EA,",
    "Bauer R, Loupe P. Population pharmacokinetic modelling and simulation",
    "of fremanezumab in healthy subjects and patients with migraine.",
    "Br J Clin Pharmacol. 2019;85(12):2721-2733.",
    "doi:10.1111/bcp.14096 (PMID 31418911).",
    sep = " "
  )
  vignette <- "Fiedler-Kelly_2019_fremanezumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for allometric scaling on CL (exponent 1.05) and on the central volume of distribution Vc (exponent 1.53; same exponent applied to both Vc,IV and Vc,SC per Table 2 footnotes d and e), normalized to the population median of 71 kg. Treated as a baseline (time-invariant) covariate in the source. Vp, Q, F, and ALAG1 are FIXED in the source and do not carry an allometric weight effect.",
      source_name        = "WT"
    ),
    ROUTE_IV = list(
      description        = "Indicator for intravenous administration of fremanezumab",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (subcutaneous)",
      notes              = "Per-subject (or per-dose-record) dosing-route indicator: 1 = IV cohort, 0 = SC cohort. Selects the route-specific central volume of distribution (Vc,IV = 2.98 L FIXED for ROUTE_IV = 1; Vc,SC = 1.88 L for ROUTE_IV = 0) and the route-specific residual-error structure (proportional only with sigma^2 = 0.0467 for IV; combined additive sigma^2 = 0.204 (ug/mL)^2 plus proportional sigma^2 = 0.0531 for SC), as reported in Fiedler-Kelly 2019 Table 2. Distinct from the rxode2 cmt event column (cmt = central for IV doses, cmt = depot for SC doses); ROUTE_IV is the additional per-subject covariate that the V and residual-error switches need. ADA was evaluated as a covariate in the source paper (only 0.7% of samples were ADA-positive) and was NOT a statistically significant predictor of fremanezumab PK, so it is not encoded as a model covariate.",
      source_name        = "Route of administration (IV vs SC; not given a NONMEM column letter in the paper text)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 2546L,
    n_observations   = 13745L,
    n_studies        = 7L,
    study            = "Pooled phase 1, 2b, and 3 fremanezumab studies (LBR-101-011, TV48125-PK-10078, LBR-101-021, LBR-101-022, TV48125-CNS-30049, TV48125-CNS-30050, TV48125-CNS-30051).",
    age_range        = "18-71 years",
    age_median       = "43 years",
    weight_range     = "43.5-131.8 kg (median 70.8 kg; 5th-95th percentiles approximately 51-101 kg)",
    weight_median    = "70.8 kg (the covariate model uses 71 kg as the median-of-population reference)",
    sex_female_pct   = 86.1,
    race_ethnicity   = c(Caucasian = 79.9, Other = 20.1),
    disease_state    = "Healthy adults (n = 74) and adults with chronic migraine or episodic migraine (n = 2474).",
    dose_range       = "Phase 1: 225 / 675 / 900 mg single IV (1-h infusion) or single SC; phase 2b: 225-900 mg SC every 28 days for 3 months (one cohort: 675 mg loading then 225 mg); phase 3 pivotal: 225 mg SC every 28 days for 3 months (with or without a 675 mg loading dose) or 675 mg SC every 28 days for 3 months; phase 3 long-term safety: 225 mg SC monthly (with 675 mg loading) or 675 mg SC quarterly for 12 months.",
    regions          = "International, multi-center; phase 1 study TV48125-PK-10078 enrolled matched Japanese and Caucasian cohorts. The full geographic breakdown is not reproduced in the main paper text.",
    ada_positive_pct = 0.7,
    sampling         = "Phase 1: full PK profiles to 90 d (LBR-101-011) or 225 d (TV48125-PK-10078). Phase 2b and phase 3 pivotal: predose troughs on d 1, 29, 57, and a follow-up sample on d 85. Phase 3 LTS: predose troughs on d 1, 85, 169, 253, 337, and a follow-up sample on d 534, plus two non-trough samples per subject collected 3-10 d or 15-20 d after a dose. 13745 fremanezumab concentrations total (2436 from phase 1 and 2b, 11309 from phase 3); 1.9% of postdose samples were below the LLOQ of 250 ng/mL and were excluded from modeling.",
    reference_subject = "Median WT = 71 kg, IV or SC route selected per ROUTE_IV; the typical-value subject reported in Fiedler-Kelly 2019 Results (CL = 0.0902 L/d, Vc,SC = 1.88 L, F = 0.658).",
    notes            = "Population-level demographics from Fiedler-Kelly 2019 Section 'Exploratory data analysis' (median age 43 y, median WT 70.8 kg, 86.1% female, 79.9% Caucasian); per-study breakdowns are in Supporting Table S1 (not on disk). Approximately 20% of PK samples were collected in the presence of preventive migraine medications, 55% in the presence of acute medications, and 9% in the presence of analgesics; ADA-positive samples comprised only 0.7% of those collected. None of the covariates examined (age, albumin, renal function, sex, race, injection site, acute / analgesic / preventive medication use, ADA status, hepatic impairment up to moderate) was a statistically significant predictor of fremanezumab PK in the forward-selection / backward-elimination analysis; only body weight was retained in the final model."
  )

  ini({
    # Structural typical-value parameters at the reference subject (median WT
    # 71 kg) per Fiedler-Kelly 2019 Table 2. Units: time in day, dose in mg,
    # volume in L -> Cc = central / vc has units mg/L = ug/mL, matching the
    # published concentration units (LLOQ 250 ng/mL = 0.25 ug/mL).
    lcl     <- log(0.0902);            label("Central clearance CL (L/day) at 71 kg")                                 # Fiedler-Kelly 2019 Table 2: CL = 0.0902 L/d
    lvc_iv  <- fixed(log(2.98));       label("Central volume of distribution for IV administration Vc,IV (L) at 71 kg") # Fiedler-Kelly 2019 Table 2: Vc,iv = 2.98 L, FIXED
    lvc_sc  <- log(1.88);              label("Central volume of distribution for SC administration Vc,SC (L) at 71 kg") # Fiedler-Kelly 2019 Table 2: Vc,SC = 1.88 L
    lka     <- log(0.180);             label("First-order SC absorption rate ka (1/day)")                             # Fiedler-Kelly 2019 Table 2: ka = 0.180 1/d
    lq      <- fixed(log(0.262));      label("Inter-compartmental clearance Q (L/day)")                               # Fiedler-Kelly 2019 Table 2: Q = 0.262 L/d, FIXED
    lvp     <- fixed(log(1.72));       label("Peripheral volume of distribution Vp (L)")                              # Fiedler-Kelly 2019 Table 2: Vp = 1.72 L, FIXED
    lfdepot <- fixed(log(0.658));      label("SC absolute bioavailability F1 (fraction)")                             # Fiedler-Kelly 2019 Table 2: F1 = 0.658, FIXED
    ltlag   <- fixed(log(0.0803));     label("SC absorption lag time ALAG1 (day)")                                    # Fiedler-Kelly 2019 Table 2: ALAG1 = 0.0803 d, FIXED

    # Allometric weight effects, normalized to the median 71 kg. The Vc
    # exponent applies to both Vc,IV and Vc,SC per Table 2 footnotes d and e
    # (the same allometric exponent is reported once and applies to both
    # route-specific Vc typical values). Vp, Q, F, and ALAG1 are FIXED with
    # no allometric weight effect in the source.
    e_wt_cl <- 1.05;                   label("Power exponent of (WT/71) on CL (unitless)")                            # Fiedler-Kelly 2019 Table 2: CL allometric exponent for weight = 1.05
    e_wt_vc <- 1.53;                   label("Power exponent of (WT/71) on Vc (unitless; same exponent on Vc,IV and Vc,SC)") # Fiedler-Kelly 2019 Table 2: Vc allometric exponent for weight = 1.53

    # Inter-individual variability. Fiedler-Kelly 2019 Table 2 reports BSV as
    # %CV on log-normal parameters (paper Methods: 'BSV in parameters was
    # modelled using an exponential form'); convert to internal variance via
    # omega^2 = log(CV^2 + 1):
    #   CL CV 23.4% -> log(1 + 0.234^2) = 0.05334
    #   Vc CV 35.1% -> log(1 + 0.351^2) = 0.11618
    #   ka CV 59.0% -> log(1 + 0.590^2) = 0.29870
    # No off-diagonal Omega elements were estimated (Results: 'scatterplots of
    # individual estimates of random effect terms ... showed no apparent
    # trends ... no off-diagonal omega matrix elements were estimated'), so a
    # diagonal Omega is used. The Vc BSV was estimated only from subjects with
    # SC administration (IV-infusion or trough-only subjects had Vc fixed at
    # the typical value during the fit); ka BSV was similarly estimated only
    # from intensive-sampling subjects. The same eta is applied here regardless
    # of route, which is the typical simulation-time interpretation; the
    # vignette flags this as a deviation from the source-paper fitting choice.
    etalcl    ~ 0.05334
    etalvc_sc ~ 0.11618
    etalka    ~ 0.29870

    # Residual error. Fiedler-Kelly 2019 Table 2 reports the residual
    # variability magnitudes as sigma^2 (NONMEM variance scale); the squared
    # CVs reported in the Results ('Estimated RV ranged from 93.3%CV ... near
    # the LLOQ to 23.1% ... >10 ug/mL') confirm this interpretation for the
    # SC combined model:
    #   At Cc = 10 ug/mL  (SC): sqrt(0.0531 + 0.204 / 10^2)  = 0.235 -> 23.5% (paper 23.1%)
    #   At Cc =  0.5 ug/mL (SC): sqrt(0.0531 + 0.204 / 0.25) = 0.932 -> 93.2% (paper 93.3%)
    # nlmixr2's prop() / add() functions take SDs, so the values below are
    # sqrt(sigma^2):
    #   IV proportional: sqrt(0.0467) = 0.21610  (CV approximately 21.6%)
    #   SC proportional: sqrt(0.0531) = 0.23043  (CV approximately 23.0%)
    #   SC additive    : sqrt(0.204)  = 0.45166  (ug/mL)
    CcpropSdIv <- 0.21610;             label("Proportional residual SD, IV cohort (fraction; = sqrt(0.0467))")        # Fiedler-Kelly 2019 Table 2: Residual variability iv = 0.0467 (sigma^2)
    CcpropSdSc <- 0.23043;             label("Proportional residual SD, SC cohort (fraction; = sqrt(0.0531))")        # Fiedler-Kelly 2019 Table 2: Residual variability SC proportional component = 0.0531 (sigma^2)
    CcaddSdSc  <- 0.45166;             label("Additive residual SD, SC cohort (ug/mL; = sqrt(0.204))")                # Fiedler-Kelly 2019 Table 2: Residual variability SC additive component = 0.204 (sigma^2, (ug/mL)^2)
  })

  model({
    # Individual structural parameters. Route-specific Vc is selected by
    # ROUTE_IV (1 = IV, 0 = SC). The same allometric weight exponent (1.53)
    # is applied to both Vc,IV and Vc,SC per Fiedler-Kelly 2019 Table 2
    # footnotes d and e.
    cl  <- exp(lcl + etalcl) * (WT / 71)^e_wt_cl
    vc_iv_typ <-          exp(lvc_iv)               * (WT / 71)^e_wt_vc
    vc_sc_typ <-          exp(lvc_sc + etalvc_sc)   * (WT / 71)^e_wt_vc
    vc  <- vc_iv_typ * ROUTE_IV + vc_sc_typ * (1 - ROUTE_IV)
    q   <- exp(lq)
    vp  <- exp(lvp)
    ka  <- exp(lka + etalka)

    fdepot <- exp(lfdepot)
    tlag   <- exp(ltlag)

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Concentration in the central compartment (mg/L = ug/mL).
    Cc <- central / vc

    # Two-compartment PK with first-order SC absorption from the depot. IV
    # doses go directly to central (cmt = central) and SC doses go to depot
    # (cmt = depot) in the event dataset; bioavailability and lag time are
    # applied to the depot compartment, so they affect SC doses only.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                               k12 * central - k21 * peripheral1
    f(depot)          <- fdepot
    alag(depot)       <- tlag

    # Route-specific combined residual error. The IV cohort uses a purely
    # proportional error model (additive component = 0); the SC cohort uses
    # a combined additive + proportional model. ROUTE_IV switches the active
    # SDs (per-subject or per-record indicator).
    CcpropSd <- CcpropSdIv * ROUTE_IV + CcpropSdSc * (1 - ROUTE_IV)
    CcaddSd  <-                          CcaddSdSc * (1 - ROUTE_IV)
    Cc ~ add(CcaddSd) + prop(CcpropSd)
  })
}
