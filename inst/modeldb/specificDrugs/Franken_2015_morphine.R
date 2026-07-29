Franken_2015_morphine <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for morphine and its two",
    "glucuronide metabolites (M3G, M6G) in 47 terminally ill adult palliative",
    "care patients (Franken 2015). Morphine: two-compartment disposition with",
    "three parallel first-order absorption routes (subcutaneous bolus,",
    "immediate-release oral liquid, controlled-release oral tablet) using",
    "route-specific fixed absorption rate constants; oral bioavailability F is",
    "estimated (SC F assumed 1). M3G and M6G are each one-compartment models",
    "fed by fixed-fraction transformation of morphine clearance (Fm1 = 0.55",
    "for M3G, Fm2 = 0.10 for M6G, both fixed from literature). Morphine",
    "clearance decreases exponentially as time-to-death (TTD, days) approaches",
    "zero. Metabolite clearance depends on estimated glomerular filtration",
    "rate (eGFR, MDRD four-variable formula) and serum albumin via shared",
    "power-form covariate exponents. Residual variability was reported as",
    "additive error on the log-transformed observation (LTBS)."
  )
  reference <- paste(
    "Franken LG, Masman AD, de Winter BCM, Koch BCP, Baar FPM,",
    "Tibboel D, van Gelder T, Mathot RAA.",
    "Pharmacokinetics of Morphine, Morphine-3-Glucuronide and",
    "Morphine-6-Glucuronide in Terminally Ill Adult Patients.",
    "Clin Pharmacokinet. 2016;55(6):697-709",
    "(published online 29 December 2015).",
    "doi:10.1007/s40262-015-0345-4.",
    sep = " "
  )
  vignette <- "Franken_2015_morphine"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ug/L"
  )

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Estimated glomerular filtration rate computed from the",
        "four-variable Modification of Diet in Renal Disease (MDRD)",
        "formula (age, sex, race, serum creatinine)."
      ),
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying renal-function covariate. Used with power scaling",
        "(CRCL / 96)^0.673 on both M3G and M6G clearance, normalised to",
        "the population median eGFR of 96 mL/min/1.73 m^2 reported in",
        "Table 1. Same coefficient applied to both metabolite clearances",
        "(paper Table 2 'Covariate effect on M3G and M6G clearance: eGFR').",
        "Source column 'eGFR' in the paper text and table; canonical",
        "register name CRCL covers MDRD-derived eGFR."
      ),
      source_name        = "eGFR (MDRD four-variable, paper Eq. 4)"
    ),
    ALB = list(
      description        = "Serum albumin concentration (g/L).",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying serum albumin. Used with power scaling",
        "(ALB / 26)^1.1 on both M3G and M6G clearance, normalised to",
        "the population median of 26 g/L (Table 1). Same coefficient",
        "applied to both metabolite clearances (paper Table 2 'Covariate",
        "effect on M3G and M6G clearance: Albumin'). g/L units; multiply",
        "g/dL values by 10."
      ),
      source_name        = "albumin (g/L)"
    ),
    TTD = list(
      description        = paste(
        "Time-to-death covariate (days remaining until the patient's",
        "recorded time of death). Required by Franken 2015 as a",
        "time-varying covariate on morphine clearance."
      ),
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Days from the observation time to the patient's recorded time",
        "of death (TTD >= 0). Used in a first-order exponential decay",
        "term on morphine CL: CL(TTD) = CL_pop - theta_D * exp(-theta_rate * TTD)",
        "with theta_D = 17.6 L/h and theta_rate = 0.13 /day (paper Eq. 3,",
        "Table 2). As TTD -> 0 (death), morphine CL drops by theta_D =",
        "17.6 L/h relative to its asymptotic value far from death.",
        "Available only retrospectively (time of death must be known);",
        "for prospective simulation, set TTD to a large value (> 50",
        "days) to recover the asymptotic CL when the time of death is",
        "unknown."
      ),
      source_name        = "TTD (paper Eq. 3 and Table 2)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 47,
    n_studies      = 1,
    age_range      = "43-93 years",
    age_median     = "71 years",
    weight_range   = "not reported",
    weight_median  = "not reported",
    sex_female_pct = 55.3,
    race_ethnicity = c(Caucasian = 95.7, AfroCaribbean = 4.3),
    disease_state  = paste(
      "Terminally ill adult palliative-care patients admitted to a",
      "Dutch hospice with prognosis of 2 days to 3 months. 95.7%",
      "had advanced malignancy (neoplasm) as the primary diagnosis;",
      "remainder had circulatory or respiratory disease. Sparse",
      "venous samples (median 2 per subject, range 1-10) drawn in",
      "both pre-terminal and terminal phases (terminal = bed-bound,",
      "semi-comatose, unable to take oral medication)."
    ),
    dose_range     = paste(
      "Morphine 15-540 mg/day administered per Dutch national",
      "palliative guidelines. Routes: oral immediate-release liquid,",
      "oral controlled-release tablet, or subcutaneous bolus / infusion.",
      "Two patients (4.2%) had concomitant codeine."
    ),
    regions        = "The Netherlands (Laurens Cadenza palliative care centre, Rotterdam)",
    notes          = paste(
      "Demographics from Franken 2015 Table 1. NONMEM 7.2 + PsN 3.7.6",
      "with FOCE-I and ADVAN5; data were log-transformed before fitting.",
      "Sampling: 152 plasma samples for morphine, M3G, and M6G; ~12% of",
      "concentrations below LOQ were discarded (M1 method). Median",
      "admission duration 33 days (range 7-457). Baseline blood",
      "chemistry medians: albumin 26 g/L (range 14-39), urea 7.2 mmol/L,",
      "bilirubin 8 umol/L, GGT 64 U/L, ALP 112 U/L, ALT 12 U/L, AST 32 U/L,",
      "CRP 67 U/L, creatinine 72 umol/L, eGFR (standard MDRD) 96",
      "mL/min/1.73 m^2 (range 27-239). Final model evaluated by 500-run",
      "bootstrap and normalised prediction distribution errors (NPDE)."
    )
  )

  ini({
    # ============================================================
    # Morphine absorption (route-specific Ka, all fixed from
    # literature per Methods 3.1: 'absorption constants could not be
    # estimated and were therefore fixed to known literature values
    # [26, 27]').
    # ============================================================
    lka_sc      <- fixed(log(10))
    label("Ka subcutaneous bolus injection (1/h, FIXED)")           # Methods 3.1: ka SC = 10 /h
    lka_oral_ir <- fixed(log(6))
    label("Ka oral immediate-release liquid (1/h, FIXED)")          # Methods 3.1: ka IR liquid = 6 /h
    lka_oral_cr <- fixed(log(0.8))
    label("Ka oral controlled-release tablet (1/h, FIXED)")         # Methods 3.1: ka CR tablet = 0.8 /h

    # Bioavailability of oral morphine; SC F is structurally fixed
    # at 1 per Methods 2.3.1.
    lfdepot <- log(0.30)
    label("Oral bioavailability F (log scale)")                     # Table 2 Final: F = 0.30 (RSE 13.6%)

    # ============================================================
    # Morphine structural disposition - Franken 2015 Table 2 Final
    # ============================================================
    lcl <- log(47.5)
    label("Morphine clearance, asymptote far from death (L/h)")     # Table 2 Final: CL = 47.5 L/h (RSE 11%); this is CL_pop at TTD -> infinity
    lvc <- log(190)
    label("Morphine central volume V1 (L)")                         # Table 2 Final: V1 = 190 L (RSE 28%)
    lq  <- log(76.1)
    label("Morphine inter-compartmental clearance Q (L/h)")         # Table 2 Final: Q = 76.1 L/h (RSE 35.7%)
    lvp <- log(243)
    label("Morphine peripheral volume V2 (L)")                      # Table 2 Final: V2 = 243 L (RSE 19%)

    # ============================================================
    # Metabolite (M3G, M6G) disposition - one-compartment each.
    # Reference values are typical-value clearances at population
    # median CRCL = 96 mL/min/1.73 m^2 and median ALB = 26 g/L.
    # ============================================================
    lcl_m3g <- log(1.44)
    label("M3G clearance at population median eGFR and albumin (L/h)") # Table 2 Final: CL_M3G = 1.44 L/h (RSE 4.8%)
    lvc_m3g <- log(8.02)
    label("M3G central volume (L)")                                    # Table 2 Final: V_M3G = 8.02 L (RSE 33.2%)
    lcl_m6g <- log(1.78)
    label("M6G clearance at population median eGFR and albumin (L/h)") # Table 2 Final: CL_M6G = 1.78 L/h (RSE 6.8%)
    lvc_m6g <- log(8.24)
    label("M6G central volume (L)")                                    # Table 2 Final: V_M6G = 8.24 L (RSE 30.7%)

    # Transformation ratios (Methods 2.3.1): fraction of morphine
    # clearance routed to each glucuronide. Fixed at literature
    # values because the study lacked mass-balance data to estimate
    # them independently.
    fm_m3g <- fixed(0.55)
    label("Fraction of morphine CL routed to M3G (FIXED)")          # Methods 2.3.1 / Table 2: Fm1 = 0.55 (fixed, literature)
    fm_m6g <- fixed(0.10)
    label("Fraction of morphine CL routed to M6G (FIXED)")          # Methods 2.3.1 / Table 2: Fm2 = 0.10 (fixed, literature)

    # ============================================================
    # TTD (time-to-death) covariate on morphine CL.
    # CL_morphine(TTD) = CL_pop - ttd_d * exp(-ttd_rate * TTD)
    # As TTD -> 0 (day of death), CL falls from 47.5 to 47.5 - 17.6
    # = 29.9 L/h, matching paper Section 3.3 simulations.
    # ============================================================
    ttd_d <- 17.6
    label("TTD: peak drop in morphine CL at TTD = 0 (L/h)")          # Table 2 Final: TTD_D = 17.6 (RSE 24.7%)
    ttd_rate <- 0.13
    label("TTD: first-order rate of CL drop approaching death (1/day)") # Table 2 Final: TTD_rate = 0.13 (RSE 32%)

    # ============================================================
    # eGFR (CRCL) and albumin (ALB) covariates on metabolite CL.
    # Same coefficient applied to both M3G and M6G clearance per
    # Table 2 row 'Covariate effect on M3G and M6G clearance'.
    # The metabolite-suffixed effect names below are duplicated
    # (m3g / m6g) with identical values to make the source-trace
    # explicit; the paper's single estimate is the source.
    # ============================================================
    e_crcl_cl_m3g <- 0.673
    label("eGFR power exponent on M3G clearance (unitless)")        # Table 2 Final: eGFR = 0.673 (RSE 16.8%); shared with M6G CL
    e_crcl_cl_m6g <- 0.673
    label("eGFR power exponent on M6G clearance (unitless)")        # Table 2 Final: eGFR = 0.673 (shared single estimate); identical to e_crcl_cl_m3g
    e_alb_cl_m3g  <- 1.1
    label("Albumin power exponent on M3G clearance (unitless)")     # Table 2 Final: Albumin = 1.1 (RSE 23.3%); shared with M6G CL
    e_alb_cl_m6g  <- 1.1
    label("Albumin power exponent on M6G clearance (unitless)")     # Table 2 Final: Albumin = 1.1 (shared single estimate); identical to e_alb_cl_m3g

    # ============================================================
    # IIV - log-normal; omega^2 = log(1 + CV^2). Reported CV% from
    # Table 2 Final 'Between-subject variability (%)' block.
    # ============================================================
    etalfdepot ~ log(1 + 0.378^2)
    # Table 2 Final IIV: F = 37.8% CV (RSE 38.3%, shrinkage 9.5%) -> omega^2 = log(1 + 0.378^2)
    etalcl ~ log(1 + 0.534^2)
    # Table 2 Final IIV: morphine CL = 53.4% CV (RSE 30.1%, shrinkage 13.3%) -> omega^2 = log(1 + 0.534^2)

    # M3G CL and M6G CL share BSV with correlation fixed to 1
    # (Section 3.1: 'The correlation between BSV of M3G and M6G
    # clearance was high and fixed to unity'). Encoded as a 2x2
    # block with off-diagonal = sqrt(var_m3g * var_m6g).
    etalcl_m3g + etalcl_m6g ~ c(
      log(1 + 0.293^2),
      sqrt(log(1 + 0.293^2) * log(1 + 0.343^2)),
      log(1 + 0.343^2)
    )
    # Table 2 Final IIV: M3G CL = 29.3% CV, M6G CL = 34.3% CV; rho fixed to 1.

    # M3G V1 and M6G V1 share BSV with correlation fixed to 1
    # (Section 3.1: 'A similar approach was used for BSV on the
    # volumes of distribution of M3G and M6G').
    etalvc_m3g + etalvc_m6g ~ c(
      log(1 + 1.517^2),
      sqrt(log(1 + 1.517^2) * log(1 + 1.430^2)),
      log(1 + 1.430^2)
    )
    # Table 2 Final IIV: M3G V1 = 151.7% CV, M6G V1 = 143.0% CV; rho fixed to 1.

    # ============================================================
    # Residual variability - 'additive error on the log scale'
    # (Methods 2.3.1). The reported values are interpreted as the
    # NONMEM $SIGMA variance estimates on the log scale; the SD
    # supplied to lnorm() is the square root.
    # ============================================================
    expSd <- sqrt(0.432)
    label("Morphine log-normal residual SD on log scale")           # Table 2 Final: morphine residual = 0.432 (variance on log scale, $SIGMA convention) -> SD = sqrt(0.432)
    expSd_m3g <- sqrt(0.246)
    label("M3G log-normal residual SD on log scale")                # Table 2 Final: M3G residual = 0.246 (variance on log scale) -> SD = sqrt(0.246)
    expSd_m6g <- sqrt(0.265)
    label("M6G log-normal residual SD on log scale")                # Table 2 Final: M6G residual = 0.265 (variance on log scale) -> SD = sqrt(0.265)
  })

  model({
    # ------------------------------------------------------------
    # Reference covariate values (Table 1 population medians).
    # ------------------------------------------------------------
    egfr_ref <- 96   # Table 1: median standard MDRD eGFR = 96 mL/min/1.73 m^2
    alb_ref  <- 26   # Table 1: median albumin = 26 g/L

    # ------------------------------------------------------------
    # Route-specific absorption rate constants (all fixed).
    # ------------------------------------------------------------
    ka_sc      <- exp(lka_sc)
    ka_oral_ir <- exp(lka_oral_ir)
    ka_oral_cr <- exp(lka_oral_cr)

    # ------------------------------------------------------------
    # Oral bioavailability with IIV. SC F is structurally fixed at 1
    # (the rxode2 default when f(depot_sc) is not assigned).
    # ------------------------------------------------------------
    f_oral <- exp(lfdepot + etalfdepot)

    # ------------------------------------------------------------
    # Morphine time-varying clearance via TTD.
    # CL_pop(TTD) = exp(lcl) - ttd_d * exp(-ttd_rate * TTD)
    # Individual CL applies log-normal IIV multiplicatively to the
    # TTD-adjusted population-typical CL.
    # ------------------------------------------------------------
    cl_pop_ttd <- exp(lcl) - ttd_d * exp(-ttd_rate * TTD)
    cl <- cl_pop_ttd * exp(etalcl)
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # ------------------------------------------------------------
    # Metabolite clearances and volumes. Shared eGFR/albumin
    # covariate factors applied multiplicatively (power form).
    # ------------------------------------------------------------
    crcl_factor_m3g <- (CRCL / egfr_ref)^e_crcl_cl_m3g
    crcl_factor_m6g <- (CRCL / egfr_ref)^e_crcl_cl_m6g
    alb_factor_m3g  <- (ALB  / alb_ref )^e_alb_cl_m3g
    alb_factor_m6g  <- (ALB  / alb_ref )^e_alb_cl_m6g
    cl_m3g <- exp(lcl_m3g + etalcl_m3g) * crcl_factor_m3g * alb_factor_m3g
    cl_m6g <- exp(lcl_m6g + etalcl_m6g) * crcl_factor_m6g * alb_factor_m6g
    vc_m3g <- exp(lvc_m3g + etalvc_m3g)
    vc_m6g <- exp(lvc_m6g + etalvc_m6g)

    # ------------------------------------------------------------
    # Morphine -> metabolite formation flux. Cleared fraction of
    # morphine CL is partitioned: Fm1 to M3G, Fm2 to M6G, and
    # (1 - Fm1 - Fm2) to remaining renal / other elimination
    # (Fig. 2 in the paper). Working in mg of morphine-equivalents
    # in every compartment; M3G and M6G plasma concentrations were
    # adjusted to their morphine equivalents using molecular weights
    # (Methods 2.3.1), so 1 mg of morphine cleared via the M3G
    # pathway becomes 1 mg of morphine-equivalent M3G in the
    # metabolite compartment.
    # ------------------------------------------------------------
    cl_form_m3g <- fm_m3g * cl
    cl_form_m6g <- fm_m6g * cl
    cl_other    <- (1 - fm_m3g - fm_m6g) * cl

    # ------------------------------------------------------------
    # ODE system. Three parallel depots feed the morphine central
    # compartment with route-specific first-order absorption. SC
    # bioavailability is the rxode2 default (1); oral routes carry
    # f_oral via f(depotN) <- f_oral.
    # ------------------------------------------------------------
    d/dt(depot)  <- -ka_sc      * depot
    d/dt(depot2) <- -ka_oral_ir * depot2
    d/dt(depot3) <- -ka_oral_cr * depot3
    d/dt(central) <-
      ka_sc      * depot      +
      ka_oral_ir * depot2     +
      ka_oral_cr * depot3     +
      q  * peripheral1 / vp   -
      q  * central     / vc   -
      (cl_form_m3g + cl_form_m6g + cl_other) * central / vc
    d/dt(peripheral1) <- q * central / vc - q * peripheral1 / vp

    d/dt(central_m3g) <- cl_form_m3g * central / vc - cl_m3g * central_m3g / vc_m3g
    d/dt(central_m6g) <- cl_form_m6g * central / vc - cl_m6g * central_m6g / vc_m6g

    # ------------------------------------------------------------
    # Oral bioavailability assignment (SC depot stays at default 1).
    # ------------------------------------------------------------
    f(depot2) <- f_oral
    f(depot3) <- f_oral

    # ------------------------------------------------------------
    # Observations. Internal states are mg of morphine-equivalents;
    # plasma concentrations are reported in ug/L (paper's assay
    # range 2-500 ug/L). 1 mg/L = 1000 ug/L.
    # ------------------------------------------------------------
    Cc     <- central     / vc     * 1000
    Cc_m3g <- central_m3g / vc_m3g * 1000
    Cc_m6g <- central_m6g / vc_m6g * 1000

    Cc     ~ lnorm(expSd)
    Cc_m3g ~ lnorm(expSd_m3g)
    Cc_m6g ~ lnorm(expSd_m6g)
  })
}
