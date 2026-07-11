Morganroth_2015_moxifloxacin <- function() {
  description <- paste(
    "Population pharmacodynamic linear concentration-effect model for",
    "moxifloxacin-induced placebo- and baseline-corrected QTcF interval",
    "prolongation (DDQTcF) in 40 healthy adult male Japanese and 40 healthy",
    "adult male Caucasian volunteers following a single 400 mg oral dose in",
    "a thorough QT study (Morganroth 2015).",
    "The linear mixed-effects PD model of Table 4 has form",
    "DDQTcF = (alpha + rho * RACE_WHITE) + (beta + gamma * RACE_WHITE) * Cc,",
    "with alpha = 1.71 ms, beta = 2.58 ms per ug/mL (Japanese reference),",
    "rho = 2.58 ms (additive intercept shift in Caucasian; USA site),",
    "gamma = -0.24 ms per ug/mL (concentration-by-country interaction;",
    "Caucasian slope = 2.34 ms per ug/mL).",
    "The source publication does not fit a popPK model; the PK driver in",
    "this file is a typical-value 1-compartment oral approximation with",
    "CL/F = 8.47 L/h, V/F = 132.6 L, and ka = 1.7 /h derived from the",
    "pooled NCA summary statistics in Morganroth 2015 Table 1 (see vignette",
    "Errata). Per operator sidecar-001 (option C) the paper's",
    "BVN(0, Sigma) subject random effects on the PD intercept and slope",
    "are OMITTED because Sigma is not numerically reported in the paper,",
    "and a small placeholder additive residual SD of 1 ms is used to satisfy",
    "rxode2's residual-error requirement; downstream users who need",
    "VPC-style simulation must add their own IIV.",
    sep = " "
  )

  reference <- paste(
    "Morganroth J, Wang Y, Thorn M, Kumagai Y, Harris S, Stockbridge N,",
    "Kleiman R, Shah R. (2015).",
    "Moxifloxacin-induced QTc interval prolongations in healthy male",
    "Japanese and Caucasian volunteers: a direct comparison in a",
    "thorough QT study.",
    "British Journal of Clinical Pharmacology 80(3):446-459.",
    "doi:10.1111/bcp.12684.",
    "Accepted Article Published Online 22 May 2015.",
    sep = " "
  )

  vignette <- "Morganroth_2015_moxifloxacin"

  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline (kg).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used by the typical-value NCA-derived PK driver to reproduce the",
        "reported by-ethnicity difference in moxifloxacin exposure: Morganroth",
        "2015 Methods 'Sample size' states 'body weight is known to affect",
        "moxifloxacin concentration given the same dose (higher body weight is",
        "associated with lower moxifloxacin concentration)'. Linear WT scaling",
        "(exponent 1) on CL/F and V/F with reference WT 71 kg (pooled mean of",
        "the two cohorts: Japanese 65.9 kg, Caucasian 76.6 kg per Results",
        "'Study population and exposure'). WT is time-fixed per subject."
      ),
      source_name        = "BWT"
    ),
    RACE_WHITE = list(
      description        = "White (Caucasian) race indicator, 1 = White, 0 = non-White (Japanese in this cohort).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Japanese; Japan study site)",
      notes              = paste(
        "Morganroth 2015 Methods 'Statistical plan' Equation 1 uses country",
        "as a binary indicator: country = 0 for Japanese subjects (Kitasato",
        "University East Hospital, Sagamihara, Japan) and country = 1 for",
        "Caucasian subjects (SeaView Research, Miami, USA). In this study",
        "cohort ethnicity is perfectly confounded with study site, so",
        "RACE_WHITE = country directly (1 = Caucasian/USA site, 0 =",
        "Japanese/Japan site); no value inversion is required. The paper",
        "frames the effect as ethnicity-based (Japanese vs Caucasian QT",
        "sensitivity) rather than site-based, so the canonical race",
        "indicator is the appropriate encoding. When simulating a subject",
        "of Caucasian ethnicity enrolled at a non-USA site (or Japanese",
        "ethnicity enrolled at a non-Japan site), the same coding applies",
        "because the paper does not disentangle ethnicity from site."
      ),
      source_name        = "country (Japan = 0, USA = 1; RACE_WHITE = country directly)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 80L,
    n_studies      = 1L,
    age_range      = "18-45 years (Japanese mean 33.8 +/- 7.9; Caucasian mean 30.9 +/- 7.2)",
    weight_range   = "Japanese 65.9 +/- 8.9 kg; Caucasian 76.6 +/- 8.3 kg",
    sex_female_pct = 0,
    race_ethnicity = c(Japanese = 50, White = 50),
    disease_state  = paste(
      "Healthy adult male volunteers meeting the ICH E14 thorough-QT-study",
      "inclusion criteria: BMI 18-28 kg/m^2 at screening; no family history",
      "of QTc prolongation or unexplainable sudden death at less than 50",
      "years of age; resting supine heart rate 50-100 beats/min; resting",
      "supine systolic BP 90-140 mmHg; resting supine diastolic BP 50-90",
      "mmHg; no baseline ECG abnormality (QTc <= 450 ms, QRS <= 110 ms,",
      "PR <= 200 ms, no second or third degree heart block)."
    ),
    dose_range     = paste(
      "Single 400 mg oral moxifloxacin (site-specific tablet formulation)",
      "administered in the fasted state (no food or liquid except water for",
      "at least 10 h prior to dosing and at least 4 h thereafter) with a",
      "minimum 3-day washout between the placebo and moxifloxacin periods",
      "of the two-period crossover design."
    ),
    regions        = "Japan (Kitasato University East Hospital, Sagamihara) for Japanese subjects; USA (SeaView Research Inc., Miami, Florida) for Caucasian subjects.",
    notes          = paste(
      "Randomized, double-blind, two-period, crossover, ICH-E14-compliant",
      "thorough QT study (Morganroth 2015 Methods 'Study design').",
      "Half of subjects at each site received placebo followed by 400 mg",
      "moxifloxacin; the other half received moxifloxacin first. Plasma",
      "moxifloxacin was quantified by validated LC-MS/MS (LLOQ 0.001",
      "ug/mL). ECG time-matched sampling at 0, 0.25, 0.5, 1, 2, 3, 4, 6,",
      "8, 12, and 23.5 h on baseline day 0 and treatment day 1 of each",
      "period; triplicate 12-lead ECGs at each time point averaged to a",
      "single value. NCA summary statistics (Table 1): Japanese Cmax 3.27",
      "+/- 0.6 ug/mL, AUC(0, inf) 52.1 +/- 9.6 ug/mL/h, t1/2 11.7 +/- 1.2",
      "h; Caucasian Cmax 2.98 +/- 0.7 ug/mL, AUC(0, inf) 42.4 +/- 7.0",
      "ug/mL/h, t1/2 10.0 +/- 1.3 h; median tmax 2 h in both groups.",
      "The observed PK differences between ethnic groups were not",
      "statistically significant and the paper attributes them primarily",
      "to body-weight differences."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Typical-value PK driver (Morganroth 2015 Table 1 pooled NCA
    # summary statistics). The paper does NOT fit a popPK model; these
    # values are fixed for simulation purposes only and are flagged as
    # such in the vignette Errata. Pooled arithmetic means across the
    # two 40-subject cohorts (Japanese, Caucasian) give:
    #   CL/F  = 400 mg / ((52.1 + 42.4)/2 ug/mL h) = 8.47 L/h
    #   t1/2  = (11.7 + 10.0)/2 = 10.85 h
    #   kel   = ln(2) / 10.85   = 0.0639 /h
    #   V/F   = CL / kel         = 132.6 L
    # ka = 1.7 /h chosen to reproduce median tmax = 2 h with the
    # first-order oral absorption / 1-compartment structure and the
    # pooled kel above.
    # ------------------------------------------------------------------
    lka <- fixed(log(1.7))
    label("Absorption rate constant ka (1/h; NCA-derived, fixed)")
    # NCA-derived to reproduce Morganroth 2015 Table 1 median tmax = 2 h
    # given the pooled 1-compartment kel = 0.0639 /h. NOT a popPK fit --
    # see vignette Errata.

    lcl <- fixed(log(8.47 * 70 / 71))
    label("Apparent clearance CL/F (L/h, reference WT 70 kg; NCA-derived, fixed)")
    # Pooled Dose / AUC(0, inf) from Morganroth 2015 Table 1: 400 mg /
    # ((52.1 + 42.4)/2 ug/mL h) = 8.47 L/h at pooled reference WT 71 kg.
    # Rescaled to canonical reference WT 70 kg via linear scaling:
    # 8.47 * 70/71 = 8.35 L/h. NOT a popPK fit -- see vignette Errata.

    lvc <- fixed(log(132.6 * 70 / 71))
    label("Apparent central volume of distribution V/F (L, reference WT 70 kg; NCA-derived, fixed)")
    # Pooled V/F = (CL/F) / kel = 8.47 / (ln(2) / 10.85) = 132.6 L at
    # pooled reference WT 71 kg; rescaled to canonical reference WT 70
    # kg via linear scaling: 132.6 * 70/71 = 130.7 L. NOT a popPK fit
    # -- see vignette Errata.

    # ------------------------------------------------------------------
    # Published PD linear-model structural parameters (Morganroth 2015
    # Table 4, "Model parameters" and "Estimates of prediction lines by
    # ethnicity" blocks). Standard errors reported to the right of each
    # estimate in Table 4.
    # ------------------------------------------------------------------
    int_ddqtcf <- 1.71
    label("DDQTcF intercept for the Japanese reference cohort (ms)")
    # Morganroth 2015 Table 4 "Intercept" row: 1.71 ms (SE 1.29, P =
    # 0.1876). Reference cohort = Japanese (country = 0 in Eq. 1).

    slope_ddqtcf <- 2.58
    label("DDQTcF slope on moxifloxacin plasma concentration, Japanese reference (ms per ug/mL)")
    # Morganroth 2015 Table 4 "Moxifloxacin plasma concentration" row:
    # 2.58 ms per ug/mL (SE 0.62, P < 0.0001). This is also reported in
    # the "Estimates of prediction lines by ethnicity" block under
    # "Japan Slope vs plasma concentration" = 2.58.

    # ------------------------------------------------------------------
    # Covariate-effect coefficients (Morganroth 2015 Table 4 "Country"
    # and "Concentration-by-country interaction" rows). Encoded as
    # additive shifts on the intercept and slope respectively.
    # ------------------------------------------------------------------
    e_race_white_int_ddqtcf <- fixed(2.58)
    label("Additive DDQTcF intercept shift in Caucasian (USA) subjects (ms)")
    # Morganroth 2015 Table 4 "Country" row: 2.58 ms (SE 1.82, P =
    # 0.1606). Applied as e_race_white_int_ddqtcf * RACE_WHITE so
    # Caucasian intercept = 1.71 + 2.58 = 4.29 ms, matching the "USA
    # Intercept" value in Table 4 "Estimates of prediction lines by
    # ethnicity" block. Encoded fixed() because the covariate effect
    # value is a reported estimate of the paper's Eq. 1 fixed-effects
    # regression, not a re-fittable IIV or residual parameter.

    e_race_white_slope_ddqtcf <- fixed(-0.24)
    label("Additive DDQTcF slope shift in Caucasian (USA) subjects (ms per ug/mL)")
    # Morganroth 2015 Table 4 "Concentration-by-country interaction"
    # row: -0.24 ms per ug/mL (SE 0.89, P = 0.7853). Applied as
    # e_race_white_slope_ddqtcf * RACE_WHITE so Caucasian slope =
    # 2.58 - 0.24 = 2.34 ms per ug/mL, matching the "USA Slope vs
    # plasma concentration" value in Table 4 "Estimates of prediction
    # lines by ethnicity" block.

    # ------------------------------------------------------------------
    # Placeholder additive residual error (per operator sidecar-001
    # option C). The paper describes an additive residual e_ij ~
    # N(0, sigma^2) but does not numerically report sigma. A 1 ms
    # placeholder is used to satisfy rxode2's residual-error machinery
    # for the typical-value simulation; downstream users who need
    # realistic within-subject variability must supply their own value.
    # The Methods "Sample size" section quotes a historical bootstrap
    # study SD of 10.5 ms for DDQTcF (used only for the power
    # calculation, not fitted in the current paper).
    # ------------------------------------------------------------------
    addSd <- fixed(1)
    label("Additive residual error on DDQTcF (ms; placeholder, not from paper)")
  })

  model({
    # ================================================================
    # 1. Typical-value 1-compartment oral PK driver: linear body-weight
    #    scaling on CL/F and V/F (exponent 1) with reference WT 70 kg.
    #    Parameters are FIXED from pooled NCA summary statistics
    #    (Morganroth 2015 Table 1); this is NOT a popPK fit.
    # ================================================================
    ka <- exp(lka)
    cl <- exp(lcl) * (WT / 70)
    vc <- exp(lvc) * (WT / 70)
    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Moxifloxacin plasma concentration in ug/mL (dose in mg, vc in L
    # gives central/vc in mg/L = ug/mL, the paper's reporting unit).
    Cc <- central / vc

    # ================================================================
    # 2. PD linear concentration-effect model (Morganroth 2015 Eq. 1
    #    fixed-effect form; random-effects terms s_ij, d_ij, e_ij
    #    omitted per operator sidecar-001 option C):
    #      DDQTcF = (alpha + rho * RACE_WHITE)
    #             + (beta + gamma * RACE_WHITE) * Cc
    # ================================================================
    intercept_ddqtcf <- int_ddqtcf   + e_race_white_int_ddqtcf   * RACE_WHITE
    slope_ind        <- slope_ddqtcf + e_race_white_slope_ddqtcf * RACE_WHITE

    DDQTcF <- intercept_ddqtcf + slope_ind * Cc

    # ================================================================
    # 3. Additive residual error on DDQTcF (placeholder per sidecar-001
    #    option C).
    # ================================================================
    DDQTcF ~ add(addSd)
  })
}
