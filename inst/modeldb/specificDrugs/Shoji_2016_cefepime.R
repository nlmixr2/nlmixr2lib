Shoji_2016_cefepime <- function() {
  description <- "Two-compartment IV population PK model for cefepime in 91 neonates, infants, and children (Shoji 2016); body-weight allometric scaling (fixed exponents 0.75 on CL and Q, 1.0 on Vss), nonlinear postmenstrual-age maturation on CL, a power effect of serum creatinine on CL, and a power effect of gestational age on Vss. Central volume of distribution enters as a fixed fraction of steady-state volume (V/Vss = 0.460)."
  reference <- paste(
    "Shoji K, Bradley JS, Reed MD, van den Anker JN, Domonoske C,",
    "Capparelli EV. Population pharmacokinetic assessment and",
    "pharmacodynamic implications of pediatric cefepime dosing for",
    "susceptible-dose-dependent organisms.",
    "Antimicrob Agents Chemother. 2016;60(4):2150-2156.",
    "doi:10.1128/AAC.02592-15.",
    sep = " "
  )
  vignette <- "Shoji_2016_cefepime"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Allometric scaling on CL and Q (exponent 0.75, fixed)",
        "and on Vss (exponent 1, fixed) per Shoji 2016 Methods 'Pharmacokinetic",
        "analysis' paragraph ('TVCL was scaled allometrically by the subject",
        "weight (weight^0.75), and TVVss was also scaled by subject weight",
        "(weight^1.0) before evaluation of other covariates'). No reference",
        "weight is centred -- the scaling enters as wt^exponent absolute."
      ),
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age in weeks / 4.35 + postnatal months)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Drives the nonlinear maturation factor on CL via",
        "[-0.09 + 1.09 * (1 - exp(-0.00958 * PMA_weeks))] (Shoji 2016",
        "Table 3). The source paper reports PMA in weeks; this model",
        "converts the canonical PAGE (months) back to weeks via",
        "pma_wk = PAGE * 4.35 so the Table-3 maturation constants",
        "(theta_6 = 1.09, theta_7 = 0.00958 per Table 5) apply unchanged.",
        "The maturation function is non-physical for PMA < ~9 weeks",
        "(where the factor turns negative); minimum observed PMA in the",
        "model-development cohort was GA 22.1 weeks (Table 1)."
      ),
      source_name        = "PMA"
    ),
    GA = list(
      description        = "Gestational age at birth",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Power-form effect on Vss:",
        "(GA / 30)^theta_9 with theta_9 = -0.548 (Shoji 2016 Table 3 / 5).",
        "Reference 30 weeks GA is the cohort median for infants under",
        "2 months. Shoji 2016 Methods states that GA was collected only",
        "for patients younger than 2 months of age; for older children",
        "the source paper does not document the imputation rule. Users",
        "must supply GA for every subject (set GA = 40 weeks for term",
        "older children to apply the formula at the term-birth reference)."
      ),
      source_name        = "GA"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Power-form effect on CL: (CREAT / 0.6)^theta_8",
        "with theta_8 = -0.392 (Shoji 2016 Table 3 / 5). Reference 0.6 mg/dL",
        "is the cohort median (Table 1). Source paper column SCR."
      ),
      source_name        = "SCR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 91L,
    n_studies      = 2L,
    age_range      = "Postnatal age 0.03-197.30 months (median 0.99, IQR 0.23-11.19); oldest patient 16 years",
    age_median     = "0.99 months",
    weight_range   = "0.58-75.00 kg (median 3.10, IQR 1.44-8.28)",
    weight_median  = "3.10 kg",
    sex_female_pct = 34.1,
    race_ethnicity = c(Caucasian = 40.7, `African American` = 34.1, Asian = 3.3, Hispanic = 22.0),
    disease_state  = "Neonates, infants, and children with suspected or proven bacterial infection receiving IV cefepime; data pooled from two previously published pediatric cefepime PK studies (Shoji 2016 refs 22 and 23)",
    dose_range     = "Varied across the pooled studies; typical dose-interval estimates assumed q12h with mean actual doses 48.8 +/- 4.7 mg/kg (all patients), 48.2 +/- 5.9 mg/kg (>= 30 days), 49.9 +/- 0.2 mg/kg (term neonates), 49.2 +/- 3.3 mg/kg (preterm neonates). 30-minute IV infusions in the source studies; intramuscular and BLQ samples excluded",
    regions        = "United States (single-center cohorts pooled from two prior pediatric cefepime PK studies, refs 22 and 23 of Shoji 2016)",
    ga_range       = "22.10-42.29 weeks (median 29, IQR 26.05-33) for infants under 2 months of age; 78.2% of those infants were preterm (GA < 36 weeks)",
    pna_range      = "0-7 days (25 subjects, 27.5%); 8-29 days (19, 20.9%); 30 days-1 year (26, 28.6%); >= 1 year (21, 23.1%)",
    creat_range    = "0.10-1.50 mg/dL (median 0.60, IQR 0.40-0.85)",
    n_observations = 664L,
    notes          = paste(
      "Baseline demographics per Shoji 2016 Table 1. 91 subjects with 664",
      "cefepime plasma concentrations (725 collected; 48 IM samples, 4 BQL,",
      "and 9 transcription / sampling / assay errors excluded). Final model",
      "fit in NONMEM 7.2 with FOCE-I; bootstrap n = 1000 with 86.9% success",
      "rate (Table 5). Covariates retained in the final model were PMA",
      "(nonlinear) and SCR on CL, and GA on Vss; sex was screened but not",
      "retained. Race was recorded (Caucasian, African American, Asian,",
      "Hispanic) but not tested as a covariate."
    )
  )

  ini({
    # Structural PK -- Shoji 2016 Table 3 final-model expressions and
    # Table 5 bootstrap estimates (rows labeled theta_1 .. theta_9).
    lcl <- log(0.395); label("Typical CL after maturation (L/h per kg^0.75)")              # Shoji 2016 Table 5: theta_1 = 0.395 +/- 0.049 (bootstrap 0.396 [0.315-0.550])
    lvss <- log(0.406); label("Typical Vss at WT = 1 kg (L/kg)")                            # Shoji 2016 Table 5: theta_2 = 0.406 +/- 0.023 (bootstrap 0.404 [0.361-0.451])
    lq <- log(0.575); label("Typical intercompartmental clearance (L/h per kg^0.75)")      # Shoji 2016 Table 5: theta_4 = 0.575 +/- 0.049 (bootstrap 0.580 [0.488-2.570])

    # V (central) / Vss ratio (estimated, Shoji 2016 Table 5 theta_3).
    # V_central = ratio * Vss; V_peripheral = (1 - ratio) * Vss.
    vc_vss_ratio <- 0.460; label("Ratio of central volume to steady-state volume (unitless)") # Shoji 2016 Table 5: theta_3 = 0.460 +/- 0.020 (bootstrap 0.452 [0.002-0.487])

    # Maturation function on CL (Shoji 2016 Table 3 equation).
    # CL_typical = theta_1 * [mat_intercept + mat_amplitude * (1 - exp(-mat_rate * PMA_weeks))] * wt^0.75 * (SCR/0.6)^theta_8
    # mat_intercept = -0.09 is hard-coded in the Table 3 equation; not listed as a
    # separate numbered theta in Table 5 (mathematically equals 1 - mat_amplitude
    # so that the maturation factor approaches 1 at PMA -> infinity, consistent
    # with theta_1 representing TVCL "after maturation"). Encoded as fixed because
    # the paper does not report an uncertainty estimate for it.
    mat_intercept <- fixed(-0.09); label("Maturation function intercept (unitless, fixed)") # Shoji 2016 Table 3 expression (hard-coded constant)
    mat_amplitude <- 1.09;         label("Maturation function amplitude (unitless)")        # Shoji 2016 Table 5: theta_6 = 1.09 +/- 0.087 (bootstrap 1.09 [0.945-1.320])
    mat_rate      <- 0.00958;      label("Maturation rate constant on PMA (1/week)")        # Shoji 2016 Table 3 expression (Table 5 theta_7 reports 0.010 +/- 0.003; Table 3 uses the more precise 0.00958)

    # Covariate exponents (Shoji 2016 Table 5).
    e_creat_cl  <- -0.392; label("Power exponent on (CREAT / 0.6) for CL (unitless)")               # Shoji 2016 Table 5: theta_8 = -0.392 +/- 0.101 (bootstrap -0.396 [-0.620 to -0.191])
    e_ga_vc_vp  <- -0.548; label("Shared power exponent on (GA / 30) for Vc and Vp (unitless)")     # Shoji 2016 Table 5: theta_9 = -0.548 +/- 0.221 (bootstrap -0.555 [-1.01 to -0.106]); source paper applies the exponent to Vss, which decomposes here into vc = ratio * Vss and vp = (1 - ratio) * Vss so the same exponent operates on both compartments

    # Allometric exponents (Shoji 2016 Methods: "TVCL was scaled allometrically
    # by the subject weight (weight^0.75), and TVVss was also scaled by subject
    # weight (weight^1.0) before evaluation of other covariates."). theta_4
    # (typical Q) is reported in Table 5 as L/h per kg^0.75, so Q also carries
    # the 0.75 allometric exponent. The Vss WT exponent applies equally to Vc
    # and Vp because both are linear shares of Vss.
    e_wt_cl_q  <- fixed(0.75); label("Shared allometric exponent on CL and Q (unitless, fixed)") # Shoji 2016 Methods (allometric scaling paragraph)
    e_wt_vc_vp <- fixed(1.00); label("Shared allometric exponent on Vc and Vp (unitless, fixed)") # Shoji 2016 Methods (allometric scaling paragraph)

    # Reference covariate values for the centred power terms (Shoji 2016
    # Table 1 cohort medians and Table 3 reference values).
    creat_ref <- 0.6; label("Reference serum creatinine (mg/dL, cohort median)")
    ga_ref    <- 30;  label("Reference gestational age (weeks)")

    # Inter-individual variability (Shoji 2016 Table 3 and Table 5 final-model CV%;
    # log-normal etas, so omega^2 = log(1 + CV^2)).
    # CL  31.8% CV -> log(1 + 0.318^2) = log(1.101124) = 0.09633
    # Vss 22.2% CV -> log(1 + 0.222^2) = log(1.049284) = 0.04811
    etalcl  ~ 0.09633 # Shoji 2016 Table 5: omega_CL 31.8% CV (bootstrap 30.6 [25.1-35.8])
    etalvss ~ 0.04811 # Shoji 2016 Table 5: omega_Vss 22.2% CV (bootstrap 21.1 [14.1-28.7])

    # Combined residual error -- Shoji 2016 Table 3 "Residual variability (%) = 66.3"
    # and Table 5 theta_5 = 0.705 ("fraction of additive and proportional error").
    # Per sidecar request-001 / response-001 (operator chose Option B), encoded
    # as a standard combined error model with sigma = 66.3% as the proportional
    # CV and theta_5 = 0.705 as the additive SD in mg/L (the source paper
    # reports plasma concentrations in ug/mL = mg/L). This matches the
    # convention used for similar antibiotic popPK extractions in this library
    # (e.g. Shi_2018_ceftazidime.R, Delattre_2010_amikacin.R) where two
    # residual-error numbers are reported without explicit units.
    propSd <- 0.663; label("Proportional residual error (fraction)") # Shoji 2016 Table 5: % sigma = 66.30 +/- 37.30 (bootstrap 66.1 [48.6-98.8])
    addSd  <- 0.705; label("Additive residual error (mg/L)")         # Shoji 2016 Table 5: theta_5 = 0.705 +/- 0.039 (bootstrap 0.705 [0.598-0.785])
  })

  model({
    # ----- Derived covariate terms -----
    # Convert canonical PAGE (months) back to source-paper PMA (weeks) so the
    # Shoji 2016 Table 3 maturation constants (theta_6, theta_7) apply unchanged.
    pma_wk <- PAGE * 4.35

    # Nonlinear postmenstrual-age maturation factor on CL (Shoji 2016 Table 3):
    #   F_PMA = mat_intercept + mat_amplitude * (1 - exp(-mat_rate * PMA_weeks))
    # Non-physical (negative) for PMA_weeks < about 9 weeks; minimum observed
    # PMA in the model-development cohort was GA 22.1 weeks (Table 1).
    fpma <- mat_intercept + mat_amplitude * (1 - exp(-mat_rate * pma_wk))

    # Power effect of serum creatinine on CL (Shoji 2016 Table 3):
    #   F_CREAT = (CREAT / 0.6)^theta_8
    creat_factor <- (CREAT / creat_ref)^e_creat_cl

    # Power effect of gestational age on Vss (Shoji 2016 Table 3):
    #   F_GA = (GA / 30)^theta_9
    ga_factor <- (GA / ga_ref)^e_ga_vc_vp

    # ----- Individual PK parameters -----
    cl  <- exp(lcl + etalcl) * fpma * WT^e_wt_cl_q * creat_factor
    vss <- exp(lvss + etalvss) * WT^e_wt_vc_vp * ga_factor
    q   <- exp(lq) * WT^e_wt_cl_q

    # Split Vss into central and peripheral volumes via the estimated ratio.
    vc <- vc_vss_ratio * vss
    vp <- (1 - vc_vss_ratio) * vss

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- ODE system (two-compartment IV) -----
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # ----- Observation -----
    # Dose in mg, vc in L -> central / vc has units mg/L (equivalent to ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
