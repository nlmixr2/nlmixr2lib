Barbolosi_2017_radioiodine <- function() {
  description <- paste(
    "Semi-mechanistic three-compartment model of radioiodine (iodine-131)",
    "therapy response in metastatic papillary thyroid cancer (Barbolosi",
    "2017), reused without re-estimation by Fusella Giuntini 2026 for",
    "protocol optimisation (RAIR-Sim). Each RAI administration is a bolus",
    "of activity (GBq) into an effective-activity state that decays at",
    "a*ln(2) per month (the delayed, not the physical, effect of iodine-131);",
    "tumor cells grow exponentially with doubling time Td under treatment",
    "and are killed at rate rho*A(t); serum thyroglobulin is produced by",
    "tumor cells at lambda per cell and eliminated at ke. A two-class",
    "Monolix mixture on Td separates non-responders (fast growth, Td 9.8",
    "months, 27.5%) from responders (Td 66.6 months, 72.5%), supplied as",
    "the binary covariate MIX_FAST_GROW. Time in months."
  )
  reference <- paste(
    "Barbolosi D, Summer I, Meille C, Serre R, Kelly A, Zerdoud S, Bournaud C,",
    "Schvartz C, Toubeau M, Toubert ME, Keller I, Taieb D. Modeling",
    "therapeutic response to radioiodine in metastatic thyroid cancer: a",
    "proof-of-concept study for individualized medicine. Oncotarget.",
    "2017;8(24):39167-39176. doi:10.18632/oncotarget.16637.",
    "Model reused for protocol simulation in: Fusella Giuntini M, Voyant C,",
    "Taieb D, Barbolosi D. A computational framework for optimizing",
    "radioiodine therapy protocols in metastatic thyroid cancer. Sci Rep.",
    "2026;16:27179. doi:10.1038/s41598-026-56267-1 (MATLAB code:",
    "https://github.com/MarieFG49/CTD_SIMU_RAI)."
  )
  vignette <- "Barbolosi_2017_radioiodine"
  units <- list(time = "month", dosing = "GBq", concentration = "ng/mL")

  # The effective-activity state is the paper's delayed-effect abstraction of
  # the administered iodine-131, not its physical activity (half-life 8 days).
  paper_specific_compartments <- "activity"

  compartmentData <- list(
    activity = list(
      analyte = "iodine-131 (effective activity)",
      units = "GBq",
      specimen = "not applicable",
      verified = TRUE
    ),
    tumor = list(analyte = "metastatic thyroid cancer cells", units = "cells", specimen = "tumor", verified = TRUE),
    thyroglobulin = list(analyte = "thyroglobulin", units = "ng/mL", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    MIX_FAST_GROW = list(
      description = paste(
        "Per-subject latent mixture-class indicator for the Td mixture.",
        "1 = non-responder (fast tumor growth under RAI; Monolix category 1,",
        "typical Td 9.8 months); 0 = responder (Monolix category 2, typical",
        "Td 66.6 months)."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (responder; majority class, 72.5% of the source cohort)",
      notes = paste(
        "Not a measured covariate: Barbolosi 2017 assigned classes with the",
        "Monolix mixture algorithm (Table 1 'Mixture': pi_cat1 = 0.275,",
        "pi_cat2 = 0.725). For population simulation draw MIX_FAST_GROW ~",
        "Bernoulli(0.275). Fusella Giuntini 2026 simulates the two classes",
        "separately ('non-responder' Td = 9.8, 'responder' Td = 66.6)."
      ),
      source_name = "Monolix mixture category (cat = 1 / 2); MIX_FAST_GROW = as.integer(cat == 1)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 50L,
    n_studies = 1L,
    age_range = "8-55 years",
    age_mean = "29 years",
    disease_state = paste(
      "Metastatic papillary thyroid carcinoma with diffuse lung metastases",
      "(all M1; pT1/pT2 7, pT3/pT4 43; N1 47, Nx 3) after total",
      "thyroidectomy, with lymph-node dissection except in 3 patients"
    ),
    dose_range = paste(
      "Iodine-131 3.7-5.5 GBq per course (one 8-year-old received 1.85 GBq",
      "for the first course), mostly every 6-8 months; mean cumulative",
      "activity 22.2 GBq; at least 3 courses; endogenous TSH stimulation"
    ),
    regions = "France (eight institutions)",
    notes = paste(
      "Barbolosi 2017 Results 'Patient characteristics' and Materials and",
      "methods (inclusion criteria); stimulated Tg measured before each RAI",
      "course over at least 3 years of follow-up. Sex distribution not",
      "reported. Fusella Giuntini 2026 Appendix B repeats the same 50-patient",
      "cohort."
    )
  )

  ini({
    # Barbolosi 2017 Results: 'parameters {rho, a, lambda} had weak
    # variability ... These parameters were subsequently fixed'; Table 1
    # lists them without S.E. and with IIV 0 for rho and a.
    lkdecay <- fixed(log(0.0169)); label("Delayed effectiveness rate a of radioiodine (1/month); effective activity decays at a*ln(2)") # Barbolosi 2017 Table 1 'a' = 0.0169 (fixed); Fusella Giuntini 2026 Appendix A Table 1 'a' = 0.0169
    lkkill <- fixed(log(0.00407)); label("Efficiency rate rho of iodine on tumor cells (1/(GBq*month))") # Barbolosi 2017 Table 1 'rho' = 0.00407 (fixed); Fusella Giuntini 2026 Table 1 = 0.00407
    lksyn <- fixed(log(3.86e-9)); label("Thyroglobulin production rate lambda per tumor cell (ng/mL/month per cell)") # Barbolosi 2017 Table 1 'lambda' = 3.86e-009 (fixed); Fusella Giuntini 2026 Table 1 = 3.86e-9

    lkdeg <- log(0.319); label("Thyroglobulin elimination rate constant ke (1/month)") # Barbolosi 2017 Table 1 'ke' = 0.319 (S.E. 0.068, R.S.E. 21%)
    lrbase_tumor <- log(1.12e9); label("Initial tumor cell count N0 at the first RAI course (cells)") # Barbolosi 2017 Table 1 'N0' = 1.12e+009 (S.E. 1.2e+003)

    # Two-class mixture on Td (Barbolosi 2017 Table 1 'Group Parameters').
    # Consistency: the Monolix category effect beta_Td = 1.92 (Table 1) gives
    # 9.8 * exp(1.92) = 66.9 months, matching the printed 66.6.
    ltdouble_pop1 <- log(9.8); label("Tumor-cell doubling time under RAI, non-responder class MIX_FAST_GROW = 1 (month)") # Barbolosi 2017 Table 1 'Td cat=1' = 9.8 (S.E. 1.6, R.S.E. 16%)
    ltdouble_pop2 <- log(66.6); label("Tumor-cell doubling time under RAI, responder class MIX_FAST_GROW = 0 (month)") # Barbolosi 2017 Table 1 'Td cat=2' = 66.6 (S.E. 29, R.S.E. 44%)

    # Initial serum Tg: Barbolosi 2017 sets Tg(0) to each patient's first
    # stimulated Tg before the first course; Fusella Giuntini 2026 fixes it at
    # 100 ng/mL for every simulation. To reproduce the Barbolosi usage, supply
    # a per-subject data column `lrbase_thyroglobulin` = log(Tg0 in ng/mL),
    # which overrides this value in rxSolve().
    lrbase_thyroglobulin <- fixed(log(100)); label("Initial serum thyroglobulin Tg0 before the first RAI course (ng/mL)") # Fusella Giuntini 2026 'Practical implementation': 'the initial thyroglobulin concentration was fixed at 100 ng/mL'

    # IIV: Monolix log-normal omegas are SDs (Table 1 'Interindividual
    # Variability'); variances below are omega^2. rho and a have omega 0 and
    # carry no eta. N0 has omega 2.05e-007 (no S.E.), numerically zero, so it
    # carries no eta either: a 4e-14 variance makes rxode2 rewrite the omega
    # matrix in every simulation. The Td omega is printed without S.E.
    etalkdeg ~ 1.3456 # Barbolosi 2017 Table 1 IIV 'ke' omega = 1.16 (S.E. 0.23); 1.16^2
    etalksyn ~ 6.1009 # Barbolosi 2017 Table 1 IIV 'lambda' omega = 2.47 (S.E. 0.34); 2.47^2
    etaltdouble ~ fixed(0.09) # Barbolosi 2017 Table 1 IIV 'Td' omega = 0.3 (no S.E.); 0.3^2, shared by both classes

    propSd <- 0.372; label("Proportional residual error on thyroglobulin (fraction)") # Barbolosi 2017 Table 1 'b' = 0.372 (S.E. 0.031); Results: 'normal distribution with proportional variance'
  })
  model({
    kdecay <- exp(lkdecay)
    kkill <- exp(lkkill)
    ksyn <- exp(lksyn + etalksyn)
    kdeg <- exp(lkdeg + etalkdeg)
    rbase_tumor <- exp(lrbase_tumor)
    rbase_thyroglobulin <- exp(lrbase_thyroglobulin)

    # Mixture class selects the typical Td; the Td eta is shared.
    tdouble <- exp(ltdouble_pop1 * MIX_FAST_GROW + ltdouble_pop2 * (1 - MIX_FAST_GROW) + etaltdouble)
    kgrow <- log(2) / tdouble

    # RAI courses are bolus doses (GBq) into `activity`, which starts at 0.
    tumor(0) <- rbase_tumor
    thyroglobulin(0) <- rbase_thyroglobulin

    # Barbolosi 2017 equations for A, N and Tg; Fusella Giuntini 2026
    # Eq. 1a prints dA/dt = -a*A, but Barbolosi 2017 and the RAIR-Sim MATLAB
    # code (diff(A) == -a*log(2)*A) both carry the ln(2) factor.
    d/dt(activity) <- -kdecay * log(2) * activity
    d/dt(tumor) <- kgrow * tumor - kkill * activity * tumor
    d/dt(thyroglobulin) <- ksyn * tumor - kdeg * thyroglobulin

    thyroglobulin ~ prop(propSd)
  })
}
