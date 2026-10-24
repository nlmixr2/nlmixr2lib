Ferl_2018_mosunetuzumab_cyno_covariates <- function() {
  description <- "Preclinical (cynomolgus monkey). Two-compartment population PK model of the anti-CD20/CD3 T-cell-dependent bispecific antibody mosunetuzumab (BTCT4465A) with first-order subcutaneous absorption, a constant linear clearance, a decaying linear clearance and a decaying Michaelis-Menten elimination; final model with covariates (cell-line source on V1, CLd and the decaying clearance; baseline CD4+ count on CL1 and the Vmax decay constant; body weight on the clearance decay constant) (Ferl 2018)."
  reference <- "Ferl GZ, Reyes A, Sun LL, Cheu M, Oldendorp A, Ramanujan S, Stefanich EG. A Preclinical Population Pharmacokinetic Model for Anti-CD20/CD3 T-Cell-Dependent Bispecific Antibodies. Clin Transl Sci. 2018;11(3):296-304. doi:10.1111/cts.12535"
  vignette <- "Ferl_2018_mosunetuzumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    FORM_MOSU_CHO = list(
      description = "Drug-substance cell-line source indicator (1 = CHO-derived material, 0 = E. coli-derived material)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (E. coli-derived)",
      notes = "Categorical covariate on V1, CLd and CL2_0 (Table S3 reports separate E. coli and CHO estimates). Only two of the seven NHP cohorts received CHO-derived material, most of it in one study with n = 44, so the paper cautions that the effect may be confounded with between-study variability.",
      source_name = "cell line source (E. coli vs CHO)"
    ),
    CD4_ABS = list(
      description = "Baseline peripheral-blood CD4+ T-cell count",
      units = "cells/uL",
      type = "continuous",
      reference_category = NULL,
      notes = "Pre-dose baseline value. Linear effect per Eq. (6), p = p_pop * (1 + theta * (CD4_ABS - 2134)), on CL1 and on the Vmax decay constant lambda1; median 2134 cells/uL (Table S3 caption). Figure S8 plots the effect over roughly 1000-7000 cells/uL.",
      source_name = "baseline CD4+ count"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effect per Eq. (6), lambda2 = lambda2_pop * (1 + theta * (WT - 3.12)), on the decay constant of the time-varying linear clearance only; median 3.12 kg (Table S3 caption). With theta = -0.562 per kg, lambda2 reaches zero at WT = 4.90 kg and turns negative (an exponentially growing clearance) above it, so the model is not valid for heavier animals. The weight range of the cohort is not reported. CL and V values are absolute for the typical 3.12 kg animal and are not weight-scaled.",
      source_name = "BWT"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "mosunetuzumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mosunetuzumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "mosunetuzumab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "cynomolgus monkey (Chinese origin)",
    n_subjects = NA,
    n_studies = 7,
    weight_median = "3.12 kg",
    disease_state = "Healthy, purpose-bred, treatment-naive cynomolgus monkeys (Charles River Laboratories, Reno, NV).",
    dose_range = "Single IV slow bolus 0.001-1 mg/kg (3-4 animals per group); repeat dosing once weekly x 4 at 0.01-1 mg/kg as IV slow bolus or 1-h IV infusion, or 1 mg/kg SC (3-8 animals per group).",
    notes = "Same seven-study NHP dataset as Ferl_2018_mosunetuzumab_cyno. Median baseline CD4+ count 2134 cells/uL; median body weight 3.12 kg. Two cohorts received CHO-derived material (one study with n = 44); the rest received E. coli-derived material. The total number of animals is not reported. BTCT4465A is the Genentech development code for mosunetuzumab (RG7828)."
  )

  ini({
    # Ferl 2018 Table S3 (final model with covariates). CL1, CL2_0, CLd, V1 and V2 were normalised
    # by the 3.12 kg median body weight in a post hoc step (Table S3 caption); the values below
    # undo that for the typical 3.12 kg animal and convert mL to L. Vmax_0 is absolute (ug/day).
    lka <- log(1.26); label("SC absorption rate constant Ka (1/day)") # Table S3 Ka = 1.26 1/day
    logitfdepot <- 1.58; label("Logit of SC bioavailability, F = 1/(1 + exp(-Z)) (unitless)") # Table S3 Z = 1.58 (F = 0.83)
    lcl_exp_inf <- log(6.95 * 3.12 / 1000); label("Constant nonsaturable clearance CL1 at the median CD4+ count (L/day)") # Table S3 CL1 = 6.95 mL/day/kg x 3.12 kg
    lcl_exp_component_ecoli <- log(855 * 3.12 / 1000); label("Initial time-varying linear clearance CL2_0, E. coli material (L/day)") # Table S3 CL2_0,E.coli = 855 mL/day/kg x 3.12 kg
    lcl_exp_component_cho <- log(120 * 3.12 / 1000); label("Initial time-varying linear clearance CL2_0, CHO material (L/day)") # Table S3 CL2_0,CHO = 120 mL/day/kg x 3.12 kg
    lcl_exp_kdes <- log(35.2); label("Decay constant lambda2 of the time-varying linear clearance at 3.12 kg (1/day)") # Table S3 lambda2 = 35.2 1/day
    lq_ecoli <- log(34.3 * 3.12 / 1000); label("Distribution clearance CLd, E. coli material (L/day)") # Table S3 CLd,E.coli = 34.3 mL/day/kg x 3.12 kg
    lq_cho <- log(61.8 * 3.12 / 1000); label("Distribution clearance CLd, CHO material (L/day)") # Table S3 CLd,CHO = 61.8 mL/day/kg x 3.12 kg
    lvc_ecoli <- log(43.6 * 3.12 / 1000); label("Central volume V1, E. coli material (L)") # Table S3 V1,E.coli = 43.6 mL/kg x 3.12 kg
    lvc_cho <- log(108 * 3.12 / 1000); label("Central volume V1, CHO material (L)") # Table S3 V1,CHO = 108 mL/kg x 3.12 kg
    lvp <- log(40.4 * 3.12 / 1000); label("Peripheral volume V2 (L)") # Table S3 V2 = 40.4 mL/kg x 3.12 kg
    lvmax <- log(1760 / 1000); label("Initial Michaelis-Menten maximum elimination rate Vmax_0 (mg/day)") # Table S3 Vmax_0 = 1760 ug/day
    lvmax_exp_kdes <- log(0.167); label("Decay constant lambda1 of the time-varying Vmax at the median CD4+ count (1/day)") # Table S3 lambda1 = 0.167 1/day
    lkm <- log(31.1); label("Michaelis-Menten constant KM (ug/mL)") # Table S3 KM = 31.1 ug/mL

    # Linear covariate slopes, Eq. (6): p = p_pop * (1 + theta * (cov - median)).
    e_cd4_abs_cl_exp_inf <- 0.000314; label("Linear effect of baseline CD4+ count on CL1 (per cell/uL, centred at 2134)") # Table S3 theta_cov CD4,CL1 = 0.000314
    e_cd4_abs_vmax_exp_kdes <- 0.000371; label("Linear effect of baseline CD4+ count on lambda1 (per cell/uL, centred at 2134)") # Table S3 theta_cov CD4,lambda1 = 0.000371
    e_wt_cl_exp_kdes <- -0.562; label("Linear effect of body weight on lambda2 (per kg, centred at 3.12 kg)") # Table S3 theta_cov BWT,lambda2 = -0.562

    # Table S2 variance-covariance matrix, order CL2_0, Vmax_0, CL1, CLd, V1, V2. Although Table S2
    # is captioned as the base Model VI matrix, its diagonal reproduces the Table S3 %CVs through
    # Eq. (5) (20, 40, 59, 39, 39, 44 vs printed 20, 40, 59, 39, 38, 44) and not the Table 2 %CVs.
    # As printed to two decimals the matrix is slightly indefinite (smallest eigenvalue -0.00095;
    # the CL2_0-V1 correlation is -0.97), so the off-diagonals are scaled by 0.99 to make it
    # positive definite; the variances are unchanged.
    etalcl_exp_component + etalvmax + etalcl_exp_inf + etalq + etalvc + etalvp ~ c(
      0.15,
      0.99 * -0.1, 0.18,
      0.99 * -0.039, 0.99 * 0.031, 0.041,
      0.99 * -0.14, 0.99 * 0.16, 0.99 * 0.052, 0.3,
      0.99 * -0.14, 0.99 * 0.086, 0.99 * 0.042, 0.99 * 0.12, 0.14,
      0.99 * -0.075, 0.99 * 0.13, 0.99 * 0.042, 0.99 * 0.19, 0.99 * 0.06, 0.14
    ) # Table S2

    # Supplement section B.9: proportional error model; magnitude not reported.
    propSd <- fixed(0); label("Proportional residual error (fraction; not reported)")
  })
  model({
    ka <- exp(lka)
    fdepot <- expit(logitfdepot)

    cl_exp_inf <- exp(lcl_exp_inf + etalcl_exp_inf) * (1 + e_cd4_abs_cl_exp_inf * (CD4_ABS - 2134))
    cl_exp_component <- exp((1 - FORM_MOSU_CHO) * lcl_exp_component_ecoli + FORM_MOSU_CHO * lcl_exp_component_cho + etalcl_exp_component)
    cl_exp_kdes <- exp(lcl_exp_kdes) * (1 + e_wt_cl_exp_kdes * (WT - 3.12))
    q <- exp((1 - FORM_MOSU_CHO) * lq_ecoli + FORM_MOSU_CHO * lq_cho + etalq)
    vc <- exp((1 - FORM_MOSU_CHO) * lvc_ecoli + FORM_MOSU_CHO * lvc_cho + etalvc)
    vp <- exp(lvp + etalvp)
    vmax <- exp(lvmax + etalvmax)
    vmax_exp_kdes <- exp(lvmax_exp_kdes) * (1 + e_cd4_abs_vmax_exp_kdes * (CD4_ABS - 2134))
    km <- exp(lkm)

    # Eq. (1) and Figure 1, t = time since the first dose
    cl <- cl_exp_inf + cl_exp_component * exp(-cl_exp_kdes * time)
    vmax_t <- vmax * exp(-vmax_exp_kdes * time)

    Cc <- central / vc

    # Eqs. (1)-(3)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (cl + q) / vc * central + q / vp * peripheral1 -
      vmax_t * Cc / (Cc + km)
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1
    f(depot) <- fdepot

    Cc ~ prop(propSd)
  })
}
