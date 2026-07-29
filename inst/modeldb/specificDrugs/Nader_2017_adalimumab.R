Nader_2017_adalimumab <- function() {
  description <- "One-compartment population PK model with first-order subcutaneous absorption and linear elimination for adalimumab in adults with moderate-to-severe hidradenitis suppurativa, with baseline body weight and baseline C-reactive protein power effects on apparent clearance, a time-delayed anti-adalimumab-antibody effect on apparent clearance, and baseline body weight on apparent volume of distribution (Nader 2017)."
  reference <- "Nader A, Beck D, Noertersheuser P, Williams D, Mostafa N. Population Pharmacokinetics and Immunogenicity of Adalimumab in Adult Patients with Moderate-to-Severe Hidradenitis Suppurativa. Clin Pharmacokinet. 2017;56(9):1091-1102. doi:10.1007/s40262-016-0502-4"
  vignette <- "Nader_2017_adalimumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on both CL/F and V/F, each normalized as (WT/94 kg)^exponent per Nader 2017 Eqs. (5) and (6). Reference weight 94 kg is the pooled cohort median (Table 1). Baseline body weight, time-fixed per subject.",
      source_name        = "WTKG"
    ),
    CRP = list(
      description        = "Baseline C-reactive protein",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline (pre-treatment) CRP measured by standard assay, as a power effect on CL/F normalized as (CRP/7.9 mg/L)^exponent per Nader 2017 Eq. (5). Reference value 7.9 mg/L is the pooled cohort median (Table 1). Time-fixed per subject.",
      source_name        = "CRP"
    ),
    ADA_POS = list(
      description        = "Anti-adalimumab antibody (AAA) positivity",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (AAA-negative; typical patient)",
      notes              = "Time-delayed multiplicative power effect on CL/F, active only 21 days (3 weeks) or more after the first adalimumab dose in AAA-positive subjects, per Nader 2017 Eq. (5) definition of COV_1 = AAA_eff. In the source data a sample was classified AAA-positive if the AAA concentration in undiluted serum was greater than 20 ng/mL AND the sample was collected within 30 days of an adalimumab dose (Methods Section 2.3). The model column represents a subject-level 'will develop AAA' indicator to be used together with the model-time gate. Model time origin (t = 0) coincides with the first adalimumab dose (subcutaneous loading). Source paper labels this covariate 'AAA'; renamed to canonical ADA_POS per inst/references/covariate-columns.md.",
      source_name        = "AAA"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 600L,
    n_studies        = 3L,
    age_range        = "18-67 years",
    age_median       = "35 years",
    weight_range     = "43-221 kg",
    weight_median    = "94 kg",
    sex_female_pct   = 66,
    race_ethnicity   = "Not reported in the population PK Table 1 summary. The pooled phase II/III trials enrolled predominantly White adults with moderate-to-severe HS.",
    disease_state    = "Moderate-to-severe hidradenitis suppurativa (Hurley stage II 53%, stage III 43%) with inadequate response to a trial of oral antibiotics; naive to anti-tumor-necrosis-factor-alpha treatment.",
    dose_range       = "Loading 160 mg SC at week 0 and 80 mg at week 2 (phase III), then 40 mg SC weekly through week 12 and thereafter, or 40 mg SC every other week (phase II open-label period).",
    regions          = "Multiregional (phase II NCT00918255; phase III PIONEER-I NCT01468207 and PIONEER-II NCT01468233).",
    crp_median       = "7.9 mg/L (range 0.1-189 mg/L)",
    aaa_positive_pct = 6.5,
    notes            = "Pooled analysis of one phase II study (n = 143) and two phase III studies PIONEER-I (n = 295) and PIONEER-II (n = 162). Baseline demographics from Nader 2017 Table 1. Of 787 randomized patients, 187 were excluded from the population PK dataset (175 randomized to placebo or discontinued before receiving adalimumab; 12 without measurable adalimumab concentrations above the LLOQ). Serum adalimumab measured by validated ELISA; assay range 3.125-50.0 ng/mL in diluted serum, LLOQ 3.125 ng/mL. AAA measured by validated ELISA; sample classified AAA-positive if AAA > 20 ng/mL in undiluted serum within 30 days of an adalimumab dose."
  )

  ini({
    # Structural PK parameters -- final-model population estimates for the reference
    # patient (94 kg, baseline CRP 7.9 mg/L, AAA-negative). Nader 2017 Table 2.
    lka <- log(0.195); label("First-order SC absorption rate constant ka (1/day)")   # Nader 2017 Table 2, ka row (0.195 1/day, 9.90% RSE)
    lcl <- log(0.667); label("Apparent clearance CL/F (L/day)")                      # Nader 2017 Table 2, CL/F row (0.667 L/day, 2.65% RSE)
    lvc <- log(13.5);  label("Apparent central volume V/F (L)")                      # Nader 2017 Table 2, Vd/F row (13.5 L, 3.40% RSE)

    # Covariate effects on CL/F and V/F -- Nader 2017 Table 2, applied per Eqs. (5)
    # and (6) with reference values WTKG_median = 94 kg and CRP_median = 7.9 mg/L.
    e_ada_cl <- 6.76;  label("Multiplier on CL/F for AAA-positive subjects after day 21 (unitless)")  # Nader 2017 Table 2, AAA on CL/F (theta_4 in Eq. 5)
    e_crp_cl <- 0.173; label("Power exponent of CRP/7.9 mg/L on CL/F (unitless)")                     # Nader 2017 Table 2, CRP on CL/F (theta_5 in Eq. 5)
    e_wt_cl  <- 0.888; label("Power exponent of WT/94 kg on CL/F (unitless)")                         # Nader 2017 Table 2, Weight on CL/F (theta_6 in Eq. 5)
    e_wt_vc  <- 0.707; label("Power exponent of WT/94 kg on V/F (unitless)")                          # Nader 2017 Table 2, Weight on V/F (theta_7 in Eq. 6)

    # Inter-individual variability -- Nader 2017 states the final model included
    # correlated exponential etas on CL/F and V/F but does not publish the CL/F-V/F
    # covariance value. Encoded here as independent (diagonal) etas; the missing
    # off-diagonal is documented in the vignette Errata.
    # NONMEM $OMEGA reported as variance (omega^2): sqrt(0.346) = 0.588 -> 58.8% CV
    # on CL/F, sqrt(0.091) = 0.302 -> 30.2% CV on V/F, matching the "% CV" column
    # in Table 2.
    etalcl ~ 0.346  # Nader 2017 Table 2: IIV on CL/F 0.346 (58.8% CV)
    etalvc ~ 0.091  # Nader 2017 Table 2: IIV on V/F  0.091 (30.2% CV)

    # Residual error -- combined additive + proportional per Methods Eq. (2).
    # NONMEM $SIGMA reported as variance: proportional sigma^2 = 0.052 -> SD 0.228;
    # additive sigma^2 = 2.06 (ug/mL)^2 -> SD 1.435 ug/mL.
    propSd <- sqrt(0.052); label("Proportional residual error (fraction)")   # Nader 2017 Table 2: proportional residual variance 0.052
    addSd  <- sqrt(2.06);  label("Additive residual error (ug/mL)")          # Nader 2017 Table 2: additive residual variance 2.06 (ug/mL)^2
  })

  model({
    # AAA effect gate: the e_ada_cl multiplier is active only for AAA-positive
    # subjects starting 21 days (3 weeks) after the first adalimumab dose, per
    # Nader 2017 Eq. (5) definition of COV_1 = AAA_eff. Model time origin (t = 0)
    # coincides with the first adalimumab dose (SC loading). The (t > 21) gate
    # evaluates to 1 after day 21 and 0 otherwise, so aaa_active is 1 only when
    # ADA_POS = 1 AND t > 21 days.
    aaa_active <- ADA_POS * (t > 21)

    # Individual PK parameters -- Nader 2017 Eqs. (5) and (6):
    #   CL/F = theta_1 * theta_4^aaa_active * (CRP/7.9)^theta_5 * (WT/94)^theta_6 * exp(eta_1)
    #   V /F = theta_3 * (WT/94)^theta_7 * exp(eta_2)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
          e_ada_cl^aaa_active *
          (CRP / 7.9)^e_crp_cl *
          (WT / 94)^e_wt_cl
    vc <- exp(lvc + etalvc) *
          (WT / 94)^e_wt_vc

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Observation: dose in mg, V/F in L -> mg/L, which is equivalent to ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
