Galluppi_2021_ulotaront <- function() {
  description <- "Two-compartment population PK model with first-order oral absorption for ulotaront (SEP-363856), a trace amine-associated receptor 1 (TAAR1) agonist with 5-HT1A agonist activity in phase III development for schizophrenia. Pooled analysis of nine studies (seven phase I, one phase II acute, one 6-month open-label extension) in 404 adult subjects (99 healthy volunteers and 305 patients with schizophrenia). Body weight was estimated as a power-form covariate on the clearance parameters (CL/F, Q/F) and the volume parameters (Vc/F, Vp/F); disease status, sex, race (Asian vs non-Asian), and age were retained as full-model covariates on CL/F only. IIV on CL/F, Vc/F, ka, Vp/F is modelled as a full 4x4 correlated BLOCK. Residual error is proportional-only per AIC/BIC (Galluppi 2021)."
  reference   <- "Galluppi GR, Polhamus DG, Fisher JM, Hopkins SC, Koblan KS. Population pharmacokinetic analysis of ulotaront in subjects with schizophrenia. CPT Pharmacometrics Syst Pharmacol. 2021;10(10):1245-1254. doi:10.1002/psp4.12692"
  vignette    <- "Galluppi_2021_ulotaront"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "ulotaront", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "ulotaront", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "ulotaront", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at study entry. Full-model power exponents estimated on the clearance parameters (CL/F, Q/F: 0.821) and the volume parameters (Vc/F, Vp/F: 0.610). Reference weight 70 kg (rounded standard; paper does not state the reference explicitly, but the base-model allometric relationship was fixed at 0.75 on clearances and 1 on volumes -- a standard 70 kg-anchored allometric parameterisation). Analysis-set weight range 45.2-135.9 kg (mean [SD] = 77.7 [15.7] kg); Galluppi 2021 Results 'Population pharmacokinetic analysis dataset'.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age at study entry",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at study entry. Full-model power exponent on CL/F = -0.154 (95% CI -0.322, 0.0147); 95% CI includes zero and Galluppi 2021 Discussion concludes no clinically meaningful age effect. Reference age 35 years (rounded to the analysis-set mean of 33.3 years; paper does not state the reference explicitly). Analysis-set age range 18-55 years (mean [SD] = 33.3 [8.7] years).",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Full-model multiplicative effect on CL/F = 0.938 (95% CI 0.843, 1.04); 95% CI includes 1.0 and the paper concludes no clinically meaningful sex effect. Analysis set: 118 women (29.2%) / 286 men.",
      source_name        = "SEX"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator (1 = Asian, 0 = other)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian)",
      notes              = "Full-model multiplicative effect on CL/F = 0.987 (95% CI 0.874, 1.12); 95% CI includes 1.0 and the paper concludes no clinically meaningful race effect. Over 80% of Asian subjects in the analysis set were from the Japanese studies (Galluppi 2021 Results 'Population pharmacokinetic analysis dataset'); the paper's 'Asian vs non-Asian' contrast should be interpreted primarily as Japanese vs non-Asian. Analysis-set race distribution: 53.7% White, 31.4% Black, 10.9% Asian, 3.9% Other/Mixed.",
      source_name        = "ASIAN"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator: 1 = healthy adult volunteer, 0 = adult patient with schizophrenia.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (schizophrenia patient)",
      notes              = "Time-fixed per subject. Galluppi 2021 Table 1 encodes the effect with a PATIENT flag (1 = schizophrenia patient) and reports Patient CL = 0.809 (95% CI 0.720, 0.908) with typical values referenced to the healthy-volunteer state. The canonical DIS_HEALTHY convention uses 0 = patient as reference; the structural typical lcl below is shifted to the patient state and the covariate coefficient is negated so the model is mathematically identical: at DIS_HEALTHY = 0 (patient reference), CL/F = 32.5 * 0.809 = 26.29 L/h; at DIS_HEALTHY = 1, exp(+log(1/0.809)) restores the paper's HV-typical 32.5 L/h. Analysis set: 99 healthy volunteers (24.5%) / 305 patients with schizophrenia (75.5%). Disease status is the only clinically meaningful non-weight covariate in the full model (Galluppi 2021 Discussion).",
      source_name        = "PATIENT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 404L,
    n_studies      = 9L,
    n_observations = 4149L,
    age_range      = "18-55 years",
    age_median     = "not reported (mean [SD] = 33.3 [8.7] years)",
    weight_range   = "45.2-135.9 kg",
    weight_median  = "not reported (mean [SD] = 77.7 [15.7] kg)",
    sex_female_pct = 29.2,
    race_ethnicity = c(White = 53.7, Black = 31.4, Asian = 10.9, `Other/Mixed` = 3.9),
    disease_state  = "Adult subjects: healthy adult volunteers (n = 99) and adult patients with schizophrenia (n = 305) meeting DSM-IV-TR or DSM-5 criteria with baseline CGI-S <= 4 and PANSS total <= 80.",
    dose_range     = "5-150 mg once-daily oral (capsule or tablet), single and multiple dose (up to 7 days in phase I; up to 12 weeks in phase II acute; 6-month open-label extension at 25/50/75 mg flexible dosing).",
    regions        = "United States, global multi-regional (phase II), Japan (Studies DA801002 and DA801004).",
    notes          = "Pooled analysis of seven phase I studies (SEP-361-101, SEP-361-103, SEP-361-105, SEP-361-106, SEP-361-111, DA801002, SEP-361-1004), one phase II acute schizophrenia study (SEP-361-201), and one 6-month open-label extension (SEP-361-202). Two bioanalytical assays: LLQ 0.02 ng/mL for the earlier phase I studies (SEP-361-101/103/105/106/111) and 0.25 ng/mL for the Japanese studies (DA801002, DA801004) and phase II. 9.4% of active-drug samples were below the LLQ, mostly at later time points (>48 h post last dose); BLQ observations were not included in the model fit. Baseline demographics from Galluppi 2021 Results 'Population pharmacokinetic analysis dataset' and supplement Table S1."
  )

  ini({
    # ----------------------------------------------------------------------
    # Structural parameters -- Galluppi 2021 Table 1 (Full model estimates).
    # Paper reports typical values with the healthy-volunteer state as the
    # implicit PATIENT reference. The canonical DIS_HEALTHY register uses
    # 0 = patient as reference; the structural typical value lcl is shifted
    # to the patient state so that:
    #   * at DIS_HEALTHY = 0, WT = 70, RACE_ASIAN = 0, SEXF = 0, AGE = 35:
    #       CL/F = exp(lcl) = 32.5 * 0.809 = 26.29 L/h
    #   * the e_dis_healthy_cl coefficient in model() applies exp(+log(1/0.809))
    #     to recover the paper's HV-typical 32.5 L/h at DIS_HEALTHY = 1.
    # No covariate effects are reported on Vc/F, Vp/F, Q/F, or ka beyond
    # weight, so lvc / lvp / lq / lka equal the paper's values directly.
    # ----------------------------------------------------------------------
    lka <- log(0.966)          ; label("First-order oral absorption rate ka (1/h)")                                              # Galluppi 2021 Table 1: ka = 0.966 1/h
    lcl <- log(32.5 * 0.809)   ; label("Apparent clearance CL/F (L/h); patient reference at WT = 70 kg, AGE = 35 y")             # Galluppi 2021 Table 1: CL/F = 32.5 L/h (HV-typical) shifted to patient reference via Patient CL = 0.809
    lvc <- log(232)            ; label("Apparent central volume Vc/F (L); typical at WT = 70 kg")                                # Galluppi 2021 Table 1: Vc/F = 232 L
    lq  <- log(0.790)          ; label("Apparent intercompartmental clearance Q/F (L/h); typical at WT = 70 kg")                 # Galluppi 2021 Table 1: Q/F = 0.790 L/h
    lvp <- log(19.3)           ; label("Apparent peripheral volume Vp/F (L); typical at WT = 70 kg")                             # Galluppi 2021 Table 1: Vp/F = 19.3 L

    # Weight-effect exponents (estimated, not fixed). Applied uniformly to
    # both clearances (CL/F, Q/F) and both volumes (Vc/F, Vp/F) per
    # Galluppi 2021 Results 'Model development results'.
    e_wt_cl_q  <- 0.821        ; label("Weight power exponent shared on CL/F and Q/F (unitless)")                                # Galluppi 2021 Table 1: Weight CL = 0.821
    e_wt_vc_vp <- 0.610        ; label("Weight power exponent shared on Vc/F and Vp/F (unitless)")                               # Galluppi 2021 Table 1: Weight V = 0.610

    # Age effect on CL (power-law about the rounded reference age).
    e_age_cl <- -0.154         ; label("Age power exponent on CL/F (unitless); reference age 35 y")                              # Galluppi 2021 Table 1: Age CL = -0.154

    # Categorical covariate effects on log-CL. The DIS_HEALTHY coefficient
    # is negated relative to the paper because DIS_HEALTHY = 1 - paper's PATIENT.
    e_dis_healthy_cl <- -log(0.809) ; label("Effect of DIS_HEALTHY = 1 on log-CL/F (unitless); = -log(Patient CL)")              # Galluppi 2021 Table 1: Patient CL = 0.809 (paper PATIENT flag = 1 - DIS_HEALTHY)
    e_race_asian_cl  <-  log(0.987) ; label("Effect of RACE_ASIAN = 1 on log-CL/F (unitless); = log(Asian CL)")                  # Galluppi 2021 Table 1: Asian CL = 0.987
    e_sexf_cl        <-  log(0.938) ; label("Effect of SEXF = 1 on log-CL/F (unitless); = log(Female CL)")                       # Galluppi 2021 Table 1: Female CL = 0.938

    # ----------------------------------------------------------------------
    # Inter-individual variability -- Galluppi 2021 Table 1. Full BLOCK(4)
    # covariance matrix on log-scale for CL/F, Vc/F, ka, Vp/F. Order
    # matches the paper's Table 1 rows. Six off-diagonal covariances all
    # estimated. Variances imply %CV via sqrt(exp(omega^2) - 1) * 100:
    # CL 40.4%, Vc 13.7%, ka 75.4%, Vp 40.9%.
    # ----------------------------------------------------------------------
    etalcl + etalvc + etalka + etalvp ~
      c(0.151,
        0.0379,   0.0187,
        -0.0646,  -0.0146,  0.450,
        -0.00568, 0.0255,   0.0573,  0.155)                                                                                     # Galluppi 2021 Table 1: BLOCK(4) IIV variance-covariance (order CL, Vc, ka, Vp)

    # ----------------------------------------------------------------------
    # Residual unexplained variability -- Galluppi 2021 Table 1.
    # Proportional-only error model selected per AIC/BIC (Galluppi 2021
    # Results 'Model development results'). sigma^2 = 0.104 gives
    # propSd = sqrt(0.104) = 0.322 (approximately 32% CV).
    # ----------------------------------------------------------------------
    propSd <- sqrt(0.104)      ; label("Proportional residual error (fraction)")                                                 # Galluppi 2021 Table 1: Residual (proportional) sigma^2 = 0.104
  })

  model({
    # Individual PK parameters. Log-additive categorical covariate effects
    # on CL, power-law allometric weight scaling on CL/F, Q/F, Vc/F, Vp/F,
    # and power-law age scaling on CL only. No covariate effects on ka.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl +
              e_dis_healthy_cl * DIS_HEALTHY +
              e_race_asian_cl  * RACE_ASIAN +
              e_sexf_cl        * SEXF) *
          (WT  / 70) ^ e_wt_cl_q *
          (AGE / 35) ^ e_age_cl
    vc <- exp(lvc + etalvc) * (WT / 70) ^ e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70) ^ e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70) ^ e_wt_vc_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL = 1000 ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
