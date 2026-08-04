Yuen_1995_ganciclovir <- function() {
  description <- paste(
    "Two-compartment population PK model for IV ganciclovir in patients with CMV",
    "retinitis, asymptomatic CMV urine shedding, or a solid-organ transplant with",
    "renal dysfunction (Yuen 1995), with an additive-intercept clearance whose",
    "renal-elimination term scales with body weight and creatinine clearance and is",
    "reduced in transplant recipients and in CMV-retinitis patients.",
    "Parameters transcribed from the Yang 2023 ganciclovir / valganciclovir",
    "population-PK model repository (Table 3), not from the primary publication;",
    "re-verify against Yuen 1995 when the primary is obtained.",
    sep = " "
  )
  reference <- paste(
    "Yuen GJ, Drusano GL, Fletcher C, Capparelli E, Connor JD, Lalezari JP, Drew L,",
    "Follansbee S, Busch D, Jacobson M, et al. Population differences in",
    "ganciclovir clearance as determined by nonlinear mixed-effects modelling.",
    "Antimicrob Agents Chemother. 1995;39(10):2350-2352. doi:10.1128/AAC.39.10.2350.",
    "Parameters transcribed from Yang W, Mak W, Gwee A, Gu M, Wu Y, Shi Y, He Q,",
    "Xiang X, Han B, Zhu X. Establishment and Evaluation of a Parametric Population",
    "Pharmacokinetic Model Repository for Ganciclovir and Valganciclovir.",
    "Pharmaceutics. 2023;15(7):1801. doi:10.3390/pharmaceutics15071801 (Table 3).",
    sep = " "
  )
  vignette <- "Yang_2023_ganciclovir_model_repository"
  units    <- list(time = "hr", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Weight enters un-normalized in three places. On clearance it multiplies",
        "the renal-elimination term (0.168 * BW * CLcr/100), so weight and",
        "creatinine clearance act jointly rather than as separate effects. On both",
        "volumes it is a simple per-kilogram coefficient (Vc = 0.381 * BW,",
        "Vp = 0.511 * BW), i.e. 0.381 and 0.511 L/kg. Weight was retained on CL,",
        "Vc and Vp; the Q effect is recorded as not reported (Yang 2023 Table 4).",
        "Yang 2023 Table 2 records the cohort weights as not reported (NR).",
        sep = " "
      ),
      source_name        = "BW"
    ),
    CRCL = list(
      description        = paste(
        "Creatinine clearance, raw (NOT BSA-normalized). Yang 2023 Table 3",
        "footnote defines CLcr as 'creatinine clearance (mL/min)'; the review does",
        "not state which estimating equation Yuen 1995 used, so the assay form is",
        "to be confirmed against the primary publication.",
        sep = " "
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Raw mL/min, not BSA-normalized -- distinct from the eGFR / CrCL / CKD-EPI /",
        "ASCC variants used by the other models in this repository, which are all",
        "in mL/min/1.73 m^2. Precedent for storing a raw-mL/min renal covariate",
        "under the canonical CRCL column: Delattre_2010_amikacin.R,",
        "Georges_2009_ceftazidime.R, Chen_2023_nemonoxacin.R. Enters the clearance",
        "equation divided by 100 and multiplied by body weight",
        "(0.168 * BW * CLcr/100), so the slope 0.168 has units of L/h per kg per",
        "(100 mL/min). CLcr was retained on CL (Yang 2023 Table 4).",
        sep = " "
      ),
      source_name        = "CLcr"
    ),
    TX_ANY = list(
      description        = "Any-solid-organ-transplant indicator (1 = solid-organ graft recipient, 0 = non-transplant patient)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-transplant patient)",
      notes              = paste(
        "Yang 2023 Table 3 footnote: 'T: T = 0 for non-transplant patients and",
        "0.76 for transplant patients'. The clearance equation applies the factor",
        "as (1 - T), so a transplant recipient's renal-elimination term is reduced",
        "by 76% relative to a non-transplant patient with the same weight and",
        "creatinine clearance. Encoded here as (1 - e_tx_any_cl * TX_ANY) with",
        "e_tx_any_cl = 0.76, so TX_ANY = 0 recovers the reference value.",
        "This is a transplant-vs-non-transplant contrast and therefore NOT",
        "expressible with the organ-specific TX_LIVER / TX_HEART / TX_LUNG",
        "canonicals, which partition within a transplant cohort. Yuen 1995's",
        "cohort was 31 CMV-retinitis patients, 17 CMV-urine-shedding patients, and",
        "5 solid-organ-transplant recipients with renal dysfunction, so only 5 of",
        "53 subjects carry TX_ANY = 1. Transplant status was tested as a",
        "demographic covariate and retained on CL (Yang 2023 Table 4). Founding",
        "example for the TX_ANY canonical.",
        sep = " "
      ),
      source_name        = "T"
    ),
    DIS_CMV_RETINITIS = list(
      description        = "CMV-retinitis indicator (1 = CMV retinitis, 0 = CMV-positive without retinitis)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CMV-positive without retinitis; in Yuen 1995, asymptomatic CMV urine shedding)",
      notes              = paste(
        "Yang 2023 Table 3 footnote: 'CMV: CMV = 0 for CMV-shedding patients and",
        "0.41 for patients with CMV retinitis'. The clearance equation applies the",
        "factor as (1 - CMV), so a retinitis patient's renal-elimination term is",
        "reduced by 41% relative to an asymptomatic urine-shedding patient with the",
        "same weight and creatinine clearance. Encoded here as",
        "(1 - e_dis_cmv_retinitis_cl * DIS_CMV_RETINITIS) with",
        "e_dis_cmv_retinitis_cl = 0.41, so DIS_CMV_RETINITIS = 0 recovers the",
        "reference value. This is a WITHIN-CMV-positive disease-presentation",
        "contrast and therefore orthogonal to DIS_CMV, which encodes",
        "transplant-recipient-with-CMV versus a non-CMV reference: every subject in",
        "the retinitis and shedding groups here would carry DIS_CMV = 1. CMV",
        "presentation was tested as a laboratory covariate and retained on CL",
        "(Yang 2023 Table 4). Founding example for the DIS_CMV_RETINITIS canonical.",
        sep = " "
      ),
      source_name        = "CMV"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 53L,
    n_studies      = 1L,
    n_observations = 558L,
    age_median     = "Not reported in Yang 2023 Table 2.",
    weight_median  = "Not reported in Yang 2023 Table 2.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported.",
    disease_state  = paste(
      "31 patients with CMV retinitis, 17 patients shedding CMV in urine, and 5",
      "solid-organ-transplant recipients with renal dysfunction.",
      sep = " "
    ),
    dose_range     = "IV ganciclovir 1.2-5.0 mg/kg as a 1 h constant-rate infusion.",
    regions        = "United States (prospective).",
    bioassay       = "HPLC, LLOQ 0.25 ug/mL.",
    notes          = paste(
      "Demographics and dosing from Yang 2023 Table 2 (53 subjects, 558",
      "observations; sex, age and weight recorded as not reported). Sampling was",
      "a mixture of sparse (0.5 and 3 h post dose) and intensive (0, 1, 2-3, 5-7",
      "and 10-12 h post dose) schedules. Estimated in NONMEM; the estimation",
      "algorithm is not reported. Covariates tested: weight, transplant status",
      "(yes/no), creatinine clearance, and CMV shedding versus CMV retinitis, with",
      "forward inclusion at p < 0.005; all four were retained on CL, weight on Vc",
      "and Vp (Yang 2023 Table 4). This is the oldest model in the repository. The",
      "primary study reported no model-evaluation results (published before 2000),",
      "but Yang 2023 Section 3.2.1 states that after careful inspection this",
      "model's performance was comparable to the others in the repository. The",
      "primary study's stated purpose was to evaluate the effect of HIV on",
      "ganciclovir clearance.",
      sep = " "
    )
  )

  ini({
    # Structural PK -- Yang 2023 Table 3, Yuen et al. (1995) row. Clearance is
    # additive-linear with a non-zero intercept, so the exponentiated intercept
    # below is the covariate-free component and e_wt_crcl_cl is an additive slope
    # (not a power exponent); the Zhou_1996_ganciclovir.R /
    # Georges_2009_ceftazidime.R pattern. Both volumes are per-kilogram
    # coefficients, so the reference subject for lvc / lvp is WT = 1 kg. Clearances
    # in L/h, volumes in L. IV-only model, so no absorption parameters.
    lcl <- log(0.382); label("Covariate-free intercept of clearance (L/h)")              # Yang 2023 Table 3 (Yuen 1995): CL = 0.382 + 0.168 * BW * CLcr/100 * (1-T) * (1-CMV)
    lvc <- log(0.381); label("Central volume coefficient at WT = 1 kg (Vc, L/kg)")       # Yang 2023 Table 3 (Yuen 1995): Vc = 0.381 * BW
    lvp <- log(0.511); label("Peripheral volume coefficient at WT = 1 kg (Vp, L/kg)")    # Yang 2023 Table 3 (Yuen 1995): Vp = 0.511 * BW
    lq  <- log(13.4) ; label("Inter-compartmental clearance (Q, L/h; not weight-scaled)") # Yang 2023 Table 3 (Yuen 1995): Q = 13.4

    # Covariate effects on clearance. The renal-elimination slope acts on the
    # PRODUCT of body weight and creatinine clearance, so weight and CLcr are not
    # separable effects here. The volume coefficients above already carry their
    # weight dependence, so they need no separate exponent parameter.
    e_wt_crcl_cl <- 0.168; label("Additive slope of WT x CLcr/100 on CL (L/h per kg per 100 mL/min)")  # Yang 2023 Table 3 (Yuen 1995): + 0.168 * BW * CLcr/100

    # Both group effects are FRACTIONAL REDUCTIONS applied as (1 - factor), per
    # the Yang 2023 Table 3 footnotes: T = 0 for non-transplant patients and 0.76
    # for transplant patients; CMV = 0 for CMV-shedding patients and 0.41 for
    # patients with CMV retinitis. The reference subject is therefore a
    # non-transplant patient with asymptomatic CMV urine shedding.
    e_tx_any_cl            <- 0.76; label("Fractional reduction of the renal CL term in transplant recipients (unitless; non-transplant reference)")     # Yang 2023 Table 3 footnote (Yuen 1995): T = 0.76 for transplant patients, applied as (1-T)
    e_dis_cmv_retinitis_cl <- 0.41; label("Fractional reduction of the renal CL term in CMV-retinitis patients (unitless; CMV-shedding reference)")      # Yang 2023 Table 3 footnote (Yuen 1995): CMV = 0.41 for CMV retinitis, applied as (1-CMV)

    # Between-subject variability. Yang 2023 Methods: %CV = sqrt(omega^2) * 100%,
    # so variance = (BSV% / 100)^2. Q and Vp carry no BSV in the source table.
    etalcl ~ 0.225625  # Yang 2023 Table 3 (Yuen 1995): BSV CL = 47.5% -> 0.475^2
    etalvc ~ 0.075625  # Yang 2023 Table 3 (Yuen 1995): BSV Vc = 27.5% -> 0.275^2

    # Residual unexplained variability: proportional.
    propSd <- 0.361; label("Proportional residual error (fraction)")  # Yang 2023 Table 3 (Yuen 1995): 36.1% proportional error
  })

  model({
    # Additive-intercept clearance: the 0.382 L/h intercept is covariate-free and
    # the renal-elimination term is scaled down multiplicatively in transplant
    # recipients and in CMV-retinitis patients. Exponential between-subject
    # variability is applied to the assembled typical value (the
    # Zhou_1996_ganciclovir.R / Georges_2009_ceftazidime.R pattern for
    # additive-intercept clearance).
    cl <- (
      exp(lcl) +
        e_wt_crcl_cl * WT * CRCL / 100 *
        (1 - e_tx_any_cl * TX_ANY) *
        (1 - e_dis_cmv_retinitis_cl * DIS_CMV_RETINITIS)
    ) * exp(etalcl)
    vc <- exp(lvc) * WT * exp(etalvc)
    vp <- exp(lvp) * WT
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> concentration in mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
