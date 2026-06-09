Wade_2008_fluconazole <- function() {
  description <- "One-compartment intravenous population PK model for fluconazole in preterm and term infants (gestational age 23-40 weeks, postnatal age <120 days) with allometric body weight on CL and V (fixed exponents 0.75 and 1.0, reference 1 kg), power effects of gestational age at birth (reference 26 weeks) and postnatal age (reference 2 weeks) on CL, and an on/off power effect of serum creatinine on CL gated when SCR > 1 mg/dL (Wade 2008)."
  reference <- "Wade KC, Wu D, Kaufman DA, Ward RM, Benjamin DK Jr, Sullivan JE, Ramey N, Jayaraman B, Hoppu K, Adamson PC, Gastonguay MR, Barrett JS. Population pharmacokinetics of fluconazole in young infants. Antimicrob Agents Chemother. 2008. doi:10.1128/AAC.00569-08"
  vignette <- "Wade_2008_fluconazole"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; missing weights carried forward up to 7 days per Wade 2008 Methods. Allometric scaling on CL (fixed exponent 0.75) and V (fixed exponent 1.0) with reference 1 kg (Wade 2008 Table 2 Base model row and abstract).",
      source_name        = "WT"
    ),
    GA = list(
      description        = "Gestational age at birth",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject; Wade 2008 uses notation BGA. Power covariate on CL with reference 26 weeks (the cohort median).",
      source_name        = "BGA"
    ),
    PNA = list(
      description        = "Postnatal age (time since birth)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Wade 2008 expresses PNA as weeks of life (day-of-life divided by 7) with reference 2 weeks; the canonical PNA carries months, so the reference is reparameterised inside model() as 2 * 7 / 30.4375 = 0.460 months (the conversion uses 30.4375 days/month per the Zhao 2018 PNA precedent). Cohort PNA at enrolment ranged from 0.14 to 12.6 weeks (Wade 2008 abstract).",
      source_name        = "PNA"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; Wade 2008 notation SCRT. Per Wade 2008 Methods, infants without a measured value were assumed to have SCR <= 1.0 mg/dL, and measured SCR > 1.0 mg/dL was carried forward for up to 7 days. Inside model() the dichotomous renal-impairment flag CR is derived from CREAT > 1, so values <= 1 (including imputed 1.0 for missing) make the SCR term collapse to 1.",
      source_name        = "SCRT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 55L,
    n_studies      = 2L,
    age_range      = "PNA at enrolment 1-88 days (cohort median 16 days); BGA 23-40 weeks",
    age_median     = "PNA 16 days; BGA 26 weeks",
    weight_range   = "0.451-7.125 kg (Wade 2008 Table 1)",
    weight_median  = "1.020 kg (Wade 2008 Table 1)",
    sex_female_pct = 44,
    race_ethnicity = c(Caucasian = 50, Black = 40, Other = 10, Hispanic = 9),
    disease_state  = "Preterm and term infants <120 days of age receiving intravenous fluconazole for prevention or treatment of invasive candidiasis (NICHD Pediatric Pharmacology Research Unit Network).",
    dose_range     = "3-12 mg/kg/dose intravenous fluconazole (clinical-care dosing, not protocolised).",
    regions        = "United States plus Helsinki, Finland (multicentre PPRU network).",
    samples        = "357 plasma fluconazole observations (217 prospectively timed, 140 scavenged from discarded clinical specimens); median 6.5 samples per infant (range 1-16).",
    notes          = "Demographics from Wade 2008 Table 1. Hispanic ethnicity overlaps with race categories (9% Hispanic of any race). The sex_female_pct value is computed as 100 - 56 (% male)."
  )

  ini({
    # Structural parameters (Wade 2008 Table 3 final-model column, reference
    # weight 1 kg, reference BGA 26 weeks, reference PNA 2 weeks, CR = 0).
    lcl <- log(0.015); label("Clearance at 1 kg, 26 weeks BGA, 2 weeks PNA, normal renal function (L/h)") # Wade 2008 Table 3: theta_CL = 0.015 L/h
    lvc <- log(1.024); label("Volume of distribution at 1 kg (L)")                                          # Wade 2008 Table 3: theta_V  = 1.024 L

    # Allometric exponents on body weight. Wade 2008 Table 2 'Base model'
    # rows fix the CL exponent at 0.75 and the V exponent at 1; only the
    # population intercepts (theta_CL, theta_V) and the covariate effects
    # below are estimated.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL (unitless)") # Wade 2008 Table 2: CL = theta_CL * (wt)^0.75
    e_wt_vc <- fixed(1.00); label("Allometric exponent on V (unitless)")  # Wade 2008 Table 2: V  = theta_V  * (wt)^1

    # Covariate effects on CL (Wade 2008 Table 3).
    e_ga_cl    <-  1.739; label("Exponent on (GA / 26 weeks) for CL (unitless)")                 # Wade 2008 Table 3: theta_CL-BGA   =  1.739
    e_pna_cl   <-  0.237; label("Exponent on (PNA / 2 weeks) for CL (unitless)")                 # Wade 2008 Table 3: theta_CL-PNA   =  0.237
    e_creat_cl <- -4.896; label("Exponent on serum creatinine for CL when CR = 1 (unitless)")    # Wade 2008 Table 3: theta_CL-SCRT  = -4.896

    # Inter-individual variability. Wade 2008 Table 3 reports the omega
    # block on the log-normal eta scale as variances (Omega 1,1 = var(CL),
    # Omega 2,2 = var(V), Omega 1,2 = cov(CL,V)). Exponential IIV model
    # per Wade 2008 Methods 'Population PK analysis'.
    etalcl + etalvc ~ c(0.11,
                        0.014, 0.057)  # Wade 2008 Table 3: var(CL) = 0.11, cov(CL,V) = 0.014, var(V) = 0.057

    # Residual unexplained variability. Wade 2008 Table 3 reports sigma
    # as variances for two strata (prospective PK samples vs scavenged
    # samples). The library encoding uses only the prospective PK
    # residual-error variances (Sigma 1 and Sigma 2), converted to SDs:
    # propSd = sqrt(0.027), addSd = sqrt(0.04). The scavenged-sample
    # residual error (Sigma 3 / Sigma 4) and the 0.953 multiplicative
    # bias factor (theta_SCAV) are documented as deviations in the
    # vignette - they are properties of the sampling protocol rather
    # than of the underlying fluconazole PK and would not apply to a
    # prospectively-designed validation study.
    propSd <- 0.1643; label("Proportional residual SD for prospective PK samples (fraction)") # Wade 2008 Table 3: Sigma1 = 0.027 (variance) -> SD = sqrt(0.027) = 0.1643
    addSd  <- 0.2;    label("Additive residual SD for prospective PK samples (ug/mL)")        # Wade 2008 Table 3: Sigma2 = 0.04  (variance) -> SD = sqrt(0.04)  = 0.2
  })

  model({
    # Renal-impairment switch CR: 1 when measured SCR exceeds 1 mg/dL,
    # 0 otherwise. Per Wade 2008 Methods, infants with no measured
    # creatinine value are assumed to have SCR <= 1.0 mg/dL, so the user
    # should supply CREAT <= 1 (e.g. 1.0) for missing measurements; CR
    # then collapses the SCR power term to 1.
    CR <- (CREAT > 1) * 1.0

    # Wade 2008 PNA reference is 2 weeks; canonical PNA is in months.
    # 2 weeks = 2 * 7 / 30.4375 = 0.460 months (using 30.4375 days/month
    # per the Zhao 2018 PNA precedent in inst/references/covariate-columns.md).
    f_pna_cl <- (PNA / (2 * 7 / 30.4375))^e_pna_cl

    # Power covariate effects (Wade 2008 Eq. for CL in abstract):
    # CL = theta_CL * (WT/1)^0.75 * (BGA/26)^theta_BGA * (PNA/2)^theta_PNA
    #      * SCRT^(theta_SCRT * CR)
    f_ga_cl   <- (GA / 26)^e_ga_cl
    f_scr_cl  <- CREAT^(e_creat_cl * CR)

    cl <- exp(lcl + etalcl) * (WT / 1)^e_wt_cl * f_ga_cl * f_pna_cl * f_scr_cl
    vc <- exp(lvc + etalvc) * (WT / 1)^e_wt_vc

    kel <- cl / vc

    # IV fluconazole dosed directly into the central compartment.
    d/dt(central) <- -kel * central

    # Plasma fluconazole concentration. Dose in mg, vc in L -> mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
