Harwood_2025_tobramycin <- function() {
  description <- "One-compartment population PK model for intravenous tobramycin in children and adolescents with cystic fibrosis (Harwood 2025), with allometric weight scaling of CL and V, a Rhodin 2009 GFR maturation function of postmenstrual age on CL, a modified-Schwartz eGFR power effect on CL, and an estimated zero-order infusion duration."
  reference <- "Harwood KH, Duffull S, Shanthikumar S, Lei A, Ranganathan S, Robinson P, Sandaradura I, Lai T, Gwee A. Personalized tobramycin dosing in children with cystic fibrosis: an AUC24-guided approach. Antimicrob Agents Chemother. 2025;69(9):e00278-25. doi:10.1128/aac.00278-25"
  vignette <- "Harwood_2025_tobramycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric power scaling on CL (estimated exponent 0.663) and V (estimated exponent 0.744), centred on 70 kg (Supplementary Equations 1-2; supplement covariate-handling item 1 states weight, unlike the other continuous covariates, was centred on 70 kg). Recorded per admission in the source data.",
      source_name = "WT"
    ),
    AGE = list(
      description = "Postnatal age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Combined with GA inside model() to form postmenstrual age in weeks, PMA = AGE * 52 + GA, exactly as written in Harwood 2025 Equation 1 ((AGE * 52) + 40, with the term-birth GA of 40 weeks substituted).",
      source_name = "AGE"
    ),
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters only through PMA = AGE * 52 + GA in the GFR maturation function. Harwood 2025 Discussion: gestational age was not routinely recorded once postnatal age exceeded 1-2 years and 40 weeks was assumed if not otherwise stated; the published dosing Equation 1 hard-codes 40. Set GA = 40 for term-born or unknown-gestation subjects.",
      source_name = "GA"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only to compute the modified (bedside) Schwartz eGFR = 41.3 * HT[m] / CREAT[mg/dL] (Supplementary Equation 5). The source equation takes height in metres; model() divides the canonical cm value by 100.",
      source_name = "HT"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only to compute the modified (bedside) Schwartz eGFR = 41.3 * HT[m] / CREAT[mg/dL] (Supplementary Equation 5; the supplement abbreviation list gives CR in mg/dL). Harwood 2025 Table 1 reports creatinine in umol/L (median 43, range 14-107); divide umol/L by 88.4 to obtain mg/dL.",
      source_name = "CR"
    )
  )

  compartmentData <- list(
    central = list(analyte = "tobramycin", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 63L,
    n_studies = 1L,
    n_observations = 450L,
    n_admissions = 94L,
    age_range = "0.2-19.3 years",
    age_median = "13.3 years",
    weight_range = "5.0-81.7 kg",
    weight_median = "43.8 kg",
    height_range = "0.56-1.85 m",
    height_median = "1.53 m",
    creatinine_range = "14-107 umol/L (median 43)",
    egfr_range = "58-265 mL/min (median 127; modified Schwartz)",
    sex_female_pct = 100 * (63 - 34) / 63,
    disease_state = "Cystic fibrosis with a respiratory (pulmonary) exacerbation treated with IV tobramycin; 86% had CF-related comorbidities (pancreatic insufficiency 71%, bronchiectasis 57%, diabetes 24%, liver disease 24%).",
    dose_range = "Initial dose median 10 mg/kg (range 0.2-14) once daily as a 30-min IV infusion; median course 12 days (range 1-26).",
    regions = "Australia (Royal Children's Hospital, Melbourne; The Children's Hospital at Westmead, Sydney)",
    co_medication = "Concomitant beta-lactam in 98% of admissions; CFTR modulators in 36%, NSAIDs 7%, diuretics 5%. None retained as covariates.",
    notes = "Multicenter prospective observational TDM study, August 2020 - August 2022 (Harwood 2025 Methods and Tables 1-2). 34/63 male. Age groups by admission: <2 years 10 (10.6%), 2-12 years 28 (29.8%), >12 years 56 (59.6%). Serum concentrations by immunoassay (LLOQ 0.6 mg/L); 29/450 (6%) BLQ handled by the M6 method."
  )

  ini({
    # Structural parameters, typical values for a 70 kg, fully renally
    # mature subject with eGFR = 127.3 mL/min (Supplementary Eqs. 1-5).
    lcl <- log(6.17); label("Clearance for a 70 kg fully mature subject at eGFR 127.3 mL/min (L/h)") # Table 3: CL 6.17 L/h (RSE 5.35%)
    lvc <- log(15.5); label("Volume of distribution for a 70 kg subject (L)") # Table 3: V 15.5 L (RSE 6.97%)
    ld1 <- log(0.801); label("Zero-order infusion duration D1 (h)") # Table 3: D1 0.801 h (RSE 10.4%)

    # Covariate exponents
    e_wt_cl <- 0.663; label("Allometric exponent on (WT/70) for CL (unitless)") # Table 3: WTCL 0.663 (RSE 7.33%)
    e_wt_vc <- 0.744; label("Allometric exponent on (WT/70) for V (unitless)") # Table 3: WTV 0.744 (RSE 9.13%)
    e_egfr_cl <- fixed(0.4); label("Power exponent on (eGFR/127.3) for CL (unitless)") # Table 3: eGFRCL 0.4 (fixed)

    # Renal maturation (Rhodin 2009 GFR parameters, fixed from the literature)
    tm50_cl <- fixed(47.7); label("Postmenstrual age at 50% GFR maturation (weeks)") # Supplementary Eq. 4 and Eq. 1: 47.7, from Rhodin 2009 Table 3
    hill <- fixed(3.4); label("Hill coefficient of the GFR maturation function (unitless)") # Supplementary Eq. 4 and Eq. 1: 3.4, from Rhodin 2009 Table 3

    # IIV: Table 3 reports the SD of each eta; variance = SD^2
    etalcl ~ 0.170^2 # Table 3: omega CL 0.170 (SD; RSE 17.2%)
    etalvc ~ 0.173^2 # Table 3: omega V 0.173 (SD; RSE 14.4%)
    etald1 ~ 0.404^2 # Table 3: omega D1 0.404 (SD; RSE 35.0%)

    # Residual error, combined additive + proportional
    addSd <- 0.173; label("Additive residual error (mg/L)") # Table 3: Additive error 0.173 (RSE 23.7%)
    propSd <- 0.103; label("Proportional residual error (fraction)") # Table 3: Proportional error 0.103 (RSE 9.09%)
  })

  model({
    # Postmenstrual age in weeks (Eq. 1: (AGE * 52) + 40, with GA = 40)
    pma <- AGE * 52 + GA

    # GFR maturation fraction (Supplementary Eq. 4)
    mf <- pma^hill / (tm50_cl^hill + pma^hill)

    # Modified (bedside) Schwartz eGFR, HT in m and CREAT in mg/dL
    # (Supplementary Eq. 5)
    egfr <- 41.3 * (HT / 100) / CREAT

    # Individual parameters (Supplementary Eqs. 1-3)
    cl <- exp(lcl + etalcl) *
      mf *
      (WT / 70)^e_wt_cl *
      (egfr / 127.3)^e_egfr_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    d1 <- exp(ld1 + etald1)

    kel <- cl / vc

    d / dt(central) <- -kel * central

    # Estimated infusion duration; dose records need rate = -2
    dur(central) <- d1

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
