Matcha_2025_amikacin <- function() {
  description <- "Two-compartment IV population PK model for amikacin in term neonates (Matcha 2025); clearance carries an uncentred exponential body-weight effect (Exp^(WT * 0.308)) and a power renal-function ratio (SCrmean / SCr)^0.397, where SCrmean is the age-typical population serum creatinine predicted from postnatal age by the Wang 2019 paediatric reference equation, so postnatal age enters clearance as renal maturation; central volume, peripheral volume, and intercompartmental clearance carry no covariates. Supports the paper's WT / SCr / PNA dosing nomogram targeting peak 24-35 mg/L and trough 2-5 mg/L."
  reference <- paste(
    "Matcha S, Dilli Batcha JS, Raju AP, Chaudhari BB, Moorkoth S,",
    "Mallayasamy S, Lewis LE. Precision dosing of amikacin in term neonates",
    "using pharmacometric approach. Pediatr Res. 2025;98:936-941.",
    "doi:10.1038/s41390-025-04044-7.",
    "The age-typical serum-creatinine reference equation (Matcha 2025 Eq. 8)",
    "is taken from Wang H, Sherwin C, Gobburu JVS, Ivaturi V. Population",
    "pharmacokinetic modeling of gentamicin in pediatrics.",
    "J Clin Pharmacol. 2019;59:1584-1596.",
    sep = " "
  )
  vignette <- "Matcha_2025_amikacin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central     = list(analyte = "amikacin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "amikacin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Current body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Enters CL through an UNCENTRED EXPONENTIAL form,",
        "Exp^(WT * 0.308) (Matcha 2025 Eq. 4) -- not the more common",
        "allometric power form. Because the form is uncentred, the lcl",
        "population value (0.064 L/h) is the clearance extrapolated to",
        "WT = 0 kg and is NOT the typical clearance of a neonate. At the",
        "cohort median weight of 2.85 kg, typical CL is 0.154 L/h when",
        "measured creatinine equals the age-typical reference, and",
        "0.128 L/h at the full cohort median (SCr 0.55 mg/dL,",
        "PNA 5.02 days, where the renal-function factor is 0.83).",
        "Cohort range 1.74-4.84 kg (Table 1); the dosing nomogram",
        "(Table 3) spans 2.00-4.50 kg.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    CREAT = list(
      description        = "Measured individual serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Paired with the internally-derived age-typical",
        "reference creatinine (see the PNA entry): CL carries the factor",
        "(SCrmean / CREAT)^0.397, so a neonate whose creatinine equals the",
        "age-typical value has factor 1 and a neonate with higher",
        "creatinine has lower clearance. Cohort median 0.55 mg/dL,",
        "range 0.16-1.49 mg/dL (Table 1); the dosing nomogram (Table 3)",
        "spans 0.15-1.50 mg/dL. Must be supplied in mg/dL, because the",
        "Wang 2019 reference equation returns mg/dL and the two values",
        "are divided.",
        sep = " "
      ),
      source_name        = "SCr"
    ),
    PNA = list(
      description        = "Postnatal age (chronological age since birth)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. The canonical PNA column is in MONTHS, but the",
        "published equation takes PNA in DAYS, so model() converts with",
        "PNA_days = PNA * 30.4375 (the same reparameterisation precedent",
        "recorded for Zhao 2018 in covariate-columns.md). PNA is not a",
        "direct covariate on any PK parameter: it drives the age-typical",
        "reference creatinine SCrmean (Matcha 2025 Eq. 8, from Wang 2019),",
        "which is the numerator of the renal-function ratio on CL. This is",
        "how the paper incorporates renal maturation -- SCrmean rises from",
        "about 0.25 mg/dL at 1 day to about 0.91 mg/dL at 28 days, so",
        "typical CL rises with postnatal age. PNA must be strictly",
        "positive; the log(PNA) term is undefined at birth. Cohort PNA at",
        "sampling 0.62-27.32 days, median 5.02 days (Table 1).",
        sep = " "
      ),
      source_name        = "PNA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 78,
    n_studies      = 1,
    n_observations = 100,
    age_range      = "PNA at start of amikacin 0.04-21.29 days; PNA at sample collection 0.62-27.32 days",
    age_median     = "PNA 1.70 days at start of amikacin; 5.02 days at sample collection",
    ga_range       = "Gestational age 37-42 weeks (term only); postmenstrual age at start of amikacin 37.11-42.22 weeks",
    ga_median      = "Gestational age 39 weeks; postmenstrual age 39.20 weeks",
    weight_range   = "Current body weight 1.74-4.84 kg; birth weight 1.75-4.50 kg",
    weight_median  = "Current body weight 2.85 kg; birth weight 2.95 kg",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported; single-centre Indian cohort",
    disease_state  = "Term neonates prescribed amikacin for suspected or confirmed sepsis; clinically unstable neonates were excluded.",
    dose_range     = "25-110 mg IV (median 40 mg), equivalent to 8.23-51.89 mg/kg (median 12.52 mg/kg)",
    regions        = "India (single-centre prospective observational study; Kasturba Medical College, Manipal)",
    renal_function = "Serum creatinine 0.16-1.49 mg/dL (median 0.55); Schwartz-estimated creatinine clearance 12.68-145.0 mL/min (median 40.98), computed with k = 0.45",
    notes          = paste(
      "Demographics from Matcha 2025 Table 1. 101 concentrations were",
      "collected from 80 neonates; after removing one outlying observation",
      "(154 mg/L) and one subject with all concentrations missing,",
      "100 concentrations from 78 subjects were used for model building",
      "(supplementary Table S1). Sampling was opportunistic at",
      "0.5-24 h after the most recent dose. Model fitted in Pumas",
      "(Julia 1.6); evaluated by VPC (n = 500) and bootstrap (n = 1000,",
      "100% success).",
      sep = " "
    ),
    sex_note       = "Sex distribution not reported in Matcha 2025."
  )

  ini({
    # Structural parameters -- Matcha 2025 Table 2 (typical estimates) and
    # Eqs. 4-7. NOTE: the WT effect on CL is uncentred and exponential, so
    # lcl is the CL intercept extrapolated to WT = 0 kg, not a neonatal
    # typical value. See covariateData$WT$notes.
    lcl <- log(0.064); label("Clearance intercept at WT = 0 kg with serum creatinine equal to the age-typical reference (L/h)")  # Matcha 2025 Table 2 TVCL; Eq. 4
    lvc <- log(1.281); label("Volume of the central compartment (L)")                                                            # Matcha 2025 Table 2 TVVC; Eq. 5
    lq  <- log(0.055); label("Intercompartmental clearance (L/h)")                                                               # Matcha 2025 Table 2 TVQ; Eq. 7
    lvp <- log(0.618); label("Volume of the peripheral compartment (L)")                                                         # Matcha 2025 Table 2 TVVP; Eq. 6

    # Covariate effects on CL
    e_wt_cl    <- 0.308; label("Exponential body-weight coefficient on CL, as Exp^(e_wt_cl * WT) (1/kg)")            # Matcha 2025 Table 2 Wt_cl; Eq. 4
    e_creat_cl <- 0.397; label("Exponent on the age-typical-to-measured serum creatinine ratio for CL (unitless)")   # Matcha 2025 Table 2 SCr_cl; Eq. 4

    # Between-subject variability. The paper reports BSV as %CV on
    # log-normally distributed parameters, converted here with
    # omega^2 = log(CV^2 + 1). Covariance between CL and VC was
    # statistically significant but was judged not clinically meaningful
    # and excluded from the final model, so the omega matrix is diagonal.
    etalcl ~ 0.106966  # CV 33.6%   # Matcha 2025 Table 2 BSV_cl
    etalvc ~ 0.156081  # CV 41.1%   # Matcha 2025 Table 2 BSV_vc

    # Residual unexplained variability
    propSd <- 0.309; label("Proportional residual error (fraction)")  # Matcha 2025 Table 2 Proportional residual error 30.9 %CV
  })

  model({
    # 1. Derived covariate terms.
    #    Age-typical population serum creatinine predicted from postnatal
    #    age (Matcha 2025 Eq. 8, adopted from Wang 2019 and also used in
    #    supplementary Table S1 to impute missing creatinine). The published
    #    equation takes PNA in days; the canonical PNA column is in months.
    pnaDays  <- PNA * 30.4375
    creatRef <- -0.02324 - 0.14545 * log(pnaDays) + 0.26964 * pnaDays^0.5

    #    Renal function: age-typical creatinine over measured creatinine.
    scrFactor <- (creatRef / CREAT)^e_creat_cl

    # 2. Individual PK parameters (Matcha 2025 Eqs. 4-7)
    cl <- exp(lcl + etalcl) * exp(e_wt_cl * WT) * scrFactor
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. Two-compartment IV disposition; amikacin is given intravenously
    #    into the central compartment.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Observation and error
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
