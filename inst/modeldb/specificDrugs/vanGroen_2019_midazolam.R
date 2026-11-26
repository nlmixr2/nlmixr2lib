vanGroen_2019_midazolam <- function() {
  description <- paste(
    "Two-compartment population PK model for an intravenous [14C]midazolam",
    "microdose (111 Bq/kg; 37.6 ng/kg, infused over 30 min) in 15 infants",
    "aged 0.6-49.1 weeks postnatal (2.6-8.9 kg) admitted to paediatric or",
    "neonatal intensive care, given either as an isolated microdose or as a",
    "microtracer mixed into a therapeutic midazolam loading dose (van Groen",
    "2019). Body weight scales the central volume by a power function",
    "referenced to 4 kg; clearance, intercompartmental clearance and",
    "peripheral volume are weight-independent. IIV on clearance only;",
    "combined additive and proportional residual error. Fitted on the",
    "radioactivity scale: doses in Bq and plasma [14C]midazolam in Bq/L."
  )
  reference <- paste(
    "van Groen BD, Vaes WH, Park BK, Krekels EHJ, van Duijn E, Korgvee LT,",
    "Maruszak W, Grynkiewicz G, Garner RC, Knibbe CAJ, Tibboel D,",
    "de Wildt SN, Turner MA (2019).",
    "Dose-linearity of the pharmacokinetics of an intravenous [14C]midazolam",
    "microdose in children.",
    "Br J Clin Pharmacol 85(10):2332-2340.",
    "doi:10.1111/bcp.14047."
  )
  vignette <- "vanGroen_2019_midazolam"
  units <- list(time = "h", dosing = "Bq", concentration = "Bq/L")

  compartmentData <- list(
    central = list(
      analyte = "[14C]midazolam",
      units = "Bq",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "[14C]midazolam",
      units = "Bq",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on the central volume, V1 = V1_4kg * (WT/4)^k1",
        "(Table 3). The 4 kg reference is close to the cohort median of",
        "3.6 kg (Table 1, range 2.6-8.9 kg). Time-fixed per subject."
      ),
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    PNA = list(
      description = "Postnatal age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on all parameters (Methods, Nonlinear mixed effects",
        "modelling) and not significant once body weight was on V1",
        "(Results). Cohort median 11.4 (0.6-49.1) weeks."
      )
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on all parameters; not retained (Results)."
    ),
    CREAT = list(
      description = "Plasma creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on all parameters; not retained. Cohort median 35 (20-51) umol/L."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on all parameters; not retained. Cohort median 17 (7-68) U/L."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on all parameters; not retained. Cohort median 42 (12-93) U/L."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on all parameters; not retained. Cohort median 9 (2-274) umol/L."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 15L,
    n_studies = 1L,
    age_range = "Postnatal age 0.6-49.1 weeks (median 11.4); gestational age at birth 23.9-41.4 weeks (median 39.4).",
    weight_range = "2.6-8.9 kg (median 3.6 kg).",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported.",
    disease_state = paste(
      "Critically ill preterm-to-2-year-old infants on the paediatric and",
      "neonatal intensive care wards, with intravenous access. Serious",
      "hepatic impairment (ASAT/ALAT > 200 U/L), renal impairment (plasma",
      "creatinine > 150 umol/L), haemofiltration, dialysis and ECMO were",
      "exclusions."
    ),
    dose_range = paste(
      "Single IV [14C]midazolam 111 Bq/kg (37.6 ng/kg) over 30 min, either",
      "alone (microdose, n = 9) or mixed into the first therapeutic",
      "midazolam loading dose (microtracer, n = 6)."
    ),
    regions = "Alder Hey Children's NHS Foundation Trust and Liverpool Women's NHS Foundation Trust, Liverpool, UK.",
    notes = paste(
      "67 plasma samples (up to 5 per infant, 0.17-36 h post-dose), of",
      "which 8 below the AMS detection limit were excluded. [14C]midazolam",
      "was quantified by UPLC-AMS (LLOQ 0.31 mBq/mL = 0.31 Bq/L). Sex was",
      "not reported. Treatment arm (microdose vs microtracer) was tested on",
      "every parameter and was not significant. NONMEM 7.4; 50-replicate",
      "bootstrap."
    )
  )

  ini({
    # Table 3 'Parameter estimates of the pharmacokinetic model for IV
    # [14C]midazolam', Estimate (RSE%) column.
    lcl <- log(2.06)
    label("Clearance (L/h)") # Table 3 CL = 2.06 L/h (RSE 24%)
    lq <- log(0.79)
    label("Intercompartmental clearance (L/h)") # Table 3 Q = 0.79 L/h (RSE 44%)
    lvc <- log(3.81)
    label("Central volume of distribution at 4 kg body weight (L)") # Table 3 V1_4kg = 3.81 L (RSE 8%)
    e_wt_vc <- 1.36
    label("Power exponent of (WT/4) on the central volume (unitless)") # Table 3 k1 = 1.36 (RSE 10%)
    lvp <- log(3.19)
    label("Peripheral volume of distribution (L)") # Table 3 V2 = 3.19 L (RSE 18%)

    etalcl ~ 0.73 # Table 3 'omega2 CL' = 0.73 (RSE 42%), a variance on the log scale

    # Table 3 residual rows read as NONMEM $SIGMA variances (see the vignette
    # for the evidence): SD = sqrt(variance).
    propSd <- sqrt(0.09)
    label("Proportional residual error (fraction)") # Table 3 'Proportional error' = 0.09 (RSE 24%), variance -> SD 0.30
    addSd <- sqrt(0.08)
    label("Additive residual error (Bq/L)") # Table 3 'Additional error' = 0.08 (RSE 50%), variance in (Bq/L)^2 -> SD 0.283 Bq/L
  })
  model({
    cl <- exp(lcl + etalcl)
    q <- exp(lq)
    vc <- exp(lvc) * (WT / 4)^e_wt_vc
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
