Zhou_2026_rituximab <- function() {
  description <- "Two-compartment population PK model for intravenous rituximab with a power effect of the serum albumin-globulin ratio on clearance, in treatment-naive Chinese adults with diffuse large B-cell lymphoma receiving R-CHOP induction therapy"
  reference <- paste(
    "Zhou WY, Ling J, Guan MM, Qu Y, Dan YH, Chen R (2026).",
    "Population Pharmacokinetics of Rituximab in Treatment-Naive Chinese Patients",
    "With Diffuse Large B-cell Lymphoma During Induction Therapy: Clinical Implications.",
    "Clin Transl Sci 19(8):e70683. doi:10.1111/cts.70683.",
    sep = " "
  )
  vignette <- "Zhou_2026_rituximab"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    ALB_GLOB_RATIO = list(
      description = "Serum albumin-to-globulin ratio (A/G ratio)",
      units = "(unitless ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Ratio of serum albumin to serum globulin, both in g/L. Reference",
        "(normalising) value 1.55, the modelling-cohort mean (Table 1: 1.55 +/-",
        "0.35). CL scales as (A/G / 1.55)^-1.49 (Eq 5), so a lower ratio gives a",
        "higher clearance. The paper does not state whether the ratio was carried",
        "per subject (baseline) or per cycle; the inclusion criteria required the",
        "biochemistry to be measured within 1 week before rituximab treatment.",
        sep = " "
      ),
      source_name = "A/G ratio"
    )
  )

  # Covariates screened in the stepwise search (forward inclusion p < 0.05,
  # backward elimination p < 0.001) but not retained. BSA passed forward
  # inclusion and failed backward elimination; the rest failed forward inclusion.
  # Terms with no canonical column (globulin alone) are listed in
  # population$notes instead.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained."
    ),
    AGE = list(description = "Age", units = "years", type = "continuous", notes = "Screened; not retained."),
    WT = list(description = "Body weight", units = "kg", type = "continuous", notes = "Screened; not retained."),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Significant in forward inclusion; removed in backward elimination (dOFV < 10.83)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened; not retained."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained."
    ),
    LDH = list(
      description = "Lactate dehydrogenase",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained."
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/L",
      type = "continuous",
      notes = "Screened; not retained."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened; not retained (the A/G ratio was retained instead)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened; not retained."
    ),
    HGB = list(description = "Hemoglobin", units = "g/L", type = "continuous", notes = "Screened; not retained."),
    PLT = list(
      description = "Platelet count",
      units = "10^9/L",
      type = "continuous",
      notes = "Screened; not retained."
    ),
    CYCLE = list(
      description = "Chemotherapy cycle number",
      units = "(count)",
      type = "count",
      notes = "Screened; not retained."
    )
  )

  compartmentData <- list(
    central = list(analyte = "rituximab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "rituximab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 45,
    n_studies = 1,
    n_observations = "115 serum rituximab concentrations (index dataset); external validation 12 patients / 28 concentrations",
    age_mean_sd = "66.44 +/- 10.85 years",
    weight_mean_sd = "61.34 +/- 9.02 kg",
    bsa_mean_sd = "1.76 +/- 0.14 m^2",
    alb_glob_ratio_mean_sd = "1.55 +/- 0.35 (interquartile range 1.33-1.66)",
    sex_female_pct = 64.4,
    race_ethnicity = c(Asian = 100),
    disease_state = "newly diagnosed, treatment-naive CD20-positive diffuse large B-cell lymphoma on R-CHOP induction therapy",
    dose_range = "375 mg/m^2 IV (500-750 mg per infusion) every 3 weeks (occasionally every 4 weeks) for 6-8 cycles; rate-escalated infusion (first infusion 50 mg/h rising by 50 mg/h every 30 min to 400 mg/h; later infusions 100 mg/h rising by 100 mg/h every 30 min to 400 mg/h)",
    regions = "China (single centre, Third Affiliated Hospital of Soochow University, Changzhou)",
    notes = paste(
      "Retrospective sparse sampling from clinical residual serum (pre-dose",
      "troughs and post-infusion samples during hospitalisation across",
      "induction cycles); assay by magnetic-particle chemiluminescence",
      "immunoassay. Baseline characteristics in Table 1 (16 male / 29 female;",
      "mean chemotherapy cycle 3.51). Covariates screened but not retained",
      "without a canonical column: serum globulin (24.15 +/- 6.03 g/L). Race is",
      "not stated; recorded as Asian for the single-centre Chinese cohort.",
      sep = " "
    )
  )

  ini({
    lcl <- log(0.0181); label("Clearance at A/G ratio 1.55 (L/h)")         # Table 2 CL 0.0181 L/h (RSE 19%); Eq 5
    lvc <- log(8.12); label("Central volume of distribution V1 (L)")       # Table 2 V1 8.12 L (RSE 8%); Eq 6
    lq <- log(0.0213); label("Intercompartmental clearance Q (L/h)")       # Table 2 Q 0.0213 L/h (RSE 29%); Eq 7
    lvp <- log(29.2); label("Peripheral volume of distribution V2 (L)")    # Table 2 V2 29.2 L (RSE 47%); Eq 8

    e_alb_glob_ratio_cl <- -1.49; label("Power exponent of the albumin-globulin ratio on CL (unitless)") # Table 2 'A/G Ratio' -1.49 (RSE 18%); Eq 5

    # IIV estimated on CL only (Results); Table 2 'IIV of CL (%)' 36.2, a
    # log-normal CV (Methods: lognormal IIV), converted with log(1 + CV^2).
    etalcl ~ log(1 + 0.362^2)

    propSd <- 0.316; label("Proportional residual error (fraction)")       # Table 2 'RV (%)' 31.6 (RSE 24%); Methods: proportional error model
  })

  model({
    # Reference A/G ratio: modelling-cohort mean, the normaliser in Eq 5.
    agRef <- 1.55

    cl <- exp(lcl + etalcl) * (ALB_GLOB_RATIO / agRef)^e_alb_glob_ratio_cl
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L == ug/mL
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
