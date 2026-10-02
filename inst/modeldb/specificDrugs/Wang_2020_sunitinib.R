Wang_2020_sunitinib <- function() {
  description <- "Two-compartment population PK model with first-order absorption and lag time for oral sunitinib (parent drug) in children and young adults (2-21 years) with refractory solid tumours, with body surface area (BSA) as a linear covariate on apparent clearance and a power covariate on apparent central volume (Wang 2020). The active metabolite SU012662 was fitted as a separate model (see Wang_2020_sunitinib_su12662)."
  reference <- paste(
    "Wang E, DuBois SG, Wetmore C, Khosravan R.",
    "Population pharmacokinetics-pharmacodynamics of sunitinib in pediatric",
    "patients with solid tumors.",
    "Cancer Chemother Pharmacol. 2020;86(2):181-192.",
    "doi:10.1007/s00280-020-04106-z.",
    sep = " "
  )
  vignette <- "Wang_2020_sunitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    BSA = list(
      description = "Baseline body surface area (DuBois and DuBois formula)",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline value; reference 1.47 m^2 (cohort median, Wang 2020 Table 2 and Table 3 footnote c). Linear effect on CL/F, power effect on Vc/F. BSA was computed as 0.20247 * height(m)^0.725 * weight(kg)^0.425 (Methods).",
      source_name = "BSA"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F and Vc/F in a separate SCM run with BSA replaced by body weight; the BSA model had the lower OFV and was selected (Wang 2020 Results).",
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F and Vc/F (Table 1); not retained.",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Tested on CL/F and Vc/F (Table 1); not retained.",
      source_name = "SEX"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = "Tested on CL/F (Table 1); not retained.",
      source_name = "RACE"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "sunitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sunitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "sunitinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 59L,
    n_studies = 2L,
    age_range = "2-21 years",
    weight_range = "16.2-100 kg (median 50.4 kg)",
    bsa_range = "0.66-2.14 m^2 (median 1.47 m^2)",
    sex_female_pct = 52.5,
    race_ethnicity = c(Asian = 5.1, NonAsian = 89.8, Unknown = 5.1),
    disease_state = "Children and young adults with refractory solid tumours, predominantly high-grade glioma, ependymoma, brain stem glioma, or sarcoma (Children's Oncology Group studies ADVL0612 and ACNS1021).",
    dose_range = "Sunitinib 15 or 20 mg/m^2 orally once daily on schedule 4/2 (4 weeks on, 2 weeks off); intact capsules or capsule contents sprinkled on yogurt or apple sauce.",
    regions = "United States and Canada (Children's Oncology Group).",
    notes = "Baseline demographics from Wang 2020 Table 2 (28 male / 31 female; 3 Asian, 53 non-Asian, 3 unknown). 365 sunitinib plasma observations; LLOQ 1 ng/mL (ADVL0612) or 0.1 ng/mL (ACNS1021)."
  )

  ini({
    # Structural parameters (Wang 2020 Table 3, sunitinib final model, typical patient with BSA 1.47 m^2)
    lka <- log(0.38); label("Absorption rate constant ka (1/h)") # Table 3 sunitinib 'ka' = 0.38 (RSE 31.5%)
    ltlag <- log(0.64); label("Absorption lag time tlag (h)") # Table 3 sunitinib 'tlag' = 0.64 (RSE 12.8%)
    lcl <- log(24.1); label("Apparent clearance CL/F at BSA 1.47 m^2 (L/h)") # Table 3 sunitinib 'CL/F' = 24.1 (RSE 6.6%); Results equation CL/F = 24.1 * [1 + 0.557 * (BSA - 1.47)]
    lvc <- log(1070); label("Apparent central volume Vc/F at BSA 1.47 m^2 (L)") # Table 3 sunitinib 'Vc/F' = 1070 (RSE 10.3%)
    lvp <- log(63.8); label("Apparent peripheral volume Vp/F (L)") # Table 3 sunitinib 'Vp/F' = 63.8 (RSE 34.3%)
    lq <- log(0.28); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 sunitinib 'Q/F' = 0.28 (RSE 177%)

    # Covariate effects
    e_bsa_cl <- 0.557; label("Linear slope of BSA on CL/F, per m^2 about 1.47 m^2 (unitless)") # Results text CL/F = 24.1 * [1 + 0.557 * (BSA - 1.47)]; Table 3 'BSA on CL/F' = 0.56 (RSE 19.9%)
    e_bsa_vc <- 1.47; label("Power exponent of BSA/1.47 on Vc/F (unitless)") # Results text Vc/F = 1070 * (BSA/1.47)^1.47; Table 3 'BSA on Vc/F' = 1.47 (RSE 19.6%)

    # IIV: Table 3 reports omega as CV%; omega^2 = log(1 + CV^2)
    etalcl ~ 0.1106 # Table 3 sunitinib 'omega (CL/F)' = 34.2%
    etalvc ~ 0.05646 # Table 3 sunitinib 'omega (Vc/F)' = 24.1%
    etalka ~ 0.5705 # Table 3 sunitinib 'omega (ka)' = 87.7%

    # Residual error
    propSd <- 0.319; label("Proportional residual error (fraction)") # Table 3 sunitinib 'sigma' = 31.9% (RSE 2.7%)
  })

  model({
    # Individual PK parameters with BSA covariate effects
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag)
    cl <- exp(lcl + etalcl) * (1 + e_bsa_cl * (BSA - 1.47))
    vc <- exp(lvc + etalvc) * (BSA / 1.47)^e_bsa_vc
    vp <- exp(lvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # mg / L * 1000 = ng/mL
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
