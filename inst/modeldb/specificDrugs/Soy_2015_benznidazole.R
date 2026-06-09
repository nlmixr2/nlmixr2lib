Soy_2015_benznidazole <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and",
    "first-order elimination for oral benznidazole in adult patients",
    "with chronic Chagas disease (Soy 2015; CINEBENZ trial, n = 39",
    "index plus n = 10 external validation). Apparent clearance CL/F =",
    "1.73 L/h, apparent volume of distribution V/F = 89.6 L, and",
    "absorption rate constant Ka = 1.15 1/h fixed from the literature",
    "(Raaflaub & Ziegler 1979). Inter-individual variability is on CL/F",
    "(33.4% CV) and V/F (68.8% CV); inter-occasion variability is on",
    "CL/F (29.5% CV), folded into the CL/F eta as BSV-equivalent for",
    "forward simulation. Residual error is combined additive (0.57",
    "mg/L) plus proportional (19.53% CV). No demographic or biological",
    "covariates were retained in the final model."
  )
  reference <- paste(
    "Soy D, Aldasoro E, Guerrero L, Posada E, Serret N, Mejia T,",
    "Urbina JA, Gascon J. Population pharmacokinetics of benznidazole",
    "in adult patients with Chagas disease.",
    "Antimicrobial Agents and Chemotherapy. 2015;59(6):3342-3349.",
    "doi:10.1128/AAC.05018-14.",
    sep = " "
  )
  vignette <- "Soy_2015_benznidazole"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  # Soy 2015 screened the covariates below and did not retain any of
  # them in the final structural model (Results: "none of the
  # covariates ... significantly influenced the BNZ CL/F"). They are
  # listed here for provenance only and are not referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL/F and V/F (Soy 2015 Results, covariate model selection paragraph). Not retained.",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Screened on CL/F and V/F (Soy 2015 Results, covariate model selection paragraph). Not retained. Index cohort was 26/39 (66.7%) female.",
      source_name        = "GENDER"
    ),
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL/F and V/F via the scm approach (Soy 2015 Results). Not retained. Separately, Soy 2015 ran an exploratory WT-popPK model with allometric exponents fixed at 0.75 (CL/F) and 1 (V/F) and reference weight = cohort mean; that model did not improve fit and the typical-value estimates were essentially unchanged from the final no-covariate model (CL/F = 1.75 L/h, V/F = 95.3 L). Reported here for completeness; not encoded in the final model().",
      source_name        = "WT"
    ),
    BMI = list(
      description        = "Body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL/F and V/F (Soy 2015 Results). Not retained.",
      source_name        = "BMI"
    ),
    CRCL = list(
      description        = "Creatinine clearance (Cockcroft-Gault)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL/F and V/F (Soy 2015 Results). Not retained. Computed via Cockcroft-Gault (Soy 2015 Methods, Data section).",
      source_name        = "CLCR"
    ),
    TPRO = list(
      description        = "Total serum protein",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL/F and V/F (Soy 2015 Results). Not retained.",
      source_name        = "TPRO"
    ),
    TBILI = list(
      description        = "Total bilirubin",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL/F and V/F (Soy 2015 Results). Not retained.",
      source_name        = "TBILI"
    ),
    HCT = list(
      description        = "Hematocrit",
      units              = "fraction",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL/F and V/F (Soy 2015 Results); reported as a fraction in Table 1 (index data set mean 0.40, SD 0.04). Not retained.",
      source_name        = "HCT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 39L,
    n_studies      = 1L,
    n_observations = 358L,
    age_range      = "19-55 years",
    age_mean       = "37.15 years (median 36 years)",
    weight_range   = "43-100 kg",
    weight_mean    = "70.55 kg",
    sex_female_pct = round(100 * 26 / 39, 1),
    race_ethnicity = c(Bolivian_pct = 96),
    disease_state  = paste(
      "Adult patients (>=18 years) with chronic Chagas disease",
      "diagnosed by two T. cruzi serologic tests; treatment-naive",
      "at enrolment. Exclusion criteria included hypersensitivity to",
      "benznidazole, immunodeficiency (including HIV), hepatic or",
      "renal impairment, pregnancy, and lactation."
    ),
    dose_range     = paste(
      "Oral benznidazole 2.5 mg/kg every 12 hours for 8 weeks",
      "(maximum 400 mg/day; Abarax, Elea Laboratory, Argentina).",
      "Median drug adherence across the cohort 99.2%."
    ),
    sampling       = paste(
      "Plasma samples on treatment days +1 (first dose) and +15",
      "(steady state) at 1 h, 3-6 h, and 6-12 h post-dose; after the",
      "final dose at 3-12 h, 12-24 h, and 24-36 h; additional trough",
      "samples at days +30 and +45. Mean 9.1 samples per subject in",
      "the index data set. HPLC-UV assay; LLOQ = 1.6 mg/L."
    ),
    validation_cohort = paste(
      "An additional external validation cohort of 10 patients (96",
      "plasma samples, mean 9.6/subject) drawn under conditions",
      "similar to the index data set. Index cohort fit (n = 39)",
      "feeds the popPK parameters; validation cohort is used only",
      "for external predictive-performance assessment, not for",
      "parameter estimation."
    ),
    regions        = paste(
      "Single center: Hospital Clinic of Barcelona, Spain.",
      "96% of patients were originally from Bolivia."
    ),
    notes          = paste(
      "CINEBENZ study, EudraCT 2011-002900-34, ClinicalTrials.gov",
      "NCT01755403. Baseline demographics from Soy 2015 Table 1",
      "(index data set column); final-model parameter estimates",
      "from Soy 2015 Table 2."
    )
  )

  ini({
    # Structural PK parameters (Soy 2015 Table 2, base/final popPK
    # column). Reported as apparent values (CL/F, V/F) because oral
    # bioavailability F could not be estimated without an intravenous
    # reference formulation (Discussion paragraph 3 and Results,
    # population PK model paragraph). Ka was fixed at 1.15 1/h from
    # Raaflaub & Ziegler 1979 (Soy 2015 reference 8); the available
    # sampling design did not support its estimation.
    lka <- fixed(log(1.15)); label("Absorption rate constant Ka (1/h; fixed from Raaflaub & Ziegler 1979)")  # Soy 2015 Table 2 ("Ka 1.15 (fixed)")
    lcl <- log(1.73);        label("Apparent oral clearance CL/F (L/h)")                                     # Soy 2015 Table 2 (CL/F)
    lvc <- log(89.6);        label("Apparent volume of distribution V/F (L)")                                # Soy 2015 Table 2 (V/F)

    # IIV on log-normal scale; omega^2 = log(1 + CV^2). Soy 2015
    # reports IIV CL/F = 33.4% CV, IIV V/F = 68.8% CV, and IOV CL/F
    # = 29.5% CV (Table 2). Following the precedent in
    # Denti_2018_levofloxacin.R, the IOV component on CL/F is folded
    # into etalcl as a BSV-equivalent for forward simulation since
    # nlmixr2lib has no idiomatic encoding for per-occasion
    # variability separate from BSV. The total variance on log-CL/F
    # is therefore log(1 + 0.334^2) + log(1 + 0.295^2) = 0.10579 +
    # 0.08344 = 0.18923 (geometric CV% ~= 45.6%).
    etalcl ~ log(1 + 0.334^2) + log(1 + 0.295^2)  # Soy 2015 Table 2: IIV CL/F = 33.4% CV plus IOV CL/F = 29.5% CV, folded as BSV-equivalent
    etalvc ~ log(1 + 0.688^2)                      # Soy 2015 Table 2: IIV V/F = 68.8% CV

    # Residual error: combined additive + proportional on linear
    # concentrations. Table 2 reports sigma^2_1 = 0.57 (additive SD,
    # mg/L) and sigma^2_2 = 19.53 (proportional, expressed as CV%).
    # The Results text (population pharmacokinetic model paragraph)
    # quotes the proportional component as 19.1% CV; Table 2 is the
    # authoritative value used here.
    addSd  <- 0.57;   label("Additive residual error (mg/L)")             # Soy 2015 Table 2 (sigma^2_1)
    propSd <- 0.1953; label("Proportional residual error (fraction CV)")  # Soy 2015 Table 2 (sigma^2_2)
  })

  model({
    # Individual PK parameters with log-normal random effects.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # One-compartment ODE system with first-order absorption from
    # the depot and first-order elimination from the central
    # compartment.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Concentration units: dose mg / volume L = mg/L, matching the
    # HPLC-UV assay (LLOQ 1.6 mg/L; analytical range 1.6-100 mg/L).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
