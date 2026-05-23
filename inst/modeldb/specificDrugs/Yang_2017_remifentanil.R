Yang_2017_remifentanil <- function() {
  description <- "One-compartment population PK model for continuous intravenous remifentanil infusion in critically ill adults receiving venoarterial extracorporeal membrane oxygenation (VA-ECMO), with sex and centrifugal-pump rotational speed as covariates on clearance (Yang 2017)."
  reference <- paste(
    "Yang S, Noh H, Hahn J, Jin BH, Min KL, Bae SK, Kim J,",
    "Park MS, Hong T, Wi J, Chang MJ.",
    "Population pharmacokinetics of remifentanil in critically ill patients",
    "receiving extracorporeal membrane oxygenation.",
    "Sci Rep 2017;7(1):16275. doi:10.1038/s41598-017-16358-6.",
    sep = " "
  )
  vignette <- "Yang_2017_remifentanil"
  units <- list(
    time          = "h",
    dosing        = "ug",
    concentration = "ng/mL"
  )
  # Dose units ug (i.e. mg/h continuous-infusion rates converted to ug/h
  # in the event table); central / vc gives ug/L = ng/mL directly,
  # matching the Yang 2017 plasma concentration units (Table 2 footnote;
  # Methods, Remifentanil concentration assay).

  covariateData <- list(
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male) under the canonical SEXF convention",
      notes              = "Yang 2017 encodes sex as a male-indicator (SEX = 0 for female, SEX = 1 for male) with female as the structural reference category (the published equation CL = 366 x 0.502^SEX yields CL = 366 L/h for SEX = 0 / female). To store under the canonical SEXF (1 = female, 0 = male) while preserving Yang's female-reference CL of 366 L/h, the effect is applied in model() as exp(e_sex_cl * (1 - SEXF)), so SEXF = 1 (female) yields factor 1 and SEXF = 0 (male) yields the paper's male-vs-female log-coefficient log(0.502) (about -0.689). This mirrors the Bajaj 2017 nivolumab pattern. The published cohort had 10 male and 5 female patients (Table 1; 67% male).",
      source_name        = "SEX"
    ),
    ECMO_PUMP_SPEED = list(
      description        = "Extracorporeal-membrane-oxygenation centrifugal-pump rotational speed",
      units              = "RPM",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL: CL_i = CL * (ECMO_PUMP_SPEED / 2350)^2.04 per Yang 2017 Results (final-model equation). Reference 2350 RPM is the cohort median (Results: 'median ECMO pump speeds of 2350 RPM', IQR 2302-2532 RPM). Treated as time-fixed per subject in Yang 2017 -- the per-subject pump speed used was the prevailing speed during the PK sampling window. The reported pump-speed range explored in the simulation analyses (Methods, Simulations) was 1700-2900 RPM.",
      source_name        = "ECMO pump speed"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 15L,
    n_studies        = 1L,
    age_range        = "19-78 years",
    age_median       = "57 years (IQR 45-69)",
    weight_range     = "40.8-94.0 kg",
    weight_median    = "65.4 kg (IQR 54.5-70.0)",
    sex_female_pct   = 33.3,
    bmi_median       = "23.8 kg/m^2 (IQR 21.2-24.2)",
    total_protein_median = "4.7 g/dL (IQR 3.9-5.3)",
    disease_state    = "Critically ill adults receiving venoarterial extracorporeal membrane oxygenation (VA-ECMO) in a cardiovascular intensive care unit. Indications for VA-ECMO included acute myocardial infarction (n = 5), non-ST-segment elevation MI (n = 4), STEMI (n = 1), ischemic cardiomyopathy (n = 1), pulmonary embolism (n = 1), coronary artery occlusive disease (n = 1), myocarditis (n = 1), atrial fibrillation with bronchiolitis (n = 1), and angina pectoris (n = 1) per Table 1.",
    dose_range       = "Continuous intravenous remifentanil infusion at a median rate of 0.35 mg/h (IQR 0.25-0.35 mg/h; full per-patient range 0.14-1.0 mg/h per Table 1). No patient received a bolus injection during the study.",
    ecmo_duration    = "Median 143 h on VA-ECMO (IQR 96-250; range 48-532 h)",
    ecmo_pump_speed  = "Median 2350 RPM (IQR 2302-2532)",
    crrt_pct         = 66.7,
    regions          = "South Korea (Seoul; Severance Cardiovascular Hospital)",
    notes            = "Single-center prospective cohort study at Severance Cardiovascular Hospital, Yonsei University, between January 2015 and December 2016 (ClinicalTrials.gov NCT02581280). 55 remifentanil plasma concentrations across 15 patients (at least 3 samples per patient). PK samples drawn from dwelling arterial lines at 8-12 h (T1), 24 h (T2), and 36-48 h (T3) of remifentanil infusion during VA-ECMO support; if remifentanil was discontinued during VA-ECMO, serial samples were collected immediately before discontinuation and at 5, 10, 15, 25, 30, and 40 min. Plasma quenched on ice with formic acid then stored at -80 C until assay. Concentrations measured by LC-MS/MS validated 0.05-500 ng/mL, LLOQ 0.05 ng/mL. ECMO circuit: Terumo Capiox SP centrifugal pump with PMEA-coated polymethylpentene Capiox EBS oxygenator and PVC tubing. Baseline demographics in Table 1. NONMEM 7.3 with FOCE INTER; bootstrap n = 5000."
  )

  ini({
    # Structural parameters from Yang 2017 Table 2 final population PK
    # estimates and Results final-model equation. Reference subject:
    # female (SEXF = 1) at the cohort-median ECMO pump speed 2350 RPM.
    # Table 2 reports CL in L/h and V in L.
    lcl <- log(366); label("Clearance for reference female at 2350 RPM (CL, L/h)")    # Table 2: TVCL = 366 L/h
    lvc <- log(41);  label("Central volume of distribution (V, L)")                   # Table 2: TVV  = 41 L

    # Covariate effects from Yang 2017 Results final-model equation:
    # CL (L/h) = 366 x 0.502^SEX x (ECMO pump speed / 2350)^2.04
    # where SEX = 0 (female) / SEX = 1 (male). Stored as a log
    # coefficient applied via exp(e_sex_cl * (1 - SEXF)) so that
    # SEXF = 1 (female) reproduces the paper's female-reference CL of
    # 366 L/h and SEXF = 0 (male) gives CL x 0.502.
    e_sex_cl       <- log(0.502); label("Log multiplier of male sex on CL (unitless; applied as exp((1 - SEXF)*e_sex_cl), reproduces paper's 0.502^SEX)") # Table 2: theta_sex on CL = 0.502
    e_pumpspeed_cl <- 2.04;       label("Power exponent of (ECMO_PUMP_SPEED / 2350) on CL (unitless)")                                                    # Table 2: theta_ECMOpumpspeed on CL = 2.04

    # Inter-individual variability. Yang 2017 Methods: "an exponential
    # variance model was evaluated and assumed to follow a log-normal
    # distribution with a mean of zero and a variance of omega^2". Table 2
    # column header is "omega CL" (no superscript), reporting the
    # square root of omega^2, i.e. the SD on the log scale; the
    # additive sigma column carries linear ng/mL units (Table 2) which
    # only resolves dimensionally if reported as SD, not variance.
    # Square the SD to obtain the variance used by nlmixr2's `~`.
    # IIV on V was not estimated (Results: "The estimate of the
    # interindividual variability on V was near zero likely because
    # of the narrow weight range in the patient population").
    etalcl ~ 0.0154  # Table 2: omega_CL = 0.124 (SD on log scale); 0.124^2 = 0.015376

    # Residual error. Yang 2017 Methods: residual variability assumed
    # normally distributed with mean zero and variance sigma^2. Table 2
    # reports "sigma proportional" and "sigma additive" with linear
    # additive units (ng/mL), consistent with SD reporting.
    propSd <- 0.387; label("Proportional residual SD (fraction)")   # Table 2: sigma proportional = 0.387
    addSd  <- 0.111; label("Additive residual SD (ng/mL)")          # Table 2: sigma additive    = 0.111 ng/mL
  })
  model({
    # Derived sex term: Yang 2017 encodes sex as a male-indicator with
    # female as the structural reference category, so (1 - SEXF)
    # reproduces the paper's SEX = 1 (male) column while keeping SEXF
    # (1 = female) as the canonical storage convention. This mirrors
    # the Bajaj 2017 nivolumab pattern.
    sex_male <- 1 - SEXF

    # Individual clearance per Yang 2017 Results final-model equation:
    #   CL (L/h) = 366 * 0.502^SEX * (ECMO pump speed / 2350)^2.04
    # Reference subject: female (sex_male = 0) at the cohort-median
    # ECMO pump speed of 2350 RPM.
    cl <- exp(lcl + etalcl) *
      exp(e_sex_cl * sex_male) *
      (ECMO_PUMP_SPEED / 2350)^e_pumpspeed_cl

    # No IIV on volume (Yang 2017 Results: not estimated due to the
    # narrow weight range in the cohort).
    vc <- exp(lvc)

    kel <- cl / vc

    # One-compartment open model with linear first-order elimination
    # and zero-order infusion input (Results: "one-compartment model
    # with a zero-order input and first-order output (linear
    # elimination)"). Continuous infusion enters the central
    # compartment directly via the event table's rate column; no depot
    # / first-order absorption is modelled because no patient received
    # a bolus injection during the study (Methods, Study procedures).
    d/dt(central) <- -kel * central

    # Plasma concentration. Dose units ug (mg/h infusion rates
    # converted to ug/h in the event table; 1 mg/h = 1000 ug/h);
    # vc units L; central / vc has units ug/L = ng/mL, matching
    # the Yang 2017 plasma units (Table 2 footnote).
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
