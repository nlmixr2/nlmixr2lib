RoubaudBaudron_2025_ceftriaxone <- function() {
  description <- paste(
    "Two-compartment population PK model for ceftriaxone given intravenously",
    "or subcutaneously (1 g every 24 h) to patients aged 65 years or older",
    "(PhASAge study, France). Subcutaneous doses enter a depot with",
    "first-order absorption and a bioavailability of 99%; intravenous doses",
    "enter the central compartment. Clearance is additive-linear in",
    "Cockcroft-Gault creatinine clearance (an intercept plus a slope on",
    "CRCL / 59.66) with an additive (normal, L/h) between-subject random",
    "effect; the absorption rate and central volume carry log-normal",
    "random effects. Proportional residual error. Parametric model fitted",
    "in Monolix 2023-R1 (SAEM)."
  )
  reference <- paste(
    "Roubaud-Baudron C, Fauchon H, Stanke-Labesque F, Paccalin M, Breilh D,",
    "Gregoire N, Forestier E, Ferry T, Gavazzi G, Goutelle S; the",
    "Infectiology and Geriatrics InterGroup (SPILF-SFGG GInGer).",
    "Pharmacokinetics of Subcutaneous and Intravenous Ceftriaxone in an Older",
    "Population: The PhASAge Study. Open Forum Infect Dis.",
    "2025;12(6):ofaf313. doi:10.1093/ofid/ofaf313."
  )
  vignette <- "RoubaudBaudron_2025_ceftriaxone"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  dosing <- c("depot", "central")

  # The between-subject random effect on clearance is additive on the linear
  # (L/h) scale: Results state 'clearance was normally distributed in the
  # final model' and Table 2 footnote c prints 'CLi = 0.35 + 0.3 x
  # (CLcri / 59.66) + etaiCL'. It acts on the linear clearance rather than on
  # the log-scale intercept lcl_int (precedents: Balice 2022 daptomycin,
  # Lin 2021 vancomycin).
  paper_specific_etas <- c("etacl")

  compartmentData <- list(
    depot = list(analyte = "ceftriaxone", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ceftriaxone", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ceftriaxone", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault formula (raw, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source symbol CLcr; Table 1 row 'GFRe, mean (SD), mL/min' with",
        "footnote b 'Cockcroft-Gault formula' (cohort mean 57.5, SD 26.5",
        "mL/min). Enters clearance additively as 0.3 x (CRCL / 59.66) L/h",
        "(Table 2 footnote c), so the intercept 0.35 L/h is the clearance at",
        "CRCL = 0 and the typical clearance at CRCL = 59.66 mL/min is",
        "0.65 L/h. The 59.66 mL/min divisor is printed in the equation without",
        "being identified as a cohort statistic. Treated as time-fixed per",
        "subject. For its dosing simulations the paper sampled CRCL from a",
        "log-normal distribution with median 51.95 mL/min and log-scale SD",
        "0.51 (Supplementary Table S3 footnote).",
        collapse = " "
      ),
      source_name = "CLcr"
    )
  )

  # Covariates screened (Methods: Monolix stepwise covariate search;
  # Results: 'No other covariate significantly improved the fit') but not
  # retained. Documentation only -- none enters model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Named in the Results as a variable that could influence ceftriaxone PK; not retained in the final model. Cohort mean 67.4 (SD 14.3) kg, Table 1."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Named in the Results as a variable that could influence ceftriaxone PK; not retained in the PK model. Cohort mean 31.1 (SD 4.8) g/L, Table 1 (9 missing). Albumin enters the paper's dosing simulations only through the Sanz-Codina 2023 unbound-fraction equation, which is not part of this model."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Not retained. Cohort mean 85.5 (SD 7.6) years, Table 1; enrolment floor 65 years."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 47L,
    n_studies = 1L,
    n_observations = "3 samples per intravenous patient (C0, Cend, C5h) and 4 per subcutaneous patient (C0, Cend, C2h, C5h) at steady state",
    age_range = "65 years or older; mean 85.5 (SD 7.6) years",
    weight_range = "mean 67.4 (SD 14.3) kg",
    sex_female_pct = 51,
    race_ethnicity = "not reported (French multicentre cohort)",
    disease_state = "Hospitalised older patients (geriatric wards) treated with ceftriaxone for infection, mainly urinary (56%) or pulmonary (23%); 77% functionally dependent (Katz ADL < 6), 96% polypharmacy",
    dose_range = "Ceftriaxone 1 g every 24 h for more than 48 h (assumed steady state); intravenous (n = 23) or subcutaneous infusion (n = 24; median 50 mL over a median 30 min, mostly abdomen or thigh)",
    regions = "France (5 hospitals), August 2019 to February 2022",
    renal_function = "Cockcroft-Gault creatinine clearance mean 57.5 (SD 26.5) mL/min (Table 1)",
    albumin = "mean 31.1 (SD 4.8) g/L (Table 1)",
    notes = "PhASAge study (NCT03583749; EudraCT 2018-001432-21), ceftriaxone substudy. 50 patients enrolled, 47 with PK samples analysed. A nonparametric Pmetrics model of the same structure is reported in Supplementary Tables S1-S2 and is not encoded here."
  )

  ini({
    # Final parametric model, Roubaud-Baudron 2025 Table 2 (Monolix 2023-R1).
    lka <- log(1.63); label("Subcutaneous absorption rate constant (1/h)") # Table 2 Ka 1.63 h-1 (RSE 15.5%)
    lvc <- log(7.13); label("Central volume of distribution (L)") # Table 2 Vc 7.13 L (RSE 8.83%)
    lvp <- log(8.29); label("Peripheral volume of distribution (L)") # Table 2 Vp 8.29 L (RSE 0.51%)
    lq <- log(9.14); label("Intercompartmental clearance (L/h)") # Table 2 Q 9.14 L/h (RSE 0.29%)
    lcl_int <- log(0.35); label("Intercept of the additive-linear CL equation: CL at CRCL = 0 (L/h)") # Table 2 CL 0.35 L/h (RSE 15.2%); footnote c equation
    e_crcl_cl <- 0.3; label("Additive CL change per unit of CRCL / 59.66 (L/h)") # Table 2 'beta Cockcroft' 0.3 (RSE 16.8%); footnote c equation
    lfdepot <- log(0.99); label("Subcutaneous bioavailability (fraction)") # Results 'Ceftriaxone Concentrations': 'bioavailability of subcutaneous ceftriaxone was estimated at 99% by population modeling' (not tabulated in Table 2)

    # Table 2 'Random effects, SD of CV, %' prints 'SD-CV%' pairs. Ka and Vc
    # are log-normal (omega SD 0.56 and 0.46 -> variances 0.3136 and 0.2116;
    # sqrt(exp(omega^2) - 1) gives 60.7% and 48.6%, matching the printed CVs).
    # CL is normally distributed: SD 0.14 L/h (variance 0.0196 (L/h)^2;
    # 0.14 / 0.35 = 40%, printed CV 40.79%).
    etalka ~ 0.3136 # Table 2 Ka random effect SD 0.56 (CV 60.86%, RSE 21.9%)
    etalvc ~ 0.2116 # Table 2 Vc random effect SD 0.46 (CV 48.7%, RSE 15.9%)
    etacl ~ 0.0196 # Table 2 CL random effect SD 0.14 L/h (CV 40.79%, RSE 11.7%)

    propSd <- 0.11; label("Proportional residual error (fraction)") # Table 2 'Residual error, b' 0.11 (RSE 9.40%)
  })
  model({
    crcl_ref <- 59.66 # mL/min, Table 2 footnote c divisor

    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q <- exp(lq)
    # Table 2 footnote c: CLi = 0.35 + 0.3 x (CLcri / 59.66) + etaiCL (L/h)
    cl <- exp(lcl_int) + e_crcl_cl * (CRCL / crcl_ref) + etacl

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot)

    # Dose in mg, vc in L: central / vc is mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
