Toyoshima_2021_peficitinib <- function() {
  description <- paste(
    "Two-compartment population PK model for oral peficitinib (a pan-Janus",
    "kinase inhibitor) in adult Asian patients with rheumatoid arthritis",
    "(Toyoshima 2021 RA patient model; 989 patients with PK data from the",
    "phase 2 RAJ1 and phase 3 RAJ3 / RAJ4 studies). Absorption is sequential",
    "zero-order then first-order with a lag time: the dose enters the depot",
    "over a zero-order duration D after the lag ALAG and is then absorbed",
    "first-order at rate Ka. Apparent clearance CL/F depends on baseline",
    "MDRD eGFR and baseline lymphocyte count through power functions",
    "centred on the RA-patient means (91.5 mL/min/1.73 m^2 and 1550 x 10^6",
    "cells/L). Interindividual variability on CL and Vc only; proportional",
    "residual error. Structural parameters were estimated with the NONMEM",
    "PRIOR NWPRI penalty informed by the companion healthy-volunteer model",
    "(Toyoshima_2021_peficitinib_healthy).",
    sep = " "
  )
  reference <- paste(
    "Toyoshima J, Shibata M, Kaibara A, Kaneko Y, Izutsu H, Nishimura T.",
    "(2021). Population pharmacokinetic analysis of peficitinib in patients",
    "with rheumatoid arthritis. Br J Clin Pharmacol 87(4):2014-2022.",
    "doi:10.1111/bcp.14605.",
    "This file encodes the final RA patient model (Table 3, Equation 5);",
    "the prior healthy-volunteer model of Supplemental Table 2 is encoded",
    "in Toyoshima_2021_peficitinib_healthy.R.",
    sep = " "
  )
  vignette <- "Toyoshima_2021_peficitinib"
  # Doses in mg and volumes in L give mg/L; the observation is scaled by 1000
  # to ng/mL, the unit of the assay LLOQ (0.25 ng/mL), Figure 2 and Table 4.
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "peficitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "peficitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "peficitinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Baseline estimated glomerular filtration rate (MDRD equation)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline eGFR calculated with the modification of diet in renal",
        "disease (MDRD) equation (Methods 2.4). Enters CL as a power",
        "function centred on the RA-patient arithmetic mean:",
        "(CRCL / 91.5)^0.213 (Equation 5; Table 3 'eGFR on CL' = 0.213).",
        "RA-patient mean 91.49 (SD 22.27), range 36.4-188.4 (Table 2).",
        "Patients with eGFR <= 40 were excluded from the phase 2/3 studies.",
        sep = " "
      ),
      source_name = "eGFR"
    ),
    LYMPH_ABS = list(
      description = "Baseline absolute peripheral-blood lymphocyte count",
      units = "cells/uL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The source reports the count in 10^6 cells/L, which is numerically",
        "identical to cells/uL, so no value transformation is needed. Enters",
        "CL as a power function centred on the RA-patient arithmetic mean:",
        "(LYMPH_ABS / 1550)^-0.104 (Equation 5; Table 3 'LYM on CL' =",
        "-0.104). RA-patient mean 1550 (SD 540), range 500-4600 (Table 2).",
        sep = " "
      ),
      source_name = "lymphocyte count (LYM)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on CL in the stepwise search (forward p < .01, backward p < .001) and not retained (Methods 2.4, Results 3.3).",
      source_name = "age"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened on CL and on Vc and not retained (Methods 2.4, Results 3.3).",
      source_name = "weight"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened on CL and not retained. RA-patient mean 40.1 g/L (Supplemental Table 1).",
      source_name = "ALB"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL and not retained.",
      source_name = "ALT"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL and not retained.",
      source_name = "AST"
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL and not retained.",
      source_name = "ALP"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened on CL and not retained.",
      source_name = "TBIL"
    ),
    NEUT = list(
      description = "Absolute neutrophil count",
      units = "cells/uL",
      type = "continuous",
      notes = "Screened on CL and on Vc and not retained. Reported in 10^6 cells/L (numerically cells/uL); RA-patient mean 5530 (Table 2).",
      source_name = "ANC"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Screened on CL and not retained.",
      source_name = "serum creatinine"
    ),
    CRP = list(
      description = "C-reactive protein",
      units = "mg/L",
      type = "continuous",
      notes = "Screened on CL and on Vc and not retained. RA-patient mean 24.45 mg/L (Table 2).",
      source_name = "CRP"
    ),
    HCT = list(
      description = "Haematocrit",
      units = "(fraction)",
      type = "continuous",
      notes = "Screened on CL and not retained. Reported as a fraction (mean 0.370; Supplemental Table 1).",
      source_name = "haematocrit"
    ),
    HGB = list(
      description = "Haemoglobin",
      units = "g/L",
      type = "continuous",
      notes = "Screened on CL and not retained.",
      source_name = "haemoglobin"
    ),
    PLT = list(
      description = "Platelet count",
      units = "10^9 cells/L",
      type = "continuous",
      notes = "Screened on CL and not retained.",
      source_name = "platelets"
    ),
    RBC = list(
      description = "Red blood cell count",
      units = "10^12 cells/L",
      type = "continuous",
      notes = "Screened on CL and not retained.",
      source_name = "red blood cells"
    ),
    URATE = list(
      description = "Serum urate",
      units = "umol/L",
      type = "continuous",
      notes = "Screened on CL and not retained. RA-patient mean 277.1 umol/L (Supplemental Table 1).",
      source_name = "urate"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL and not retained. 74.2% of RA patients were female (Table 2).",
      source_name = "sex"
    ),
    REGION = list(
      description = "Region of enrolment (Japan, Korea or Taiwan)",
      units = "(categorical)",
      type = "categorical",
      notes = "Screened on CL and not retained. 94.7% Japan, 3.0% Korea, 2.3% Taiwan (Table 2).",
      source_name = "region"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 989,
    n_studies = 3,
    age_range = "20-86 years (mean 55.3, SD 11.8; Table 2, n = 1011)",
    weight_range = "29.9-117.4 kg (mean 58.1, SD 12.4; Table 2, n = 1011)",
    sex_female_pct = 74.2,
    race_ethnicity = "Asian (enrolled in Japan, Korea and Taiwan)",
    disease_state = "Rheumatoid arthritis with inadequate response to conventional DMARDs (RAJ3) or methotrexate (RAJ4), or on monotherapy (RAJ1)",
    dose_range = "25, 50, 100 or 150 mg orally once daily, fed, in the morning",
    regions = "Japan 94.7%, Korea 3.0%, Taiwan 2.3%",
    renal_function = "MDRD eGFR 36.4-188.4 mL/min/1.73 m^2 (mean 91.49)",
    notes = paste(
      "4919 plasma concentrations from 989 patients in the phase 2 RAJ1",
      "(12 weeks) and phase 3 RAJ3 / RAJ4 (52 weeks) studies; mostly trough",
      "samples plus a post-dose sample at week 4 or 8 (Table 1). Table 2",
      "summarises 1011 RA patients. Patients receiving etanercept were",
      "excluded. Samples > 48 h after the last dose (all studies) or < 19 h",
      "after the last dose in RAJ1 were excluded, as were BLQ values",
      "(LLOQ 0.25 ng/mL) and outliers above mean + 3 SD of log trough",
      "concentration per arm (Methods 2.3). NONMEM 7.3, FOCE-I.",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- Toyoshima 2021 Table 3 'Estimate' column.
    # Apparent (oral) parameters; bioavailability is not separately estimated.
    lcl <- log(91.7)
    label("Apparent clearance CL/F at eGFR 91.5 and lymphocytes 1550 (L/h)") # Table 3 'CL (L/h)' 91.7 (RSE 2.3%); Equation 5
    lvc <- log(280)
    label("Apparent central volume Vc/F (L)") # Table 3 'Vc (L)' 280 (RSE 2.3%)
    lvp <- log(122)
    label("Apparent peripheral volume Vp/F (L)") # Table 3 'Vp (L)' 122 (RSE 8.5%)
    lq <- log(10.2)
    label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 'Q(L/h)' 10.2 (RSE 4.9%)
    lka <- log(5.83)
    label("First-order absorption rate constant (1/h)") # Table 3 'Ka (L/h)' 5.83 (RSE 7.6%); unit printed as L/h, a typo for 1/h (Supplemental Table 2 prints 1/h)
    ltlag <- log(0.132)
    label("Absorption lag time (h)") # Table 3 'ALAG (h)' 0.132 (RSE 2%)
    ld1 <- log(1.37)
    label("Duration of zero-order input into the depot (h)") # Table 3 'D (h)' 1.37 (RSE 4.3%)

    # Covariate effects -- power functions centred on the RA-patient means
    e_crcl_cl <- 0.213
    label("Power exponent of eGFR (CRCL / 91.5) on CL (unitless)") # Table 3 'eGFR on CL' 0.213 (RSE 19.4%); Equation 5
    e_lymph_abs_cl <- -0.104
    label("Power exponent of lymphocyte count (LYMPH_ABS / 1550) on CL (unitless)") # Table 3 'LYM on CL' -0.104 (RSE 26.9%); Equation 5

    # IIV -- Table 3 'Random effect for IIV (omega^2)'; footnote a: CV% = sqrt(exp(omega^2) - 1) x 100
    etalcl ~ 0.0639 # Table 3 omega^2 CL 0.0639 (25.7% CV, shrinkage 12.9%)
    etalvc ~ 0.143 # Table 3 omega^2 Vc 0.143 (39.2% CV, shrinkage 48.8%)

    # Residual error -- Table 3 'Residual error Proportional' 0.496; footnote c:
    # variability (49.6%) = estimate x 100, i.e. the estimate is on the SD scale.
    propSd <- 0.496
    label("Proportional residual error (fraction)") # Table 3 'Proportional' 0.496 (49.6%, RSE 1.6%)
  })
  model({
    # Individual parameters (Equation 1: P_i = theta x exp(eta_i); Equation 5 for CL)
    cl <- exp(lcl + etalcl) * (CRCL / 91.5)^e_crcl_cl * (LYMPH_ABS / 1550)^e_lymph_abs_cl
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q <- exp(lq)
    ka <- exp(lka)
    tlag <- exp(ltlag)
    d1 <- exp(ld1)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Sequential zero- then first-order absorption: the dose is released into
    # the depot over d1 after the lag (dose records need rate = -2).
    alag(depot) <- tlag
    dur(depot) <- d1

    # mg / L -> ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
