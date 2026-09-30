Toyoshima_2021_peficitinib_healthy <- function() {
  description <- paste(
    "Two-compartment population PK model for oral peficitinib (a pan-Janus",
    "kinase inhibitor) in healthy Japanese adult volunteers (Toyoshima 2021",
    "prior healthy-volunteer model; 98 subjects from five phase 1 studies",
    "at 150 mg). Absorption is sequential zero-order then first-order with",
    "a lag time: the dose enters the depot over a zero-order duration D",
    "after the lag ALAG and is then absorbed first-order at rate Ka.",
    "Interindividual variability on CL, Vp, Q, Ka, ALAG, D and relative",
    "bioavailability F (no IIV on Vc); proportional residual error; no",
    "covariates. The authors used this model as the NWPRI prior for the",
    "RA patient model (Toyoshima_2021_peficitinib).",
    sep = " "
  )
  reference <- paste(
    "Toyoshima J, Shibata M, Kaibara A, Kaneko Y, Izutsu H, Nishimura T.",
    "(2021). Population pharmacokinetic analysis of peficitinib in patients",
    "with rheumatoid arthritis. Br J Clin Pharmacol 87(4):2014-2022.",
    "doi:10.1111/bcp.14605.",
    "This file encodes the prior healthy-volunteer model of Supplemental",
    "Table 2; the final RA patient model of Table 3 is encoded in",
    "Toyoshima_2021_peficitinib.R.",
    sep = " "
  )
  vignette <- "Toyoshima_2021_peficitinib"
  # Doses in mg and volumes in L give mg/L; the observation is scaled by 1000
  # to ng/mL, the unit of the assay LLOQ (0.25 ng/mL) and Supplemental Figure 2.
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "peficitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "peficitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "peficitinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on CL; did not meet the forward-addition criterion (p < .01) so no covariate was retained in the healthy-volunteer model (Results 3.2).",
      source_name = "age"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened on CL and on Vp and not retained (Methods 2.4, Results 3.2).",
      source_name = "weight"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (MDRD)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "Screened on CL and not retained. Healthy-volunteer mean 93.85 (SD 13.95; Table 2).",
      source_name = "eGFR"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened on CL and not retained. Healthy-volunteer mean 44.2 g/L (Supplemental Table 1).",
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
    )
  )

  population <- list(
    species = "human",
    n_subjects = 98,
    n_studies = 5,
    age_range = "20-69 years (mean 34.7, SD 12.0; Table 2)",
    weight_range = "48.8-77.9 kg (mean 64.1, SD 7.4; Table 2)",
    sex_female_pct = 4.1,
    race_ethnicity = "Japanese",
    disease_state = "Healthy volunteers, including normal-function control subjects of the hepatic (PK10) and renal (PK11) impairment studies",
    dose_range = "150 mg orally, single dose (PK10, PK11, PK12, PK27) and single then once-daily multiple dose for 7 days (PK20)",
    regions = "Japan",
    renal_function = "MDRD eGFR 60.7-130.8 mL/min/1.73 m^2 (mean 93.85)",
    notes = paste(
      "2464 plasma concentrations from 98 subjects in five clinical",
      "pharmacology studies with rich sampling to 48-72 h (Table 1).",
      "Subjects with renal or hepatic impairment were excluded. Food",
      "condition (fasted, high-fat or normal fed) differed between phase 1",
      "studies but was not included in the model (Discussion).",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- Toyoshima 2021 Supplemental Table 2 'Estimate' column.
    # Apparent (oral) parameters. Typical bioavailability is not tabulated;
    # only its IIV is, so the typical value is the reference 1.
    lcl <- log(83.4)
    label("Apparent clearance CL/F (L/h)") # Supplemental Table 2 'CL (L/h)' 83.4 (RSE 3%); Discussion 'estimated CL was 83.4 ... L/h in the prior healthy volunteer ... model'
    lvc <- log(253)
    label("Apparent central volume Vc/F (L)") # Supplemental Table 2 'Vc (L)' 253 (RSE 3.3%)
    lvp <- log(124)
    label("Apparent peripheral volume Vp/F (L)") # Supplemental Table 2 'Vp (L)' 124 (RSE 11.4%)
    lq <- log(9.19)
    label("Apparent intercompartmental clearance Q/F (L/h)") # Supplemental Table 2 'Q (L/h)' 9.19 (RSE 9.2%)
    lka <- log(5.04)
    label("First-order absorption rate constant (1/h)") # Supplemental Table 2 'Ka (1/h)' 5.04 (RSE 12.3%)
    ltlag <- log(0.133)
    label("Absorption lag time (h)") # Supplemental Table 2 'ALAG (h)' 0.133 (RSE 6.6%)
    ld1 <- log(1.3)
    label("Duration of zero-order input into the depot (h)") # Supplemental Table 2 'D (h)' 1.3 (RSE 4.7%)
    lfdepot <- fixed(log(1))
    label("Relative bioavailability (fraction)") # not estimated: Supplemental Table 2 lists IIV on F only; typical F is the reference 1

    # IIV -- Supplemental Table 2 'Random effect for IIV (omega^2)';
    # footnote a: CV% = sqrt(exp(omega^2) - 1) x 100
    etalcl ~ 0.0068 # Supplemental Table 2 omega^2 CL 0.0068 (8.3% CV)
    etalvp ~ 0.427 # Supplemental Table 2 omega^2 Vp 0.427 (73% CV)
    etalq ~ 0.27 # Supplemental Table 2 omega^2 Q 0.27 (55.7% CV)
    etalka ~ 1.21 # Supplemental Table 2 omega^2 Ka 1.21 (153.4% CV)
    etaltlag ~ 0.525 # Supplemental Table 2 omega^2 ALAG 0.525 (83.1% CV)
    etald1 ~ 0.249 # Supplemental Table 2 omega^2 D 0.249 (53.2% CV)
    etalfdepot ~ 0.0626 # Supplemental Table 2 omega^2 F 0.0626 (25.4% CV)

    # Residual error -- footnote c: variability (33.8%) = estimate x 100 (SD scale)
    propSd <- 0.338
    label("Proportional residual error (fraction)") # Supplemental Table 2 'Proportional' 0.338 (33.8%, RSE 3.4%)
  })
  model({
    # Individual parameters (Equation 1: P_i = theta x exp(eta_i))
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)
    d1 <- exp(ld1 + etald1)
    fdepot <- exp(lfdepot + etalfdepot)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Sequential zero- then first-order absorption: the dose is released into
    # the depot over d1 after the lag (dose records need rate = -2).
    f(depot) <- fdepot
    alag(depot) <- tlag
    dur(depot) <- d1

    # mg / L -> ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
