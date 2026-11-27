Bae_2019_vancomycin <- function() {
  description <- "Two-compartment IV population PK model for vancomycin in Korean adult inpatients receiving therapeutic drug monitoring (Bae 2019). Clearance follows a power function of Cockcroft-Gault creatinine clearance (reference 72 mL/min) and is replaced by a separate absolute clearance in patients on hemodialysis or continuous renal replacement therapy; peripheral volume scales linearly with body weight (reference 60 kg)."
  reference <- "Bae SH, Yim DS, Lee H, Park AR, Kwon JE, Sumiko H, Han S. Application of Pharmacometrics in Pharmacotherapy: Open-Source Software for Vancomycin Therapeutic Drug Management. Pharmaceutics. 2019;11(5):224. doi:10.3390/pharmaceutics11050224"
  vignette <- "Bae_2019_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Bae 2019 Table 1: mean 61.6 kg (range 30.0-126.7). Enters the peripheral volume linearly as V2 = theta3 * (WT/60) (Table 2), with no estimated exponent; 60 kg is the reference weight the Results text quotes for the typical V2.",
      source_name = "Weight"
    ),
    CRCL = list(
      description = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Bae 2019 Table 1 (footnote 1: calculated by the Cockcroft-Gault equation): mean 77.0 mL/min (range 4.57-279); the Discussion gives a median of 72.4 mL/min. Enters CL as (CRCL/72)^theta2 for subjects on neither hemodialysis nor CRRT (Table 2). No truncation or flooring of CRCL is described. Stored under the canonical CRCL column, which accepts raw mL/min when the source does not BSA-normalize, with the assay form recorded here. Not read for subjects with RRT_HEMODIAL_STATUS = 1 or RRT_CRRT_STATUS = 1.",
      source_name = "CLCR"
    ),
    RRT_HEMODIAL_STATUS = list(
      description = "Hemodialysis treatment-status indicator (1 = patient receiving hemodialysis at the time of vancomycin therapy, 0 = not)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no hemodialysis)",
      notes = "Bae 2019 Table 1: 20 of 220 patients received hemodialysis. Per-subject treatment status during vancomycin therapy; the paper does not model dialysis-session timing. A value of 1 REPLACES the CRCL-driven clearance with the absolute clearance CL_HD = 0.334 L/h (Table 2). Mutually exclusive with RRT_CRRT_STATUS: set at most one of the two to 1.",
      source_name = "HD"
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous renal replacement therapy treatment-status indicator (1 = patient receiving CRRT at the time of vancomycin therapy, 0 = not)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no CRRT)",
      notes = "Bae 2019 Table 1: 9 of 220 patients received CRRT (modality not specified). A value of 1 REPLACES the CRCL-driven clearance with the absolute clearance CL_CRRT = 0.716 L/h (Table 2). Mutually exclusive with RRT_HEMODIAL_STATUS: set at most one of the two to 1.",
      source_name = "CRRT"
    )
  )

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 220L,
    n_studies = 1L,
    n_concentrations = 1020L,
    age_range = "21-98 years",
    age_mean = "63 years",
    weight_range = "30.0-126.7 kg",
    weight_mean = "61.6 kg",
    sex_female_pct = 36.8,
    race_ethnicity = "Korean (single-centre cohort, Seoul St. Mary's Hospital)",
    disease_state = "Adult inpatients receiving intravenous vancomycin with therapeutic drug monitoring; 20 patients on hemodialysis and 9 on continuous renal replacement therapy.",
    dose_range = "Intravenous vancomycin per routine clinical care; individual regimens varied.",
    regions = "Republic of Korea",
    renal_function = "Cockcroft-Gault creatinine clearance mean 77.0 mL/min (range 4.57-279), median 72.4 mL/min, 25th percentile 40 mL/min; serum creatinine mean 1.7 mg/dL (range 0.20-13.3).",
    notes = "Bae 2019 Table 1 and Methods 2.2: retrospective TDM data, June 2015 to August 2017, 139 male / 81 female. Fit in NONMEM 7.4 with FOCE-I. Concentrations above 100 (printed as ng/mL, evidently mg/L) were excluded as outliers."
  )

  ini({
    # Structural parameters (Bae 2019 Table 2, 'Estimate' column).
    lcl <- log(2.82); label("Clearance at CRCL = 72 mL/min, no hemodialysis or CRRT (L/h)") # Table 2 theta1 = 2.82 (RSE 4.18%)
    e_crcl_cl <- 0.836; label("Power exponent on (CRCL/72) for clearance (unitless)") # Table 2 theta2 = 0.836 (RSE 6.89%)
    lcl_crrt <- log(0.716); label("Clearance in patients on continuous renal replacement therapy (L/h)") # Table 2 CL_CRRT = 0.716 (RSE 11.0%)
    lcl_dial <- log(0.334); label("Clearance in patients on hemodialysis (L/h)") # Table 2 CL_HD = 0.334 (RSE 11.9%)
    lvc <- log(31.8); label("Central volume of distribution V1 (L)") # Table 2 V1 = 31.8 (RSE 7.01%)
    lq <- log(11.7); label("Intercompartmental clearance Q (L/h)") # Table 2 Q = 11.7 (RSE 7.42%)
    lvp <- log(75.4); label("Peripheral volume of distribution V2 at WT = 60 kg (L)") # Table 2 theta3 = 75.4 (RSE 7.91%)

    # Inter-individual variability (Table 2, reported as %CV with
    # CV = sqrt(exp(omega^2) - 1) x 100 per Methods Eq. 3), so
    # omega^2 = log(1 + CV^2). No IIV is reported on V1 or Q.
    etalcl ~ 0.685147 # Table 2 omega_CL = 99.2% CV -> log(1 + 0.992^2)
    etalvp ~ 0.216775 # Table 2 omega_V2 = 49.2% CV -> log(1 + 0.492^2)

    # Residual error: Table 2 reports only a proportional component.
    propSd <- 0.253; label("Proportional residual error (fraction)") # Table 2 sigma_prop = 0.253 (RSE 2.91%)
  })
  model({
    # Clearance: CRCL power model off renal replacement therapy; an absolute
    # clearance replaces it on hemodialysis or CRRT (Table 2 rows CL, CL_CRRT,
    # CL_HD). The two RRT indicators are mutually exclusive.
    cl_renal <- exp(lcl) * (CRCL / 72)^e_crcl_cl
    no_rrt <- (1 - RRT_HEMODIAL_STATUS) * (1 - RRT_CRRT_STATUS)
    cl_typ <- cl_renal * no_rrt +
      exp(lcl_dial) * RRT_HEMODIAL_STATUS +
      exp(lcl_crrt) * RRT_CRRT_STATUS
    cl <- cl_typ * exp(etalcl)

    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp) * (WT / 60)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
