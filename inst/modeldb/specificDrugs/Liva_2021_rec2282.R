Liva_2021_rec2282 <- function() {
  description <- "Two-compartment population PK model for the oral pan-histone-deacetylase inhibitor REC-2282 (AR-42) in 55 adults with relapsed or refractory solid tumours or haematologic malignancies pooled from the first-in-human trial OSU09102 (NCT01129193) and the phase 1 acute myeloid leukaemia trial OSU11130 (NCT01798901) (Liva 2021). Oral doses of 20-80 mg enter a depot after an absorption lag time and pass through one transit compartment to the central compartment, with the same first-order rate constant ka for both steps; elimination is first-order from the central compartment and all volumes and clearances are apparent (CL/F, Vc/F, Q/F, Vp/F). Fat-free mass (Janmahasatian) is a power covariate on CL/F centred on the cohort median of 50.8 kg, and tumour type (haematologic vs solid) and formulation (tablet vs capsule) are linear fractional covariates on the lag time. IIV on CL/F and ka; residual error is additive on log-transformed concentrations."
  reference <- paste(
    "Liva S, Chen M, Mortazavi A, Walker A, Wang J, Dittmar K, Hofmeister C, Coss CC, Phelps MA.",
    "Population Pharmacokinetic Analysis from First-in-Human Data for HDAC Inhibitor, REC-2282 (AR-42),",
    "in Patients with Solid Tumors and Hematologic Malignancies: A Case Study for Evaluating Flat vs.",
    "Body Size Normalized Dosing.",
    "Eur J Drug Metab Pharmacokinet. 2021;46:807-816.",
    "doi:10.1007/s13318-021-00722-z.",
    sep = " "
  )
  vignette <- "Liva_2021_rec2282"

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
  )

  compartmentData <- list(
    depot = list(analyte = "REC-2282", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "REC-2282", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "REC-2282", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "REC-2282", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass by the Janmahasatian equation",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed at baseline. Liva 2021 Section 2.4 computes FFM from total body weight and",
        "BMI as FFM = 9.27e3 * TBW / (6.68e3 + 216 * BMI) in males and",
        "9.27e3 * TBW / (8.78e3 + 244 * BMI) in females (Janmahasatian et al. 2005). Continuous",
        "covariates were normalised by the population median (Section 2.3, Equation 1), and",
        "Table 2 gives the median FFM of all patients as 50.8 kg (range 29.9-76.0). FFM was the",
        "only covariate retained on CL/F; it beat weight, BSA, height, LBW and sex in the",
        "univariate screen (Table 4, dOFV -9.02) but explains only about 2.6% of the IIV on CL/F",
        "(Discussion).",
        sep = " "
      ),
      source_name = "FFM (Liva 2021 Table 3 'Covariates FFM'; Table 2; Section 2.4)"
    ),
    TUMTP_HEME = list(
      description = "Haematologic-malignancy tumour-type indicator; 1 = haematologic malignancy, 0 = solid tumour",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (solid tumour)",
      notes = paste(
        "Time-fixed. Liva 2021 Section 2.3 evaluates tumour type 'as dichotomous variables for",
        "solid vs. heme tumor'. Supplemental Table 2 lists the malignancies: haematologic =",
        "relapsed multiple myeloma (17), relapsed lymphoma (10) and relapsed/refractory AML (13);",
        "solid = transitional cell carcinoma (4), breast (2), neurofibromatosis type 2 (5),",
        "meningioma (2), Sertoli cell (1), lung (1), melanoma (1) and carcinoma of unknown",
        "primary (1). All OSU11130 patients have AML, so the authors note that tumour type 'may",
        "be attributed to a study effect between these two trials' (Discussion). The paper does",
        "not state which level is coded 1; see the vignette Assumptions and deviations section.",
        "Applied as ALAG = TV * (1 + 0.942 * TUMTP_HEME) (Equation 2).",
        sep = " "
      ),
      source_name = "TMR (Liva 2021 Table 3 'Covariates TMR')"
    ),
    FORM_CAPSULE = list(
      description = "Capsule formulation indicator; 1 = capsule, 0 = tablet",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (capsule) carries the typical lag time",
      notes = paste(
        "Time-fixed per dosing occasion. Liva 2021 Table 3 footnote defines FORM as",
        "'formulation (capsule vs. tablet)'. The paper does not state which level is coded 1",
        "nor which patients received which formulation; see the vignette Assumptions and",
        "deviations section. Applied as ALAG = TV * (1 + 0.74 * (1 - FORM_CAPSULE)), i.e. the",
        "tablet carries the 0.74 fractional increase in the lag time.",
        sep = " "
      ),
      source_name = "FORM (Liva 2021 Table 3 'Covariates FORM')"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Significant on CL/F in the univariate screen (Table 4, dOFV -5.76) but not retained in the final model; FFM was kept instead. Median 76.1 kg (42.9-122.4), Table 2.",
      source_name = "Weight (Liva 2021 Table 4)"
    ),
    BSA = list(
      description = "Body surface area (Mosteller)",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Significant on CL/F in the univariate screen (Table 4, dOFV -7.01) but not retained. Median 1.9 m^2 (1.4-2.4), Table 2.",
      source_name = "BSA (Liva 2021 Table 4)"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Significant on CL/F in the univariate screen (Table 4, dOFV -7.87) but not retained.",
      source_name = "Height (Liva 2021 Table 4)"
    ),
    LBM = list(
      description = "Lean body weight (James equation from total body weight and height)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "The paper labels this lean body weight (LBW), computed by the James equation from total body weight and height (Section 2.4). Significant on CL/F in the univariate screen (Table 4, dOFV -8.46) but not retained. Median 53.1 kg (34.2-79.8), Table 2.",
      source_name = "LBW (Liva 2021 Table 4)"
    ),
    SEXF = list(
      description = "Biological sex indicator; 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Significant on CL/F in the univariate screen (Table 4, dOFV -8.72) but not retained; 51.8% female (Table 2).",
      source_name = "SEX (Liva 2021 Table 4)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 55L,
    n_studies = 2L,
    age_range = "20-80 years (mean (SD) 60.0 (14.1))",
    age_median = "63 years",
    weight_range = "42.9-122.4 kg (mean (SD) 76.5 (17.8))",
    weight_median = "76.1 kg",
    ffm_range = "29.9-76.0 kg (median 50.8, mean (SD) 51.4 (12.3))",
    bmi_range = "18.5-43.6 kg/m^2 (median 26.0)",
    sex_female_pct = 51.8,
    disease_state = "Relapsed or refractory malignancies: multiple myeloma (17), lymphoma (10), acute myeloid leukaemia (13), and solid tumours including neurofibromatosis type 2, meningioma, transitional cell carcinoma and breast cancer (16)",
    dose_range = "20-80 mg orally three times weekly (Mon/Wed/Fri; one OSU11130 cohort 40 mg four times weekly) for 3 weeks of each 28-day cycle",
    regions = "USA (The Ohio State University Comprehensive Cancer Center)",
    notes = paste(
      "57 patients were enrolled (44 in OSU09102, 13 in OSU11130; Table 1) and 882 plasma",
      "concentrations from 55 of them were modelled (Results 3.2). Table 2 summarises the 56",
      "patients with anthropometric data. Patients were fasted before dosing. Sampling at",
      "pre-dose, (0.25), 0.5, 1, 1.5, 2, 4, 8, 10 and 24 h on day 1, day 5 (OSU11130) and/or",
      "day 19 (OSU09102) (Section 2.1; Supplemental Table 1).",
      sep = " "
    )
  )

  ini({
    # Final covariate model estimates, Liva 2021 Table 3 'Covariate model' columns.
    lcl <- log(11.6)
    label("Apparent clearance CL/F at the median FFM of 50.8 kg (L/h)") # Table 3 covariate model: CL = 11.6 L/h (RSE 9.0%)
    lvc <- log(105)
    label("Apparent central volume of distribution Vc/F (L)") # Table 3 covariate model: Vc = 105 L (RSE 8.0%)
    lq <- log(7.1)
    label("Apparent inter-compartmental clearance Q/F (L/h)") # Table 3 covariate model: Q = 7.1 L/h (RSE 12.5%)
    lvp <- log(76.5)
    label("Apparent peripheral volume of distribution Vp/F (L)") # Table 3 covariate model: Vp = 76.5 L (RSE 12.3%)
    lka <- log(1.28)
    label("Absorption rate constant ka, depot to transit1 and transit1 to central (1/h)") # Table 3 covariate model: ka = 1.28 1/h (RSE 10.6%); Fig. 2 shows the same ka on both steps
    ltlag <- log(0.0453)
    label("Absorption lag time ALAG for a solid-tumour patient on capsules (h)") # Table 3 covariate model: ALAG = 0.0453 h (RSE 36.2%); no IIV in the final model (Results 3.2)

    # Covariate effects. Equation 1 (continuous, power on the median-normalised
    # covariate) and Equation 2 (categorical, TV * (1 + theta * COV)).
    e_ffm_cl <- 0.493
    label("Power exponent on (FFM / 50.8 kg) for CL/F (unitless)") # Table 3 covariate model: FFM = 0.493 (RSE 33.9%)
    e_heme_tlag <- 0.942
    label("Fractional increase in ALAG for haematologic malignancy (unitless)") # Table 3 covariate model: TMR = 0.942 (RSE 60.6%)
    e_tablet_tlag <- 0.74
    label("Fractional increase in ALAG for the tablet formulation (unitless)") # Table 3 covariate model: FORM = 0.74 (RSE 50.0%)

    # IIV. Exponential model (Section 2.3); Table 3 reports IIV as CV%, converted
    # with omega^2 = log(CV^2 + 1).
    etalcl ~ 0.05068 # Table 3 covariate model: IIV CL = 22.8 CV% (RSE 20.6%, shrinkage 15.8%); log(0.228^2 + 1)
    etalka ~ 0.20500 # Table 3 covariate model: IIV ka = 47.7 CV% (RSE 11.5%, shrinkage 5.3%); log(0.477^2 + 1)

    # Residual error. Section 2.3: 'additive error model for log-transformed
    # data'; the Table 3 footnote states the epsilon estimate is an SD.
    expSd <- 0.545
    label("Additive residual SD on log-transformed REC-2282 concentrations") # Table 3 covariate model: epsilon (proportional) = 0.545 (RSE 4.5%)
  })

  model({
    # Individual parameters (Equations 1 and 2)
    cl <- exp(lcl + etalcl) * (FFM / 50.8)^e_ffm_cl
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag) * (1 + e_heme_tlag * TUMTP_HEME) * (1 + e_tablet_tlag * (1 - FORM_CAPSULE))

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Fig. 2: oral dose -> (ALAG) -> depot -(ka)-> transit -(ka)-> Vc <-(Q)-> Vp
    d/dt(depot) <- -ka * depot
    d/dt(transit1) <- ka * depot - ka * transit1
    d/dt(central) <- ka * transit1 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # central in mg, vc in L -> mg/L; x 1000 gives ng/mL (= ug/L, the Fig. 1 axis unit)
    Cc <- central / vc * 1000

    Cc ~ lnorm(expSd)
  })
}
