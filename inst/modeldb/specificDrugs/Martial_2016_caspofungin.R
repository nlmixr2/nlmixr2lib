Martial_2016_caspofungin <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous caspofungin in adult",
    "intensive-care-unit patients (all Child-Pugh B through hypoalbuminaemia) treated",
    "for suspected or proven invasive fungal infection or prophylaxis (Martial 2016).",
    "Body weight enters a priori by allometry on CL (exponent 0.75) and V1 (exponent 1),",
    "standardised to 70 kg; Q and V2 carry no weight scaling. Log-normal IIV on CL, V1",
    "and V2 (encoded as diagonal; the reported CL-V1 correlation has no published",
    "value) and a combined proportional + fixed additive residual error. Stepwise",
    "covariate screening, including the Child-Pugh score, retained no further covariates.",
    sep = " "
  )
  reference <- paste(
    "Martial LC, Bruggemann RJM, Schouten JA, van Leeuwen HJ, van Zanten AR,",
    "de Lange DW, Muilwijk EW, Verweij PE, Burger DM, Aarnoutse RE, Pickkers P,",
    "Dorlo TPC. Dose Reduction of Caspofungin in Intensive Care Unit Patients with",
    "Child Pugh B Will Result in Suboptimal Exposure. Clin Pharmacokinet.",
    "2016;55:723-733. doi:10.1007/s40262-015-0347-2. PMCID: PMC4875935.",
    sep = " "
  )
  vignette <- "Martial_2016_caspofungin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "caspofungin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "caspofungin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Included a priori (Methods, 'Pharmacokinetic Model'): allometric exponent 0.75",
        "on CL and 1 on V1, standardised to a typical 70 kg patient (Table 2 footnote b;",
        "Discussion). Baseline value only; ICU cohort median 75 kg (range 50-99 kg, Table 1).",
        sep = " "
      ),
      source_name = "BW"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 21L,
    n_studies = 1L,
    n_observations = 419L,
    age_range = "45-80 years (median 71)",
    age_median = "71 years",
    weight_range = "50-99 kg (median 75)",
    weight_median = "75 kg",
    bmi_range = "19.0-36.4 kg/m^2 (median 24.9)",
    sex_female_pct = 38,
    race_ethnicity = "Not reported (Dutch ICU cohort).",
    disease_state = paste(
      "Adults admitted to the intensive care unit receiving caspofungin for suspected or",
      "proven fungal infection (16 normally sterile site, 4 blood) or prophylaxis (1).",
      "All 21 patients were classified Child-Pugh B and all had hypoalbuminaemia",
      "(<= 34 g/L); none were neutropenic (Table 1).",
      sep = " "
    ),
    hepatic_function = "Child-Pugh B in 100% of patients, driven by hypoalbuminaemia (non-cirrhotic).",
    dose_range = paste(
      "70 mg IV loading dose on day 1, then 50 mg daily (BW <= 80 kg) or 70 mg daily",
      "(BW > 80 kg), each infused over about 1 h; one patient received 35 mg daily for",
      "3 days. Sampling for up to 14 days of treatment.",
      sep = " "
    ),
    regions = "Netherlands (NCT01533558).",
    notes = paste(
      "Baseline demographics from Martial 2016 Table 1. Intensive sampling on days 3 and",
      "7 (pre-dose to 24 h), plus trough samples on other days; 419 concentrations.",
      sep = " "
    )
  )

  ini({
    # Structural parameters: Martial 2016 Table 2 (final NONMEM FOCE-I estimates).
    lcl <- log(0.55); label("Clearance for a 70 kg patient (L/h)") # Table 2: CL = 0.55 L/h [RSE 7%]
    lvc <- log(8.93); label("Central volume V1 for a 70 kg patient (L)") # Table 2: V1 = 8.93 L [RSE 8%]
    lq <- log(0.707); label("Intercompartmental clearance Q (L/h)") # Table 2: Q = 0.707 L/h [RSE 3%]
    lvp <- log(4.98); label("Peripheral volume V2 (L)") # Table 2: V2 = 4.98 L [RSE 17%]

    # A priori allometric exponents (Methods, 'Pharmacokinetic Model'; Discussion).
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL (unitless)") # Methods: power exponent 0.75 on CL
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V1 (unitless)") # Methods: power exponent 1 on V1

    # Log-normal IIV reported as %CV; omega^2 = log(CV^2 + 1). IIV on Q held
    # at 0 in the source (omitted here). The CL-V1 correlation (dOFV 8.6,
    # Results) has no reported value, so the block is diagonal.
    etalcl ~ log(1 + 0.307^2) # Table 2: IIV CL 30.7% CV (shrinkage 0%)
    etalvc ~ log(1 + 0.256^2) # Table 2: IIV V1 25.6% CV (shrinkage 6%)
    etalvp ~ log(1 + 0.758^2) # Table 2: IIV V2 75.8% CV (shrinkage 16%)

    # Combined residual error (Results, 'Population Pharmacokinetic Model').
    propSd <- 0.148; label("Proportional residual error (fraction)") # Table 2: proportional residual error 0.148 [RSE 17%]
    addSd <- fixed(0.01); label("Additive residual error (mg/L)") # Table 2: additional residual error 0.01 mg/L, not estimated
  })

  model({
    # Weight enters CL and V1 only; Table 2 footnote c places Q and V2 at the
    # population typical value without weight standardisation.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Caspofungin is infused over about 1 h into the central compartment; the
    # infusion duration is set on the dose record.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L = mg/L (equivalent to ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
