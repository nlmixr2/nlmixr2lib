Smit_2019_gentamicin <- function() {
  description <- "Two-compartment IV population PK model for gentamicin in morbidly obese adults undergoing bariatric surgery and non-obese healthy volunteers, all with normal renal function (Smit 2019). Clearance and central volume scale as power functions of total body weight (reference 70 kg; estimated exponents 0.729 and 1.25); peripheral volume and intercompartmental clearance carry no covariate. Correlated log-normal IIV on clearance and central volume; combined additive and proportional residual error."
  reference <- "Smit C, Wasmann RE, Goulooze SC, Hazebroek EJ, Van Dongen EPA, Burgers DMT, Mouton JW, Bruggemann RJM, Knibbe CAJ. A Prospective Clinical Study Characterizing the Influence of Morbid Obesity on the Pharmacokinetics of Gentamicin: Towards Individualized Dosing in Obese Patients. Clin Pharmacokinet. 2019;58(10):1333. doi:10.1007/s40262-019-00762-4"
  vignette <- "Smit_2019_gentamicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Gentamicin was infused intravenously over 30 min (Methods, Study
  # Procedures), so the dose enters `central` directly. Total gentamicin was
  # measured in plasma (Methods, Drug Assay).
  compartmentData <- list(
    central = list(analyte = "gentamicin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "gentamicin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference weight 70 kg (Table 2 Vc = Vc70kg x (TBW/70)^X and CL = CL70kg x (TBW/70)^Z; Eqs. 4-5). Table 1: morbidly obese median 148.8 kg (range 109-221), non-obese median 72.9 kg (range 53-86). Enters CL and Vc as power functions with estimated exponents.",
      source_name = "TBW"
    )
  )

  # Screened in the Smit 2019 covariate analysis (Methods and Results; ESM
  # Fig. S1) but not retained in the final model; documentation only.
  covariatesDataExcluded <- list(
    LBM = list(
      description = "Lean body weight (Janmahasatian formula)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "LBW on Vc and CL dropped the base-model OFV by 19.1 and 17.3 points, less than TBW (44.4 and 30.2), with poorer goodness of fit (Results); not retained.",
      source_name = "LBW"
    ),
    ABW = list(
      description = "Adjusted body weight, IBW + 0.4 x (TBW - IBW), IBW by the Devine formula",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "ABW on Vc and CL dropped the base-model OFV by 18.8 and 21.2 points, less than TBW (Results); not retained.",
      source_name = "ABW"
    ),
    CRCL = list(
      description = "Measured glomerular filtration rate from a 24-h urine collection on the study day (not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Absolute mL/min, not normalized to 1.73 m^2. GFR seemed to influence clearance in the base model, but the correlation disappeared once TBW was on clearance (Results); not retained. BMI, age, sex, and MDRD / Cockcroft-Gault eGFR were also screened without effect.",
      source_name = "GFR"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 28L,
    n_studies = 1L,
    n_concentrations = 280L,
    age_range = "19-54 years",
    age_median = "40.5 years (morbidly obese), 22.0 years (non-obese)",
    weight_range = "53-221 kg",
    weight_median = "148.8 kg (morbidly obese), 72.9 kg (non-obese)",
    bmi_range = "18-65 kg/m^2",
    sex_female_pct = 50,
    disease_state = "20 morbidly obese adults (BMI > 40 kg/m^2, or > 35 kg/m^2 with co-morbidities) undergoing laparoscopic gastric bypass or sleeve gastrectomy, and 8 non-obese healthy volunteers (BMI 18-25 kg/m^2); all with eGFR > 60 mL/min.",
    dose_range = "Single 30-min IV infusion: 5 mg/kg lean body weight (morbidly obese) or 5 mg/kg total body weight (non-obese); 240-480 mg.",
    regions = "The Netherlands (St. Antonius Hospital, Nieuwegein)",
    renal_function = "Measured 24-h urine GFR median 171.9 mL/min (range 110-230) in the morbidly obese and 123.7 mL/min (91-170) in the non-obese.",
    notes = "Prospective study (Dutch Trial Registry NTR6058; Smit 2019 Table 1). Ten samples per subject from 5 min to 24 h after the end of infusion; 24 of 280 concentrations (8.6%) below the 0.4 mg/L LLOQ were handled with the M3 method. NONMEM 7.3. Eta shrinkage 4% (CL) and 7% (Vc)."
  )

  ini({
    # Structural parameters (Smit 2019 Table 2, final model). The paper
    # reports clearances in L/min; model() converts them to L/h.
    lcl <- log(0.0892); label("Clearance at WT = 70 kg (L/min)") # Table 2 CL70kg 0.0892 L/min (RSE 5.6%)
    lvc <- log(11.9); label("Central volume at WT = 70 kg (L)") # Table 2 Vc70kg 11.9 L (RSE 8.8%)
    lq <- log(0.0848); label("Intercompartmental clearance (L/min)") # Table 2 Q 0.0848 L/min (RSE 8.2%)
    lvp <- log(7.29); label("Peripheral volume (L)") # Table 2 Vp 7.29 L (RSE 5.7%)

    # Covariate exponents on (WT/70)
    e_wt_cl <- 0.729; label("Power exponent on (WT/70) for CL (unitless)") # Table 2 Z 0.729 (RSE 9.6%)
    e_wt_vc <- 1.25; label("Power exponent on (WT/70) for Vc (unitless)") # Table 2 X 1.25 (RSE 10.8%)

    # Between-subject variability. Table 2 prints IIV as % (CL 18.1%, Vc 19.2%);
    # ESM Eq. 1 defines log-normal etas, so omega^2 = log(1 + CV^2). The
    # covariance 0.0316 is printed on the omega scale.
    etalcl + etalvc ~ c(
      log(1 + 0.181^2),
      0.0316, log(1 + 0.192^2)
    ) # Table 2: IIV CL 18.1%; covariance IIV Vc-CL 0.0316; IIV Vc 19.2%

    # Residual variability: ESM Eq. 2, Y = Cpred + Cpred x eps1 + eps2;
    # Table 2 footnote b: estimates reported as standard deviations.
    propSd <- 0.159; label("Proportional residual error (fraction)") # Table 2 proportional error 0.159 (RSE 8.2%)
    addSd <- 0.206; label("Additive residual error (mg/L)") # Table 2 additive error 0.206 mg/L (RSE 8.4%)
  })
  model({
    # Individual PK parameters (Table 2; Eqs. 4-5). Clearances x 60 converts
    # L/min to L/h so that time is in hours.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * 60
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq) * 60
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volume in L, so central / vc is mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
