Smit_2019_tobramycin <- function() {
  description <- "Two-compartment IV population PK model for tobramycin in morbidly obese adults undergoing bariatric surgery and non-obese healthy volunteers with normal renal function, total body weight 57-194 kg (Smit 2019). Clearance is a linear function of the MDRD eGFR de-indexed for body surface area (mL/min, centred at 115 mL/min); central volume scales linearly with total body weight (reference 70 kg). Combined additive and proportional residual error."
  reference <- "Smit C, Wasmann RE, Wiezer MJ, van Dongen HPA, Mouton JW, Bruggemann RJM, Knibbe CAJ. Tobramycin Clearance Is Best Described by Renal Function Estimates in Obese and Non-obese Individuals: Results of a Prospective Rich Sampling Pharmacokinetic Study. Pharm Res. 2019;36(8):112. doi:10.1007/s11095-019-2651-2. Correction: Pharm Res. 2020;37(2):27. doi:10.1007/s11095-020-2756-7 (corrected Table 2)"
  vignette <- "Smit_2019_tobramycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Tobramycin is given as a 0.5-h IV infusion (Smit 2019 Methods, Study
  # Procedures), so the dose enters `central` directly. `central` is verified:
  # total tobramycin PLASMA concentrations were measured (Methods, Cobas TOBR2
  # immunoassay). The paper does not state what the peripheral compartment
  # represents.
  compartmentData <- list(
    central = list(analyte = "tobramycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tobramycin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference weight 70 kg (Smit 2019 Table II, Vc = Vc70kg * (TBW/70)). Linear (exponent fixed to 1; the 'TBW on Vc' model in Table III has 8 parameters, i.e. no estimated exponent). Table I medians 137.8 kg (obese, range 103-194) and 66.3 kg (non-obese, range 57-91).",
      source_name = "TBW"
    ),
    CRCL = list(
      description = "MDRD estimated glomerular filtration rate DE-INDEXED for body surface area: MDRD (mL/min/1.73 m^2) * BSA / 1.73, with the 4-variable MDRD = 186.3 * (creatinine umol/L / 88.4)^-1.154 * age^-0.203 * 0.742 (if female) * 1.210 (if black) and the Du Bois BSA = 0.007184 * TBW^0.425 * height(cm)^0.725",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "RAW, NOT BSA-normalized (mL/min, not mL/min/1.73 m^2). Smit 2019 Methods ('de-indexed for body surface area (BSA) by multiplying with individual BSA/1.73') and supplemental equations; the Discussion names the Du Bois BSA formula. Enters CL as a CENTRED linear term CL = CL115 * (1 + Z * (CRCL - 115)) (corrected Table 2). Table I medians 127.5 mL/min (obese, range 77-171) and 102.6 mL/min (non-obese, range 91-120). The authors warn against extrapolating outside the observed 77-171 mL/min range.",
      source_name = "MDRD"
    )
  )

  # Screened in the Smit 2019 covariate analysis (Methods, Covariate Analysis;
  # Results, Table III) but not retained in the final model.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Listed among the screened covariates (Methods, Covariate Analysis); not retained. Table I medians 43.0 years (obese) and 22.5 years (non-obese).",
      source_name = "age"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female)",
      units = "(binary)",
      type = "binary",
      reference_category = "male",
      notes = "Listed among the screened covariates (Methods, Covariate Analysis); not retained.",
      source_name = "sex"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Listed among the screened covariates; not retained.",
      source_name = "BMI"
    ),
    LBW = list(
      description = "Lean body weight (Janmahasatian)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear covariate on Vc gave OFV -23.9, worse than TBW (power) on Vc at -25.9 (Results); not retained.",
      source_name = "LBW"
    ),
    ABW = list(
      description = "Adjusted body weight, IBW + 0.4 * (TBW - IBW)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear covariate on Vc gave OFV -29.2 but TBW was selected for its better goodness of fit at high concentrations (p > 0.05 between them; Results); not retained.",
      source_name = "ABW"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 28L,
    n_studies = 1L,
    n_concentrations = 280L,
    age_range = "20-54 years",
    age_median = "43.0 years (obese), 22.5 years (non-obese)",
    weight_range = "57-194 kg",
    weight_median = "137.8 kg (obese), 66.3 kg (non-obese)",
    bmi_range = "19-53 kg/m^2",
    sex_female_pct = 53.6,
    disease_state = "20 morbidly obese adults (BMI > 40 kg/m^2, or > 35 kg/m^2 with comorbidities) undergoing laparoscopic gastric sleeve or bypass surgery, and 8 non-obese healthy volunteers (BMI 18-25 kg/m^2); all with eGFR >= 60 mL/min and no nephrotoxic co-medication.",
    dose_range = "Single 0.5-h IV infusion: 5 mg/kg lean body weight (Janmahasatian) in obese patients, 5 mg/kg total body weight in non-obese volunteers; median 340 mg (range 240-480) and 320 mg (280-440)",
    regions = "Netherlands (St. Antonius Hospital, Nieuwegein)",
    renal_function = "De-indexed MDRD median 127.5 mL/min (obese, range 77-171) and 102.6 mL/min (non-obese, range 91-120); 24-h urine GFR median 163.3 and 124.7 mL/min",
    notes = "Prospective rich-sampling study (Dutch Trial Registry NTR6058), 10 samples per subject over 24 h; 23 of 280 concentrations (8.2%) below the 0.3 mg/L LLOQ, handled with the M3 method. NONMEM 7.3. Female percentage is 15/28 from Table I (male/female 9/11 obese, 4/4 non-obese)."
  )

  ini({
    # Structural parameters: Smit 2019 Table II final model, as corrected by the
    # 2020 correction (doi:10.1007/s11095-020-2756-7), which fixes the units of CL
    # and Q to L/h (originally printed L/min) and Vp to 4.34 L.
    lcl <- log(6.33); label("Clearance at de-indexed MDRD = 115 mL/min (L/h)") # Corrected Table 2 CL_MDRD115 6.33 L/h (RSE 2.6%)
    lvc <- log(10.6); label("Central volume of distribution at WT = 70 kg (L)") # Corrected Table 2 Vc70kg 10.6 L (RSE 11%)
    lq <- log(6.69); label("Intercompartmental clearance (L/h)") # Corrected Table 2 Q 6.69 L/h (RSE 12%)
    lvp <- log(4.34); label("Peripheral volume of distribution (L)") # Corrected Table 2 Vp 4.34 L (RSE 18%); original Table II printed 4.35

    e_crcl_cl <- 0.0099; label("Linear slope of CL on (CRCL - 115) (per mL/min)") # Corrected Table 2 Z 0.00990 (RSE 10%); Fig. 5 nomogram 0.0099

    # IIV reported as CV% (Table 2), converted with omega^2 = log(1 + CV^2) for
    # the exponential model theta_i = theta * exp(eta_i) of Eq. 1.
    etalcl ~ 0.014297 # Corrected Table 2 IIV CL 12.0% -> log(1 + 0.120^2)
    etalvc ~ 0.060142 # Corrected Table 2 IIV Vc 24.9% -> log(1 + 0.249^2)

    # Residual error, Eq. 2: Y = Cpred + Cpred * eps1 + eps2. Additive error is
    # printed in mg/L, so both entries are read as standard deviations.
    propSd <- 0.116; label("Proportional residual error (fraction)") # Corrected Table 2 proportional error 0.116 (RSE 12%)
    addSd <- 0.346; label("Additive residual error (mg/L)") # Corrected Table 2 additive error 0.346 mg/L (RSE 12%)
  })
  model({
    # Individual PK parameters (corrected Table 2 final-model equations)
    cl <- exp(lcl + etalcl) * (1 + e_crcl_cl * (CRCL - 115))
    vc <- exp(lvc + etalvc) * (WT / 70)
    q <- exp(lq)
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
