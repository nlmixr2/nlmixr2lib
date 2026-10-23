Downes_2020_vancomycin <- function() {
  description <- "One-compartment IV population PK model for vancomycin in critically ill children aged 2-18 years in the PICU (Downes 2020). Clearance scales allometrically with body weight (exponent 0.75, reference 34 kg), as a power function of the combined creatinine-cystatin C full age spectrum (FAScomb) estimated GFR (exponent 0.903, reference 145 mL/min/1.73 m^2), and is 0.699-fold lower in females; volume scales linearly with body weight. Additive plus proportional residual error."
  reference <- "Downes KJ, Zane NR, Zuppa AF. Effect of Cystatin C on Vancomycin Clearance Estimation in Critically Ill Children Using a Population Pharmacokinetic Modeling Approach. Ther Drug Monit. 2020;42(6):848-855. doi:10.1097/FTD.0000000000000796. Preliminary results presented as Downes K, Zuppa A. Cystatin C Improves Estimation of Vancomycin Clearance in Critically Ill Children Using a Population Pharmacokinetic Modeling Approach. Open Forum Infect Dis. 2018;5(Suppl 1):S431 (IDWeek 2018 abstract 1402)."
  vignette <- "Downes_2020_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Vancomycin was given as an intermittent IV infusion, so the dose enters
  # `central` directly. Methods ('Dosing, PK sampling, and biomarker
  # measurement') names the vancomycin assay (Abbott Architect chemiluminescent
  # microparticle immunoassay) but not the matrix -- only cystatin C is stated
  # to be measured in plasma -- so "plasma" is the repository default, not a
  # paper-sourced claim.
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling on CL (exponent 0.75, fixed) and V (exponent 1, fixed), normalized to 34 kg (Table 4 footnote a). Table 2 cohort median 33.9 kg (IQR 25.9-45.9).",
      source_name = "WT"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate from the combined creatinine-cystatin C full age spectrum equation (FAScomb), BSA-normalized",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 1: FAScomb eGFR = 107.3 / [(Cr/Q_Cr + CysC/Q_CysC) / 2], with serum creatinine Cr in mg/dL, plasma cystatin C CysC in mg/L, Q_CysC = 0.82 mg/L, and Q_Cr the mean or median creatinine of the 1-year age-specific distribution in healthy children (the age-specific Q_Cr table is not reproduced in the paper; it comes from the FAS publication, Pottel H et al.). Cystatin C values were rounded to the nearest 0.1 mg/L before entering the equation. Stored under canonical CRCL, which covers BSA-normalized GFR estimates regardless of filtration marker (the register lists the combined creatinine-cystatin C CKD-EPI variant under the same canonical). Reference value 145 mL/min/1.73 m^2 equals the Table 2 cohort median at PK sampling (IQR 109-175). The baseline value was carried for doses before the PK-sampling interval and the PK-sampling value for that interval, so the covariate is time-varying in the source data.",
      source_name = "eGFR (FAScomb)"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "male",
      notes = "Multiplicative factor 0.699 on CL if female (Table 4 footnote a: '(0.699 if female)'), encoded as e_sexf_cl^SEXF. 6 of 20 subjects (30%) were female (Table 2). The Discussion notes this sex effect on vancomycin CL had not been previously reported and may reflect the small sample.",
      source_name = "female sex"
    )
  )

  # Tested in the covariate search but not retained (Table 3 and its
  # footnote a), so documentation only and not referenced in model().
  covariatesDataExcluded <- list(
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested directly as a covariate on CL (Table 3 step 2, OFV change -8.933); not retained. It enters the final model only through the FAScomb eGFR. Table 2 median 0.4 mg/dL (IQR 0.3-0.5) on the PK sampling day.",
      source_name = "Cr"
    ),
    CYSC = list(
      description = "Plasma cystatin C",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested directly as a covariate on CL (Table 3 step 2, OFV change -18.289); not retained. It enters the final model only through the FAScomb eGFR. Table 2 median 0.5 mg/L (IQR 0.4-0.8) on the PK sampling day.",
      source_name = "CysC"
    ),
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL and V at each covariate step (Table 3 footnote a); did not significantly improve fit.",
      source_name = "age"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL and V at each covariate step (Table 3 footnote a); did not significantly improve fit.",
      source_name = "height"
    ),
    CONMED_INOTROPE = list(
      description = "Vasoactive medication use indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "no vasoactive medication",
      notes = "Tested on CL and V at each covariate step (Table 3 footnote a); did not significantly improve fit.",
      source_name = "vasoactive medications"
    ),
    CONMED_STEROID = list(
      description = "Systemic corticosteroid use indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "no steroid use",
      notes = "Tested on CL at each covariate step (Table 3 footnote a); did not significantly improve fit. 8 of 20 subjects (40%) received steroids (Table 2).",
      source_name = "steroids"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20L,
    n_studies = 1L,
    n_sites = 1L,
    n_concentrations = 83L,
    age_range = "3.9-18.2 years (eligibility 2-18 years)",
    age_median = "12.7 years (IQR 8.8-15.5)",
    weight_median = "33.9 kg (IQR 25.9-45.9)",
    sex_female_pct = 30,
    disease_state = "Critically ill children in the PICU receiving intravenous vancomycin for any indication (sepsis or rule-out sepsis 65%, pneumonia 20%); 55% mechanically ventilated; chronic lung disease 40%, neurologic disorder 40%, immunocompromising condition 25%. Excluded: hypo- or hyperthyroidism, renal replacement therapy, ECMO.",
    dose_range = "Clinician-chosen intermittent IV regimens, typically 10-15 mg/kg every 6-8 h infused over 1-2 h; initial dosing median 54.8 mg/kg/day (IQR 37.4-59.4), every 6 h in 70% and every 8 h in 30% (Table 2).",
    regions = "United States (Children's Hospital of Philadelphia)",
    renal_function = "FAScomb eGFR median 145 mL/min/1.73 m^2 (IQR 109-175) and bedside Schwartz eGFR median 138 (IQR 97-179) at PK sampling; renal function was largely normal or augmented (Discussion).",
    notes = "Prospective observational study, December 2015 to November 2017. 83 vancomycin concentrations (58 research samples from one steady-state dosing interval -- 1 h and 3.5 h after the end of infusion and a pre-dose trough -- plus 25 clinical samples); one below-quantification clinical trough excluded. Assay range 3.0-50.0 mg/L. Fit in NONMEM 7.4 with FOCE-I."
  )

  ini({
    # Structural parameters: Downes 2020 Table 4, for a typical 34-kg male with
    # an FAScomb eGFR of 145 mL/min/1.73 m^2.
    lcl <- log(4.02); label("Clearance at WT=34 kg, CRCL=145 mL/min/1.73 m^2, male (L/h)") # Table 4: CL 4.02 L/hr (95% CI 3.6-4.2), RSE 5.93
    lvc <- log(16.4); label("Volume of distribution at WT=34 kg (L)")                       # Table 4: V 16.4 L (95% CI 13.8-19.0), RSE 8.05

    # Allometric exponents. Results 'Population PK model': the CL exponent was
    # tested as a random-effect parameter, estimated 0.74, 'thus we included
    # 0.75 as a fixed value'. Methods: 'V_d was scaled to 1 in the base model',
    # and the Table 4 footnote equation carries (WT/34) with no exponent.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on (WT/34) for CL (unitless)") # Results; Table 4 footnote a (WT/34)^0.75
    e_wt_vc <- fixed(1);    label("Allometric exponent on (WT/34) for V (unitless)")  # Methods; Table 4 footnote a (WT/34)

    # Covariate effects on CL (Table 4 footnote a:
    # CL = 4.02 * (WT/34)^0.75 * (eGFR/145)^0.903 * (0.699 if female)).
    e_crcl_cl <- 0.903; label("Power exponent on (CRCL/145) for CL (unitless)")       # Table 4: eGFR from FAScomb 0.903 (95% CI 0.48-1.33), RSE 24.1
    e_sexf_cl <- 0.699; label("Multiplicative factor on CL for females (unitless)")   # Table 4: Female sex 0.699 (95% CI 0.48-0.92), RSE 16.0

    # IIV. Table 4 footnote b: IIV 'presented as percent coefficient of
    # variation calculated by the square root of the variance x 100', so
    # omega^2 = (23.4/100)^2. IIV on V is reported as NA (not estimated).
    etalcl ~ 0.054756 # Table 4: IIV clearance 23.4 (CV), RSE 34.9; omega^2 = 0.234^2

    # Residual error. Results: 'additive plus proportional error structure'.
    # Footnote c states the proportional term is reported as a standard
    # deviation; the additive row has no footnote and no units and is read as an
    # SD in mg/L by maintainers' ruling (2026-09-22). The alternative
    # reading, a NONMEM SIGMA variance, would give
    # addSd = sqrt(2.44) = 1.562 mg/L.
    propSd <- 0.179; label("Proportional residual error (fraction)") # Table 4: residual proportional 0.179, RSE 65.1
    addSd <- 2.44;   label("Additive residual error (mg/L)")         # Table 4: residual additive 2.44, RSE 106.0
  })
  model({
    # Downes 2020 Table 4 footnote a.
    cl <- exp(lcl + etalcl) * (WT / 34)^e_wt_cl * (CRCL / 145)^e_crcl_cl * e_sexf_cl^SEXF
    vc <- exp(lvc) * (WT / 34)^e_wt_vc

    kel <- cl / vc
    d/dt(central) <- -kel * central

    # Dose in mg and volume in L, so central/vc is mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
