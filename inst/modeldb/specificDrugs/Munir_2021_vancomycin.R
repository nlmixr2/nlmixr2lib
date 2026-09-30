Munir_2021_vancomycin <- function() {
  description <- "One-compartment IV population PK model for vancomycin in 58 adult surgical patients in Pakistan (Munir 2021). Clearance carries linear, median-centred effects of raw Cockcroft-Gault creatinine clearance (increasing, reference 101.15 mL/min) and total body weight (decreasing as printed in Eq. 3, reference 75 kg); volume of distribution has no retained covariate. Additive residual error."
  reference <- "Munir MM, Rasheed H, Khokhar MI, Khan RR, Saeed HA, Abbas M, Ali M, Bilal R, Nawaz HA, Khan AM, Qamar S, Anjum SM, Usman M. Dose Tailoring of Vancomycin Through Population Pharmacokinetic Modeling Among Surgical Patients in Pakistan. Front Pharmacol. 2021;12:721819. doi:10.3389/fphar.2021.721819"
  vignette <- "Munir_2021_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Cockcroft-Gault creatinine clearance (raw, NOT BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Munir 2021 Methods (Patient selection and sampling): 'creatinine clearance (CRCL) was calculated using the Cockcroft and Gault equation'; no BSA normalization is mentioned, so values are raw mL/min (per inst/references/covariate-columns.md, CRCL accepts raw Cockcroft-Gault when the source does not BSA-normalize; precedents Zhou_2019_vancomycin.R, Delattre_2010_amikacin.R). Table 1 cohort median 101.15 mL/min (range 15.9-177.2); the same 101.15 is the centring value in Eq. 2.",
      source_name = "CRCL"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Munir 2021 Table 1 cohort median 75 kg (range 53-129); the same 75 kg is the centring value in Eq. 3. The printed linear effect decreases CL with weight (factor 1.24 at 53 kg, 0.41 at 129 kg) and reaches zero at 165.9 kg, so the model must not be used above the observed weight range.",
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Munir 2021 Covariate analysis: screened by stepwise covariate modelling but not retained (Results: 'the age of the patient was not a significant covariate for vancomycin CL'), although Figures 1A and 2A show CL falling with age. Table 1 median 54 years (range 25-86)."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "categorical",
      notes = "Munir 2021 Covariate analysis: sex screened by stepwise covariate modelling; not retained. Table 1: 39 male / 19 female."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Munir 2021 Covariate analysis: SeCR screened; not retained -- its information enters the final model through the Cockcroft-Gault CRCL. Table 1 median 0.935 mg/dL (range 0.4-4.7)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 58L,
    n_studies = 1L,
    age_range = "25-86 years",
    age_median = "54 years",
    weight_range = "53-129 kg",
    weight_median = "75 kg",
    sex_female_pct = 32.8,
    race_ethnicity = "Pakistani (single-centre study, Lahore General Hospital)",
    disease_state = "Adults (> 18 years) admitted to a surgical unit after a major surgical procedure and treated with IV vancomycin: grossly contaminated wounds, gas gangrene, and severe peritonitis (Table 1).",
    dose_range = "500-1000 mg IV over a 0.5 h infusion (most patients received 1000 mg); samples after the first dose",
    regions = "Pakistan (Lahore)",
    renal_function = "Cockcroft-Gault CRCL median 101.15 mL/min (range 15.9-177.2); serum creatinine median 0.935 mg/dL (range 0.4-4.7)",
    n_concentrations = 176L,
    notes = "Prospective single-centre sampling, August-December 2018. 176 plasma concentrations from 58 patients (average 3 per patient, range 1-7), median concentration 23.5 mg/L (range 1.9-52.6). HPLC-UV assay, LLOQ 0.25 mg/L. Fit in NONMEM 7.4.4 with FOCE-I, ADVAN1 TRANS2; a two-compartment model was unstable. Table 1 prints 39 male (73.7%) / 19 female (26.3%) and patient-type percentages that sum over 38 patients, not 58; sex_female_pct here is 19/58. Evaluated by goodness-of-fit plots, a 500-replicate VPC and a 1000-sample bootstrap (Table 2)."
  )

  ini({
    # Structural parameters (Munir 2021 Table 2 'Final estimate' column).
    # Footnote b: CL is at CRCL = 101.15 mL/min and weight 75 kg.
    lcl <- log(2.45)
    label("Clearance at CRCL = 101.15 mL/min and WT = 75 kg (L/h)") # Munir 2021 Table 2: CL = 2.45 L/h (RSE 2%); bootstrap 2.46 (2.35-2.59); Eq. 1 constant
    lvc <- log(22.6)
    label("Volume of distribution (L)") # Munir 2021 Table 2: V = 22.6 L (RSE 5%); bootstrap 22.7 (20.3-26.1)

    # Linear median-centred covariate effects on CL (Munir 2021 Eqs. 1-3):
    #   CL = 2.45 * CLCRCL * CLWT * exp(eta1)
    #   CLCRCL = 1 + 0.0046 * (CRCL - 101.15)
    #   CLWT   = 1 - 0.011  * (WT - 75)
    # Table 2 prints the weight coefficient as 0.011 ('proportional change in
    # CL with weight', footnote d); the minus sign lives in Eq. 3 (confirmed in
    # the JATS MathML), so the signed slope is -0.011 per kg.
    e_crcl_cl <- 0.0046
    label("Linear coefficient for centred CRCL on CL (per mL/min)") # Munir 2021 Table 2: CL-CRCL = 0.0046 (RSE 13%); bootstrap 0.0045 (0.0033-0.0063); Eq. 2
    e_wt_cl <- -0.011
    label("Linear coefficient for centred WT on CL (per kg)") # Munir 2021 Table 2: CL-WT = 0.011 (RSE 10%); bootstrap 0.0107 (0.0087-0.013); sign from Eq. 3 '1 - 0.011 x (WT - 75)'

    # Between-subject variability, exponential (Methods). Table 2 footnote e:
    # BSV 'expressed in percentage coefficient of variation', converted to
    # log-normal variance as omega^2 = log(1 + CV^2). The CV reading is
    # confirmed by Table 3's simulated trough SDs (see the vignette).
    etalcl ~ 0.012688 # Munir 2021 Table 2: BSV-CL = 11.3% CV (RSE 38%); log(1 + 0.113^2)
    etalvc ~ 0.050668 # Munir 2021 Table 2: BSV-Vd = 22.8% CV (RSE 51%); log(1 + 0.228^2)

    # Residual error. Table 2 row 'Additive error' = 3.07 and Results 'the
    # additive error was 3.07'; the same Results paragraph also calls the
    # model 'proportional', but 3.07 is only plausible as an additive SD in
    # mg/L (see the vignette).
    addSd <- 3.07
    label("Additive residual error (mg/L)") # Munir 2021 Table 2: Additive error = 3.07 (RSE 10%); bootstrap 3.0 (2.43-3.67)
  })
  model({
    # Individual PK parameters (Munir 2021 Eqs. 1-3). Both covariate factors
    # are linear and centred on the cohort medians.
    crcl_factor <- 1 + e_crcl_cl * (CRCL - 101.15)
    wt_factor <- 1 + e_wt_cl * (WT - 75)
    cl <- exp(lcl + etalcl) * crcl_factor * wt_factor
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, volume in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
