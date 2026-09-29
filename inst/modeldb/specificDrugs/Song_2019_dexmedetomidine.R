Song_2019_dexmedetomidine <- function() {
  description <- "Two-compartment population PK model for intravenous dexmedetomidine in mechanically ventilated children aged 2-12 years in the ICU after neurosurgery, with fixed allometric body-weight scaling to a 70 kg reference"
  reference <- paste(
    "Song IK, Yi S, Lim HS, Lee JH, Kim EH, Cho JY, Kim MC, Kim JT, Kim HS.",
    "A Population Pharmacokinetic Model of Intravenous Dexmedetomidine for",
    "Mechanically Ventilated Children after Neurosurgery.",
    "J Clin Med. 2019;8(10):1563. doi:10.3390/jcm8101563",
    sep = " "
  )
  vignette <- "Song_2019_dexmedetomidine"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline body weight. The only covariate retained in the final model.",
        "Enters as a fixed allometric power model standardised to 70 kg",
        "(Song 2019 Methods 2.5 and Table 2 footnote c):",
        "CL = CLpop * (WT/70)^0.75, V1 = V1pop * (WT/70),",
        "Q = Qpop * (WT/70)^0.75, V2 = V2pop * (WT/70).",
        "Study cohort median weight 22.0 kg (low-dose) and 23.0 kg (high-dose);",
        "no child weighed anywhere near the 70 kg reference, so the reference",
        "is an extrapolation anchor rather than an observed size.",
        sep = " "
      ),
      source_name = "BW"
    )
  )

  # Screened by forward selection (Song 2019 Methods 2.5) but not retained:
  # "No other covariates had a significant influence on the weight-adjusted PK
  # parameters" (Results 3.3). Recorded here for provenance only; none is
  # referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened in the Song 2019 forward-selection covariate step; not retained."
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      notes = "Screened in the Song 2019 forward-selection covariate step; not retained."
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Screened in the Song 2019 forward-selection covariate step; not retained."
    ),
    LBM = list(
      description = "Lean body mass",
      units = "kg",
      type = "continuous",
      notes = "Screened in the Song 2019 forward-selection covariate step; not retained."
    ),
    IBW = list(
      description = "Ideal body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened in the Song 2019 forward-selection covariate step; not retained."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened in the Song 2019 forward-selection covariate step; not retained."
    ),
    BODYFAT_PCT = list(
      description = "Body fat percentage",
      units = "%",
      type = "continuous",
      notes = "Screened in the Song 2019 forward-selection covariate step; not retained."
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "dexmedetomidine",
      units = "ug",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "dexmedetomidine",
      units = "ug",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 29,
    n_studies = 1,
    age_range = "2-12 years",
    age_median = "8.0 years (low-dose, IQR 5.0-10.0); 7.0 years (high-dose, IQR 3.3-10.3)",
    weight_median = "22.0 kg (low-dose, IQR 19.5-30.5); 23.0 kg (high-dose, IQR 16.5-37.8)",
    sex_female_pct = 51.7,
    disease_state = paste(
      "Mechanically ventilated in the ICU after elective neurosurgery",
      "(mostly craniotomy and tumour removal); ASA physical status 1-2;",
      "no cardiovascular, hepatic or renal disease.",
      sep = " "
    ),
    dose_range = paste(
      "0.25 ug/kg IV loading over 10 min then 0.25 ug/kg/h for 50 min (low-dose, n = 15);",
      "0.5 ug/kg IV loading over 10 min then 0.5 ug/kg/h for 50 min (high-dose, n = 14)",
      sep = " "
    ),
    regions = "Republic of Korea (single centre, Seoul National University Hospital)",
    notes = paste(
      "Baseline demographics in Song 2019 Table 1; 264 plasma samples, 19 below the",
      "0.005 ug/L LLOQ. Sampling before infusion, at 10/30/60 min after the start of",
      "infusion, and at 15/30/60/120/240/480 min after the end of infusion.",
      "Race/ethnicity is not tabulated in Table 1; the Discussion refers to the",
      "cohort as Korean children. Two of 31 enrolled patients were excluded before",
      "allocation. Trial registration KCT0001150.",
      sep = " "
    )
  )

  ini({
    # Structural parameters, standardised to a 70 kg body weight
    # (Song 2019 Table 2, 'Structural model' rows; bootstrap 95% CI in
    # parentheses below is the paper's, not encoded here).
    lcl <- log(81.0)
    label("Clearance standardised to 70 kg (L/h)") # Table 2: CL_pop 81.0, RSE 5.5%, bootstrap 81.1 (72.9-90.9)
    lvc <- log(64.2)
    label("Central volume of distribution standardised to 70 kg (L)") # Table 2: V1_pop 64.2, RSE 12.6%, bootstrap 63.7 (50.6-81.0)
    lq <- log(116.4)
    label("Intercompartmental clearance standardised to 70 kg (L/h)") # Table 2: Q_pop 116.4, RSE 13.1%, bootstrap 119.2 (90.6-156.0)
    lvp <- log(167)
    label("Peripheral volume of distribution standardised to 70 kg (L)") # Table 2: V2_pop 167, RSE 12.5%, bootstrap 167 (132-217)

    # Allometric exponents. Methods 2.5 states they were held at the
    # theory-based values (0.75 for flows, 1 for volumes) citing refs 23-24,
    # and Table 2 reports no uncertainty for them, so both are fixed.
    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent of body weight on CL and Q (unitless)") # Methods 2.5; Table 2 footnote c
    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent of body weight on V1 and V2 (unitless)") # Methods 2.5; Table 2 footnote c

    # Inter-individual variability. Methods 2.5 specifies an exponential
    # (log-normal) model Pi = theta * exp(eta_i) with Var(eta) = omega^2, and
    # Table 2 reports each IIV as a CV%. The paper does not state which CV
    # definition it used, so the exact log-normal relation is applied:
    # omega^2 = log(1 + CV^2). See the vignette 'Assumptions and deviations'
    # section for the alternative reading (CV% = omega * 100) and its size.
    etalcl ~ 0.0708694 # Table 2, row 'CL (CV%)' = 27.1; log(1 + 0.271^2)
    etalvc ~ 0.3074847 # Table 2, row 'V1 (CV%)' = 60.0; log(1 + 0.600^2)
    etalq ~ 0.1972832 # Table 2, row 'Q (CV%)' = 46.7; log(1 + 0.467^2)
    etalvp ~ 0.3136780 # Table 2, row 'V2 (CV%)' = 60.7; log(1 + 0.607^2)

    # Combined additive + proportional residual error (Methods 2.5 lists this
    # among the tested forms; Table 2 reports both components for the final model).
    addSd <- 0.0227
    label("Additive residual error (ug/L)") # Table 2: additive error 0.0227 ug/L, RSE 81.5%
    propSd <- 0.427
    label("Proportional residual error (fraction)") # Table 2: proportional error 42.7%, RSE 14.2%
  })

  model({
    # Individual parameters: fixed allometric power model on body weight,
    # standardised to 70 kg (Song 2019 Table 2 footnote c).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q <- exp(lq + etalq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment disposition with first-order elimination, IV input into
    # the central compartment (Results 3.3).
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
