Sima_2022_ciprofloxacin <- function() {
  description <- "One-compartment parent + one-compartment metabolite population PK model for intravenous ciprofloxacin and its active metabolite desethylene ciprofloxacin in 29 critically ill adults (Sima 2022). Both analytes share one volume of distribution Vd and are modelled in MOLAR units (dose in umol, concentrations in umol/L = nmol/mL). Ciprofloxacin is lost by first-order elimination K and by first-order conversion Kpm to desethylene ciprofloxacin, which is eliminated by first-order Km. All covariate effects are exponential in the RAW, UNCENTERED covariate (Monolix log-linear form): age lowers Vd (-0.022 per year) and Km (-0.035 per year), measured creatinine clearance raises Kpm (+0.81 per mL/s), and carriage of at least one CYP1A2 rs762551 variant allele raises Km (+0.6 on the log scale). Proportional residual error on both analytes."
  reference <- "Sima M, Bobek D, Cihlarova P, Rysanek P, Rousarova J, Berousek J, Kuchar M, Vymazal T, Slanar O. Factors Affecting the Metabolic Conversion of Ciprofloxacin and Exposure to Its Main Active Metabolites in Critically Ill Patients: Population Pharmacokinetic Analysis of Desethylene Ciprofloxacin. Pharmaceutics. 2022;14(8):1627. doi:10.3390/pharmaceutics14081627"
  vignette <- "Sima_2022_ciprofloxacin"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  compartmentData <- list(
    central = list(analyte = "ciprofloxacin", units = "umol", specimen = "serum", verified = TRUE),
    central_desethylenecip = list(
      analyte = "desethylene ciprofloxacin",
      units = "umol",
      specimen = "serum",
      verified = TRUE
    )
  )

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 2: median 57 years (IQR 49-71). Enters Vd and Km as exp(beta x AGE), UNCENTERED, so",
        "lvc and lkel_desethylenecip are the values extrapolated to age 0; at the cohort median 57",
        "years the age factor is exp(-0.022 x 57) = 0.285 on Vd and exp(-0.035 x 57) = 0.136 on Km."
      ),
      source_name = "age"
    ),
    CRCL = list(
      description = "Measured creatinine clearance from a 24-h urine collection, NOT normalized to body surface area",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Sima 2022 Methods 2.1: CLCR = UCR x V / SCR, with UCR the urine creatinine (umol/L), V the",
        "urinary flow rate (mL/s) over a 24-h urine collection and SCR the enzymatic serum creatinine",
        "(umol/L). The paper reports CLCR in mL/s (the Czech SI convention): Table 2 median 1.29",
        "mL/s (IQR 0.74-1.91), i.e. 77.4 mL/min (IQR 44.4-114.6). This column carries mL/min and",
        "model() divides it by 60 before applying the published coefficient 0.81, which is per mL/s.",
        "A user who supplies mL/s here will understate Kpm by a large factor. The effect is",
        "exponential and UNCENTERED: exp(0.81 x CLCR_mL_s) is 1 at CLCR = 0 and 2.84 at the cohort",
        "median 1.29 mL/s."
      ),
      source_name = "CLCR"
    ),
    SNP_CYP1A2_RS762551 = list(
      description = "CYP1A2 rs762551 carrier of at least one variant allele (1 = wt/v or v/v; 0 = wt/wt)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wt/wt)",
      notes = paste(
        "Table 4 footnote: 'CYP1A2_v is at least one variant allele in the CYP1A2 genotype'. Table 2:",
        "wt/wt 6 (20.7%), wt/v 19 (65.5%), v/v 4 (13.8%), so 23 of 29 (79.3%) patients are carriers.",
        "The paper does not say which nucleotide (A or C) it calls the variant allele at rs762551",
        "(c.-163C>A, CYP1A2*1F), so this column is NOT interchangeable with",
        "SNP_CYP1A2_RS762551_C_CARRIER: a C-carrier indicator groups AC with CC, whereas an",
        "A-carrier indicator groups AC with AA. Genotyped by TaqMan allele-specific RT-PCR."
      ),
      source_name = "CYP1A2_v"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight. Tested as a continuous covariate; not retained.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 2: median 90 kg (IQR 70-100). Listed among the continuous covariates tested in Methods 2.5 step (2); no point estimate published."
    ),
    HT = list(
      description = "Height. Tested as a continuous covariate; not retained.",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 2: median 175 cm (IQR 168-182). Positively associated with the desethylene ciprofloxacin / ciprofloxacin AUC12 ratio in univariate regression (r2 = 0.2277, p = 0.0089) but not retained in the popPK model; the Discussion attributes the height association to its correlation with age."
    ),
    CREAT = list(
      description = "Serum creatinine. Tested as a continuous covariate; not retained.",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 2: median 65 umol/L (IQR 53-103). No point estimate published."
    ),
    BILI = list(
      description = "Serum total bilirubin. Tested as a continuous covariate; not retained.",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 2: median 11.3 umol/L (IQR 7.6-17.6). No point estimate published."
    ),
    SEXF = list(
      description = "Female sex indicator. Tested as a categorical covariate; not retained.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Results 3: 20 males and 9 females. No point estimate published."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 29L,
    n_studies = 1L,
    n_observations = 174L,
    age_range = "IQR 49-71 years",
    age_median = "57 years",
    weight_range = "IQR 70-100 kg",
    weight_median = "90 kg",
    sex_female_pct = 31.0,
    race_ethnicity = NULL,
    disease_state = "Critically ill adult ICU patients receiving intravenous ciprofloxacin as part of routine care",
    dose_range = "400 or 600 mg every 12 h by intravenous infusion (physician's choice)",
    regions = "Czechia (single centre: Department of Anesthesiology and ICM, Motol University Hospital, Prague)",
    renal_function = "Measured (24-h urine) creatinine clearance median 1.29 mL/s (IQR 0.74-1.91), i.e. 77 mL/min (IQR 44-115); serum creatinine median 65 umol/L (IQR 53-103).",
    pharmacogenomics = "CYP1A2 rs762551 wt/wt 6 (20.7%), wt/v 19 (65.5%), v/v 4 (13.8%); ABCB1 rs2032582, ABCB1 rs1045642 and SLCO1A2 rs4148977 also genotyped but not retained (Table 2).",
    notes = paste(
      "Prospective low-intervention PK study, February 2019 to June 2020 (EudraCT 2019-003732-24).",
      "Three serum samples per patient at 1, 4 and 11.5 h after the end of an infusion; 87",
      "ciprofloxacin and 87 desethylene ciprofloxacin concentrations (Results 3). UHPLC-MS/MS with",
      "LOQ 30 ng/mL (ciprofloxacin) and 3.0 ng/mL (desethylene ciprofloxacin) (Table 1). Mass",
      "concentrations were converted to molar concentrations (nmol/mL) for analysis (Methods 2.4).",
      "Estimation by SAEM in Monolix 2021R1. Formyl ciprofloxacin and oxociprofloxacin were",
      "measured but only desethylene ciprofloxacin was modelled.",
      collapse = " "
    )
  )

  ini({
    # Structural parameters -- Sima 2022 Table 4 'Fixed Effects'. Results 3
    # equation block (Monolix log-linear covariate model):
    #   log(Vd)  = log(Vd_pop)  + beta_Vd_age * age + eta_Vd
    #   log(K)   = log(K_pop)   + eta_K
    #   log(Km)  = log(Km_pop)  + beta_Km_age * age + beta_Km_CYP1A2_v * (CYP1A2 = v) + eta_Km
    #   log(Kpm) = log(Kpm_pop) + beta_Kpm_CLCR * CLCR + eta_Kpm
    # The covariates are NOT centred, so Vd_pop, Km_pop and Kpm_pop are the
    # values extrapolated to age 0 / CLCR 0 and are not typical-patient values.
    lvc <- log(565.62); label("Volume of distribution of ciprofloxacin and desethylene ciprofloxacin at age 0 (L)") # Table 4, Vd_pop = 565.62 L (R.S.E. 25.7%)
    lkel <- log(0.07); label("Ciprofloxacin elimination rate constant K (1/h)") # Table 4, K_pop = 0.07 1/h (R.S.E. 21.3%)
    lkel_desethylenecip <- log(3.81); label("Desethylene ciprofloxacin elimination rate constant Km at age 0, CYP1A2 wt/wt (1/h)") # Table 4, Km_pop = 3.81 1/h (R.S.E. 34.6%)
    lkmet_desethylenecip <- log(0.017); label("Ciprofloxacin-to-desethylene ciprofloxacin transfer rate constant Kpm at CLCR 0 (1/h)") # Table 4, Kpm_pop = 0.017 1/h (R.S.E. 30.1%)

    # Covariate effects -- exponential in the raw covariate.
    e_age_vc <- -0.022; label("Exponential age effect on Vd (per year)") # Table 4, beta_Vd_age = -0.022 (R.S.E. 20.5%)
    e_age_kel_desethylenecip <- -0.035; label("Exponential age effect on Km (per year)") # Table 4, beta_Km_age = -0.035 (R.S.E. 18.3%)
    e_cyp1a2_kel_desethylenecip <- 0.6; label("Log-scale shift in Km for CYP1A2 rs762551 variant-allele carriers (unitless)") # Table 4, beta_Km_CYP1A2_v = 0.6 (R.S.E. 34.5%)
    e_crcl_kmet_desethylenecip <- 0.81; label("Exponential creatinine-clearance effect on Kpm (per mL/s)") # Table 4, beta_Kpm_CLCR = 0.81 (R.S.E. 19.2%)

    # IIV. Table 4 reports these under 'Standard deviation of the random
    # effects', so each published value is an SD on the log scale and is
    # squared to the variance nlmixr2 expects. No correlations were reported.
    etalvc ~ 0.0841 # Table 4, omega_Vd = 0.29 SD -> 0.29^2 = 0.0841 (R.S.E. 16.7%)
    etalkel ~ 0.0576 # Table 4, omega_K = 0.24 SD -> 0.24^2 = 0.0576 (R.S.E. 35.7%)
    etalkel_desethylenecip ~ 0.0961 # Table 4, omega_Km = 0.31 SD -> 0.31^2 = 0.0961 (R.S.E. 30.3%)
    etalkmet_desethylenecip ~ 0.3481 # Table 4, omega_Kpm = 0.59 SD -> 0.59^2 = 0.3481 (R.S.E. 16.6%)

    # Residual error. Results 3: 'A proportional error model was the most
    # accurate for the residual ... variability'; Monolix's proportional
    # error parameter b is the SD of the proportional error.
    propSd <- 0.22; label("Proportional residual error, ciprofloxacin (fraction)") # Table 4, b1_parent drug = 0.22 (R.S.E. 11.0%)
    propSd_desethylenecip <- 0.25; label("Proportional residual error, desethylene ciprofloxacin (fraction)") # Table 4, b2_metabolite = 0.25 (R.S.E. 11.5%)
  })

  model({
    # CRCL is supplied in mL/min; the published coefficient is per mL/s.
    crclSi <- CRCL / 60

    vc <- exp(lvc + e_age_vc * AGE + etalvc)
    kel <- exp(lkel + etalkel)
    kel_desethylenecip <- exp(
      lkel_desethylenecip +
        e_age_kel_desethylenecip * AGE +
        e_cyp1a2_kel_desethylenecip * SNP_CYP1A2_RS762551 +
        etalkel_desethylenecip
    )
    kmet_desethylenecip <- exp(lkmet_desethylenecip + e_crcl_kmet_desethylenecip * crclSi + etalkmet_desethylenecip)

    # Parent-metabolite structure (Results 3): one compartment for each
    # analyte, first-order elimination of both, and unidirectional transfer
    # from parent to metabolite. The transfer drains the parent, as in
    # Monolix's transfer() macro, so total ciprofloxacin loss is
    # (K + Kpm) * central. Both analytes share the single volume Vd, and the
    # model is molar, so no molecular-weight correction is applied to the
    # transferred amount. Dose (umol of ciprofloxacin) goes into `central` as
    # an IV infusion.
    d/dt(central) <- -kel * central - kmet_desethylenecip * central
    d/dt(central_desethylenecip) <- kmet_desethylenecip * central - kel_desethylenecip * central_desethylenecip

    Cc <- central / vc
    Cc_desethylenecip <- central_desethylenecip / vc

    Cc ~ prop(propSd)
    Cc_desethylenecip ~ prop(propSd_desethylenecip)
  })
}
