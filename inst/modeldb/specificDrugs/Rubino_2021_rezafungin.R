Rubino_2021_rezafungin <- function() {
  description <- "Four-compartment population PK model for rezafungin after IV infusion in healthy subjects (three phase 1 studies) and in patients with candidemia and/or invasive candidiasis (phase 2 STRIVE) (Rubino 2021), with serum-albumin power effects on CL and the first peripheral volume, a female-sex shift on CL, a linear body-surface-area slope on the central volume, BSA power effects on both first and second peripheral volumes, and infection-status shifts on the central and second peripheral volumes."
  reference <- "Rubino CM, Flanagan S. Population pharmacokinetics of rezafungin in patients with fungal infections. Antimicrob Agents Chemother. 2021;65(11):e00842-21. doi:10.1128/AAC.00842-21"
  vignette <- "Rubino_2021_rezafungin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "rezafungin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "rezafungin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "rezafungin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral3 = list(analyte = "rezafungin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on CL (exponent -0.844) and on the first peripheral volume Vp1 (exponent -0.829), both with reference 4.2 g/dL (Rubino 2021 Table 2 footnote a). Rubino 2021 reports albumin in US-convention g/dL; the canonical register unit is SI g/L, so model() converts with alb_gdL <- ALB * 0.1 before the power terms. Reference 4.2 g/dL = 42 g/L. The Methods state that albumin was 'normalized for variable interlaboratory reference ranges' before the covariate analysis; the normalization formula is not given, so supply albumin as measured (or on the local laboratory's scale mapped to a standard reference range, if that is available). Median albumin was 4.76-4.87 g/dL in the three phase 1 studies and 2.81 g/dL (range 1.62-4.97) in STRIVE (Table 1). Lower albumin predicts faster CL and a larger Vp1, consistent with a higher unbound fraction of this highly protein-bound drug (Discussion).",
      source_name = "albumin"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear slope on the central volume, Vc = [11.1 + 7.54 x (BSA - 1.83)] x (...), and power effects on Vp1 (exponent 1.14) and Vp2 (exponent 2.47) with reference 1.83 m^2 (Rubino 2021 Table 2 footnote a). 1.83 m^2 is the median BSA in STRIVE and the MAD and QTcF studies (Table 1). The BSA computation formula is not stated in the source. BSA did not affect CL.",
      source_name = "BSA"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Rubino 2021 Table 2 footnote a defines SEXF = 1 for females and 0 for males, identical to the canonical orientation. Linear fractional effect on CL: CL x (1 - 0.133 x SEXF), i.e. 13.3% lower clearance in females.",
      source_name = "SEXF"
    ),
    DIS_INFECT_ACTIVE = list(
      description = "Infection-status indicator (patient enrolled with a systemic fungal infection)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy phase 1 subject)",
      notes = "Rubino 2021 Methods define infection status as 'healthy' (phase 1 subjects) or 'infected' (phase 2 STRIVE patients enrolled with a systemic fungal infection -- candidemia and/or invasive candidiasis). The paper's 'infected' indicator maps directly onto DIS_INFECT_ACTIVE = 1 with no re-expression; it is time-fixed per subject because every STRIVE patient was enrolled with an active infection. Linear fractional effects: Vc x (1 + 0.478 x DIS_INFECT_ACTIVE) and Vp2 x (1 + 1.69 x DIS_INFECT_ACTIVE) (Table 2 footnote a). The typical values 11.1 L and 6.69 L are therefore the healthy-subject reference. The indicator is confounded with albumin: every healthy subject had albumin above the 95th percentile of the STRIVE patients (Results, Fig. 1).",
      source_name = "infected"
    )
  )

  # Covariates screened in the stepwise forward-selection / backward-elimination
  # analysis (Rubino 2021 Materials and Methods, 'Population pharmacokinetic
  # modeling') but NOT retained in the final model. Documented for provenance
  # only; these names are not referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened but not retained."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened but not retained. The earlier phase 1 model (Lakota 2018) scaled every parameter to body weight allometrically; that scaling was removed before the covariate analysis, and BSA was selected instead."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened but not retained."
    ),
    IBW = list(
      description = "Ideal body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened but not retained."
    ),
    CRCL = list(
      description = "Creatinine clearance normalized to body surface area",
      units = "mL/min/1.73m^2",
      type = "continuous",
      notes = "Screened but not retained. STRIVE enrolled patients with CrCL as low as 8.57 mL/min/1.73 m^2 (Table 1); the absence of a renal-function effect on CL is the paper's basis for concluding that no renal dose adjustment is needed."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened but not retained."
    ),
    RACE_BLACK = list(
      description = "Black race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Race was screened but not retained; the paper does not state how race was coded for the screen."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 135L,
    n_studies = 4L,
    n_healthy = 66L,
    n_patients = 69L,
    n_observations = 1518L,
    age_range = "20-88 years",
    age_median = "36-44.5 years (phase 1 studies); 57.5 years (STRIVE)",
    weight_median = "71-78 kg by study",
    bsa_range = "1.22-2.37 m^2",
    bsa_median = "1.81-1.83 m^2 by study",
    albumin_median = "4.76-4.87 g/dL (phase 1 studies); 2.81 g/dL (STRIVE, range 1.62-4.97)",
    sex_female_pct = 44.9,
    race_ethnicity = c(White = 86.0, Black = 11.0, Other = 2.2, `Not reported` = 0.7),
    disease_state = "Pooled cohort of 66 healthy subjects from three phase 1 studies and 69 patients with candidemia and/or invasive candidiasis from the phase 2 STRIVE trial (NCT02734862). STRIVE enrolled patients with severe renal impairment (CrCL down to 8.57 mL/min/1.73 m^2).",
    dose_range = "Phase 1: single 50, 100, 200 or 400 mg as a 60-min IV infusion (SAD); 100 or 200 mg on days 1 and 8, or 400 mg on days 1, 8 and 15, as 60-min IV infusions (MAD); single 600 or 1,400 mg infused at 400 mg/h (QTcF). Phase 2: 400 mg IV once weekly, or 400 mg in week 1 followed by 200 mg weekly.",
    studies = "NCT02516904 (SAD), NCT02551549 (MAD), CD101.IV.1.06 (QTcF); phase 2 STRIVE NCT02734862",
    notes = "Demographics from Rubino 2021 Table 1 (the STRIVE column counts 70 patients for race and sex, one of whom was later excluded as an outlier; the sex and race percentages above pool 136 subjects). Table 1 prints the STRIVE weight range as '34-55' kg against a mean of 74.5 kg (SD 21.3), which is a typesetting error in the source; the weight range is therefore not recorded here. 1,518 plasma concentrations (416 from STRIVE, six per patient) were modelled after excluding 15 outlier observations and one STRIVE patient whose profile was inconsistent with rezafungin PK; only two concentrations were below the lower limit of quantification. NONMEM 7.2, FOCE-I."
  )

  ini({
    # Structural parameters (Rubino 2021 Table 2, 'Population mean' column).
    # Typical values are for a healthy (DIS_INFECT_ACTIVE = 0) male at
    # albumin = 4.2 g/dL and BSA = 1.83 m^2.
    lcl <- log(0.254); label("Clearance at albumin 4.2 g/dL, male (L/h)") # Rubino 2021 Table 2 'CL, liter/h' = 0.254 (3.74 %SEM)
    lvc <- log(11.1); label("Central volume at BSA 1.83 m^2, healthy (L)") # Rubino 2021 Table 2 'Vc, liter' = 11.1 (8.87 %SEM)
    lq <- log(18.2); label("Intercompartmental clearance to peripheral1 (L/h)") # Rubino 2021 Table 2 'CLd1, liter/h' = 18.2 (1.09 %SEM)
    lvp <- log(14.6); label("First peripheral volume at albumin 4.2 g/dL and BSA 1.83 m^2 (L)") # Rubino 2021 Table 2 'Vp1, liter' = 14.6 (4.42 %SEM)
    lq2 <- log(0.541); label("Intercompartmental clearance to peripheral2 (L/h)") # Rubino 2021 Table 2 'CLd2, liter/h' = 0.541 (8.11 %SEM)
    lvp2 <- log(6.69); label("Second peripheral volume at BSA 1.83 m^2, healthy (L)") # Rubino 2021 Table 2 'Vp2, liter' = 6.69 (11.8 %SEM)
    lq3 <- log(0.0743); label("Intercompartmental clearance to peripheral3 (L/h)") # Rubino 2021 Table 2 'CLd3, liter/h' = 0.0743 (7.62 %SEM)
    lvp3 <- log(13.6); label("Third peripheral volume (L)") # Rubino 2021 Table 2 'Vp3, liter' = 13.6 (8.38 %SEM)

    # Covariate effects (Rubino 2021 Table 2 and footnote a equations).
    e_alb_cl <- -0.844; label("Exponent of (albumin / 4.2 g/dL) on CL (unitless)") # Rubino 2021 Table 2 'CL-albumin power' = -0.844 (9.34 %SEM)
    e_sexf_cl <- -0.133; label("Proportional change in CL for females (fraction)") # Rubino 2021 Table 2 'Proportional change in females' = -0.133 (24.3 %SEM)
    e_bsa_vc <- 7.54; label("Linear slope of Vc on (BSA - 1.83 m^2) (L per m^2)") # Rubino 2021 Table 2 'Vc-BSA slope' = 7.54 (23.6 %SEM)
    e_infect_vc <- 0.478; label("Proportional change in Vc for infected patients (fraction)") # Rubino 2021 Table 2 'Proportional change in infected patients' (Vc block) = 0.478 (27.7 %SEM)
    e_alb_vp <- -0.829; label("Exponent of (albumin / 4.2 g/dL) on Vp1 (unitless)") # Rubino 2021 Table 2 'Vp1-albumin power' = -0.829 (10.6 %SEM)
    e_bsa_vp <- 1.14; label("Exponent of (BSA / 1.83 m^2) on Vp1 (unitless)") # Rubino 2021 Table 2 'Vp1-BSA power' = 1.14 (15.6 %SEM)
    e_bsa_vp2 <- 2.47; label("Exponent of (BSA / 1.83 m^2) on Vp2 (unitless)") # Rubino 2021 Table 2 'Vp2-BSA power' = 2.47 (18.4 %SEM)
    e_infect_vp2 <- 1.69; label("Proportional change in Vp2 for infected patients (fraction)") # Rubino 2021 Table 2 'Proportional change in infected patients' (Vp2 block) = 1.69 (19.5 %SEM)

    # IIV on Vp2 is not a separate random effect: the Vp1 eta is reused,
    # multiplied by this estimated scaling term.
    vp2_eta_scale <- 1.71; label("Scaling of the Vp1 eta applied to Vp2 (unitless)") # Rubino 2021 Table 2 'IIV scaling term relative to Vp1 IIV' = 1.71 (18.8 %SEM)

    # Interindividual variability (Rubino 2021 Table 2). The paper prints the
    # variance with the %CV in parentheses, and its %CV is sqrt(omega^2)
    # (sqrt(0.0562) = 0.237, sqrt(0.148) = 0.385). The Vp1 variance is
    # printed as 0.0107, the same number as the CL-Vp1 covariance on the row
    # above; that value contradicts both its own %CV (sqrt(0.0107) = 10.3%,
    # not 21.7%) and the printed CL-Vp1 r^2 (0.0107^2 / (0.0562 x 0.0107)
    # = 0.19, not 0.043). Both of those give omega^2 near 0.047
    # (0.217^2 = 0.0471; 0.0107^2 / (0.0562 x 0.043) = 0.0474), so 0.0471 is
    # used. Vc and Vp1 are uncorrelated (no covariance was estimated).
    etalcl + etalvc + etalvp ~ c(
      0.0562,
      0.0604, 0.148,
      0.0107, 0, 0.0471
    ) # Rubino 2021 Table 2: omega^2 CL 0.0562 (23.7%), Vc 0.148 (38.5%), Vp1 0.0471 from the printed 21.7% CV; covariances CL-Vc 0.0604 and CL-Vp1 0.0107

    # Residual error: the additive component was removed from the full
    # multivariable model (Results), leaving a proportional error.
    propSd <- 0.0891; label("Proportional residual error (fraction)") # Rubino 2021 Table 2 'Proportional error, %CV' = 8.91 (2.57 %SEM)
  })

  model({
    # 1. Derived covariate terms. The canonical ALB column is SI g/L; the
    #    albumin exponents were estimated against g/dL with a 4.2 g/dL
    #    reference, so convert first.
    alb_gdL <- ALB * 0.1

    # 2. Individual parameters, following Rubino 2021 Table 2 footnote a:
    #    CL  = 0.254 (1 - 0.133 SEXF) (albumin/4.2)^-0.844
    #    Vc  = [11.1 + 7.54 (BSA - 1.83)] (1 + 0.478 infected)
    #    Vp1 = 14.6 (albumin/4.2)^-0.829 (BSA/1.83)^1.14
    #    Vp2 = 6.69 (1 + 1.69 infected) (BSA/1.83)^2.47
    #    The Vc BSA effect is an additive slope in litres, so it is added to
    #    the typical value exp(lvc) before the exponential eta is applied.
    cl <- exp(lcl + etalcl) * (1 + e_sexf_cl * SEXF) * (alb_gdL / 4.2)^e_alb_cl
    vc <- (exp(lvc) + e_bsa_vc * (BSA - 1.83)) * exp(etalvc) * (1 + e_infect_vc * DIS_INFECT_ACTIVE)
    q <- exp(lq)
    vp <- exp(lvp + etalvp) * (alb_gdL / 4.2)^e_alb_vp * (BSA / 1.83)^e_bsa_vp
    q2 <- exp(lq2)
    vp2 <- exp(lvp2 + vp2_eta_scale * etalvp) * (1 + e_infect_vp2 * DIS_INFECT_ACTIVE) * (BSA / 1.83)^e_bsa_vp2
    q3 <- exp(lq3)
    vp3 <- exp(lvp3)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2
    k14 <- q3 / vc
    k41 <- q3 / vp3

    # 4. Four-compartment disposition with first-order elimination; doses
    #    enter central as a zero-order IV infusion from the event table.
    d/dt(central) <- -kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2 -
      k14 * central + k41 * peripheral3
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(peripheral3) <- k14 * central - k41 * peripheral3

    # 5. Observation: mg / L = mg/L (plasma rezafungin).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
