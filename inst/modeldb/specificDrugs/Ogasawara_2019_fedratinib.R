Ogasawara_2019_fedratinib <- function() {
  description <- "Two compartment oral PK model of fedratinib with first-order absorption and a lag time in patients with myelofibrosis, polycythemia vera or essential thrombocythemia (Ogasawara 2019)"
  reference <- "Ogasawara K, Zhou S, Krishna G, Palmisano M, Li Y. Population pharmacokinetics of fedratinib in patients with myelofibrosis, polycythemia vera, and essential thrombocythemia. Cancer Chemother Pharmacol. 2019;84(4):707-718. doi:10.1007/s00280-019-03929-9"
  vignette <- "Ogasawara_2019_fedratinib"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "fedratinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "fedratinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "fedratinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    DIS_PV = list(
      description = "Polycythemia vera disease-state indicator (1 = polycythemia vera)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (myelofibrosis or essential thrombocythemia)",
      notes = paste(
        "Multiplicative factors of 1.54 on CL/F and 1.87 on V2/F for polycythemia vera patients,",
        "per Ogasawara 2019 Table 2 footnotes b and c. The reference group pools primary myelofibrosis,",
        "post-polycythemia-vera myelofibrosis, post-essential-thrombocythemia myelofibrosis and essential",
        "thrombocythemia (90.0 percent of the 452 subjects); only the 45 subjects with active polycythemia vera",
        "carry a value of 1. Patients whose polycythemia vera has progressed to myelofibrosis take the value 0:",
        "the paper's Discussion reports that post-PV myelofibrosis CL/F returns to the myelofibrosis level,",
        "so the indicator marks the active polycythemia vera phenotype rather than a polycythemia vera history.",
        "Time-fixed per subject."
      ),
      source_name = "PV"
    ),
    CRCL = list(
      description = "Creatinine clearance, raw (not body-surface-area normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL/F with exponent 0.294, normalized to 78.3 mL/min per the Ogasawara 2019",
        "Table 2 footnote b equation. Raw un-normalized creatinine clearance in mL/min, NOT the",
        "mL/min/1.73 m^2 body-surface-area-normalized default of this canonical - supply values on",
        "the raw scale. The analysis dataset median was 78.5 mL/min (range 20.1 to 181), so the",
        "78.3 normalization constant printed in the equation differs very slightly from the Table 1",
        "median; the printed equation governs (a 0.075 percent difference in the resulting CL/F).",
        "Missing baseline values were imputed at the study-population median by the authors."
      ),
      source_name = "CLcr"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on V2/F with exponent 0.727, normalized to the study-population median of 70.1 kg (Ogasawara 2019 Table 2 footnote c). Baseline weight; observed range 39.5 to 135 kg.",
      source_name = "Weight"
    ),
    DOSE = list(
      description = "Administered fedratinib dose level",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Use case (a): per-subject assigned once-daily dose level, entering V2/F as the power term",
        "(DOSE / 400)^-0.279 per Ogasawara 2019 Table 2 footnote c. The final model was fit to doses of",
        "100 mg and above only; the authors excluded the 30 to 60 mg data because CL/F and V2/F were",
        "dose-dependent below about 120 mg. Supply DOSE in mg matching the amt of the dose records."
      ),
      source_name = "Dose"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "year",
      type = "continuous",
      notes = "Screened in the stepwise covariate search; no clinically meaningful effect on fedratinib PK over 20 to 95 years and not retained in the final model (Ogasawara 2019 Discussion)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened in the stepwise covariate search; no clinically meaningful effect and not retained in the final model (Ogasawara 2019 Discussion)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Component of the NCI-ODWG liver-function classification screened as a covariate; neither mild nor moderate hepatic impairment had a clinically meaningful effect and no hepatic term was retained (Ogasawara 2019 Discussion)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Component of the NCI-ODWG liver-function classification screened as a covariate; no hepatic term was retained in the final model (Ogasawara 2019 Discussion)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 452,
    n_studies = 6,
    n_observations = 3442,
    age_range = "20-95 years",
    age_median = "65 years",
    weight_range = "39.5-135 kg",
    weight_median = "70.1 kg",
    sex_female_pct = 44.9,
    race_ethnicity = c(Caucasian = 88.3, `African-American` = 1.5, Asian = 9.7, Other = 0.4),
    disease_state = "Myeloproliferative neoplasms: primary myelofibrosis (51.3 percent), post-polycythemia-vera myelofibrosis (20.1 percent), post-essential-thrombocythemia myelofibrosis (11.3 percent), polycythemia vera (10.0 percent) and essential thrombocythemia (7.3 percent)",
    dose_range = "100-800 mg orally once daily; the final model was fit to doses of 100 mg and above, and 400 mg QD was the modal regimen (51.1 percent of subjects)",
    renal_function = "Median creatinine clearance 78.5 mL/min (range 20.1 to 181); no subjects with latent renal impairment were enrolled",
    hepatic_function = "NCI-ODWG classification: 320 normal, 115 mild, 17 moderate; no severe hepatic impairment",
    studies = "TED12037 (phase 1, NCT00631462), ARD11936 (NCT01420770), ARD12042 (NCT01420783), ARD12181 (NCT01523171) and ARD12188 (Japanese, NCT01692366) phase 2 studies, and EFC12153 (phase 3, NCT01437787)",
    notes = "Demographics from Ogasawara 2019 Table 1; per-study designs, dosing regimens and PK sampling schedules from Supplementary Table 1. Estimated in NONMEM 7.3.0 by FOCE with interaction on natural-log-transformed concentrations."
  )

  ini({
    # Structural parameters - apparent (oral) parameters, so all volumes and
    # clearances carry the /F of the published parameterization.
    lka <- log(1.57); label("Absorption rate constant (1/h)") # Table 2, TVKa
    lcl <- log(13.0); label("Apparent clearance CL/F (L/h)") # Table 2, TVCL/F
    lvc <- log(311); label("Apparent central volume of distribution V2/F (L)") # Table 2, TVV2/F
    lq <- log(45.2); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2, TVQ/F
    lvp <- log(1460); label("Apparent peripheral volume of distribution V3/F (L)") # Table 2, TVV3/F
    ltlag <- log(0.265); label("Absorption lag time (h)") # Table 2, TVALAG1

    # Covariate effects on CL/F - Table 2 footnote b:
    #   CL/F (L/h) = 13.0 * 1.54(if PV) * (CLcr/78.3)^0.294
    e_dis_pv_cl <- 1.54; label("Multiplicative factor on CL/F for polycythemia vera (unitless)") # Table 2, PV on CL/F
    e_crcl_cl <- 0.294; label("Power exponent of creatinine clearance on CL/F (unitless)") # Table 2, CLcr on CL/F

    # Covariate effects on V2/F - Table 2 footnote c:
    #   V2/F (L) = 311 * 1.87(if PV) * (Weight/70.1)^0.727 * (Dose/400)^-0.279
    e_dis_pv_vc <- 1.87; label("Multiplicative factor on V2/F for polycythemia vera (unitless)") # Table 2, PV on V2/F
    e_wt_vc <- 0.727; label("Power exponent of body weight on V2/F (unitless)") # Table 2, Weight on V2/F
    e_dose_vc <- -0.279; label("Power exponent of dose level on V2/F (unitless)") # Table 2, Dose on V2/F

    # Inter-individual variability. nlmixr2 stores the variance (omega^2), which
    # is what Table 2 reports directly. CL/F and V2/F share an estimated
    # covariance, so they are declared as a block; ka is independent.
    # Implied coefficients of variation, sqrt(exp(omega^2) - 1): CL/F 53.9,
    # V2/F 67.7 and ka 133.9 per cent.
    etalcl + etalvc ~ c(0.255, 0.197, 0.383) # Table 2 random effects: omega2 CL/F, COV CL/F-V2/F, omega2 V2/F
    etalka ~ 1.07 # Table 2 random effects: omega2 Ka

    # Residual error. Concentrations were natural-log transformed and the
    # residual was additive on that log scale, i.e. log-normal in linear space.
    # Table 2 reports the variance sigma^2 = 0.201, so the SD is sqrt(0.201).
    expSd <- 0.448; label("Log-scale additive residual error SD (unitless)") # Table 2, sigma2 Log additive 0.201; sqrt(0.201) = 0.448
  })

  model({
    # Individual parameters. Reference values are the normalization constants
    # printed in the Table 2 footnote equations: CLcr 78.3 mL/min, weight
    # 70.1 kg and dose 400 mg. The binary polycythemia vera effects are applied
    # as factor^DIS_PV so they contribute 1 for the myelofibrosis / essential
    # thrombocythemia reference group.
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag)

    cl <- exp(lcl + etalcl) *
      e_dis_pv_cl^DIS_PV *
      (CRCL / 78.3)^e_crcl_cl

    vc <- exp(lvc + etalvc) *
      e_dis_pv_vc^DIS_PV *
      (WT / 70.1)^e_wt_vc *
      (DOSE / 400)^e_dose_vc

    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    Cc <- central / vc # mg/L
    Cc ~ lnorm(expSd)
  })
}
