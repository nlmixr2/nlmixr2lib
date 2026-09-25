Hahn_2019_sufentanil <- function() {
  description <- paste(
    "Two-compartment IV population PK model for sufentanil in critically ill",
    "adults supported with venoarterial extracorporeal membrane oxygenation",
    "(VA-ECMO) after myocardial infarction. First-order elimination, no",
    "absorption (continuous IV infusion into the central compartment). Two",
    "covariates were retained after forward selection and backward",
    "elimination: an exponential effect of tympanic body temperature on",
    "clearance (+23.0% per degC, centered at the cohort median 36.9 degC) and",
    "a power effect of total plasma protein on the peripheral volume",
    "(exponent 2.46, centered at the cohort median 4.5 g/dL = 45 g/L).",
    "Inter-individual variability was estimated on CL and V2 only.",
    "Proportional residual error 29.0%. Relative to non-ECMO reference data",
    "the volumes are markedly increased and clearance decreased, consistent",
    "with circuit sequestration of this lipophilic, highly protein-bound drug",
    "and with reduced hepatic blood flow in critical illness. Parameter",
    "values from Hahn 2019 Table 2 (final model column)."
  )
  reference <- paste(
    "Hahn J, Yang S, Min KL, Kim D, Jin BH, Park C, Park MS, Wi J, Chang MJ.",
    "Population pharmacokinetics of intravenous sufentanil in critically ill",
    "patients supported with extracorporeal membrane oxygenation therapy.",
    "Critical Care 2019;23:248. doi:10.1186/s13054-019-2508-4.",
    sep = " "
  )
  vignette <- "Hahn_2019_sufentanil"
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "ug/L"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Blood was drawn from an existing arterial line into
  # EDTA tubes, centrifuged, and the plasma assayed for sufentanil by
  # HPLC-MS/MS (Methods "Sample collection and plasma concentration assay"),
  # so the observed matrix is plasma.
  compartmentData <- list(
    central = list(analyte = "sufentanil", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "sufentanil", units = "ug", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    BODYTEMP = list(
      description = "Tympanic body temperature",
      units = "degC",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying within subject in principle -- body temperature is an",
        "ICU vital sign recorded serially from the electronic medical record",
        "(Methods 'Dosing, administration, and data collection') -- and the",
        "paper does not state whether the fitted covariate was an admission",
        "value or a per-record value. Hahn 2019 cohort median (range) 36.9",
        "(33-38.7) degC (Table 1). Applied as an exponential effect on",
        "systemic clearance centered at the cohort median:",
        "CL = 37.8 * exp(e_bodytemp_cl * (BODYTEMP - 36.9)) L/h with",
        "e_bodytemp_cl = 0.207 per degC (Results 'Population PK model",
        "building' final-model equation; Table 2 row Theta Temp). Clearance",
        "rises about 23.0% per degC, i.e. it falls in hypothermia. The",
        "Discussion attributes the direction to reduced total hepatic blood",
        "flow at low temperature (sufentanil hepatic extraction ratio ~0.7)",
        "and to slowed CYP3A4 metabolic activity. The paper's Monte Carlo",
        "simulations span 33, 35, 36.7, 38 and 39 degC (Methods 'Model",
        "evaluation and simulations'), which is the range over which the",
        "effect was intended to be used."
      ),
      source_name = "temperature"
    ),
    TPRO = list(
      description = "Total plasma protein",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Hahn 2019 reports total plasma protein in the US convention g/dL:",
        "cohort median (range) 4.5 (2.1-6) g/dL (Table 1). The canonical",
        "column is SI g/L, so the reference is 45 g/L and user data must be",
        "supplied in g/L (1 g/dL = 10 g/L). Because the effect is a ratio to",
        "the reference, (TPRO/45) in g/L equals the paper's",
        "(total plasma protein/4.5) in g/dL exactly and no inline unit",
        "conversion is needed. Applied as a power effect on the peripheral",
        "volume of distribution:",
        "V2 = 1640 * (TPRO/45)^e_tpro_vp L with e_tpro_vp = 2.46 (Results",
        "'Population PK model building' final-model equation; Table 2 row",
        "Theta T.Prot). The positive exponent means higher total protein",
        "predicts a larger peripheral volume -- the opposite direction to",
        "several other drugs, which the Discussion acknowledges explicitly,",
        "proposing that sufentanil binding is driven mainly by alpha-1 acid",
        "glycoprotein rather than total protein and that low total protein",
        "may instead mark impaired hepatic function. The Discussion also",
        "cautions that V2 shrinkage was relatively high (33%) and that the",
        "relationship is observational. The paper's Monte Carlo simulations",
        "span 2, 4, 6 and 8 g/dL (= 20, 40, 60, 80 g/L), which extends above",
        "the observed cohort maximum of 6 g/dL."
      ),
      source_name = "total plasma protein"
    )
  )

  covariatesDataExcluded <- list(
    TBIL = list(
      description = "Total bilirubin",
      units = "mg/dL",
      type = "continuous",
      notes = paste(
        "Significant on CL in the univariate covariate screen (drop in OFV",
        "of 6.0 points) but removed during backward elimination and not",
        "retained in the final model (Results 'Population PK model",
        "building'). No point estimate is reported for the effect. Cohort",
        "median (range) 1.9 (0.3-6.6) mg/dL (Table 1)."
      )
    ),
    LBW = list(
      description = "Lean body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Significant on V2 in the univariate covariate screen (drop in OFV",
        "of 3.4 points) but removed during backward elimination and not",
        "retained in the final model (Results 'Population PK model",
        "building'). No point estimate is reported for the effect. Cohort",
        "median (range) 55.3 (36.8-58.6) kg (Table 1). The Discussion notes",
        "that weight-related covariates were not retained and attributes",
        "this to the large impact of ECMO on sufentanil PK, consistent with",
        "prior cardiopulmonary-bypass sufentanil studies."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20L,
    n_studies = 1L,
    n_observations = 106L,
    age_range = "23-88 years",
    age_median = "55 years",
    weight_range = "52.9-92.5 kg",
    weight_median = "69.4 kg",
    lean_body_weight_range = "36.8-58.6 kg",
    bmi_range = "20.5-31.8 kg/m^2",
    sex_female_pct = 20,
    race_ethnicity = c(Asian = 100),
    disease_state = paste(
      "Critically ill adults (>=19 years) receiving sufentanil-based",
      "analgesia and sedation during venoarterial ECMO for myocardial",
      "infarction: ST-elevation MI (12), acute MI (4) and non-ST-elevation",
      "MI (4). All received mechanical ventilation and started ECMO within",
      "12 h of MI onset. Median (range) APACHE II 29 (15-36); median VA-ECMO",
      "duration 138 (52.9-263) h; median ECMO flow rate 3 (0.6-4.1) L/min.",
      "Nine of 20 concurrently received continuous venovenous",
      "hemodiafiltration. Patients with known sufentanil allergy or taking",
      "interacting medications were excluded."
    ),
    baseline_labs = paste(
      "Median (range): total plasma protein 4.5 (2.1-6) g/dL, total",
      "bilirubin 1.9 (0.3-6.6) mg/dL, blood urea nitrogen 21.3 (7.5-58)",
      "mg/dL, serum creatinine 1.4 (0.4-4.9) mg/dL, partial pressure of",
      "carbon dioxide 29.1 (13.5-46.7) mmHg, tympanic body temperature 36.9",
      "(33-38.7) degC (Table 1)."
    ),
    dose_range = paste(
      "Continuous IV infusion of sufentanil, initial rate 12.5 ug/h for",
      "patients under 60 kg (5 of 20) or 17.5 ug/h for patients at or above",
      "60 kg (15 of 20), titrated to a target Richmond Agitation Sedation",
      "Scale score. Median (range) infusion duration 110 (34-260) h.",
      "Midazolam was co-administered as needed (3 or 5 mg IV bolus then",
      "4.5 mg/h infusion)."
    ),
    co_medication = "Midazolam for supplemental sedation in all patients.",
    regions = "South Korea (single-centre, Severance Cardiovascular Hospital, Seoul)",
    notes = paste(
      "Prospective cohort PK study conducted January 2016 to June 2017",
      "(IRB 4-20140919; ClinicalTrials.gov NCT02581280). Sampling began",
      "within the first 48 h of ECMO: after 3 and 12 h of infusion then",
      "every 24 h to 96 h; after infusion cessation at 0, 0.5, 1, 2, 6 and",
      "12 h then every 24 h to 72 h. 106 plasma samples from 20 patients.",
      "Concentrations below the 0.02 ug/L lower limit of quantification were",
      "excluded from the analysis. Estimation by FOCE+I in NONMEM 7.4."
    )
  )

  ini({
    # Fixed effects -- Hahn 2019 Table 2, "Final model (RSE%)" column, and the
    # final-model equations printed in Results "Population PK model building".
    # Typical values apply at the covariate reference point: body temperature
    # 36.9 degC and total plasma protein 4.5 g/dL (45 g/L), both cohort medians.
    lcl <- log(37.8); label("Systemic clearance CL (L/h)")                         # Hahn 2019 Table 2, Theta CL row: 37.8 (RSE 3%)
    lvc <- log(229); label("Central volume of distribution V1 (L)")                # Hahn 2019 Table 2, Theta V1 row: 229 (RSE 10%)
    lq <- log(41); label("Intercompartmental clearance Q (L/h)")                   # Hahn 2019 Table 2, Theta Q row: 41 (RSE 12%)
    lvp <- log(1640); label("Peripheral volume of distribution V2 (L)")            # Hahn 2019 Table 2, Theta V2 row: 1640 (RSE 9%)

    # Covariate effects -- Hahn 2019 Table 2 and the final-model equations.
    e_bodytemp_cl <- 0.207; label("Exponential coefficient of (BODYTEMP - 36.9) on CL (per degC)") # Hahn 2019 Table 2, Theta Temp row: 0.207 (RSE 5%)
    e_tpro_vp <- 2.46; label("Power exponent of TPRO/45 on Vp (unitless)")         # Hahn 2019 Table 2, Theta T.Prot row: 2.46 (RSE 7%)

    # Inter-individual variability -- Hahn 2019 Table 2, 'Inter-individual
    # variability' section. IIV was retained on CL and V2 only (Results:
    # 'IIVs were included for CL and V2, since they significantly improved
    # model performance'). The rows are labelled omega^2, i.e. the reported
    # numbers are NONMEM log-normal VARIANCES on theta_i = theta_pop * exp(eta_i)
    # (Methods 'Population PK model development'), and are used here unchanged.
    etalcl ~ 0.167 # Hahn 2019 Table 2, row 'omega CL 2' final model = 0.167 (RSE 57%, shrinkage 20%)
    etalvp ~ 1.13 # Hahn 2019 Table 2, row 'omega V2 2' final model = 1.13 (RSE 48%, shrinkage 22%)

    # Residual variability -- proportional only: c_ij = cp_ij * (1 + eps_ij)
    # with var(eps) = sigma^2 (Methods 'Population PK model development').
    # Table 2 reports sigma^2 = 0.0841, so the SD nlmixr2 expects is
    # sqrt(0.0841) = 0.29 exactly.
    propSd <- 0.29; label("Proportional residual error (fraction)")                # Hahn 2019 Table 2, row 'sigma 2 proportional' final model = 0.0841; sqrt(0.0841) = 0.29
  })
  model({
    # 1. Individual PK parameters with the Hahn 2019 final-model covariate
    # equations (Results "Population PK model building"):
    #   CL = 37.8 * EXP(0.207 * (temperature - 36.9)) L/h
    #   V1 = 229 L
    #   V2 = 1640 * (total plasma protein / 4.5)^2.46 L
    #   Q  = 41 L/h
    # TPRO is carried in SI g/L, so the reference 45 g/L is the paper's
    # 4.5 g/dL; the ratio is identical in either unit.
    cl <- exp(lcl + etalcl) * exp(e_bodytemp_cl * (BODYTEMP - 36.9))
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp) * (TPRO / 45)^e_tpro_vp

    # 2. Micro-constants (1/h with cl, q in L/h and vc, vp in L).
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. Two-compartment disposition with IV dosing only. Sufentanil was given
    # as a continuous IV infusion into the central compartment; there is no
    # absorption step (Methods "Dosing, administration, and data collection").
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Observation and error. central in ug, vc in L -> Cc in ug/L, matching
    # the paper's 0.3-0.6 ug/L target concentration band and the 0.02 ug/L
    # lower limit of quantification.
    Cc <- central / vc

    Cc ~ prop(propSd)
  })
}
