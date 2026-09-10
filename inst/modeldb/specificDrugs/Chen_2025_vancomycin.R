Chen_2025_vancomycin <- function() {
  description <- "One-compartment IV population PK model for vancomycin in adult intensive-care patients with sepsis, developed from 11,046 routine therapeutic-drug-monitoring concentrations in 4,006 MIMIC-IV patients (Chen 2025). Clearance is a near-linear power function of Cockcroft-Gault creatinine clearance (exponent 0.997, reference 93 mL/min) multiplied by an exponential Charlson-Comorbidity-Index term; volume of distribution is a shallow power function of body weight (exponent 0.205, reference 84 kg). This is the population-PK arm of a four-way comparison (PPK, Bayesian, random forest, hybrid PPK-ML) of AUC24 prediction; only the PPK structural model is a pharmacokinetic model and only it is packaged here. Both between-subject variances and both residual-error magnitudes are absent from the article and its supplement and are encoded as zero."
  reference <- "Chen K, Wang C, Wei Y, Ma S, Huang W, Dong Y, Wang Y. Machine learning and population pharmacokinetics: a hybrid approach for optimizing vancomycin therapy in sepsis patients. Microbiol Spectr. 2025;13(5):e00499-25. doi:10.1128/spectrum.00499-25"
  vignette <- "Chen_2025_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Vancomycin was given intravenously (an inclusion
  # criterion was "documented intravenous vancomycin use during their ICU
  # stay"), so the dose enters `central` directly and there is no depot
  # state. Chen 2025 Methods "Data collection" says only "vancomycin
  # treatment regimen and blood concentrations"; the MIMIC-IV vancomycin
  # assay is a serum measurement, but the paper never states serum versus
  # plasma, so `specimen` is left unverified.
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = FALSE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CLCR. Chen 2025 equation 15 enters it as the power term (CLCR/93)^0.997. The reference 93 mL/min is the Table 1 cohort mean creatinine clearance (overall 93.25, SD 66.15; training set 93.28, SD 66.42; testing set 91.19, SD 41.74). Chen 2025 does not name the estimating equation or state any BSA normalization -- Methods 'Data collection' lists only 'laboratory data, including creatinine, creatinine clearance, and red blood cell count' as extracted from MIMIC-IV 2.2, where the derived creatinine-clearance field is Cockcroft-Gault in raw mL/min. Stored under the canonical CRCL column in raw mL/min, following the Alqahtani_2018_vancomycin.R, Buelga_2005_vancomycin.R and Delattre_2010_amikacin.R precedents for raw Cockcroft-Gault (see inst/references/covariate-columns.md, CRCL entry). The exponent is 0.997 (bootstrap 95% CI 0.94-1.04), i.e. statistically indistinguishable from direct proportionality, so clearance is effectively linear in creatinine clearance over the observed range.",
      source_name        = "CLCR"
    ),
    SCORE_CCI = list(
      description        = "Charlson Comorbidity Index total score",
      units              = "(SCORE_CCI units, weighted comorbidity count)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CCI. Chen 2025 equation 15 enters it as exp(-0.151 * (CCI/5.62)); the whole product -0.151 * (CCI/5.62) sits in the exponent, confirmed against the publisher's typeset equation image (spectrum.00499-25.m015.jpg in the EuropePMC supplementary bundle for PMC12054080), because pdftotext flattens the superscript and makes the nesting ambiguous. The scaling constant 5.62 is the Table 1 TRAINING-set mean CCI (5.62, SD 3.10; overall cohort mean 5.63, SD 3.10; testing set 5.91, SD 2.84). Note that the term is a plain ratio and is NOT centred: at the cohort mean CCI of 5.62 the term evaluates to exp(-0.151) = 0.860, not 1, so the typical clearance of an average-comorbidity patient is 0.860 * 3.35 = 2.88 L/h rather than the 3.35 L/h that Chen 2025 Results calls 'the typical value of the CL population'. The model is encoded exactly as the equation prints; see the vignette Errata. Higher comorbidity burden lowers vancomycin clearance (coefficient -0.151, bootstrap 95% CI -0.21 to -0.09).",
      source_name        = "CCI"
    ),
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column WT, defined by Chen 2025 Results as 'patient weight in kilograms'. Chen 2025 equation 16 enters it as the power term (WT/84)^0.205 on volume of distribution. Table 1 cohort mean weight is 82.56 kg (SD 25.35; training set 82.49, SD 25.34; testing set 87.49, SD 26.05); the paper does not state where the reference 84 kg comes from, and it matches neither the reported overall mean nor either subset mean exactly -- most likely the training-set median, which Chen 2025 does not report. Used as printed. Chen 2025 Limitations note that height was missing for 61.7% of MIMIC-IV patients, so BMI could not be computed and no adjusted body weight was available for obese patients; total body weight is therefore the only size descriptor in the model. NOTE: Table S1 labels the 0.205 exponent row 'WT on CL', which contradicts equation 16 placing it on V; the equation is taken as authoritative (see vignette Errata).",
      source_name        = "WT"
    )
  )

  # Collected by Chen 2025 (Methods 'Data collection', tabulated in Table 1)
  # but NOT retained in the final PPK model, whose only covariates are CLCR
  # and CCI on CL and WT on V (equations 15-16). Chen 2025 states that "the
  # stepwise regression method was used to screen covariates to establish the
  # final PPK model" but never enumerates the candidate set or reports the
  # dropped terms, so these are recorded as collected-and-not-retained rather
  # than as confirmed screening failures. Several of them WERE retained by
  # the paper's random-forest and hybrid machine-learning models (Results
  # 'Machine learning model': TAD, DOSE24, WT, age, hematocrit, red blood
  # cell count, hemoglobin and CLCR for the random forest; TAD, DOSE24, CL,
  # V, age, WT, hematocrit and hemoglobin for the hybrid), but those are
  # regression predictors of concentration, not pharmacokinetic covariates,
  # and the ML models are not packaged here.
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Subject age at ICU admission",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 1 overall mean 65.39 years (SD 15.81). Inclusion required age over 18 years at ICU admission. Not retained in the final PPK model; selected as a predictor in both the random-forest and hybrid ML models.",
      source_name        = "Age"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "male",
      notes              = "Table 1 reports male percentage: overall 58.81% male, so 41.19% female. Not retained in the final PPK model.",
      source_name        = "Male"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 1 overall mean 1.36 mg/dL (SD 1.27). Not retained in the final PPK model; it enters indirectly as the input to the CLCR field that is retained.",
      source_name        = "Creatinine"
    ),
    HCT = list(
      description        = "Hematocrit",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 1 overall mean 29.57% (SD 5.17). Not retained in the final PPK model; selected as a predictor in both the random-forest and hybrid ML models.",
      source_name        = "Hematocrit"
    ),
    RBC = list(
      description        = "Red blood cell count",
      units              = "10^6/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 1 overall mean 3.24 (SD 0.62). Not retained in the final PPK model; selected as a predictor in the random-forest ML model.",
      source_name        = "Red blood cell"
    ),
    HGB = list(
      description        = "Hemoglobin",
      units              = "g/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 1 overall mean 9.51 g/dL (SD 1.73). Not retained in the final PPK model; selected as a predictor in both the random-forest and hybrid ML models.",
      source_name        = "Hemoglobin"
    ),
    WBC = list(
      description        = "White blood cell count",
      units              = "10^3/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 1 overall mean 13.68 (SD 9.29). Not retained in the final PPK model.",
      source_name        = "White blood cell"
    ),
    PLT = list(
      description        = "Platelet count",
      units              = "10^3/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 1 overall mean 215.36 (SD 117.81). Not retained in the final PPK model.",
      source_name        = "Platelet"
    ),
    LACT = list(
      description        = "Serum lactate",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 1 overall mean 2.16 mmol/L (SD 1.47). Not retained in the final PPK model.",
      source_name        = "Lactate"
    ),
    PT = list(
      description        = "Prothrombin time",
      units              = "s",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 1 overall mean 17.47 s (SD 8.72). Not retained in the final PPK model.",
      source_name        = "Prothrombin time"
    ),
    SOFA = list(
      description        = "Sequential Organ Failure Assessment score",
      units              = "(score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 1 overall mean 6.21 (SD 3.43). The only baseline variable other than hospital length of stay that differed significantly between the training and testing sets (6.22 vs 5.36, P < 0.05). Not retained in the final PPK model.",
      source_name        = "SOFA"
    ),
    APSIII = list(
      description        = "Acute Physiology Score III",
      units              = "(score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 1 overall mean 51.07 (SD 19.47). Not retained in the final PPK model.",
      source_name        = "APSIII"
    ),
    SAPS_II = list(
      description        = "Simplified Acute Physiology Score II",
      units              = "(score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 1 overall mean 40.17 (SD 13.43); the Table 1 row label is misspelled 'SQPSII'. Not retained in the final PPK model.",
      source_name        = "SQPSII"
    ),
    GCS = list(
      description        = "Glasgow Coma Scale score",
      units              = "(score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 1 overall mean 13.38 (SD 2.88). Not retained in the final PPK model.",
      source_name        = "GCS"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 4059L,
    n_studies        = 1L,
    n_sites          = 1L,
    n_concentrations = 11046L,
    age_range        = "adults over 18 years at ICU admission",
    age_mean         = "65.39 years (SD 15.81)",
    weight_mean      = "82.56 kg (SD 25.35)",
    sex_female_pct   = 41.19,
    race_ethnicity   = "White 65.31%, Other 21.11%, Black 10.89%, Asian 2.69% (Table 1, overall cohort)",
    disease_state    = "Adults admitted to intensive care and meeting the Sepsis-3 definition, receiving intravenous vancomycin, with an ICU stay longer than 24 hours and at least one vancomycin concentration measurement. Excluded: patients receiving renal replacement therapy, pregnant patients, patients with no dosing information before a concentration measurement, and patients who died within 48 hours. Mean SOFA 6.21, mean Charlson Comorbidity Index 5.63, mean APS III 51.07; ICU mortality 11.92%, in-hospital mortality 18.48%; mean ICU stay 7.17 days and mean hospital stay 15.62 days.",
    dose_range       = "Not reported. The 24 hour vancomycin dose (DOSE24) was a predictor in the machine-learning arms of the paper, but Chen 2025 reports neither its distribution nor any protocolized regimen; dosing was routine clinical care recorded in MIMIC-IV.",
    regions          = "United States (Beth Israel Deaconess Medical Center, Boston, via the MIMIC-IV 2.2 database)",
    renal_function   = "Creatinine clearance mean 93.25 mL/min (SD 66.15); serum creatinine mean 1.36 mg/dL (SD 1.27). Patients on renal replacement therapy were excluded.",
    notes            = "Retrospective analysis of routine electronic-health-record data from MIMIC-IV 2.2. Of 4,059 eligible patients, the 53 who had both a peak (1-2 h after dosing) and a trough (30 min to 1 h before the next dose) within the SAME dosing interval were held out as the testing set and the remaining 4,006 patients, contributing 11,046 concentrations, were the PPK training set; the parameter estimates packaged here come from that 4,006-patient training fit. Variables with more than 20% missing data were dropped and the remainder imputed by classification-and-regression-trees multiple imputation (R MICE package). Fit in NONMEM 7.5; covariates selected by stepwise regression; evaluated by goodness-of-fit diagnostics (Figure S1) and a bootstrap (Table S1). Chen 2025's wider purpose is a four-way comparison of AUC24 prediction between this PPK model, a Bayesian posterior using it, a random forest, and a hybrid model feeding the PPK individual CL and V into the random forest; in the 53-patient testing set the PPK model alone performed worst (MAPE 68.17%, F30 34.6%) and the Bayesian posterior best (MAPE 13.37%, F30 94.2%). Only the PPK structural model is a pharmacokinetic model and only it is packaged here."
  )

  ini({
    # Structural parameters, Chen 2025 equations 15-16, with estimates and
    # bootstrap results in supplementary Table S1 (file
    # spectrum.00499-25-s0001.pdf).
    #
    # UNITS: Chen 2025 Results says "the typical value of the CL population
    # was 3.35 mL/min", but supplementary Table S1 labels the same row
    # "CL (L*h^-1) 3.35 (3%)". L/h is correct and mL/min is a units error in
    # the main text: 3.35 mL/min = 0.201 L/h would give vancomycin a terminal
    # half-life of 340 h, and Figure 4 shows the model's own predicted AUC24
    # spanning roughly 250-1150 ug*h/mL, which a ~3 L/h clearance reproduces
    # (a 2,000 mg daily dose over a clearance of 2.88 L/h gives 694
    # ug*h/mL) and a 0.2 L/h clearance does not. See vignette Errata.
    #
    # lcl is the clearance of a patient with CRCL = 93 mL/min and a Charlson
    # Comorbidity Index of ZERO, not of a typical patient: the CCI term in
    # equation 15 is an uncentred ratio, so it contributes exp(-0.151) =
    # 0.860 at the cohort mean CCI of 5.62.
    lcl <- log(3.35);  label("Clearance at CRCL=93 mL/min and SCORE_CCI=0 (CL, L/h)")  # Chen 2025 Table S1 "CL (L*h^-1)" = 3.35 (RSE 3%); bootstrap median 3.33, 95% CI 3.15-3.54, bias -1.0%
    lvc <- log(98.5);  label("Volume of distribution at WT=84 kg (V, L)")              # Chen 2025 Table S1 "V (L)" = 98.5 (RSE 1%); bootstrap median 98.80, 95% CI 95.46-101.59, bias 2.7%

    # Covariate effects. All three are estimated (Table S1 reports an RSE and
    # a bootstrap confidence interval for each), so none is wrapped in
    # fixed().
    e_crcl_cl <- 0.997;  label("Power exponent on (CRCL/93 mL/min) for CL (unitless)")            # Chen 2025 Table S1 "CLCR on CL" = 0.997 (RSE 2%); bootstrap median 0.99, 95% CI 0.94-1.04, bias -0.6%
    e_cci_cl  <- -0.151; label("Exponential coefficient on (SCORE_CCI/5.62) for CL (unitless)")   # Chen 2025 Table S1 "Charlson on CL" = -0.151 (RSE 17%); bootstrap median -0.15, 95% CI -0.21 to -0.09, bias -0.2%
    e_wt_vc   <- 0.205;  label("Power exponent on (WT/84 kg) for V (unitless)")                   # Chen 2025 equation 16; Table S1 row (mislabelled "WT on CL") = 0.205 (RSE 23%); bootstrap median 0.20, 95% CI 0.11-0.30, bias -0.1%

    # Between-subject variability. Chen 2025 equation 1 declares an
    # exponential IIV model with variance omega^2, and equations 15-16 carry
    # exp(eta1) on CL and exp(eta2) on V, so the STRUCTURE of the random
    # effects is fully specified. No omega VALUE appears anywhere in the
    # article or in its supplement: Table S1 lists only the five fixed
    # effects above, and the paper reports no omega, no CV%, and no shrinkage.
    # Both variances are therefore encoded as zero rather than invented, so
    # simulations from this model are population-typical unless the user
    # supplies their own omega.
    etalcl ~ fixed(0)  # Chen 2025 equations 1 and 15 (exponential IIV on CL declared; magnitude not reported)
    etalvc ~ fixed(0)  # Chen 2025 equations 1 and 16 (exponential IIV on V declared; magnitude not reported)

    # Residual variability. Chen 2025 Results "PPK model" states that "the
    # hybrid model described the residual variation in equation 4", i.e. the
    # combined error model Cobs = Cpred * (1 + eps) + eps1, which maps onto
    # nlmixr2's prop(propSd) + add(addSd). ("Hybrid" here names the combined
    # RESIDUAL-error model of equation 4, not the paper's hybrid PPK-ML
    # model.) As with the omegas, neither delta^2 nor delta1^2 is reported
    # anywhere in the article or supplement, so both magnitudes are encoded
    # as zero rather than invented.
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- not reported in the source)")  # Chen 2025 equation 4 selected; no residual-error estimate published
    addSd  <- fixed(0); label("Additive residual SD (ug/mL; 0 -- not reported in the source)")         # Chen 2025 equation 4 selected; no residual-error estimate published
  })
  model({
    # Individual PK parameters. Chen 2025 equations 15-16:
    #   CL (L/h) = 3.35 * (CLCR/93)^0.997 * exp(-0.151 * (CCI/5.62)) * exp(eta1)
    #   V  (L)   = 98.5 * (WT/84)^0.205 * exp(eta2)
    # The CCI term is an uncentred ratio inside the exponential, exactly as
    # the published equation prints it (verified against the publisher's
    # typeset equation image spectrum.00499-25.m015.jpg).
    cl <- exp(lcl + etalcl) * (CRCL / 93)^e_crcl_cl * exp(e_cci_cl * (SCORE_CCI / 5.62))
    vc <- exp(lvc + etalvc) * (WT / 84)^e_wt_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, volume in L, so central/vc is mg/L == ug/mL, the units in
    # which vancomycin concentrations and the paper's AUC24 (ug*h/mL) are
    # conventionally reported.
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
