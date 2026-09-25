Zhou_2019_vancomycin <- function() {
  description <- "One-compartment IV population PK model for vancomycin in 70 Chinese geriatric patients (age >= 65 years) with hospital- or community-acquired pulmonary infection (Zhou 2019). Clearance scales by power exponent with raw Cockcroft-Gault creatinine clearance (mL/min, reference 56.28); volume of distribution carries no retained covariate. Estimated from 125 steady-state trough concentrations collected by routine therapeutic drug monitoring."
  reference <- "Zhou Y, Gao F, Chen C, Ma L, Yang T, Liu X, Liu Y, Wang X, Zhao X, Que C, Li S, Lv J, Cui Y, Yang L. Development of a Population Pharmacokinetic Model of Vancomycin and its Application in Chinese Geriatric Patients with Pulmonary Infections. Eur J Drug Metab Pharmacokinet. 2019;44(3):361-370. doi:10.1007/s13318-018-0534-2"
  vignette <- "Zhou_2019_vancomycin"
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
      notes = "Zhou 2019 Sect. 3.3: 'CLCR was calculated using the Cockcroft-Gault formula'; no BSA normalization is mentioned anywhere in the paper, so values are raw mL/min (per inst/references/covariate-columns.md, CRCL accepts raw Cockcroft-Gault when the source does not BSA-normalize). Table 1 cohort mean 56.3 mL/min (SD 22.1). The reference value 56.28 mL/min in the Eq. (9) power term is that same cohort MEAN carried to four significant figures (Table 1 prints it rounded to 56.3), not a median -- Sect. 4 describes normally distributed continuous variables as mean +/- SD. The paper stratifies its dosing recommendations at CLCR = 50 mL/min.",
      source_name = "CLCR"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Zhou 2019 Table 2 model 4: WT on CL alone lowered OFV by 8.559 (p < 0.05) and Fig. 1a shows WT positively related to CL, but in forward inclusion on top of CLCR (model 6) it gave only dOFV -3.781 (p > 0.05) and was not retained in the final model. Table 1 cohort mean 60.7 kg (SD 10.2)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Zhou 2019 Table 2 model 3: AGE on CL alone lowered OFV by 8.782 (p < 0.05) and Fig. 1b shows AGE inversely related to CL. It survived forward inclusion (model 5, dOFV -4.696) but was dropped in backward elimination (model 8, dOFV +4.696 < 6.64, p > 0.01). Table 1 cohort mean 78.3 years (SD 6.96)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Zhou 2019 Sect. 3.3 covariate list; Fig. 1c shows SCR inversely related to CL. Not retained -- its information enters the final model through the Cockcroft-Gault CLCR. Table 1 cohort mean 90.6 umol/L (SD 31)."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "categorical",
      notes = "Zhou 2019 Sect. 3.3 screened SEX via Eqs. (7)-(8); not retained in the final model. Table 1: 49 male / 21 female."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Zhou 2019 Sect. 3.3 covariate list; not retained. Table 1 reports 29.3 +/- 4.12 under the header 'ALB/g/dL', but 29.3 g/dL is physiologically impossible and the companion 'TP/g/dL' row reads 59.8 -- both are g/L mislabelled as g/dL."
    ),
    TP = list(
      description = "Total serum protein",
      units = "g/L",
      type = "continuous",
      notes = "Zhou 2019 Sect. 3.3 covariate list; not retained. Table 1 reports 59.8 +/- 8.92 under a 'g/dL' header that is g/L (see ALB note)."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mmol/L",
      type = "continuous",
      notes = "Zhou 2019 Sect. 3.3 covariate list; not retained. Table 1 median 10.5 mmol/L (range 3.18-86)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Zhou 2019 Sect. 3.3 covariate list; not retained. No cohort summary printed in Table 1."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Zhou 2019 Sect. 3.3 covariate list; not retained. No cohort summary printed in Table 1."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Zhou 2019 Sect. 3.3 covariate list (TBIL); not retained. No cohort summary printed in Table 1."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 70L,
    n_studies = 1L,
    age_range = ">= 65 years (inclusion criterion)",
    age_median = "78.3 years (mean, SD 6.96)",
    weight_median = "60.7 kg (mean, SD 10.2)",
    height_median = "161 cm (mean, SD 10)",
    sex_female_pct = 30,
    race_ethnicity = "Chinese (single-center study, Peking University First Hospital, Beijing)",
    disease_state = "Geriatric inpatients aged >= 65 years with pulmonary infection -- 57 (81.4%) hospital-acquired pneumonia, 13 (18.6%) community-acquired pneumonia. Comorbid respiratory failure in 46 (65.7%) and hypertension in 38 (54.3%). Patients with multiple organ failure, renal replacement therapy or low-volume shock were excluded.",
    dose_range = "500 mg IV (q6h, q8h, q12h, q24h, q48h) or 1000 mg IV (q8h, q12h), each given as a 1.5-2 h infusion; mean daily dose 1.55 g/day (SD 0.770)",
    regions = "China (Beijing)",
    renal_function = "Cockcroft-Gault CLCR mean 56.3 mL/min (SD 22.1); serum creatinine mean 90.6 umol/L (SD 31)",
    n_concentrations = 125L,
    notes = "Retrospective single-center therapeutic-drug-monitoring study, January 2012 to December 2016 (Sect. 3.1; the Abstract Methods says January 2011, an internal inconsistency in the paper). 125 observations from 70 patients, 1-5 per patient (mean 1.79); observed concentration mean 17 mg/L (SD 8.03). ALL observations are steady-state TROUGH samples drawn 0.5-2 h before the fourth or fifth dose, so no peak or distribution-phase data informed the fit -- the volume of distribution is therefore only weakly identified and its 154 L point estimate is high for vancomycin (about 2.5 L/kg at the cohort mean weight), implying a typical terminal half-life near 44 h. Assay: chemiluminescent microparticle immunoassay, ARCHITECT i1000, CV < 10%, linear range 3-100 mg/L. Fit in NONMEM 7.3.0 with FOCE-I; ADVAN1 TRANS2. Evaluated by goodness-of-fit plots, NPDE (mean 0.248, variance 1.28) and 1000-sample bootstrap (92.1% success rate; Table 4)."
  )

  ini({
    # Structural parameters (Zhou 2019 Table 3 'Final model' column). The
    # reference subject has CLCR = 56.28 mL/min, the cohort mean.
    lcl <- log(2.45)
    label("Clearance at CLCR=56.28 mL/min (L/h)") # Zhou 2019 Table 3: theta_1CL = 2.45 L/h (RSE 6.9%); bootstrap median 2.43, 95% CI (2.09, 2.81) in Table 4
    lvc <- log(154)
    label("Volume of distribution (L)") # Zhou 2019 Table 3: theta_2Vd = 154 L (RSE 9.2%); bootstrap median 154, 95% CI (117, 191) in Table 4

    # Covariate effect (Zhou 2019 Eq. 9):
    #   CL (L/h) = 2.45 * (CLCR / 56.28)^0.542
    e_crcl_cl <- 0.542
    label("Power exponent on (CRCL/56.28) for CL") # Zhou 2019 Table 3: theta_3CLCR on CL = 0.542 (RSE 35.1%); bootstrap median 0.538, 95% CI (0.206, 0.878) in Table 4

    # Inter-individual variability, exponential per Eqs. (1)-(2):
    #   CL_i = theta_CL * exp(eta_i), Vd_i = theta_Vd * exp(eta_i).
    # Table 3's omega entries are NONMEM $OMEGA variances of eta, entered here
    # directly (no CV back-transformation). Equivalent CVs are
    # sqrt(exp(0.174)-1) = 43.6% on CL and sqrt(exp(0.339)-1) = 63.4% on Vd.
    # The variance reading is confirmed by Table 5: the paper's own 1000-run
    # simulation prints between-subject CVs of 45.6-52.1% across the seven
    # regimens, which this model reproduces at 41-48% on the individual
    # prediction and 48-55% once the proportional residual term is added in
    # quadrature; reading the same numbers as SDs would give only about 22%.
    etalcl ~ 0.174 # Zhou 2019 Table 3, row 'omega CL' = 0.174 (RSE 21.2%); bootstrap median 0.162, 95% CI (0.092, 0.256) in Table 4
    etalvc ~ 0.339 # Zhou 2019 Table 3, row 'omega V' = 0.339 (RSE 37.8%); bootstrap median 0.289, 95% CI (0.121, 0.557) in Table 4

    # Residual error, combined per Eq. (3): C_obs = C_pre * (1 + eps1) + eps2.
    # Table 3 gives sigma_1 = 0.0657 as a NONMEM $SIGMA variance, so the
    # proportional SD is sqrt(0.0657) = 0.2563. sigma_2 was held at zero, so
    # the additive term is present in the paper's error model but contributes
    # nothing.
    propSd <- 0.2563
    label("Proportional residual error (fraction)") # Zhou 2019 Table 3: sigma_1 = 0.0657 variance (RSE 34.2%) -> sqrt(0.0657) = 0.2563
    addSd <- fixed(0)
    label("Additive residual error (mg/L)") # Zhou 2019 Table 3: sigma_2 = 0 FIX
  })
  model({
    # Individual PK parameters. CL scales with raw Cockcroft-Gault creatinine
    # clearance as a power function referenced to the cohort mean 56.28 mL/min
    # (Eq. 9); Vd has no retained covariate (Eq. 10).
    cl <- exp(lcl + etalcl) * (CRCL / 56.28)^e_crcl_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, volume in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
