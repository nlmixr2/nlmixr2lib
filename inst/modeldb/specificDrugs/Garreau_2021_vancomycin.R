Garreau_2021_vancomycin <- function() {
  description <- "Two-compartment IV population PK model for vancomycin given by continuous infusion to 78 critically ill adults in intensive care, 22 of them on continuous renal replacement therapy (Garreau 2021). Clearance depends on Cockcroft-Gault creatinine clearance (ideal-body-weight based) in patients off CRRT and on the CRRT effluent flow rate in patients on CRRT; the central volume depends on ideal body weight. As printed in the paper's final-model equations, each covariate enters as the exponential of a centred power term, CL = CLpop * exp((CRCL/41.4)^0.5), rather than as the plain power form of the Methods. Proportional residual error."
  reference <- "Garreau R, Falquet B, Mioux L, Bourguignon L, Ferry T, Tod M, Wallet F, Friggeri A, Richard JC, Goutelle S. Population Pharmacokinetics and Dosing Simulation of Vancomycin Administered by Continuous Injection in Critically Ill Patient. Antibiotics (Basel). 2021;10(10):1228. doi:10.3390/antibiotics10101228"
  vignette <- "Garreau_2021_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Cockcroft-Gault creatinine clearance computed with ideal body weight (raw, NOT BSA-normalized); 0 in patients on CRRT or anuric",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Garreau 2021 Sect. 4.1: CRCL was estimated by the Cockcroft-Gault equation with total, ideal and adjusted body weight (CG_TBW, CG_IBW, CG_AjBW); Sect. 3 states the retained descriptor is 'CRCL based on IBW' (CG_IBW). No BSA normalization. 'If the patient had anuria or was undergoing CRRT, the value of CRCL was set to 0 mL/min.' Reference 41.4 mL/min is the centring value of the Table 2 final-model equation (Eq. 1 defines the centring value as the population median). Learning-set Table 1 mean 50.4 +/- 29 mL/min. Used only when RRT_CRRT_STATUS = 0.",
      source_name = "CRCL (CG_IBW)"
    ),
    IBW = list(
      description = "Ideal body weight (Devine formula)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Garreau 2021 Sect. 4.1: IBW derived with the Devine formula (men 50 + 2.3 * (height_in - 60); women 45.5 + 2.3 * (height_in - 60)). Reference 64.1 kg is the centring value of the Table 2 final-model equation for V1. Learning-set Table 1 mean 63.5 +/- 9 kg.",
      source_name = "IBW"
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous renal replacement therapy indicator (1 = on CRRT, 0 = not)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not on CRRT)",
      notes = "Garreau 2021 Table 2 final-model equation switches the CL covariate model on 'CRRT = 0' / 'CRRT = 1'. Sect. 4.1: all dialysed patients were on continuous RRT (no intermittent haemodialysis). 22/78 (28.2%) of the learning set were on CRRT. The paper does not say whether the flag varied within a patient; it is encoded as a subject-level status and can be supplied time-varying in the data if needed.",
      source_name = "CRRT"
    ),
    RRT_CRRT_EFFLUENT_FLOW = list(
      description = "CRRT effluent flow rate",
      units = "mL/h",
      type = "continuous",
      reference_category = NULL,
      notes = "Garreau 2021 reports CRRT_EFR in mL/min (Table 1: learning set 35.6 +/- 18.7 mL/min) and centres it at 20 mL/min in the Table 2 final-model equation. The canonical column is in mL/h, so the model divides by 1200 mL/h (= 20 mL/min * 60). Multiply a mL/min value by 60 on ingestion. 'Calculated according to the CRRT technique' (Sect. 4.1); the paper does not print the formula. Used only when RRT_CRRT_STATUS = 1; its value is ignored otherwise.",
      source_name = "CRRTEFR"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Garreau 2021 Sect. 4.1 collected daily TBW and Sect. 3 reports that IBW was a better descriptor of V1 than TBW; not retained. Table 1 mean 77.9 +/- 20.4 kg."
    ),
    ABW = list(
      description = "Adjusted body weight",
      units = "kg",
      type = "continuous",
      notes = "Garreau 2021 Sect. 4.1 derived AjBW from TBW and screened it (and CG_AjBW); not retained."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Garreau 2021 Sect. 4.1 covariate list; not retained. Table 1 mean 68.9 +/- 12.3 years."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Garreau 2021 Sect. 4.1 covariate list; not retained. Table 1: 57 men / 21 women."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Garreau 2021 Sect. 4.1 covariate list; not retained. Table 1 mean 27.6 +/- 7.7 kg/m^2."
    ),
    BODYTEMP = list(
      description = "Body temperature",
      units = "degC",
      type = "continuous",
      notes = "Garreau 2021 Sect. 4.1 covariate list; not retained. Table 1 mean 37.8 +/- 1.1 degC."
    ),
    SOFA = list(
      description = "Sequential organ failure assessment score",
      units = "points",
      type = "continuous",
      notes = "Garreau 2021 Sect. 4.1 covariate list; not retained. Table 1 mean 11 +/- 4."
    ),
    SAPS_II = list(
      description = "Simplified Acute Physiology Score II (French IGS-II)",
      units = "points",
      type = "continuous",
      notes = "Garreau 2021 Table 1 reports IGS-II 55.9 +/- 17.3; Sect. 4.1 states that all collected covariates were tested; not retained."
    ),
    DIS_SEPTIC_SHOCK = list(
      description = "Septic shock at admission (Sepsis-3 definition)",
      units = "(binary)",
      type = "binary",
      notes = "Garreau 2021 Sect. 4.1 covariate list; not retained. Table 1: 83.3% of the learning set."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Garreau 2021 Sect. 4.1 covariate list; enters the final model only through the Cockcroft-Gault CRCL. Table 1 mean 139.3 +/- 75.7 umol/L."
    ),
    TPRO = list(
      description = "Serum total protein",
      units = "g/L",
      type = "continuous",
      notes = "Garreau 2021 Sect. 4.1 covariate list; not retained. Table 1 mean 60 +/- 11.3 g/L."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Garreau 2021 Sect. 4.1 covariate list; not retained. Table 1 mean 23.1 +/- 5.9 g/L."
    ),
    EVLWI = list(
      description = "Extravascular lung water index (PiCCO)",
      units = "mL/kg",
      type = "continuous",
      notes = "Garreau 2021 Sect. 3: PiCCO variables were assessed as covariates but none showed a significant influence. Table 1 mean 9.2 +/- 5.5 mL/kg."
    ),
    PVPI = list(
      description = "Pulmonary vascular permeability index (PiCCO)",
      units = "(unitless)",
      type = "continuous",
      notes = "Garreau 2021 Sect. 3: PiCCO variables were assessed as covariates but none showed a significant influence. Table 1 mean 1.9 +/- 1."
    ),
    CARDIAC_INDEX = list(
      description = "Cardiac index (PiCCO)",
      units = "L/min/m^2",
      type = "continuous",
      notes = "Garreau 2021 Sect. 3: PiCCO variables were assessed as covariates but none showed a significant influence. Table 1 mean 3.0 +/- 1.5 (unit misprinted there as umol/L)."
    ),
    CONMED_INOTROPE = list(
      description = "Concomitant inotropic agent",
      units = "(binary)",
      type = "binary",
      notes = "Garreau 2021 Sect. 3: an effect of inotropic agents on vancomycin clearance was found during model building but not retained because of borderline significance (p = 0.067)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 78L,
    n_studies = 1L,
    age_median = "68.9 years (mean, SD 12.3)",
    weight_median = "77.9 kg (mean, SD 20.4); IBW 63.5 kg (mean, SD 9)",
    height_median = "168.7 cm (mean, SD 8.5)",
    sex_female_pct = 26.9,
    race_ethnicity = "Not reported (single-centre French cohort)",
    disease_state = "Critically ill adults in intensive care receiving vancomycin by continuous infusion with invasive PiCCO haemodynamic monitoring; septic shock in 83.3%; SOFA 11 +/- 4; IGS-II 55.9 +/- 17.3. Patients with myeloma, cystic fibrosis or burns over 20% of body surface were excluded.",
    dose_range = "Loading dose 22.7 +/- 7.5 mg/kg then continuous infusion 28.6 +/- 9.4 mg/kg/day (mean +/- SD)",
    regions = "France (Croix-Rousse Hospital, Lyon)",
    renal_function = "Cockcroft-Gault CRCL mean 50.4 +/- 29 mL/min; serum creatinine 139.3 +/- 75.7 umol/L; 22/78 (28.2%) on continuous renal replacement therapy (effluent flow 35.6 +/- 18.7 mL/min)",
    n_concentrations = 335L,
    notes = "Learning dataset collected December 2013 to April 2015; 335 concentrations from 78 patients by routine morning TDM sampling over 4.1 +/- 2 days. External validation dataset (Centre Hospitalier Lyon Sud, February 2019 to September 2020): 417 concentrations from 84 patients, 4 (4.8%) on CRRT. Assay: immunoturbidimetry on Abbott Architect C8000, LLOQ 1.1 mg/L, linear range 1.1-100 mg/L. Fit in Monolix 2020R1 with SAEM."
  )

  ini({
    # Structural parameters (Garreau 2021 Table 2, fixed effects). CLpop is
    # clearance with the covariate term removed: the final-model equation is
    # CL0 = CLpop * exp((CRCL/41.4)^0.5), so the clearance of a patient off
    # CRRT at CRCL = 41.4 mL/min is 0.79 * e = 2.15 L/h, and V1pop is scaled by
    # exp((IBW/64.1)^alpha), so V1 at IBW = 64.1 kg is 27.3 * e = 74.2 L.
    lcl <- log(0.79)
    label("Clearance with the covariate term removed (L/h)") # Garreau 2021 Table 2: CLpop = 0.79 L/h (RSE 12.7%)
    lvc <- log(27.3)
    label("Central volume with the IBW term removed (L)") # Garreau 2021 Table 2: V1pop = 27.3 L (RSE 45.1%)
    lq <- log(6.08)
    label("Intercompartmental clearance (L/h)") # Garreau 2021 Table 2: Qpop = 6.08 L/h (RSE 41.8%); Sect. 2.2 rounds it to 6.1
    lvp <- log(61.3)
    label("Peripheral volume (L)") # Garreau 2021 Table 2: V2pop = 61.3 L (RSE 9.7%); Sect. 2.2 prose prints 63.1, a digit transposition -- the table value is used

    # Covariate exponents. Each sits INSIDE an exponential per the Table 2
    # final-model equations (see model()). The CRCL and effluent-flow exponents
    # appear only in those equations and are absent from the Table 2 list of
    # estimated fixed effects, so they are encoded as fixed.
    e_crcl_cl <- fixed(0.5)
    label("Exponent on (CRCL/41.4) inside exp() for CL off CRRT (unitless)") # Garreau 2021 Table 2 final-model equation: CL0 = CLpop * e^((CRCL/41.4)^0.5) if CRRT = 0
    e_rrt_crrt_effluent_flow_cl <- fixed(0.69)
    label("Exponent on (CRRT effluent flow/20 mL/min) inside exp() for CL on CRRT (unitless)") # Garreau 2021 Table 2 final-model equation: CL0 = CLpop * e^((CRRTEFR/20)^0.69) if CRRT = 1
    e_ibw_vc <- 1.88
    label("Exponent on (IBW/64.1) inside exp() for V1 (unitless)") # Garreau 2021 Table 2: alpha = 1.88 (RSE 67.1%); final-model equation V1 = V1pop * e^((IBW/64.1)^alpha)

    # Inter-individual variability. Table 2 reports omega as the STANDARD
    # DEVIATION of a log-normal random effect ('Random effects (standard
    # deviation)'); the variances entered here are omega^2.
    etalcl ~ 0.5776 # Garreau 2021 Table 2: omega_CL = 0.76 SD (RSE 11.1%) -> 0.76^2
    etalvc ~ 0.3721 # Garreau 2021 Table 2: omega_V = 0.61 SD (RSE 51.1%) -> 0.61^2
    etalq ~ 0.2401 # Garreau 2021 Table 2: omega_Q = 0.49 SD (RSE 33.8%) -> 0.49^2
    etalvp ~ 0.2304 # Garreau 2021 Table 2: omega_V2 = 0.48 SD (RSE 16.5%) -> 0.48^2

    # Residual error: proportional (Sect. 2.2, 'a two-compartment model with a
    # proportional residual error'), Monolix y = f + b * f * e.
    propSd <- 0.13
    label("Proportional residual error (fraction)") # Garreau 2021 Table 2: b = 0.13 (RSE 6.3%)
  })
  model({
    # Clearance covariate model (Garreau 2021 Table 2 final-model equations):
    #   CL0 = CLpop * exp((CRCL / 41.4)^0.5)       if CRRT = 0
    #   CL0 = CLpop * exp((CRRTEFR / 20)^0.69)     if CRRT = 1   (CRRTEFR in mL/min)
    # RRT_CRRT_EFFLUENT_FLOW is in mL/h, so the 20 mL/min centring value is
    # 1200 mL/h. Both bases are >= 0 (CRCL is set to 0 in CRRT patients), so the
    # arithmetic switch below never evaluates a negative power.
    crcl_term <- exp((CRCL / 41.4)^e_crcl_cl)
    crrt_term <- exp((RRT_CRRT_EFFLUENT_FLOW / 1200)^e_rrt_crrt_effluent_flow_cl)
    cl_cov <- (1 - RRT_CRRT_STATUS) * crcl_term + RRT_CRRT_STATUS * crrt_term

    cl <- exp(lcl + etalcl) * cl_cov
    # V1 = V1pop * exp((IBW / 64.1)^alpha) * exp(eta_V1)
    vc <- exp(lvc + etalvc) * exp((IBW / 64.1)^e_ibw_vc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
