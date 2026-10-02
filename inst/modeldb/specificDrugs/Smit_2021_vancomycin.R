Smit_2021_vancomycin <- function() {
  description <- "Two-compartment IV population PK model for vancomycin in normal-weight, overweight and obese children and adolescents aged 1-18 years with a wide range of renal function (Smit 2021). Clearance scales as a power function of total body weight (estimated exponent 0.745, reference 22.1 kg) and linearly with bedside-Schwartz creatinine clearance capped at 120 mL/min/1.73 m^2 (reference 100 mL/min/1.73 m^2); central and peripheral volumes scale linearly with total body weight and intercompartmental clearance as a power function of total body weight (exponent 0.599). Residual variability is additive on log-transformed concentrations."
  reference <- "Smit C, Goulooze SC, Bruggemann RJM, Sherwin CM, Knibbe CAJ. Dosing Recommendations for Vancomycin in Children and Adolescents with Varying Levels of Obesity and Renal Dysfunction: a Population Pharmacokinetic Study in 1892 Children Aged 1-18 Years. AAPS J. 2021;23(3):53. doi:10.1208/s12248-021-00577-x"
  vignette <- "Smit_2021_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Vancomycin is given as an IV infusion (60-min infusions per Smit 2021
  # Methods, Patients and Setting), so the dose enters `central` directly.
  # `central` is verified: the Abstract and Methods state vancomycin SERUM
  # concentrations were measured (Abbott Architect immunoassay). The paper
  # does not state what matrix the peripheral compartment represents.
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference weight 22.1 kg (Smit 2021 Table II footnote: TVCL/TVV1/TVQ/TVV2 are typical values for an individual weighing 22.1 kg; supplement control stream divides WT by 22.1). Table I medians 20.6 kg (normal weight), 25.0 kg (overweight) and 30.0 kg (obese); overall range 5.8-188 kg. Enters CL (power, estimated exponent), V1 and V2 (linear) and Q (power, estimated exponent).",
      source_name = "WT"
    ),
    CRCL = list(
      description = "Creatinine clearance from the revised bedside Schwartz equation, BSA-normalized: CLcr (mL/min/1.73 m^2) = 0.413 * length (cm) / serum creatinine (mg/dL)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying in the source analysis (serum creatinine within 168 h of a dose, next-observation-carried-backward within an individual; supplement Methods). CAPPED at 120 mL/min/1.73 m^2 inside the model (Table II footnote a; control stream IF(SCHW.GT.120) SCHW_MAX=120) -- supply the uncapped value, the model applies the cap. Enters CL linearly as (min(CRCL, 120) / 100) (control stream THETA(5) exponent 1 FIX). Main-text Eq. 1 prints the constant as 0.41; the supplement Results give 0.413. Table I medians 121.2 / 114.7 / 111.7 mL/min/1.73 m^2 by weight group; range 8.6-963.5.",
      source_name = "SCHW"
    )
  )

  # Screened in the Smit 2021 covariate analysis (supplement Methods and
  # Results) but not retained in the final model; documentation only.
  covariatesDataExcluded <- list(
    NEUTROPENIA = list(
      description = "Neutropenia indicator (absolute neutrophil count < 1.5e9 cells/L)",
      units = "(binary)",
      type = "binary",
      reference_category = "not neutropenic",
      notes = "Tested as a binary covariate on CL: OFV +0.8 versus the model with TBW and CLcr on CL (supplement Results); not retained.",
      source_name = "NPEN"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Serum creatinine (power) with TBW on CL gave OFV -2126.9 versus the structural model, worse than the -2356.1 of bedside Schwartz CLcr with TBW (supplement Results); not retained as a standalone covariate.",
      source_name = "CREAT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1892L,
    n_studies = 1L,
    n_sites = 21L,
    n_concentrations = 5524L,
    age_range = "1-18 years",
    age_median = "6.9 years (normal weight; IQR 2.9-13.2), 7.2 years (overweight), 6.9 years (obese)",
    weight_range = "5.8-188 kg",
    weight_median = "20.6 kg (normal weight), 25.0 kg (overweight), 30.0 kg (obese)",
    sex_female_pct = 43.7,
    race_ethnicity = c(Caucasian = 87.9, Asian = 0.7, Hispanic = 1.4, `African American` = 2.3, Other = 7.8),
    disease_state = "Hospitalized children and adolescents receiving intravenous vancomycin at 21 Intermountain Healthcare hospitals (Utah, USA); 1344 normal weight, 247 overweight and 301 obese (WHO/CDC BMI-for-age > 85th / > 95th percentile). Patients on renal replacement therapy or ECMO were excluded. 670 (35%) admitted to ICU; 17% neutropenic.",
    dose_range = "Clinician-chosen regimens, generally 15-20 mg/kg two to four times daily as 60-min IV infusions",
    regions = "United States (Utah)",
    renal_function = "Bedside Schwartz CLcr median ~115 mL/min/1.73 m^2, range 8.6-963.5; 12 patients below 30 mL/min/1.73 m^2",
    notes = "Retrospective multicenter TDM data collected 2006-2012 (Smit 2021 Methods; Table I). Female percentage is the size-weighted complement of Table I % male (57.3, 53.4, 54.1). Race percentages are pooled over the three weight groups of Table I. Log-transformed concentrations fitted in NONMEM 7.4 with FOCE-I; eta shrinkage 24% (CL) and 57% (V2), epsilon shrinkage 16%."
  )

  ini({
    # Structural parameters (Smit 2021 Table II, final model). Reference
    # individual: 22.1 kg total body weight, bedside Schwartz CLcr 100 mL/min/1.73 m^2.
    lcl <- log(2.12); label("Clearance at WT = 22.1 kg and CRCL = 100 mL/min/1.73 m^2 (L/h)") # Table II TVCL 2.12 L/h (RSE 1%)
    lvc <- log(8.90); label("Central volume at WT = 22.1 kg (L)") # Table II TVV1 8.90 L (RSE 3%)
    lq <- log(1.55); label("Intercompartmental clearance at WT = 22.1 kg (L/h)") # Table II TVQ 1.55 (RSE 5%)
    lvp <- log(12.3); label("Peripheral volume at WT = 22.1 kg (L)") # Table II TVV2 12.3 L (RSE 6%)

    # Covariate exponents
    e_wt_cl <- 0.745; label("Power exponent on (WT/22.1) for CL (unitless)") # Table II theta1 0.745 (RSE 2%)
    e_wt_q <- 0.599; label("Power exponent on (WT/22.1) for Q (unitless)") # Table II theta2 0.599 (RSE 9%)
    e_wt_vc_vp <- fixed(1); label("Exponent on (WT/22.1) for V1 and V2 (unitless)") # Table II TVV1/TVV2 x (TBW/22.1), linear; supplement control stream THETA(7) (1) FIX
    e_crcl_cl <- fixed(1); label("Exponent on (min(CRCL,120)/100) for CL (unitless)") # Table II TVCL x (SCHW/100), linear; supplement control stream THETA(5) (1) FIX

    # Between-subject variability (Table II, as CV%; footnote c: CV = sqrt(exp(omega^2) - 1),
    # so omega^2 = log(CV^2 + 1)). Covariance printed directly on the omega scale.
    # IIV on V1 and Q is 0 FIX in the supplement control stream and is omitted.
    etalcl + etalvp ~ c(
      0.079152,
      -0.085, 0.792993
    ) # Table II: CL 28.7% CV -> log(1 + 0.287^2); cov CL-V2 -0.085; V2 110% CV -> log(1 + 1.10^2)

    # Residual variability: NONMEM Y = log(F) + ERR(1) on log-transformed
    # concentrations (supplement control stream $ERROR). Table II prints 0.0789,
    # read as the $SIGMA variance (the control stream $SIGMA 0.0788 is labelled
    # 'PROP ERR IN LOGDOMAIN'), so the log-scale SD is sqrt(0.0789).
    expSd <- 0.28089; label("Log-scale residual SD (unitless)") # Table II proportional error 0.0789 (RSE 6%) = sigma^2; sqrt(0.0789) = 0.2809
  })
  model({
    # Bedside Schwartz CLcr capped at 120 mL/min/1.73 m^2 (Table II footnote a;
    # control stream IF(SCHW.GT.120) SCHW_MAX=120).
    crcl_cap <- min(CRCL, 120)

    # Individual PK parameters (Table II; supplement control stream $PK)
    cl <- exp(lcl + etalcl) * (WT / 22.1)^e_wt_cl * (crcl_cap / 100)^e_crcl_cl
    vc <- exp(lvc) * (WT / 22.1)^e_wt_vc_vp
    q <- exp(lq) * (WT / 22.1)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 22.1)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volume in L, so central / vc is mg/L.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
