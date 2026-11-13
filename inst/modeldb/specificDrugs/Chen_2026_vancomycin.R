Chen_2026_vancomycin <- function() {
  description <- "One-compartment IV population PK model for vancomycin in critically ill Chinese children (Chen 2026), the Chen 2023 paediatric model re-estimated on a pooled two-centre dataset (821 children, 1767 concentrations). Clearance carries fixed allometric body weight (exponent 0.75, reference 70 kg), a sigmoid postmenstrual-age maturation function (TM50 37.1 weeks), a near-proportional power effect of BUN-and-creatinine-based estimated GFR capped at 120 mL/min/1.73 m^2 (exponent 1.02, reference 109 mL/min/1.73 m^2), and a 21.7% reduction in patients after cardiothoracic surgery on vasoactive agents. Central volume scales linearly with body weight (reference 70 kg) and rises as serum albumin falls (power exponent 0.255 on the 36.2 g/L over albumin ratio). The paper used this model as the Bayesian prior that generated the training targets for a CatBoost pre-dose CL/Vd predictor; only the population PK model is packaged. Residual variability is combined exponential-plus-additive, encoded with the exponential arm replaced by its first-order-equivalent proportional arm."
  reference <- "Chen J, Dai L, Xu H, Yao Y, Huang X, Ma J, Liu X, Yao H, Li J, Wu D, Sun J, Yu Z, Li Y, Yang H, Bu S. Machine-learning prediction of pre-dose pharmacokinetics optimizes initial vancomycin dosing in critically ill children. Front Pharmacol. 2026;17:1839423. doi:10.3389/fphar.2026.1839423"
  vignette <- "Chen_2026_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Vancomycin was given by intermittent IV infusion over 1 or 3 h (Chen 2026
  # Methods, 'Data collection and preprocessing'), so the dose enters
  # `central` directly. The same paragraph states that SERUM concentrations
  # were analysed by HPLC.
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 70 kg. Chen 2026 Table 1 median 11.00 kg (IQR 6.00-20.00). The allometric exponents are fixed as in the Chen 2023 framework that Chen 2026 retained: 0.75 on CL and 1 on Vd, both printed as constants without an RSE in the Results equations and in the Supplementary Table S4 formula notes.",
      source_name = "WT"
    ),
    PAGE = list(
      description = "Postmenstrual age, computed as gestational age plus postnatal age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "WEEKS, not the register-default months: TM50 = 37.1 weeks (Supplementary Table S4) and the sigmoid maturation function is only meaningful on the week scale. Chen 2026 Table 1 median 131.14 weeks (IQR 61.43-374.43); 80 children (9.74%) had a PMA below 44 weeks.",
      source_name = "PMA"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate from the BUN-and-creatinine-based paediatric equation, BSA-normalized: eGFR (mL/min/1.73 m^2) = 40.7 * (height (m) / SCr (mg/dL))^0.64 * (30 / BUN (mg/dL))^0.202",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Chen 2026 Methods equation (the same equation as Chen 2023 Equation 11, restated in the Supplementary Table S4 formula notes). Supply the UNCAPPED value: the model applies the 120 mL/min/1.73 m^2 ceiling of the Chen 2023 framework internally. Evidence that the ceiling was retained in the update: Chen 2026 states that 'the original structural and covariate framework was retained', and the Table 1 eGFR summary (median 108.40, IQR 85.17-120.00) has its upper quartile exactly at 120.00, whereas the per-centre uncapped values in Supplementary Table S3 reach an upper quartile of 140.5. Reference value 109 mL/min/1.73 m^2 as printed in the Results equation. Unit conversions from the Table 1 reporting units: height cm -> m, SCr umol/L -> mg/dL (divide by 88.4), BUN mmol/L -> mg/dL (multiply by 2.8).",
      source_name = "eGFR"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 36.2 g/L as printed in the Chen 2026 Results Vd equation. Albumin is in the DENOMINATOR, (36.2 / ALB)^0.255, so hypoalbuminaemia raises Vd. Chen 2026 Table 1 median 36.30 g/L (IQR 31.60-40.50).",
      source_name = "ALB"
    ),
    SURG_CARD_VASOACTIVE = list(
      description = "Cardiothoracic surgery with concomitant vasoactive support",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no cardiothoracic surgery, or cardiothoracic surgery without concomitant vasoactive agents)",
      notes = "Chen 2026 Methods, 'Data collection and preprocessing': 'We defined CTS patients as those who underwent CTS and received concomitant vasoactive agents (including milrinone, dopamine, epinephrine or norepinephrine) during vancomycin therapy ... and CTS status was handled as a time-varying covariate in the PPK model.' Supply it per record. Table 1: 86 of 821 children (10.48%).",
      source_name = "CTS"
    )
  )

  # Recorded in Chen 2026 Table 1 but not part of the population PK model.
  # Sex, ICU admission and the PMA class entered only the CatBoost predictor.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "One of the six CatBoost predictors, not a covariate of the population PK model ('no new covariates ... materially improved the model during updating'). Table 1: 343 of 821 female (41.78%).",
      source_name = "Sex"
    ),
    ICU = list(
      description = "Intensive care unit admission",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not in ICU)",
      notes = "CatBoost predictor only. Table 1: 675 of 821 (82.22%).",
      source_name = "ICU"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Input (converted to metres) to the eGFR equation carried by CRCL. Table 1 median 81.00 cm (IQR 62.00-115.00).",
      source_name = "Height"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Input to the eGFR equation (divide by 88.4 for mg/dL). Table 1 median 23.80 umol/L (IQR 17.60-33.00).",
      source_name = "Scr"
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Input to the eGFR equation (multiply by 2.8 for mg/dL). Table 1 median 3.60 mmol/L (IQR 2.45-5.00).",
      source_name = "BUN"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 821L,
    n_studies = 1L,
    n_sites = 2L,
    n_concentrations = 1767L,
    age_range = "under 18 years (eligibility criterion)",
    age_median = "1.75 years postnatal (IQR 0.41-6.41); postmenstrual age 131.14 weeks (IQR 61.43-374.43)",
    weight_median = "11.00 kg (IQR 6.00-20.00)",
    sex_female_pct = 41.78,
    race_ethnicity = c(Chinese = 100),
    disease_state = "Critically ill children (suspected or proven infection plus pSOFA >= 2) who received intravenous vancomycin for at least 3 consecutive days and had at least one serum concentration. 82.22% were in an ICU and 10.48% were after cardiothoracic surgery on vasoactive agents. Excluded: renal replacement therapy and extracorporeal membrane oxygenation.",
    dose_range = "Intermittent IV infusion over 1 or 3 h (dose amounts not tabulated in Chen 2026)",
    regions = "China (Xinhua Hospital, Shanghai Jiao Tong University School of Medicine, 753 children; People's Hospital of the Inner Mongolia Autonomous Region, 68 children)",
    renal_function = "eGFR median 108.40 mL/min/1.73 m^2 (IQR 85.17-120.00, after the 120 ceiling)",
    notes = "Retrospective therapeutic-drug-monitoring data, 2013-06 to 2023-03. 1767 concentrations (840 peaks, 927 troughs); 89.2% of children contributed at least two. The Chen 2023 one-compartment model was re-estimated in NONMEM 7.6 with FOCE-I on the pooled data with the structure and covariates unchanged (Supplementary Table S4 lists original and updated estimates side by side). A 53-child external cohort from Shanghai Children's Medical Center was used only to validate the CatBoost predictor."
  )

  ini({
    # Structural parameters (Chen 2026 Supplementary Table S4, 'Updated model'
    # column; the same values are printed in the Results equations). The
    # reference subject weighs 70 kg, is fully mature, has an eGFR of 109
    # mL/min/1.73 m^2 and an albumin of 36.2 g/L, and is not after
    # cardiothoracic surgery.
    lcl <- log(7.84); label("Clearance standardized to WT=70 kg, mature, eGFR=109 mL/min/1.73 m^2, non-CTS (L/h)")  # Chen 2026 Table S4: theta_CL = 7.84 L/h (RSE 1.9%)
    lvc <- log(36.6); label("Volume of distribution standardized to WT=70 kg and ALB=36.2 g/L (L)")                 # Chen 2026 Table S4: theta_V = 36.6 L (RSE 1.6%)

    # Allometric exponents, printed as constants without uncertainty in the
    # Results equations and the Table S4 formula notes.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on (WT/70) for CL (unitless)")  # Chen 2026 Results CL equation: (WT/70)^0.75
    e_wt_vc <- fixed(1);    label("Allometric exponent on (WT/70) for Vd (unitless)")  # Chen 2026 Results Vd equation: linear (WT/70)

    # Sigmoid postmenstrual-age maturation on CL, MF = 1 / (1 + (PMA/TM50)^Hill).
    pma_tm50 <- 37.1;  label("Postmenstrual age at 50% mature CL (weeks)")                       # Chen 2026 Table S4: TM50 = 37.1 (RSE 3.8%)
    pma_hill <- -1.64; label("Hill coefficient of the sigmoid CL maturation function (unitless)")  # Chen 2026 Table S4: Hill = -1.64 (RSE 9.1%)

    # Covariate effects on CL.
    e_crcl_cl <- 1.02;  label("Power exponent on (capped eGFR / 109) for CL (unitless)")                                # Chen 2026 Table S4: theta_eGFR = 1.02 (RSE 3.3%)
    e_cts_cl  <- 0.783; label("Multiplicative CL factor for cardiothoracic surgery with vasoactive support (unitless)")  # Chen 2026 Table S4: theta_1 = 0.783 (RSE 3.2%)

    # Covariate effect on Vd.
    e_alb_vc <- 0.255; label("Power exponent on (36.2 / ALB) for Vd (unitless)")  # Chen 2026 Table S4: theta_ALB = 0.255 (RSE 25.9%)

    # Correlated between-subject variability on CL and Vd. Table S4 reports the
    # diagonals as %CV and the off-diagonal as a covariance on the OMEGA scale,
    # so the diagonals are read on the same scale: omega^2 = (%CV / 100)^2,
    # following Chen_2023_vancomycin.R for the same table layout.
    #   CL: 0.283^2 = 0.080089     Vd: 0.272^2 = 0.073984
    #   implied correlation 0.0548 / sqrt(0.080089 * 0.073984) = 0.712
    etalcl + etalvc ~ c(0.080089,
                        0.054800, 0.073984)  # Chen 2026 Table S4: CL 28.3 %CV (RSE 5.8%), Cov(CL-Vd) 0.0548 (RSE 7.6%), Vd 27.2 %CV (RSE 16.6%)

    # Residual variability: combined exponential and additive (Table S4 rows
    # 'Exponential (%CV)' and 'Additive (mg/L)'). rxode2 cannot solve a
    # log-normal arm combined with an additive arm, so the exponential arm is
    # encoded as its first-order equivalent proportional arm,
    # F * exp(eps1) ~= F * (1 + eps1), as in Chen_2023_vancomycin.R.
    propSd <- 0.177; label("Proportional residual error (fraction)")  # Chen 2026 Table S4: Exponential 17.7 %CV (RSE 10.7%)
    addSd  <- 0.154; label("Additive residual SD (ug/mL)")            # Chen 2026 Table S4: Additive 0.154 mg/L (RSE 37.7%); 1 mg/L == 1 ug/mL
  })
  model({
    # Renal-function ceiling of the Chen 2023 covariate model, retained in the
    # update (see covariateData$CRCL$notes for the Table 1 evidence).
    crcl_capped <- min(CRCL, 120)

    # Sigmoid postmenstrual-age maturation; Hill is negative, so MF rises from
    # 0 towards 1 as PMA increases.
    maturation_cl <- 1 / (1 + (PAGE / pma_tm50)^pma_hill)

    # Chen 2026 Results equations:
    #   CL = 7.84 * (WT/70)^0.75 * MF * (eGFR/109)^1.02 * 0.783^CTS * exp(eta1)
    #   Vd = 36.6 * (WT/70) * (36.2/ALB)^0.255 * exp(eta2)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * maturation_cl *
      (crcl_capped / 109)^e_crcl_cl * e_cts_cl^SURG_CARD_VASOACTIVE
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc * (36.2 / ALB)^e_alb_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg and volume in L give mg/L, i.e. ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
