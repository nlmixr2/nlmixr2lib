Han_2016_bevacizumab <- function() {
  description <- "Two-compartment population PK model for IV bevacizumab in adult cancer patients (Han 2016) with allometric body-weight scaling and covariate effects of sex, baseline albumin, baseline alkaline phosphatase, and concomitant interferon alpha on clearance."
  reference <- "Han K, Peyret T, Marchand M, Quartino A, Gosselin NH, Girish S, Allison DE, Jin J. Population pharmacokinetics of bevacizumab in cancer patients with external validation. Cancer Chemother Pharmacol. 2016;78(2):341-351. doi:10.1007/s00280-016-3079-6"
  vignette <- "Han_2016_bevacizumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power scaling on CL and Q (shared exponent) and on V1 and V2 (shared exponent), centered at 70 kg. Model-building cohort median 74.8 kg (range 38.6-195 kg; Han 2016 Table 2 row BWT).",
      source_name        = "BWT"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 = female (reference); CL and V1 are larger in males.",
      notes              = "Han 2016 reports the categorical effect as 'Male on CL = 1.14' and 'Male on V1 = 1.18' (Table 3). Inverted from source-paper notation (source 'Male = 1') to the canonical SEXF (1 = female) by computing the male indicator as (1 - SEXF) in model().",
      source_name        = "Male (= 1 - SEXF)"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on CL with paper-reported reference value 39 g/L (model-building median, Han 2016 Table 2; Fig. 2 typical-patient definition). For patients with missing ALB in the source fit, Han 2016 imputed ALB = 41.8 g/L (Table 3 'Missing ALBU on CL'); users who wish to mirror that imputation should pre-fill missing ALB with 41.8 g/L on the data table.",
      source_name        = "ALBU"
    ),
    ALP = list(
      description        = "Baseline serum alkaline phosphatase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on CL with paper-reported reference value 109 U/L (model-building median, Han 2016 Table 2; Fig. 2 typical-patient definition). For patients with missing ALP in the source fit, Han 2016 imputed ALP = 76.3 U/L (Table 3 'Missing BALP on CL'); users who wish to mirror that imputation should pre-fill missing ALP with 76.3 U/L on the data table.",
      source_name        = "BALP"
    ),
    CONMED_IFNALPHA = list(
      description        = "Concomitant interferon alpha treatment indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = no concomitant IFN-alpha",
      notes              = "Han 2016 Table 3 reports 'IFNa on CL = 0.844' as the multiplicative factor when IFN-alpha is coadministered (102 of 1792 patients in the model-building cohort, all in the renal cell carcinoma RCC study BO17705). New canonical registered alongside this extraction as a sibling of CONMED_IFNB1A (the existing CONMED_IFNB1A entry explicitly invites the CONMED_IFNALPHA name for the alpha-interferon species).",
      source_name        = "IFNa"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 1792L,
    n_studies        = 15L,
    age_range        = "20-88 years",
    age_median       = "59 years",
    weight_range     = "38.6-195 kg",
    weight_median    = "74.8 kg",
    sex_female_pct   = 47,
    race_ethnicity   = "Caucasian 51.8%, Asian 3.7%, Black 2.7%, Hispanic 0.9%, Other 2.8%, unknown 38.1% (race recorded for 1113 of 1792 patients; Han 2016 Table 2).",
    disease_state    = "Adult cancer patients: colon/colorectal, non-small cell lung, kidney, pancreatic, breast, prostate (HRPC), and brain (glioblastoma) cancers; Phases I-IV studies in early and metastatic disease.",
    dose_range       = "Bevacizumab 1-20 mg/kg intravenous infusion (30-90 min) given once every 1, 2, or 3 weeks, as a single agent, in combination with chemotherapy, or with interferon alpha.",
    regions          = "Multinational (15 Genentech / Roche studies, predominantly North America and Europe).",
    renal_function   = "Normal 49.7%, mild impairment 37.2%, moderate impairment 10.9% (Han 2016 Table 2 row 'Renal function').",
    co_medication    = "Single agent 5.8%, chemotherapy 88.5%, interferon alpha 5.7% (Han 2016 Table 2 row 'Concomitant treatment').",
    n_observations   = "8943 bevacizumab serum concentrations (model-building) plus 1670 concentrations from 146 Japanese patients in three additional studies (external validation; not used for fitting).",
    notes            = "Pooled model-building dataset comprises 15 Genentech / Roche studies (Han 2016 Table 1); external validation uses three Japanese studies (JO18157, JO19901, JO19907). Bevacizumab serum concentrations measured by ELISA with LLOQ 78 ng/mL; concentrations below LLOQ were omitted before model-building. <5% of samples were below LLOQ and all were pre-dose."
  )

  ini({
    # ---- Structural parameters (reference 70 kg adult, female, no IFN-alpha) ----
    # Han 2016 Table 3 final-model estimates. Units in the paper are mL/h and mL;
    # converted here to L/day for consistency with antibody dosing once every 1-3 weeks
    # (CL_mL_h * 24/1000 -> L/day; V_mL / 1000 -> L).
    lcl <- log(8.6 * 24 / 1000);  label("Clearance (CL, L/day)")                                # Table 3 row CL = 8.6 mL/h
    lvc <- log(2678 / 1000);      label("Central volume of distribution (V1, L)")               # Table 3 row V1 = 2678 mL
    lq  <- log(18.6 * 24 / 1000); label("Inter-compartmental clearance (Q, L/day)")             # Table 3 row Q  = 18.6 mL/h
    lvp <- log(2423 / 1000);      label("Peripheral volume of distribution (V2, L)")            # Table 3 row V2 = 2423 mL

    # ---- Allometric exponents (estimated, shared between CL+Q and V1+V2) ----
    e_wt_cl_q  <- 0.589; label("Shared allometric exponent on log(WT/70) for CL and Q (unitless)")  # Table 3 row 'BWT on CL and Q'
    e_wt_vc_vp <- 0.470; label("Shared allometric exponent on log(WT/70) for V1 and V2 (unitless)") # Table 3 row 'BWT on V1 and V2'

    # ---- Categorical covariate effects on log-scale ----
    # Han 2016 Methods page 343: "if categorical covariate is not 0, then Effect = exp(theta_eff)".
    # Table 3 reports the e^theta_eff value; we encode theta_eff = log(reported_value) so
    # the model() block can write "effect = exp(theta * indicator)".
    e_male_lcl <- log(1.14);  label("Log-effect on CL for male vs female (unitless; Table 3 e^theta = 1.14)")  # Table 3 row 'Male on CL' = 1.14
    e_male_lvc <- log(1.18);  label("Log-effect on V1 for male vs female (unitless; Table 3 e^theta = 1.18)")  # Table 3 row 'Male on V1' = 1.18
    e_ifna_lcl <- log(0.844); label("Log-effect on CL for concomitant interferon alpha (unitless; Table 3 e^theta = 0.844)")  # Table 3 row 'IFNa on CL' = 0.844

    # ---- Continuous covariate effects (power-form exponents on CL) ----
    # Reference values are the model-building cohort medians (Han 2016 Table 2; Fig. 2 caption
    # defines the typical patient as 70 kg, female, ALB = 39 g/L, ALP = 109 U/L, no IFN-alpha).
    e_alb_lcl <- -0.473; label("Power exponent of (ALB / 39) on CL (unitless)")  # Table 3 row 'ALBU on CL' = -0.473
    e_alp_lcl <-  0.312; label("Power exponent of (ALP / 109) on CL (unitless)") # Table 3 row 'BALP on CL' = 0.312

    # ---- Inter-individual variability (log-normal; omega^2 = log(CV^2 + 1)) ----
    # Han 2016 Table 3 reports IIV as %CV. Full block IIV on CL, V1, V2 was declared in
    # the Methods (page 344) but Table 3 reports only the diagonal CVs; the off-diagonals
    # are unreported in the article and supplement. Encoded here as independent etas with
    # the reported diagonal variances; deviation documented in the vignette Errata section.
    etalcl ~ 0.08184  # IIV CL  29.2% -> omega^2 = log(1 + 0.292^2)
    etalvc ~ 0.03293  # IIV V1  18.3% -> omega^2 = log(1 + 0.183^2)
    etalvp ~ 0.15819  # IIV V2  41.4% -> omega^2 = log(1 + 0.414^2)

    # ---- Residual error (combined proportional + additive on linear scale) ----
    propSd <- 0.218;   label("Proportional residual error (fraction)")              # Table 3 row 'Prop. error (%)' = 21.8
    addSd  <- 0.0553;  label("Additive residual error (mg/L = microgram/mL)")       # Table 3 row 'Add. error (microgram/mL)' = 0.0553
  })
  model({
    # ---- Indicator derivations from canonical covariate columns ----
    # SEXF = 1 for female (reference) and 0 for male; paper's 'Male' indicator is the complement.
    male <- 1 - SEXF

    # ---- Covariate multiplicative factors ----
    cl_cov_male <- exp(e_male_lcl * male)
    vc_cov_male <- exp(e_male_lvc * male)
    cl_cov_ifna <- exp(e_ifna_lcl * CONMED_IFNALPHA)
    cl_cov_alb  <- (ALB / 39)^e_alb_lcl
    cl_cov_alp  <- (ALP / 109)^e_alp_lcl

    # ---- Individual PK parameters ----
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q  * cl_cov_male * cl_cov_ifna * cl_cov_alb * cl_cov_alp
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp * vc_cov_male
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    # ---- Two-compartment ODE system (IV infusion to central) ----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # ---- Observation and error model ----
    # Dose in mg, V in L -> Cc in mg/L = microgram/mL (matches Han 2016 Table 3 additive-error units).
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
