Wu_2024_gfr_maturation <- function() {
  description <- paste(
    "Glomerular filtration rate (GFR) maturation model for preterm",
    "and term-born individuals from birth to 18 years of age",
    "(Wu 2024 Pharm Res). The model simultaneously characterises",
    "inulin clearance values (assumed equal to GFR; Eq. 1) and",
    "serum creatinine concentrations (assumed at quasi-steady state",
    "with synthesis = production / GFR; Eq. 2). GFR at birth",
    "GFRbirth is a linear function of birthweight and postnatal",
    "maturation follows a sigmoidal Emax (Hill) function of postnatal",
    "age (PNA, days; Eq. 3 and Eq. 4) with GFRmax allometrically",
    "scaled to current weight (exponent fCW) and PNA50 power-scaled",
    "to gestational age (exponent GAPNA50) so the rate of postnatal",
    "maturation depends on the degree of prematurity. Creatinine",
    "synthesis rate uses the Pierce 2021 (Kidney Int) age- and",
    "sex-dependent k with height and BSA scaling (Table S1 of the",
    "supplement). Time axis is PNA in days; the model has no drug",
    "compartments and no dosing -- it predicts GFR (mL/min, observed",
    "as inulin clearance) and Scr (mg/dL) at each user-supplied",
    "observation time given the subject covariates.",
    sep = " "
  )
  reference <- paste(
    "Wu Y, Allegaert K, Flint RB, Goulooze SC, Valitalo PAJ,",
    "de Hoog M, Mulla H, Sherwin CMT, Simons SHP, Krekels EHJ,",
    "Knibbe CAJ, Voller S.",
    "When will the Glomerular Filtration Rate in Former Preterm",
    "Neonates Catch up with Their Term Peers?",
    "Pharm Res. 2024;41(4):637-649.",
    "doi:10.1007/s11095-024-03677-3.",
    sep = " "
  )
  vignette <- "Wu_2024_gfr_maturation"
  units <- list(
    time          = "day (postnatal age, PNA)",
    dosing        = "n/a (no exogenous dosing; endogenous GFR / Scr model)",
    concentration = "mL/min (inulin CL = GFR); mg/dL (serum creatinine Scr)"
  )

  covariateData <- list(
    WT = list(
      description        = "Current body weight (CW)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Drives the allometric scaling of GFRmax with",
        "reference 1.75 kg (= 1,750 g) per Wu 2024 Table II. The paper",
        "expresses the reference as '1,750 g'; this model uses canonical",
        "WT in kg with the equivalent reference 1.75 kg so the form is",
        "(WT / 1.75)^fCW. Cohort range across the inulin / Scr datasets:",
        "0.4 kg to 83.3 kg (Wu 2024 Table I)."
      ),
      source_name        = "CW"
    ),
    WT_BIRTH = list(
      description        = "Birth weight (Bwb), time-fixed per subject",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at birth. Drives a linear effect on GFRbirth with",
        "reference 1.75 kg (= 1,750 g) per Wu 2024 Eq. 3 and Table II.",
        "Cohort range: 0.43 kg to 5.24 kg (Wu 2024 Table I)."
      ),
      source_name        = "Bwb"
    ),
    GA = list(
      description        = "Gestational age at birth",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at birth. Drives a power effect on PNA50",
        "(PNA50 = TVPNA50 * (GA/34)^GAPNA50) with reference 34 weeks",
        "per Wu 2024 Table II. Higher GA accelerates maturation",
        "(GAPNA50 negative). Cohort range: 23 to 44 weeks (Table I)."
      ),
      source_name        = "GA"
    ),
    HT = list(
      description        = "Subject height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Enters the Pierce 2021 creatinine synthesis",
        "rate equation linearly: syn_rate ~ k * HT / 88.4 * (BSA/1.73)",
        "(Wu 2024 supplement Table S1). When measured height was",
        "missing, the source data assembly used interpolation,",
        "extrapolation or growth-chart imputation (Methods, Imputation",
        "of Missing Demographic Information)."
      ),
      source_name        = "Height"
    ),
    BSA = list(
      description        = "Body surface area, Haycock-equation",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Computed from WT and HT via the Haycock",
        "equation in the source data assembly (Wu 2024 Methods). Enters",
        "the Pierce 2021 creatinine synthesis rate via (BSA/1.73)",
        "(reference 1.73 m^2 = the standardised adult body surface)."
      ),
      source_name        = "BSA"
    ),
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Time-fixed. Selects between the male and female branches of",
        "the Pierce 2021 k(age) constant in the creatinine synthesis",
        "rate (Wu 2024 supplement Table S1). Female:male ratio in the",
        "Wu 2024 cohort: inulin CL 229/154 (60% female), Scr datasets",
        "201/193 (51% female) (Table I).",
        "Missing sex (n=7 in the inulin CL dataset) was imputed by",
        "random binomial draw with probability 0.5 (Methods)."
      ),
      source_name        = "Sex"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 383L + 71L + 98L + 125L + 100L,
    n_studies      = 5L,
    age_range      = "PNA 0 to 6,570 days (0 days to 18 years) across the combined inulin CL and Scr datasets",
    weight_range   = "Current weight 0.4 to 83.3 kg; birthweight 430 to 5,000 g",
    gestational_age_range = "23 to 44 weeks (Wu 2024 Table I)",
    sex_female_pct = round(100 * (229 + 36 + 46 + 59 + 60) / (383 + 71 + 98 + 125 + 100), 1),
    race_ethnicity = "Not reported (pooled European / North American cohorts: Belgium, Netherlands, Poland, USA)",
    disease_state  = paste(
      "Preterm and term-born neonates, infants, children and",
      "adolescents. The inulin clearance subset (n=383) covers the",
      "infancy-to-adolescence range with the bulk of measurements in",
      "the first 90 days of life; the four serum creatinine datasets",
      "(n=394) span birth to 17.7 years and include intensive-care",
      "neonates (Erasmus MC Sophia, dataset 1; DINO study preterm",
      "<32 weeks GA, dataset 2) and adolescent ex-preterm survivors",
      "vs term controls at 11 years (PREMATCH-Leuven dataset 3;",
      "Krakow / Polish-American Children's Hospital dataset 4).",
      "Subjects undergoing acute kidney injury (AKI: Scr increase",
      ">= 0.3 mg/dL or >= 150 % over the previous observation; KDIGO",
      "definition) were excluded for the duration of the AKI episode",
      "(Methods, Exclusion of Scr During AKI)."
    ),
    dose_range     = "n/a (no exogenous drug PK modelled)",
    regions        = "Europe (Netherlands, Belgium, Poland) and United States (Dayton OH); pooled literature inulin CL data (Methods reference 12 and 39-44)",
    samples_inulin_cl = "431 inulin clearance values from 383 subjects",
    samples_scr    = "2,181 serum creatinine concentrations from 394 subjects (1,391 + 565 + 125 + 100 across datasets 1-4 per Table I)",
    notes          = paste(
      "Demographics from Wu 2024 Table I. The pooled inulin CL data",
      "include 333 subjects below 90 days PNA from the earlier Wu et",
      "al. analysis (paper reference 12) plus 50 additional subjects",
      "above 90 days PNA from a January 2023 literature search",
      "(references 39 to 44). Five clinical / literature sources",
      "contributed: literature-pooled inulin CL (1), Sophia ICU Scr",
      "(2), DINO study Scr (3), PREMATCH Leuven Scr (4) and Krakow",
      "outpatient Scr (5).",
      "Missing demographics in the inulin CL dataset (GA, birthweight,",
      "current weight, height, sex) were imputed per Methods 'Imputation",
      "of Missing Demographic Information' using the Fenton growth chart,",
      "the Haycock BSA equation, growth-curve interpolation, and a",
      "random binomial sex draw."
    )
  )

  ini({
    # ---- Structural fixed effects (Wu 2024 Table II) ----
    # The final-model equation in Table II is, with paper symbols on
    # natural scale and Bwb / CW in grams, PNA in days, GA in weeks:
    #   GFR = TVGFRbirth * (Bwb/1750)
    #       + ( TVGFRmax * (CW/1750)^fCW - TVGFRbirth * (Bwb/1750) )
    #         * PNA^gamma / ( (TVPNA50 * (GA/34)^GAPNA50)^gamma + PNA^gamma )
    # The model below uses canonical WT and WT_BIRTH in kg with the
    # equivalent reference 1.75 kg (= 1750 g), so the form is identical
    # numerically.
    lgfrbirth <- log(1.26)
    label("Typical GFR at birth TVGFRbirth at WT_BIRTH = 1.75 kg reference (mL/min)")  # Wu 2024 Table II: TVGFRbirth = 1.26 mL/min (RSE 4 %)
    lgfrmax   <- log(8.98)
    label("Typical GFR maximum TVGFRmax at WT = 1.75 kg reference (mL/min)")  # Wu 2024 Table II: TVGFRmax = 8.98 mL/min (RSE 12 %)
    lpna50    <- log(34)
    label("Typical PNA at half-max GFR maturation TVPNA50 at GA = 34 weeks reference (days)")  # Wu 2024 Table II: TVPNA50 = 34 days (RSE 31 %)
    lhill     <- log(1.03)
    label("Hill coefficient gamma on PNA in the GFR maturation function (unitless)")  # Wu 2024 Table II: gamma = 1.03 (RSE 7 %)

    # ---- Covariate-effect parameters (Wu 2024 Table II) ----
    e_wt_gfrmax <- 0.738
    label("Allometric exponent of current weight (WT / 1.75 kg) on GFRmax (unitless)")  # Wu 2024 Table II: fCW = 0.738 (RSE 5 %)
    e_ga_pna50  <- -3.61
    label("Power exponent of (GA / 34 weeks) on PNA50 (unitless; negative = higher GA accelerates maturation)")  # Wu 2024 Table II: GAPNA50 = -3.61 (RSE 18 %)

    # ---- Inter-individual variability (Wu 2024 Table II) ----
    # CV% to omega^2 conversion: omega^2 = log(1 + CV^2).
    etalgfrbirth ~ 0.09229  # CV 31.1 % -> omega^2 = log(1 + 0.311^2) = 0.09229 (Wu 2024 Table II, RSE 16 %, shrinkage 57 %)
    etalgfrmax   ~ 0.02850  # CV 17.0 % -> omega^2 = log(1 + 0.17^2)  = 0.02850 (Wu 2024 Table II, RSE 10 %, shrinkage 47 %)
    etalpna50    ~ 0.27184  # CV 55.9 % -> omega^2 = log(1 + 0.559^2) = 0.27184 (Wu 2024 Table II, RSE 43 %, shrinkage 57 %)
    etalhill     ~ 0.11862  # CV 35.5 % -> omega^2 = log(1 + 0.355^2) = 0.11862 (Wu 2024 Table II, RSE 20 %, shrinkage 55 %)

    # ---- Residual error (Wu 2024 Table II 'Residual errors' rows) ----
    # Inulin CL observations were log-transformed with constant
    # additive residual error in NONMEM (Methods, Stochastic Model).
    # That is, log(Cinulin_obs) = log(GFR) + eps with eps ~ N(0, sigma^2);
    # equivalently Cinulin ~ lnorm(SD) on the linear scale. The reported
    # 0.033 is the variance on the log scale; SD = sqrt(0.033).
    expSd_Cinulin <- sqrt(0.033)
    label("Log-scale residual SD for inulin clearance observations (unitless; lognormal form)")  # Wu 2024 Table II: sigma^2 inulin CL-log additive = 0.033 (RSE 55 %)
    # Scr observations carry a combined additive + proportional residual.
    # Reported values are variances; nlmixr2 prop()/add() expect SDs.
    addSd_Scr <- sqrt(0.00187)
    label("Additive residual SD on serum creatinine (mg/dL)")  # Wu 2024 Table II: sigma^2 Scr-additive = 0.00187 (RSE 36 %)
    propSd_Scr <- sqrt(0.0198)
    label("Proportional residual SD on serum creatinine (fraction)")  # Wu 2024 Table II: sigma^2 Scr-proportional = 0.0198 (RSE 40 %)
  })

  model({
    # ---- Time axis and unit harmonisation ----
    # The model treats the rxode2 time variable as PNA in days,
    # matching units$time. AGE in years is derived from PNA for use
    # in the Pierce 2021 creatinine-synthesis function.
    pna_d  <- time
    age_yr <- pna_d / 365.25

    # ---- Individual GFR-maturation parameters ----
    # GFRbirth = TVGFRbirth * (WT_BIRTH / 1.75 kg) per Eq. 3.
    # 1.75 kg matches the paper's 1,750 g reference; the WT_BIRTH
    # canonical column is in kg so (WT_BIRTH / 1.75) reproduces (Bwb_g / 1750).
    gfrbirth <- exp(lgfrbirth + etalgfrbirth) * (WT_BIRTH / 1.75)
    # GFRmax = TVGFRmax * (WT / 1.75 kg)^fCW per Table II.
    gfrmax   <- exp(lgfrmax + etalgfrmax) * (WT / 1.75)^e_wt_gfrmax
    # PNA50 = TVPNA50 * (GA / 34 weeks)^GAPNA50 per Table II.
    pna50_d  <- exp(lpna50 + etalpna50) * (GA / 34)^e_ga_pna50
    hill     <- exp(lhill + etalhill)

    # ---- GFR maturation function (Eq. 4 with covariate substitutions) ----
    GFR <- gfrbirth +
      (gfrmax - gfrbirth) * pna_d^hill /
        (pna50_d^hill + pna_d^hill)

    # ---- Pierce 2021 creatinine synthesis rate (Wu 2024 supplement Table S1) ----
    # k(age, sex) constants (Pierce et al. 2021 Kidney Int 99(4):948-956):
    #   Males:   k = 39.0 * 1.008^(age-12) for 1 to <12 y
    #            k = 39.0 * 1.045^(age-12) for 12 to <18 y
    #            k = 50.8                  for 18 to 25 y
    #   Females: k = 36.1 * 1.008^(age-12) for 1 to <12 y
    #            k = 36.1 * 1.023^(age-12) for 12 to <18 y
    #            k = 41.4                  for 18 to 25 y
    # Wu 2024 applies Pierce across the full neonatal-to-adolescent
    # range, extrapolating the 1 to <12 y formula below 1 year of age
    # (Methods: Pierce 2021 selected because it improved 12-18 y Scr
    # prediction relative to Schwartz; cf. Smeets 2018 cross-validation
    # cited in the Discussion).
    k_m_lt12 <- 39.0 * 1.008^(age_yr - 12)
    k_m_lt18 <- 39.0 * 1.045^(age_yr - 12)
    k_m_ge18 <- 50.8
    k_male   <- (age_yr < 12) * k_m_lt12 +
                (age_yr >= 12) * (age_yr < 18) * k_m_lt18 +
                (age_yr >= 18) * k_m_ge18
    k_f_lt12 <- 36.1 * 1.008^(age_yr - 12)
    k_f_lt18 <- 36.1 * 1.023^(age_yr - 12)
    k_f_ge18 <- 41.4
    k_female <- (age_yr < 12) * k_f_lt12 +
                (age_yr >= 12) * (age_yr < 18) * k_f_lt18 +
                (age_yr >= 18) * k_f_ge18
    k_pierce <- SEXF * k_female + (1 - SEXF) * k_male

    # Synthesis rate (mg/min * 100) per Pierce 2021 / Wu 2024 Table S1.
    # The 88.4 and 1.73 constants are the Pierce normalisation
    # references (88.4 cm height and 1.73 m^2 BSA).
    syn_rate <- k_pierce * HT / 88.4 * (BSA / 1.73)

    # ---- Observations ----
    # Inulin CL observation (mL/min) is GFR by paper assumption (Eq. 1).
    Cinulin <- GFR
    # Steady-state creatinine: Scr (mg/dL) = synth_rate / GFR with
    # synth_rate in (mg/min * 100) and GFR in mL/min; the factor-of-100
    # cancels with the mL -> dL conversion (1 dL = 100 mL).
    Scr <- syn_rate / GFR

    Cinulin ~ lnorm(expSd_Cinulin)
    Scr     ~ add(addSd_Scr) + prop(propSd_Scr)
  })
}
