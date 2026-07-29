Hwang_2022_tremelimumab <- function() {
  description <- "Two-compartment population PK model for tremelimumab (anti-CTLA-4 IgG2 kappa) with regimen-dependent sigmoidal time-varying clearance in adults with advanced solid tumours, dosed as monotherapy or in combination with durvalumab (Hwang 2022)"
  reference <- "Hwang M, Chia YL, Zheng Y, Chen CC-K, He J, Song X, Zhou D, Goldberg SB, Siu LL, Planchard D, Peters S, Mann H, Krug L, Even C. Population pharmacokinetic modelling of tremelimumab in patients with advanced solid tumours and the impact of disease status on time-varying clearance. Br J Clin Pharmacol. 2023;89(5):1601-1616. doi:10.1111/bcp.15622"
  vignette <- "Hwang_2022_tremelimumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL and Vc with reference weight 71 kg (Hwang 2022 Table 3 covariate rows and equations on page 1609).",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL with reference 39 g/L (Hwang 2022 Table 3 covariate row and equation on page 1609). Hwang 2022 imputed physiologically infeasible baseline values to the population median (Table 2 footnote e). Hypoalbuminaemia clinically defined as serum albumin <25 g/L.",
      source_name        = "ALB"
    ),
    COMBO_DURVA = list(
      description        = "Indicator for tremelimumab co-administered with durvalumab",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tremelimumab monotherapy)",
      notes              = "Hwang 2022 NONMEM control-stream column COMB selects between monotherapy and combination-with-durvalumab values of the time-varying-CL Tmax and lambda parameters via IF(COMB.EQ.0) / IF(COMB.EQ.1) (supplement Supplementary Results, Final-model NONMEM control stream). Stored under the canonical COMBO_DURVA per inst/references/covariate-columns.md.",
      source_name        = "COMB"
    )
  )

  population <- list(
    n_subjects     = 956L,
    n_studies      = 5L,
    age_range      = "22-87 years",
    age_median     = "65 years",
    weight_range   = "35.5-149.0 kg",
    weight_median  = "71.6 kg",
    sex_female_pct = 36.8,
    race_ethnicity = c(White = 75.5, Asian = 19.3, Black = 1.4, Other = 3.9),
    disease_state  = "Advanced solid tumours (pleural or peritoneal malignant mesothelioma 32.6%, lung cancer 22.1%, urothelial bladder cancer 17.1%, other 12.1%, biliary tract carcinoma 6.2%, oesophagus carcinoma 5.3%, breast cancer 3.1%)",
    dose_range     = "1-10 mg/kg or 75-750 mg flat IV every 4 weeks (Q4W); monotherapy 10 mg/kg Q4W or 750 mg flat, combination with durvalumab 1 mg/kg, 3 mg/kg, or 75 mg flat Q4W",
    regions        = "Global; Japan 6.7%, Korea 7.7%, Other 85.6%",
    treatment_regimen = "Monotherapy 38.8%, in combination with durvalumab 61.2%",
    ecog_distribution = "ECOG 0 38.9%, ECOG 1 60.9%, ECOG 2 0.1%",
    notes          = "Baseline demographics per Hwang 2022 Table 2 (model development N = 958 patients; the modelling dataset after sample exclusions was 4,043 PK samples from 956 patients). Trials pooled: Study 02 (NCT01938612, biliary tract / oesophageal cancer / SCCHN), Study 06 (NCT02000947, NSCLC), Study 10 (NCT02261220, urothelial cancer), DETERMINE (NCT01843374, pleural / peritoneal malignant mesothelioma), and D4884C00001 (NCT02527434, urothelial / triple-negative breast / pancreatic ductal adenocarcinoma). Median tumour size 73 mm (range 10-668). Median baseline albumin 39 g/L (range 15-52); median baseline LDH 209.5 U/L; postbaseline ADA-positive 6.1% (combination 7.8%, monotherapy 2.8%). External validation (n = 554, 4 phase 2/3 studies: ARCTIC, EAGLE, CONDOR, MYSTIC) not used for parameter estimation."
  )

  ini({
    # Structural parameters - typical-value baseline (no regimen / covariate
    # effects yet). Reference covariates (Hwang 2022 page 1609 final-model
    # equations): WT = 71 kg, ALB = 39 g/L. CL is the baseline value at
    # time = 0 since the time-varying multiplier exp(EMPIR) -> exp(0) = 1
    # at t = 0 (NONMEM control-stream form: EMPIR = TMAX * TIME^LAM /
    # (TC50^LAM + TIME^LAM); EMPIR(0) = 0). The half-life implied by these
    # baseline values is ~18 days, matching Hwang 2022 (page 1609).
    lcl  <- log(0.276); label("Baseline clearance CL_BASE for the reference covariates (L/day)")    # Hwang 2022 Table 3: CL = 0.276 L/day
    lvc  <- log(3.64);  label("Central volume of distribution Vc (L)")                                # Hwang 2022 Table 3: Vc = 3.64 L
    lq   <- log(0.357); label("Intercompartmental clearance Q (L/day)")                               # Hwang 2022 Table 3: Q  = 0.357 L/day
    lvp  <- log(2.40);  label("Peripheral volume of distribution Vp (L)")                             # Hwang 2022 Table 3: Vp = 2.40 L

    # Covariate effects (Hwang 2022 Table 3 covariate rows and the final-model
    # CL and Vc equations on page 1609). Power-style on continuous covariates.
    e_alb_cl <- -0.996; label("Power exponent of ALB on CL (unitless)")                                # Hwang 2022 Table 3: covariate 1 ALB on CL = -0.996
    e_wt_vc  <-  0.606; label("Power exponent of WT on Vc (unitless)")                                 # Hwang 2022 Table 3: covariate 2 WT on Vc  =  0.606
    e_wt_cl  <-  0.638; label("Power exponent of WT on CL (unitless)")                                 # Hwang 2022 Table 3: covariate 3 WT on CL  =  0.638

    # Time-varying clearance: a sigmoidal-Hill multiplier on log-CL,
    # parameterized separately for monotherapy (mono) and combination-with-
    # durvalumab (combo). The control-stream form (supplement
    # Supplementary Results) is
    #   EMPIR = (Tmax + ETA_Tmax) * TIME^lambda / (TC50^lambda + TIME^lambda)
    #   CL    = CL_BASE * exp(EMPIR) * covariate effects * exp(eta_CL)
    # so EMPIR(0) = 0 and EMPIR(t -> Inf) = (Tmax + eta_Tmax). A common TC50
    # is shared between the two regimens, with regimen-specific Tmax and lambda.
    cl_tmax_mono   <-  0.151; label("Asymptotic log-change in CL on monotherapy (unitless)")          # Hwang 2022 Table 3: covariate 4 Tmax (monotherapy)            =  0.151
    cl_tc50        <- 95.1;   label("Time at which the change in CL is 50%% of the asymptote (days)") # Hwang 2022 Table 3: covariate 5 TC50 (common)                = 95.1
    cl_lambda_mono <- 14.5;   label("Sigmoidicity exponent of time on CL on monotherapy (unitless)")  # Hwang 2022 Table 3: covariate 6 Lambda (monotherapy)         = 14.5
    cl_tmax_combo  <- -0.187; label("Asymptotic log-change in CL on combination with durvalumab (unitless)") # Hwang 2022 Table 3: covariate 7 Tmax (combination therapy)   = -0.187
    cl_lambda_combo <- 3.20;  label("Sigmoidicity exponent of time on CL on combination with durvalumab (unitless)") # Hwang 2022 Table 3: covariate 8 Lambda (combination therapy) =  3.20

    # IIV (Hwang 2022 Table 3). Reported values are variances on the internal
    # NONMEM scale (omega^2); the published "CV%" column equals
    # sqrt(omega^2) * 100 (a common log-normal short-hand). The supplement
    # control stream declares ETA(1) on CL and ETA(2) on Vc as $OMEGA BLOCK(2)
    # (correlated), but Hwang 2022 Table 3 does not report a final off-diagonal
    # estimate. Diagonal IIV is used here and the omitted off-diagonal is
    # called out in the vignette's Assumptions and deviations.
    etalcl    ~ 0.113   # Hwang 2022 Table 3: omega^2_CL   = 0.113 (33.6%% CV)
    etalvc    ~ 0.0536  # Hwang 2022 Table 3: omega^2_Vc   = 0.0536 (23.2%% CV)

    # Additive (not log-normal) IIV on Tmax. NONMEM control stream
    # (supplement Supplementary Results) defines Tmax_i = THETA + ETA(3) for
    # both regimens — a single shared eta is added to whichever regimen-
    # specific Tmax theta is active. Reported as omega^2 = 0.151
    # (38.9%% sqrt-based "CV%").
    etacl_tmax ~ 0.151  # Hwang 2022 Table 3: omega^2_Tmax = 0.151 (38.9%% CV)

    # Residual error (Hwang 2022 Table 3). Combined proportional + additive
    # in linear concentration space. NONMEM control stream
    #   Y = THETA(5)*F*EPS(1) + F + THETA(6)*EPS(2)
    # with EPS variances fixed to 1, so the THETA values (0.306 and 0.119)
    # are the proportional and additive residual SDs directly.
    propSd <- 0.306;  label("Proportional residual error (fraction)")              # Hwang 2022 Table 3: proportional residual error = 0.306
    addSd  <- 0.119;  label("Additive residual error (ug/mL)")                     # Hwang 2022 Table 3: additive residual error    = 0.119 ug/mL
  })
  model({
    # 1. Time-varying-CL multiplier (Hwang 2022 page 1609 final equation +
    # supplement NONMEM EMPIR = TMAX * TIME^LAM / (TC50^LAM + TIME^LAM)).
    # COMBO_DURVA = 0 picks the monotherapy Tmax and lambda; COMBO_DURVA = 1
    # picks the combination-with-durvalumab values. The single additive eta
    # etacl_tmax is added to whichever Tmax is active.
    cl_tmax_active   <- cl_tmax_mono   * (1 - COMBO_DURVA) + cl_tmax_combo   * COMBO_DURVA
    cl_lambda_active <- cl_lambda_mono * (1 - COMBO_DURVA) + cl_lambda_combo * COMBO_DURVA
    cl_tmax_i        <- cl_tmax_active + etacl_tmax
    cl_tv_mult       <- cl_tmax_i * t^cl_lambda_active /
      (cl_tc50^cl_lambda_active + t^cl_lambda_active)

    # 2. Individual structural parameters (Hwang 2022 page 1609). cl_base is
    # the t = 0 individual clearance after covariate adjustment; cl is the
    # time-varying instantaneous clearance.
    cl_base <- exp(lcl + etalcl) *
      (ALB / 39)^e_alb_cl *
      (WT  / 71)^e_wt_cl
    cl <- cl_base * exp(cl_tv_mult)

    vc <- exp(lvc + etalvc) * (WT / 71)^e_wt_vc
    q  <- exp(lq)
    vp <- exp(lvp)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. ODE system
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # 5. Observation and combined error
    # Dose in mg and volumes in L give central/vc in mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
