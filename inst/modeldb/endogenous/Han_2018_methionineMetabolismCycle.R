Han_2018_methionineMetabolismCycle <- function() {
  description <- "Preclinical (rat). Seven-compartment mechanistic methionine metabolism cycle (MMC) model in Zucker Diabetic Fatty (ZDF) rats vs non-diabetic controls; predicts plasma methionine and homocysteine after IV methionine."
  reference   <- "Han N, Chae JW, Jeon J, Lee J, Back HM, Song B, Kwon KI, Kim SK, Yun HY. Prediction of Methionine and Homocysteine levels in Zucker diabetic fatty (ZDF) rats as a DIS_DIAB animal model after consumption of a Methionine-rich diet. Nutr Metab (Lond). 2018;15:14. doi:10.1186/s12986-018-0247-1"
  vignette    <- "Han_2018_methionineMetabolismCycle"
  units       <- list(time = "hour", dosing = "mmol/kg", concentration = "mmol/L")

  covariateData <- list(
    DIS_DIAB = list(
      description        = "Type-2 diabetes mellitus indicator: 1 = Zucker Diabetic Fatty (ZDF/Gmi fa/fa) rat as DIS_DIAB model; 0 = non-diabetic control (ZDF/Gmi fa/?).",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed cohort indicator from Han 2018 Methods (Study design). Multiplicative-on-log-scale effects (Han 2018 Eqs 8-11) on K_SH, K_HM, K_HC, K_HP, K_PH; coefficients derived from the Table 1 control vs ZDF point estimates as theta = log(K_ZDF / K_control). K_CM, K_MC, K_MS, K_SS, K_EL, and Vc are common to both cohorts (Han 2018 Results, paragraph 1: 'predicted without consideration of disease status').",
      source_name        = "ZDF"
    )
  )

  population <- list(
    species       = "rat (Zucker Diabetic Fatty: ZDF/Gmi fa/fa diabetic + ZDF/Gmi fa/? control)",
    n_studies     = 1,
    disease_state = "Two cohorts: ZDF/Gmi fa/fa rats as a Type-2-diabetes mellitus model and ZDF/Gmi fa/? non-diabetic littermate controls (Han 2018 Methods, Study design).",
    dose_range    = "Single intravenous bolus 0.8 mmol/kg methionine (mixture of L-methionine 0.6 mmol/kg + L-methionine-d4 0.2 mmol/kg). Blood sampled at 0, 10, 30, 60, 120, 210, 300, and 420 minutes after administration.",
    regions       = "Republic of Korea (Chungnam National University, Daejeon; Seoul National University).",
    notes         = "NONMEM 7.3.0 with Perl-Speaks-NONMEM 4.3.0, FOCE-INTER estimation. The Han 2018 paper does not numerically specify the residual-error variances ('Residual variability was explained by a combined error model using both additive and proportional error (equation was not shown)'); the packaged model therefore omits residual error and is intended for typical-value + IIV simulation rather than refitting. The paper reports IIV (%) as 70.1 for K_HM and 45.4 for K_HP in Table 1, interpreted here as %CV with variance = log((CV/100)^2 + 1). Per-group subject counts are not given in the trimmed text. Vc = 0.15 L/kg is reported as the methionine plasma volume; the same Vc is reused for the homocysteine plasma compartment in the absence of a separately-reported homocysteine volume (see vignette Assumptions)."
  )

  ini({
    # Methionine compartment-transfer and elimination rate constants (Han 2018 Table 1; common to both cohorts -- see Results paragraph 1)
    lkcm <- log(0.13);          label("Rate constant: methionine plasma -> hepatic (K_CM, 1/h)")     # Han 2018 Table 1 (point estimate 0.13, RSE 16.8 %)
    lkmc <- fixed(log(0.0629)); label("Rate constant: methionine hepatic -> plasma (K_MC, 1/h)")     # Han 2018 Table 1 (no RSE reported; fixed during estimation)
    lkms <- log(0.605);         label("Rate constant: methionine hepatic -> SAM via MAT (K_MS, 1/h)")# Han 2018 Table 1 (point estimate 0.605, RSE 9.5 %)
    lkss <- fixed(log(3220));   label("Rate constant: SAM -> SAH via methyltransferases (K_SS, 1/h)")# Han 2018 Table 1 (no RSE reported; fixed during estimation)
    lkel <- fixed(log(0.245));  label("Methionine elimination rate constant (K_EL, 1/h)")            # Han 2018 Table 1 (no RSE reported; fixed during estimation)
    lvc  <- fixed(log(0.15));   label("Apparent volume of distribution for plasma compartments (V_c, L/kg)")  # Han 2018 Table 1 (no RSE reported; fixed during estimation)

    # SAH -> homocysteine and homocysteine-fate rate constants (Han 2018 Table 1 control-rat point estimates)
    lksh <- log(11.6);          label("Rate constant: SAH -> homocysteine via SAHH (K_SH, 1/h)")      # Han 2018 Table 1 control (point estimate 11.6, RSE 20.5 %)
    lkhm <- log(30.7);          label("Rate constant: homocysteine -> methionine via BHMT (K_HM, 1/h)")# Han 2018 Table 1 control (point estimate 30.7, RSE 31.1 %)
    lkhc <- log(11.1);          label("Rate constant: homocysteine -> cysteine via CBS (K_HC, 1/h)") # Han 2018 Table 1 control (point estimate 11.1, RSE 4.8 %)
    lkhp <- log(142);           label("Rate constant: homocysteine hepatic -> plasma (K_HP, 1/h)")    # Han 2018 Table 1 control (point estimate 142, RSE 18.5 %)
    lkph <- log(5.13);          label("Rate constant: homocysteine plasma -> hepatic (K_PH, 1/h)")    # Han 2018 Table 1 control (point estimate 5.13, RSE 10.4 %)

    # DIS_DIAB (ZDF) covariate effects (multiplicative on log scale, Han 2018 Eqs 8-11).
    # Coefficient = log(K_ZDF / K_control); ZDF cohort recovers when DIS_DIAB = 1.
    # Han 2018 Eq 8 (K_HM) and Eq 9 (K_HP) combine the etas with the ZDF effect;
    # Eq 10 covers one of the remaining rate constants (printed as K_HC); the source's
    # Eq 11 is reproduced for the last cohort-different rate constant in Table 1.
    # The five coefficients below correspond, one-for-one, to the five Table 1 rate
    # constants whose point estimates differ between cohorts.
    e_t2dm_ksh <- fixed(log(13.5 / 11.6));  label("DIS_DIAB (ZDF) log-scale effect on K_SH")  # log(13.5/11.6) = +0.1519; Han 2018 Table 1
    e_t2dm_khm <- fixed(log(2.54 / 30.7));  label("DIS_DIAB (ZDF) log-scale effect on K_HM (Eq 8)")  # log(2.54/30.7) = -2.4910; Han 2018 Table 1
    e_t2dm_khc <- fixed(log(0.52 / 11.1));  label("DIS_DIAB (ZDF) log-scale effect on K_HC (Eq 10)") # log(0.52/11.1) = -3.0606; Han 2018 Table 1
    e_t2dm_khp <- fixed(log(20.4 / 142));   label("DIS_DIAB (ZDF) log-scale effect on K_HP (Eq 9)")  # log(20.4/142)  = -1.9405; Han 2018 Table 1
    e_t2dm_kph <- fixed(log(0.032 / 5.13)); label("DIS_DIAB (ZDF) log-scale effect on K_PH")         # log(0.032/5.13)= -5.0760; Han 2018 Table 1

    # Inter-individual variability (Han 2018 Table 1 reports IIV as %; interpreted as CV%).
    # variance = log((CV/100)^2 + 1) for a log-normal random effect on the underlying rate constant.
    etalkhm ~ log((70.1/100)^2 + 1)  # Han 2018 Table 1 IIV 70.1 %
    etalkhp ~ log((45.4/100)^2 + 1)  # Han 2018 Table 1 IIV 45.4 %

    # Residual variability: Han 2018 ('Residual variability was explained by a combined
    # error model using both additive and proportional error (equation was not shown)')
    # does not give numeric values, so residual error is intentionally omitted; the model
    # is intended for typical-value + IIV simulation. See population$notes.
  })

  model({
    # Individual rate constants (1/h)
    kcm <- exp(lkcm)
    kmc <- exp(lkmc)
    kms <- exp(lkms)
    kss <- exp(lkss)
    kel <- exp(lkel)
    vc  <- exp(lvc)

    # Cohort-stratified rate constants (DIS_DIAB = 0 -> control; DIS_DIAB = 1 -> ZDF).
    ksh <- exp(lksh           + e_t2dm_ksh * DIS_DIAB)
    khm <- exp(lkhm + etalkhm + e_t2dm_khm * DIS_DIAB)
    khc <- exp(lkhc           + e_t2dm_khc * DIS_DIAB)
    khp <- exp(lkhp + etalkhp + e_t2dm_khp * DIS_DIAB)
    kph <- exp(lkph           + e_t2dm_kph * DIS_DIAB)

    # ODE system (Han 2018 Eqs 1-7). Compartment paper-mapping:
    #   met    = A(1) methionine systemic circulation (plasma; observed)
    #   methep = A(2) methionine hepatic / metabolically-active compartment
    #   sam    = A(3) S-adenosylmethionine
    #   sah    = A(4) S-adenosyl-L-homocysteine
    #   hcyhep = A(5) homocysteine hepatic / metabolically-active compartment
    #   cys    = A(6) cysteine (terminal sink)
    #   hcy    = A(7) homocysteine systemic circulation (plasma; observed)
    # Eq 6 source-side prints `dA(6)/dt = K_HC * A(4)`; corrected here to `K_HC * A(5)`
    # to match Fig 1 (homocysteine -> cysteine via CBS) and the mass-balance of Eq 5
    # (which already debits K_HC * A(5) from the homocysteine hepatic pool). The cysteine
    # compartment is a terminal sink with no observation and no feedback, so the
    # correction does not affect predictions for the observed methionine / homocysteine
    # compartments. Documented in the vignette Errata.
    d/dt(met)    <- -kcm * met    + kmc * methep
    d/dt(methep) <-  kcm * met    - kmc * methep - kms * methep - kel * methep + khm * hcyhep
    d/dt(sam)    <-  kms * methep - kss * sam
    d/dt(sah)    <-  kss * sam    - ksh * sah
    d/dt(hcyhep) <-  ksh * sah    - khm * hcyhep - khc * hcyhep - khp * hcyhep + kph * hcy
    d/dt(cys)    <-  khc * hcyhep
    d/dt(hcy)    <-  khp * hcyhep - kph * hcy

    # Observation variables (plasma concentrations in mmol/L = mM).
    Cmet <- met / vc
    Chcy <- hcy / vc
  })
}
