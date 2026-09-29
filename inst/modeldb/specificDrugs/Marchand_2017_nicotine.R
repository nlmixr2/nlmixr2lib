Marchand_2017_nicotine <- function() {
  description <- "Two-compartment population PK model for nicotine with zero-order absorption and an additive mono-exponentially decaying background-exposure component, in healthy adult smokers using heated tobacco, conventional cigarettes, nasal spray or gum. Solve with rxSolve(useLinCmt = FALSE): rxode2's default ODE-to-linCmt conversion silently discards peripheral1 from an explicit k12 / k21 parameterisation, leaving AUC correct but the terminal phase mono-exponential."
  reference <- "Marchand M, Brossard P, Merdjan H, Lama N, Weitkunat R, Ludicke F. Nicotine Population Pharmacokinetics in Healthy Adult Smokers: A Retrospective Analysis. Eur J Drug Metab Pharmacokinet. 2017;42(6):943-954. doi:10.1007/s13318-017-0405-2"
  vignette <- "Marchand_2017_nicotine"
  units <- list(time = "min", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on relative bioavailability, referenced to the typical subject weight of 69.1 kg (Marchand 2017 Table 4 and Sect. 3.4). Time-fixed baseline weight.",
      source_name = "WT"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Exponential effect on apparent clearance; Cl/F is 26 pct higher in females than in males (Marchand 2017 Table 4).",
      source_name = "SEX"
    ),
    CYP2A6 = list(
      description = "Baseline CYP2A6 metabolic activity, expressed as a percentage",
      units = "percent",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on Cl/F and on the background baseline concentration C0, referenced to the typical subject value of 29.2 percent (Marchand 2017 Table 4). Learning-dataset mean 31.5 (SD 18.2) percent, Table 2. CYP2A6 is the major enzyme metabolizing nicotine; 2 subjects with a missing value were excluded from the covariate analysis.",
      source_name = "CYP2A6 activity"
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Black)",
      notes = "Exponential effect on the background baseline concentration C0; being Black increased C0 by 50 pct (Marchand 2017 Table 4). The authors note this effect could not be distinguished from a study effect because all Black subjects in the learning dataset came from a single study (Sect. 4).",
      source_name = "RACE (black versus non-black)"
    ),
    DOSE_NICOTINE_MG = list(
      description = "Nominal nicotine dose of the product used on this occasion",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Nominal nicotine ISO yield for inhaled products (0.5 mg for THS, 0.1-1.5 mg for conventional cigarettes depending on brand); 1 mg for the nasal spray and 2 mg for the gum (Marchand 2017 Sect. 2.2). Power effect on relative bioavailability referenced to 0.5 mg. Must equal the dose amount on the corresponding dose record.",
      source_name = "nicotine ISO yield / nicotine dose"
    ),
    FORM_NICOTINE_CC = list(
      description = "Conventional cigarette product indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Tobacco Heating System, the reference product)",
      notes = "Per-occasion indicator. Exponential effect on relative bioavailability (Marchand 2017 Table 4).",
      source_name = "nature of product = CC"
    ),
    FORM_NICOTINE_NNS = list(
      description = "Nicotine nasal spray product indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Tobacco Heating System, the reference product)",
      notes = "Per-occasion indicator. Exponential effect on relative bioavailability (Marchand 2017 Table 4).",
      source_name = "nature of product = NNS"
    ),
    FORM_NICOTINE_GUM = list(
      description = "Nicotine gum product indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Tobacco Heating System, the reference product)",
      notes = "Per-occasion indicator. Exponential effects on relative bioavailability and on the zero-order absorption duration (Marchand 2017 Table 4). All gum in the analysis was the mentholated variant, so FORM_NICOTINE_MENTHOL is also 1 on gum occasions.",
      source_name = "nature of product = GUM"
    ),
    FORM_NICOTINE_MENTHOL = list(
      description = "Mentholated product variant indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (regular, non-mentholated variant)",
      notes = "Per-occasion indicator, orthogonal to the product-type indicators (mentholated THS, mentholated cigarette and the gum are all menthol variants). Exponential effects on the apparent central volume and on the zero-order absorption duration (Marchand 2017 Table 4). The authors caution that menthol may be confounded with region, product and product-use behaviour (Sect. 4).",
      source_name = "menthol"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "year",
      type = "continuous",
      notes = "Screened on V1/F as a secondary covariate; not retained in the final model (Marchand 2017 Sect. 3.3 and 3.4)."
    ),
    HT = list(
      description = "Body height",
      units = "m",
      type = "continuous",
      notes = "Screened on Tdur as a secondary covariate; not retained (Marchand 2017 Sect. 3.3)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m2",
      type = "continuous",
      notes = "Screened as a secondary covariate; not retained (Marchand 2017 Sect. 3.4)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened as a secondary covariate; not retained (Marchand 2017 Sect. 3.4)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened as a secondary covariate; not retained (Marchand 2017 Sect. 3.4)."
    ),
    BILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened as a secondary covariate; not retained (Marchand 2017 Sect. 3.4)."
    ),
    CRCL = list(
      description = "Creatinine clearance, Cockcroft-Gault",
      units = "mL/min",
      type = "continuous",
      notes = "Screened as a primary covariate on nicotine clearance; not retained (Marchand 2017 Suppl. Table 1 and Sect. 3.4)."
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened on Tdur (white versus non-white); not retained (Marchand 2017 Sect. 3.3)."
    ),
    SMOKE_CPD_SCORE = list(
      description = "Fagerstrom cigarettes-per-day score",
      units = "(ordinal score 0-3)",
      type = "count",
      notes = "Daily cigarette use at baseline screened as a secondary covariate; not retained (Marchand 2017 Sect. 3.4)."
    ),
    SMOKE_TTFC_SCORE = list(
      description = "Fagerstrom time-to-first-cigarette score",
      units = "(ordinal score 0-3)",
      type = "count",
      notes = "First FTND item screened as a secondary covariate; not retained (Marchand 2017 Sect. 3.4)."
    )
  )

  compartmentData <- list(
    central = list(analyte = "nicotine", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "nicotine", units = "ug", specimen = "plasma", verified = FALSE)
  )

  population <- list(
    species = "human",
    n_subjects = 246,
    n_studies = 4,
    age_range = "21-66 years",
    age_median = "33.5 years (mean, SD 9.23)",
    weight_median = "70.1 kg (mean, SD 13.9)",
    sex_female_pct = 45.5,
    race_ethnicity = c(White = 35, Black = 13.8, Asian = 50.8, Other = 0.4),
    disease_state = "healthy adult smokers of at least 10 cigarettes per day",
    dose_range = "single use of Tobacco Heating System or conventional cigarette (nicotine ISO yield 0.1-1.5 mg), 1 mg regular nicotine nasal spray, or 2 mg mentholated nicotine gum",
    regions = "USA (25.2 pct), EU (24.4 pct), Japan (50.4 pct)",
    notes = "The learning dataset used to build this model comprised 246 subjects and 6843 measurable concentrations from 4 randomized two-period crossover single-product-use trials (Marchand 2017 Table 1 and Sect. 3.1-3.2); baseline covariates are Table 2. The covariate model was developed in the 244 subjects with complete covariate information. A further 456 subjects from 4 ad libitum use studies formed an external validation dataset (702 subjects overall); in that dataset every parameter was held at the learning-dataset estimate except C0, which was re-estimated at 2.10 ng/mL because the smoking abstinence period was much shorter (Sect. 3.5.2)."
  )

  ini({
    # Structural parameters. Typical subject (Marchand 2017 Sect. 3.4): male,
    # not Black, 69.1 kg, using the regular variant of THS at a 0.5 mg nicotine
    # ISO yield, baseline CYP2A6 activity 29.2 percent.
    lvc <- log(70.0); label("Apparent central volume of distribution V1/F (L)") # Table 3: 70.0 L, RSE 2.8
    lcl <- log(0.407); label("Apparent clearance Cl/F (L/min)") # Table 3: 0.407 L/min, RSE 3.0
    lvp <- log(171); label("Apparent peripheral volume of distribution V2/F (L)") # Table 3: 171 L, RSE 3.3
    lq <- log(0.171); label("Apparent inter-compartmental clearance Cl2/F (L/min)") # Table 3: 0.171 L/min, RSE 3.5
    ld1 <- log(5.30); label("Duration of zero-order absorption Tdur (min)") # Table 3: 5.30 min, RSE 1.2
    lfrel <- fixed(log(1)); label("Relative bioavailability of the reference product THS (fraction)") # Table 3 row Frel-THS: 1, reference product
    lrbase <- log(0.358); label("Background nicotine concentration at first product use C0 (ng/mL)") # Table 3: 0.358 ng/mL, RSE 1.5

    # Covariate effects. Every coefficient below is on the log scale exactly as
    # Table 3 reports it; the corresponding closed forms are Table 4.
    e_cyp2a6_cl <- 0.322; label("Power exponent on (CYP2A6 / 29.2) for Cl/F (unitless)") # Table 3 dCldCYP2A6: 0.322, RSE 1.6
    e_sexf_cl <- 0.235; label("Log-scale shift in Cl/F for female subjects (unitless)") # Table 3 dCldSEX (female): 0.235, RSE 1.6
    e_dose_nicotine_mg_frel <- -0.573; label("Power exponent on (DOSE_NICOTINE_MG / 0.5) for Frel (unitless)") # Table 3 dFreldDOSE: -0.573, RSE 1.6
    e_wt_frel <- -0.715; label("Power exponent on (WT / 69.1) for Frel (unitless)") # Table 3 dFreldWT: -0.715, RSE 1.6
    e_form_nicotine_cc_frel <- 0.0189; label("Log-scale shift in Frel for conventional cigarette versus THS (unitless)") # Table 3 dFreld-CC: 0.0189, RSE 1.6
    e_form_nicotine_nns_frel <- -1.42; label("Log-scale shift in Frel for nicotine nasal spray versus THS (unitless)") # Table 3 dFreld-NNS: -1.42, RSE 1.5
    e_form_nicotine_gum_frel <- -0.489; label("Log-scale shift in Frel for nicotine gum versus THS (unitless)") # Table 3 dFreld-GUM: -0.489, RSE 1.6
    e_form_nicotine_menthol_vc <- 0.0912; label("Log-scale shift in V1/F for mentholated variants (unitless)") # Table 3 dVdMENTH: 0.0912, RSE 1.6
    e_cyp2a6_rbase <- -0.401; label("Power exponent on (CYP2A6 / 29.2) for C0 (unitless)") # Table 3 dC0dCYP2A6: -0.401, RSE 1.6
    e_race_black_rbase <- 0.408; label("Log-scale shift in C0 for Black subjects (unitless)") # Table 3 dC0dBLACK: 0.408, RSE 1.6
    e_form_nicotine_menthol_d1 <- 0.0530; label("Log-scale shift in Tdur for mentholated variants (unitless)") # Table 3 dTdurdMENTH: 0.0530, RSE 1.6
    e_form_nicotine_gum_d1 <- 2.14; label("Log-scale shift in Tdur for nicotine gum (unitless)") # Table 3 dTdurd-GUM: 2.14, RSE 1.2

    # Inter-individual variability. Table 3 reports the log-scale variance in
    # its 'Variance' column and a companion 'IIV' column that is its square
    # root expressed as a percentage, so the variances below are used as
    # printed without any CV-to-variance conversion.
    etalvc ~ 0.641 # Table 3 V1/F variance 0.641, RSE 4.1, shrinkage 1.9
    etalcl ~ 0.467 # Table 3 Cl/F variance 0.467, RSE 4.2, shrinkage 2.5
    etalvp ~ 0.715 # Table 3 V2/F variance 0.715, RSE 4.4, shrinkage 29.4
    etalq ~ 1.93 # Table 3 Cl2/F variance 1.93, RSE 4.4, shrinkage 14.3
    etald1 ~ 0.141 # Table 3 Tdur variance 0.141, RSE 4.1, shrinkage 12.4
    etalrbase ~ 0.233 # Table 3 C0 variance 0.233, RSE 4.5, shrinkage 22.6
    etalfrel ~ 0.489 # Table 3 prints one shared variance 0.489 on each of the three dFreld rows and none on the Frel-THS row

    # Residual error. The base model comparison in Supplementary Table 2 shows a
    # log-additive residual beating both proportional and combined models by a
    # large margin, which is a log-normal error in nlmixr2 terms. Phoenix NLME
    # reports this term in the fixed-effect estimate column, i.e. as a standard
    # deviation rather than a variance.
    expSd <- 0.289; label("Residual standard deviation on the natural-log concentration scale (log ng/mL)") # Table 3 row Residual (log domain): 0.289, RSE 0.9
  })

  model({
    # Doses are supplied in mg while the states are carried in ug, so that
    # amount / volume with the volume in L is directly in ng/mL.
    ug_per_mg <- 1000

    # Covariate reference values, Marchand 2017 Table 4 and Sect. 3.4.
    wt_ref <- 69.1 # kg
    cyp2a6_ref <- 29.2 # percent CYP2A6 activity
    dose_ref <- 0.5 # mg nicotine ISO yield of THS

    # 1. Individual parameters
    vc <- exp(lvc + etalvc) * exp(e_form_nicotine_menthol_vc * FORM_NICOTINE_MENTHOL)
    cl <- exp(lcl + etalcl) * (CYP2A6 / cyp2a6_ref)^e_cyp2a6_cl * exp(e_sexf_cl * SEXF)
    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)
    d1 <- exp(ld1 + etald1) *
      exp(
        e_form_nicotine_menthol_d1 * FORM_NICOTINE_MENTHOL +
          e_form_nicotine_gum_d1 * FORM_NICOTINE_GUM
      )
    rbase <- exp(lrbase + etalrbase) *
      (CYP2A6 / cyp2a6_ref)^e_cyp2a6_rbase *
      exp(e_race_black_rbase * RACE_BLACK)

    # Relative bioavailability. THS is the reference product with Frel fixed to
    # 1 and, per Table 3, carries no random effect: the single shared Frel
    # variance is printed on the three product-difference rows only, so the
    # random effect is gated to the non-THS products.
    nonths <- FORM_NICOTINE_CC + FORM_NICOTINE_NNS + FORM_NICOTINE_GUM
    frel <- exp(
      lfrel +
        e_form_nicotine_cc_frel * FORM_NICOTINE_CC +
        e_form_nicotine_nns_frel * FORM_NICOTINE_NNS +
        e_form_nicotine_gum_frel * FORM_NICOTINE_GUM +
        etalfrel * nonths
    ) *
      (DOSE_NICOTINE_MG / dose_ref)^e_dose_nicotine_mg_frel *
      (WT / wt_ref)^e_wt_frel

    # 2. Micro-constants, and the macroscopic terminal rate constant that drives
    # the background sub-model. Marchand 2017 Fig. 2 derives beta from the
    # disposition parameters rather than estimating it.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    ksum <- k12 + k21 + kel
    beta <- 0.5 * (ksum - sqrt(ksum * ksum - 4 * k21 * kel))

    # 3. Two-compartment linear disposition, Marchand 2017 Fig. 2
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Zero-order absorption of the product-use dose directly into the central
    # compartment over Tdur. Dose records must carry rate = -2 so that rxode2
    # uses the modelled duration dur(central) = d1.
    dur(central) <- d1
    f(central) <- frel * ug_per_mg

    # 5. Observation. The background sub-model is a mono-exponential decay from
    # C0 at the time of first product use, added to the product-use
    # concentration to give the total plasma nicotine that was fitted
    # (Marchand 2017 Fig. 2). Time zero is the time of first product use.
    Cprod <- central / vc
    Cbkgrd <- rbase * exp(-beta * t)
    Cc <- Cprod + Cbkgrd
    Cc ~ lnorm(expSd)
  })
}
