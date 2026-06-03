Lakota_2018_ertapenem <- function() {
  description <- "Three-compartment population PK model for ertapenem in adults across a wide range of body sizes (Lakota 2018)"
  reference <- "Lakota EA, Landersdorfer CB, Zhang L, Nafziger AN, Bertino JS Jr, Bhavnani SM, Forrest A. Population Pharmacokinetic Analyses for Ertapenem in Subjects with a Wide Range of Body Sizes. Antimicrob Agents Chemother. 2018;62(10):e00784-18. doi:10.1128/AAC.00784-18"
  vignette <- "Lakota_2018_ertapenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on CL with reference 95.90 kg (dataset median); exponent 0.278 (Lakota 2018 Table 3, footnote b).",
      source_name        = "WT"
    ),
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on Vc with reference 2.06 m^2 (dataset median); exponent 1.86 (Lakota 2018 Table 3, footnote c). Strongly correlated with WT in the dataset (r = 0.966).",
      source_name        = "BSA"
    )
  )

  # Covariates screened in the formal forward-selection / backward-elimination
  # procedure (Lakota 2018 Methods) but not retained in the final model. Listed
  # here so the source-paper covariate screen is preserved for provenance.
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL (Lakota 2018 Table 2 Step 1); P = 0.05823, did not meet the alpha = 0.01 forward-selection threshold.",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = male",
      notes              = "Screened as a binary effect on PK parameters (Lakota 2018 Methods); not retained. Renamed from source column SEX to the canonical SEXF per covariate-columns.md.",
      source_name        = "SEX"
    ),
    HT = list(
      description        = "Subject height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a body-size measure (Lakota 2018 Methods); not retained.",
      source_name        = "HT"
    ),
    BMI = list(
      description        = "Body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a body-size measure (Lakota 2018 Table 2); P < 0.001 on Vc but BSA was preferred in step 1, and P = 0.00046 on CL but WT was preferred in step 2.",
      source_name        = "BMI"
    ),
    IBW = list(
      description        = "Ideal body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a body-size measure (Lakota 2018 Table 2); P = 0.04472 on Vc, did not meet alpha = 0.01 forward-selection threshold.",
      source_name        = "IBW"
    ),
    FFM = list(
      description        = "Fat-free mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a body-size measure (Lakota 2018 Table 2); inferior to BSA on Vc and to WT on CL.",
      source_name        = "FFM"
    ),
    CRCL = list(
      description        = "Creatinine clearance (Cockcroft-Gault) and BSA-normalized creatinine clearance (CLcrn) were both screened",
      units              = "mL/min and mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as renal-function covariates on CL (Lakota 2018 Methods); not retained because the cohort lacked any subjects with renal impairment (CLcr range 83-267 mL/min; CLcrn range 59.9-210). Discussion notes other ertapenem popPK models (Burkhardt et al., Dailly et al.) do include a CLcr-CL relationship.",
      source_name        = "CLCR, CLCRN"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 30,
    n_studies       = 1,
    age_range       = "19-50 years",
    age_median      = "37 years",
    weight_range    = "49.8-140 kg",
    weight_median   = "95.9 kg",
    bsa_range       = "1.47-2.58 m^2",
    bsa_median      = "2.06 m^2",
    bmi_range       = "18.8-49.9 kg/m^2",
    sex_female_pct  = 50,
    race_ethnicity  = "Not reported in baseline demographics (Lakota 2018 Table 1)",
    disease_state   = "Healthy adult volunteers; 10 normal-weight (BMI 18.5-24.9), 10 obese (BMI 30-39.9), and 10 morbidly obese (BMI >= 40 kg/m^2)",
    dose_range      = "Single 1 g IV infusion over 30 min",
    regions         = "United States (single-site clinical study)",
    notes           = "Subjects re-analysed from the Chen et al. NCA study (Lakota 2018 reference 11). 356 ertapenem serum samples; all above LLOQ; no outliers excluded. See Lakota 2018 Table 1 for full baseline demographics."
  )

  ini({
    # Structural parameters - reference values for a subject at the dataset
    # median (WT = 95.90 kg, BSA = 2.06 m^2). All from Lakota 2018 Table 3
    # ("Population mean / Final estimate" column).
    lcl  <- log(1.79);  label("Clearance for the reference subject (CL, L/h)")     # Table 3
    lvc  <- log(4.76);  label("Central volume of distribution (Vc, L)")            # Table 3
    lq   <- log(6.71);  label("Inter-compartmental clearance 1 (CLd1, L/h)")       # Table 3
    lvp  <- log(2.96);  label("First peripheral volume of distribution (Vp1, L)")  # Table 3
    lq2  <- log(0.296); label("Inter-compartmental clearance 2 (CLd2, L/h)")       # Table 3
    lvp2 <- log(1.1);   label("Second peripheral volume of distribution (Vp2, L)") # Table 3

    # Covariate effects - power functions, reference values at dataset medians.
    # CL  = 1.79 * (WT/95.90)^0.278   (Table 3, footnote b)
    # Vc  = 4.76 * (BSA/2.06)^1.86    (Table 3, footnote c)
    e_wt_cl  <- 0.278; label("Power exponent of WT on CL (unitless, reference 95.90 kg)")  # Table 3
    e_bsa_vc <- 1.86;  label("Power exponent of BSA on Vc (unitless, reference 2.06 m^2)") # Table 3

    # IIV: omega^2 = log(CV^2 + 1).
    # Lakota 2018 Table 3 reports CV% (with %SEM as a relative SE on the CV
    # estimate, not on the variance):
    #   CL:        CV = 9.57%  -> omega^2 = log(0.0957^2 + 1) ~= 0.00912
    #   Vc, Vp1:   CV = 7.92%  shared between Vc and Vp1 (Table 3 footnote d;
    #              correlation 1.00 collapsed to a single OMEGA term during
    #              full-multivariable model evaluation - Results paragraph 4).
    # The shared OMEGA is applied to both lvc and lvp in model() below; only
    # one eta (etalvc) is declared in ini().
    etalcl ~ 0.00912
    etalvc ~ 0.00625

    # Residual error - proportional only. Lakota 2018 Table 3 reports
    # "Proportional error 0.0166 (12.9 % CV)"; the 0.0166 is sigma^2 and
    # sqrt(0.0166) = 0.1288 is the CV / SD on the linear scale.
    # The additive residual term was dropped from the final model
    # (estimated < 1e-8 during development, Results paragraph "Final
    # population pharmacokinetic model").
    propSd <- 0.129; label("Proportional residual error (fraction)") # Table 3
  })
  model({
    # Individual PK parameters with reference-centered power covariates.
    # The same etalvc is added to lvp because Lakota 2018 retained a single
    # OMEGA term shared between Vc and Vp1.
    cl  <- exp(lcl + etalcl) * (WT  / 95.90)^e_wt_cl
    vc  <- exp(lvc + etalvc) * (BSA / 2.06)^e_bsa_vc
    q   <- exp(lq)
    vp  <- exp(lvp + etalvc)
    q2  <- exp(lq2)
    vp2 <- exp(lvp2)

    # Micro-constants
    kel <- cl  / vc
    k12 <- q   / vc
    k21 <- q   / vp
    k13 <- q2  / vc
    k31 <- q2  / vp2

    # 3-compartment linear PK; dose is IV infusion to central (rate set on
    # the dose record; no depot compartment).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-                                                      k13 * central - k31 * peripheral2

    # Concentration in mg/L (dose in mg, volume in L)
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
