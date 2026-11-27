Krause_2026_warfarin <- function() {
  description <- "Two-compartment population PK model for oral warfarin in adults (Krause 2026), with first-order absorption after a lag time, linear elimination, and body weight on clearance, intercompartmental clearance and both volumes through fixed allometric exponents (0.75 on flows, 1 on volumes; reference 75 kg). Fitted in Monolix to the single-dose O'Reilly warfarin dataset and used in the paper as the adult reference model for deriving pediatric doses by steady-state exposure matching."
  reference <- paste(
    "Krause A, Celliere G. Strategies for pediatric dose derivation from",
    "population pharmacokinetic models. J Pharmacokinet Pharmacodyn.",
    "2026;53(5):46. doi:10.1007/s10928-026-10053-8.",
    "Parameter values are the full-precision estimates of the deposited",
    "Monolix run Warfarin-TwoCompLinPKAllomScalingBW",
    "(Supplementary Material 2, populationParameters.txt), which Table 2",
    "of the article prints rounded to two significant figures.",
    "Warfarin data from O'Reilly RA, Aggeler PM. Circulation.",
    "1968;38(1):169-177. doi:10.1161/01.cir.38.1.169.",
    sep = " "
  )
  vignette <- "Krause_2026_warfarin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline body weight; power model (WT/75)^beta on CL, Q, Vc and Vp with the exponents fixed to 0.75 (flows) and 1 (volumes). Reference weight 75 kg (Supplementary Material 1 covariate equations; mlxtran logtWt = log(wt/75)).",
      source_name = "wt"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex, 1 = female",
      units = "(binary)",
      type = "binary",
      reference_category = "male",
      notes = "Present in the deposited dataset (column 'sex', 5 of 32 subjects coded 0, matching the 5 females of the O'Reilly cohort) and declared as a covariate in the Monolix project, but not used in any parameter model.",
      source_name = "sex"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Present in the deposited dataset and declared as a covariate in the Monolix project, but not used in any parameter model; the paper scales by body size only and assumes maturation is complete.",
      source_name = "age"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "warfarin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "warfarin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "warfarin", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 32L,
    n_studies = 1L,
    age_range = "21-63 years",
    age_median = "27.5 years (deposited dataset)",
    weight_range = "40-102 kg (deposited dataset)",
    weight_median = "71.65 kg (deposited dataset)",
    sex_female_pct = 15.6,
    disease_state = "Adults receiving a single oral warfarin loading dose (O'Reilly warfarin dataset).",
    dose_range = "Single oral dose of 1.5 mg/kg, 60-153 mg.",
    observations = "247 warfarin plasma concentrations: 1, 2, 3, 6, 9, 12 h post-dose in 14 subjects and 24, 36, 48, 72, 96, 120 h in all subjects (Application section).",
    notes = "Age range, dose and sampling scheme from the Application section of Krause 2026. Weight range, median age and weight, and the female percentage (5 of 32) are computed from the deposited dataset (Supplementary Material 2, Monolix/data/warfarin_data.csv); the article prints no demographics table."
  )

  ini({
    # Structural parameters for a 75 kg adult. Values are the deposited
    # populationParameters.txt estimates; Table 2 prints them rounded.
    ltlag <- log(0.735187266103877); label("Absorption lag time (h)") # Table 2 tlag 0.74 h; populationParameters.txt Tlag_pop
    lka <- log(0.959138951210313); label("Absorption rate constant (1/h)") # Table 2 ka 0.96 1/h; populationParameters.txt ka_pop
    lcl <- log(0.139322460947242); label("Apparent clearance for a 75 kg adult (L/h)") # Table 2 CL 0.14 L/h; populationParameters.txt Cl_pop
    lvc <- log(6.56729256615457); label("Apparent central volume for a 75 kg adult (L)") # Table 2 Vc 6.57 L; populationParameters.txt V1_pop
    lq <- log(3.8387377539771); label("Apparent intercompartmental clearance for a 75 kg adult (L/h)") # Table 2 Q 3.84 L/h; populationParameters.txt Q_pop
    lvp <- log(2.12301525481965); label("Apparent peripheral volume for a 75 kg adult (L)") # Table 2 Vp 2.12 L; populationParameters.txt V2_pop

    # Allometric exponents, fixed (mlxtran method=FIXED)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL (unitless)") # Table 2 beta CL, body weight 0.75
    e_wt_q <- fixed(0.75); label("Allometric exponent of body weight on Q (unitless)") # Table 2 beta Q, body weight 0.75
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on Vc (unitless)") # Table 2 beta Vc, body weight 1.00
    e_wt_vp <- fixed(1); label("Allometric exponent of body weight on Vp (unitless)") # Table 2 beta Vp, body weight 1.00

    # IIV: Monolix reports omega as the SD of the log-normal random effect;
    # the variances below are omega^2 and reproduce the Table 2 %CV column.
    etaltlag ~ 0.310774 # Table 2 tlag omega 0.56 (60.4 %CV); populationParameters.txt omega_Tlag 0.557471
    etalka ~ 0.545827 # Table 2 ka omega 0.74 (85.2 %CV); populationParameters.txt omega_ka 0.738801
    etalcl ~ 0.0682129 # Table 2 CL omega 0.26 (26.6 %CV); populationParameters.txt omega_Cl 0.261176
    etalvc ~ 0.00482925 # Table 2 Vc omega 0.07 (7.00 %CV); populationParameters.txt omega_V1 0.0694928
    etalq ~ 6.21021 # Table 2 Q omega 2.49 (2228.9 %CV); populationParameters.txt omega_Q 2.492029
    etalvp ~ 0.202592 # Table 2 Vp omega 0.45 (47.4 %CV); populationParameters.txt omega_V2 0.450102

    # Residual error: Monolix combined1, y = Cc + (a + b * Cc) * e
    addSd <- 0.207331015689362; label("Additive residual error (mg/L)") # Table 2 a 0.21; populationParameters.txt a1
    propSd <- 0.0485207465519773; label("Proportional residual error (fraction)") # Table 2 b 0.05; populationParameters.txt b1
  })

  model({
    # Individual parameters; body weight normalised to 75 kg
    # (Supplementary Material 1 covariate equations)
    tlag <- exp(ltlag + etaltlag)
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 75)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc
    q <- exp(lq + etalq) * (WT / 75)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 75)^e_wt_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Monolix pkmodel(Tlag, ka, V, Cl, k12, k21), written out as in
    # Supplementary Material 1
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag

    # Dose in mg and volumes in L give mg/L, the unit of the dataset
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
