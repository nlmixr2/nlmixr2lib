Wattanakul_2020_piperaquine <- function() {
  description <- paste(
    "Population PK model for oral piperaquine in 1,000 African patients",
    "(mostly children) with uncomplicated falciparum malaria from a",
    "multicentre post-licensing pharmacovigilance study in Burkina Faso,",
    "Ghana, Mozambique and Tanzania (Wattanakul 2020). The structure is",
    "the Hoglund 2017 pooled meta-analysis model refitted to these data",
    "with a NONMEM frequentist prior: two transit absorption compartments",
    "(ka = ktr) feeding a three-compartment disposition model, fixed",
    "allometric body-weight scaling (exponent 0.75 on clearances, 1 on",
    "volumes, reference 54 kg), a sigmoid age-maturation function on",
    "elimination clearance, and a fixed dose-occasion increment in",
    "relative bioavailability. Between-subject and between-occasion",
    "variability on bioavailability and mean transit time. Doses are",
    "piperaquine base in mg; predictions are plasma piperaquine in ng/mL.",
    "Sister model files from the same paper:",
    "modellib('Wattanakul_2020_piperaquine_qtc') (absolute QTc Emax model)",
    "and modellib('Wattanakul_2020_piperaquine_dqtc') (change-from-baseline",
    "QTc Emax model).",
    sep = " "
  )
  reference <- paste(
    "Wattanakul T, Ogutu B, Kabanywanyi AM, Asante K-P, Oduro A, Adjei A,",
    "Sie A, Sevene E, Macete E, Compaore G, Valea I, Osei I, Winterberg M,",
    "Gyapong M, Adjuik M, Abdulla S, Owusu-Agyei S, White NJ, Day NPJ,",
    "Tinto H, Baiden R, Binka F, Tarning J. Pooled multicenter analysis of",
    "cardiovascular safety and population pharmacokinetic properties of",
    "piperaquine in African patients with uncomplicated falciparum malaria.",
    "Antimicrob Agents Chemother. 2020;64(7):e01848-19.",
    "doi:10.1128/AAC.01848-19. PMC7318010.",
    "Open Access under CC BY 4.0. The PK parameters are in Table 2; the",
    "structural prior is Hoglund 2017 (doi:10.1371/journal.pmed.1002212,",
    "modellib('Hoglund_2017_piperaquine')).",
    sep = " "
  )
  vignette <- "Wattanakul_2020_piperaquine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Fixed allometric scaling with exponent 0.75 on CL/F, Q1/F and Q2/F",
        "and 1 on Vc/F, Vp1/F and Vp2/F (Methods, 'Population",
        "pharmacokinetic modeling of piperaquine'), referenced to 54 kg",
        "(Table 2 footnote b: 'the typical individual in the prior",
        "population with a body weight of 54 kg').",
        sep = " "
      ),
      source_name = "body weight"
    ),
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Sigmoid enzyme maturation on elimination clearance, Equation 3:",
        "CL_i = TVCL * AGE^Hill / (MF50^Hill + AGE^Hill) with MF50 = 0.575",
        "years and Hill = 5.51, both carried fixed from the prior model",
        "(Table 2). The maturation factor is above 0.99 from about 1.3",
        "years of age.",
        sep = " "
      ),
      source_name = "AGE"
    ),
    OCC = list(
      description = "Dose-occasion index within one three-day course (1 = first daily dose, 2 = second, 3 = third)",
      units = "(count)",
      type = "count",
      reference_category = NULL,
      notes = paste(
        "Serves two roles. (1) Dose-occasion effect on relative",
        "bioavailability, fixed at 0.237 per consecutive dose (Table 2;",
        "Methods: '24% increased relative bioavailability between each",
        "consecutive dose'), encoded additively F_OCC = 1 + 0.237 * (OCC - 1)",
        "as in modellib('Hoglund_2017_piperaquine'). (2) Between-occasion",
        "variability on F and MTT (Equation 2), multiplexed onto",
        "etaiov_*_1 .. etaiov_*_3 sharing one variance per parameter;",
        "OCC outside 1..3 gives no occasion eta. Carry the OCC of the most",
        "recent dose on every observation row. For repeated monthly",
        "courses restart OCC at 1 for each course.",
        sep = " "
      ),
      source_name = "dose occasion"
    )
  )

  covariatesDataExcluded <- list(
    BODYTEMP = list(
      description = "Body temperature at enrolment",
      units = "degC",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search (Methods) and not retained: 'No additional significant covariate relationships were found for the patients studied here' (Results)."
    ),
    HGB = list(
      description = "Hemoglobin",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search and not retained. No coefficient reported."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search and not retained. No coefficient reported."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search and not retained. No coefficient reported."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search and not retained. No coefficient reported."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search and not retained. No coefficient reported."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search and not retained. No coefficient reported."
    ),
    POT = list(
      description = "Serum potassium",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search and not retained. No coefficient reported. Malaria parasite count and serum chloride were screened too (Methods) and not retained; they have no register column and are not listed here."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "piperaquine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "piperaquine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "piperaquine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 1000L,
    n_studies = 1L,
    n_observations = 2989L,
    age_range = "6 months and older; 0.6% < 1 y, 23.8% 1 to < 5 y, 45.0% 5 to < 12 y, 12.7% 12 to < 18 y, 17.9% >= 18 y (Table 1)",
    age_median = "7.5 years (IQR 5-12)",
    weight_range = "5 kg and above (inclusion criterion)",
    weight_median = "21 kg (IQR 15-38)",
    sex_female_pct = 51.8,
    race_ethnicity = "Black African (all participants enrolled at African sites; not tabulated)",
    disease_state = "Uncomplicated Plasmodium falciparum malaria. Severe malaria, pregnancy, lactation, QT-prolonging co-medication, family history of sudden death and baseline QTcF or QTcB above 450 ms were exclusion criteria.",
    dose_range = "Dihydroartemisinin-piperaquine (Eurartesim) once daily for 3 days by body-weight band (Table 5, old WHO regimen): 160 to 1,280 mg piperaquine phosphate per dose, directly observed, fasted 3 h before and after dosing.",
    regions = "Burkina Faso (n = 299), Ghana (n = 442), Mozambique (n = 89), Tanzania (n = 170); 10 sites",
    notes = paste(
      "Nested pharmacokinetic-ECG cohort of an 11,028-patient",
      "pharmacovigilance study (NCT02199951). 1 to 5 plasma samples per",
      "patient at about 0, 48, 52, 120, 144 and 168 h after the first dose;",
      "1.44% of samples below the 1.50 ng/mL LLOQ were omitted. Because",
      "sampling stopped at day 7 the Hoglund 2017 meta-analysis model",
      "(8,776 samples, 728 individuals) was used as a frequentist prior",
      "($PRIOR, NONMEM 7.3, FOCE-I).",
      sep = " "
    )
  )

  ini({
    # Structural parameters: Wattanakul 2020 Table 2 'Population estimate'
    # column, typical values at 54 kg.
    lmtt <- log(2.13)
    label("Mean transit time of the two-transit-compartment absorption chain MTT (h)")
    # Table 2: MTT (h) = 2.13 (%RSE 1.11; 95% CI 2.09-2.18)
    lcl <- log(53.1)
    label("Apparent elimination clearance CL/F at WT = 54 kg, fully mature (L/h)")
    # Table 2: CL/F (liter/h) = 53.1 (%RSE 2.77; 95% CI 50.2-56.1)
    lvc <- log(1730)
    label("Apparent central volume of distribution Vc/F at WT = 54 kg (L)")
    # Table 2: VC/F (liter) = 1,730 (%RSE 8.04; 95% CI 1,441-1,991)
    lq <- log(282)
    label("Apparent intercompartmental clearance to peripheral1 Q1/F at WT = 54 kg (L/h)")
    # Table 2: Q1/F (liter/h) = 282 (%RSE 5.60; 95% CI 249-310)
    lvp <- log(3290)
    label("Apparent first peripheral volume Vp1/F at WT = 54 kg (L)")
    # Table 2: VP1/F (liter) = 3,290 (%RSE 5.10; 95% CI 2,949-3,595)
    lq2 <- log(82.9)
    label("Apparent intercompartmental clearance to peripheral2 Q2/F at WT = 54 kg (L/h)")
    # Table 2: Q2/F (liter/h) = 82.9 (%RSE 2.42; 95% CI 78.9-86.6)
    lvp2 <- log(25100)
    label("Apparent second peripheral volume Vp2/F at WT = 54 kg (L)")
    # Table 2: VP2/F (liter) = 25,100 (%RSE 1.77; 95% CI 24,170-25,925)
    lfdepot <- fixed(log(1))
    label("Relative bioavailability F at the first dose occasion (unitless)")
    # Table 2: F = '1 fixed'

    e_doseocc_f <- fixed(0.237)
    label("Increment in relative bioavailability per consecutive dose occasion (fraction)")
    # Table 2: 'Dose occasion effect on F' = 0.237, held at the prior value
    mat_mf50 <- fixed(0.575)
    label("Age at 50 percent maturation of elimination clearance (years)")
    # Table 2: AGE50 (yr) = 0.575, held at the prior value
    mat_hill <- fixed(5.51)
    label("Hill coefficient of the clearance maturation function (unitless)")
    # Table 2: Hill = 5.51, held at the prior value

    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL/F, Q1/F and Q2/F (unitless)")
    # Methods: 'allometric function on all clearance (exponent of 0.75)'
    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on Vc/F, Vp1/F and Vp2/F (unitless)")
    # Methods: 'and volume of distribution (exponent of 1) parameters'

    # IIV / IOV: Table 2 footnote b, %CV = 100 * sqrt(exp(omega^2) - 1), so
    # omega^2 = log(CV^2 + 1). No IIV was reported on CL/F or Q1/F.
    etalfdepot ~ 0.136211
    # Table 2: F IIV 38.2% CV (%RSE 3.90); log(0.382^2 + 1) = 0.136211
    etalmtt ~ 0.131576
    # Table 2: MTT IIV 37.5% CV (%RSE 9.83); log(0.375^2 + 1) = 0.131576
    etalvc ~ 0.598301
    # Table 2: VC/F IIV 90.5% CV (%RSE 19.8); log(0.905^2 + 1) = 0.598301
    etalvp ~ 0.0533095
    # Table 2: VP1/F IIV 23.4% CV (%RSE 24.1); log(0.234^2 + 1) = 0.0533095
    etalq2 ~ 0.0708694
    # Table 2: Q2/F IIV 27.1% CV (%RSE 11.9); log(0.271^2 + 1) = 0.0708694
    etalvp2 ~ 0.0963315
    # Table 2: VP2/F IIV 31.8% CV (%RSE 1.01); log(0.318^2 + 1) = 0.0963315

    etaiov_fdepot_1 ~ 0.168209
    # Table 2: F IOV 42.8% CV (%RSE 2.07); log(0.428^2 + 1) = 0.168209
    etaiov_fdepot_2 ~ fixed(0.168209)
    # Same IOV variance as occasion 1
    etaiov_fdepot_3 ~ fixed(0.168209)
    # Same IOV variance as occasion 1
    etaiov_mtt_1 ~ 0.182162
    # Table 2: MTT IOV 44.7% CV (%RSE 1.20); log(0.447^2 + 1) = 0.182162
    etaiov_mtt_2 ~ fixed(0.182162)
    # Same IOV variance as occasion 1
    etaiov_mtt_3 ~ fixed(0.182162)
    # Same IOV variance as occasion 1

    # Residual error additive on the log scale (Methods), i.e. proportional
    # on the linear scale. Table 2 footnote a: sigma is a variance.
    propSd <- sqrt(0.198)
    label("Proportional residual SD of plasma piperaquine (log-scale SD)")
    # Table 2: sigma = 0.198 (variance, %RSE 4.13; 95% CI 0.167-0.232); SD = sqrt(0.198) = 0.445
  })

  model({
    # Between-occasion variability multiplexed on the dose-occasion index.
    iov_fdepot <- (OCC == 1) * etaiov_fdepot_1 + (OCC == 2) * etaiov_fdepot_2 + (OCC == 3) * etaiov_fdepot_3
    iov_mtt <- (OCC == 1) * etaiov_mtt_1 + (OCC == 2) * etaiov_mtt_2 + (OCC == 3) * etaiov_mtt_3

    # Equation 3: sigmoid maturation of elimination clearance.
    maturation_cl <- AGE^mat_hill / (mat_mf50^mat_hill + AGE^mat_hill)

    allom_cl <- (WT / 54)^e_wt_cl
    allom_v <- (WT / 54)^e_wt_vc

    cl <- exp(lcl) * allom_cl * maturation_cl
    vc <- exp(lvc + etalvc) * allom_v
    q <- exp(lq) * allom_cl
    vp <- exp(lvp + etalvp) * allom_v
    q2 <- exp(lq2 + etalq2) * allom_cl
    vp2 <- exp(lvp2 + etalvp2) * allom_v

    # Two transit compartments with ka = ktr: depot -> transit1 -> transit2
    # -> central has three equal transitions, ktr = 3 / MTT (same convention
    # as the Hoglund 2017 prior).
    mtt <- exp(lmtt + etalmtt + iov_mtt)
    ktr <- 3 / mtt

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(central) <- ktr * transit2 - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Relative bioavailability with the fixed dose-occasion increment.
    f_occ <- 1 + e_doseocc_f * (OCC - 1)
    f(depot) <- f_occ * exp(lfdepot + etalfdepot + iov_fdepot)

    # Dose in mg piperaquine base, volume in L -> mg/L; x1000 -> ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
