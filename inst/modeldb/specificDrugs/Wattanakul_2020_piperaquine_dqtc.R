Wattanakul_2020_piperaquine_dqtc <- function() {
  description <- paste(
    "Sigmoid Emax model of the change from baseline in the heart-rate-",
    "corrected QT interval (DeltaQTcSSB, study-specific correction",
    "QT / RR^0.476) as a function of plasma piperaquine in 1,000 African",
    "patients (mostly children) with uncomplicated falciparum malaria",
    "(Wattanakul 2020, supplementary Table S2): DeltaQTc = (0 + eta) +",
    "47.5 * Cc^1.22 / (Cc^1.22 + EC50^1.22), EC50 = 319 ng/mL at the",
    "study-median age of 7.5 years with a linear age effect on EC50",
    "(2.87 percent per year). The piperaquine population PK model of the",
    "same paper is embedded so the prediction is driven by simulated",
    "concentration; doses are piperaquine base in mg. Sister model files:",
    "modellib('Wattanakul_2020_piperaquine') (population PK) and",
    "modellib('Wattanakul_2020_piperaquine_qtc') (absolute QTc Emax",
    "model).",
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
    "Open Access under CC BY 4.0.",
    "The DeltaQTc parameters are in supplementary Table S2",
    "(AAC.01848-19-s0001.pdf) and the model is supplementary Equation 2;",
    "the embedded PK parameters are in Table 2.",
    sep = " "
  )
  vignette <- "Wattanakul_2020_piperaquine"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "(the observation QTcS is the change from pre-treatment baseline in the study-specific heart-rate-corrected QT interval, DeltaQTcSSB, ms; the driving plasma piperaquine concentration Cc is in ng/mL)"
  )

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
        "Two roles. (1) Embedded PK: sigmoid enzyme maturation on",
        "elimination clearance (Equation 3; MF50 = 0.575 years, Hill = 5.51).",
        "(2) PD: linear effect on EC50, EC50_i = EC50 * (1 + e_age_ec50 *",
        "(AGE - 7.5)) * exp(eta). The paper reports the effect only as",
        "'Effect of age on EC50 (%)' per year without printing the",
        "equation; the linear form centred on the study-median age of 7.5",
        "years (Table 1) is the maintainers' reading, and it reproduces the",
        "Discussion statement that EC50 is lower in young children than in",
        "adults. See the vignette Assumptions and deviations.",
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
    POT = list(
      description = "Serum potassium",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Significant on the QTc baseline in the stepwise search (a 1 mmol/L increase lowered baseline by 1.06 ms, 0.25 percent) but removed from the final model as clinically negligible given the narrow potassium range (IQR 3.70-4.48 mmol/L) (Results, 'Relationship between piperaquine concentration and QTc interval')."
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
    # ==================================================================
    # EMBEDDED PIPERAQUINE PK -- identical to
    # modellib('Wattanakul_2020_piperaquine'), Table 2. The PK residual
    # error is deliberately absent: the PD model was fitted sequentially
    # to individually predicted concentrations, so concentration is a
    # latent driver here.
    # ==================================================================
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


    # ==================================================================
    # CONCENTRATION-DeltaQTcSSB MODEL -- supplementary Table S2 and
    # supplementary Equation 2:
    #   DeltaQTc = (DeltaQTcBaseline + eta) + Emax * Cp^gamma / (Cp^gamma + EC50^gamma) + eps
    # DeltaQTc is the change from each patient's pre-treatment QTcSSB.
    # ==================================================================
    e0 <- fixed(0)
    label("Typical baseline change in QTcSSB (ms)")
    # Table S2: DeltaQTcBaseline = '0 fixed'
    lemax <- log(47.5)
    label("Maximum piperaquine-induced change from baseline in QTcSSB (ms)")
    # Table S2: Emax (ms) = 47.5 (%RSE 9.48; 95% CI 42.1-53.4)
    lec50 <- log(319)
    label("Piperaquine concentration giving half-maximal DeltaQTcSSB at AGE = 7.5 y (ng/mL)")
    # Table S2: EC50 (ng/ml) = 319 (%RSE 16.8; 95% CI 260-409)
    lhill <- log(1.22)
    label("Hill coefficient of the concentration-DeltaQTcSSB Emax function (unitless)")
    # Table S2: gamma = 1.22 (%RSE 6.41; 95% CI 1.11-1.35)
    e_age_ec50 <- 0.0287
    label("Fractional change in EC50 per year of age above 7.5 y (1/year)")
    # Table S2: 'Effect of age on EC50 (%)' = 2.87 (%RSE 35.2; 95% CI 1.31-4.80); linear form centred on the Table 1 median age is the maintainers' reading

    # Table S2 footnote c: baseline IIV additive, absolute SD in ms.
    etae0 ~ 112.36
    # Table S2: DeltaQTcBaseline IIV 10.6 ms (%RSE 15.2; 95% CI 8.05-14.1); variance 10.6^2 = 112.36
    # Table S2 footnote d: exponential IIV as %CV, omega^2 = log(CV^2 + 1).
    etalemax ~ 0.0839881
    # Table S2: Emax IIV 29.6% CV (%RSE 22.9); log(0.296^2 + 1) = 0.0839881
    etalec50 ~ 1.74203
    # Table S2: EC50 IIV 217% CV (%RSE 7.82); log(2.17^2 + 1) = 1.74203

    # Footnote calls sigma an 'additive residual error (variance)' but the
    # row is in ms; read as an SD, as for Table 4.
    addSd <- 12.6
    label("Additive residual SD of DeltaQTcSSB (ms)")
    # Table S2: sigma (ms) = 12.6 (%RSE 3.72; 95% CI 11.7-13.5), read as an SD
  })

  model({
    # ---- embedded piperaquine PK ---------------------------------------
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

    # Plasma piperaquine in ng/mL (dose mg base, volume L), the latent
    # driver of the QTc effect.
    Cc <- 1000 * central / vc

    # ---- sigmoid Emax concentration-DeltaQTcSSB relationship (Table S2 Eq. 2) ----
    ec50 <- exp(lec50 + etalec50) * (1 + e_age_ec50 * (AGE - 7.5))
    emax <- exp(lemax + etalemax)
    hill <- exp(lhill)
    cp <- max(Cc, 0)
    QTcS <- e0 + etae0 + emax * cp^hill / (cp^hill + ec50^hill)
    QTcS ~ add(addSd)
  })
}
