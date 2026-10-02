# Population pharmacokinetic IPD meta-analysis model for oral lumefantrine in
# non-pregnant adults, quantifying drug-drug interactions with commonly used
# antiretroviral and antituberculosis treatment (Francis 2020, Antimicrob
# Agents Chemother 64(5):e02394-19; doi:10.1128/AAC.02394-19).

Francis_2020_lumefantrine <- function() {
  description <- paste(
    "Population PK model for oral lumefantrine from an individual",
    "participant data meta-analysis of 10 studies (793 non-pregnant adults,",
    "6,100 concentrations; HIV-malaria coinfected, malaria-infected,",
    "HIV-infected and healthy volunteers from sub-Saharan Africa and the",
    "USA) treated with artemether-lumefantrine (Francis 2020 Antimicrob",
    "Agents Chemother). Savic transit-compartment absorption (MTT 2.86 h,",
    "7.58 transit compartments, separate first-order ka) into",
    "three-compartment disposition with first-order elimination; allometric",
    "body-weight scaling (exponent 0.75 on clearances, 1 on volumes) centred",
    "at 57 kg. Drug-drug interactions: lopinavir-ritonavir-based ART",
    "(-50.1% CL/F, +67.2% F, -47.6% ka), efavirenz-based ART (+89.9% CL/F)",
    "and rifampicin-based antituberculosis treatment (+142% CL/F).",
    "Study- and dose-occasion-specific relative bioavailability (evening",
    "doses of InterACT/SEACAT = reference), a 2.28-fold dried-blood-spot",
    "matrix factor and an estimated 4.3 h delay of the unobserved 5th dose",
    "in two studies. BSV on CL/F, F and Q1/F, between-visit variability on",
    "CL/F and between-occasion (per-dose) variability on F, MTT and ka;",
    "combined additive and proportional residual error.",
    sep = " "
  )
  reference <- paste(
    "Francis J, Barnes KI, Workman L, Kredo T, Vestergaard LS, Hoglund RM,",
    "Byakika-Kibwika P, Lamorde M, Walimbwa SI, Chijioke-Nwauche I,",
    "Sutherland CJ, Merry C, Scarsi KK, Nyagonde N, Lemnge MM, Khoo SH,",
    "Bygbjerg IC, Parikh S, Aweeka FT, Tarning J, Denti P (2020). An",
    "individual participant data population pharmacokinetic meta-analysis",
    "of drug-drug interactions between lumefantrine and commonly used",
    "antiretroviral treatment. Antimicrobial Agents and Chemotherapy",
    "64(5):e02394-19. doi:10.1128/AAC.02394-19.",
    sep = " "
  )
  vignette <- "Francis_2020_lumefantrine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Between-visit variability on CL/F (Table 3 'BVV') is carried by a
  # second subject-level eta on the CL/F line, not a separate typical value.
  paper_specific_etas <- c("etabvv_cl")

  compartmentData <- list(
    depot = list(analyte = "lumefantrine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lumefantrine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "lumefantrine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "lumefantrine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Total body weight. Allometric scaling of all clearances (exponent",
        "0.75, fixed) and volumes (exponent 1, fixed) centred at 57 kg, the",
        "median body weight of the pooled population (Table 3 footnote c;",
        "Methods 'Population pharmacokinetic modeling'). Fat-free mass and",
        "normal fat mass were tested as alternatives and did not improve the",
        "fit."
      ),
      source_name = "WT"
    ),
    CONMED_LPV = list(
      description = "Lopinavir-ritonavir-based antiretroviral therapy (1 = on LPV/r-based ART)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = coadministered lopinavir-ritonavir-based ART (the ritonavir",
        "boost is part of the regimen and is the CYP3A4 inhibitor the paper",
        "credits for the effect); 0 = not on LPV/r. Fractional effects",
        "(Table 3): CL/F x (1 - 0.501), F x (1 + 0.672), ka x (1 - 0.476).",
        "Together these give the ~3.4-fold higher lumefantrine AUC reported",
        "in the Results."
      ),
      source_name = "LPV-r"
    ),
    CONMED_EFV = list(
      description = "Efavirenz-based antiretroviral therapy (1 = on EFV-based ART)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = coadministered efavirenz-based ART; 0 = not on efavirenz",
        "(reference pools ART-naive, nevirapine- and dolutegravir-based ART,",
        "which had no significant effect). Fractional effect on CL/F",
        "x (1 + 0.899) (Table 3), i.e. 47% lower AUC (Results)."
      ),
      source_name = "EFV"
    ),
    CONMED_RIF = list(
      description = "Rifampicin-based antituberculosis treatment (1 = on rifampicin)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = on rifampicin-based antituberculosis treatment at the time of",
        "artemether-lumefantrine dosing (chronic, fully induced); 0 = not.",
        "Fractional effect on CL/F x (1 + 1.42) (Table 3), i.e. 59% lower",
        "AUC (Results). Four of the 13 rifampicin-treated participants were",
        "also on efavirenz; a trend toward an even higher clearance with",
        "the combination was not statistically significant and was not",
        "retained, so the two effects combine multiplicatively here as in",
        "the final model."
      ),
      source_name = "RIF"
    ),
    STUDY_SEACAT = list(
      description = "SEACAT 2.4.1 or 2.4.2 study participant (South Africa)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = participant of the SEACAT 2.4.1 or SEACAT 2.4.2 study. Selects",
        "the dose-occasion bioavailability effects of Table 3: the first",
        "(morning) dose, OCC = 1, has F x (1 - 0.486); every subsequent",
        "morning dose, OCC = 3, 5, ..., has F x (1 - 0.772); evening doses",
        "(OCC = 2, 4, 6) are at the reference. The SEACAT protocol doses at",
        "0, 8 and 24 h and every 12 h thereafter starting in the morning,",
        "so odd OCC values are the morning doses (Table 2)."
      ),
      source_name = "study (SEACAT)"
    ),
    STUDY_UGANDA = list(
      description = "Uganda study 1, 2, 3 or 4 participant",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = participant of any of the four Ugandan studies (doses taken",
        "with a standard Ugandan breakfast). All doses have F x (1 - 0.269)",
        "(Table 3; Results '(a) Diurnal variation': 'The value of relative",
        "bioavailability in the Ugandan studies was similar and was found",
        "to be 26.9% lower than the reference')."
      ),
      source_name = "study (Uganda studies)"
    ),
    STUDY_NIGERIA1 = list(
      description = "Nigeria study 1 participant",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = participant of Nigeria study 1 (Parikh et al.). Two effects:",
        "the 6th (morning, observed) dose, OCC = 6, has F x (1 - 0.608)",
        "(Table 3), and the unobserved 5th dose, OCC = 5, is delayed by",
        "tlag_unobs = 4.30 h relative to its imputed time (Table 3 footnote",
        "d)."
      ),
      source_name = "study (Nigeria study 1)"
    ),
    STUDY_USHV = list(
      description = "US healthy-volunteer study (German et al. 2009) participant",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = participant of the U.S. healthy-volunteer study. Bioavailability",
        "was not different from the reference, but the unobserved 5th dose,",
        "OCC = 5, is delayed by tlag_unobs = 4.30 h relative to its imputed",
        "time (Table 3 footnote d)."
      ),
      source_name = "study (U.S. healthy volunteer study)"
    ),
    OCC = list(
      description = "Dose number within the artemether-lumefantrine course (1 = first dose)",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Integer dose index, set on every dose record and carried forward on",
        "the observation records that follow it (OCC = k from the k-th dose",
        "until the (k+1)-th dose). It plays two roles: (1) the per-dose",
        "between-occasion variability on F, MTT and ka (each dose is one",
        "occasion; 12 occasion slots are provided so that the paper's",
        "extended 4-, 5- and 6-day regimens can be simulated; OCC values",
        "outside 1-12 carry no BOV), and (2) selecting the study-specific",
        "dose-occasion effects (see STUDY_SEACAT, STUDY_NIGERIA1,",
        "STUDY_USHV). The source does not state the occasion definition;",
        "one occasion per dose is assumed (see vignette Assumptions)."
      ),
      source_name = "dose number"
    ),
    SAMPLE_DBS = list(
      description = "Dried-blood-spot (capillary whole blood) sample indicator (per observation)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = concentration measured in a capillary dried blood spot (Nigeria",
        "study 2, the only DBS study); 0 = venous plasma (all other",
        "studies). The predicted concentration is multiplied by the",
        "estimated matrix scaling factor 2.28 (Table 3 footnote e)."
      ),
      source_name = "sample matrix"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_NVP = list(
      description = "Nevirapine-based antiretroviral therapy",
      units = "(binary)",
      type = "binary",
      notes = "Tested; discordant per-study trends and no significant effect in the pooled data after adjusting for other factors (Results '(ii) Drug-drug interactions'). Not retained."
    ),
    CONMED_DOLUTEGRAVIR = list(
      description = "Dolutegravir-based antiretroviral therapy",
      units = "(binary)",
      type = "binary",
      notes = "Tested; 'Dolutegravir-based ART did not alter lumefantrine exposure' (Results). Not retained."
    ),
    HIV_POS = list(
      description = "HIV infection",
      units = "(binary)",
      type = "binary",
      notes = "Tested; a trend toward higher CL/F in ART-naive HIV+ SEACAT participants was not consistent across studies and was not retained (Results '(iii) HIV and malaria disease effects')."
    ),
    DIS_MALARIA_ACUTE = list(
      description = "Malaria infection at the time of dosing",
      units = "(binary)",
      type = "binary",
      notes = "Tested; no significant consistent difference in PK parameters attributable to malaria (Results '(iii)'). Not retained."
    ),
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      notes = "Tested (with normal fat mass) as an alternative body-size descriptor to total body weight; did not improve the fit (Results '(i)')."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 793L,
    n_studies = 10L,
    n_observations = 6100L,
    age_range = "adults; study-group median ages 26-43.5 years, overall range about 14-70 years (Table 1)",
    weight_range = "median 57 kg (allometric reference); study-group ranges spanning 25-104.7 kg (Table 1)",
    disease_state = paste(
      "Non-pregnant adults: 41% HIV-malaria coinfected, 36% malaria-infected",
      "(HIV-uninfected), 20% HIV-infected (malaria-uninfected) and 3%",
      "healthy volunteers (Abstract). Concomitant treatment: none,",
      "nevirapine-, efavirenz-, lopinavir-ritonavir- or dolutegravir-based",
      "ART, and rifampicin-based antituberculosis treatment (13",
      "participants)."
    ),
    dose_range = paste(
      "Artemether-lumefantrine (Coartem, 20/120 mg tablets) 4 tablets",
      "(480 mg lumefantrine) per dose, twice daily for 3 days (6 doses at",
      "0, 8, 24, 36, 48 and 60 h), or a single 480 mg dose in two study",
      "phases (Table 2)."
    ),
    regions = "South Africa, Tanzania, Uganda, Nigeria, USA",
    notes = paste(
      "WWARN individual participant data meta-analysis. Studies: SEACAT",
      "2.4.1 and 2.4.2 (South Africa), InterACT (Tanzania), Uganda studies",
      "1-4, Nigeria studies 1-2 and a U.S. healthy-volunteer study. 341",
      "(5.59%) samples were below the LLOQ and handled with the M6 method.",
      "Nigeria study 2 measured capillary dried blood spots; all other",
      "studies measured venous plasma."
    )
  )

  ini({
    # ---- Disposition (typical 57-kg adult; Table 3 footnote c) --------------
    lcl <- log(3.28); label("Apparent clearance CL/F (L/h)") # Table 3: CL = 3.28 L/h (95% CI 3.14-3.46)
    lvc <- log(60); label("Apparent central volume Vc/F (L)") # Table 3: central volume = 60 L (95% CI 56.3-63.9)
    lq <- log(0.63); label("Apparent intercompartmental clearance to first peripheral Q1/F (L/h)") # Table 3: Q1 = 0.63 L/h (95% CI 0.60-0.67)
    lvp <- log(182); label("Apparent first peripheral volume Vp1/F (L)") # Table 3: V first peripheral = 182 L (95% CI 171-195)
    lq2 <- log(1.55); label("Apparent intercompartmental clearance to second peripheral Q2/F (L/h)") # Table 3: Q2 = 1.55 L/h (95% CI 1.43-1.72)
    lvp2 <- log(39.1); label("Apparent second peripheral volume Vp2/F (L)") # Table 3: V second peripheral = 39.1 L (95% CI 37.1-41.2)

    # ---- Allometric exponents (fixed; Methods) ------------------------------
    e_wt_cl <- fixed(0.75); label("Allometric exponent on clearances (unitless)") # Methods: exponent fixed to 0.75 for clearance parameters
    e_wt_vc <- fixed(1); label("Allometric exponent on volumes (unitless)") # Methods: exponent fixed to 1 for volumes of distribution

    # ---- Absorption ----------------------------------------------------------
    lmtt <- log(2.86); label("Mean absorption transit time MTT (h)") # Table 3: MTT = 2.86 h (95% CI 2.74-2.94)
    lnn <- log(7.58); label("Number of hypothetical transit compartments NN (unitless)") # Table 3: NN = 7.58 (95% CI 7.06-8.10)
    lka <- log(0.727); label("First-order absorption rate constant ka (1/h)") # Table 3: Ka = 0.727 1/h (95% CI 0.62-0.83)
    lfdepot <- fixed(log(1)); label("Relative oral bioavailability F (reference: InterACT/SEACAT evening doses)") # Table 3: F = 1 FIXED

    # ---- Drug-drug interactions (fractional effects; Table 3) ---------------
    e_conmed_lpv_cl <- -0.501; label("Fractional effect of LPV/r-based ART on CL/F (unitless)") # Table 3: Lopinavir-ritonavir on CL = -50.1% (95% CI -53.0 to -46.4)
    e_conmed_lpv_fdepot <- 0.672; label("Fractional effect of LPV/r-based ART on F (unitless)") # Table 3: Lopinavir-ritonavir on F = 67.2% (95% CI 49.1-88.9)
    e_conmed_lpv_ka <- -0.476; label("Fractional effect of LPV/r-based ART on ka (unitless)") # Table 3: Lopinavir-ritonavir on Ka = -47.6% (95% CI -56.5 to -37.4)
    e_conmed_efv_cl <- 0.899; label("Fractional effect of efavirenz-based ART on CL/F (unitless)") # Table 3: Efavirenz on CL = 89.9% (95% CI 81.1-99.7)
    e_conmed_rif_cl <- 1.42; label("Fractional effect of rifampicin-based TB treatment on CL/F (unitless)") # Table 3: Rifampin-based TB treatment on CL = 142% (95% CI 111-180)

    # ---- Study and dose-occasion effects on F (fractional; Table 3) ---------
    e_study_seacat_first_fdepot <- -0.486; label("Fractional effect on F of the first (morning) dose in SEACAT (unitless)") # Table 3: First dose in SECAT on F = -48.6% (95% CI -54.9 to -41.7)
    e_study_seacat_am_fdepot <- -0.772; label("Fractional effect on F of consecutive morning doses in SEACAT (unitless)") # Table 3: Consecutive morning doses in SEACAT on F = -77.2% (95% CI -80.7 to -73.8)
    e_study_uganda_fdepot <- -0.269; label("Fractional effect on F in the Uganda studies (unitless)") # Table 3: Uganda studies on F = -26.9% (95% CI -32.3 to -20.7)
    e_study_nigeria1_fdepot <- -0.608; label("Fractional effect on F of the 6th (morning) dose in Nigeria study 1 (unitless)") # Table 3: Nigeria study 1 on F = -60.8% (95% CI -73.2 to -47.3); Results text says 60.1%, see vignette
    tlag_unobs <- 4.30; label("Delay of the unobserved 5th dose in Nigeria study 1 and the US healthy-volunteer study (h)") # Table 3: Delay for unobserved dose = 4.30 h (95% CI 2.84-5.73)
    e_sample_dbs_cc <- 2.28; label("Dried-blood-spot to plasma concentration scaling factor (fold)") # Table 3: Scaling factor for DBS concn = 2.28-fold (95% CI 2.05-2.55)

    # ---- Between-subject / between-visit variability ------------------------
    # Table 3 reports 'approximate CV%' (footnote b); for log-normal etas the
    # approximate CV is sqrt(omega^2), so omega^2 = (CV/100)^2.
    etalcl ~ 0.043264 # Table 3: BSV on CL 20.8% (approximate CV) -> omega^2 = 0.208^2
    etabvv_cl ~ 0.023716 # Table 3: BVV on CL 15.4% (approximate CV) -> omega^2 = 0.154^2
    etalfdepot ~ 0.091204 # Table 3: BSV on F 30.2% (approximate CV) -> omega^2 = 0.302^2
    etalq ~ 0.084681 # Table 3: BSV on Q1 29.1% (approximate CV) -> omega^2 = 0.291^2

    # ---- Between-occasion (per-dose) variability ----------------------------
    etaiov_fdepot_1 ~ 0.319225  # Table 3: BOV on F 56.5% (approximate CV) -> omega^2 = (56.5/100)^2
    etaiov_fdepot_2 ~ fixed(0.319225)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_fdepot_3 ~ fixed(0.319225)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_fdepot_4 ~ fixed(0.319225)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_fdepot_5 ~ fixed(0.319225)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_fdepot_6 ~ fixed(0.319225)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_fdepot_7 ~ fixed(0.319225)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_fdepot_8 ~ fixed(0.319225)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_fdepot_9 ~ fixed(0.319225)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_fdepot_10 ~ fixed(0.319225)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_fdepot_11 ~ fixed(0.319225)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_fdepot_12 ~ fixed(0.319225)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_mtt_1 ~ 0.101761  # Table 3: BOV on MTT 31.9% (approximate CV) -> omega^2 = (31.9/100)^2
    etaiov_mtt_2 ~ fixed(0.101761)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_mtt_3 ~ fixed(0.101761)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_mtt_4 ~ fixed(0.101761)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_mtt_5 ~ fixed(0.101761)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_mtt_6 ~ fixed(0.101761)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_mtt_7 ~ fixed(0.101761)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_mtt_8 ~ fixed(0.101761)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_mtt_9 ~ fixed(0.101761)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_mtt_10 ~ fixed(0.101761)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_mtt_11 ~ fixed(0.101761)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_mtt_12 ~ fixed(0.101761)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_ka_1 ~ 0.544644  # Table 3: BOV on ka 73.8% (approximate CV) -> omega^2 = (73.8/100)^2
    etaiov_ka_2 ~ fixed(0.544644)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_ka_3 ~ fixed(0.544644)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_ka_4 ~ fixed(0.544644)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_ka_5 ~ fixed(0.544644)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_ka_6 ~ fixed(0.544644)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_ka_7 ~ fixed(0.544644)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_ka_8 ~ fixed(0.544644)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_ka_9 ~ fixed(0.544644)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_ka_10 ~ fixed(0.544644)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_ka_11 ~ fixed(0.544644)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)
    etaiov_ka_12 ~ fixed(0.544644)  # same BOV variance as occasion 1 (NONMEM BLOCK(1) SAME)

    # ---- Residual error (Table 3) ---------------------------------------------
    addSd <- 32.9; label("Additive residual error (ng/mL)") # Table 3: Additive error = 32.9 ng/mL (95% CI 32.2-33.5)
    propSd <- 0.142; label("Proportional residual error (fraction)") # Table 3: Proportional error = 14.2% (95% CI 13.9-14.4)
  })

  model({
    # Occasion (dose-number) indicators.
    oc1 <- OCC == 1
    oc2 <- OCC == 2
    oc3 <- OCC == 3
    oc4 <- OCC == 4
    oc5 <- OCC == 5
    oc6 <- OCC == 6
    oc7 <- OCC == 7
    oc8 <- OCC == 8
    oc9 <- OCC == 9
    oc10 <- OCC == 10
    oc11 <- OCC == 11
    oc12 <- OCC == 12

    iov_fdepot <- oc1 * etaiov_fdepot_1 +
      oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 +
      oc4 * etaiov_fdepot_4 +
      oc5 * etaiov_fdepot_5 +
      oc6 * etaiov_fdepot_6 +
      oc7 * etaiov_fdepot_7 +
      oc8 * etaiov_fdepot_8 +
      oc9 * etaiov_fdepot_9 +
      oc10 * etaiov_fdepot_10 +
      oc11 * etaiov_fdepot_11 +
      oc12 * etaiov_fdepot_12
    iov_mtt <- oc1 * etaiov_mtt_1 +
      oc2 * etaiov_mtt_2 +
      oc3 * etaiov_mtt_3 +
      oc4 * etaiov_mtt_4 +
      oc5 * etaiov_mtt_5 +
      oc6 * etaiov_mtt_6 +
      oc7 * etaiov_mtt_7 +
      oc8 * etaiov_mtt_8 +
      oc9 * etaiov_mtt_9 +
      oc10 * etaiov_mtt_10 +
      oc11 * etaiov_mtt_11 +
      oc12 * etaiov_mtt_12
    iov_ka <- oc1 * etaiov_ka_1 +
      oc2 * etaiov_ka_2 +
      oc3 * etaiov_ka_3 +
      oc4 * etaiov_ka_4 +
      oc5 * etaiov_ka_5 +
      oc6 * etaiov_ka_6 +
      oc7 * etaiov_ka_7 +
      oc8 * etaiov_ka_8 +
      oc9 * etaiov_ka_9 +
      oc10 * etaiov_ka_10 +
      oc11 * etaiov_ka_11 +
      oc12 * etaiov_ka_12

    # Disposition: allometric scaling at 57 kg and fractional DDI effects on
    # CL/F. Between-visit variability on CL/F is a single subject-level draw
    # (one visit); redraw etabvv_cl for a second visit of the same subject.
    cl_bsv <- exp(lcl + etalcl)
    cl <- cl_bsv * exp(etabvv_cl) * (WT / 57)^e_wt_cl *
      (1 + e_conmed_lpv_cl * CONMED_LPV) *
      (1 + e_conmed_efv_cl * CONMED_EFV) *
      (1 + e_conmed_rif_cl * CONMED_RIF)
    vc <- exp(lvc) * (WT / 57)^e_wt_vc
    q <- exp(lq + etalq) * (WT / 57)^e_wt_cl
    vp <- exp(lvp) * (WT / 57)^e_wt_vc
    q2 <- exp(lq2) * (WT / 57)^e_wt_cl
    vp2 <- exp(lvp2) * (WT / 57)^e_wt_vc

    # Absorption. Savic 2007 transit model: ktr = (NN + 1) / MTT, the
    # separately estimated ka empties the absorption compartment (depot).
    mtt <- exp(lmtt + iov_mtt)
    nn <- exp(lnn)
    ktr <- (nn + 1) / mtt
    ka <- exp(lka + iov_ka) * (1 + e_conmed_lpv_ka * CONMED_LPV)

    # Relative bioavailability of the current dose: reference = evening doses
    # of InterACT and SEACAT; study- and dose-occasion-specific fractional
    # effects (Results '(a) Diurnal variation'; Table 3).
    seacat_first <- STUDY_SEACAT * oc1
    seacat_am <- STUDY_SEACAT * (oc3 + oc5 + oc7 + oc9 + oc11)
    nigeria1_d6 <- STUDY_NIGERIA1 * oc6
    fdepot <- exp(lfdepot + etalfdepot + iov_fdepot) *
      (1 + e_conmed_lpv_fdepot * CONMED_LPV) *
      (1 + e_study_seacat_first_fdepot * seacat_first) *
      (1 + e_study_seacat_am_fdepot * seacat_am) *
      (1 + e_study_uganda_fdepot * STUDY_UGANDA) *
      (1 + e_study_nigeria1_fdepot * nigeria1_d6)

    # Estimated delay of the unobserved 5th dose, whose time was imputed as
    # exactly 12 h before the observed 6th dose (Results '(c) Dosing time').
    dly <- tlag_unobs * (STUDY_NIGERIA1 + STUDY_USHV) * oc5

    # Transit input (Savic gamma density) written out explicitly rather than
    # through rxode2's transit() macro, which rescales its internal dose lookup
    # by bioavailability and delivers zero dose under f(depot) <- 0 in nlmixr2
    # UI form. podo(depot) is the un-scaled amount of the most recent dose, so
    # fdepot multiplies it here. MTT (2.86 h) is far shorter than the 8-12 h
    # dosing interval, so the most-recent-dose semantics of tad()/podo() lose
    # no mass (asserted in the vignette).
    tdose <- tad(depot) - dly
    if (tdose > 0) {
      transit_in <- fdepot * exp(log(podo(depot)) + log(ktr) + nn * log(ktr * tdose) -
        ktr * tdose - lgamma(nn + 1))
    } else {
      transit_in <- 0
    }

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(depot) <- transit_in - ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # The dose record only arms tad()/podo(); the gamma input delivers it.
    f(depot) <- 0

    # Venous plasma lumefantrine (mg/L x 1000 = ng/mL); dried-blood-spot
    # samples are scaled by the estimated 2.28-fold matrix factor.
    Cc <- 1000 * central / vc * (1 + (e_sample_dbs_cc - 1) * SAMPLE_DBS)
    Cc ~ add(addSd) + prop(propSd)
  })
}
