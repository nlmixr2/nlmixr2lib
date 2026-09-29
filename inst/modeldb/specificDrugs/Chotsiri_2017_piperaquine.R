Chotsiri_2017_piperaquine <- function() {
  description <- paste(
    "Three-compartment population pharmacokinetic model for piperaquine in 16",
    "healthy Thai adult volunteers, from Chotsiri 2017. A single oral dose of",
    "three co-formulated dihydroartemisinin-piperaquine tablets was given on",
    "two separate occasions, with and without a concomitant single low dose of",
    "primaquine, in an open-label randomized three-way crossover. Absorption is",
    "a two-transit-compartment chain whose transit rate constant and final",
    "absorption rate constant were not distinguishable and were set equal;",
    "body weight enters as a fixed allometric function on all clearance and",
    "volume terms centred on the study-median 64 kg. Between-subject",
    "variability is carried on relative bioavailability, clearance, central",
    "volume and the second inter-compartmental clearance, with between-occasion",
    "variability on relative bioavailability and mean transit time.",
    "Primaquine coadministration was",
    "screened on every pharmacokinetic parameter and produced no clinically",
    "relevant interaction, so no primaquine covariate appears in the final",
    "model. Doses are expressed as piperaquine BASE. Sister model files from",
    "the same paper: modellib('Chotsiri_2017_dihydroartemisinin')",
    "(dihydroartemisinin population PK) and",
    "modellib('Chotsiri_2017_piperaquine_qtc') (this same PK model driving the",
    "linear concentration-QTc prolongation model).",
    sep = " "
  )
  reference <- paste(
    "Chotsiri P, Wattanakul T, Hoglund RM, Hanboonkunupakarn B,",
    "Pukrittayakamee S, Blessborn D, Jittamala P, White NJ, Day NPJ,",
    "Tarning J.",
    "Population pharmacokinetics and electrocardiographic effects of",
    "dihydroartemisinin-piperaquine in healthy volunteers.",
    "Br J Clin Pharmacol. 2017;83(12):2752-2766.",
    "doi:10.1111/bcp.13372. PMC5698590.",
    "Open Access under CC BY-NC 4.0.",
    "Parameter estimates are in Table 2 ('Pharmacokinetic parameters of",
    "piperaquine'); secondary exposure parameters used for validation are in",
    "Table 3. Sister model files from the same paper:",
    "modellib('Chotsiri_2017_dihydroartemisinin') and",
    "modellib('Chotsiri_2017_piperaquine_qtc').",
    sep = " "
  )
  vignette <- "Chotsiri_2017_dihydroartemisinin_piperaquine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Fixed allometric function on all clearance parameters (CL/F, Q1/F,",
        "Q2/F; exponent 0.75) and all volume parameters (Vc/F, Vp1/F, Vp2/F;",
        "exponent 1), centred on the study-median body weight of 64 kg",
        "(Methods, 'Population pharmacokinetic analysis', Equations 3 and 4).",
        "Unlike the dihydroartemisinin model of the same paper it improved the",
        "piperaquine fit (dOFV = -5.95; Results, 'Population pharmacokinetic",
        "properties of piperaquine'). Set WT = 64 to recover the tabulated",
        "typical values.",
        sep = " "
      ),
      source_name = "BW"
    ),
    OCC = list(
      description = paste(
        "Integer-valued occasion index identifying which of the two",
        "dihydroartemisinin-piperaquine dosing occasions a dose belongs to,",
        "for between-occasion variability on relative bioavailability and",
        "the mean transit time.",
        sep = " "
      ),
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Each volunteer received dihydroartemisinin-piperaquine twice -- once",
        "alone and once with primaquine -- in randomised order, separated by",
        "an 8-week washout (Methods, 'Study design'). Values 1 and 2 are",
        "multiplexed inside model() onto etaiov_fdepot_1 / etaiov_fdepot_2 and",
        "etaiov_mtt_1 / etaiov_mtt_2, each pair sharing one variance.",
        "OCC = 0, or any value outside 1..2, yields the occasion-free",
        "typical value.",
        sep = " "
      ),
      source_name = "OCC"
    )
  )

  # Covariates the source screened and did not retain in the final model.
  # Documentation only -- none is referenced in model().
  covariatesDataExcluded <- list(
    CONMED_PRIMAQUINE = list(
      description = "Concomitant single low-dose primaquine administration indicator.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (dihydroartemisinin-piperaquine alone)",
      notes = paste(
        "Screened two ways and not retained. In the stepwise approach it was",
        "not significant. The full covariate approach (a categorical primaquine",
        "effect on every pharmacokinetic parameter simultaneously, bootstrapped",
        "n = 1000) did show 'a median 37.3% (95% CI -67.6%, 33.7%) decrease in",
        "central volume of distribution and a median 26.8% (95% CI -21.2%,",
        "62.5%) increase in mean transit absorption time', but 'the 95% CI for",
        "these covariate effects included a zero effect, so a lack of effect",
        "could not be excluded' (Results, 'Drug-drug interactions'; Figure 3B).",
        "Those two medians are reported for the exploratory full-covariate run,",
        "not for the final model, and are therefore not encoded here.",
        sep = " "
      )
    ),
    AGE = list(
      description = "Subject age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened by stepwise forward inclusion / backward elimination and not retained: 'No other covariates were significant in the stepwise covariate approach' (Results). No coefficient is reported."
    ),
    AST = list(
      description = "Aspartate aminotransferase.",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate approach (Methods) and not retained. No coefficient is reported."
    ),
    ALT = list(
      description = "Alanine aminotransferase.",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate approach (Methods) and not retained. No coefficient is reported."
    ),
    ALP = list(
      description = "Alkaline phosphatase.",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate approach (Methods) and not retained. No coefficient is reported."
    ),
    HGB = list(
      description = "Haemoglobin.",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate approach (Methods) and not retained. No coefficient is reported."
    ),
    BUN = list(
      description = "Blood urea nitrogen.",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate approach (Methods) and not retained. No coefficient is reported."
    ),
    SCR = list(
      description = "Serum creatinine.",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate approach (Methods) and not retained. No coefficient is reported."
    ),
    ALB = list(
      description = "Serum albumin.",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate approach (Methods) and not retained. Table 1 reports the cohort median as 4.25 g/dL. No coefficient is reported."
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Explicitly not evaluated: 'Gender was not evaluated as a covariate owing to the substantial imbalance between male and female subjects' (Methods, 'Population pharmacokinetic analysis')."
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
    n_subjects = 16L,
    n_studies = 1L,
    age_range = "22-53 years",
    age_median = "40 years",
    weight_range = "54.0-71.4 kg",
    weight_median = "64.1 kg",
    sex_female_pct = 68.75,
    race_ethnicity = c(Asian = 100),
    disease_state = "Healthy volunteers. Participants with malaria or glucose-6-phosphate dehydrogenase deficiency, and pregnant or lactating women, were excluded.",
    dose_range = "Single oral dose of three co-formulated dihydroartemisinin-piperaquine tablets (40 mg dihydroartemisinin + 320 mg piperaquine phosphate each), i.e. 960 mg piperaquine phosphate = 514 mg piperaquine base, given 30 min after a light meal (~200 kcal, 8 g fat).",
    regions = "Bangkok, Thailand",
    n_observations = 623L,
    notes = paste(
      "Open-label, randomized, three-way crossover (NCT01525511, TMEC 12-004,",
      "OXTREC 58-11) conducted 18 June to 2 November 2012 at the Faculty of",
      "Tropical Medicine, Mahidol University. Every volunteer received",
      "primaquine alone first (1-week washout), then dihydroartemisinin-",
      "piperaquine alone and dihydroartemisinin-piperaquine plus primaquine in",
      "random order with an 8-week washout. Plasma sampling at 0, 0.25, 0.5, 1,",
      "1.5, 2, 3, 4, 6, 8, 10, 12 and 24 h post dose plus days 3, 4, 7, 11, 15,",
      "22 and 36; LC-MS/MS with a lower limit of quantification of 1.50 ng/mL.",
      "Only 2.3% of piperaquine samples were below the limit of quantification",
      "and these were omitted. Estimated in NONMEM 7.3 with FOCE-I. Baseline",
      "demographics are in Table 1. Note the Discussion states 'three males and",
      "13 female' whereas the Methods state 'five males out of 16 subjects';",
      "sex_female_pct here uses the Methods figure (11 of 16). See the vignette",
      "Errata.",
      sep = " "
    )
  )

  ini({
    # ================================================================
    # Structural parameters -- Chotsiri 2017 Table 2, 'Pharmacokinetic
    # parameters of piperaquine'. Table 2 footnote a: 'Parameter
    # estimates are based on the typical individual in the population
    # with a body weight of 64 kg.' All disposition parameters are
    # apparent (relative to F, which is fixed to unity).
    #
    # DOSE UNITS. The tablets are labelled in piperaquine phosphate but
    # the assay measures piperaquine base, and the fitted CL/F is on the
    # base scale: AUC_inf = Dose / (CL/F) = 514 / 27.4 = 18 800 ng h/mL
    # sits between the paper's own two published median AUCs of 17 700
    # and 19 600 h ng/mL (Table 3), whereas the unconverted 960 mg would
    # give 35 000. Dose this model in mg of piperaquine BASE
    # (320 mg piperaquine phosphate tetrahydrate, MW 999.56, carries
    # 320 * 535.51 / 999.56 = 171.4 mg base). See vignette Errata.
    # ================================================================

    # ---- Absorption -------------------------------------------------
    # Two transit compartments (Results, 'Population pharmacokinetic
    # properties of piperaquine': 'A transit compartment absorption
    # model with two transit compartments was superior to all other
    # models evaluated'). Unlike the dihydroartemisinin model of the
    # same paper, the transit rate constant and the rate constant out of
    # the last transit compartment were NOT distinguishable and were set
    # equal ('There was no significant change in model fit when the
    # transit rate between transit compartments and the absorption rate
    # from the last transit compartment to the central compartment were
    # set to be equal, dOFV = 0.564'), which is why Table 2 carries no
    # k_a row for piperaquine.
    #
    # With ka = ktr all three transfers (depot -> transit1 -> transit2
    # -> central) share one rate, so MTT spans n + 1 = 3 transfers and
    # ktr = 3 / MTT (Savic 2007), the same convention the sibling
    # Chotsiri_2019_piperaquine.R and Hoglund_2017_piperaquine.R use.
    # Confirmed against the paper's own secondary parameters: it
    # reproduces the published median Tmax of 3.76-3.98 h (4.14 h
    # simulated for the typical 64-kg subject) whereas ktr = 2 / MTT
    # gives 5.74 h. See the vignette Source trace.
    lmtt <- log(3.13)
    label("Mean transit time through the two-transit-compartment absorption chain MTT (h)")
    # Chotsiri 2017 Table 2: MTT = 3.13 h (%RSE 9.42; bootstrap 95% CI 2.66-3.84)

    # ---- Disposition ------------------------------------------------
    # Three-compartment disposition (Results: a three-compartment model
    # beat a two-compartment model by dOFV = -297, and a fourth
    # compartment added only dOFV = -0.500).
    lcl <- log(27.4)
    label("Apparent elimination clearance CL/F at WT = 64 kg (L/h)")
    # Chotsiri 2017 Table 2: CL/F = 27.4 L/h (%RSE 5.50; bootstrap 95% CI 24.6-30.4)

    lvc <- log(751)
    label("Apparent central volume of distribution Vc/F at WT = 64 kg (L)")
    # Chotsiri 2017 Table 2: V_C/F = 751 L (%RSE 23.5; bootstrap 95% CI 470-1160)

    lq <- log(206)
    label("Apparent inter-compartmental clearance to peripheral1, Q1/F at WT = 64 kg (L/h)")
    # Chotsiri 2017 Table 2: Q_P1/F = 206 L/h (%RSE 9.56; bootstrap 95% CI 166-242)

    lvp <- log(1900)
    label("Apparent first peripheral volume of distribution Vp1/F at WT = 64 kg (L)")
    # Chotsiri 2017 Table 2: V_P1/F = 1900 L (%RSE 8.23; bootstrap 95% CI 1660-2260)

    lq2 <- log(71.5)
    label("Apparent inter-compartmental clearance to peripheral2, Q2/F at WT = 64 kg (L/h)")
    # Chotsiri 2017 Table 2: Q_P2/F = 71.5 L/h (%RSE 9.01; bootstrap 95% CI 58.5-84.4)

    lvp2 <- log(13500)
    label("Apparent second peripheral volume of distribution Vp2/F at WT = 64 kg (L)")
    # Chotsiri 2017 Table 2: V_P2/F = 13 500 L (%RSE 8.95; bootstrap 95% CI 11 400-16 000)

    # ---- Relative bioavailability -----------------------------------
    lfdepot <- fixed(log(1))
    label("Relative bioavailability F for the typical subject (unitless)")
    # Chotsiri 2017 Table 2: F (%) = '100 Fixed'. Methods: 'Between-subject
    # and between-occasion variability was also evaluated on the relative
    # bioavailability, fixed to unity for the population'.

    # ---- Allometric scaling -----------------------------------------
    # Methods, 'Population pharmacokinetic analysis': 'Body weight was
    # introduced into the pharmacokinetic model as a fixed allometric
    # function on all volume, clearance and distribution parameters,
    # centred on the median body weight (64 kg) of the study population'
    # (Equations 3 and 4). Fixed, not estimated, and reported without
    # uncertainty.
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL/F, Q1/F and Q2/F, referenced to 64 kg (unitless)")
    # Chotsiri 2017 Methods, Equation 3: clearance exponent 3/4

    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on Vc/F, Vp1/F and Vp2/F, referenced to 64 kg (unitless)")
    # Chotsiri 2017 Methods, Equation 4: volume exponent 1

    # ================================================================
    # Variability. Table 2's variability column is headed '%CV of
    # BSV/BOV' and footnote a gives the transform explicitly: the %CV
    # is '100 * sqrt(exp(estimate) - 1)'. Converting back to the
    # internal variance scale is therefore omega^2 = log(CV^2 + 1).
    # For piperaquine only the MTT row carries the '*' that marks
    # between-OCCASION variability; the F, CL/F, Vc/F and Q_P2/F rows
    # are between-SUBJECT variability (Results: 'substantial
    # between-subject and between-occasion variability in the
    # absorption of piperaquine, with additional between-subject
    # variability in the elimination clearance, the inter-compartmental
    # clearance and the central volume of distribution').
    # ================================================================

    etalfdepot ~ 0.0315384
    # Chotsiri 2017 Table 2: F BSV = 17.9% CV (%RSE 34.0; 95% CI 0.178%-26.1%); log(0.179^2 + 1) = 0.0315384

    etalcl ~ 0.0118110
    # Chotsiri 2017 Table 2: CL/F BSV = 10.9% CV (%RSE 37.2; 95% CI 0.109%-15.72%); log(0.109^2 + 1) = 0.0118110

    etalvc ~ 0.1653246
    # Chotsiri 2017 Table 2: V_C/F BSV = 42.4% CV (%RSE 40.9; 95% CI 0.406%-62.9%); log(0.424^2 + 1) = 0.1653246

    etalq2 ~ 0.0564569
    # Chotsiri 2017 Table 2: Q_P2/F BSV = 24.1% CV (%RSE 36.3; 95% CI 0.203%-37.3%); log(0.241^2 + 1) = 0.0564569

    # Between-occasion variability, two occasions (Methods, 'Study
    # design'). Each occasion's eta shares the single published
    # variance; OCC selects which one is active inside model().
    #
    # NOTE on the F row. Piperaquine is the one parameter in Table 2
    # carrying BOTH a between-subject and a between-occasion entry: the
    # printed F row has two variability lines, '17.9% (34.0%)' with no
    # asterisk and '19.1% (13.3%)*' with one. That second line is lost
    # in text extractions that collapse the multi-line table cell, and
    # is the direct tabular evidence for the Results sentence
    # 'substantial between-subject AND between-occasion variability in
    # the absorption of piperaquine'.
    etaiov_fdepot_1 ~ 0.0358313
    # Chotsiri 2017 Table 2: F BOV = 19.1% CV (%RSE 13.3; 95% CI 13.5%-23.3%); log(0.191^2 + 1) = 0.0358313
    etaiov_fdepot_2 ~ fixed(0.0358313)
    # Same published BOV variance as occasion 1

    etaiov_mtt_1 ~ 0.0986537
    # Chotsiri 2017 Table 2: MTT BOV = 32.2% CV (%RSE 13.4; 95% CI 21.1%-37.8%); log(0.322^2 + 1) = 0.0986537
    etaiov_mtt_2 ~ fixed(0.0986537)
    # Same published BOV variance as occasion 1

    # ================================================================
    # Residual error. Methods, 'Population pharmacokinetic analysis':
    # 'Residual unexplained variability was modelled as an additive
    # error on the log-transformed observed concentrations (equivalent
    # to an exponential error on an arithmetic scale).' That is exactly
    # a log-normal residual, so the canonical `expSd` / `lnorm()` pair
    # is used rather than the small-sigma proportional approximation.
    #
    # Table 2's sigma row is a VARIANCE: its footnote defines sigma_PK
    # as 'residual exponential error variance of drug measurements',
    # and the companion sigma_PD of 146 for the QTc model is only
    # dimensionally coherent as a variance (146 ms^2 -> SD 12.1 ms; an
    # SD of 146 ms is impossible for an ECG interval residual). The
    # log-scale SD is therefore sqrt(0.137) = 0.370.
    # ================================================================
    expSd <- sqrt(0.137)
    label("Log-scale residual standard deviation for plasma piperaquine (unitless)")
    # Chotsiri 2017 Table 2: sigma_PK = 0.137 (%RSE 9.22; bootstrap 95% CI 0.111-0.161), a variance -> SD = sqrt(0.137)
  })

  model({
    # ---- Occasion multiplexing --------------------------------------
    # Two dosing occasions of dihydroartemisinin-piperaquine per
    # subject. OCC outside 1..2 (e.g. 0) collapses the term to zero and
    # returns the occasion-free typical prediction.
    iov_fdepot <- (OCC == 1) * etaiov_fdepot_1 + (OCC == 2) * etaiov_fdepot_2
    iov_mtt <- (OCC == 1) * etaiov_mtt_1 + (OCC == 2) * etaiov_mtt_2

    # ---- Allometric scaling on the study-median 64 kg ---------------
    allom_cl <- (WT / 64)^e_wt_cl
    allom_v <- (WT / 64)^e_wt_vc

    # ---- Individual parameters --------------------------------------
    mtt <- exp(lmtt + iov_mtt)

    cl <- exp(lcl + etalcl) * allom_cl
    vc <- exp(lvc + etalvc) * allom_v
    q <- exp(lq) * allom_cl
    vp <- exp(lvp) * allom_v
    q2 <- exp(lq2 + etalq2) * allom_cl
    vp2 <- exp(lvp2) * allom_v

    fdepot <- exp(lfdepot + etalfdepot + iov_fdepot)

    # ka = ktr, so the mean transit time spans all three transfers.
    ktr <- 3 / mtt

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ---- ODE system --------------------------------------------------
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(central) <- ktr * transit2 - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    f(depot) <- fdepot

    # ---- Observation -------------------------------------------------
    # Doses are in mg of piperaquine base and volumes in L, so
    # central / vc is mg/L; the factor 1000 converts to the assay's
    # ng/mL.
    Cc <- 1000 * central / vc

    Cc ~ lnorm(expSd)
  })
}
