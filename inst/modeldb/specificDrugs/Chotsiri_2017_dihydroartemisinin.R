Chotsiri_2017_dihydroartemisinin <- function() {
  description <- paste(
    "Two-compartment population pharmacokinetic model for dihydroartemisinin",
    "(DHA) in 16 healthy Thai adult volunteers, from Chotsiri 2017. A single",
    "oral dose of three co-formulated dihydroartemisinin-piperaquine tablets",
    "(120 mg DHA) was given on two separate occasions, with and without a",
    "concomitant single low dose of primaquine, in an open-label randomized",
    "three-way crossover. Absorption is a six-transit-compartment chain in",
    "which the transit rate constant and the rate constant out of the last",
    "transit compartment into central were estimated separately; body weight",
    "enters as a fixed allometric function on all clearance and volume terms",
    "centred on the study-median 64 kg. Variability is dominated by absorption:",
    "between-occasion variability on relative bioavailability, mean transit",
    "time and the absorption rate constant, with between-subject variability",
    "only on elimination clearance. Primaquine coadministration was screened",
    "on every pharmacokinetic parameter and produced no clinically relevant",
    "interaction, so no primaquine covariate appears in the final model.",
    "Sister model files from the same paper:",
    "modellib('Chotsiri_2017_piperaquine') (piperaquine population PK) and",
    "modellib('Chotsiri_2017_piperaquine_qtc') (piperaquine PK plus the linear",
    "concentration-QTc prolongation model).",
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
    "dihydroartemisinin'); secondary exposure parameters used for validation",
    "are in Table 3. Sister model files from the same paper:",
    "modellib('Chotsiri_2017_piperaquine') and",
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
        "Fixed allometric function on all clearance parameters (CL/F, Q/F;",
        "exponent 0.75) and all volume parameters (Vc/F, Vp/F; exponent 1),",
        "centred on the study-median body weight of 64 kg (Methods,",
        "'Population pharmacokinetic analysis', Equations 3 and 4). Adding it",
        "did not improve the DHA fit (dOFV = 0.819) but it was retained in the",
        "final model 'based on the strong biological prior and previously",
        "published results' (Results, 'Population pharmacokinetic properties",
        "of DHA'). Set WT = 64 to recover the tabulated typical values.",
        sep = " "
      ),
      source_name = "BW"
    ),
    OCC = list(
      description = paste(
        "Integer-valued occasion index identifying which of the two",
        "dihydroartemisinin-piperaquine dosing occasions a dose belongs to,",
        "for between-occasion variability on relative bioavailability, mean",
        "transit time and the absorption rate constant.",
        sep = " "
      ),
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Each volunteer received dihydroartemisinin-piperaquine twice -- once",
        "alone and once with primaquine -- in randomised order, separated by",
        "an 8-week washout (Methods, 'Study design'). The paper reports eta",
        "shrinkage separately 'on study occasions 1 and 2' (Results,",
        "'Population pharmacokinetic properties of DHA'), confirming two",
        "occasions. Values 1 and 2 are multiplexed inside model() onto",
        "etaiov_fdepot_1 / etaiov_fdepot_2, etaiov_mtt_1 / etaiov_mtt_2 and",
        "etaiov_ka_1 / etaiov_ka_2, each pair sharing one variance. OCC = 0,",
        "or any value outside 1..2, yields the occasion-free typical value.",
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
        "not significant; in the full covariate approach (a categorical",
        "primaquine effect placed on every pharmacokinetic parameter",
        "simultaneously, bootstrapped n = 1000) 'the impact of primaquine",
        "coadministration was less than +/-25% on primary pharmacokinetic",
        "parameters' for DHA (Results, 'Drug-drug interactions'; Figure 3A).",
        "The paper reports the bootstrap densities graphically only, with no",
        "point estimate, so no coefficient can be encoded.",
        sep = " "
      )
    ),
    AGE = list(
      description = "Subject age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened by stepwise forward inclusion / backward elimination and not retained: 'No significant covariates were identified in the stepwise covariate approach' (Results). No coefficient is reported."
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
    depot = list(analyte = "dihydroartemisinin", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "dihydroartemisinin", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "dihydroartemisinin", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "dihydroartemisinin", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "dihydroartemisinin", units = "mg", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "dihydroartemisinin", units = "mg", specimen = "administration site", verified = TRUE),
    transit6 = list(analyte = "dihydroartemisinin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "dihydroartemisinin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dihydroartemisinin", units = "mg", specimen = "plasma", verified = TRUE)
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
    dose_range = "Single oral dose of three co-formulated dihydroartemisinin-piperaquine tablets (40 mg dihydroartemisinin + 320 mg piperaquine phosphate each), i.e. 120 mg dihydroartemisinin, given 30 min after a light meal (~200 kcal, 8 g fat).",
    regions = "Bangkok, Thailand",
    n_observations = 384L,
    notes = paste(
      "Open-label, randomized, three-way crossover (NCT01525511, TMEC 12-004,",
      "OXTREC 58-11) conducted 18 June to 2 November 2012 at the Faculty of",
      "Tropical Medicine, Mahidol University. Every volunteer received",
      "primaquine alone first (1-week washout), then dihydroartemisinin-",
      "piperaquine alone and dihydroartemisinin-piperaquine plus primaquine in",
      "random order with an 8-week washout. Plasma sampling at 0, 0.25, 0.5, 1,",
      "1.5, 2, 3, 4, 6, 8, 10, 12 and 24 h post dose; LC-MS/MS with a lower",
      "limit of quantification of 2.00 ng/mL. 15% of DHA concentrations were",
      "below the limit of quantification and were omitted (M1); an M3",
      "categorical treatment gave similar performance. Estimated in NONMEM",
      "7.3 with FOCE-I. Baseline demographics are in Table 1. Note the",
      "Discussion states 'three males and 13 female' whereas the Methods",
      "state 'five males out of 16 subjects'; sex_female_pct here uses the",
      "Methods figure (11 of 16). See the vignette Errata.",
      sep = " "
    )
  )

  ini({
    # ================================================================
    # Structural parameters -- Chotsiri 2017 Table 2, 'Pharmacokinetic
    # parameters of dihydroartemisinin'. Table 2 footnote a: 'Parameter
    # estimates are based on the typical individual in the population
    # with a body weight of 64 kg.' All disposition parameters are
    # apparent (relative to F, which is fixed to unity).
    # ================================================================

    # ---- Absorption -------------------------------------------------
    # Six transit compartments (Results, 'Population pharmacokinetic
    # properties of DHA': 'A transit compartment absorption model with
    # six transit compartments was superior to all other absorption
    # models evaluated'). Unlike the piperaquine model of the same
    # paper, the transit rate constant and the rate constant out of the
    # last transit compartment were estimated SEPARATELY here
    # ('Estimating both the transit rate between transit compartments
    # and the absorption rate from the last transit compartment to the
    # central compartment resulted in a significantly improved model
    # fit compared with when setting them to be equal, dOFV = -17.6').
    #
    # MTT therefore spans only the six ktr-governed transfers, so
    # ktr = 6 / MTT and the seventh transfer runs at ka. That reading
    # is confirmed against the paper's own secondary parameters: it
    # reproduces the published median Tmax of 1.27-1.30 h (1.28 h
    # simulated for the typical 64-kg subject), whereas the alternative
    # ktr = 7 / MTT gives 1.18 h. See the vignette Source trace.
    lmtt <- log(0.567)
    label("Mean transit time through the six-transit-compartment absorption chain MTT (h)")
    # Chotsiri 2017 Table 2: MTT = 0.567 h (%RSE 11.4; bootstrap 95% CI 0.527-0.818)

    lka <- log(2.89)
    label("Absorption rate constant from the last transit compartment into central (1/h)")
    # Chotsiri 2017 Table 2: k_a = 2.89 /h (%RSE 37.1; bootstrap 95% CI 1.88-6.99)

    # ---- Disposition ------------------------------------------------
    # Two-compartment disposition (Results: a two-compartment model beat
    # a one-compartment model both with BLQ omitted, dOFV = -26.0, and
    # under the M3 method, dOFV = -12.3; a third compartment added only
    # dOFV = -6.55 and gave an implausibly long terminal half-life).
    lcl <- log(148)
    label("Apparent elimination clearance CL/F at WT = 64 kg (L/h)")
    # Chotsiri 2017 Table 2: CL/F = 148 L/h (%RSE 10.6; bootstrap 95% CI 121-183)

    lvc <- log(214)
    label("Apparent central volume of distribution Vc/F at WT = 64 kg (L)")
    # Chotsiri 2017 Table 2: V_C/F = 214 L (%RSE 16.9; bootstrap 95% CI 148-287)

    lq <- log(28.5)
    label("Apparent inter-compartmental clearance Q/F at WT = 64 kg (L/h)")
    # Chotsiri 2017 Table 2: Q_P/F = 28.5 L/h (%RSE 26.0; bootstrap 95% CI 15.5-44.1)

    lvp <- log(65.9)
    label("Apparent peripheral volume of distribution Vp/F at WT = 64 kg (L)")
    # Chotsiri 2017 Table 2: V_P/F = 65.9 L (%RSE 19.1; bootstrap 95% CI 42.5-91.5)

    # ---- Relative bioavailability -----------------------------------
    lfdepot <- fixed(log(1))
    label("Relative bioavailability F on an occasion-free reference occasion (unitless)")
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
    label("Allometric exponent of body weight on CL/F and Q/F, referenced to 64 kg (unitless)")
    # Chotsiri 2017 Methods, Equation 3: clearance exponent 3/4

    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on Vc/F and Vp/F, referenced to 64 kg (unitless)")
    # Chotsiri 2017 Methods, Equation 4: volume exponent 1

    # ================================================================
    # Variability. Table 2's variability column is headed '%CV of
    # BSV/BOV' and footnote a gives the transform explicitly: the %CV
    # is '100 * sqrt(exp(estimate) - 1)'. Converting back to the
    # internal variance scale is therefore omega^2 = log(CV^2 + 1).
    # Rows marked '*' in Table 2 are between-OCCASION variability; the
    # unmarked row is between-SUBJECT variability.
    # ================================================================

    # Between-subject variability on elimination clearance only
    # (Results: 'substantial between-occasion variability in the
    # absorption of DHA, with additional between-subject variability in
    # the elimination clearance of DHA').
    etalcl ~ 0.0519860
    # Chotsiri 2017 Table 2: CL/F BSV = 23.1% CV (%RSE 14.2; 95% CI 15.2-27.6); log(0.231^2 + 1) = 0.0519860

    # Between-occasion variability, two occasions (Methods, 'Study
    # design'; Results reports shrinkage separately 'on study occasions
    # 1 and 2'). Each occasion's eta shares the single published
    # variance; OCC selects which one is active inside model().
    etaiov_fdepot_1 ~ 0.1212269
    # Chotsiri 2017 Table 2: F BOV = 35.9% CV (%RSE 20.1; 95% CI 21.4%-50.4%); log(0.359^2 + 1) = 0.1212269
    etaiov_fdepot_2 ~ fixed(0.1212269)
    # Same published BOV variance as occasion 1

    etaiov_mtt_1 ~ 0.2442598
    # Chotsiri 2017 Table 2: MTT BOV = 52.6% CV (%RSE 14.2; 95% CI 36.0%-67.6%); log(0.526^2 + 1) = 0.2442598
    etaiov_mtt_2 ~ fixed(0.2442598)
    # Same published BOV variance as occasion 1

    etaiov_ka_1 ~ 0.5833881
    # Chotsiri 2017 Table 2: k_a BOV = 89.0% CV (%RSE 23.7; 95% CI 46.0%-169%); log(0.890^2 + 1) = 0.5833881
    etaiov_ka_2 ~ fixed(0.5833881)
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
    # log-scale SD is therefore sqrt(0.358) = 0.598.
    # ================================================================
    expSd <- sqrt(0.358)
    label("Log-scale residual standard deviation for plasma dihydroartemisinin (unitless)")
    # Chotsiri 2017 Table 2: sigma_PK = 0.358 (%RSE 9.07; bootstrap 95% CI 0.292-0.418), a variance -> SD = sqrt(0.358)
  })

  model({
    # ---- Occasion multiplexing --------------------------------------
    # Two dosing occasions of dihydroartemisinin-piperaquine per
    # subject. OCC outside 1..2 (e.g. 0) collapses every term to zero
    # and returns the occasion-free typical prediction.
    iov_fdepot <- (OCC == 1) * etaiov_fdepot_1 + (OCC == 2) * etaiov_fdepot_2
    iov_mtt <- (OCC == 1) * etaiov_mtt_1 + (OCC == 2) * etaiov_mtt_2
    iov_ka <- (OCC == 1) * etaiov_ka_1 + (OCC == 2) * etaiov_ka_2

    # ---- Allometric scaling on the study-median 64 kg ---------------
    allom_cl <- (WT / 64)^e_wt_cl
    allom_v <- (WT / 64)^e_wt_vc

    # ---- Individual parameters --------------------------------------
    mtt <- exp(lmtt + iov_mtt)
    ka <- exp(lka + iov_ka)

    cl <- exp(lcl + etalcl) * allom_cl
    vc <- exp(lvc) * allom_v
    q <- exp(lq) * allom_cl
    vp <- exp(lvp) * allom_v

    fdepot <- exp(lfdepot + iov_fdepot)

    # Six transit transfers governed by ktr; the seventh, out of
    # transit6 into central, runs at the separately estimated ka.
    ktr <- 6 / mtt

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- ODE system --------------------------------------------------
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(transit6) <- ktr * transit5 - ka * transit6
    d/dt(central) <- ka * transit6 - kel * central -
      k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # ---- Observation -------------------------------------------------
    # Doses are in mg and volumes in L, so central / vc is mg/L; the
    # factor 1000 converts to the assay's ng/mL.
    Cc <- 1000 * central / vc

    Cc ~ lnorm(expSd)
  })
}
