Lallemand_2023_benzylpenicillin_horse <- function() {
  description <- "Preclinical (horse). Three-compartment population pharmacokinetic model for benzylpenicillin (BP) in horses, developed as an international meta-analysis to derive PK/PD cutoff values supporting equine clinical breakpoints (VetCAST approach). Seven parallel absorption depots reproduce the five modalities the authors modelled simultaneously: (1) IV sodium or potassium BP dosed directly into the central compartment; (2) IM sodium BP with two SEQUENTIAL first-order rate constants, the rapid Ka being replaced by the slower Ka2 once time-after-dose exceeds a lag time; (3) IM procaine BP from four single-ingredient products (France Depocilline, Sweden Penovet, USA1 Norocillin, Japan Procaine BP G sol 'KS'), sharing one Ka and one logit-scale bioavailability that are each modified by a 'source of dataset' covariate; (4) the Duplocilline fixed combination, whose total BP dose is split by SPC composition into a procaine BP fraction and a benzathine BP fraction absorbed in PARALLEL, each with its own Ka and bioavailability; and (5) IM penethamate hydriodide given at three separate injection sites, each site carrying its own sequential rapid-then-slow Ka pair and its own switch delay. All seven depots and the IV route share one set of disposition parameters (Vc, Vp, Vp2, CL, Q, Q2), which were FROZEN in the full model at the values estimated from the 23 IV profiles alone; plasma clearance carries a 'source of dataset' covariate with France as reference. The PK/PD layer reproduces the paper's Monte Carlo machinery: a fixed unbound fraction of 0.4 converts total to free BP, and the free AUC and the cumulative time above MIC are integrated as states so that fAUC/MIC and fT>MIC can be read directly. Pooled meta-analysis of 40 horses providing 63 rich profiles and 1022 BP plasma concentrations from France, Sweden, Japan and the USA (Lallemand 2023)."
  reference <- paste(
    "Lallemand EA, Bousquet-Melou A, Chapuis L, Davis J, Ferran AA, Kukanich B,",
    "Kuroda T, Lacroix MZ, Minamijima Y, Olsen L, Pelligand L, Portugal FR,",
    "Roques BB, Santschi EM, Wilson KE, Toutain P-L.",
    "Pharmacokinetic-pharmacodynamic cutoff values for benzylpenicillin in horses",
    "to support the establishment of clinical breakpoints for benzylpenicillin",
    "antimicrobial susceptibility testing in horses.",
    "Front Microbiol. 2023;14:1282949. doi:10.3389/fmicb.2023.1282949",
    sep = " "
  )
  vignette <- "Lallemand_2023_benzylpenicillin_horse"
  units <- list(time = "h", dosing = "ug/kg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot1      = list(analyte = "benzylpenicillin", units = NA_character_, specimen = "administration site", verified = FALSE),
    depot2      = list(analyte = "benzylpenicillin", units = NA_character_, specimen = "administration site", verified = FALSE),
    depot3      = list(analyte = "benzylpenicillin", units = NA_character_, specimen = "administration site", verified = FALSE),
    depot4      = list(analyte = "benzylpenicillin", units = NA_character_, specimen = "administration site", verified = FALSE),
    depot5      = list(analyte = "benzylpenicillin", units = NA_character_, specimen = "administration site", verified = FALSE),
    depot6      = list(analyte = "benzylpenicillin", units = NA_character_, specimen = "administration site", verified = FALSE),
    depot7      = list(analyte = "benzylpenicillin", units = NA_character_, specimen = "administration site", verified = FALSE),
    central     = list(analyte = "benzylpenicillin", units = NA_character_, specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "benzylpenicillin", units = NA_character_, specimen = "plasma", verified = FALSE),
    peripheral2 = list(analyte = "benzylpenicillin", units = NA_character_, specimen = "plasma", verified = FALSE),
    auc_free    = list(analyte = "benzylpenicillin", units = NA_character_, specimen = "not applicable", verified = FALSE),
    t_above_mic = list(analyte = "benzylpenicillin", units = NA_character_, specimen = "not applicable", verified = FALSE)
  )

  covariateData <- list(
    STUDY_SWEDEN_IV = list(
      description        = "Swedish IV cohort indicator (Lallemand 2023 'source of dataset' level Nationcode = 1)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (French cohort, Nationcode = 0, the reference source of dataset)",
      notes              = "1 = one of the 4 Swedish horses that received a single IV dose of sodium BP (Olsen 2018) and no IM administration. Modifies plasma clearance only. Set to 0 for every other cohort. Mutually exclusive with the other four STUDY_* indicators; when all five are 0 the horse is French. Time-fixed per subject. The coefficient was estimated in the IV-only model (Lallemand 2023 Table 4) and then FROZEN in the full model.",
      source_name        = "Nationcode == 1"
    ),
    STUDY_SWEDEN_IM = list(
      description        = "Swedish IM cohort indicator (Lallemand 2023 'source of dataset' level Nationcode = 11)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (French cohort, Nationcode = 0, the reference source of dataset)",
      notes              = "1 = one of the 8 Swedish horses enrolled in the 2x2 crossover comparing four IM administrations of procaine BP against seven IM administrations of sodium BP (Olsen 2013). These horses received no IV dose, so their plasma clearance was estimated indirectly (Bayesian) in the full model rather than frozen from the IV analysis. Modifies plasma clearance, the procaine BP absorption rate constant, and the procaine BP logit bioavailability. A SEPARATE level from STUDY_SWEDEN_IV: Lallemand 2023 Section 2.2 explains there was no reason to assume these 8 horses shared the clearance of the 4 Swedish IV horses (243 vs 351 mL/kg/h). Time-fixed per subject.",
      source_name        = "Nationcode == 11"
    ),
    STUDY_USA1 = list(
      description        = "USA1 cohort indicator (Lallemand 2023 'source of dataset' level Nationcode = 2)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (French cohort, Nationcode = 0, the reference source of dataset)",
      notes              = "1 = one of the 7 USA1 horses (Younkin 2019) that received two IV administrations of potassium BP at a 6 h interval followed by a single IM dose of procaine BP (Norocillin) 6 h after the second IV dose. Modifies plasma clearance, the procaine BP absorption rate constant, and the procaine BP logit bioavailability. Only 3 sampling points were available so this dataset was not analysable by NCA, and sampling stopped 12 h after the IM dose; the resulting very slow absorption (MAT 230 h) should be treated with caution per Lallemand 2023 Discussion. Time-fixed per subject. The clearance coefficient was estimated in the IV-only model (Table 4) and then FROZEN in the full model.",
      source_name        = "Nationcode == 2"
    ),
    STUDY_USA2 = list(
      description        = "USA2 cohort indicator (Lallemand 2023 'source of dataset' level Nationcode = 4)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (French cohort, Nationcode = 0, the reference source of dataset)",
      notes              = "1 = one of the 6 USA2 horses (Wilson 2022) that received a single IV administration of potassium BP together with gentamicin 6.6 mg/kg. IV data only; no IM administration. Modifies plasma clearance only. Time-fixed per subject. The coefficient was estimated in the IV-only model (Table 4) and then FROZEN in the full model.",
      source_name        = "Nationcode == 4"
    ),
    STUDY_JAPAN = list(
      description        = "Japanese cohort indicator (Lallemand 2023 'source of dataset' level Nationcode = 3)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (French cohort, Nationcode = 0, the reference source of dataset)",
      notes              = "1 = one of the 6 Japanese horses that received a single IM dose of procaine BP (Procaine BP G sol for Animals 'KS'). No IV data were collected in Japan, so plasma clearance was estimated indirectly (Bayesian) in the full model. Modifies plasma clearance, the procaine BP absorption rate constant, and the procaine BP logit bioavailability. Time-fixed per subject.",
      source_name        = "Nationcode == 3"
    ),
    FORM_BP_NA_IM = list(
      description        = "Sodium benzylpenicillin solution given by the intramuscular route (Geepenil) formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (IV administration of soluble sodium or potassium BP, the reference observation block)",
      notes              = "1 = the observation was recorded after IM sodium BP. Lallemand 2023 estimated a SEPARATE additive residual-error standard deviation for each of the five observation blocks (Supplementary Table S7 rows stdev0-stdev4) while sharing one multiplicative CV across all blocks; this indicator selects the IM sodium BP additive SD (0.162677 ug/mL). Observation-level, not subject-level: a horse in the Swedish crossover contributes rows with this indicator set to 1 on its sodium BP occasion and to 0 on its procaine BP occasion. Mutually exclusive with the other three FORM_BP_* indicators.",
      source_name        = "observation block CObsPeniNa in Supplementary Data 2"
    ),
    FORM_BP_PROC = list(
      description        = "Procaine benzylpenicillin single-ingredient IM suspension formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (IV administration of soluble sodium or potassium BP, the reference observation block)",
      notes              = "1 = the observation was recorded after IM procaine BP from one of the four single-ingredient products (Depocilline, Penovet, Norocillin, or the Japanese 'KS' product). Selects the procaine BP additive residual SD (0.004579 ug/mL) per Lallemand 2023 Supplementary Table S7 row stdev3. Does NOT cover the procaine BP contained in the Duplocilline combination, which forms its own observation block (see FORM_BP_DUPLO). Observation-level. Mutually exclusive with the other three FORM_BP_* indicators.",
      source_name        = "observation block CObsPROC in Supplementary Data 2"
    ),
    FORM_BP_DUPLO = list(
      description        = "Duplocilline (procaine benzylpenicillin + benzathine benzylpenicillin fixed combination) IM suspension formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (IV administration of soluble sodium or potassium BP, the reference observation block)",
      notes              = "1 = the observation was recorded after IM Duplocilline. Selects the Duplocilline additive residual SD (0.006537 ug/mL) per Lallemand 2023 Supplementary Table S7 row stdev2. The measured concentration is the SUM of the BP released by both prodrug components, which the model reproduces by feeding both the procaine and the benzathine depot into the same central compartment. Observation-level. Mutually exclusive with the other three FORM_BP_* indicators.",
      source_name        = "observation block CObsBenza in Supplementary Data 2"
    ),
    FORM_BP_PENETH = list(
      description        = "Penethamate hydriodide (Penetavet) IM suspension formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (IV administration of soluble sodium or potassium BP, the reference observation block)",
      notes              = "1 = the observation was recorded after IM penethamate hydriodide. Selects the penethamate additive residual SD (0.008068 ug/mL) per Lallemand 2023 Supplementary Table S7 row stdev4. The measured concentration is the SUM of the BP released from all three injection sites, which the model reproduces by feeding the three site depots into the same central compartment. Observation-level. Mutually exclusive with the other three FORM_BP_* indicators.",
      source_name        = "observation block CObsPenethamate in Supplementary Data 2"
    )
  )

  population <- list(
    species        = "horse (Equus caballus); Standardbred trotters (Sweden), Thoroughbreds (Japan, USA2), Quarter Horses / Quarter Horse crosses / Paint / Thoroughbred (USA1), French Trotters, Arabian and Spanish half-bloods and a Saddlebred (France), and one mixed pony (Sweden)",
    n_subjects     = 40L,
    n_studies      = 5L,
    n_profiles     = 63L,
    n_observations = 1022L,
    age_range      = "4-24 years (Supplementary Table S3; age not reported for 2 of the 40 horses)",
    weight_range   = "308-592 kg (Supplementary Table S3)",
    sex_female_pct = 50,
    sex_notes      = "20 mares/females, 15 geldings, 4 stallions/males, 1 horse with sex not reported (Supplementary Table S3)",
    disease_state  = "healthy adult horses enrolled in experimental pharmacokinetic trials; no clinical infection",
    dose_range     = "5.617-25.11 mg/kg expressed as benzylpenicillin base (Supplementary Table S2); the reference regimen throughout the paper is 12.36 mg/kg = 22,000 IU/kg",
    regions        = "France, Sweden, Japan, United States (two independent US datasets, USA1 and USA2)",
    formulations   = "IV sodium BP (France, Sweden), IV potassium BP (USA1, USA2), IM procaine BP (France Depocilline, Sweden Penovet, USA1 Norocillin, Japan 'KS'), IM procaine BP + benzathine BP (France Duplocilline), IM penethamate hydriodide (France Penetavet), IM sodium BP (Sweden Geepenil)",
    notes          = "Meta-analysis assembled for the VetCAST clinical-breakpoint project. The French data were generated specifically for this project and the Japanese data were previously unpublished; the Swedish (Olsen 2013, 2018), USA1 (Younkin 2019) and USA2 (Wilson 2022) data had been published before. Designs were deliberately unbalanced: no IV data in Japan, IV data only in USA2, and limits of quantification differing by site (5.5 ng/mL Sweden, 10 ng/mL France and USA1, 30 ng/mL Japan, 40 ng/mL USA2). BLOQ values (<5% of measurements) were discarded. All concentrations measured by LC-MS/MS. Because the goal was a maximally generic PK/PD cutoff, NO biological covariates (body weight, age, sex, breed) were carried into the model; 'source of dataset' is the only covariate. Extent of BP plasma protein binding measured at 62.8% in Swedish horses (Olsen 2013), giving the unbound fraction of 0.4 used for the PK/PD indices."
  )

  paper_specific_compartments <- c("auc_free", "t_above_mic")

  ini({
    # ---------------------------------------------------------------------
    # Disposition, shared by every route and formulation.
    # FROZEN in the full model at the values estimated from the 23 IV
    # profiles alone (Lallemand 2023 Section 3.3: "the disposition
    # parameters for French horses (CL, CL2, CL3, Vc, V2, and V3) were
    # used as a reference"; Supplementary Data 2 marks each of these six
    # fixef() declarations "(freeze)").
    # Volumes in mL/kg and clearances in mL/kg/h, so a dose in ug/kg
    # gives a concentration in ug/mL.
    # ---------------------------------------------------------------------
    lvc <- fixed(log(106.184336892764))
    label("Central volume of distribution Vc (mL/kg; frozen from the IV analysis)")  # Supplementary Data 2 fixef(tvV(freeze)); rounded to 106 in Table 4 and Supplementary Table S7
    lvp <- fixed(log(46.2952044536849))
    label("First peripheral volume of distribution V2 (mL/kg; frozen from the IV analysis)")  # Supplementary Data 2 fixef(tvV2(freeze)); rounded to 46.3 in Table 4
    lvp2 <- fixed(log(50.4297016344536))
    label("Second peripheral volume of distribution V3 (mL/kg; frozen from the IV analysis)")  # Supplementary Data 2 fixef(tvV3(freeze)); rounded to 50.4 in Table 4
    lcl <- fixed(log(480.899353392651))
    label("Plasma clearance CL for French horses, the reference source of dataset (mL/kg/h; frozen from the IV analysis)")  # Supplementary Data 2 fixef(tvCl(freeze)); rounded to 481 in Table 4 and Supplementary Table S7
    lq <- fixed(log(125.961026710481))
    label("Inter-compartmental clearance CL2 to the first peripheral compartment (mL/kg/h; frozen from the IV analysis)")  # Supplementary Data 2 fixef(tvCl2(freeze)); rounded to 126 in Table 4
    lq2 <- fixed(log(24.9818543761259))
    label("Inter-compartmental clearance CL3 to the second peripheral compartment (mL/kg/h; frozen from the IV analysis)")  # Supplementary Data 2 fixef(tvCl3(freeze)); rounded to 25.0 in Table 4

    # ---------------------------------------------------------------------
    # 'Source of dataset' covariate on plasma clearance, log-additive
    # (Phoenix: Cl = tvCl * exp(dCldNationcode_k * (Nationcode == k))).
    # France (all indicators 0) is the reference. The three coefficients
    # that came from the IV-only analysis are FIXED; the two that could
    # only be estimated indirectly in the full model (Japan and the
    # Swedish IM cohort, neither of which received an IV dose) are
    # estimated. See the vignette Errata for why this split follows the
    # paper prose rather than the (freeze) flags in Supplementary Data 2.
    # ---------------------------------------------------------------------
    e_study_sweden_iv_cl <- fixed(-0.683254)
    label("Effect of the Swedish IV cohort on log CL (scalar; frozen from the IV analysis)")  # Table 4 dCldSource_of_data1 = -0.68325471; Supplementary Table S7 dCldNationcode1 = -0.6833. Gives CL = 243 mL/kg/h
    e_study_usa1_cl <- fixed(-0.125766)
    label("Effect of the USA1 cohort on log CL (scalar; frozen from the IV analysis)")  # Table 4 dCldSource_of_data2 = -0.12576635; Supplementary Table S7 dCldNationcode2 = -0.1258. Gives CL = 424 mL/kg/h. Not significant (95% CI -0.329 to 0.077) but retained
    e_study_usa2_cl <- fixed(-0.136097)
    label("Effect of the USA2 cohort on log CL (scalar; frozen from the IV analysis)")  # Table 4 dCldSource_of_data3 = -0.13609778; Supplementary Table S7 dCldNationcode4 = -0.1361. Gives CL = 420 mL/kg/h, i.e. 87.3% of the French value
    e_study_japan_cl <- -0.156691638818043
    label("Effect of the Japanese cohort on log CL (scalar; Bayesian estimate from the full model)")  # Supplementary Data 2 fixef(dCldNationcode3); Supplementary Table S7 = -0.1567. Gives CL = 411 mL/kg/h (Table 6, Supplementary Table S8 ClearanceJP)
    e_study_sweden_im_cl <- -0.314630693388979
    label("Effect of the Swedish IM cohort on log CL (scalar; Bayesian estimate from the full model)")  # Supplementary Data 2 fixef(dCldNationcode11); Supplementary Table S7 = -0.3146. Gives CL = 351 mL/kg/h (Table 6, Supplementary Table S8 ClearanceSW11)

    # ---------------------------------------------------------------------
    # Depot 1: IM sodium BP (Geepenil, Swedish horses).
    # Two SEQUENTIAL first-order rate constants: the rapid Ka acts from
    # each dose until time-after-dose reaches the lag time, after which
    # the slower Ka2 takes over. In Supplementary Data 2 this is the
    # SW1PENI / SW2PENI switch driven by a sequence{} block.
    # Bioavailability is multiplicative here, NOT ilogit-transformed
    # (Supplementary Data 2: stparm(FPeniNa = tvFPeniNa * exp(nFPeniNa))).
    # ---------------------------------------------------------------------
    lka_im_na <- log(1.02419485674906)
    label("Rapid first-order absorption rate constant Ka1 after IM sodium BP (1/h)")  # Supplementary Data 2 fixef(tvKa1PeniNa); Supplementary Table S7 tvKa1BP-Na = 1.024. MAT1 = 1/Ka1 = 0.976 h (Table 6)
    lka2_im_na <- log(0.248123902641013)
    label("Slow first-order absorption rate constant Ka2 after IM sodium BP, active after the lag time (1/h)")  # Supplementary Data 2 fixef(tvKa2PeniNa); Supplementary Table S7 tvKa2BP-Na = 0.248. MAT2 = 1/Ka2 = 4.03 h (Table 6)
    ltlag_im_na <- log(0.470979905259535)
    label("Time after each IM sodium BP dose at which Ka1 is replaced by Ka2 (h)")  # Supplementary Data 2 fixef(tvtlagPeni); Supplementary Table S7 tvTlagBP-Na = 0.471 h (Table 6)
    lfdepot_im_na <- log(0.891328516980401)
    label("Bioavailability of BP after IM sodium BP (fraction; multiplicative scale)")  # Supplementary Data 2 fixef(tvFPeniNa) = 0.8913 -> 89.1% (Table 6, Supplementary Table S7)

    # ---------------------------------------------------------------------
    # Depot 2: IM procaine BP, single-ingredient products.
    # One first-order Ka and one logit-scale bioavailability, each
    # modified by the 'source of dataset' covariate with France as the
    # reference. Note the unusual covariate form the authors used for F:
    # the covariate MULTIPLIES the logit-scale typical value before the
    # inverse-logit is taken (Supplementary Data 2:
    # FPROC = ilogit(tvFPROC * exp(dFNationcode_k * (Nationcode == k)) + nFPROC)).
    # ---------------------------------------------------------------------
    lka_proc <- log(0.046594811352403)
    label("First-order absorption rate constant Ka after IM procaine BP, French formulation Depocilline (1/h)")  # Supplementary Data 2 fixef(tvKaPROC); Supplementary Table S7 tvKaPROC = 0.0466. MAT = 21.46 h (Table 6)
    e_study_sweden_im_ka_proc <- fixed(0.0390078308660437)
    label("Effect of the Swedish IM cohort on log Ka for procaine BP, Penovet (scalar)")  # Supplementary Data 2 fixef(dKadNationcode11(freeze)); Supplementary Table S7 = 0.0390. Gives Ka = 0.048 /h, MAT = 20.64 h (Table 6)
    e_study_usa1_ka_proc <- fixed(-2.37155662564783)
    label("Effect of the USA1 cohort on log Ka for procaine BP, Norocillin (scalar)")  # Supplementary Data 2 fixef(dKadNationcode2(freeze)); Supplementary Table S7 = -2.372. Gives Ka = 0.00435 /h, MAT = 230 h (Table 6)
    e_study_japan_ka_proc <- fixed(-0.00427442230486608)
    label("Effect of the Japanese cohort on log Ka for procaine BP (scalar)")  # Supplementary Data 2 fixef(dKadNationcode3(freeze)); Supplementary Table S7 = -0.0043. Gives Ka = 0.04640 /h, MAT = 21.55 h (Table 6)
    logitfdepot_proc <- 7.83103159739662
    label("Logit of the bioavailability of BP after IM procaine BP, French formulation Depocilline (scalar)")  # Supplementary Data 2 fixef(tvFPROC); Supplementary Table S7 tvFBP-PROC = 7.8310 -> ilogit = 99.960% (Supplementary Table S8)
    e_study_sweden_im_fdepot_proc <- -0.00177886964815447
    label("Multiplicative effect of the Swedish IM cohort on the logit bioavailability of procaine BP (scalar)")  # Supplementary Data 2 fixef(dFNationcode11); Supplementary Table S7 = -0.0018. Gives F = 99.960% (Supplementary Table S8)
    e_study_usa1_fdepot_proc <- fixed(0.0517205803762012)
    label("Multiplicative effect of the USA1 cohort on the logit bioavailability of procaine BP (scalar)")  # Supplementary Data 2 fixef(dFNationcode2(freeze)); Supplementary Table S7 = 0.052. Gives F = 99.974% (Supplementary Table S8)
    e_study_japan_fdepot_proc <- fixed(0.422679139787995)
    label("Multiplicative effect of the Japanese cohort on the logit bioavailability of procaine BP (scalar)")  # Supplementary Data 2 fixef(dFNationcode3(freeze)); Supplementary Table S7 = 0.423. Gives F = 99.999% (Supplementary Table S8)

    # ---------------------------------------------------------------------
    # Depots 3 and 4: Duplocilline, a fixed combination of procaine BP
    # and benzathine BP given to French horses. The total BP dose is
    # split by SPC composition and the two prodrug fractions are absorbed
    # in PARALLEL, each with its own Ka and its own ilogit bioavailability.
    # Both feed the same central compartment because the assay measures
    # the sum of the BP released by the two ingredients.
    # ---------------------------------------------------------------------
    lka_duplo_pro <- log(0.0827110143014943)
    label("First-order absorption rate constant Ka of the procaine BP fraction of Duplocilline (1/h)")  # Supplementary Data 2 fixef(tvKaBenza1); Supplementary Table S7 = 0.0827. MAT = 12.1 h (Table 6)
    logitfdepot_duplo_pro <- 1.53153195755327
    label("Logit of the bioavailability of BP from the procaine BP fraction of Duplocilline (scalar)")  # Supplementary Data 2 fixef(tvFBenza1); Supplementary Table S7 = 1.532 -> ilogit = 82.2% (Table 6)
    lka_duplo_benza <- log(0.00935419980483117)
    label("First-order absorption rate constant Ka of the benzathine BP fraction of Duplocilline (1/h)")  # Supplementary Data 2 fixef(tvKaBenza2); Supplementary Table S7 = 0.0094. MAT = 107 h (Table 6)
    logitfdepot_duplo_benza <- 10.9421855761952
    label("Logit of the bioavailability of BP from the benzathine BP fraction of Duplocilline (scalar)")  # Supplementary Data 2 fixef(tvFBenza2); Supplementary Table S7 = 10.9422 -> ilogit = 100.0% (Table 6)

    # ---------------------------------------------------------------------
    # Depots 5-7: IM penethamate hydriodide (Penetavet) given at three
    # separate injection sites, one per daily administration. Each site
    # has its own sequential rapid-then-slow Ka pair and its own switch
    # delay measured from that site's dose (Supplementary Data 2 uses one
    # sequence{} block per site precisely so each site's clock is
    # independent). All three sites share one ilogit bioavailability.
    # ---------------------------------------------------------------------
    lka_peneth1 <- log(0.0350071439943249)
    label("Rapid absorption rate constant of penethamate at injection site 1 (1/h)")  # Supplementary Data 2 fixef(tvKaPenethamate1); Supplementary Table S7 = 0.0350. MAT = 28.57 h (Supplementary Table S8)
    lka2_peneth1 <- log(0.0141803796444364)
    label("Slow absorption rate constant of penethamate at injection site 1, active after its lag time (1/h)")  # Supplementary Data 2 fixef(tvKaPenethamate_slow1); Supplementary Table S7 = 0.0142
    ltlag_peneth1 <- log(70.6482526955672)
    label("Time after the site 1 penethamate dose at which the rapid Ka is replaced by the slow Ka (h)")  # Supplementary Data 2 fixef(tvtlag1); Supplementary Table S7 tvTlag1 = 70.6 h
    lka_peneth2 <- log(0.0223223860881139)
    label("Rapid absorption rate constant of penethamate at injection site 2 (1/h)")  # Supplementary Data 2 fixef(tvKaPenethamate2); Supplementary Table S7 = 0.0223. MAT = 44.80 h (Supplementary Table S8)
    lka2_peneth2 <- log(0.108333675324903)
    label("Slow absorption rate constant of penethamate at injection site 2, active after its lag time (1/h)")  # Supplementary Data 2 fixef(tvKaPenethamate_slow2); Supplementary Table S7 = 0.1083
    ltlag_peneth2 <- log(27.2167070491644)
    label("Time after the site 2 penethamate dose at which the rapid Ka is replaced by the slow Ka (h)")  # Supplementary Data 2 fixef(tvtlag2); Supplementary Table S7 tvTlag2 = 27.2 h
    lka_peneth3 <- log(0.0151705494676296)
    label("Rapid absorption rate constant of penethamate at injection site 3 (1/h)")  # Supplementary Data 2 fixef(tvKaPenethamate3); Supplementary Table S7 = 0.0152. MAT = 65.92 h (Supplementary Table S8)
    lka2_peneth3 <- log(0.00162233822089588)
    label("Slow absorption rate constant of penethamate at injection site 3, active after its lag time (1/h)")  # Supplementary Data 2 fixef(tvKaPenethamate_slow3); Supplementary Table S7 = 0.0016
    ltlag_peneth3 <- log(42.3487955074074)
    label("Time after the site 3 penethamate dose at which the rapid Ka is replaced by the slow Ka (h)")  # Supplementary Data 2 fixef(tvtlag3); Supplementary Table S7 tvTlag3 = 42.3 h
    logitfdepot_peneth <- 0.788807915876975
    label("Logit of the bioavailability of BP from penethamate, shared by all three injection sites (scalar)")  # Supplementary Data 2 fixef(tvFPenethamate); Supplementary Table S7 = 0.7888 -> ilogit = 68.8% (Table 6)

    # ---------------------------------------------------------------------
    # PK/PD integration constants (Lallemand 2023 Section 2.3).
    # ---------------------------------------------------------------------
    fu <- fixed(0.4)
    label("Fraction of BP unbound in equine plasma (unitless)")  # Section 2.3: plasma protein binding measured at 62.8% in Swedish horses (Olsen 2013); "For our simulations, fu was introduced into the model, with a typical value of 0.4"
    mic <- fixed(0.25)
    label("Minimum inhibitory concentration used for the fT>MIC integration (mg/L, numerically equal to ug/mL)")  # Section 2.3 explored 0.0625, 0.125, 0.25, 0.375, 0.5, 1 and 2 mg/L; 0.25 mg/L is the PK/PD cutoff the paper concludes with. Change this to explore another MIC

    # ---------------------------------------------------------------------
    # Between-subject variability. Lallemand 2023 used an exponential
    # (log-normal) random-effect model throughout, with full OMEGA blocks
    # per formulation. Values are the full-model estimates from
    # Supplementary Data 2 ranef() blocks, read as row-major lower
    # triangles. BSV% = 100 * sqrt(exp(omega^2) - 1) (paper equation 2).
    # ---------------------------------------------------------------------
    # Disposition block, FROZEN in the full model
    # (Supplementary Data 2: ranef(block(nV, nCl, nV2, nV3, nCl2, nCl3)(freeze))).
    # Diagonals reproduce the BSV% quoted in the Discussion: Vc 174%,
    # CL 26.4%, V2 56.7%, V3 42%.
    etalvc + etalcl + etalvp + etalvp2 + etalq + etalq2 ~
      fixed(1.400552,
            0.20942975, 0.067363321,
            0.40218913, 0.087954509, 0.27982625,
            0.34658421, 0.06450593, 0.14584899, 0.16522015,
            0.31712857, 0.049439529, 0.1098892, 0.054298124, 0.23712655,
            0.33452693, 0.065632766, 0.077679041, 0.10881838, 0.12799386, 0.2143448)

    # IM sodium BP block, Phoenix order (nFPeniNa, ntlagPeni, nKa1PeniNA, nKa2PeniNA)
    etalfdepot_im_na + etaltlag_im_na + etalka_im_na + etalka2_im_na ~
      c(0.030733235,
        0.00072792607, 0.044429534,
        0.0021984035, -0.073040728, 0.80538569,
        0.0059777827, 0.0016714216, -0.16919682, 0.044433223)

    # Duplocilline block, Phoenix order (nKaBenza1, nKaBenza2, nFBenza1, nFBenza2)
    etalka_duplo_pro + etalka_duplo_benza + etalogitfdepot_duplo_pro + etalogitfdepot_duplo_benza ~
      c(0.12663874,
        0.1798759, 0.61633102,
        0.081958913, 0.77856969, 2.1630137,
        0.20130889, 0.58867596, 1.3440868, 1.160584)

    # Procaine BP block, Phoenix order (nKaPROC, nFPROC).
    # Supplementary Table S9 note: BSV of Ka is 34% with low shrinkage;
    # the large logit-scale BSV on F corresponds to a very small BSV on F
    # itself because bioavailability is essentially complete.
    etalka_proc + etalogitfdepot_proc ~ c(0.10963699, 0.52361095, 3.6624274)

    # Penethamate blocks, one per injection site, all FROZEN
    # (Supplementary Data 2 ranef(...)(freeze)). Each pairs the site's
    # rapid and slow Ka; the correlations are near 1 (0.999), so these
    # blocks are positive definite but ill-conditioned.
    etalka_peneth1 + etalka2_peneth1 ~ fixed(1.2160073, 0.24768147, 0.050500768)
    etalka_peneth2 + etalka2_peneth2 ~ fixed(0.094452959, 0.044713983, 0.021186634)
    etalka_peneth3 + etalka2_peneth3 ~ fixed(3.1459704, 0.30982059, 0.030569386)
    etalogitfdepot_peneth ~ 0.0047176976

    # ---------------------------------------------------------------------
    # Residual error: one multiplicative (proportional) CV shared by every
    # observation block, plus a block-specific additive standard deviation.
    # Phoenix wrote this as
    #   observe(CObs = C + CEps * sqrt(1 + C^2 * (CMultStdev/sigma())^2))
    # with error(CEps = <block additive SD>), which expands to a residual
    # variance of addSd^2 + (propSd * C)^2 -- i.e. exactly nlmixr2's
    # add() + prop() combination.
    # ---------------------------------------------------------------------
    propSd <- 0.283453226691478
    label("Proportional (multiplicative) residual error CV, shared by all observation blocks (fraction)")  # Supplementary Data 2 fixef(tvCMultStdev); Table 6 and Supplementary Table S7 = 0.283
    addSd_iv <- 0.00792812036590013
    label("Additive residual error SD for the IV observation block (ug/mL)")  # Supplementary Data 2 error(CEpsIV); Supplementary Table S7 stdev0 = 0.007928
    addSd_na_im <- 0.162677174710096
    label("Additive residual error SD for the IM sodium BP observation block (ug/mL)")  # Supplementary Data 2 error(CEpsPeniNa); Supplementary Table S7 stdev1 = 0.162677
    addSd_duplo <- 0.00653731719464946
    label("Additive residual error SD for the Duplocilline observation block (ug/mL)")  # Supplementary Data 2 error(CEpsBenza); Supplementary Table S7 stdev2 = 0.006537
    addSd_proc <- 0.00457902506048033
    label("Additive residual error SD for the IM procaine BP observation block (ug/mL)")  # Supplementary Data 2 error(CEpsPROC); Supplementary Table S7 stdev3 = 0.004579
    addSd_peneth <- 0.00806824116414487
    label("Additive residual error SD for the IM penethamate observation block (ug/mL)")  # Supplementary Data 2 error(CEpsPenethamate); Supplementary Table S7 stdev4 = 0.008068
  })

  model({
    # -------------------------------------------------------------------
    # Formulation composition constant.
    # Duplocilline delivers 12.4 mg/kg of BP split into 5.96376 mg/kg as
    # procaine BP and 6.43624 mg/kg as benzathine BP per the SPC
    # (Supplementary Table S2). The fraction is folded into the
    # bioavailability of each Duplocilline depot so that the user doses
    # the TOTAL BP amount into both depots.
    # -------------------------------------------------------------------
    frac_duplo_pro <- 5.96376 / 12.4  # = 0.4809484, unitless

    # -------------------------------------------------------------------
    # Individual disposition parameters
    # -------------------------------------------------------------------
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    vp2 <- exp(lvp2 + etalvp2)
    q <- exp(lq + etalq)
    q2 <- exp(lq2 + etalq2)
    cl <- exp(lcl + etalcl +
      e_study_sweden_iv_cl * STUDY_SWEDEN_IV +
      e_study_sweden_im_cl * STUDY_SWEDEN_IM +
      e_study_usa1_cl * STUDY_USA1 +
      e_study_usa2_cl * STUDY_USA2 +
      e_study_japan_cl * STUDY_JAPAN)

    # -------------------------------------------------------------------
    # Individual absorption parameters
    # -------------------------------------------------------------------
    # Depot 1, IM sodium BP: sequential Ka1 -> Ka2 at tlag
    ka_im_na <- exp(lka_im_na + etalka_im_na)
    ka2_im_na <- exp(lka2_im_na + etalka2_im_na)
    tlag_im_na <- exp(ltlag_im_na + etaltlag_im_na)
    f_im_na <- exp(lfdepot_im_na + etalfdepot_im_na)

    # Depot 2, IM procaine BP: Ka and logit F both carry the
    # 'source of dataset' covariate. For F the covariate multiplies the
    # logit-scale typical value, matching Supplementary Data 2.
    ka_proc <- exp(lka_proc + etalka_proc +
      e_study_sweden_im_ka_proc * STUDY_SWEDEN_IM +
      e_study_usa1_ka_proc * STUDY_USA1 +
      e_study_japan_ka_proc * STUDY_JAPAN)
    logit_f_proc <- logitfdepot_proc * exp(
      e_study_sweden_im_fdepot_proc * STUDY_SWEDEN_IM +
        e_study_usa1_fdepot_proc * STUDY_USA1 +
        e_study_japan_fdepot_proc * STUDY_JAPAN) + etalogitfdepot_proc
    f_proc <- 1 / (1 + exp(-logit_f_proc))

    # Depots 3 and 4, Duplocilline: two parallel prodrug fractions
    ka_duplo_pro <- exp(lka_duplo_pro + etalka_duplo_pro)
    ka_duplo_benza <- exp(lka_duplo_benza + etalka_duplo_benza)
    logit_f_duplo_pro <- logitfdepot_duplo_pro + etalogitfdepot_duplo_pro
    logit_f_duplo_benza <- logitfdepot_duplo_benza + etalogitfdepot_duplo_benza
    f_duplo_pro <- 1 / (1 + exp(-logit_f_duplo_pro))
    f_duplo_benza <- 1 / (1 + exp(-logit_f_duplo_benza))

    # Depots 5-7, penethamate: one sequential Ka pair per injection site
    ka_peneth1 <- exp(lka_peneth1 + etalka_peneth1)
    ka2_peneth1 <- exp(lka2_peneth1 + etalka2_peneth1)
    tlag_peneth1 <- exp(ltlag_peneth1)
    ka_peneth2 <- exp(lka_peneth2 + etalka_peneth2)
    ka2_peneth2 <- exp(lka2_peneth2 + etalka2_peneth2)
    tlag_peneth2 <- exp(ltlag_peneth2)
    ka_peneth3 <- exp(lka_peneth3 + etalka_peneth3)
    ka2_peneth3 <- exp(lka2_peneth3 + etalka2_peneth3)
    tlag_peneth3 <- exp(ltlag_peneth3)
    logit_f_peneth <- logitfdepot_peneth + etalogitfdepot_peneth
    f_peneth <- 1 / (1 + exp(-logit_f_peneth))

    # -------------------------------------------------------------------
    # Sequential-absorption switches. tad(<depot>) is the time since that
    # depot's most recent dose, so each administration restarts its own
    # rapid-then-slow sequence exactly as the Phoenix sequence{} blocks
    # did, without hard-coding the dose times.
    # -------------------------------------------------------------------
    ka_na_active <- ifelse(tad(depot1) < tlag_im_na, ka_im_na, ka2_im_na)
    ka_pen1_active <- ifelse(tad(depot5) < tlag_peneth1, ka_peneth1, ka2_peneth1)
    ka_pen2_active <- ifelse(tad(depot6) < tlag_peneth2, ka_peneth2, ka2_peneth2)
    ka_pen3_active <- ifelse(tad(depot7) < tlag_peneth3, ka_peneth3, ka2_peneth3)

    # -------------------------------------------------------------------
    # Concentrations and ODE system. The clearance parameterisation
    # mirrors Supplementary Data 2 exactly:
    #   deriv(A1 = -Cl*C1 - Cl2*(C1 - C2) - Cl3*(C1 - C3))
    # -------------------------------------------------------------------
    Cc <- central / vc
    Cp <- peripheral1 / vp
    Cp2 <- peripheral2 / vp2

    d/dt(depot1) <- -ka_na_active * depot1
    d/dt(depot2) <- -ka_proc * depot2
    d/dt(depot3) <- -ka_duplo_pro * depot3
    d/dt(depot4) <- -ka_duplo_benza * depot4
    d/dt(depot5) <- -ka_pen1_active * depot5
    d/dt(depot6) <- -ka_pen2_active * depot6
    d/dt(depot7) <- -ka_pen3_active * depot7
    d/dt(central) <- ka_na_active * depot1 +
      ka_proc * depot2 +
      ka_duplo_pro * depot3 +
      ka_duplo_benza * depot4 +
      ka_pen1_active * depot5 +
      ka_pen2_active * depot6 +
      ka_pen3_active * depot7 -
      cl * Cc - q * (Cc - Cp) - q2 * (Cc - Cp2)
    d/dt(peripheral1) <- q * (Cc - Cp)
    d/dt(peripheral2) <- q2 * (Cc - Cp2)

    # -------------------------------------------------------------------
    # Bioavailability. For the two Duplocilline depots the SPC composition
    # fraction is folded in, so the user doses the total BP amount into
    # both depot3 and depot4.
    # -------------------------------------------------------------------
    f(depot1) <- f_im_na
    f(depot2) <- f_proc
    f(depot3) <- frac_duplo_pro * f_duplo_pro
    f(depot4) <- (1 - frac_duplo_pro) * f_duplo_benza
    f(depot5) <- f_peneth
    f(depot6) <- f_peneth
    f(depot7) <- f_peneth

    # -------------------------------------------------------------------
    # PK/PD integration (Lallemand 2023 Section 2.3). The free AUC and the
    # cumulative time above MIC are carried as states so that
    # fAUC/MIC = auc_free / mic and fT>MIC = t_above_mic can be read
    # directly off the solution, reproducing the Phoenix
    # deriv(AUC_...) / deriv(T_above_...) constructs.
    # -------------------------------------------------------------------
    Cu <- fu * Cc
    d/dt(auc_free) <- Cu
    d/dt(t_above_mic) <- (Cu >= mic)

    # -------------------------------------------------------------------
    # Observation. One analyte (BP in plasma) with a shared proportional
    # CV and an observation-block-specific additive SD; the indicator form
    # below returns addSd_iv when all four IM formulation flags are 0.
    # -------------------------------------------------------------------
    addSd <- addSd_iv +
      (addSd_na_im - addSd_iv) * FORM_BP_NA_IM +
      (addSd_proc - addSd_iv) * FORM_BP_PROC +
      (addSd_duplo - addSd_iv) * FORM_BP_DUPLO +
      (addSd_peneth - addSd_iv) * FORM_BP_PENETH
    Cc ~ add(addSd) + prop(propSd)
  })
}
