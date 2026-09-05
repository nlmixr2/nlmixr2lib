RodriguezFernandez_2024_warfarin <- function() {
  description <- "Warfarin population PK/PD model for INR in Caribbean Hispanic (Puerto Rican) patients: a one-compartment first-order-oral PK layer fixed from Reyes-Gonzalez 2020 (CYP2C9-genotype-specific elimination rate constant, weight-proportional volume) driving an inhibitory sigmoid-Imax effect on the zero-order synthesis rate of two parallel three-compartment transit chains, combined into INR. VKORC1 -1639G>A is a predictor of both the baseline INR and the warfarin IC50, the latter modelled as a sum of per-allele contributions. Caribbean Hispanic IC50 values are 3-5x higher than previously published European estimates."
  reference <- paste(
    "Rodriguez-Fernandez K, Reynaldo-Fernandez G, Reyes-Gonzalez S, de las Barreras C,",
    "Rodriguez-Vera L, Vlaar C, Monbaliu J-CM, Stelzer T, Duconge J, Mangas-Sanjuan V.",
    "New insights into the role of VKORC1 polymorphisms for optimal warfarin dose selection",
    "in Caribbean Hispanic patients through an external validation of a population PK/PD model.",
    "Biomed Pharmacother. 2024;170:115977. doi:10.1016/j.biopha.2023.115977.",
    "PMID: 38056237. PMCID: PMC10853672.",
    "PD parameters from Table 2; Imax, INRmax, transit-chain length and the fixed MTT values",
    "from the Results ('Base population PK/PD model' and 'Final population PK/PD model').",
    "The PK layer is fixed (no IIV) from the upstream publication",
    "Reyes-Gonzalez S, de las Barreras C, Reynaldo G, Rodriguez-Vera L, Vlaar C, Lopez Mejias V,",
    "Monbaliu J-CM, Stelzer T, Mangas V, Duconge J. Genotype-driven pharmacokinetic simulations",
    "of warfarin levels in Puerto Ricans. Drug Metab Pers Ther. 2020;35(4).",
    "doi:10.1515/dmpt-2020-0135. PMID: 34704696. PMCID: PMC7892629 (Methods: ka, F, Vd/kg and",
    "the six CYP2C9-genotype elimination rate constants).",
    "The transit-chain / INR structural form follows Hamberg 2007 and Hamberg 2010;",
    "see modellib('Hamberg_2007_warfarin_s')."
  )
  vignette <- "RodriguezFernandez_2024_warfarin"
  paper_specific_compartments <- c(
    "coag_s1", "coag_s2", "coag_s3",
    "coag_l1", "coag_l2", "coag_l3"
  )

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot   = list(analyte = "warfarin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "warfarin", units = "mg", specimen = "plasma", verified = TRUE),
    coag_s1 = list(analyte = "anticoagulant response", units = "fraction", specimen = "not applicable", verified = TRUE),
    coag_s2 = list(analyte = "anticoagulant response", units = "fraction", specimen = "not applicable", verified = TRUE),
    coag_s3 = list(analyte = "anticoagulant response", units = "fraction", specimen = "not applicable", verified = TRUE),
    coag_l1 = list(analyte = "anticoagulant response", units = "fraction", specimen = "not applicable", verified = TRUE),
    coag_l2 = list(analyte = "anticoagulant response", units = "fraction", specimen = "not applicable", verified = TRUE),
    coag_l3 = list(analyte = "anticoagulant response", units = "fraction", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Warfarin apparent volume of distribution is strictly proportional to body weight at 0.14 L/kg (Reyes-Gonzalez 2020 Methods, taken from the literature and not estimated), which is the 'body weight on V' covariate retained in the final model of Rodriguez-Fernandez 2024. Cohort median 83 kg, range 51-159 kg (Table 1).",
      source_name        = "body weight"
    ),
    CYP2C9_S1_COUNT = list(
      description        = "Count of CYP2C9*1 (wild-type) alleles per subject (0, 1, or 2)",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = "*1/*1 (CYP2C9_S1_COUNT == 2) is the reference diplotype; kel = 0.0189 1/h.",
      notes              = "CYP2C9_S1_COUNT + CYP2C9_S2_COUNT + CYP2C9_S3_COUNT + CYP2C9_S5_COUNT + CYP2C9_S6_COUNT + CYP2C9_S8_COUNT = 2 per subject. Reyes-Gonzalez 2020 groups *3/*5/*6/*8 together as the pooled reduced-function allele 'n', so the model derives a pooled n-count from the four columns. Cohort distribution (Rodriguez-Fernandez 2024 Table 1, n=138): *1/*1 71%, *1/*2 15%, *1/*3 5.1%, *1/*5 0.72%, *1/*8 1.44%, *2/*2 1.44%, *2/*3 4.34%, *2/*5 0.72%.",
      source_name        = "CYP2C9 genotype"
    ),
    CYP2C9_S2_COUNT = list(
      description        = "Count of CYP2C9*2 reduced-function alleles per subject (0, 1, or 2)",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "See CYP2C9_S1_COUNT. *2 carries its own elimination rate constants in Reyes-Gonzalez 2020 (*1/*2 0.0158, *2/*2 0.0130, *2/*n 0.009 1/h) and is not pooled into 'n'.",
      source_name        = "CYP2C9 genotype"
    ),
    CYP2C9_S3_COUNT = list(
      description        = "Count of CYP2C9*3 reduced-function alleles per subject (0, 1, or 2)",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "See CYP2C9_S1_COUNT. Pooled with *5, *6 and *8 into the Reyes-Gonzalez 2020 'n' allele class.",
      source_name        = "CYP2C9 genotype"
    ),
    CYP2C9_S5_COUNT = list(
      description        = "Count of CYP2C9*5 reduced-function alleles per subject (0, 1, or 2)",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "See CYP2C9_S1_COUNT. Pooled with *3, *6 and *8 into the Reyes-Gonzalez 2020 'n' allele class. Present in the Rodriguez-Fernandez 2024 cohort as *1/*5 (n=1) and *2/*5 (n=1) (Table 1).",
      source_name        = "CYP2C9 genotype"
    ),
    CYP2C9_S6_COUNT = list(
      description        = "Count of CYP2C9*6 reduced-function alleles per subject (0, 1, or 2)",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "See CYP2C9_S1_COUNT. Pooled with *3, *5 and *8 into the Reyes-Gonzalez 2020 'n' allele class. No *6 carriers were present in the Rodriguez-Fernandez 2024 cohort (Table 1); the column is retained because Reyes-Gonzalez 2020 explicitly names *6 as a member of the 'n' class, so the model is defined for *6 carriers.",
      source_name        = "CYP2C9 genotype"
    ),
    CYP2C9_S8_COUNT = list(
      description        = "Count of CYP2C9*8 reduced-function alleles per subject (0, 1, or 2)",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "See CYP2C9_S1_COUNT. Pooled with *3, *5 and *6 into the Reyes-Gonzalez 2020 'n' allele class. Present in the Rodriguez-Fernandez 2024 cohort as *1/*8 (n=2) (Table 1); *8 is an African-ancestry-enriched allele and its presence is one of the features distinguishing this Caribbean Hispanic cohort from the European cohorts of the Hamberg models.",
      source_name        = "CYP2C9 genotype"
    ),
    VKORC1_1639G_COUNT = list(
      description        = "Count of VKORC1 -1639G alleles per subject (0, 1, or 2). The complementary -1639A count is 2 - VKORC1_1639G_COUNT.",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = "No single reference level: IC50 is a per-allele sum (5.88 mg/L per G allele + 4.61 mg/L per A allele) and the baseline INR is a three-level categorical (G/A 1.78, G/G 1.84, A/A 2.18).",
      notes              = "Cohort distribution (Rodriguez-Fernandez 2024 Table 1, n=138): G/A 45%, G/G 41.3%, A/A 13.7%. The per-allele IC50 sum reproduces the genotype IC50 values quoted in the Results exactly: G/G 2*5.88 = 11.76, G/A 5.88+4.61 = 10.49, A/A 2*4.61 = 9.22 mg/L. The baseline INR effect is categorical rather than per-allele and is NOT monotone in the G count (G/A 1.78 < G/G 1.84 < A/A 2.18), so it cannot be collapsed to an allele-dosage term.",
      source_name        = "VKORC1 haplotype"
    )
  )

  # Covariates screened during the univariate covariate analysis but not
  # retained in the final model ("Other covariates, binary or (normalized)
  # continuous variables, were investigated during the modeling process, but no
  # statistical improvement (p < 0.01) after their inclusion was observed" --
  # Results, Final population PK/PD model). No point estimates are reported for
  # these, so they are documentation only and are not referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age at baseline",
      units       = "years",
      type        = "continuous",
      notes       = "Screened univariately on the PD parameters; not retained (dOFV did not reach the 6.63-unit / p<0.01 threshold). Cohort median 68 years, range 31-90 (Table 1). Note that age IS a covariate on clearance in the Hamberg warfarin models, but the PK layer here is fixed from Reyes-Gonzalez 2020, which carries no age term."
    ),
    SMOKE = list(
      description = "Current-smoker indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened; not retained. Cohort: 59 (43%) yes, 79 (57%) no (Table 1)."
    ),
    DIS_DIAB = list(
      description = "Type 2 diabetes mellitus comorbidity indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened; not retained. Cohort: 39 (28%) yes (Table 1)."
    ),
    RACE_WHITE = list(
      description = "Self-reported White Hispanic race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened; not retained. Cohort self-reported race (Table 1): White Hispanic 33 (24%), Black Hispanic 26 (19%), Admixed 37 (27%), other/not reported 42 (30%)."
    ),
    RACE_BLACK_HISPANIC = list(
      description = "Self-reported Black Hispanic race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened; not retained. See RACE_WHITE for the cohort distribution."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 138L,
    n_studies      = 1L,
    age_range      = "31-90 years",
    age_median     = "68 years",
    weight_range   = "51-159 kg",
    weight_median  = "83 kg",
    sex_female_pct = NA_real_,
    race_ethnicity = "Caribbean Hispanic (Puerto Rican), self-reported: White Hispanic 24%, Black Hispanic 19%, Admixed 27%, other/not reported 30%",
    disease_state  = "Adults on long-term warfarin anticoagulation for thromboembolic disorders",
    dose_range     = "7-82 mg total weekly warfarin (24- or 48-hour dosing interval)",
    regions        = "Puerto Rico (Veterans Affairs Caribbean Healthcare System anticoagulation clinic, San Juan)",
    cyp2c9_freq    = "*1/*1 71%, *1/*2 15%, *1/*3 5.1%, *1/*5 0.72%, *1/*8 1.44%, *2/*2 1.44%, *2/*3 4.34%, *2/*5 0.72% (Table 1)",
    vkorc1_freq    = "G/A 45%, G/G 41.3%, A/A 13.7% (Table 1)",
    n_pd_records   = "1033 INR observations",
    notes          = paste(
      "Elderly, predominantly male cohort; the paper states 'mostly males' but does not report a",
      "sex breakdown, so sex_female_pct is NA. No warfarin concentrations were measured -- the PK",
      "layer is a deterministic genotype-driven simulation carried from Reyes-Gonzalez 2020, so",
      "every subject sharing a CYP2C9 genotype and body weight has an identical concentration-time",
      "profile and no PK IIV exists. Observation period exceeded 2000 days per subject (Discussion).",
      "See Rodriguez-Fernandez 2024 Table 1."
    )
  )

  ini({
    # ==================================================================
    # PK layer -- Reyes-Gonzalez 2020 Methods (Study Design). Fixed with
    # no IIV: "Given the absence of experimental PK values in the
    # recruited patients, inter-individual variability (IIV) was not
    # associated to the PK parameters" (Rodriguez-Fernandez 2024,
    # Base population PK/PD model).
    # ==================================================================
    lka  <- fixed(log(1.19))    ; label("Warfarin first-order absorption rate constant ka (1/h)")                    # Reyes-Gonzalez 2020 Methods (fixed from the literature, their ref [14])
    lvc  <- fixed(log(0.14))    ; label("Warfarin apparent volume of distribution per kg body weight (L/kg)")        # Reyes-Gonzalez 2020 Methods (population average from the literature, their ref [5]); F assumed 100%

    # CYP2C9-diplotype-specific elimination rate constants. Reyes-Gonzalez
    # 2020 pools *3, *5, *6 and *8 into a single reduced-function allele
    # class "n", giving six diplotype classes. All six are fixed literature
    # values (their ref [15]), not estimated in either publication.
    lkel_11 <- fixed(log(0.0189)); label("Warfarin elimination rate constant kel, CYP2C9 *1/*1 (1/h)")               # Reyes-Gonzalez 2020 Methods
    lkel_12 <- fixed(log(0.0158)); label("Warfarin elimination rate constant kel, CYP2C9 *1/*2 (1/h)")               # Reyes-Gonzalez 2020 Methods
    lkel_1n <- fixed(log(0.0132)); label("Warfarin elimination rate constant kel, CYP2C9 *1/*n, n in {*3,*5,*6,*8} (1/h)")  # Reyes-Gonzalez 2020 Methods
    lkel_22 <- fixed(log(0.0130)); label("Warfarin elimination rate constant kel, CYP2C9 *2/*2 (1/h)")               # Reyes-Gonzalez 2020 Methods
    lkel_2n <- fixed(log(0.0090)); label("Warfarin elimination rate constant kel, CYP2C9 *2/*n (1/h)")               # Reyes-Gonzalez 2020 Methods
    lkel_nn <- fixed(log(0.0075)); label("Warfarin elimination rate constant kel, CYP2C9 *n/*n (1/h)")               # Reyes-Gonzalez 2020 Methods

    # ==================================================================
    # PD layer -- Rodriguez-Fernandez 2024 Table 2 and Results
    # ==================================================================
    # Mean transit times, fixed to the Hamberg 2010 values. Table 2 reports
    # them in days (1.13 and 4.62); the Results give the hour equivalents
    # (27.2 and 110.9 h), which are used here because units$time is "h".
    lmtt1 <- fixed(log(27.2))  ; label("Mean transit time MTT1 of the fast coagulation-factor chain (h; 1.13 d)")    # Table 2 (FIX); hour value from Results, Base population PK/PD model
    lmtt2 <- fixed(log(110.9)) ; label("Mean transit time MTT2 of the slow coagulation-factor chain (h; 4.62 d)")    # Table 2 (FIX); hour value from Results, Base population PK/PD model

    # Baseline INR by VKORC1 haplotype -- a three-level categorical effect.
    lrbase_ga <- log(1.78)     ; label("Typical baseline INR, VKORC1 G/A (unitless)")                                 # Table 2 (bootstrap median 1.76, 95%CI 1.48-1.91)
    lrbase_gg <- log(1.84)     ; label("Typical baseline INR, VKORC1 G/G (unitless)")                                 # Table 2 (bootstrap median 1.85, 95%CI 1.53-2.05)
    lrbase_aa <- log(2.18)     ; label("Typical baseline INR, VKORC1 A/A (unitless)")                                 # Table 2 (bootstrap median 2.14, 95%CI 1.97-2.23)

    # IC50 as a sum of per-allele contributions. Reproduces the genotype
    # values quoted in the Results: G/G 11.76, G/A 10.49, A/A 9.22 mg/L.
    lec50_g   <- log(5.88)     ; label("Warfarin IC50 contribution per VKORC1 -1639G allele (mg/L)")                  # Table 2 (bootstrap median 5.91, 95%CI 5.61-6.37)
    lec50_a   <- log(4.61)     ; label("Warfarin IC50 contribution per VKORC1 -1639A allele (mg/L)")                  # Table 2 (bootstrap median 4.63, 95%CI 4.41-4.98)

    lhill     <- log(1.47)     ; label("Sigmoidicity (Hill) coefficient gamma of the inhibitory Imax function (unitless)")  # Table 2, paper symbol gamma (bootstrap median 1.48, 95%CI 1.27-1.61)
    limax     <- fixed(log(1)) ; label("Maximum fractional inhibition Imax of vitamin K epoxide reductase (1 = 100%)")      # Results, Base population PK/PD model ("complete (100%) inhibition ... (Imax)")
    linrmax   <- fixed(log(20)); label("Maximum INR increment above baseline (unitless)")                                   # Results, Base population PK/PD model ("the maximum INR was set to 20 as previously reported")

    # IIV -- exponential, reported in Table 2 as CV%. Variances on the log
    # scale via omega^2 = log(1 + CV^2): 23% -> 0.05154, 34% -> 0.10938.
    etalrbase ~ 0.05154        # Table 2, Inter-individual variability, Baseline 23% (shrinkage 13%; bootstrap median 24%, 95%CI 19-30%)
    etalec50  ~ 0.10938        # Table 2, Inter-individual variability, IC50 34% (shrinkage 44%; bootstrap median 34%, 95%CI 28-39%)

    # Residual unexplained variability -- proportional model on INR.
    propSd    <- 0.27          ; label("Proportional residual error on INR (fraction)")                                     # Table 2, Residual unexplained variability, Proportional 27% (shrinkage 6%; bootstrap median 26%, 95%CI 22-31%)
  })

  model({
    # ---- 1. CYP2C9 diplotype class from per-allele counts ----------------
    # Reyes-Gonzalez 2020 pools *3, *5, *6 and *8 into one reduced-function
    # allele class "n"; the six indicators below are mutually exclusive and
    # sum to 1 for any subject whose six count columns sum to 2.
    cyp2c9_n_count <- CYP2C9_S3_COUNT + CYP2C9_S5_COUNT + CYP2C9_S6_COUNT + CYP2C9_S8_COUNT

    is_cyp2c9_11 <- (CYP2C9_S1_COUNT == 2)
    is_cyp2c9_12 <- (CYP2C9_S1_COUNT == 1) * (CYP2C9_S2_COUNT == 1)
    is_cyp2c9_1n <- (CYP2C9_S1_COUNT == 1) * (cyp2c9_n_count == 1)
    is_cyp2c9_22 <- (CYP2C9_S2_COUNT == 2)
    is_cyp2c9_2n <- (CYP2C9_S2_COUNT == 1) * (cyp2c9_n_count == 1)
    is_cyp2c9_nn <- (cyp2c9_n_count == 2)

    lkel_geno <- lkel_11 * is_cyp2c9_11 +
                 lkel_12 * is_cyp2c9_12 +
                 lkel_1n * is_cyp2c9_1n +
                 lkel_22 * is_cyp2c9_22 +
                 lkel_2n * is_cyp2c9_2n +
                 lkel_nn * is_cyp2c9_nn

    # ---- 2. VKORC1 haplotype indicators ----------------------------------
    is_vkorc1_gg <- (VKORC1_1639G_COUNT == 2)
    is_vkorc1_ga <- (VKORC1_1639G_COUNT == 1)
    is_vkorc1_aa <- (VKORC1_1639G_COUNT == 0)

    # ---- 3. Individual parameters ----------------------------------------
    ka  <- exp(lka)
    # Volume is strictly proportional to body weight (0.14 L/kg), so no
    # reference weight and no allometric exponent are involved.
    vc  <- exp(lvc) * WT
    kel <- exp(lkel_geno)
    cl  <- kel * vc

    mtt1   <- exp(lmtt1)
    mtt2   <- exp(lmtt2)
    hill   <- exp(lhill)
    imax   <- exp(limax)
    inrmax <- exp(linrmax)

    rbase <- exp(lrbase_gg * is_vkorc1_gg +
                 lrbase_ga * is_vkorc1_ga +
                 lrbase_aa * is_vkorc1_aa +
                 etalrbase)

    # Per-allele IC50 contributions summed on the LINEAR scale (the paper's
    # genotype IC50 values are exact sums of the two allele contributions),
    # then given exponential IIV.
    ec50 <- (exp(lec50_g) * VKORC1_1639G_COUNT +
             exp(lec50_a) * (2 - VKORC1_1639G_COUNT)) * exp(etalec50)

    # ---- 4. One-compartment first-order oral PK --------------------------
    # F is assumed complete (100%), so CL and V are apparent values
    # (Reyes-Gonzalez 2020 Methods; Rodriguez-Fernandez 2024 Results).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    Cc <- central / vc

    # ---- 5. Inhibitory sigmoid Imax on clotting-factor synthesis ---------
    inhibition <- imax * Cc^hill / (ec50^hill + Cc^hill)

    # ---- 6. Two parallel transit chains, three compartments each ---------
    # "two transit compartment chains with three compartments each to account
    # for the delay between exposure and INR response" (Results, Base
    # population PK/PD model). Each chain's zero-order input is reduced by
    # the fractional inhibition, following the Hamberg transit-chain form
    # ktr = 1 / MTT (see modellib('Hamberg_2007_warfarin_s')).
    # "_s" is the SHORT-MTT (fast) chain, "_l" the LONG-MTT (slow) chain.
    ktr1 <- 1 / mtt1
    d/dt(coag_s1) <- ktr1 * (1 - inhibition - coag_s1)
    d/dt(coag_s2) <- ktr1 * (coag_s1 - coag_s2)
    d/dt(coag_s3) <- ktr1 * (coag_s2 - coag_s3)

    ktr2 <- 1 / mtt2
    d/dt(coag_l1) <- ktr2 * (1 - inhibition - coag_l1)
    d/dt(coag_l2) <- ktr2 * (coag_l1 - coag_l2)
    d/dt(coag_l3) <- ktr2 * (coag_l2 - coag_l3)

    # Drug-free steady state of every chain compartment is unity.
    coag_s1(0) <- 1
    coag_s2(0) <- 1
    coag_s3(0) <- 1
    coag_l1(0) <- 1
    coag_l2(0) <- 1
    coag_l3(0) <- 1

    # ---- 7. INR observation ----------------------------------------------
    # Hamberg-family form: the two chain end-compartments combine as a
    # PRODUCT, and INR rises from the subject's baseline towards
    # baseline + INRmax as inhibition deepens. Rodriguez-Fernandez 2024
    # reports no INR-response exponent (lambda) in Table 2, and its Table 3
    # dosing recommendations are only reproducible with lambda = 1 (see the
    # validation vignette); the exponent is therefore absent here.
    INR <- rbase + inrmax * (1 - coag_s3 * coag_l3)

    INR ~ prop(propSd)
  })
}
