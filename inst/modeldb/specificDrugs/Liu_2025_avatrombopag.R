Liu_2025_avatrombopag <- function() {
  description <- paste(
    "Population PK model for avatrombopag in 92 healthy Chinese adults",
    "given a single 20 mg oral dose under fasting or high-fat-fed",
    "conditions in a two-sequence, four-period replicate bioequivalence",
    "trial. One-compartment disposition with first-order elimination",
    "preceded by a nine-transit-compartment absorption chain (dose ->",
    "depot -> transit1 .. transit9, transfer rate ktr = 9 / MTT)",
    "followed by first-order absorption (ka) from transit9 into the",
    "central compartment. A food-by-genotype co-effect acts on the",
    "absorption fraction Fa (relative bioavailability): Fa = 1 +",
    "COV1 * FED_HIGHFAT + COV2 * (1 - ABCB1_C1236T_HET) + COV3 *",
    "FED_HIGHFAT * (1 - ABCB1_C1236T_HET), so the reference Fa = 1 is a",
    "fasting ABCB1 (C1236T) TC heterozygote. Between-subject variability",
    "is estimated on CL/F, Vd/F, ka and MTT; inter-occasion variability",
    "is estimated on ka across the four study periods. The source paper",
    "additionally simulated platelet counts by linking this PK model to",
    "a previously published chronic-liver-disease platelet-count PK/PD",
    "life-span model (Nomoto 2018, J Clin Pharmacol 58:1629-1638); no",
    "PD parameter of that model is reported in Liu 2025, so only the PK",
    "layer is implemented here."
  )
  reference <- paste(
    "Liu X, Chen L, Ju G, Li C, Liu B, Fei Y, Wang X, Gao Y, He Q,",
    "Zhu X, Ouyang D. (2025). Investigation of the ABCB1 Gene",
    "Polymorphism and Food Effects on the Avatrombopag Pharmacokinetics",
    "in Chinese Individuals: A Population Pharmacokinetic/",
    "Pharmacodynamic Analysis. Pharmaceuticals 18(6):903.",
    "doi:10.3390/ph18060903"
  )
  vignette <- "Liu_2025_avatrombopag"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list(
    FED_HIGHFAT = list(
      description        = paste(
        "1 = the 20 mg avatrombopag dose was taken about 30 min after a",
        "high-fat, high-calorie meal; 0 = the dose was taken after an",
        "overnight fast."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (overnight fast)",
      notes              = paste(
        "Liu 2025 Methods 'Study Design and Pharmacokinetic Assessments':",
        "the fed condition is the FDA-standard high-fat, high-calorie",
        "breakfast providing 800-1000 kcal total (about 150 kcal protein,",
        "250 kcal carbohydrate, 500-600 kcal fat) consumed about 30 min",
        "before dosing. Time-fixed per subject in this study: the fasting",
        "arm (n = 47) and the fed arm (n = 45) are disjoint subject groups",
        "(Table 1), so FED_HIGHFAT does not vary within a subject even",
        "though the design is a four-period replicate crossover.",
        "Enters the model only through the absorption fraction Fa",
        "(Table 3 footnote, COV1 and COV3); the paper found no food effect",
        "on CL/F or Vd/F typical values."
      ),
      source_name        = "FOOD"
    ),
    ABCB1_C1236T_HET = list(
      description        = paste(
        "1 = subject is an ABCB1 (C1236T, rs1128503) TC heterozygote;",
        "0 = subject is homozygous at that locus (CC wild-type or TT",
        "variant, pooled)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = paste(
        "The Liu 2025 Fa reference group is ABCB1_C1236T_HET = 1 (TC",
        "heterozygote) combined with FED_HIGHFAT = 0, for which Fa = 1."
      ),
      notes              = paste(
        "Liu 2025 Table 3 footnote defines its model covariate as",
        "'ABCB1 = 0 if the genotype is heterozygote, 1 if homozygote',",
        "i.e. the exact complement of the canonical heterozygote",
        "indicator. The model body therefore uses (1 - ABCB1_C1236T_HET)",
        "wherever the paper writes ABCB1. Pooling CC and TT into a single",
        "'non-TC' stratum is the paper's own finding (Results 2.2.2: CC",
        "and TT showed similar exposure, Cmax ratio 102.4% fasting /",
        "113.7% fed and AUC0-t ratio 108.7% fasting / 108.4% fed, while",
        "TC differed markedly from both). Genotype frequencies in the",
        "92-subject cohort were 13.0% CC, 46.7% TC, 40.2% TT",
        "(Discussion paragraph 4). Time-fixed per subject (germline",
        "genotype)."
      ),
      source_name        = "ABCB1"
    ),
    OCC = list(
      description        = paste(
        "Integer study period (1-4) of the two-sequence, four-period",
        "replicate bioequivalence design; indexes the inter-occasion",
        "variability on ka."
      ),
      units              = "(integer)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Liu 2025 Methods: 'This open-label, crossover study consisted of",
        "four periods with two sequences, separated by a 14-day washout.'",
        "Each subject therefore contributes four single-dose occasions.",
        "Table 3 reports a single OCC [KA] variance that is shared across",
        "occasions (the NONMEM $OMEGA BLOCK(1) SAME idiom); nlmixr2 has no",
        "SAME shortcut, so four etas are declared and the variance is",
        "fixed to the shared value for occasions 2-4. For single-occasion",
        "simulations pass OCC = 1 so the first IOV eta applies."
      ),
      source_name        = "OCC"
    )
  )

  # Screened in the Liu 2025 stepwise covariate search but not retained in
  # the final model (Methods 4.4: 'Covariate analysis included demographic
  # data (e.g., age, weight, sex), lab test results (e.g., ALT, AST),
  # clinical trial conditions (e.g., food intake), and genotype
  # information'; only the food-by-ABCB1(C1236T) co-effect on Fa survived
  # backward elimination at p < 0.01). Documentation only -- none of these
  # is referenced in model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened but not retained. Liu 2025 Limitations explicitly",
        "attributes this to the homogeneous bioequivalence-study",
        "population: 'the analysis population was derived from a",
        "bioequivalence study involving relatively homogeneous",
        "participants, which may limit the ability to identify certain",
        "covariates, such as body weight, that were identified in another",
        "study'. Cohort means were 63.22 +/- 7.88 kg (fasting) and",
        "60.95 +/- 6.32 kg (fed) (Table 1)."
      )
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened but not retained; median 28 years, range 18-45 (Results 2.1, Table 1)."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened but not retained. Liu 2025 Results 2.2.1 and Table S2:",
        "'no statistically significant variation in AVA exposure was",
        "observed between males and females'. The cohort was 83 male /",
        "9 female."
      )
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained (Methods 4.4); Table 1 means 23.04 +/- 12.35 (fasting) and 19.46 +/- 12.31 U/L (fed)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained (Methods 4.4); Table 1 means 21.2 +/- 4.71 (fasting) and 18.97 +/- 4.5 U/L (fed)."
    ),
    CYP2C9_PM_IM = list(
      description = "Pooled CYP2C9 poor-or-intermediate-metabolizer phenotype indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Genotyped (*1/*3 = intermediate metabolizer, n = 6; no *3/*3",
        "poor metabolizers) and shown by NCA to raise fed-state exposure",
        "1.70-fold (Abstract; Results 2.2.2), but NOT retained as a",
        "covariate in the final population PK model (Table 3 carries only",
        "COV1-COV3, all food-by-ABCB1 terms). Liu 2025 Limitations",
        "attributes the non-retention to the small IM sample size and the",
        "directionally inconsistent fasting-state signal."
      )
    ),
    ABCB1_C3435T_HET = list(
      description = "ABCB1 C3435T (rs1045642) heterozygote indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Genotyped (40% CC / 47% TC / 13% TT) and screened; Results 2.2.2",
        "reports 'no significant differences were observed in Cmax and AUC",
        "values across the ABCB1 (C3435T) genotypes' and it was not",
        "retained in the final model."
      )
    ),
    ABCB1_G2677TA = list(
      description = "ABCB1 G2677T/A (rs2032582) genotype",
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste(
        "Genotyped (GA/GG/TA/TG/TT strata, Table 1) and screened;",
        "Results 2.2.2 reports 'no significant associations were found",
        "between the ABCB1 (G2677T/A) polymorphism and Cmax or AUC",
        "parameters' and it was not retained in the final model. No",
        "canonical covariate column is registered for this SNP because no",
        "nlmixr2lib model uses it in model()."
      )
    )
  )

  compartmentData <- list(
    depot    = list(analyte = "avatrombopag", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "avatrombopag", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "avatrombopag", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "avatrombopag", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "avatrombopag", units = "mg", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "avatrombopag", units = "mg", specimen = "administration site", verified = TRUE),
    transit6 = list(analyte = "avatrombopag", units = "mg", specimen = "administration site", verified = TRUE),
    transit7 = list(analyte = "avatrombopag", units = "mg", specimen = "administration site", verified = TRUE),
    transit8 = list(analyte = "avatrombopag", units = "mg", specimen = "administration site", verified = TRUE),
    transit9 = list(analyte = "avatrombopag", units = "mg", specimen = "administration site", verified = TRUE),
    central  = list(analyte = "avatrombopag", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species           = "human",
    n_subjects        = 92L,
    n_studies         = 1L,
    n_pk_observations = 5923L,
    age_range         = "18-45 years",
    age_median        = "28 years",
    weight_range      = "mean 63.22 +/- 7.88 kg (fasting arm) and 60.95 +/- 6.32 kg (fed arm)",
    sex_female_pct    = 9.8,
    race_ethnicity    = "Chinese (single-centre Chinese cohort; no further ethnic breakdown reported)",
    disease_state     = "healthy volunteers",
    dose_range        = paste(
      "single 20 mg oral avatrombopag tablet with 200 mL water, repeated",
      "over four periods separated by a 14-day washout; 47 subjects dosed",
      "under an overnight fast and 45 dosed about 30 min after a high-fat,",
      "high-calorie meal"
    ),
    regions           = "China (Clinical Trials Unit of Hunan Province People's Hospital, Changsha, Hunan)",
    genotypes         = paste(
      "CYP2C9: 90% *1/*1 normal metabolizer, 10% *1/*3 intermediate",
      "metabolizer, no *3/*3. ABCB1 (C1236T): 14% CC, 51% TC, 35% TT.",
      "ABCB1 (C3435T): 40% CC, 47% TC, 13% TT. ABCB1 (G2677T/A): 20%",
      "GG, 67% heterozygous, 13% homozygous variant. All variants were in",
      "Hardy-Weinberg equilibrium (Results 2.1, 2.2.2 and Table 1)."
    ),
    notes             = paste(
      "Open-label, single-centre, two-sequence, four-period replicate",
      "bioequivalence trial (IRB approval [2023]-32.1, Hunan Province",
      "People's Hospital). Sampling at predose and 1, 2, 3, 4, 5, 6, 7,",
      "8, 9, 10, 11, 12, 16, 24, 48, 72 and 96 h post dose; plasma",
      "avatrombopag by UPLC-MS/MS with an AVA-D8 internal standard,",
      "linear range 1.000-400.000 ng/mL. Estimation in NONMEM 7.5 with",
      "FOCE-I; stepwise covariate modelling (forward p < 0.05, backward",
      "p < 0.01); evaluated by GOF plots, VPC and non-parametric",
      "bootstrap (Table 3). Baseline demographics are Liu 2025 Table 1."
    )
  )

  ini({
    # ---- Structural PK parameters (Liu 2025 Table 3, 'Final Model' estimates) ----
    lcl  <- log(7.85);  label("Apparent clearance CL/F (L/h)")                                      # Table 3 'CL/F (L/h) 7.85' (RSE 14%; bootstrap median 7.81 [5.92-10.62])
    lvc  <- log(199);   label("Apparent central volume of distribution Vd/F (L)")                   # Table 3 'Vd/F (L) 199' (RSE 13%; bootstrap median 198.01 [149.77-267.64])
    lka  <- log(0.582); label("First-order absorption rate constant ka from transit9 to central (1/h)")  # Table 3 'ka (/h) 0.582' (RSE 13%; bootstrap median 0.581 [0.4723-0.785])
    lmtt <- log(1.33);  label("Mean transit time through the nine-compartment absorption chain (h); ktr = 9 / MTT")  # Table 3 'MTT (h) 1.33' (RSE 5%; bootstrap median 1.33 [1.22-1.45])

    # Reference absorption fraction. Liu 2025 Table 3 footnote writes the
    # covariate model as "Fa = 1 x (1 + COV1 x FOOD + COV2 x ABCB1 +
    # COV3 x FOOD x ABCB1)", so the leading Fa is structurally 1, not an
    # estimated theta.
    lfdepot <- fixed(log(1)); label("Reference absorption fraction Fa for a fasting ABCB1 (C1236T) TC heterozygote (unitless)")  # Table 3 footnote: "Fa = 1 x (1 + ...)"

    # ---- Food-by-genotype co-effect on the absorption fraction Fa ----
    # Table 3 footnote defines FOOD = 1 fed / 0 fasting and ABCB1 = 0
    # heterozygote / 1 homozygote, i.e. ABCB1 = 1 - ABCB1_C1236T_HET.
    # Resulting Fa by stratum: fasting TC 1.000, fed TC 0.9347,
    # fasting non-TC 0.7280, fed non-TC 0.8397.
    e_fed_highfat_fdepot  <- -0.0653; label("COV1: linear effect of the high-fat fed state on Fa in ABCB1 (C1236T) TC heterozygotes (fraction)")            # Table 3 'COV1 * -0.0653' (RSE 111%; bootstrap median -0.0702 [-0.3069-0.2698])
    e_abcb1hom_fdepot     <- -0.272;  label("COV2: linear effect of ABCB1 (C1236T) homozygosity (CC or TT) on Fa under fasting (fraction)")                 # Table 3 'COV2 * -0.272' (RSE 15%; bootstrap median -0.28 [-0.47-0.01])
    e_fed_abcb1hom_fdepot <-  0.177;  label("COV3: linear food-by-genotype interaction on Fa (high-fat fed x ABCB1 (C1236T) homozygous, fraction)")         # Table 3 'COV3 * 0.177' (RSE 56%; bootstrap median 0.186 [-0.216-0.431])

    # ---- Between-subject variability (Liu 2025 Table 3, 'Between-subject variability') ----
    # Reported on the variance scale for an exponential IIV model
    # (Methods 4.4: "Inter-individual variability (IIV) was modeled using
    # an exponential approach").
    etalcl  ~ 0.182   # Table 3 BSV 'CL/F 0.182'  (RSE 35%); sqrt(exp(0.182)-1) = 44.7% CV
    etalvc  ~ 0.159   # Table 3 BSV 'Vd/F 0.159'  (RSE 32%); sqrt(exp(0.159)-1) = 41.6% CV
    etalka  ~ 0.504   # Table 3 BSV 'KA 0.504'    (RSE 35%); sqrt(exp(0.504)-1) = 81.6% CV
    etalmtt ~ 0.131   # Table 3 BSV 'MTT 0.131'   (RSE 23%); sqrt(exp(0.131)-1) = 37.5% CV

    # ---- Inter-occasion variability on ka across the four study periods ----
    # Table 3 reports a single 'OCC [KA] 0.913' variance shared by every
    # occasion; nlmixr2 has no NONMEM SAME shortcut, so occasions 2-4 are
    # fixed to the occasion-1 value.
    etaiov_ka_1 ~ 0.913        # Table 3 'OCC [KA] 0.913' (RSE 12%); sqrt(exp(0.913)-1) = 122% CV
    etaiov_ka_2 ~ fix(0.913)   # shared with occasion 1
    etaiov_ka_3 ~ fix(0.913)   # shared with occasion 1
    etaiov_ka_4 ~ fix(0.913)   # shared with occasion 1

    # ---- Residual unexplained variability (Liu 2025 Table 3) ----
    propSd <- 0.239; label("Proportional residual error (fraction)")   # Table 3 'sigma prop (%) 23.90%' (RSE 4%; bootstrap median 23.79% [22.28-25.74%])
    addSd  <- 1.18;  label("Additive residual error (ng/mL)")          # Table 3 'sigma add (ng/mL) 1.18' (RSE 3%; bootstrap median 1.17 [0.86-1.48])
  })

  model({
    # 1. Occasion indicators and the inter-occasion variability on ka.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 +
              oc3 * etaiov_ka_3 + oc4 * etaiov_ka_4

    # 2. Individual parameters.
    cl  <- exp(lcl  + etalcl)
    vc  <- exp(lvc  + etalvc)
    ka  <- exp(lka  + etalka + iov_ka)
    mtt <- exp(lmtt + etalmtt)

    # Absorption fraction. Liu 2025 Table 3 footnote:
    #   Fa = 1 x (1 + COV1 x FOOD + COV2 x ABCB1 + COV3 x FOOD x ABCB1)
    # with FOOD = 1 fed / 0 fasting and ABCB1 = 1 homozygote / 0
    # heterozygote. The canonical column ABCB1_C1236T_HET has the opposite
    # orientation (1 = TC heterozygote), hence the (1 - ...) complement.
    fdepot <- exp(lfdepot) *
      (1 + e_fed_highfat_fdepot  * FED_HIGHFAT +
           e_abcb1hom_fdepot     * (1 - ABCB1_C1236T_HET) +
           e_fed_abcb1hom_fdepot * FED_HIGHFAT * (1 - ABCB1_C1236T_HET))

    # 3. Micro-constants. ktr = n_transit / MTT for the nine-compartment
    #    chain (dose -> depot -> transit1 .. transit9 is nine ktr steps).
    ktr <- 9 / mtt
    kel <- cl / vc

    # 4. ODE system: nine-transit absorption chain, first-order absorption
    #    from transit9 into central, first-order elimination.
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2
    d/dt(transit3) <-  ktr * transit2 - ktr * transit3
    d/dt(transit4) <-  ktr * transit3 - ktr * transit4
    d/dt(transit5) <-  ktr * transit4 - ktr * transit5
    d/dt(transit6) <-  ktr * transit5 - ktr * transit6
    d/dt(transit7) <-  ktr * transit6 - ktr * transit7
    d/dt(transit8) <-  ktr * transit7 - ktr * transit8
    d/dt(transit9) <-  ktr * transit8 - ka  * transit9
    d/dt(central)  <-  ka  * transit9 - kel * central

    # 5. Bioavailability.
    f(depot) <- fdepot

    # 6. Observation and error. Doses are in mg and vc is in L, so
    #    central / vc is mg/L; x 1000 converts to the ng/mL scale on which
    #    the assay (linear range 1.000-400.000 ng/mL) and the additive
    #    residual error (1.18 ng/mL) are reported.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
