Yao_2025_flurbiprofen_r <- function() {
  description <- "Two-compartment population PK model for R(-)-flurbiprofen, the less anti-inflammatory enantiomer liberated from intravenous flurbiprofen axetil, in 67 Chinese adults undergoing elective unilateral joint replacement under spinal anaesthesia (Yao 2025 Table 3). Plasma is the central compartment and cerebrospinal fluid (CSF) is the peripheral compartment, so the paper's Vp and Q are the CSF distribution volume and the plasma-CSF intercompartmental clearance; both matrices were assayed enantioselectively. Typical values Vc = 17.0 L, CL = 11.8 L/h, Vcsf = 79.1 L, Q = 0.45 L/h. Body surface area scales Vc as a power function normalised to the cohort median 2.6 m^2, and POR rs1057868 genotype acts on CL as two genotype indicators relative to the paper's AA reference group. Parameters are apparent values relative to the nominal 100 mg flurbiprofen axetil dose."
  reference <- paste(
    "Yao H, Luo X, Yuan J, Zhang H, An H, Feng Y.",
    "Exploring the population pharmacokinetic and pharmacogenetics characteristics of",
    "flurbiprofen isomers in selective joint replacement patients with postoperative pain.",
    "Drug Des Devel Ther. 2025;19:9169-9183.",
    "doi:10.2147/DDDT.S542722.",
    "Companion model for the S(+) enantiomer: modellib('Yao_2025_flurbiprofen_s').",
    sep = " "
  )
  vignette <- "Yao_2025_flurbiprofen"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE because Yao 2025 Results
  # ("Plasma and CSF were conceptualized as the central and peripheral
  # compartments, respectively (Figure S1)") states the matrix of both
  # states explicitly, and the assay section confirms the analyte is the
  # R(-) enantiomer of flurbiprofen in both plasma and CSF.
  compartmentData <- list(
    central = list(analyte = "R(-)-flurbiprofen", units = "mg", specimen = "plasma", verified = TRUE),
    csf     = list(analyte = "R(-)-flurbiprofen", units = "mg", specimen = "CSF", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description        = "Body surface area at baseline",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters Vc as a power function normalised to the cohort median,",
        "Vc = tvVc * (BSA / 2.6)^1.37, per Yao 2025 Methods Equation I ('Continuous covariates",
        "(age, weight, height, and BMI) were incorporated using a power function after",
        "normalization to their median values') and the Table 3 row 'BSA on Vc'. The",
        "normalisation constant 2.6 m^2 is the cohort median from Table 1 (IQR 2.5-2.7, range",
        "2.0-3.2). NOTE: the paper's tabulated BSA values are not reproducible from its own",
        "height and weight (median 1.61 m and 70 kg give about 1.76 m^2 by DuBois or 1.79 m^2 by",
        "Mosteller, not 2.6 m^2), and the BSA formula is never stated. Because the exponent was",
        "estimated against the paper's own BSA scale, this model must be driven with BSA on that",
        "same scale; the vignette Errata documents the discrepancy and the virtual cohort samples",
        "BSA from the reported distribution rather than recomputing it from height and weight."
      ),
      source_name        = "BSA"
    ),
    SNP_POR_RS1057868_GA = list(
      description        = "POR rs1057868 (POR*28) heterozygous GA genotype indicator (1 = GA, 0 = otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the paper's AA reference genotype group, n = 17 of 67)",
      notes              = paste(
        "Time-fixed per subject (germline genotype). Yao 2025 Table 1 reports rs1057868 as",
        "AA/GA/GG = 17/41/9 and Table 3 estimates a separate CL coefficient for GA and for GG,",
        "so AA is the reference category. Effect enters CL multiplicatively as",
        "exp(e_snp_por_rs1057868_ga_cl) = exp(-0.29) = 0.748, i.e. GA subjects have about 25%",
        "lower apparent CL than AA subjects. POR encodes cytochrome P450 oxidoreductase, which",
        "supplies electrons to microsomal CYP enzymes. This locus deviates from Hardy-Weinberg",
        "equilibrium in this cohort (Results), which the authors attribute to the small sample",
        "size. The paper does NOT state which genotype is the wild type, so this indicator is",
        "named by the reported genotype letters rather than by wild-type/variant status; see the",
        "covariate register entry for the unresolved strand orientation."
      ),
      source_name        = "POR (rs1057868) GA"
    ),
    SNP_POR_RS1057868_GG = list(
      description        = "POR rs1057868 (POR*28) homozygous GG genotype indicator (1 = GG, 0 = otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the paper's AA reference genotype group, n = 17 of 67)",
      notes              = paste(
        "Time-fixed per subject (germline genotype). Paired with SNP_POR_RS1057868_GA; the two",
        "indicators are mutually exclusive and both are 0 for the AA reference group. Effect",
        "enters CL multiplicatively as exp(e_snp_por_rs1057868_gg_cl) = exp(-2.01) = 0.134, i.e.",
        "GG subjects have about 87% lower apparent CL than AA subjects. Yao 2025 Discussion:",
        "'the GA and GG genotypes of the POR gene significantly inhibit the clearance of",
        "R-flurbiprofen compared to the AA genotype ... may be attributed to compromised electron",
        "transfer from POR to CYP450 enzymes'. Only 9 of 67 subjects were GG, so this coefficient",
        "is imprecisely estimated (bootstrap 95% CI -5.2 to -1.2)."
      ),
      source_name        = "POR (rs1057868) GG"
    )
  )

  # Covariates screened by Yao 2025 but NOT retained in the final R(-)
  # model. Documentation only -- checkModelConventions() does not require
  # these to be referenced in model(). BMI and a second POR-on-Vp effect
  # entered during forward inclusion but were dropped in backward
  # elimination (Results, R(-) section: "the removal of BMI and POR
  # (rs1057868) genotype covariates did not result in a significant
  # increase in -2LL"), and Table 3 carries neither a BMI row nor a
  # POR-on-Vp row.
  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Median 27.1 kg/m^2 (IQR 25.1-29.4, range 18.4-40.1; Yao 2025 Table 1). Significant on the peripheral (CSF) volume of distribution during forward inclusion but removed in backward elimination; absent from Table 3."
    ),
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Median 70 kg (IQR 64-78, range 47-96; Yao 2025 Table 1). Screened; Discussion states weight showed no substantial impact on either enantiomer."
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Median 161 cm (IQR 158-164, range 140-180; Yao 2025 Table 1, reported in m). Screened; not retained."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Median 70 years (IQR 67-75, range 49-83; Yao 2025 Table 1). Screened; not retained."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "57 of 67 subjects female (85.1%; Yao 2025 Table 1). Screened; Discussion states gender showed no substantial impact."
    ),
    SNP_PXR_RS3814055 = list(
      description = "PXR (NR1I2) rs3814055 genotype",
      units       = "(binary)",
      type        = "categorical",
      notes       = "CC/CT/TT = 28/32/7 (Yao 2025 Table 1). Screened as a candidate covariate on CL and Vd; not retained."
    ),
    SNP_PXR_RS1523127 = list(
      description = "PXR (NR1I2) rs1523127 genotype",
      units       = "(binary)",
      type        = "categorical",
      notes       = "AA/CA/CC = 3/60/4 (Yao 2025 Table 1); deviates from Hardy-Weinberg equilibrium in this cohort. Screened; not retained."
    ),
    SNP_POR_RS1135612 = list(
      description = "POR rs1135612 genotype",
      units       = "(binary)",
      type        = "categorical",
      notes       = "CC/CT/TT = 7/33/27 (Yao 2025 Table 1). Screened; not retained."
    ),
    SNP_ABCB1_RS1045642 = list(
      description = "ABCB1 rs1045642 (c.3435C>T) genotype",
      units       = "(binary)",
      type        = "categorical",
      notes       = "AA/GA/GG = 12/40/15 (Yao 2025 Table 1). Retained on CL in the companion S(+) model but not in the R(-) model."
    ),
    SNP_ABCB1_RS4148738 = list(
      description = "ABCB1 rs4148738 genotype",
      units       = "(binary)",
      type        = "categorical",
      notes       = "AA/GA/GG = 15/37/15 (Yao 2025 Table 1). Screened; not retained."
    ),
    SNP_CYP2C9_MONOMORPHIC = list(
      description = "CYP2C9 rs1057910 / rs1799853 / rs182132442 genotypes",
      units       = "(binary)",
      type        = "categorical",
      notes       = "All 67 subjects were homozygous wild type at every assayed CYP2C9 locus (Yao 2025 Table 1 and Results: 'All patients were homozygous for the wild-type CYP2C9 allele'), so the locus was monomorphic and could not be evaluated even though CYP2C9 is the principal flurbiprofen-metabolising enzyme."
    ),
    SNP_UGT1A9_RS28898617 = list(
      description = "UGT1A9 rs28898617 genotype",
      units       = "(binary)",
      type        = "categorical",
      notes       = "All 67 subjects TT (Yao 2025 Table 1); monomorphic, excluded from analysis per Results."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 67L,
    n_studies      = 1L,
    age_range      = "49-83 years",
    age_median     = "70 years (IQR 67-75; mean 71)",
    weight_range   = "47-96 kg",
    weight_median  = "70 kg (IQR 64-78; mean 71)",
    height_range   = "1.40-1.80 m",
    height_median  = "1.61 m (IQR 1.58-1.64; mean 1.62)",
    bmi_range      = "18.4-40.1 kg/m^2",
    bmi_median     = "27.1 kg/m^2 (IQR 25.1-29.4; mean 27.3)",
    bsa_range      = "2.0-3.2 m^2",
    bsa_median     = "2.6 m^2 (IQR 2.5-2.7; mean 2.6)",
    sex_female_pct = 85.1,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Adults aged 18-85 years, ASA physical status I-III, scheduled for elective unilateral joint replacement surgery under spinal (subarachnoid) anaesthesia, treated for postoperative pain. Patients with asthma, hepatic or renal dysfunction, peptic ulcer, NSAID allergy, recent NSAID or CYP2C9 inhibitor/inducer exposure, or subnormal plasma total protein or albumin were excluded. Ethics approval #2019PHB169-01 (Peking University People's Hospital); ClinicalTrials.gov NCT04128410.",
    dose_range     = "Single 100 mg intravenous injection of flurbiprofen axetil (FEX, 5050E; Beijing Tide Pharmaceutical) infused at 2 mL/min, given after 1 mg intravenous midazolam and before subarachnoid ropivacaine. Flurbiprofen axetil is a lipid-microsphere prodrug hydrolysed by carboxylesterase to flurbiprofen within about 5 minutes.",
    regions        = "China (single centre, Peking University People's Hospital, Beijing).",
    n_observations = "134 concentrations (67 plasma + 67 CSF) from 67 patients. Medical-ethics constraints allowed only ONE CSF sample per participant, so patients were block-randomised into 10 groups of about 7 and each group was sampled at a single nominal post-dose time (5, 10, 15, 20, 25, 30, 35, 40, 45 or 50 min); the paired venous plasma sample was drawn simultaneously from the contralateral median cubital vein. Assay LLOQ 0.1 ug/mL (plasma, linear 0.1-10 ug/mL) and 1 ng/mL (CSF, linear 1-100 ng/mL) by enantioselective LC-MS/MS on a CHIRALPAK-IG3 column.",
    notes          = "70 patients enrolled October 2019 to June 2020; 67 analysed (3 excluded for undetectable genotypes). Phoenix NLME 8.3 with FOCE-ELS. Validated by 1000-replicate bootstrap and visual predictive check. The extremely sparse one-sample-per-subject design means the reported inter-individual variances are weakly identified. Concomitant midazolam and ropivacaine were not modelled."
  )

  ini({
    # ================================================================
    # Structural PK parameters -- Yao 2025 Table 3 ("Pharmacokinetic
    # Parameter Estimates and Bootstrap Results of Final PPK Model for
    # R(-)-Flurbiprofen"), Estimate column. Plasma is the central
    # compartment and CSF the peripheral compartment (Results,
    # Pharmacokinetic Study), so the paper's Vp is the CSF volume
    # (encoded lvcsf, following the Kumpulainen_2010_flurbiprofen
    # precedent) and the paper's Q is the plasma-CSF intercompartmental
    # clearance.
    #
    # These are APPARENT values relative to the nominal 100 mg
    # flurbiprofen axetil dose -- Table 3's own abbreviation footnote
    # reads "CL, apparent clearance; VC, apparent volume of
    # distribution in the central compartment" even though dosing was
    # intravenous. The paper never states the amount entered in the
    # dataset; the vignette Errata documents the 100 mg basis and the
    # cross-check against Zhang 2018, which used the same product and
    # study design at the same centre.
    # ================================================================
    lvc   <- log(17.0); label("Apparent central (plasma) volume of distribution Vc (L)")      # Table 3: Vc = 17.0 L (bootstrap median 17.3, 95% CI 15.0-20.1). The abstract quotes 17.1 L for the same parameter; Table 3 is the designated final-model parameter table and is used here (see vignette Errata).
    lcl   <- log(11.8); label("Apparent plasma clearance CL (L/h)")                           # Table 3: CL = 11.8 L/h (bootstrap median 8.7, 95% CI 2.92-15.54); Discussion reports the same estimate to more digits as 11.76 L/h
    lvcsf <- log(79.1); label("Apparent peripheral (CSF) volume of distribution Vp (L)")      # Table 3: Vp = 79.1 L (bootstrap median 55.3, 95% CI 33.5-58.3 -- the published CI does not contain the point estimate; see vignette Errata); abstract: "the Vd of R(-)-flurbiprofen (CSF 79.1 L VS plasma 17.1 L)"
    lq    <- log(0.45); label("Plasma-CSF intercompartmental clearance Q (L/h)")              # Table 3: Q = 0.45 L/h (bootstrap median 0.32, 95% CI 0.19 to 0.92); abstract restates this as the "CSF CL" of R(-)-flurbiprofen (0.45 L/h)

    # ================================================================
    # Covariate effect on Vc -- body surface area, power function
    # normalised to the cohort median.
    # Methods, Population Covariate Analysis: "Continuous covariates
    # (age, weight, height, and BMI) were incorporated using a power
    # function after normalization to their median values, as specified
    # in Equation I", i.e. Vc = tvVc * (BSA / median(BSA))^exponent.
    # Median BSA = 2.6 m^2 (Table 1).
    # Sign check against the Discussion: "A larger BSA is associated
    # with a greater apparent volume of distribution for
    # R-flurbiprofen" -- a positive exponent. Consistent.
    # ================================================================
    e_bsa_vc <- 1.37; label("Power exponent on (BSA / 2.6 m^2) for Vc (unitless)")            # Table 3: BSA on Vc = 1.37 (bootstrap median 0.59, 95% CI -0.07 to 1.96 -- the bootstrap CI spans zero; see vignette Errata)

    # ================================================================
    # Covariate effects on CL -- POR rs1057868 (POR*28) genotype.
    # Methods, Population Covariate Analysis: "Categorical covariates,
    # including sex and genetic polymorphisms, were modeled using
    # indicator variables within an exponential function (Equation II),
    # with the reference category coded as 0 and other categories
    # assigned values of 1 or higher." Table 3 lists one coefficient
    # for GA and one for GG, so AA (n = 17) is the reference group and
    # the two indicators enter a single exponential:
    #   CL = tvCL * exp(e_GA * I(GA) + e_GG * I(GG)) * exp(eta)
    # Sign check against the Discussion: "the GA and GG genotypes of
    # the POR gene significantly inhibit the clearance of
    # R-flurbiprofen compared to the AA genotype" -- exp(-0.29) = 0.748
    # (25% lower) and exp(-2.01) = 0.134 (87% lower). Both negative,
    # monotonic in G-allele count. Consistent.
    # ================================================================
    e_snp_por_rs1057868_ga_cl <- -0.29; label("Exponential coefficient on CL for POR rs1057868 GA vs AA (unitless)") # Table 3: POR (rs1057868) on CL, GA = -0.29 (bootstrap median -0.15, 95% CI -0.42 to -0.05)
    e_snp_por_rs1057868_gg_cl <- -2.01; label("Exponential coefficient on CL for POR rs1057868 GG vs AA (unitless)") # Table 3: POR (rs1057868) on CL, GG = -2.01 (bootstrap median -1.42, 95% CI -5.2 to -1.2)

    # ================================================================
    # Inter-individual variability -- Yao 2025 Table 3, "Interindividual
    # variability (%)" block. The table's own abbreviation footnote
    # defines the reported quantity as "omega, square root of
    # interindividual variance for parameters", so each tabulated value
    # is the SD of eta on the log scale and the internal nlmixr2
    # variance is omega^2. (The "(%)" in the block header is
    # inconsistent with the footnote; taking the values as CV fractions
    # instead would give omega^2 = log(1 + CV^2), which differs by less
    # than 3% for every value here, so the distinction is immaterial.)
    # The exponential IIV model is stated in Methods, Structural Model.
    # No correlations between etas were reported.
    # ================================================================
    etalvc   ~ 0.0009  # Table 3: omega_Vc = 0.03 -> variance 0.03^2 = 0.0009
    etalcl   ~ 0.0144  # Table 3: omega_CL = 0.12 -> variance 0.12^2 = 0.0144
    etalvcsf ~ 0.0484  # Table 3: omega_Vp = 0.22 -> variance 0.22^2 = 0.0484
    etalq    ~ 0.0196  # Table 3: omega_Q  = 0.14 -> variance 0.14^2 = 0.0196

    # ================================================================
    # Residual error -- Yao 2025 Table 3, "Random residual variability"
    # block: one additive error on the plasma observations and one
    # multiplicative (proportional) error on the CSF observations,
    # identical to the S(+) model's Table 2 values.
    #
    # The R(-) Results section describes the structural model as
    # adopting "a proportional residual error model", which conflicts
    # with Table 3's "Plasma, additive error" row. Table 3 is the
    # designated final-model parameter table and is followed here.
    #
    # Both magnitudes are implausibly small relative to the observed
    # concentrations (plasma is a few mg/L; an additive SD of 0.003 mg/L
    # = 3 ng/mL sits far below the 0.1 ug/mL plasma LLOQ, and a 0.1%
    # proportional CSF error is far tighter than any bioanalytical
    # assay). They are transcribed exactly as published. The same
    # pattern appears in Zhang 2018 (sigma = 0.0023 mg/L additive) from
    # the same group and centre, so it reflects this group's Phoenix
    # reporting convention rather than a transcription error here.
    # ================================================================
    addSd        <- 0.003; label("Additive residual error SD on plasma concentration Cc (mg/L)")   # Table 3: "Plasma, additive error, sigma" = 0.003
    propSd_Ccsf  <- 0.001; label("Proportional residual error SD on CSF concentration Ccsf (fraction)") # Table 3: "CSF, multiplicative error, sigma" = 0.001
  })

  model({
    # ---- Individual PK parameters ---------------------------------
    # BSA scales Vc as a power function normalised to the cohort median
    # 2.6 m^2 (Yao 2025 Methods Equation I; Table 3 "BSA on Vc").
    # POR rs1057868 genotype acts on CL through two mutually exclusive
    # indicators, both 0 for the AA reference group (Methods Equation
    # II; Table 3).
    cl   <- exp(lcl + etalcl) *
            exp(e_snp_por_rs1057868_ga_cl * SNP_POR_RS1057868_GA +
                e_snp_por_rs1057868_gg_cl * SNP_POR_RS1057868_GG)
    vc   <- exp(lvc + etalvc) * (BSA / 2.6)^e_bsa_vc
    vcsf <- exp(lvcsf + etalvcsf)
    q    <- exp(lq    + etalq)

    # ---- Concentrations -------------------------------------------
    Cc   <- central / vc    # R(-)-flurbiprofen in plasma (mg/L)
    Ccsf <- csf / vcsf      # R(-)-flurbiprofen in CSF (mg/L)

    # ---- ODE system ------------------------------------------------
    # Standard two-compartment IV disposition with plasma as the
    # central compartment and CSF as the peripheral compartment
    # (Results: "Plasma and CSF were conceptualized as the central and
    # peripheral compartments, respectively (Figure S1)"). Elimination
    # is first-order from plasma only; the paper reports no separate
    # CSF elimination pathway -- the abstract's "CSF CL" is the
    # intercompartmental clearance Q of Table 3.
    d/dt(central) <- -cl * Cc - q * (Cc - Ccsf)
    d/dt(csf)     <-            q * (Cc - Ccsf)

    # ---- Residual error --------------------------------------------
    Cc   ~ add(addSd)
    Ccsf ~ prop(propSd_Ccsf)
  })
}
