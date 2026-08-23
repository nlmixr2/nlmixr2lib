Yao_2025_flurbiprofen_s <- function() {
  description <- "Two-compartment population PK model for S(+)-flurbiprofen, the pharmacologically more active enantiomer liberated from intravenous flurbiprofen axetil, in 67 Chinese adults undergoing elective unilateral joint replacement under spinal anaesthesia (Yao 2025 Table 2). Plasma is the central compartment and cerebrospinal fluid (CSF) is the peripheral compartment, so the paper's Vp and Q are the CSF distribution volume and the plasma-CSF intercompartmental clearance; both matrices were assayed enantioselectively. Typical values Vc = 25.6 L, CL = 16.7 L/h, Vcsf = 32.6 L, Q = 0.39 L/h. ABCB1 rs1045642 genotype is the only retained covariate and acts on CL as two genotype indicators relative to the paper's AA reference group. Parameters are apparent values relative to the nominal 100 mg flurbiprofen axetil dose."
  reference <- paste(
    "Yao H, Luo X, Yuan J, Zhang H, An H, Feng Y.",
    "Exploring the population pharmacokinetic and pharmacogenetics characteristics of",
    "flurbiprofen isomers in selective joint replacement patients with postoperative pain.",
    "Drug Des Devel Ther. 2025;19:9169-9183.",
    "doi:10.2147/DDDT.S542722.",
    "Companion model for the R(-) enantiomer: modellib('Yao_2025_flurbiprofen_r').",
    sep = " "
  )
  vignette <- "Yao_2025_flurbiprofen"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE because Yao 2025 Results
  # ("Plasma and CSF were conceptualized as the central and peripheral
  # compartments, respectively (Figure S1)") states the matrix of both
  # states explicitly, and the assay section confirms the analyte is the
  # S(+) enantiomer of flurbiprofen in both plasma and CSF.
  compartmentData <- list(
    central = list(analyte = "S(+)-flurbiprofen", units = "mg", specimen = "plasma", verified = TRUE),
    csf     = list(analyte = "S(+)-flurbiprofen", units = "mg", specimen = "CSF", verified = TRUE)
  )

  covariateData <- list(
    SNP_ABCB1_RS1045642_GA = list(
      description        = "ABCB1 rs1045642 heterozygous GA genotype indicator (1 = GA, 0 = otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the paper's AA reference genotype group, n = 12 of 67)",
      notes              = paste(
        "Time-fixed per subject (germline genotype). Yao 2025 Table 1 reports rs1045642 as",
        "AA/GA/GG = 12/40/15 and Table 2 estimates a separate CL coefficient for GA and for GG,",
        "so AA is the reference category. Effect enters CL multiplicatively as",
        "exp(e_snp_abcb1_rs1045642_ga_cl) = exp(-1.52) = 0.219, i.e. GA subjects have about 78%",
        "lower apparent CL than AA subjects. The paper does NOT state which genotype is the",
        "wild type, so this indicator is named by the reported genotype letters rather than by",
        "wild-type/variant status; see the covariate register entry for the unresolved strand",
        "orientation. Not directly composable with the pooled carrier indicator",
        "SNP_ABCB1_RS1045642, whose reference is the c.3435CC wild type."
      ),
      source_name        = "ABCB1 (rs1045642) GA"
    ),
    SNP_ABCB1_RS1045642_GG = list(
      description        = "ABCB1 rs1045642 homozygous GG genotype indicator (1 = GG, 0 = otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the paper's AA reference genotype group, n = 12 of 67)",
      notes              = paste(
        "Time-fixed per subject (germline genotype). Paired with SNP_ABCB1_RS1045642_GA; the two",
        "indicators are mutually exclusive and both are 0 for the AA reference group. Effect enters",
        "CL multiplicatively as exp(e_snp_abcb1_rs1045642_gg_cl) = exp(0.19) = 1.21, i.e. GG",
        "subjects have about 21% higher apparent CL than AA subjects. Yao 2025 Discussion:",
        "'the GA genotype of the ABCB1 gene inhibits the clearance of S-flurbiprofen, while the",
        "GG genotype increases its clearance rate'. The non-monotonic GA/GG pattern is reported",
        "as-is; the AA reference group had only 12 subjects."
      ),
      source_name        = "ABCB1 (rs1045642) GG"
    )
  )

  # Covariates screened by Yao 2025 but NOT retained in the final S(+) model.
  # Documentation only -- checkModelConventions() does not require these to be
  # referenced in model(). BSA was significant on the volume of distribution
  # during forward inclusion but was dropped in backward elimination
  # (Results, S(+) section: "Backward elimination confirmed the necessity of
  # ABCB1 genotype (rs1045642) retention whereas BSA elimination failed to
  # demonstrate statistical significance"), and Table 2 carries no BSA row.
  # The Discussion confirms it: "body surface area (BSA) significantly
  # influences the apparent volume of distribution of R(-)-flurbiprofen, but
  # not S(+)-flurbiprofen".
  covariatesDataExcluded <- list(
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Median 2.6 m^2 (IQR 2.5-2.7, range 2.0-3.2; Yao 2025 Table 1). Entered the S(+) model on the volume of distribution during forward inclusion but was removed in backward elimination and does not appear in Table 2. Retained in the companion R(-) model on Vc."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Median 27.1 kg/m^2 (IQR 25.1-29.4, range 18.4-40.1; Yao 2025 Table 1). Screened as a candidate covariate; not retained."
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
      notes       = "CC/CT/TT = 28/32/7 (Yao 2025 Table 1). Screened as a candidate covariate on CL and Vd; not retained in the final S(+) model."
    ),
    SNP_PXR_RS1523127 = list(
      description = "PXR (NR1I2) rs1523127 genotype",
      units       = "(binary)",
      type        = "categorical",
      notes       = "AA/CA/CC = 3/60/4 (Yao 2025 Table 1); deviates from Hardy-Weinberg equilibrium in this cohort. Screened; not retained."
    ),
    SNP_POR_RS1057868 = list(
      description = "POR rs1057868 (POR*28) genotype",
      units       = "(binary)",
      type        = "categorical",
      notes       = "AA/GA/GG = 17/41/9 (Yao 2025 Table 1); deviates from Hardy-Weinberg equilibrium in this cohort. Retained on CL in the companion R(-) model but not in the S(+) model."
    ),
    SNP_POR_RS1135612 = list(
      description = "POR rs1135612 genotype",
      units       = "(binary)",
      type        = "categorical",
      notes       = "CC/CT/TT = 7/33/27 (Yao 2025 Table 1). Screened; not retained."
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
    # Structural PK parameters -- Yao 2025 Table 2 ("Pharmacokinetic
    # Parameter Estimates and Bootstrap Results of Final PPK Model for
    # S(+)-Flurbiprofen"), Estimate column. Plasma is the central
    # compartment and CSF the peripheral compartment (Results,
    # Pharmacokinetic Study: "Plasma and CSF were conceptualized as the
    # central and peripheral compartments, respectively"), so the
    # paper's Vp is the CSF volume (encoded lvcsf, following the
    # Kumpulainen_2010_flurbiprofen precedent) and the paper's Q is the
    # plasma-CSF intercompartmental clearance.
    #
    # These are APPARENT values relative to the nominal 100 mg
    # flurbiprofen axetil dose -- Table 2's own abbreviation footnote
    # reads "CL, apparent clearance; VC, apparent volume of
    # distribution in the central compartment" even though dosing was
    # intravenous. The paper never states the amount entered in the
    # dataset; the vignette Errata documents the 100 mg basis and the
    # cross-check against Zhang 2018, which used the same product and
    # study design at the same centre.
    # ================================================================
    lvc   <- log(25.6); label("Apparent central (plasma) volume of distribution Vc (L)")      # Table 2: Vc = 25.6 L (bootstrap median 25.3, 95% CI 21.5-28.6)
    lcl   <- log(16.7); label("Apparent plasma clearance CL (L/h)")                           # Table 2: CL = 16.7 L/h (bootstrap median 13.6, 95% CI 7.2-23.6); Discussion reports the same estimate to more digits as 16.67 L/h
    lvcsf <- log(32.6); label("Apparent peripheral (CSF) volume of distribution Vp (L)")      # Table 2: Vp = 32.6 L (bootstrap median 37.0, 95% CI 21.5-79.1); abstract: "the Vd of ... S(+)-flurbiprofen (CSF 32.6 L VS plasma 25.6 L)"
    lq    <- log(0.39); label("Plasma-CSF intercompartmental clearance Q (L/h)")              # Table 2: Q = 0.39 L/h (bootstrap median 0.44, 95% CI 0.25-0.93); abstract restates this as the "CSF CL" of S(+)-flurbiprofen (0.39 L/h)

    # ================================================================
    # Covariate effects on CL -- ABCB1 rs1045642 genotype.
    # Methods, Population Covariate Analysis: "Categorical covariates,
    # including sex and genetic polymorphisms, were modeled using
    # indicator variables within an exponential function (Equation II),
    # with the reference category coded as 0 and other categories
    # assigned values of 1 or higher." Table 2 lists one coefficient
    # for GA and one for GG, so AA (n = 12) is the reference group and
    # the two indicators enter a single exponential:
    #   CL = tvCL * exp(e_GA * I(GA) + e_GG * I(GG)) * exp(eta)
    # Sign check against the Discussion: "the GA genotype of the ABCB1
    # gene inhibits the clearance of S-flurbiprofen, while the GG
    # genotype increases its clearance rate" -- exp(-1.52) = 0.219
    # (78% lower) and exp(0.19) = 1.21 (21% higher). Consistent.
    # ================================================================
    e_snp_abcb1_rs1045642_ga_cl <- -1.52; label("Exponential coefficient on CL for ABCB1 rs1045642 GA vs AA (unitless)") # Table 2: ABCB1 (rs1045642) on CL, GA = -1.52 (bootstrap median -0.96, 95% CI -3.85 to -0.10)
    e_snp_abcb1_rs1045642_gg_cl <-  0.19; label("Exponential coefficient on CL for ABCB1 rs1045642 GG vs AA (unitless)") # Table 2: ABCB1 (rs1045642) on CL, GG = 0.19 (bootstrap median 0.29, 95% CI -0.16 to 0.93)

    # ================================================================
    # Inter-individual variability -- Yao 2025 Table 2, "Interindividual
    # variability (%)" block. The table's own abbreviation footnote
    # defines the reported quantity as "omega, square root of
    # interindividual variance for parameters", so each tabulated value
    # is the SD of eta on the log scale and the internal nlmixr2
    # variance is omega^2. (The "(%)" in the block header is
    # inconsistent with the footnote; taking the values as CV fractions
    # instead would give omega^2 = log(1 + CV^2), which differs by less
    # than 3% for every value here, so the distinction is immaterial.)
    # The exponential IIV model is stated in Methods, Structural Model:
    # "This study used the exponential model to describe the
    # inter-individual variability in pharmacokinetic (PK) parameters."
    # No correlations between etas were reported.
    # ================================================================
    etalvc   ~ 0.0169  # Table 2: omega_Vc = 0.13 -> variance 0.13^2 = 0.0169
    etalcl   ~ 0.0256  # Table 2: omega_CL = 0.16 -> variance 0.16^2 = 0.0256
    etalvcsf ~ 0.0625  # Table 2: omega_Vp = 0.25 -> variance 0.25^2 = 0.0625
    etalq    ~ 0.0225  # Table 2: omega_Q  = 0.15 -> variance 0.15^2 = 0.0225

    # ================================================================
    # Residual error -- Yao 2025 Table 2, "Random residual variability"
    # block: one additive error on the plasma observations and one
    # multiplicative (proportional) error on the CSF observations. This
    # is what the S(+) Results section calls "a combined additive and
    # proportional error model": the combination is across matrices,
    # not within one matrix.
    #
    # Both magnitudes are implausibly small relative to the observed
    # concentrations (plasma is a few mg/L; an additive SD of 0.003 mg/L
    # = 3 ng/mL sits far below the 0.1 ug/mL plasma LLOQ, and a 0.1%
    # proportional CSF error is far tighter than any bioanalytical
    # assay). They are transcribed exactly as published. The same
    # pattern appears in Zhang 2018 (sigma = 0.0023 mg/L additive) from
    # the same group and centre, so it reflects this group's Phoenix
    # reporting convention rather than a transcription error here. The
    # vignette Errata records the implication: simulated profiles are
    # effectively noise-free apart from the IIV.
    # ================================================================
    addSd        <- 0.003; label("Additive residual error SD on plasma concentration Cc (mg/L)")   # Table 2: "Plasma, additive error, sigma" = 0.003
    propSd_Ccsf  <- 0.001; label("Proportional residual error SD on CSF concentration Ccsf (fraction)") # Table 2: "CSF, multiplicative error, sigma" = 0.001
  })

  model({
    # ---- Individual PK parameters ---------------------------------
    # ABCB1 rs1045642 genotype acts on CL through two mutually
    # exclusive indicators, both 0 for the AA reference group
    # (Yao 2025 Methods Equation II; Table 2).
    cl   <- exp(lcl + etalcl) *
            exp(e_snp_abcb1_rs1045642_ga_cl * SNP_ABCB1_RS1045642_GA +
                e_snp_abcb1_rs1045642_gg_cl * SNP_ABCB1_RS1045642_GG)
    vc   <- exp(lvc   + etalvc)
    vcsf <- exp(lvcsf + etalvcsf)
    q    <- exp(lq    + etalq)

    # ---- Concentrations -------------------------------------------
    Cc   <- central / vc    # S(+)-flurbiprofen in plasma (mg/L)
    Ccsf <- csf / vcsf      # S(+)-flurbiprofen in CSF (mg/L)

    # ---- ODE system ------------------------------------------------
    # Standard two-compartment IV disposition with plasma as the
    # central compartment and CSF as the peripheral compartment
    # (Results: "Plasma and CSF were conceptualized as the central and
    # peripheral compartments, respectively (Figure S1)"). Elimination
    # is first-order from plasma only; the paper reports no separate
    # CSF elimination pathway -- the abstract's "CSF CL" is the
    # intercompartmental clearance Q of Table 2.
    d/dt(central) <- -cl * Cc - q * (Cc - Ccsf)
    d/dt(csf)     <-            q * (Cc - Ccsf)

    # ---- Residual error --------------------------------------------
    Cc   ~ add(addSd)
    Ccsf ~ prop(propSd_Ccsf)
  })
}
