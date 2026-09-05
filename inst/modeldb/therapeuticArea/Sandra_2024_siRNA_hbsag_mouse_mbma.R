Sandra_2024_siRNA_hbsag_mouse_mbma <- function() {
  description <- paste(
    "MBMA. Preclinical (mouse).",
    "Kinetic-pharmacodynamic (K-PD) model coupled to an inhibitory indirect",
    "response model (IRM) for the baseline- and placebo-corrected serum",
    "hepatitis B surface antigen (HBsAg) time course after a single dose of an",
    "anti-HBV short interfering RNA (siRNA) in HBV-infected mice. Fitted by",
    "model-based meta-analysis to 237 study-arm-mean HBsAg observations",
    "digitised from 25 active treatment arms in 9 placebo- or",
    "vehicle-controlled non-clinical studies, covering 13 siRNAs in 3 delivery",
    "classes: N-acetylgalactosamine-conjugated (GalNAc), lipid-nanoparticle",
    "formulated (LNP), and cholesterol-conjugated (chol, all co-administered",
    "with a GalNAc-conjugated melittin-like endolytic peptide). The dose enters",
    "a virtual biophase compartment depot_kpd as an absolute amount in ng and",
    "decays at the first-order rate kel; the amount in depot_kpd inhibits",
    "synthesis of the response through a sigmoid Imax term with drug-specific",
    "id50 and hill. kel is siRNA-CLASS specific (GalNAc 0.03291, LNP 0.1368,",
    "chol 0.2403 /day, i.e. biophase half-lives of 21.06, 5.07 and 2.89 days);",
    "id50 is siRNA-SPECIFIC and spans 41.46 ng (chol NAG-MLP chol-siHBV-74+77)",
    "to 7451 ng (GalNAc VIR-2218); hill was estimated for the 3 siRNAs with",
    "more than one dose level and FIXED to 1 for the other 10. The single",
    "system parameter kdeg (0.7196 /day) is shared across all siRNAs and",
    "carries the only random effect, a BETWEEN-STUDY variance (CV 31.4%).",
    "The output Cc is the double-corrected HBsAg ratio: it starts at exactly 1",
    "(ksyn = kdeg, so the response is at steady state before treatment) and",
    "falls toward 0 as the siRNA suppresses HBsAg synthesis; it is NOT a drug",
    "concentration. Suitable simulation scope is the study-arm-mean HBsAg",
    "ratio time course, NOT individual-mouse HBsAg and NOT plasma or liver",
    "siRNA pharmacokinetics -- the paper notes that a K-PD biophase amount",
    "cannot be extrapolated back to plasma PK because siRNA plasma and liver",
    "kinetics are completely temporally disconnected."
  )

  reference <- paste(
    "Sandra L, T'jollyn H, Vermeulen A, Ackaert O, Perez-Ruixo JJ.",
    "Model-based meta-analysis to quantify the effects of short interfering",
    "RNA therapeutics on hepatitis B surface antigen turnover in hepatitis",
    "B-infected mice.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(5):729-742.",
    "doi:10.1002/psp4.13129"
  )

  vignette <- "Sandra_2024_siRNA_hbsag_mouse_mbma"

  units <- list(
    time          = "day",
    dosing        = "ng",
    concentration = "fraction/fraction"
  )

  # The dose is a bolus into the virtual biophase compartment (Eq 5:
  # K(t = 0) = Dose), which buildModelDb() cannot auto-detect because it is
  # neither `depot` nor `central`.
  dosing <- c("depot_kpd")

  covariateData <- list(
    # ------------------------------------------------------------------
    # siRNA identity. Exactly one of the 13 TRT_* indicators is 1 on any
    # given study arm (Sandra 2024 eligibility criterion 2: siRNA given as
    # monotherapy). The indicator selects the drug-specific id50 and hill,
    # and -- because each siRNA belongs to exactly one delivery class --
    # also determines the class-specific kel. The class is DERIVED inside
    # model() from these indicators rather than carried as three further
    # columns, so it is impossible to specify a drug/class pair that
    # Sandra 2024 Table 1 does not contain.
    #
    # These are MBMA study-arm-level covariates (properties of a
    # published trial arm, not of an individual mouse). They are
    # registered as a single grouped entry in
    # inst/references/covariate-columns.md, following the FLARE /
    # NAPROXEN / TRAMADOL MBMA-arm precedent, and are well-formed members
    # of the auto-approved TRT_<drug> family alongside TRT_EPHEDRINE and
    # TRT_CAFEDRINE_THEODRENALINE.
    # ------------------------------------------------------------------
    TRT_OLX703A = list(
      description        = "Study-arm indicator that the administered siRNA is OLX703A (GalNAc-conjugated).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (a different siRNA was administered on this arm)",
      notes              = "MBMA study-arm-level treatment indicator. GalNAc class, so it selects kel = exp(lkel_galnac). Sandra 2024 Table 1 (Olix business report 2022, ref 41; 9 mg/kg SC, AAV-HBV mouse) and Table 3 (ID50 = 2581 ng, RSE 37.09%). Exactly one TRT_* indicator is 1 per arm.",
      source_name        = "OLX703A (Sandra 2024 Tables 1 and 3)"
    ),
    TRT_AB729 = list(
      description        = "Study-arm indicator that the administered siRNA is AB-729 (imdusiran; GalNAc-conjugated).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (a different siRNA was administered on this arm)",
      notes              = "MBMA study-arm-level treatment indicator. GalNAc class. One of the three siRNAs for which the Hill exponent was ESTIMATED rather than fixed to 1 (gamma = 2.547), because more than one dose level was available. Sandra 2024 Table 1 (Thi et al. 2022 ref 37, 1 and 3 mg/kg SC; Amy C.H. Lee 2018 ref 40, 1/3/9 mg/kg SC and 3 mg/kg SC) and Table 3 (ID50 = 4435 ng, RSE 30.50%).",
      source_name        = "AB-729 (Sandra 2024 Tables 1 and 3)"
    ),
    TRT_ALG125755 = list(
      description        = "Study-arm indicator that the administered siRNA is ALG-125755 (GalNAc-conjugated).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (a different siRNA was administered on this arm)",
      notes              = "MBMA study-arm-level treatment indicator. GalNAc class. Sandra 2024 Table 1 (Hong et al. 2022, ref 38; 5 mg/kg SC, AAV-HBV mouse) and Table 3 (ID50 = 6287 ng, RSE 16.52%).",
      source_name        = "ALG-125755 (Sandra 2024 Tables 1 and 3)"
    ),
    TRT_ALG125819 = list(
      description        = "Study-arm indicator that the administered siRNA is ALG-125819 (GalNAc-conjugated).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (a different siRNA was administered on this arm)",
      notes              = "MBMA study-arm-level treatment indicator. GalNAc class. Sandra 2024 Table 1 (Hong et al. 2022, ref 38; 5 mg/kg SC, AAV-HBV mouse) and Table 3 (ID50 = 7426 ng, RSE 15.60%).",
      source_name        = "ALG-125819 (Sandra 2024 Tables 1 and 3)"
    ),
    TRT_VIR2218 = list(
      description        = "Study-arm indicator that the administered siRNA is VIR-2218 (elebsiran; GalNAc-conjugated).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (a different siRNA was administered on this arm)",
      notes              = "MBMA study-arm-level treatment indicator. GalNAc class; the least potent siRNA in the analysis. Sandra 2024 Table 1 (Noack et al. 2022, ref 39; 3 mg/kg SC, AAV-HBV C57BL/6 mouse) and Table 3 (ID50 = 7451 ng, RSE 28.69%). Sandra 2024 Discussion notes that the derived biophase half-life of 21.06 days agrees with the 21.4 days Boianelli et al. (ref 43) estimated for the structurally very similar ALN-HBV02.",
      source_name        = "VIR-2218 (Sandra 2024 Tables 1 and 3)"
    ),
    TRT_ARB1740 = list(
      description        = "Study-arm indicator that the administered siRNA is ARB-1740 (lipid-nanoparticle formulated).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (a different siRNA was administered on this arm)",
      notes              = "MBMA study-arm-level treatment indicator. LNP class. One of the three siRNAs for which the Hill exponent was ESTIMATED rather than fixed to 1 (gamma = 2.03), because four dose levels were available. Sandra 2024 Table 1 (Thi et al. 2019, ref 36; 0.03/0.1/0.3/1 mg/kg IV, AAV-HBV C57BL/6 mouse) and Table 3 (ID50 = 72.01 ng, RSE 41.37%). Sandra 2024 Results/Discussion identify this study as the single driver of the retained between-study variability on kdeg (kdeg about twofold larger than in the other studies).",
      source_name        = "ARB-1740 (Sandra 2024 Tables 1 and 3)"
    ),
    TRT_ARB1467 = list(
      description        = "Study-arm indicator that the administered siRNA is ARB-1467 (lipid-nanoparticle formulated).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (a different siRNA was administered on this arm)",
      notes              = "MBMA study-arm-level treatment indicator. LNP class. Sandra 2024 Table 1 (Amy C.H. Lee 2018, ref 40; 0.3 and 3 mg/kg IV, AAV-HBV mouse) and Table 3 (ID50 = 124.0 ng, RSE 53.57%).",
      source_name        = "ARB-1467 (Sandra 2024 Tables 1 and 3)"
    ),
    TRT_ARC520 = list(
      description        = "Study-arm indicator that the administered siRNA is ARC-520 (cholesterol-conjugated, co-administered with the GalNAc-conjugated melittin-like peptide NAG-MLP).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (a different siRNA was administered on this arm)",
      notes              = "MBMA study-arm-level treatment indicator. chol class. Sandra 2024 Table 1 (Trubetskoy et al. 2017, ref 35; 8 mg/kg IV, HDI NOD-SCID mouse) and Table 3 (ID50 = 44.01 ng, RSE 34.59%). Sandra 2024 Discussion records that clinical development of ARC-520 was discontinued for toxicity attributed to the NAG-MLP excipient ARC-EX1, which this model does not represent.",
      source_name        = "ARC-520 (Sandra 2024 Tables 1 and 3)"
    ),
    TRT_SIHBV74 = list(
      description        = "Study-arm indicator that the administered siRNA is NAG-MLP chol-siHBV-74 (cholesterol-conjugated, co-administered with the GalNAc-conjugated melittin-like peptide NAG-MLP).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (a different siRNA was administered on this arm)",
      notes              = "MBMA study-arm-level treatment indicator. chol class. Sandra 2024 Table 1 (Wooddell et al. 2013 ref 34, 0.25 and 6 mg/kg IV; Sebastyen et al. 2015 ref 33, 3 and 6 mg/kg IV; HDI NOD-SCID mouse) and Table 3 (ID50 = 48.97 ng, RSE 26.37%).",
      source_name        = "NAG-MLP chol-siHBV-74 (Sandra 2024 Tables 1 and 3)"
    ),
    TRT_SIHBV75 = list(
      description        = "Study-arm indicator that the administered siRNA is NAG-MLP chol-siHBV-75 (cholesterol-conjugated, co-administered with the GalNAc-conjugated melittin-like peptide NAG-MLP).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (a different siRNA was administered on this arm)",
      notes              = "MBMA study-arm-level treatment indicator. chol class. Sandra 2024 Table 1 (Wooddell et al. 2013, ref 34; 0.25 and 6 mg/kg IV, HDI NOD-SCID mouse) and Table 3 (ID50 = 166.9 ng, RSE 37.15%).",
      source_name        = "NAG-MLP chol-siHBV-75 (Sandra 2024 Tables 1 and 3)"
    ),
    TRT_SIHBV76 = list(
      description        = "Study-arm indicator that the administered siRNA is NAG-MLP chol-siHBV-76 (cholesterol-conjugated, co-administered with the GalNAc-conjugated melittin-like peptide NAG-MLP).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (a different siRNA was administered on this arm)",
      notes              = "MBMA study-arm-level treatment indicator. chol class; the least potent of the cholesterol-conjugated siRNAs. Sandra 2024 Table 1 (Wooddell et al. 2013, ref 34; 0.25 and 6 mg/kg IV, HDI NOD-SCID mouse) and Table 3 (ID50 = 275.4 ng, RSE 33.97%).",
      source_name        = "NAG-MLP chol-siHBV-76 (Sandra 2024 Tables 1 and 3)"
    ),
    TRT_SIHBV77 = list(
      description        = "Study-arm indicator that the administered siRNA is NAG-MLP chol-siHBV-77 (cholesterol-conjugated, co-administered with the GalNAc-conjugated melittin-like peptide NAG-MLP).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (a different siRNA was administered on this arm)",
      notes              = "MBMA study-arm-level treatment indicator. chol class. One of the three siRNAs for which the Hill exponent was ESTIMATED rather than fixed to 1 (gamma = 1.313), because more than one dose level was available. Sandra 2024 Table 1 (Wooddell et al. 2013 ref 34, 0.25 and 6 mg/kg IV; Sebastyen et al. 2015 ref 33, 3 and 6 mg/kg IV; HDI NOD-SCID mouse) and Table 3 (ID50 = 274.0 ng, RSE 32.58%).",
      source_name        = "NAG-MLP chol-siHBV-77 (Sandra 2024 Tables 1 and 3)"
    ),
    TRT_SIHBV74_77 = list(
      description        = "Study-arm indicator that the administered siRNA is the fixed combination NAG-MLP chol-siHBV-74 + 77 (cholesterol-conjugated, co-administered with the GalNAc-conjugated melittin-like peptide NAG-MLP).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (a different siRNA was administered on this arm)",
      notes              = "MBMA study-arm-level treatment indicator. chol class; the most potent entity in the analysis. Sandra 2024 treats the 74 + 77 co-formulation as a single entity with its own ID50 rather than as a combination of two separately-parameterised siRNAs, so it gets its own indicator and must NOT be encoded by setting both TRT_SIHBV74 and TRT_SIHBV77 to 1. Sandra 2024 Table 1 (Sebastyen et al. 2015, ref 33; 3 and 6 mg/kg IV, HDI NOD-SCID mouse) and Table 3 (ID50 = 41.46 ng, RSE 34.05%).",
      source_name        = "NAG-MLP chol-siHBV-74 + 77 (Sandra 2024 Tables 1 and 3)"
    )
  )

  compartmentData <- list(
    depot_kpd = list(
      analyte  = "short interfering RNA (siRNA)",
      units    = "ng",
      specimen = "administration site",
      verified = TRUE
    ),
    effect = list(
      analyte  = "hepatitis B surface antigen (HBsAg), baseline- and placebo-corrected ratio",
      units    = "fraction",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species       = "mouse (HBV-infected: hydrodynamic-injection NOD-SCID and AAV-HBV C57BL/6)",
    n_subjects    = 147L,
    n_studies     = 9L,
    disease_state = "HBV-infected mice: hydrodynamic injection (HDI) of pHBV1.3 plasmid (10-13.5 ug, NOD-SCID, 21 days post-infection) or adeno-associated viral vector AAV-HBV (AAV2/8 or AAV8, 1e11 vg, mostly C57BL/6, 28 days post-infection)",
    age_range     = "6-8 weeks where reported (Sandra 2024 Table 1; not available for the AAV-HBV studies)",
    weight_median = "20 g where reported; a mean total bodyweight of 20 g was ASSUMED where the study did not report one, in order to convert mg/kg doses to the absolute ng doses this model consumes (Sandra 2024 'Data extraction and curation')",
    sex           = "female (HDI studies where reported), male (Thi et al. 2019), not reported for most AAV-HBV studies",
    dose_range    = "0.03-9 mg/kg as published (Sandra 2024 Table 1); IV for the chol- and LNP-siRNAs, SC for the GalNAc-siRNAs. Converted to absolute ng amounts for this model, which spans roughly 600 ng to 180,000 ng at a 20 g bodyweight.",
    notes         = paste(
      "MBMA at the study-arm level: each modelled observation is the",
      "baseline- and placebo-corrected MEAN HBsAg of one published treatment",
      "arm at one timepoint (237 such observations from 25 active treatment",
      "arms). n_subjects = 147 is the sum of the n_study column of Sandra",
      "2024 Table 1 over its ten rows; Table 1 has ten rows but nine distinct",
      "studies, because Amy C.H. Lee 2018 (ref 40) contributed two separate",
      "experiments. The model predicts study-arm-mean responses and is NOT",
      "suitable for individual-mouse simulation. Route of administration is",
      "not a model covariate: Sandra 2024 Discussion states that SC",
      "absorption of GalNAc-siRNAs is complete and that plasma and liver",
      "kinetics are completely temporally disconnected, so no adaptation was",
      "made for IV vs SC. The type of murine HBV model (HDI vs AAV-HBV) was",
      "tested on kdeg and did not significantly improve the fit or reduce the",
      "between-study variability, so it is not a covariate either."
    )
  )

  ini({
    # ================================================================
    # System parameter (Sandra 2024 Eq 6, Table 3 'Fixed effects').
    # kdeg is the ONLY system-specific parameter and is shared across all
    # 13 siRNAs and all 3 classes. Because the dependent variable is both
    # baseline- and placebo-corrected, Sandra 2024 assumes the response is
    # at steady state before treatment, so ksyn / kdeg = 1 and therefore
    # ksyn = kdeg; ksyn is NOT a separate estimated parameter.
    # ================================================================
    lkdeg <- log(0.7196)
    label("Log HBsAg degradation rate kdeg (1/day)")  # Sandra 2024 Table 3: kdeg = 0.7196 /day (RSE 12.50%)

    # ================================================================
    # Class-specific biophase elimination rate (Sandra 2024 Eq 5, Table 3
    # 'KDE'). The paper's covariate analysis (Table 2) retained siRNA
    # CLASS on KDE; adding it dropped the objective function from 178.15
    # to 155.19 and reduced the between-study variability on KDE to a
    # negligible value that was then fixed to 0.
    # ================================================================
    lkel_galnac <- log(0.03291)
    label("Log biophase elimination rate kel for GalNAc-conjugated siRNAs (1/day; t1/2 = 21.06 days)")  # Sandra 2024 Table 3: KDE_GalNAc = 0.03291 /day (RSE 8.80%); Results: ln(2)/0.03291 = 21.06 days

    lkel_lnp <- log(0.1368)
    label("Log biophase elimination rate kel for lipid-nanoparticle siRNAs (1/day; t1/2 = 5.07 days)")  # Sandra 2024 Table 3: KDE_LNP = 0.1368 /day (RSE 9.53%); Results: ln(2)/0.1368 = 5.07 days

    lkel_chol <- log(0.2403)
    label("Log biophase elimination rate kel for cholesterol-conjugated siRNAs (1/day; t1/2 = 2.89 days)")  # Sandra 2024 Table 3: KDE_chol = 0.2403 /day (RSE 4.86%); Results: ln(2)/0.2403 = 2.89 days

    # ================================================================
    # siRNA-specific potency (Sandra 2024 Eq 7, Table 3 'ID50').
    # id50 is the AMOUNT in the biophase compartment (ng, not a
    # concentration) giving half-maximal inhibition of ksyn. Adding a
    # siRNA-specific ID50 on top of the class-specific KDE dropped the
    # objective function from 139.02 to 53.76 and reduced the
    # between-study variability on ID50 to a negligible value fixed to 0
    # (Sandra 2024 Table 2).
    # Listed below in the paper's Table 3 order: GalNAc, then LNP, then
    # chol, each block sorted by increasing ID50.
    # ================================================================
    lid50_olx703a <- log(2581)
    label("Log ID50 for OLX703A, GalNAc (ng in biophase)")  # Sandra 2024 Table 3: ID50 OLX703A = 2581 ng (RSE 37.09%)

    lid50_ab729 <- log(4435)
    label("Log ID50 for AB-729 / imdusiran, GalNAc (ng in biophase)")  # Sandra 2024 Table 3: ID50 AB-729 = 4435 ng (RSE 30.50%)

    lid50_alg125755 <- log(6287)
    label("Log ID50 for ALG-125755, GalNAc (ng in biophase)")  # Sandra 2024 Table 3: ID50 ALG-125755 = 6287 ng (RSE 16.52%)

    lid50_alg125819 <- log(7426)
    label("Log ID50 for ALG-125819, GalNAc (ng in biophase)")  # Sandra 2024 Table 3: ID50 ALG-125819 = 7426 ng (RSE 15.60%)

    lid50_vir2218 <- log(7451)
    label("Log ID50 for VIR-2218 / elebsiran, GalNAc (ng in biophase)")  # Sandra 2024 Table 3: ID50 VIR-2218 = 7451 ng (RSE 28.69%)

    lid50_arb1740 <- log(72.01)
    label("Log ID50 for ARB-1740, LNP (ng in biophase)")  # Sandra 2024 Table 3: ID50 ARB-1740 = 72.01 ng (RSE 41.37%)

    lid50_arb1467 <- log(124.0)
    label("Log ID50 for ARB-1467, LNP (ng in biophase)")  # Sandra 2024 Table 3: ID50 ARB-1467 = 124.0 ng (RSE 53.57%)

    lid50_sihbv74_77 <- log(41.46)
    label("Log ID50 for NAG-MLP chol-siHBV-74 + 77, chol (ng in biophase)")  # Sandra 2024 Table 3: ID50 NAG-MLP chol-siHBV-74 + 77 = 41.46 ng (RSE 34.05%)

    lid50_arc520 <- log(44.01)
    label("Log ID50 for ARC-520, chol (ng in biophase)")  # Sandra 2024 Table 3: ID50 ARC-520 = 44.01 ng (RSE 34.59%)

    lid50_sihbv74 <- log(48.97)
    label("Log ID50 for NAG-MLP chol-siHBV-74, chol (ng in biophase)")  # Sandra 2024 Table 3: ID50 NAG-MLP chol-siHBV-74 = 48.97 ng (RSE 26.37%)

    lid50_sihbv75 <- log(166.9)
    label("Log ID50 for NAG-MLP chol-siHBV-75, chol (ng in biophase)")  # Sandra 2024 Table 3: ID50 NAG-MLP chol-siHBV-75 = 166.9 ng (RSE 37.15%)

    lid50_sihbv77 <- log(274.0)
    label("Log ID50 for NAG-MLP chol-siHBV-77, chol (ng in biophase)")  # Sandra 2024 Table 3: ID50 NAG-MLP chol-siHBV-77 = 274.0 ng (RSE 32.58%)

    lid50_sihbv76 <- log(275.4)
    label("Log ID50 for NAG-MLP chol-siHBV-76, chol (ng in biophase)")  # Sandra 2024 Table 3: ID50 NAG-MLP chol-siHBV-76 = 275.4 ng (RSE 33.97%)

    # ================================================================
    # Hill / sigmoidicity exponent (Sandra 2024 Eq 7, Table 3 'gamma').
    # Table 3 heads this block "gamma (1 FIX all except)": gamma was
    # FIXED to 1 for every siRNA and only estimated for the three with
    # more than one dose level in the dataset, where estimating it
    # significantly improved the fit (Table 2 last row, OFV 53.76 ->
    # 4.85). lhill_other therefore carries the fixed value shared by the
    # remaining ten siRNAs.
    # ================================================================
    lhill_ab729 <- log(2.547)
    label("Log Hill exponent gamma for AB-729 (unitless)")  # Sandra 2024 Table 3: gamma AB-729 = 2.547 (RSE 17.53%)

    lhill_arb1740 <- log(2.03)
    label("Log Hill exponent gamma for ARB-1740 (unitless)")  # Sandra 2024 Table 3: gamma ARB-1740 = 2.03 (RSE 13.50%)

    lhill_sihbv77 <- log(1.313)
    label("Log Hill exponent gamma for NAG-MLP chol-siHBV-77 (unitless)")  # Sandra 2024 Table 3: gamma NAG-MLP chol-siHBV-77 = 1.313 (RSE 10.31%)

    lhill_other <- fixed(log(1))
    label("Log Hill exponent gamma for the ten remaining siRNAs (unitless)")  # Sandra 2024 Table 3 heading: "gamma (1 FIX all except)" -- gamma fixed to 1 for every siRNA other than AB-729, ARB-1740 and NAG-MLP chol-siHBV-77

    # ================================================================
    # BETWEEN-STUDY variability (Sandra 2024 Eq 8, Table 3 'Random
    # effects'). This is a STUDY-level random effect, NOT individual
    # between-subject variability -- one draw per published study, not
    # per mouse (see references/pbpk-qsp-mbma.md).
    #
    # Table 3 reports "CV%" = 31.4 and Table 2 the unrounded 31.39. That
    # column is the omega STANDARD DEVIATION times 100, not the
    # log-normal CV: the Results text separately converts it with a cited
    # formula (ref 42) to "a CV% of 32.31%", and there would be nothing
    # to convert if the table already held the log-normal CV. Confirming
    # the direction: sqrt(exp(0.3139^2) - 1) = 0.3218, which rounds to
    # the ~32.3% quoted in the text. The ini() value is the VARIANCE.
    #
    # The between-study variabilities on KDE and on ID50 are deliberately
    # ABSENT: Sandra 2024 Table 2 shows both shrank to negligible values
    # once the class effect on KDE and the siRNA effect on ID50 were
    # added, and its note records "Negligible BSV estimate values were
    # fixed to 0". They are omitted rather than written as
    # `~ fixed(0)` because a zero-variance diagonal makes the omega
    # matrix singular and breaks rxSolve with "chol(): decomposition
    # failed".
    # ================================================================
    eta_study_lkdeg ~ 0.09853321  # Sandra 2024 Table 3 BSV kdeg = 31.4% (RSE 30.8%), unrounded 31.39% in Table 2; variance = 0.3139^2 = 0.09853321

    # ================================================================
    # Residual variability (Sandra 2024 Eq 9): an additive error on the
    # log scale, i.e. log-normal on the natural scale, whose SD is NOT
    # constant across observations. Eq 9 is
    #   log(d2HBsAg_ijt) = log(pred_ijt) + CV%_ijt / sqrt(N_ijt - 1) * eps
    # with eps ~ N(0, sigma^2). expSd below is that sigma. Because the
    # weight CV%_ijt / sqrt(N_ijt - 1) is built from the DIGITISED
    # per-arm summary statistics of the fitted dataset -- the observed
    # within-arm CV% at that timepoint (Eq 1, on the PERCENT scale) and
    # that arm's animal count -- it is a property of the published data,
    # not of the model, and is left to downstream code exactly as in
    # `Boucher_2018_naproxen_mbma.R` and `Vargo_2014_statins_ezetimibe_mbma.R`.
    #
    # CONSEQUENCE, and it is a large one: applying expSd on its own
    # UNDERSTATES the residual badly. Eq 1 puts CV% on the percent scale,
    # so sigma comes out of order 1/100; for a typical arm in this
    # dataset (CV% about 30, n_arm 5) the operative log-scale SD is
    # 0.0342 * 30 / sqrt(4) = 0.51, roughly 15 times the bare 0.0342.
    # To simulate a study arm of size n with observed within-arm CV% c,
    # scale this SD by c / sqrt(n - 1). The paper's own Figure 4-6
    # simulations are DETERMINISTIC, so they need no residual at all.
    # ================================================================
    expSd <- 0.0342
    label("Residual error multiplier sigma (Eq 9); operative log-scale SD per observation is expSd * (arm CV% / sqrt(arm n - 1))")  # Sandra 2024 Table 3: RV = 3.42% (RSE 4.68%)
  })

  model({
    # ================================================================
    # siRNA delivery class, derived from the treatment indicator.
    # Each of the 13 siRNAs belongs to exactly one class (Sandra 2024
    # Table 1), so the class is a function of the drug and is derived
    # here rather than supplied, which makes an inconsistent drug/class
    # pair unrepresentable.
    # ================================================================
    is_galnac <- TRT_OLX703A + TRT_AB729 + TRT_ALG125755 +
      TRT_ALG125819 + TRT_VIR2218
    is_lnp <- TRT_ARB1740 + TRT_ARB1467
    is_chol <- TRT_ARC520 + TRT_SIHBV74 + TRT_SIHBV75 +
      TRT_SIHBV76 + TRT_SIHBV77 + TRT_SIHBV74_77

    # Class-specific biophase elimination rate. Exactly one of the three
    # class indicators is 1, so this selects one of the three estimates.
    kel <- exp(
      lkel_galnac * is_galnac +
        lkel_lnp * is_lnp +
        lkel_chol * is_chol
    )

    # siRNA-specific potency. Exactly one TRT_* indicator is 1, so this
    # selects one of the thirteen estimates.
    id50 <- exp(
      lid50_olx703a * TRT_OLX703A +
        lid50_ab729 * TRT_AB729 +
        lid50_alg125755 * TRT_ALG125755 +
        lid50_alg125819 * TRT_ALG125819 +
        lid50_vir2218 * TRT_VIR2218 +
        lid50_arb1740 * TRT_ARB1740 +
        lid50_arb1467 * TRT_ARB1467 +
        lid50_sihbv74_77 * TRT_SIHBV74_77 +
        lid50_arc520 * TRT_ARC520 +
        lid50_sihbv74 * TRT_SIHBV74 +
        lid50_sihbv75 * TRT_SIHBV75 +
        lid50_sihbv77 * TRT_SIHBV77 +
        lid50_sihbv76 * TRT_SIHBV76
    )

    # Hill exponent: the three estimated values, otherwise the fixed 1.
    # The complement term is 1 for the ten siRNAs whose gamma was FIXED.
    hill <- exp(
      lhill_ab729 * TRT_AB729 +
        lhill_arb1740 * TRT_ARB1740 +
        lhill_sihbv77 * TRT_SIHBV77 +
        lhill_other * (1 - TRT_AB729 - TRT_ARB1740 - TRT_SIHBV77)
    )

    # System parameter with its between-STUDY random effect (Eq 8).
    kdeg <- exp(lkdeg + eta_study_lkdeg)
    # Eq 6: the response is at steady state before treatment because the
    # data are baseline- and placebo-corrected, so ksyn / kdeg = 1.
    ksyn <- kdeg

    # Eq 7: fractional inhibition of HBsAg synthesis by the amount of
    # siRNA currently in the biophase. E = 1/2 exactly when
    # depot_kpd == id50, for any hill.
    inh <- depot_kpd^hill / (id50^hill + depot_kpd^hill)

    # Eq 5: the dose is a bolus into the virtual biophase compartment,
    # which then decays first order. There is no measured siRNA
    # concentration anywhere in this model.
    d/dt(depot_kpd) <- -kel * depot_kpd

    # Eq 6: inhibitory indirect response on the synthesis of the
    # double-corrected HBsAg ratio.
    d/dt(effect) <- ksyn * (1 - inh) - kdeg * effect
    # Steady state before treatment: effect(0) = ksyn / kdeg = 1.
    effect(0) <- 1

    # Cc is the baseline- and placebo-corrected HBsAg ratio
    # (Sandra 2024 Eq 4, "double delta"), NOT a drug concentration.
    Cc <- effect

    # Eq 9: additive residual on the log scale (log-normal on the
    # natural scale). expSd is the Eq 9 sigma only; the per-observation
    # reweighting by the published arm's CV% and size is left to
    # downstream code -- see the ini() note above, which explains why
    # applying expSd unscaled understates the residual by roughly 15x.
    Cc ~ lnorm(expSd)
  })
}
