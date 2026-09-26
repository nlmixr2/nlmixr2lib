Sokolov_2019_antipcsk9_qsp <- function() {
  description <- paste(
    "QSP. Sokolov 2019 lipoprotein-homeostasis model benchmarking two",
    "anti-PCSK9 modalities, monoclonal antibodies and small interfering RNA,",
    "in healthy subjects and hypercholesterolemia patients on background",
    "statins. Seventeen ODEs and 51 parameters. A four-state endogenous core",
    "carries plasma PCSK9 (nmol), LDL-C, the hepatic LDL-receptor pool as a",
    "ratio to baseline, and Lp(a); LDL-C production is proportional to",
    "VLDL-C, LDL-C clearance is LDL-receptor mediated, and LDL-receptor",
    "turnover is driven by a power function of PCSK9 relative to baseline",
    "with a negative-feedback power term on LDL-C. Lp(a) carries parallel",
    "LDL-receptor-dependent and independent clearance. Three mAbs",
    "(alirocumab, evolocumab, RG-7652) each get a one-compartment",
    "first-order-absorption PK model with explicit 1:1 bimolecular PCSK9",
    "binding and clearance of the complex; two siRNAs (inclisiran, ALN-PCS)",
    "each get a one-compartment lumped liver model whose liver amount drives",
    "fractional Imax inhibition of PCSK9 synthesis. Triglycerides, HDL-C,",
    "total cholesterol, non-HDL-C and apoB are algebraic readouts. All five",
    "drugs share one parameter set except the PCSK9-on-LDL-receptor exponent",
    "n1, which the authors fitted separately by modality and which the",
    "TRT_ANTIPCSK9_SIRNA covariate selects. Deterministic: the source fitted",
    "trial-level aggregate data by nonlinear fixed effects and reports no",
    "IIV and no residual-error magnitude. Five printed equations are",
    "internally inconsistent and are corrected here; the baseline PCSK9",
    "amount is not tabulated and was back-solved from Table 1. See the",
    "vignette Errata.",
    sep = " "
  )
  reference <- paste(
    "Sokolov V, Helmlinger G, Nilsson C, Zhudenkov K, Skrtic S, Hamren B,",
    "Peskov K, Hurt-Camejo E, Jansson-Lofmark R (2019). Comparative",
    "quantitative systems pharmacology modeling of anti-PCSK9 therapeutic",
    "modalities in hypercholesterolemia. J Lipid Res 60(9):1610-1621.",
    "doi:10.1194/jlr.M092486. PMCID: PMC6718444.",
    "Model equations from main-text equations 1-25; parameter values from",
    "supplemental Table S2; the TG and apoB partition constants from the",
    "Methods 'Structure of the mathematical model' narrative.",
    sep = " "
  )
  vignette <- "Sokolov_2019_antipcsk9_qsp"

  # Every state is a paper-mechanistic lipoprotein / PCSK9 / per-drug species.
  # The endogenous core states (pcsk9, ldlr, lpa) and the five per-drug
  # chains do not map onto canonical PK compartment roles: the model carries
  # five drugs simultaneously, so `depot` / `central` / `complex` each need a
  # drug qualifier, and none of the five qualifiers is a registered
  # metabolite suffix. `ldl` is canonical and is used unqualified. The siRNA
  # circulating state is named `liver_<drug>` because the source calls it
  # "the amount of drug in the liver" (main text following equation 24) and
  # figure 2A labels its axis "Plasma or liver PK".
  paper_specific_compartments <- c(
    "pcsk9",
    "ldlr",
    "lpa",
    "depot_aliro",
    "central_aliro",
    "complex_aliro",
    "depot_evolo",
    "central_evolo",
    "complex_evolo",
    "depot_rg",
    "central_rg",
    "complex_rg",
    "depot_inc",
    "liver_inc",
    "depot_aln",
    "liver_aln"
  )

  units <- list(
    time = "day",
    dosing = "mg",
    concentration = paste(
      "mg/dL for LDL-C, VLDL-C, HDL-C, Lp(a), total cholesterol, non-HDL-C",
      "and apoB; nmol/L for plasma PCSK9 and for the three mAbs; triglyceride",
      "is a unitless ratio to baseline",
      sep = " "
    )
  )

  covariateData <- list(
    TRT_ANTIPCSK9_SIRNA = list(
      description = "Anti-PCSK9 siRNA treatment-arm indicator",
      units = "unitless",
      type = "binary",
      reference_category = "0 (monoclonal-antibody arm: alirocumab, evolocumab or RG-7652)",
      notes = paste(
        "Selects which of the two fitted values of the PCSK9-on-LDL-receptor",
        "exponent n1 applies. Supplemental Table S2 lists n1 = 0.14 in each of",
        "the alirocumab, evolocumab and RG-7652 blocks and n1 = 0.26 in each of",
        "the inclisiran and ALN-PCS blocks; the Methods explain that",
        "'parameters of PCSK9 effects on LDL-R turnover were fitted separately:",
        "for mAbs (alirocumab and evolocumab) and for siRNA (inclisiran and",
        "ALN-PCS)'. Every other parameter in the model is shared across the two",
        "modalities, so this single indicator is the whole difference between",
        "the paper's 'two previously identified sets of parameters' (Results,",
        "'PCSK9-LDL-C relationship'). The model is not defined for an arm that",
        "co-administers an mAb and an siRNA; the source never simulates one.",
        sep = " "
      ),
      source_name = "not a source data column; the source encodes modality by which parameter set it runs"
    )
  )

  compartmentData <- list(
    pcsk9 = list(analyte = "PCSK9", units = "nmol", specimen = "plasma", verified = TRUE),
    ldl = list(analyte = "LDL cholesterol", units = "mg/dL", specimen = "plasma", verified = TRUE),
    ldlr = list(
      analyte = "hepatic LDL receptor",
      units = "fraction of baseline receptor number",
      specimen = "tissue",
      verified = TRUE
    ),
    lpa = list(analyte = "lipoprotein(a)", units = "mg/dL", specimen = "plasma", verified = TRUE),
    depot_aliro = list(
      analyte = "alirocumab",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central_aliro = list(analyte = "alirocumab", units = "nmol", specimen = "plasma", verified = TRUE),
    complex_aliro = list(
      analyte = "alirocumab-PCSK9 complex",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    depot_evolo = list(
      analyte = "evolocumab",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central_evolo = list(analyte = "evolocumab", units = "nmol", specimen = "plasma", verified = TRUE),
    complex_evolo = list(
      analyte = "evolocumab-PCSK9 complex",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    depot_rg = list(analyte = "RG-7652", units = "nmol", specimen = "administration site", verified = TRUE),
    central_rg = list(analyte = "RG-7652", units = "nmol", specimen = "plasma", verified = TRUE),
    complex_rg = list(
      analyte = "RG-7652-PCSK9 complex",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    depot_inc = list(analyte = "inclisiran", units = "mg", specimen = "administration site", verified = TRUE),
    liver_inc = list(analyte = "inclisiran", units = "mg", specimen = "tissue", verified = TRUE),
    depot_aln = list(analyte = "ALN-PCS", units = "mg", specimen = "administration site", verified = TRUE),
    liver_aln = list(analyte = "ALN-PCS", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 17L,
    disease_state = paste(
      "Healthy subjects and patients with familial or nonfamilial",
      "hypercholesterolemia. Four trials enrolled healthy subjects and 16",
      "primarily enrolled hypercholesterolemia patients; roughly 70% of the",
      "hypercholesterolemia subjects on the mAb trials also received statins",
      "as standard of care. For the siRNA trials only 13% of healthy subjects",
      "received concomitant statins, while all hypercholesterolemia subjects",
      "did.",
      sep = " "
    ),
    dose_range = paste(
      "Alirocumab 50-300 mg SC; evolocumab 7-420 mg SC; RG-7652 10-800 mg SC",
      "single or 40-150 mg on days 1 and 14; inclisiran 25-800 mg SC;",
      "ALN-PCS 0.015-0.4 mg/kg IV. Single and multiple dosing, treatment",
      "periods 0.5 to 24 weeks for the mAbs and up to 8 months for inclisiran.",
      sep = " "
    ),
    notes = paste(
      "Study-level aggregated (per dosing arm) mean data from 17 published",
      "clinical trials comprising 68 dosing arms, enumerated in supplemental",
      "Table S1; subject-level data were not available, so the source could",
      "not apply mixed-effects modeling and reports no between-subject",
      "variability. n_subjects is not reported as a pooled total. Per-arm",
      "baseline PCSK9, LDL-C, HDL-C and Lp(a) were set from each trial's own",
      "published baselines and are not tabulated; see the vignette Errata for",
      "the shipped defaults.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Supplemental Table S2 tabulates all 51 parameters. Time is in DAYS.
    #
    # FIXED vs ESTIMATED follows Table S2's "Estimation method and
    # References" column exactly: rows marked 'fitted' (30 of them, counting
    # the n1 row repeated in each of the five per-drug blocks) are estimated
    # and carry the tabulated 95% confidence interval in their comment; rows
    # marked 'taken from the literature', 'taken from the FDA report',
    # 'assumed equal to ...' or 'calculated' (21 rows) are wrapped in
    # fixed(). 30 + 21 = 51, which reconciles with the Methods statement
    # that 21 of 51 parameters came from the literature.
    #
    # Parameters are NOT log-transformed: the source fitted them on the
    # natural scale by nonlinear fixed effects and reports symmetric
    # confidence intervals from the Fisher information matrix.
    # =====================================================================

    # --- Endogenous lipoprotein core -----------------------------------
    kPCSK9deg <- 1.5; label("PCSK9 degradation rate constant (1/day)") # Table S2, fitted, 95% CI [1.4; 1.61]
    kLDLcdeg <- fixed(0.231); label("LDL-C clearance rate constant (1/day)") # Table S2, taken from the literature (ref 22)
    kLDLrturn <- 3.37; label("LDL-receptor turnover rate constant (1/day)") # Table S2, fitted, 95% CI [1.58; 7.21]
    kLpAdeg <- 0.09; label("LDL-receptor-independent Lp(a) degradation rate constant (1/day)") # Table S2, fitted, 95% CI [0.06; 0.13]
    kLpAdeg2 <- 0.04; label("LDL-receptor-dependent Lp(a) degradation rate constant (1/day)") # Table S2, fitted, 95% CI [0.03; 0.06]
    n2 <- 0.52; label("Power of LDL-C feedback on LDL-receptor degradation (unitless)") # Table S2, fitted, 95% CI [0.44; 0.6]

    # n1 takes one of two fitted values, selected by TRT_ANTIPCSK9_SIRNA.
    n1mab <- 0.14; label("Power of PCSK9 on LDL-receptor degradation, mAb arms (unitless)") # Table S2, fitted, 95% CI [0.12; 0.17]; identical in the alirocumab, evolocumab and RG-7652 blocks
    n1sirna <- 0.26; label("Power of PCSK9 on LDL-receptor degradation, siRNA arms (unitless)") # Table S2, fitted, 95% CI [0.25; 0.27]; identical in the inclisiran and ALN-PCS blocks

    # --- Plasma volume and lipid partition constants --------------------
    Vpl <- fixed(2.75); label("Plasma volume (L)") # Table S2, taken from the literature (ref 20)
    lamtg <- 0.34; label("Influence of triglyceride on HDL-C (unitless)") # Table S2, fitted, 95% CI [0.23; 0.47]
    lamapoB <- fixed(0.654); label("Non-HDL-C to apoB conversion coefficient (unitless)") # Table S2, 'calculated'; printed there as -0.654, sign corrected -- see vignette Errata
    lamtgVLDL <- fixed(0.78); label("Fraction of the plasma triglyceride pool carried by VLDL (unitless)") # Methods: 'primarily divided between VLDL (78%) and LDL (22%) particles'
    lamtgLDL <- fixed(0.22); label("Fraction of the plasma triglyceride pool carried by LDL (unitless)") # Methods: 'primarily divided between VLDL (78%) and LDL (22%) particles'
    lamapoBLDL <- fixed(0.9); label("Fraction of plasma apoB carried by LDL particles (unitless)") # Methods: 'LDL contributes approximately 90% to plasma apoB content in healthy subjects'
    lamapoBVLDL <- fixed(0.1); label("Fraction of plasma apoB carried by VLDL particles (unitless)") # Methods: complement of the 90% LDL apoB contribution

    # --- Baselines. Table S2 marks PCSK9, LDL-C, Lp(a) and HDL-C baselines
    #     'taken from each arm of each trial' and tabulates no pooled value;
    #     only the VLDL-C baseline carries a number. See vignette Errata.
    BaselinePCSK9 <- fixed(6.3832); label("Baseline plasma PCSK9 amount (nmol)") # back-solved from Table 1 -- not printed in any source; see vignette Errata
    BaselineVLDLc <- fixed(23.166); label("Baseline VLDL-C (mg/dL)") # Table S2, median across the evolocumab and inclisiran trials
    BaselineLDLc <- fixed(100); label("Baseline LDL-C (mg/dL)") # not printed in any source; rounded placeholder, see vignette Errata
    BaselineHDLc <- fixed(50); label("Baseline HDL-C (mg/dL)") # not printed in any source; rounded placeholder, see vignette Errata
    BaselineLpA <- fixed(30); label("Baseline Lp(a) (mg/dL)") # not printed in any source; rounded placeholder, see vignette Errata

    # --- Alirocumab ------------------------------------------------------
    MWaliro <- fixed(146000); label("Alirocumab molecular weight (g/mol)") # Table S2, taken from the FDA report
    kabsaliro <- 0.1; label("Alirocumab first-order absorption rate constant (1/day)") # Table S2, fitted, 95% CI [0.1; 0.11]
    CLaliro <- 0.52; label("Alirocumab clearance (L/day)") # Table S2, fitted, 95% CI [0.48; 0.56]
    Vdaliro <- 1.37; label("Alirocumab volume of distribution (L)") # Table S2, fitted, 95% CI [1.17; 1.6]
    konaliro <- 0.94; label("Alirocumab-PCSK9 association rate constant (L/nmol/day)") # Table S2, fitted, 95% CI [0.71; 1.24]
    Kdaliro <- 0.52; label("Alirocumab-PCSK9 dissociation constant (nmol/L)") # Table S2, fitted, 95% CI [0.38; 0.71]

    # --- Evolocumab ------------------------------------------------------
    MWevolo <- fixed(141800); label("Evolocumab molecular weight (g/mol)") # Table S2, taken from the FDA report
    kabsevolo <- 0.095; label("Evolocumab first-order absorption rate constant (1/day)") # Table S2, fitted, 95% CI [0.09; 0.101]
    CLevolo <- 0.454; label("Evolocumab clearance (L/day)") # Table S2, fitted, 95% CI [0.433; 0.475]
    Vdevolo <- 1.34; label("Evolocumab volume of distribution (L)") # Table S2, fitted, 95% CI [1.233; 1.447]
    konevolo <- fixed(0.94); label("Evolocumab-PCSK9 association rate constant (L/nmol/day)") # Table S2, assumed equal to konaliro
    Kdevolo <- fixed(0.016); label("Evolocumab-PCSK9 dissociation constant (nmol/L)") # Table S2, taken from the FDA report

    # --- RG-7652 ---------------------------------------------------------
    MWrg <- fixed(141800); label("RG-7652 molecular weight (g/mol)") # Table S2, assumed equal to MWaliro
    kabsrg <- 0.05; label("RG-7652 first-order absorption rate constant (1/day)") # Table S2, fitted, 95% CI [0.04; 0.05]
    CLrg <- 0.35; label("RG-7652 clearance (L/day)") # Table S2, fitted, 95% CI [0.33; 0.37]
    Vdrg <- 0.71; label("RG-7652 volume of distribution (L)") # Table S2, fitted, 95% CI [0.61; 0.82]
    konrg <- fixed(0.94); label("RG-7652-PCSK9 association rate constant (L/nmol/day)") # Table S2, assumed equal to konaliro
    Kdrg <- fixed(0.52); label("RG-7652-PCSK9 dissociation constant (nmol/L)") # Table S2, assumed equal to Kdaliro

    # --- Inclisiran ------------------------------------------------------
    kabsinc <- 0.04; label("Inclisiran first-order liver-uptake rate constant (1/day)") # Table S2, fitted, 95% CI [0.03; 0.05]
    kelinc <- 0.01; label("Inclisiran liver elimination rate constant (1/day)") # Table S2, fitted, 95% CI [0.01; 0.02]
    Imaxinc <- 0.77; label("Maximum fractional inhibition of PCSK9 synthesis by inclisiran (unitless)") # Table S2, fitted, 95% CI [0.75; 0.79]
    ID50inc <- 21.74; label("Liver inclisiran amount giving half-maximal PCSK9 synthesis inhibition (mg)") # Table S2, fitted, 95% CI [15.32; 30.85]

    # --- ALN-PCS ---------------------------------------------------------
    kabsaln <- 2.59; label("ALN-PCS first-order liver-uptake rate constant (1/day)") # Table S2, fitted, 95% CI [0.7; 9.55]
    kelaln <- 0.13; label("ALN-PCS liver elimination rate constant (1/day)") # Table S2, fitted, 95% CI [0.11; 0.15]
    Imaxaln <- 0.78; label("Maximum fractional inhibition of PCSK9 synthesis by ALN-PCS (unitless)") # Table S2, fitted, 95% CI [0.7; 0.84]
    ID50aln <- 2.55; label("Liver ALN-PCS amount giving half-maximal PCSK9 synthesis inhibition (mg)") # Table S2, fitted, 95% CI [1.74; 3.75]
  })

  model({
    # ===================================================================
    # 1. Modality-specific LDL-receptor exponent. Supplemental Table S2
    #    gives n1 = 0.14 for all three mAbs and n1 = 0.26 for both siRNAs.
    # ===================================================================
    n1 <- n1mab * (1 - TRT_ANTIPCSK9_SIRNA) + n1sirna * TRT_ANTIPCSK9_SIRNA

    # ===================================================================
    # 2. VLDL-C is held at its baseline. The Discussion is explicit:
    #    "Because PCSK9 affected the clearance of LDL particles only, the
    #    amount of VLDL particles, and therefore their cholesterol
    #    fraction, were not affected by anti-PCSK9 treatment. Hence, the
    #    VLDL-C concentration was fixed at baseline level and did not
    #    change in response to treatment."
    # ===================================================================
    VLDLc <- BaselineVLDLc

    # ===================================================================
    # 3. Dose unit conversion, equation 10:
    #      dose[nmol] = 10^6 * dose[mg] / MW[g/mol]
    #    Applied as a bioavailability scalar so that all five drugs are
    #    dosed in mg. The two siRNA chains are already in mg (equations
    #    21-24) and need no conversion.
    # ===================================================================
    f(depot_aliro) <- 1e6 / MWaliro
    f(depot_evolo) <- 1e6 / MWevolo
    f(depot_rg) <- 1e6 / MWrg

    # ===================================================================
    # 4. Plasma PCSK9, equation 25 (equation 20 without the siRNA terms,
    #    equation 1 without any drug). Synthesis is the baseline amount
    #    times the degradation constant, multiplicatively inhibited by
    #    each siRNA; loss is first-order degradation plus 1:1 binding to
    #    each mAb, with the reverse reaction returning free PCSK9.
    #
    #    CORRECTED. As printed, the evolocumab terms of equations 20 and
    #    25 are malformed: the association term
    #    -konevolo*(Ac_evolo/Vd_evolo)*PCSK9 is missing, the dissociation
    #    term carries a minus rather than a plus, and a complex-clearance
    #    term -(CL_evolo/Vd_evolo)*comE appears that belongs to the comE
    #    balance (equation 16), not to the PCSK9 balance. The alirocumab
    #    and RG-7652 terms of the same equations are well formed and
    #    mutually symmetric, and the repair below simply makes evolocumab
    #    match them. The vignette shows that this reading reproduces the
    #    printed evolocumab PCSK9 statistics of Table 1 to within 0.4
    #    percentage points; see the vignette Errata.
    # ===================================================================
    d/dt(pcsk9) <- kPCSK9deg * BaselinePCSK9 *
      (1 - Imaxaln * liver_aln / (liver_aln + ID50aln)) *
      (1 - Imaxinc * liver_inc / (liver_inc + ID50inc)) -
      kPCSK9deg * pcsk9 -
      konaliro * (central_aliro / Vdaliro) * pcsk9 + konaliro * Kdaliro * complex_aliro -
      konevolo * (central_evolo / Vdevolo) * pcsk9 + konevolo * Kdevolo * complex_evolo -
      konrg * (central_rg / Vdrg) * pcsk9 + konrg * Kdrg * complex_rg

    # ===================================================================
    # 5. LDL-C, equation 2. Production is proportional to VLDL-C relative
    #    to its baseline; clearance is first-order and scaled by the
    #    relative free LDL-receptor number.
    # ===================================================================
    d/dt(ldl) <- kLDLcdeg * BaselineLDLc * (VLDLc / BaselineVLDLc) -
      kLDLcdeg * ldl * ldlr

    # ===================================================================
    # 6. Free LDL-receptor ratio, equation 3. Zero-order appearance at
    #    the turnover constant, degradation accelerated by PCSK9 above
    #    baseline (power n1) and by LDL-C above baseline (power n2, the
    #    endocytic-internalization feedback). Starts at 1 by construction.
    # ===================================================================
    d/dt(ldlr) <- kLDLrturn -
      kLDLrturn * ldlr * (pcsk9 / BaselinePCSK9)^n1 * (ldl / BaselineLDLc)^n2

    # ===================================================================
    # 7. Lp(a), equation 4. Parallel LDL-receptor-independent and
    #    LDL-receptor-dependent clearance; synthesis is set so that the
    #    state holds at baseline when ldlr is 1.
    # ===================================================================
    d/dt(lpa) <- (kLpAdeg + kLpAdeg2) * BaselineLpA -
      (kLpAdeg + kLpAdeg2 * ldlr) * lpa

    # ===================================================================
    # 8. Alirocumab, equations 11-13.
    # ===================================================================
    d/dt(depot_aliro) <- -kabsaliro * depot_aliro
    d/dt(central_aliro) <- kabsaliro * depot_aliro -
      CLaliro * (central_aliro / Vdaliro) -
      konaliro * (central_aliro / Vdaliro) * pcsk9 +
      konaliro * Kdaliro * complex_aliro
    d/dt(complex_aliro) <- konaliro * (central_aliro / Vdaliro) * pcsk9 -
      konaliro * Kdaliro * complex_aliro -
      (CLaliro / Vdaliro) * complex_aliro

    # ===================================================================
    # 9. Evolocumab, equations 14-16. CORRECTED: equation 16 is printed
    #    with a '+' where its '=' belongs and misspells the drug
    #    subscript as 'evolvo'. Written here symmetric with the
    #    alirocumab complex balance of equation 13.
    # ===================================================================
    d/dt(depot_evolo) <- -kabsevolo * depot_evolo
    d/dt(central_evolo) <- kabsevolo * depot_evolo -
      CLevolo * (central_evolo / Vdevolo) -
      konevolo * (central_evolo / Vdevolo) * pcsk9 +
      konevolo * Kdevolo * complex_evolo
    d/dt(complex_evolo) <- konevolo * (central_evolo / Vdevolo) * pcsk9 -
      konevolo * Kdevolo * complex_evolo -
      (CLevolo / Vdevolo) * complex_evolo

    # ===================================================================
    # 10. RG-7652, equations 17-19.
    # ===================================================================
    d/dt(depot_rg) <- -kabsrg * depot_rg
    d/dt(central_rg) <- kabsrg * depot_rg -
      CLrg * (central_rg / Vdrg) -
      konrg * (central_rg / Vdrg) * pcsk9 +
      konrg * Kdrg * complex_rg
    d/dt(complex_rg) <- konrg * (central_rg / Vdrg) * pcsk9 -
      konrg * Kdrg * complex_rg -
      (CLrg / Vdrg) * complex_rg

    # ===================================================================
    # 11. Inclisiran, equations 21-22, and ALN-PCS, equations 23-24.
    #     CORRECTED: equation 24's left-hand side is printed as
    #     dAc_inc/dt while every term on its right-hand side carries the
    #     'aln' subscript; it is the ALN-PCS circulating balance.
    # ===================================================================
    d/dt(depot_inc) <- -kabsinc * depot_inc
    d/dt(liver_inc) <- kabsinc * depot_inc - kelinc * liver_inc
    d/dt(depot_aln) <- -kabsaln * depot_aln
    d/dt(liver_aln) <- kabsaln * depot_aln - kelaln * liver_aln

    # ===================================================================
    # 12. Initial conditions. Every endogenous state starts at its
    #     baseline, which the Methods justify: "because under steady-state
    #     conditions system variables do not change over time, and
    #     corresponding values can thus be fixed according to respective
    #     baseline levels". The LDL-receptor state is a ratio to baseline
    #     and therefore starts at 1.
    # ===================================================================
    pcsk9(0) <- BaselinePCSK9
    ldl(0) <- BaselineLDLc
    ldlr(0) <- 1
    lpa(0) <- BaselineLpA

    # ===================================================================
    # 13. Algebraic readouts.
    #
    #     Triglyceride, equation 5, is a ratio to baseline and equals 1
    #     at baseline because lamtgVLDL + lamtgLDL = 1.
    #
    #     HDL-C, equation 6, is CORRECTED. As printed,
    #       HDLc = HDLc_bl * (1 - lamtg * LDLc/Baseline_LDLc),
    #     which returns 0.66 * HDLc_bl at baseline rather than HDLc_bl,
    #     and moves HDL-C DOWN when LDL-C falls. The Methods describe
    #     HDL-C as "inversely related to plasma TGs", equation 5 defines
    #     exactly such a baseline-referenced TG ratio, and the Discussion
    #     states "HDL-C increases by 5% to 10% in virtually every
    #     anti-PCSK9 trial". The form below is the reading that satisfies
    #     all three: it returns HDLc_bl at baseline and yields +5.8% at
    #     the deepest simulated LDL-C reduction. See the vignette Errata.
    #
    #     apoB, equation 9, uses the sign-corrected lamapoB so that apoB
    #     at baseline is lamapoB * (baseline non-HDL-C), a positive mass
    #     concentration.
    # ===================================================================
    TG <- (VLDLc / BaselineVLDLc) * lamtgVLDL + (ldl / BaselineLDLc) * lamtgLDL
    HDLc <- BaselineHDLc * (1 - lamtg * (TG - 1))
    TC <- HDLc + ldl + VLDLc
    NonHDLc <- TC - HDLc
    ApoB <- lamapoB * (BaselineLDLc + BaselineVLDLc) *
      (lamapoBLDL * (ldl / BaselineLDLc) + lamapoBVLDL * (VLDLc / BaselineVLDLc))

    # Plasma concentrations. The Methods state that "PCSK9 plasma
    # concentrations were subsequently calculated by dividing the derived
    # PCSK9 quantities by plasma volume (2.75 l)"; the same division gives
    # the free-mAb concentrations that drive the binding terms above.
    Cpcsk9 <- pcsk9 / Vpl
    LDLc <- ldl
    LpAc <- lpa
    Caliro <- central_aliro / Vdaliro
    Cevolo <- central_evolo / Vdevolo
    Crg <- central_rg / Vdrg

    # No residual-error model and no IIV. The source fitted study-level
    # aggregated means by a nonlinear FIXED-effects method and says so
    # explicitly: "such clinical data are often reported at an aggregated
    # trial level ... This prevents us from applying mixed-effects modeling
    # and evaluating between-subject variability." The 95% confidence
    # intervals in supplemental Table S2 are parameter-uncertainty bands
    # from the Fisher information matrix and likelihood profiling, not
    # random-effect magnitudes.
  })
}
