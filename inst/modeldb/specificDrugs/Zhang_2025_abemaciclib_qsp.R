Zhang_2025_abemaciclib_qsp <- function() {
  description <- paste(
    "QSP. CDK4/6 target-engagement and pRB / TOPO-IIa cell-cycle biomarker model",
    "for abemaciclib (ABE) and its three active metabolites M2, M18 and M20 in",
    "patients with metastatic breast cancer, including brain metastases. Four",
    "analytes compete reversibly for a shared CDK4 pool and a shared CDK6 pool",
    "(0.1 uM each) in plasma and in cerebrospinal fluid, and the resulting CDK",
    "occupancy fraction drives a four-compartment cell-cycle transit chain",
    "(early G1 precursor, late G1 pRB, S-phase TOPO-IIa, G2M) whose framework",
    "comes from Tate 2014. The abemaciclib and metabolite concentrations are NOT",
    "fitted here: Zhang 2025 generated them with a whole-body PK-Sim 9.1 PBPK",
    "model that is a platform port (system-specific parameters were taken from",
    "the built-in PK-Sim database, distribution used the built-in Rodgers and",
    "Rowland method, no organ ODEs / volumes / blood flows are published and no",
    ".pksim5 project was deposited), so that layer is not reproducible from the",
    "on-disk sources and is deliberately NOT extracted. Total plasma",
    "concentrations are instead supplied per record as the canonical",
    "time-varying covariates CP_ABE_NGML, CP_M2_NGML, CP_M18_NGML and",
    "CP_M20_NGML, and the albumin-scaled fraction unbound is applied inside",
    "model(), following the Liang_2024_osimertinib_qsp precedent. Deterministic",
    "mechanism model: Zhang 2025 propagated variability through PK-Sim virtual",
    "population physiology rather than fitting IIV or residual error, so no etas",
    "and no error model are encoded. The CSF limb applies the fixed Table 1",
    "CSF-to-plasma ratio K_CSF,p to the free plasma concentration, exactly as",
    "Methods 2.1 describes. Validated against the paper's own prose claims, the",
    "model reproduces five of the six numeric statements in Section 3.6,",
    "including both 150 mg BID claims (CSF CDK4 above 90%, CSF CDK6 above 80%);",
    "the one it misses is CSF CDK6 above 90% at 200 mg BID, where it returns",
    "about 86%. That shortfall is small and is not necessarily structural: the",
    "analyte concentrations used to score it are reconstructed from Table 2 plus",
    "the published metabolite mass split rather than taken from the",
    "un-extracted PBPK layer, and that reconstruction is itself about 17% low",
    "against the paper's own quoted total-analyte trough. See the vignette",
    "'Assumptions and deviations' section, which quantifies both candidate",
    "causes.",
    sep = " "
  )
  reference <- paste(
    "Zhang C, Li S, Ren J, Lang R. (2025). Physiologically Based Pharmacokinetic",
    "Model of Plasma and Intracranial Pharmacokinetics and CDK4/6 Occupancy of",
    "Abemaciclib to Optimizing Dosing Regimen for Brain Metastatic Patients.",
    "ACS Omega 10(9):9245-9256. doi:10.1021/acsomega.4c09472. PMCID PMC11904693.",
    "Target-engagement equation from Methods 2.1 (Eq 1, after Wong 2019);",
    "fraction-unbound albumin scaling from Methods 2.1 (Eq 2, after Alsmadi",
    "2021); biomarker chain from Methods 2.2 (Eqs 4-7, framework from Tate",
    "2014). All parameter values are from main-text Table 1. Note that the",
    "trimmed markdown companion of this paper retains Table 1 and Section 3.6",
    "but drops every display equation (eight 'formula-not-decoded' markers), so",
    "Eqs 1-7 were recovered from the PDF.",
    sep = " "
  )
  vignette <- "Zhang_2025_abemaciclib"

  units <- list(
    time          = "h",
    dosing        = "(not applicable; the analyte concentrations are supplied as the time-varying covariates CP_ABE_NGML, CP_M2_NGML, CP_M18_NGML and CP_M20_NGML)",
    concentration = "nmol/L"
  )

  # The bound drug-CDK complexes are indexed by kinase isoform (CDK4 / CDK6),
  # by analyte (parent + three metabolites) and by matrix (plasma / CSF), which
  # the canonical bare `complex` name cannot express. Declared here rather than
  # silently warned on, following the Liang_2024_osimertinib_qsp precedent.
  # `precursor`, `prb` and `topo2a` are the cell-cycle transit states of Eqs 4-6.
  paper_specific_compartments <- c(
    "complex_cdk4_abe",     "complex_cdk4_m2",     "complex_cdk4_m18",     "complex_cdk4_m20",
    "complex_cdk6_abe",     "complex_cdk6_m2",     "complex_cdk6_m18",     "complex_cdk6_m20",
    "complex_cdk4_abe_csf", "complex_cdk4_m2_csf", "complex_cdk4_m18_csf", "complex_cdk4_m20_csf",
    "complex_cdk6_abe_csf", "complex_cdk6_m2_csf", "complex_cdk6_m18_csf", "complex_cdk6_m20_csf",
    "precursor", "prb", "topo2a"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The occupancy states are molar concentrations of the
  # drug-CDK complex (the CDK pools are declared as 0.1 uM concentrations, not
  # amounts). The three biomarker states are dimensionless expression levels
  # relative to the untreated baseline (see the `prb0` note in ini()).
  compartmentData <- list(
    complex_cdk4_abe     = list(analyte = "abemaciclib-CDK4 complex",  units = "nmol/L", specimen = "plasma", verified = TRUE),
    complex_cdk4_m2      = list(analyte = "M2-CDK4 complex",           units = "nmol/L", specimen = "plasma", verified = TRUE),
    complex_cdk4_m18     = list(analyte = "M18-CDK4 complex",          units = "nmol/L", specimen = "plasma", verified = TRUE),
    complex_cdk4_m20     = list(analyte = "M20-CDK4 complex",          units = "nmol/L", specimen = "plasma", verified = TRUE),
    complex_cdk6_abe     = list(analyte = "abemaciclib-CDK6 complex",  units = "nmol/L", specimen = "plasma", verified = TRUE),
    complex_cdk6_m2      = list(analyte = "M2-CDK6 complex",           units = "nmol/L", specimen = "plasma", verified = TRUE),
    complex_cdk6_m18     = list(analyte = "M18-CDK6 complex",          units = "nmol/L", specimen = "plasma", verified = TRUE),
    complex_cdk6_m20     = list(analyte = "M20-CDK6 complex",          units = "nmol/L", specimen = "plasma", verified = TRUE),
    complex_cdk4_abe_csf = list(analyte = "abemaciclib-CDK4 complex",  units = "nmol/L", specimen = "CSF", verified = TRUE),
    complex_cdk4_m2_csf  = list(analyte = "M2-CDK4 complex",           units = "nmol/L", specimen = "CSF", verified = TRUE),
    complex_cdk4_m18_csf = list(analyte = "M18-CDK4 complex",          units = "nmol/L", specimen = "CSF", verified = TRUE),
    complex_cdk4_m20_csf = list(analyte = "M20-CDK4 complex",          units = "nmol/L", specimen = "CSF", verified = TRUE),
    complex_cdk6_abe_csf = list(analyte = "abemaciclib-CDK6 complex",  units = "nmol/L", specimen = "CSF", verified = TRUE),
    complex_cdk6_m2_csf  = list(analyte = "M2-CDK6 complex",           units = "nmol/L", specimen = "CSF", verified = TRUE),
    complex_cdk6_m18_csf = list(analyte = "M18-CDK6 complex",          units = "nmol/L", specimen = "CSF", verified = TRUE),
    complex_cdk6_m20_csf = list(analyte = "M20-CDK6 complex",          units = "nmol/L", specimen = "CSF", verified = TRUE),
    # The biomarker matrix is SKIN, not tumor -- Zhang 2025 says so explicitly in
    # the Discussion: "it is noteworthy that the biomarker model predicts changes
    # in pRB and TOPO-IIa occurring in skin tissue, rather than in tumor tissue.
    # Because the only alterations in pRB and TOPO-IIa within skin tissue in
    # humans has been experimentally determined, it is assumed that the
    # differences between the pRB and TOPO-IIa biomarker expression in tumors and
    # skin tissue are negligible." The controlled specimen vocabulary
    # (R/conventions.R) has no `skin` term, so these carry the generic `tissue`
    # and name skin in `analyte`; recording `tumor` here would contradict the
    # source.
    precursor            = list(analyte = "cells accumulated in early G1 (precursor compartment PC), skin", units = "fraction of untreated baseline", specimen = "tissue", verified = TRUE),
    prb                  = list(analyte = "phosphorylated retinoblastoma protein (pRB) expression, skin",   units = "fraction of untreated baseline", specimen = "tissue", verified = TRUE),
    topo2a               = list(analyte = "topoisomerase-II alpha (TOPO-IIa) expression, skin",             units = "fraction of untreated baseline", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    CP_ABE_NGML = list(
      description        = "Time-varying TOTAL abemaciclib plasma concentration driving CDK4/6 engagement",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL (not free) plasma concentration. The free concentration that",
        "actually drives binding is computed inside model() as",
        "CP_ABE_NGML * fup_abe, with fup_abe scaled to the subject's albumin by",
        "the paper's Eq 2. Zhang 2025 Eq 1 defines its driver as 'the free drug",
        "concentration in plasma or the CSF'. Supplying a total concentration is",
        "the more useful convention because that is what a PK model returns.",
        "In the source study this concentration is produced by a whole-body",
        "PK-Sim 9.1 PBPK model that is not reproducible from the published",
        "sources and is not extracted (see description).",
        "Published steady-state anchors in patients (Zhang 2025 Table 2,",
        "predicted column, Patnaik 2016 regimens): Cmin 144 ng/mL and Cmax",
        "198 ng/mL at 100 mg BID; Cmin 174 and Cmax 250 at 150 mg BID; Cmin 231",
        "and Cmax 333 at 200 mg BID. Total active analytes (ABE + M2 + M20)",
        "reach a mean plasma Cmin of 498 ng/mL at the standard 200 mg BID",
        "regimen, against a stated safety ceiling of 693 ng/mL (Zhang 2025",
        "Sections 2.6 and 3.6)."
      ),
      source_name        = "(none; computed by the PK-Sim PBPK model, not a named NONMEM data column)"
    ),
    CP_M2_NGML = list(
      description        = "Time-varying TOTAL M2 (active abemaciclib metabolite) plasma concentration driving CDK4/6 engagement",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL plasma concentration of the active metabolite M2; the free",
        "concentration is CP_M2_NGML * fup_m2 inside model(). Zhang 2025",
        "Introduction (citing the FDA multidiscipline review, ref 14) reports",
        "that M2 contributes 13% of total plasma mass in vivo. Set to 0 to",
        "study the parent alone."
      ),
      source_name        = "(none; computed by the PK-Sim PBPK model, not a named NONMEM data column)"
    ),
    CP_M18_NGML = list(
      description        = "Time-varying TOTAL M18 (active abemaciclib metabolite) plasma concentration driving CDK4/6 engagement",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL plasma concentration of the active metabolite M18; the free",
        "concentration is CP_M18_NGML * fup_m18 inside model(). Zhang 2025",
        "Introduction (ref 14) reports that M18 contributes 5% of total plasma",
        "mass in vivo. M18 is the metabolite for which no CSF data existed, so",
        "its K_CSF,p was assumed to be 1.0 (Zhang 2025 Methods 2.1). Set to 0",
        "to study the parent alone."
      ),
      source_name        = "(none; computed by the PK-Sim PBPK model, not a named NONMEM data column)"
    ),
    CP_M20_NGML = list(
      description        = "Time-varying TOTAL M20 (active abemaciclib metabolite) plasma concentration driving CDK4/6 engagement",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL plasma concentration of the active metabolite M20; the free",
        "concentration is CP_M20_NGML * fup_m20 inside model(). Zhang 2025",
        "Introduction (ref 14) reports that M20 contributes 26% of total plasma",
        "mass in vivo, making it the largest metabolite contributor. Set to 0",
        "to study the parent alone."
      ),
      source_name        = "(none; computed by the PK-Sim PBPK model, not a named NONMEM data column)"
    ),
    ALB = list(
      description        = "Plasma albumin concentration; scales the fraction unbound of all four analytes through the paper's Eq 2",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Canonical SI units g/L. Zhang 2025 reports albumin in g/dL throughout",
        "(divide by 10): 3.1 g/dL = 31 g/L in breast-cancer patients versus",
        "4.5 g/dL = 45 g/L in healthy subjects (Table 1 and Methods 2.1, after",
        "Alsmadi 2021). The Table 1 fup values are the PATIENT values, i.e.",
        "they already correspond to ALB = 31 g/L, so setting ALB = 31 recovers",
        "them exactly. Zhang 2025 Methods 2.5 swept albumin over 1.0-6.0 g/dL",
        "(10-60 g/L) to reflect the physiological range in cancer patients;",
        "note that in the source paper that sweep changes exposure through the",
        "PBPK layer as well as through binding, whereas here only the binding",
        "limb is reproduced."
      ),
      source_name        = "albumin (g/dl) (Zhang 2025 Table 1, physiological block)"
    )
  )

  population <- list(
    species        = "human (in silico; PK-Sim virtual populations matched to published clinical cohorts)",
    n_subjects     = NA_integer_,
    n_studies      = 4L,
    age_range      = "51-63 years (means of the contributing cohorts; Patnaik 2016 median 60, range 44-73)",
    sex_female_pct = NA_real_,
    disease_state  = "metastatic breast cancer, including hormone-receptor-positive / HER2-negative brain metastases; validation cohorts also included non-small cell lung cancer and other solid tumors",
    dose_range     = "abemaciclib 100-400 mg once daily or twice daily; validation at 100 and 150 mg OD and 100, 150 and 200 mg BID; 200 mg BID is the standard regimen",
    regions        = "United States, Europe and Japan (the contributing clinical studies)",
    notes          = paste(
      "Virtual-population demographics were taken per scenario from the four",
      "clinical plasma-PK studies and two CSF-PK studies tabulated in Zhang",
      "2025 Table 2 (Patnaik 2016, Tolaney 2020, Fujiwara 2016, Kim 2018),",
      "with PK-Sim mean values substituted where a study did not report a",
      "characteristic. Cohort sizes per arm ranged from 3 to 72 subjects;",
      "female proportion ranged from 33% to 100% across arms, so no single",
      "value is recorded here.",
      "Disease physiology differs from healthy subjects in three respects that",
      "Zhang 2025 states explicitly (Methods 2.1, after Alsmadi 2021): CYP3A4",
      "expression 3.02 versus 4.32 uM, albumin 3.1 versus 4.5 g/dL, and",
      "hematocrit 0.33 versus 0.43. Only the albumin difference enters the",
      "layer extracted here, through Eq 2; the other two act on the",
      "un-extracted PBPK layer.",
      "The pRB and TOPO-IIa observations the biomarker model was verified",
      "against are SKIN-biopsy measurements from Patnaik 2016, not tumor",
      "biopsies: Zhang 2025 states in the Discussion that the biomarker model",
      "'predicts changes in pRB and TOPO-IIa occurring in skin tissue, rather",
      "than in tumor tissue', because only the skin-tissue alterations have",
      "been experimentally determined in humans, and that the tumor-versus-skin",
      "difference is assumed negligible. The same Discussion paragraph notes",
      "that the paper does NOT ultimately treat these two biomarkers as",
      "predictive of clinical response (inhibition was of similar magnitude at",
      "150 and 200 mg BID), and uses a CDK4/6 occupancy above 90% as the",
      "efficacy predictor instead."
    )
  )

  ini({
    # ======================================================================
    # Molecular weights -- used only to convert the ng/mL covariates to the
    # nmol/L scale on which Kd and the CDK pools are reported.
    # Zhang 2025 Table 1, physicochemical block (source: ref 31, Posada 2020).
    # ======================================================================
    mw_abe <- fixed(506.6)
    label("Molecular weight of abemaciclib (g/mol)")
    # Zhang 2025 Table 1, MW row, ABE column.

    mw_m2 <- fixed(478.5)
    label("Molecular weight of metabolite M2 (g/mol)")
    # Zhang 2025 Table 1, MW row, M2 column.

    mw_m18 <- fixed(494.5)
    label("Molecular weight of metabolite M18 (g/mol)")
    # Zhang 2025 Table 1, MW row, M18 column.

    mw_m20 <- fixed(522.6)
    label("Molecular weight of metabolite M20 (g/mol)")
    # Zhang 2025 Table 1, MW row, M20 column.

    # ======================================================================
    # Fraction unbound in plasma, at the reference (patient) albumin level.
    # Zhang 2025 Table 1 fup row: "calculated using eq 2, original data was
    # from ref14". These are therefore the PATIENT values, i.e. the values at
    # albumin 3.1 g/dL, and model() re-scales them to any other albumin with
    # the same Eq 2 (see the algebra in the model() comment).
    # ======================================================================
    fup_abe <- fixed(0.039)
    label("Fraction of abemaciclib unbound in plasma at the reference albumin of 3.1 g/dL (unitless)")
    # Zhang 2025 Table 1, f_up row, ABE column.

    fup_m2 <- fixed(0.093)
    label("Fraction of M2 unbound in plasma at the reference albumin of 3.1 g/dL (unitless)")
    # Zhang 2025 Table 1, f_up row, M2 column.

    fup_m18 <- fixed(0.046)
    label("Fraction of M18 unbound in plasma at the reference albumin of 3.1 g/dL (unitless)")
    # Zhang 2025 Table 1, f_up row, M18 column.

    fup_m20 <- fixed(0.030)
    label("Fraction of M20 unbound in plasma at the reference albumin of 3.1 g/dL (unitless)")
    # Zhang 2025 Table 1, f_up row, M20 column.

    alb_ref <- fixed(31)
    label("Reference plasma albumin at which the fup values above apply (g/L; 3.1 g/dL in the source)")
    # Zhang 2025 Table 1, physiological block: "albumin (g/dl) 3.1", described
    # in Methods 2.1 as the mean value in patients (ref 30, Alsmadi 2021).

    # ======================================================================
    # CSF-to-plasma ratio of the FREE concentration. Zhang 2025 Methods 2.1:
    # 0.68 for ABE ("the mean ratio of ABE free concentration between the CSF
    # and plasma", ref 16), 0.14 for M2 and 0.91 for M20 (ref 33), and 1.0
    # assumed for M18 because no data were available.
    # ======================================================================
    kcsf_abe <- fixed(0.68)
    label("CSF-to-plasma free-concentration ratio for abemaciclib (unitless)")
    # Zhang 2025 Table 1, K_CSF,P row, ABE column; Methods 2.1.

    kcsf_m2 <- fixed(0.14)
    label("CSF-to-plasma free-concentration ratio for M2 (unitless)")
    # Zhang 2025 Table 1, K_CSF,P row, M2 column; Methods 2.1.

    kcsf_m18 <- fixed(1.0)
    label("CSF-to-plasma free-concentration ratio for M18 (unitless; assumed, no data available)")
    # Zhang 2025 Methods 2.1: "as no data were available for M18, KCSF,p was
    # assumed to be 1.0 for this metabolite." Table 1, K_CSF,P row, M18 column.

    kcsf_m20 <- fixed(0.91)
    label("CSF-to-plasma free-concentration ratio for M20 (unitless)")
    # Zhang 2025 Table 1, K_CSF,P row, M20 column; Methods 2.1.

    # ======================================================================
    # Equilibrium dissociation constants, one per analyte per kinase isoform.
    # Zhang 2025 Table 1 reports them as paired "CDK4/CDK6" entries; source
    # column: "ref15 for three metabolites and ref35 for ABE".
    # ======================================================================
    kd4_abe <- fixed(0.6)
    label("Equilibrium dissociation constant of abemaciclib for CDK4 (nmol/L)")
    # Zhang 2025 Table 1, Kd for CDK4/6 row, ABE column, first of the pair 0.6/8.2.

    kd6_abe <- fixed(8.2)
    label("Equilibrium dissociation constant of abemaciclib for CDK6 (nmol/L)")
    # Zhang 2025 Table 1, Kd for CDK4/6 row, ABE column, second of the pair 0.6/8.2.

    kd4_m2 <- fixed(1.2)
    label("Equilibrium dissociation constant of M2 for CDK4 (nmol/L)")
    # Zhang 2025 Table 1, Kd for CDK4/6 row, M2 column, pair 1.2/1.3.

    kd6_m2 <- fixed(1.3)
    label("Equilibrium dissociation constant of M2 for CDK6 (nmol/L)")
    # Zhang 2025 Table 1, Kd for CDK4/6 row, M2 column, pair 1.2/1.3.

    kd4_m18 <- fixed(1.2)
    label("Equilibrium dissociation constant of M18 for CDK4 (nmol/L)")
    # Zhang 2025 Table 1, Kd for CDK4/6 row, M18 column, pair 1.2/2.7.

    kd6_m18 <- fixed(2.7)
    label("Equilibrium dissociation constant of M18 for CDK6 (nmol/L)")
    # Zhang 2025 Table 1, Kd for CDK4/6 row, M18 column, pair 1.2/2.7.

    kd4_m20 <- fixed(1.5)
    label("Equilibrium dissociation constant of M20 for CDK4 (nmol/L)")
    # Zhang 2025 Table 1, Kd for CDK4/6 row, M20 column, pair 1.5/1.9.

    kd6_m20 <- fixed(1.9)
    label("Equilibrium dissociation constant of M20 for CDK6 (nmol/L)")
    # Zhang 2025 Table 1, Kd for CDK4/6 row, M20 column, pair 1.5/1.9.

    # ======================================================================
    # Dissociation rate constant and starting CDK expression.
    # ======================================================================
    koff <- fixed(6.0)
    label("First-order dissociation rate constant of the drug-CDK complex (1/h)")
    # Zhang 2025 Table 1, koff row: 0.10 /min, "fitted through nonlinear
    # analysis using the data from ref 36 for ABE, three metabolites were
    # assumed to be identical". Converted to the model time base of hours:
    # 0.10 /min * 60 min/h = 6.0 /h. A single koff applies to all four
    # analytes and to both isoforms.

    cdk0 <- fixed(100)
    label("Starting expression of each CDK isoform, CDK4/6_0 (nmol/L)")
    # Zhang 2025 Methods 2.1: "The starting expression of CDK4/6 (CDK4/60) was
    # set at default 0.1 uM due to not reported data in human." 0.1 uM = 100
    # nmol/L. Applied to the CDK4 pool and to the CDK6 pool separately, in
    # plasma and in CSF, because the paper reports the four occupancies
    # separately throughout (Figure 4 legend; every panel of Figure 5).

    # ======================================================================
    # Biomarker (cell-cycle transit) rate constants. Zhang 2025 Table 1,
    # biomarker block, source column: "kTo and kTo were optimized based on the
    # clinical study data, ref16; The remaining rate constants were taken from
    # the ref37" (Tate 2014).
    # ======================================================================
    krb <- fixed(0.16)
    label("Rate constant from early G1 to late G1, k_pRB (1/h)")
    # Zhang 2025 Table 1, kRB row.

    kto <- fixed(0.14)
    label("Rate constant from late G1 to S phase, k_TOPO-IIa (1/h)")
    # Zhang 2025 Table 1, kTO row.

    khh <- fixed(0.21)
    label("Rate constant from S phase to G2M phase, k_HH (1/h)")
    # Zhang 2025 Table 1, kHH row.

    kel_pc <- fixed(0.05)
    label("Elimination rate constant from the early-G1 precursor compartment (1/h)")
    # Zhang 2025 Table 1, kel row (units not printed in Table 1; 1/h by
    # consistency with the other three rate constants in the same block and
    # with Eq 7, in which kel is summed with k_pRB).

    prb0 <- fixed(1)
    label("Baseline pRB expression scale factor (fraction of untreated baseline)")
    # Zhang 2025 Eq 7 defines k_in = pRB0 * (k_pRB + k_el) and says only that
    # "pRB0 is the baseline expression in patients" (ref 16); no numeric value
    # is printed anywhere in the paper or Table 1. Setting it to 1 is EXACT
    # rather than an assumption: Eqs 4-6 are linear in k_in, so PC, pRB and
    # TOPO-IIa all scale proportionally to pRB0, and every quantity the paper
    # reports is a percent-of-baseline ratio (Figure 3 y-axes are "Change in
    # pRB (% of baseline)" and "Change in TOPO-IIa (% of baseline)"), which is
    # invariant to pRB0. Fixing it at 1 makes the states read directly as
    # fractions of the untreated baseline.
  })

  model({
    # =====================================================================
    # 1. Albumin-scaled fraction unbound (Zhang 2025 Eq 2).
    #
    #    As printed, Eq 2 is
    #        fup = 1 / (1 + ((1 - fup_h) * [Alb_p]) / ([Alb_h] * fup_h))
    #    which needs the HEALTHY-subject fup_h, a value the paper does not
    #    print. Rearranging, Eq 2 is equivalent to fup(Alb) = 1/(1 + B*Alb)
    #    with a single analyte-specific binding constant
    #        B = (1 - fup_h) / (fup_h * [Alb_h]),
    #    and B is fully determined by the Table 1 PATIENT value together with
    #    the patient albumin at which it applies:
    #        B = (1/fup_ref - 1) / alb_ref.
    #    So the scaling below uses only printed numbers, and reduces exactly
    #    to the Table 1 fup values when ALB = alb_ref = 31 g/L.
    # =====================================================================
    fu_abe <- 1 / (1 + (1 / fup_abe - 1) / alb_ref * ALB)
    fu_m2  <- 1 / (1 + (1 / fup_m2  - 1) / alb_ref * ALB)
    fu_m18 <- 1 / (1 + (1 / fup_m18 - 1) / alb_ref * ALB)
    fu_m20 <- 1 / (1 + (1 / fup_m20 - 1) / alb_ref * ALB)

    # =====================================================================
    # 2. Free plasma concentrations, converted from the ng/mL covariates to
    #    the nmol/L scale of Kd and cdk0:
    #        C[nmol/L] = C[ng/mL] * fu / MW[g/mol] * 1000
    # =====================================================================
    cf_abe <- CP_ABE_NGML * fu_abe / mw_abe * 1000
    cf_m2  <- CP_M2_NGML  * fu_m2  / mw_m2  * 1000
    cf_m18 <- CP_M18_NGML * fu_m18 / mw_m18 * 1000
    cf_m20 <- CP_M20_NGML * fu_m20 / mw_m20 * 1000

    # =====================================================================
    # 3. Free CSF concentrations. Zhang 2025 Methods 2.1 equates the unbound
    #    brain-interstitial concentration with the free CSF concentration and
    #    uses the Table 1 K_CSF,p ratios. See the DEVIATION note in
    #    `description`: this fixed ratio reproduces the paper at 100 mg BID
    #    but under-predicts CSF occupancy at higher doses.
    # =====================================================================
    ccsf_abe <- kcsf_abe * cf_abe
    ccsf_m2  <- kcsf_m2  * cf_m2
    ccsf_m18 <- kcsf_m18 * cf_m18
    ccsf_m20 <- kcsf_m20 * cf_m20

    # =====================================================================
    # 4. Free (unoccupied) CDK available in each pool. Zhang 2025 Methods 2.1
    #    states that the four analytes are "simulated in the same manner
    #    simultaneously, and the sum of the occupancy by ABE and metabolites
    #    gave the total CDK4/6 occupancy", i.e. they COMPETE for one shared
    #    pool per isoform per matrix. Independent per-analyte pools are ruled
    #    out arithmetically: summing four independently-computed occupancies
    #    exceeds 100% by a wide margin at the paper's own reported exposures,
    #    whereas the shared-pool form below is bounded by construction and
    #    reproduces both of the Section 3.6 prose claims at 150 mg BID.
    # =====================================================================
    free_cdk4 <- cdk0 - (complex_cdk4_abe + complex_cdk4_m2 + complex_cdk4_m18 + complex_cdk4_m20)
    free_cdk6 <- cdk0 - (complex_cdk6_abe + complex_cdk6_m2 + complex_cdk6_m18 + complex_cdk6_m20)

    free_cdk4_csf <- cdk0 - (complex_cdk4_abe_csf + complex_cdk4_m2_csf + complex_cdk4_m18_csf + complex_cdk4_m20_csf)
    free_cdk6_csf <- cdk0 - (complex_cdk6_abe_csf + complex_cdk6_m2_csf + complex_cdk6_m18_csf + complex_cdk6_m20_csf)

    # =====================================================================
    # 5. Target engagement. Zhang 2025 Eq 1, recovered from the PDF (the
    #    trimmed markdown drops it as `formula-not-decoded`):
    #        dN/dt = (koff / Kd) * CDK_unbound * C_drug - koff * CDK_bound
    #    where N is the drug-CDK complex. koff/Kd is the association rate
    #    constant. Applied per analyte, per isoform, in plasma and in CSF.
    # =====================================================================
    d/dt(complex_cdk4_abe) <- koff / kd4_abe * free_cdk4 * cf_abe - koff * complex_cdk4_abe
    d/dt(complex_cdk4_m2)  <- koff / kd4_m2  * free_cdk4 * cf_m2  - koff * complex_cdk4_m2
    d/dt(complex_cdk4_m18) <- koff / kd4_m18 * free_cdk4 * cf_m18 - koff * complex_cdk4_m18
    d/dt(complex_cdk4_m20) <- koff / kd4_m20 * free_cdk4 * cf_m20 - koff * complex_cdk4_m20

    d/dt(complex_cdk6_abe) <- koff / kd6_abe * free_cdk6 * cf_abe - koff * complex_cdk6_abe
    d/dt(complex_cdk6_m2)  <- koff / kd6_m2  * free_cdk6 * cf_m2  - koff * complex_cdk6_m2
    d/dt(complex_cdk6_m18) <- koff / kd6_m18 * free_cdk6 * cf_m18 - koff * complex_cdk6_m18
    d/dt(complex_cdk6_m20) <- koff / kd6_m20 * free_cdk6 * cf_m20 - koff * complex_cdk6_m20

    d/dt(complex_cdk4_abe_csf) <- koff / kd4_abe * free_cdk4_csf * ccsf_abe - koff * complex_cdk4_abe_csf
    d/dt(complex_cdk4_m2_csf)  <- koff / kd4_m2  * free_cdk4_csf * ccsf_m2  - koff * complex_cdk4_m2_csf
    d/dt(complex_cdk4_m18_csf) <- koff / kd4_m18 * free_cdk4_csf * ccsf_m18 - koff * complex_cdk4_m18_csf
    d/dt(complex_cdk4_m20_csf) <- koff / kd4_m20 * free_cdk4_csf * ccsf_m20 - koff * complex_cdk4_m20_csf

    d/dt(complex_cdk6_abe_csf) <- koff / kd6_abe * free_cdk6_csf * ccsf_abe - koff * complex_cdk6_abe_csf
    d/dt(complex_cdk6_m2_csf)  <- koff / kd6_m2  * free_cdk6_csf * ccsf_m2  - koff * complex_cdk6_m2_csf
    d/dt(complex_cdk6_m18_csf) <- koff / kd6_m18 * free_cdk6_csf * ccsf_m18 - koff * complex_cdk6_m18_csf
    d/dt(complex_cdk6_m20_csf) <- koff / kd6_m20 * free_cdk6_csf * ccsf_m20 - koff * complex_cdk6_m20_csf

    # =====================================================================
    # 6. Total occupancy per isoform per matrix, as a percentage.
    #    Zhang 2025 Methods 2.1: "CDK4/6 occupancy was calculated as the
    #    percentage of the CDK4/6 that was binding to the ABE and metabolites
    #    at a particular time point ... the sum of the occupancy by ABE and
    #    metabolites gave the total CDK4/6 occupancy."
    #    The paper's efficacy threshold is a TROUGH occupancy of at least 90%
    #    (Methods 2.6).
    # =====================================================================
    occCdk4Plasma <- 100 * (cdk0 - free_cdk4)     / cdk0
    occCdk6Plasma <- 100 * (cdk0 - free_cdk6)     / cdk0
    occCdk4Csf    <- 100 * (cdk0 - free_cdk4_csf) / cdk0
    occCdk6Csf    <- 100 * (cdk0 - free_cdk6_csf) / cdk0

    # =====================================================================
    # 7. Biomarker (cell-cycle transit) model. Zhang 2025 Eqs 4-6, with the
    #    zero-order input of Eq 7:
    #        dPC/dt   = k_in - k_el*PC - k_pRB*PC*(1 - CO)
    #        dpRB/dt  = k_pRB*PC*(1 - CO) - k_TOPO*pRB
    #        dTOPO/dt = k_TOPO*pRB - k_HH*TOPO
    #        k_in     = pRB0*(k_pRB + k_el)
    #
    #    WHICH occupancy is CO? The paper's only definition is the verbatim
    #    sentence after Eq 6, "CO is the CDK occupancy fraction", and it never
    #    says which of the four occupancy traces it computes (CDK4 / CDK6,
    #    plasma / CSF) is meant -- the quantity is under-determined by a factor
    #    of four. Resolved to CSF CDK6 (operator decision, task
    #    oare_PMC11904693 sidecar request-001 answer D, 2026-08-20) on a
    #    reproduction check rather than on preference:
    #
    #      * Eqs 4-6 have the exact drug-free-vs-constant-CO steady state
    #            pRB/pRB_baseline = (1-CO)*(kel + kpRB) / (kel + kpRB*(1-CO))
    #        Back-solving Figure 3's predicted plateaus through it gives the CO
    #        the authors actually used: 88.2% at 150 mg BID, 91.1% at 200 mg
    #        BID.
    #      * Section 3.6 states in words that at 150 mg BID "CDK4 occupancy
    #        exceeds 90% and CDK6 occupancy exceeds 80% (Figure 5E)", and the
    #        Figure 5 caption assigns panel E to CSF. 88.2% is above 80% but
    #        below 90%, so of the four traces only CSF CDK6 is consistent with
    #        the back-solve -- and this argument rests on the paper's own prose,
    #        not on digitising a figure.
    #      * The rival reading (CO is a single pooled "CDK4/6" quantity, since
    #        the paper writes that phrase throughout and Table 1 gives Kd as one
    #        paired "CDK4/6" row) predicts ~96%, which would put the pRB plateau
    #        near 11-15% of baseline against Figure 3's ~36%.
    #      * Site: Section 3.6 opens "Intracranial CDK4/6 occupancy serves as
    #        predictors of clinical efficacy", so CSF is the paper's own
    #        efficacy driver.
    #
    #    See the vignette "Assumptions and deviations" section, which tabulates
    #    all four candidate traces against the back-solve so a reader can redo
    #    the check.
    # =====================================================================
    co <- occCdk6Csf / 100

    kin <- prb0 * (krb + kel_pc)

    d/dt(precursor) <- kin - kel_pc * precursor - krb * precursor * (1 - co)
    d/dt(prb)       <- krb * precursor * (1 - co) - kto * prb
    d/dt(topo2a)    <- kto * prb - khh * topo2a

    # Untreated (CO = 0) steady state of Eqs 4-6, so the system starts at the
    # patient's drug-free baseline:
    #     PC   = k_in/(k_el + k_pRB) = pRB0
    #     pRB  = k_pRB*PC/k_TOPO
    #     TOPO = k_TOPO*pRB/k_HH
    precursor(0) <- prb0
    prb(0)       <- krb * prb0 / kto
    topo2a(0)    <- krb * prb0 / khh

    # =====================================================================
    # 8. Reported outputs: expression as a percentage of the untreated
    #    baseline, matching the y-axes of Zhang 2025 Figure 3. Deterministic
    #    mechanism model -- the paper reports no residual error and no IIV for
    #    this layer, so none is encoded (see `description`).
    # =====================================================================
    prbPct  <- 100 * prb    / (krb * prb0 / kto)
    topoPct <- 100 * topo2a / (krb * prb0 / khh)

    # Mass-balance check quantities: each equals cdk0 at all times.
    totalCdk4Plasma <- free_cdk4     + complex_cdk4_abe     + complex_cdk4_m2     + complex_cdk4_m18     + complex_cdk4_m20
    totalCdk6Plasma <- free_cdk6     + complex_cdk6_abe     + complex_cdk6_m2     + complex_cdk6_m18     + complex_cdk6_m20
    totalCdk4Csf    <- free_cdk4_csf + complex_cdk4_abe_csf + complex_cdk4_m2_csf + complex_cdk4_m18_csf + complex_cdk4_m20_csf
    totalCdk6Csf    <- free_cdk6_csf + complex_cdk6_abe_csf + complex_cdk6_m2_csf + complex_cdk6_m18_csf + complex_cdk6_m20_csf
  })
}
