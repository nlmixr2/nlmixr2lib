Zhang_2026_ribociclib_qsp <- function() {
  description <- paste(
    "QSP. CDK4/6 target-engagement model for ribociclib (RIB) in the",
    "cerebrospinal fluid of patients with breast cancer and brain metastasis",
    "(BCBM). Ribociclib binds reversibly to a CDK4 pool and to a CDK6 pool",
    "(0.1 uM each) in CSF, and the resulting occupancy is the paper's efficacy",
    "readout: a CDK4/6 occupancy of at least 90% is the stated therapeutic",
    "target. Zhang 2026 includes NO tumor compartment -- Methods states",
    "explicitly that 'we used the free RIB concentration in CSF to calculate",
    "CDK4/6 occupancy' -- so occupancy is computed in CSF only, unlike the",
    "plasma-plus-CSF structure of the same team's abemaciclib model",
    "(Zhang_2025_abemaciclib_qsp).",
    "The ribociclib concentration is NOT fitted here: Zhang 2026 generated it",
    "with a whole-body PK-Sim 9.1 PBPK model that is a platform port",
    "(system-specific parameters came from the built-in PK-Sim database,",
    "distribution used the built-in Rodgers and Rowland method, and no organ",
    "ODEs, organ volumes, blood flows, per-organ partition coefficients or",
    ".pksim5 project are published -- the paper reports no volume term of any",
    "kind), so that layer is not reproducible from the on-disk sources and is",
    "deliberately NOT extracted. Total plasma concentration is instead supplied",
    "per record as the canonical time-varying covariate CP_RIB_NGML, and the",
    "fraction unbound and the CSF-to-plasma unbound ratio are applied inside",
    "model(), following the Zhang_2025_abemaciclib_qsp and",
    "Liang_2024_osimertinib_qsp precedents.",
    "Deterministic mechanism model: Zhang 2026 propagated variability through",
    "PK-Sim virtual-population physiology rather than fitting IIV or residual",
    "error, so no etas and no error model are encoded.",
    "Validated against the paper's own published predicted free-CSF",
    "concentrations (Table 2) the model reproduces both dosing claims that the",
    "paper states unambiguously: at 600 mg OD it returns CDK4 95.0% and CDK6",
    "82.9%, against the paper's 'CDK4 occupancy exceeded 90%, CDK6 stayed below",
    "90% but above 70%'; and at 900 mg OD it returns CDK4 98.1% and CDK6 92.8%,",
    "against the paper's 'both CDK4 and CDK6 occupancies surpassed 90%'. At",
    "400 mg OD it returns CDK6 73.9%, consistent with the paper's 'failed to",
    "reach 90%', but CDK4 91.7%, which marginally exceeds 90% where the paper",
    "reports a failure; see the vignette 'Assumptions and deviations' section,",
    "which quantifies that one discrepancy and notes that the anchor",
    "concentration comes from a glioblastoma validation cohort rather than the",
    "BCBM simulation cohort of Figure 3.",
    sep = " "
  )
  reference <- paste(
    "Zhang C, Li S, Ren J, Wen X. (2026). Optimizing ribociclib dosing in",
    "breast cancer with brain metastasis patients using a physiologically",
    "based pharmacokinetic model. BMC Cancer 26:204.",
    "doi:10.1186/s12885-026-15561-x. PMCID PMC12888661.",
    "Target-engagement equation from Methods, Eq 1 (attributed there to refs 17",
    "and 31). All parameter values are from main-text Table 1 except the",
    "starting CDK4/6 expression, which is given in the Methods prose beneath",
    "Eq 1, and the CSF-to-plasma unbound ratio Kp,uu, which is reported in the",
    "Results as a predicted quantity per dose level. Note that the trimmed",
    "markdown companion of this paper renders Eq 1 as a",
    "'formula-not-decoded' marker, so Eq 1 was recovered from the PDF.",
    "The single supplementary file (MOESM1) is one figure (S1, simulated",
    "intracranial CDK4/6 occupancy by regimen) and contains no parameter",
    "table.",
    sep = " "
  )
  vignette <- "Zhang_2026_ribociclib"

  units <- list(
    time          = "h",
    dosing        = "(not applicable; the ribociclib concentration is supplied as the time-varying covariate CP_RIB_NGML)",
    concentration = "nmol/L"
  )

  # The bound drug-CDK complexes are indexed by kinase isoform (CDK4 / CDK6)
  # and are specific to the CSF matrix, which the canonical bare `complex` name
  # cannot express. Declared here rather than silently warned on, following the
  # Zhang_2025_abemaciclib_qsp precedent.
  paper_specific_compartments <- c(
    "complex_cdk4_csf",
    "complex_cdk6_csf"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Both states are molar concentrations of the drug-CDK
  # complex (the CDK pools are declared as 0.1 uM concentrations, not amounts).
  compartmentData <- list(
    complex_cdk4_csf = list(analyte = "ribociclib-CDK4 complex", units = "nmol/L", specimen = "CSF", verified = TRUE),
    complex_cdk6_csf = list(analyte = "ribociclib-CDK6 complex", units = "nmol/L", specimen = "CSF", verified = TRUE)
  )

  covariateData <- list(
    CP_RIB_NGML = list(
      description        = "Time-varying TOTAL ribociclib plasma concentration driving CDK4/6 engagement in CSF",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL (not free) plasma concentration. The free CSF concentration that",
        "actually drives binding is computed inside model() as",
        "CP_RIB_NGML * fup * kpuu, because Zhang 2026 defines Kp,uu as the",
        "'CSF-to-plasma unbound RIB concentration ratio'. Eq 1 defines its own",
        "driver as 'the concentration of free RIB present in CSF'; supplying a",
        "total plasma concentration is the more useful convention because that",
        "is what a PK model returns, and it matches the",
        "Zhang_2025_abemaciclib_qsp precedent.",
        "In the source study this concentration is produced by a whole-body",
        "PK-Sim 9.1 PBPK model that is not reproducible from the published",
        "sources and is not extracted (see description).",
        "Published steady-state anchors in patients (Zhang 2026 Table 2,",
        "PREDICTED column, Infante 2016 regimens, Day 18): Cmin 524 ng/mL and",
        "Cmax 1291 ng/mL at 400 mg OD; Cmin 797 and Cmax 2237 at 600 mg OD;",
        "Cmin 1456 and Cmax 3720 at 900 mg OD. Day 1 predicted Cmax is also",
        "tabulated for 750 mg (1240) and 1200 mg (2298) OD. Zhang 2026 reports",
        "no plasma profile for any twice-daily regimen, so the 200 / 300 /",
        "400 mg BID arms of its Figure 3 cannot be driven from published",
        "numbers alone."
      ),
      source_name        = "(none; computed by the PK-Sim PBPK model, not a named NONMEM data column)"
    )
  )

  population <- list(
    species        = "human (in silico; PK-Sim virtual populations matched to published clinical cohorts)",
    n_subjects     = NA_integer_,
    n_studies      = 4L,
    age_range      = "median 49 to 65 years across the contributing cohorts (Curigliano 57 and 65; Infante 60; Tien 49; Johnson 53)",
    sex_female_pct = NA_real_,
    disease_state  = "breast cancer with brain metastasis (BCBM) for the simulation population; the validation cohorts were HR-positive / HER2-negative breast cancer, Rb-positive breast cancer plus liposarcoma and colon cancer, recurrent glioblastoma, and glioblastoma",
    dose_range     = "ribociclib 200 to 1200 mg once daily or twice daily; plasma validation at 400, 600, 750, 900 and 1200 mg OD; CSF validation at 400, 600 and 900 mg OD; 600 mg OD is the standard approved regimen and 300 mg BID is the paper's recommended optimum for BCBM",
    regions        = NA_character_,
    notes          = paste(
      "Virtual-population demographics (sample size, age distribution, female",
      "proportion) were taken per scenario from the four clinical studies",
      "tabulated in Zhang 2026 Table 2 (Curigliano, Infante, Tien, Johnson),",
      "with PK-Sim database mean values substituted where a study did not",
      "report a characteristic. Cohort sizes per arm ranged from 3 to 24",
      "subjects. All simulations were performed in the fasting state.",
      "The virtual population used for the dose-optimisation simulations of",
      "Figure 3 was demographically aligned with Johnson's cohort.",
      "Disease physiology differs from healthy subjects in three respects that",
      "Zhang 2026 states explicitly (Methods): hepatic CYP3A4 expression 3.02",
      "versus roughly 5.5 uM (an approximately 45% downregulation), albumin",
      "3.1 versus 4.5 g/dL, and hematocrit 0.33 versus 0.43. Gastric emptying",
      "time was also set to a patient value of 190 min against 57-72 min in",
      "healthy individuals. NONE of these four enters the layer extracted",
      "here: all act on the un-extracted PBPK layer. The albumin difference",
      "reaches this layer only indirectly, through the Table 1 fraction",
      "unbound, which is itself a measured breast-cancer-patient value.",
      "There is no experimentally determined CDK4/6 occupancy in humans;",
      "Zhang 2026 states this as a limitation, and the starting CDK4/6",
      "expression is a default rather than a measurement."
    )
  )

  ini({
    # ======================================================================
    # Molecular weight -- used only to convert the ng/mL covariate to the
    # nmol/L scale on which Ki and the CDK pools are reported.
    # ======================================================================
    mw_rib <- fixed(434.6)
    label("Molecular weight of ribociclib (g/mol)")
    # Zhang 2026 Table 1, physicochemical block, MW row (source: ref 18).

    # ======================================================================
    # Plasma protein binding and CSF partitioning.
    # ======================================================================
    fup <- fixed(0.12)
    label("Fraction of ribociclib unbound in plasma (unitless)")
    # Zhang 2026 Table 1, physicochemical block, f_up row (source: ref 19,
    # Tien et al.). This is a measured value in breast-cancer patients, so it
    # already corresponds to the patient albumin of 3.1 g/dL recorded in the
    # Table 1 physiological block. Unlike the same team's abemaciclib paper,
    # Zhang 2026 publishes NO albumin-scaling equation for f_up, so none is
    # encoded here; the paper's 2.0-6.0 g/dL albumin sweep acts through the
    # un-extracted PBPK layer.

    kpuu <- fixed(1.03)
    label("CSF-to-plasma UNBOUND concentration ratio Kp,uu (unitless)")
    # Zhang 2026 Results, "PBPK model validation for plasma PK": the
    # PBPK-predicted mean Kp,uu values were 0.99 at 400 mg, 1.03 at 600 mg and
    # 1.10 at 900 mg (against observed 1.22, 1.29 and 1.63 respectively).
    # Set to the 600 mg value because 600 mg OD is the standard approved
    # regimen the paper benchmarks against; override in ini() for the other
    # dose levels. NOTE this is the paper's own PREDICTED output rather than a
    # Table 1 input parameter -- it is the one quantity in this file that the
    # un-extracted PBPK layer produced. It is close to unity at every reported
    # dose, so the free CSF concentration is nearly the free plasma
    # concentration.

    # ======================================================================
    # Equilibrium dissociation constants, one per kinase isoform.
    # ======================================================================
    kd4 <- fixed(10)
    label("Equilibrium dissociation constant of ribociclib for CDK4 (nmol/L)")
    # Zhang 2026 Table 1, CDK4/6 occupancy block, "K i (nM) 10/39" (source:
    # ref 24), first of the pair. The paper labels this Ki, but its own gloss
    # beneath Eq 1 states that "koff/Ki is equivalent to the association rate
    # constant (kon) of RIB to CDK4/6", which makes Ki an equilibrium
    # DISSOCIATION constant (Ki = koff/kon = Kd); it is therefore recorded
    # under the canonical `kd` name, matching Zhang_2025_abemaciclib_qsp.
    # Distinct from the Table 1 "auto-inhibition against CYP3A4, Ki 8.6 uM",
    # which is a true inhibition constant and belongs to the un-extracted
    # PBPK layer.

    kd6 <- fixed(39)
    label("Equilibrium dissociation constant of ribociclib for CDK6 (nmol/L)")
    # Zhang 2026 Table 1, CDK4/6 occupancy block, "K i (nM) 10/39", second of
    # the pair. Ribociclib is roughly 4-fold less potent against CDK6 than
    # against CDK4, which is the whole reason the paper's CDK6 occupancy
    # lags its CDK4 occupancy at every dose.

    # ======================================================================
    # Dissociation rate constant and starting CDK expression.
    # ======================================================================
    koff <- fixed(3.78)
    label("First-order dissociation rate constant of the drug-CDK complex (1/h)")
    # Zhang 2026 Table 1, CDK4/6 occupancy block, koff row: 0.063 /min,
    # "the fitting was performed using nonlinear analysis and the data were
    # sourced from Ref. 25". Converted to the model time base of hours:
    # 0.063 /min * 60 min/h = 3.78 /h. A single koff applies to both isoforms.
    # UNITS CONFLICT, resolved in favour of the table: Table 1 gives koff in
    # min^-1 while the prose gloss beneath Eq 1 calls it h^-1. The table is
    # taken as authoritative, and the same team's abemaciclib paper likewise
    # tabulates koff in min^-1 (0.10 /min). The choice does not affect any
    # steady-state result: koff cancels out of the equilibrium occupancy
    # C/(C + Kd) and sets only the rate of approach to it.

    cdk0 <- fixed(100)
    label("Starting expression of each CDK isoform, CDK4/6_0 (nmol/L)")
    # Zhang 2026 Methods, prose beneath Eq 1: "The initial expression level of
    # CDK4/6 was set to a value of 0.1 uM, given the absence of reported human
    # data." 0.1 uM = 100 nmol/L. Applied to the CDK4 pool and to the CDK6
    # pool separately, because the paper reports the two occupancies
    # separately throughout (Table 3 columns "CDK/4" and "CDK/6"; every panel
    # of Figure 3).
  })

  model({
    # =====================================================================
    # 1. Free ribociclib concentration in CSF.
    #
    #    Zhang 2026 defines Kp,uu as the "CSF-to-plasma unbound RIB
    #    concentration ratio", so
    #        C_CSF,free = Kp,uu * f_up * C_plasma,total
    #    and the covariate carries the TOTAL plasma concentration.
    #
    #    Cross-check against the paper's own numbers: the Table 2 predicted
    #    steady-state plasma Cmin at 600 mg OD is 797 ng/mL, giving
    #    797 * 0.12 * 1.03 = 98.5 ng/mL free in CSF, against the separately
    #    predicted Johnson-cohort CSF trough of 82.1 ng/mL at the same dose
    #    (20% apart, different virtual populations).
    # =====================================================================
    ccsf_ngml <- CP_RIB_NGML * fup * kpuu

    # Converted to the nmol/L scale of kd4, kd6 and cdk0:
    #     C[nmol/L] = C[ng/mL] / MW[g/mol] * 1000
    ccsf <- ccsf_ngml / mw_rib * 1000

    # =====================================================================
    # 2. Free (unoccupied) CDK available in each pool. One analyte, so each
    #    pool has a single binding partner and no competition term.
    # =====================================================================
    free_cdk4 <- cdk0 - complex_cdk4_csf
    free_cdk6 <- cdk0 - complex_cdk6_csf

    # =====================================================================
    # 3. Target engagement. Zhang 2026 Eq 1, recovered from the PDF (the
    #    trimmed markdown drops it as `formula-not-decoded`):
    #        dN/dt = (koff / Ki) * CDK_unbound * C_RIB - koff * CDK_bound
    #    where N is the drug-CDK complex and koff/Ki is the association rate
    #    constant. Applied per isoform, in CSF only -- Zhang 2026 Methods:
    #    "no tumor compartment was included in the model. Instead, we used the
    #    free RIB concentration in CSF to calculate CDK4/6 occupancy."
    # =====================================================================
    d/dt(complex_cdk4_csf) <- koff / kd4 * free_cdk4 * ccsf - koff * complex_cdk4_csf
    d/dt(complex_cdk6_csf) <- koff / kd6 * free_cdk6 * ccsf - koff * complex_cdk6_csf

    # =====================================================================
    # 4. Occupancy per isoform, as a percentage. Zhang 2026 uses a CDK4/6
    #    occupancy of at least 90% as the efficacy threshold (Methods,
    #    "Simulations for optimum dosing regimen for BCBM").
    #
    #    NOTE the paper is internally inconsistent about WHICH occupancy
    #    statistic the 90% threshold applies to: the Methods define the target
    #    as "maximal target occupancy (i.e., peak occupancy levels) at or
    #    above 90%", while the Table 3 sensitivity analysis reports "Tough
    #    [trough] CDK4/6 occupancy in CSF". Both are computable from these
    #    outputs; the vignette reports both and tabulates the consequence.
    # =====================================================================
    occCdk4Csf <- 100 * complex_cdk4_csf / cdk0
    occCdk6Csf <- 100 * complex_cdk6_csf / cdk0

    # Mass-balance check quantities: each equals cdk0 at all times.
    totalCdk4Csf <- free_cdk4 + complex_cdk4_csf
    totalCdk6Csf <- free_cdk6 + complex_cdk6_csf
  })
}
