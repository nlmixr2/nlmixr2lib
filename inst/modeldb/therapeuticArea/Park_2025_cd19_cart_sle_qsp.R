Park_2025_cd19_cart_sle_qsp <- function() {
  description <- paste(
    "QSP. Integrated quantitative systems pharmacology model of FMC63-based",
    "anti-CD19 CAR-T cell therapy repurposed for refractory systemic lupus",
    "erythematosus (Park 2025). 37 ODE states covering (1) cellular kinetics of",
    "CD4+ and CD8+ CAR-T cells, each resolved into naive / central-memory /",
    "effector-memory / terminally-differentiated-effector immunophenotypes",
    "distributing between peripheral blood and bone marrow, (2) explicit",
    "CAR-target complex formation from CAR density, CD19 antigen density and",
    "the FMC63 kon / koff, driving Michaelis-Menten CAR-T expansion and CD19+",
    "B cell killing in both compartments, (3) threshold-gated zero-order",
    "reconstitution of naive (non-autoreactive) CD19+ B cells once total CAR-T",
    "cells contract below CARTrev.min, (4) a K-PD model of fludarabine +",
    "cyclophosphamide lymphodepletion driving precursor-dependent indirect",
    "response models for host CD4+ and CD8+ T cells and a Friberg",
    "myelosuppression chain for neutrophils, and (5) SLE disease biomarkers",
    "(anti-dsDNA autoantibodies, interferon-alpha, complement C3, proteinuria)",
    "coupled to the reactive B cell pool and to the host immune cells. Three",
    "univariate literature regressions translate C3, IFN-alpha and anti-dsDNA",
    "into the SLE disease activity index (SLE-DAI). Deterministic dosing: the",
    "CAR-T product is delivered as compartment initial conditions computed from",
    "dose, body weight and product composition, and the lymphodepletion",
    "chemotherapy amount is the initial condition of the K-PD compartment.",
    sep = " "
  )

  reference <- paste(
    "Park H, Mugundu GM, Singh AP. Mechanistic Evaluation of Anti-CD19 CAR-T",
    "Cell Therapy Repurposed in Systemic Lupus Erythematosus Using a",
    "Quantitative Systems Pharmacology Model. Clin Transl Sci.",
    "2025;18(2):e70146. doi:10.1111/cts.70146",
    sep = " "
  )

  vignette <- "Park_2025_cd19_cart_sle_qsp"

  units <- list(
    time          = "day",
    dosing        = "CAR-positive cells/kg (CAR-T product); mg (lymphodepletion chemotherapy)",
    concentration = "cells/L (CAR-T and B cell states); cells/uL (host lymphocytes and neutrophils)"
  )

  # The 37 mechanistic states are paper-specific: CAR-T immunophenotype pools in
  # two physiological spaces, CAR-target complexes, reactive and reconstituting
  # CD19+ B cell pools, a lymphodepletion K-PD compartment, precursor-dependent
  # host T cell pools, a Friberg neutrophil chain, and four disease biomarkers.
  # None map onto the canonical PK compartment roles.
  paper_specific_compartments <- c(
    "cd8_tn_pb", "cd8_cm_pb", "cd8_em_pb", "cd8_ef_pb",
    "cd4_tn_pb", "cd4_cm_pb", "cd4_em_pb", "cd4_ef_pb",
    "cd8_tn_bm", "cd8_cm_bm", "cd8_em_bm", "cd8_ef_bm",
    "cd4_tn_bm", "cd4_cm_bm", "cd4_em_bm", "cd4_ef_bm",
    "cplx_cd8_pb", "cplx_cd4_pb", "cplx_cd8_bm", "cplx_cd4_bm",
    "cd19b_pb", "cd19b_rev_pb", "cd19b_bm",
    "ldc", "cd4t_pre", "cd4t", "cd8t_pre", "cd8t",
    "nt_prol", "nt_tr1", "nt_tr2", "nt_tr3", "nt",
    "ifna", "aab", "c3", "ptu"
  )

  compartmentData <- list(
    cd8_tn_pb    = list(analyte = "naive CD8+ CAR-T cells", units = "cells/L", specimen = "whole blood", verified = FALSE),
    cd8_cm_pb    = list(analyte = "central memory CD8+ CAR-T cells", units = "cells/L", specimen = "whole blood", verified = FALSE),
    cd8_em_pb    = list(analyte = "effector memory CD8+ CAR-T cells", units = "cells/L", specimen = "whole blood", verified = FALSE),
    cd8_ef_pb    = list(analyte = "effector CD8+ CAR-T cells", units = "cells/L", specimen = "whole blood", verified = FALSE),
    cd4_tn_pb    = list(analyte = "naive CD4+ CAR-T cells", units = "cells/L", specimen = "whole blood", verified = FALSE),
    cd4_cm_pb    = list(analyte = "central memory CD4+ CAR-T cells", units = "cells/L", specimen = "whole blood", verified = FALSE),
    cd4_em_pb    = list(analyte = "effector memory CD4+ CAR-T cells", units = "cells/L", specimen = "whole blood", verified = FALSE),
    cd4_ef_pb    = list(analyte = "effector CD4+ CAR-T cells", units = "cells/L", specimen = "whole blood", verified = FALSE),
    cd8_tn_bm    = list(analyte = "naive CD8+ CAR-T cells", units = "cells/L", specimen = "tissue", verified = FALSE),
    cd8_cm_bm    = list(analyte = "central memory CD8+ CAR-T cells", units = "cells/L", specimen = "tissue", verified = FALSE),
    cd8_em_bm    = list(analyte = "effector memory CD8+ CAR-T cells", units = "cells/L", specimen = "tissue", verified = FALSE),
    cd8_ef_bm    = list(analyte = "effector CD8+ CAR-T cells", units = "cells/L", specimen = "tissue", verified = FALSE),
    cd4_tn_bm    = list(analyte = "naive CD4+ CAR-T cells", units = "cells/L", specimen = "tissue", verified = FALSE),
    cd4_cm_bm    = list(analyte = "central memory CD4+ CAR-T cells", units = "cells/L", specimen = "tissue", verified = FALSE),
    cd4_em_bm    = list(analyte = "effector memory CD4+ CAR-T cells", units = "cells/L", specimen = "tissue", verified = FALSE),
    cd4_ef_bm    = list(analyte = "effector CD4+ CAR-T cells", units = "cells/L", specimen = "tissue", verified = FALSE),
    cplx_cd8_pb  = list(analyte = "CAR-CD19 complexes on CD8+ CAR-T cells", units = "number/L", specimen = "whole blood", verified = FALSE),
    cplx_cd4_pb  = list(analyte = "CAR-CD19 complexes on CD4+ CAR-T cells", units = "number/L", specimen = "whole blood", verified = FALSE),
    cplx_cd8_bm  = list(analyte = "CAR-CD19 complexes on CD8+ CAR-T cells", units = "number/L", specimen = "tissue", verified = FALSE),
    cplx_cd4_bm  = list(analyte = "CAR-CD19 complexes on CD4+ CAR-T cells", units = "number/L", specimen = "tissue", verified = FALSE),
    cd19b_pb     = list(analyte = "autoreactive CD19+ B cells", units = "cells/L", specimen = "whole blood", verified = FALSE),
    cd19b_rev_pb = list(analyte = "reconstituted naive CD19+ B cells", units = "cells/L", specimen = "whole blood", verified = FALSE),
    cd19b_bm     = list(analyte = "CD19+ B cells", units = "cells/L", specimen = "tissue", verified = FALSE),
    ldc          = list(analyte = "fludarabine plus cyclophosphamide K-PD signal", units = "mg", specimen = "plasma", verified = FALSE),
    cd4t_pre     = list(analyte = "host CD4+ T cell precursors", units = "cells/uL", specimen = "tissue", verified = FALSE),
    cd4t         = list(analyte = "host CD4+ T cells", units = "cells/uL", specimen = "whole blood", verified = FALSE),
    cd8t_pre     = list(analyte = "host CD8+ T cell precursors", units = "cells/uL", specimen = "tissue", verified = FALSE),
    cd8t         = list(analyte = "host CD8+ T cells", units = "cells/uL", specimen = "whole blood", verified = FALSE),
    nt_prol      = list(analyte = "proliferating neutrophil precursors", units = "cells/uL", specimen = "tissue", verified = FALSE),
    nt_tr1       = list(analyte = "maturing neutrophils, transit 1", units = "cells/uL", specimen = "tissue", verified = FALSE),
    nt_tr2       = list(analyte = "maturing neutrophils, transit 2", units = "cells/uL", specimen = "tissue", verified = FALSE),
    nt_tr3       = list(analyte = "maturing neutrophils, transit 3", units = "cells/uL", specimen = "tissue", verified = FALSE),
    nt           = list(analyte = "circulating neutrophils", units = "cells/uL", specimen = "whole blood", verified = FALSE),
    ifna         = list(analyte = "interferon-alpha", units = "pg/L", specimen = "serum", verified = FALSE),
    aab          = list(analyte = "anti-double-stranded-DNA autoantibodies", units = "IU/L", specimen = "serum", verified = FALSE),
    c3           = list(analyte = "complement C3 protein", units = "mg/L", specimen = "serum", verified = FALSE),
    ptu          = list(analyte = "urinary protein", units = "mg/g creatinine", specimen = "urine", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      notes              = paste(
        "Scales the blood volume (Vb = 5 * WT/70 L) and bone marrow volume",
        "(Vbm = 1.6 * WT/70 L), the total administered CAR-T cell number",
        "(dose_cart * WT), the C3 synthesis rate (1.5 mg/kg/h * 24 * WT) and",
        "the daily urine volume (Vu = 1 cc/kg/h = WT * 24 / 100 dL/day).",
        "Park 2025 does not report individual body weights, and the clinical",
        "source of the cellular-kinetic data (Mackensen 2022, Nat Med",
        "28:2124-2132) is not open access, so no per-patient weights could be",
        "verified on disk. Every allometric term is normalized to 70 kg, which",
        "is the reference used throughout the validation vignette."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 5L,
    n_studies      = 1L,
    age_range      = "18-24 years",
    weight_range   = "not reported by Park 2025; 70 kg reference used",
    sex_female_pct = 80,
    disease_state  = paste(
      "Severe, treatment-refractory systemic lupus erythematosus. Model",
      "qualification used an external cohort of 13 additional patients",
      "(7 SLE, 6 other B cell-mediated autoimmune diseases) followed for",
      "up to 2 years."
    ),
    dose_range     = paste(
      "1 x 10^6 FMC63-based anti-CD19 CAR-positive T cells/kg as a single",
      "intravenous infusion. Lymphodepletion:",
      "fludarabine 25 mg/m2/day on days -5 to -3 plus cyclophosphamide",
      "1000 mg/m2 on day -3. Model-based dose-ranging covered 1, 5, 10, 50",
      "and 100% of the clinical dose."
    ),
    regions        = "Germany (Erlangen); Mackensen 2022 / Muller 2024 case series.",
    notes          = paste(
      "Four female (18-24 years) and one male (23 years) SLE patient.",
      "Infusion product composition 75.9% CD4+ and 23.0% CD8+ CAR-T cells,",
      "with 7.29% naive, 22.1% central memory, 60.9% effector memory and",
      "9.77% effector immunophenotypes (Methods 2.1 and Appendix S1).",
      "Cellular kinetics and biomarker profiles were digitized by the authors",
      "from Mackensen 2022 (Nat Med 28:2124-2132) and Muller 2024."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # CAR-T cellular kinetics and CD19+ B cell depletion (Park 2025 Table 1,
    # "Cellular kinetics of CAR-T cells and CD19+ B cell depletion" block).
    # Parameter symbols in the trailing comments are the control-stream names
    # from Appendix S2 (the Monolix MLXTRAN source), which is the authority for
    # every equation encoded below. Appendix S2 is distributed in the PubMed
    # Central open-access package for PMC11815715 as CTS-18-e70146-s004.docx;
    # Appendix S1 (the typeset equations 1-53 and the SLE-DAI regression output)
    # is CTS-18-e70146-s003.docx. Where the two disagree the control stream is
    # used, and each such case is listed in the vignette Errata.
    # ---------------------------------------------------------------------
    lkexp_max <- log(1.8)
    label("Maximum first-order CAR-T expansion rate constant, CD8+ subset (KexpMax, 1/day)")   # Park 2025 Table 1 (Kexp1_max)
    lkkill_max <- log(2.76)
    label("Maximum first-order CD19+ B cell killing rate constant (KkillMax, 1/day)")          # Park 2025 Table 1 (Kkill_max)
    lec50_exp <- fixed(log(0.000019))
    label("CAR-target complexes per effector cell giving 50% of maximum expansion (EC50exp, number/cell)")  # Park 2025 Table 1; estimated first, then held constant for identifiability (Appendix S1)
    lrm1 <- log(0.15)
    label("Differentiation rate constant, naive to central memory (Rm1, 1/day)")               # Park 2025 Table 1 (rm1)
    lrm2 <- log(0.11)
    label("Differentiation rate constant, central memory to effector memory (Rm2, 1/day)")     # Park 2025 Table 1 (rm2)
    lrm3 <- log(0.72)
    label("Differentiation rate constant, effector memory to effector (Rm3, 1/day)")           # Park 2025 Table 1 (rm3)
    lkel_e <- log(47.38)
    label("Elimination rate constant of effector CAR-T cells (Kel.e, 1/day)")                  # Park 2025 Table 1 (Kel_e)
    lkel_m <- log(0.69)
    label("Elimination rate constant of naive and memory CAR-T cells (Kel.m, 1/day)")          # Park 2025 Table 1 (Kel_m)
    lk12 <- log(1.41)
    label("Distribution rate constant, blood to bone marrow (K12, 1/day)")                     # Park 2025 Table 1 (K12)
    lk21 <- log(0.0017)
    label("Redistribution rate constant, bone marrow to blood (K21, 1/day)")                   # Park 2025 Table 1 (K21)
    lkc50_cart <- log(7.71)
    label("CAR-target complexes per CD19+ B cell giving 50% of maximum killing (KC50, number/cell)")  # Park 2025 Table 1 (IC50)
    lcart_rev_min <- log(2.38)
    label("Total CAR-T threshold below which naive B cell reconstitution starts (CARTrev.min, cells/uL)")  # Park 2025 Table 1 (TtoR); Table 1 prints cells/L but Appendix S2 compares against the cells/uL total
    lk0_rev_cd19b <- log(1.58e6)
    label("Zero-order generation rate constant of naive CD19+ B cells (K0.rev, cells/L/day)")  # Park 2025 Table 1 (K0_rev_CD19B)
    lcd19b_pb_max <- log(449e6)
    label("Maximum CD19+ B cell level in peripheral blood (CD19B.PB.max, cells/L)")            # Park 2025 Table 1 (CD19B_PB_MAX)

    # System constants held at literature values (Park 2025 Table 1 "Note"
    # column gives the source reference for each).
    lrbase_cd19b <- fixed(log(1.43e8))
    label("Baseline reactive CD19+ B cells in peripheral blood (cells/L)")                     # Park 2025 Table 1 reference 21 (regressor CD19B_PB_int * 1e6)
    cd19b_bm_0 <- fixed(1.50e9)
    label("Baseline CD19+ B cells in bone marrow (cells/L)")                                   # Park 2025 Table 1 reference 22 (CD19B_T_0)
    kg_cd19b <- fixed(0.0193)
    label("First-order proliferation rate constant of CD19+ B cells (1/day)")                  # Park 2025 Table 1 reference 23 (Kg_CD19B)
    density_car <- fixed(15000)
    label("CAR density on the CAR-T cell surface (receptors/cell)")                            # Park 2025 Table 1 references 17 and 24 (Density_CAR)
    density_cd19b <- fixed(14000)
    label("CD19 antigen density on the CD19+ B cell surface (antigens/cell)")                  # Park 2025 Table 1 reference 18 (Density_CD19B)
    kon_car <- fixed(2.1e5)
    label("Association rate constant of the FMC63 CAR for CD19 (1/M/s)")                       # Park 2025 Table 1 reference 25 (Kon_orig)
    koff_car <- fixed(0.0000681)
    label("Dissociation rate constant of the FMC63 CAR from CD19 (1/s)")                       # Park 2025 Table 1 reference 25 (Koff_orig)
    vb_ref <- fixed(5)
    label("Peripheral blood volume at the 70 kg reference weight (L)")                         # Park 2025 Table 1 reference 22 (Vb)
    vbm_ref <- fixed(1.6)
    label("Bone marrow volume at the 70 kg reference weight (L)")                              # Park 2025 Table 1 reference 22 (Vt)
    ratio_kexp <- fixed(1.12)
    label("Ratio of the CD8+ to CD4+ maximum expansion rate constants (unitless)")             # Appendix S2 (Kexp2_max = Kexp1_max/1.12), from Salem 2023 as cited in Methods 2.4

    # ---------------------------------------------------------------------
    # Lymphodepletion chemotherapy (LDC) effect on host immune cells
    # (Park 2025 Table 1, "LDC effect on host immune cells" block).
    # ---------------------------------------------------------------------
    lkel_ldc <- log(1.39)
    label("First-order elimination rate constant of the K-PD lymphodepletion signal (1/day)")  # Park 2025 Table 1 (Kel)
    kprol_cd4pre <- fixed(0.0007)
    label("Differentiation rate constant, precursor to CD4+ T cells (1/day)")                  # Park 2025 Table 1 reference 26 (Pc_CD4_Kp)
    lk0_cd4pre <- log(210.18)
    label("Zero-order proliferation rate constant of the CD4+ T cell precursor (cells/uL/day)")  # Park 2025 Table 1 (Pc_CD4_K0)
    lkout_cd4 <- log(1.02)
    label("First-order disappearance rate constant of host CD4+ T cells (1/day)")              # Park 2025 Table 1 (CD4_Kout)
    lemax_cd4 <- log(0.85)
    label("Maximum LDC effect on host CD4+ T cells (unitless)")                                # Park 2025 Table 1 (Emax_CD4)
    lec50_cd4 <- fixed(log(0.15))
    label("LDC concentration producing 50% of the maximum CD4+ T cell effect (mg/L)")          # Park 2025 Table 1; estimated first, then held constant (Results 3.1)
    kprol_cd8pre <- fixed(0.00045)
    label("Differentiation rate constant, precursor to CD8+ T cells (1/day)")                  # Park 2025 Table 1 reference 26 (Pc_CD8_Kp)
    lk0_cd8pre <- log(525.59)
    label("Zero-order proliferation rate constant of the CD8+ T cell precursor (cells/uL/day)")  # Park 2025 Table 1 (Pc_CD8_K0)
    lkout_cd8 <- log(1.49)
    label("First-order disappearance rate constant of host CD8+ T cells (1/day)")              # Park 2025 Table 1 (CD8_Kout)
    lemax_cd8 <- log(0.96)
    label("Maximum LDC effect on host CD8+ T cells (unitless)")                                # Park 2025 Table 1 (Emax_CD8)
    lec50_cd8 <- fixed(log(0.067))
    label("LDC concentration producing 50% of the maximum CD8+ T cell effect (mg/L)")          # Park 2025 Table 1; estimated first, then held constant (Results 3.1)
    lrbase_nt <- log(5617.22)
    label("Steady-state circulating neutrophil count (cells/uL)")                              # Park 2025 Table 1 (Circ0_NT)
    lmtt_nt <- fixed(log(4.35))
    label("Mean transit time of the neutrophil proliferation chain (day)")                     # Park 2025 Table 1 reference 27 (MTT_NT)
    lgamma_nt <- log(0.14)
    label("Feedback exponent on the circulating neutrophil count (unitless)")                  # Park 2025 Table 1 (gam_NT)
    lemax_nt <- log(27.32)
    label("Maximum LDC effect on host neutrophils (unitless)")                                 # Park 2025 Table 1 (Emax_NT)
    lec50_nt <- fixed(log(1203.42))
    label("LDC concentration producing 50% of the maximum neutrophil effect (mg/L)")           # Park 2025 Table 1; estimated first, then held constant (Results 3.1)

    # ---------------------------------------------------------------------
    # SLE disease biomarkers (Park 2025 Table 1, "Changes in biomarkers").
    # ---------------------------------------------------------------------
    lksyn_aab <- log(15875.77)
    label("Release rate of anti-dsDNA autoantibodies from reactive B cells (IU/cell/day)")     # Park 2025 Table 1 (Ksyn_AAB)
    lkdeg_ifna <- fixed(log(23.04))
    label("Degradation rate constant of interferon-alpha in peripheral blood (1/day)")         # Park 2025 Table 1 reference 28 (Kdeg_IFN)
    n_pdc <- fixed(1146788.991)
    label("Circulating plasmacytoid dendritic cell count (cells/L)")                           # Appendix S2 control stream (N_pDC); equals 1.147 cells/uL
    lgamma_c3 <- log(0.1)
    label("Shaping exponent for the anti-dsDNA drive on C3 catabolism (unitless)")             # Park 2025 Table 1 (C3gam)
    ksyn_c3 <- fixed(1.5)
    label("Zero-order synthesis rate of complement C3 (mg/kg/h)")                              # Park 2025 Table 1 reference 29 (C3_BP_Ksyn = 1.5 * 24 * WT)
    lgamma_ptu <- log(1.09)
    label("Shaping exponent for the C3 drive on proteinuria (unitless)")                       # Park 2025 Table 1 (PTUgam)
    ksyn_ptu <- fixed(1000)
    label("Zero-order urinary protein excretion rate (mg/day)")                                # Appendix S2 control stream (PTU_Ksyn = 1000/1); Park 2025 Table 1 prints 100, see vignette Errata
    lgamma_ifna_cd4 <- log(0.09)
    label("Shaping exponent for interferon-alpha activation of anti-dsDNA release (unitless)") # Park 2025 Table 1 (gam_IFNA_4)
    lgamma_ifna_cd8 <- log(0.8)
    label("Shaping exponent for interferon-alpha activation of proteinuria (unitless)")        # Park 2025 Table 1 (gam_IFNA_8)

    # ---------------------------------------------------------------------
    # Individual baselines. In the source these are per-subject regressors in
    # the Monolix dataset (Appendix S2 input block), not estimated parameters.
    # The defaults below are the medians of the five modelled patients, read
    # off Figure S1 panel E; see the vignette source-trace table.
    # ---------------------------------------------------------------------
    lrbase_aab <- fixed(log(470000))
    label("Baseline anti-dsDNA autoantibody concentration (IU/L)")                             # digitized from Park 2025 Figure S1E, median of 5 patients (regressor dsDNA0)
    lrbase_ifna <- fixed(log(21))
    label("Baseline interferon-alpha concentration (pg/mL)")                                   # digitized from Park 2025 Figure S1E, median of 5 patients (regressor INFA0)
    lrbase_c3 <- fixed(log(555))
    label("Baseline complement C3 concentration (mg/L)")                                       # digitized from Park 2025 Figure S1E, median of 5 patients (regressor C30)
    lrbase_ptu <- fixed(log(3060))
    label("Baseline urinary protein excretion (mg/g creatinine)")                              # digitized from Park 2025 Figure S1E, median of 5 patients (regressor PTU0)

    # ---------------------------------------------------------------------
    # Dosing inputs. The CAR-T product is delivered through the compartment
    # initial conditions (Appendix S2 gives the infusion as regressor-derived
    # initial amounts and a mock zero bolus), so the dose level and the product
    # composition are model parameters rather than event-table records.
    # ---------------------------------------------------------------------
    dose_cart <- fixed(1e6)
    label("Administered CAR-positive T cell dose (cells/kg)")                                  # Park 2025 Methods 2.1
    pct_cd4 <- fixed(75.9)
    label("CD4+ fraction of the infused CAR-T product (percent)")                              # Park 2025 Methods 2.1
    pct_cd8 <- fixed(23.0)
    label("CD8+ fraction of the infused CAR-T product (percent)")                              # Park 2025 Methods 2.1
    pct_tn <- fixed(7.29)
    label("Naive fraction of the infused CAR-T product (percent)")                             # Park 2025 Methods 2.1
    pct_tcm <- fixed(22.1)
    label("Central memory fraction of the infused CAR-T product (percent)")                    # Park 2025 Methods 2.1
    pct_tem <- fixed(60.9)
    label("Effector memory fraction of the infused CAR-T product (percent)")                   # Park 2025 Methods 2.1
    pct_tef <- fixed(9.77)
    label("Effector fraction of the infused CAR-T product (percent)")                          # Park 2025 Methods 2.1
    dose_ldc <- fixed(1935)
    label("Total lymphodepletion chemotherapy amount entering the K-PD compartment (mg)")      # not published; reconstructed from Park 2025 Discussion (fludarabine 25 mg/m2/day x 3 days + cyclophosphamide 1000 mg/m2) at a 1.8 m2 reference body surface area, see vignette Errata

    # ---------------------------------------------------------------------
    # SLE-DAI translation (Park 2025 Section 2.7 and Figure 4A). Coefficients
    # are the literature-based (SLEDAI.Ref) univariate linear regressions in
    # the Appendix S1 R regression output; the paper uses these, not the
    # in-study fits, to a priori predict clinical scores (Results 3.3).
    # ---------------------------------------------------------------------
    sledai_c3_int <- fixed(16.99)
    label("SLE-DAI intercept for the complement C3 regression (score)")                        # Park 2025 Appendix S1 regression output, SLEDAI.Ref column
    sledai_c3_slope <- fixed(-0.01)
    label("SLE-DAI slope on complement C3 (score per mg/L)")                                   # Park 2025 Appendix S1 regression output, SLEDAI.Ref column
    sledai_ifna_int <- fixed(4.21)
    label("SLE-DAI intercept for the interferon-alpha regression (score)")                     # Park 2025 Appendix S1 regression output, SLEDAI.Ref column
    sledai_ifna_slope <- fixed(0.08)
    label("SLE-DAI slope on interferon-alpha (score per pg/mL)")                               # Park 2025 Appendix S1 regression output, SLEDAI.Ref column
    sledai_aab_int <- fixed(6.96)
    label("SLE-DAI intercept for the anti-dsDNA regression (score)")                           # Park 2025 Appendix S1 regression output, SLEDAI.Ref column
    sledai_aab_slope <- fixed(0.05)
    label("SLE-DAI slope on anti-dsDNA autoantibodies (score per IU/mL)")                      # Park 2025 Appendix S1 regression output, SLEDAI.Ref column; see vignette Errata on the Figure 4A3 slope

    # ---------------------------------------------------------------------
    # Inter-individual variability. Park 2025 Table 1 reports the Monolix
    # random effect as omega (the standard deviation of the normally
    # distributed eta on the log-transformed parameter), so the variances
    # entered here are omega^2. Appendix S1 confirms a log-normal IIV
    # distribution fitted by SAEM in Monolix 2023R1.
    # ---------------------------------------------------------------------
    etalkexp_max ~ 0.0025          # Park 2025 Table 1, omega = 0.05
    etalkkill_max ~ 0.007569       # Park 2025 Table 1, omega = 0.087
    etalec50_exp ~ fixed(0.04)     # Park 2025 Table 1, omega = 0.2 held constant
    etalcart_rev_min ~ 0.6561      # Park 2025 Table 1, omega = 0.81
    etalk0_rev_cd19b ~ 1.0816      # Park 2025 Table 1, omega = 1.04
    etalcd19b_pb_max ~ 0.6561      # Park 2025 Table 1, omega = 0.81
    etalk0_cd4pre ~ 0.4624         # Park 2025 Table 1, omega = 0.68
    etalkout_cd4 ~ 0.0225          # Park 2025 Table 1, omega = 0.15
    etalkout_cd8 ~ 0.25            # Park 2025 Table 1, omega = 0.5
    etalrbase_nt ~ 0.0841          # Park 2025 Table 1, omega = 0.29
    etalmtt_nt ~ fixed(0.04)       # Park 2025 Table 1, omega = 0.2 held constant
    etalemax_nt ~ 0.0676           # Park 2025 Table 1, omega = 0.26
    etalksyn_aab ~ 0.16            # Park 2025 Table 1, omega = 0.4
    etalkdeg_ifna ~ fixed(0.04)    # Park 2025 Table 1, omega = 0.2 held constant
  })

  model({
    # -------------------------------------------------------------------
    # Physical constants and unit conversions (Appendix S2 control stream).
    # -------------------------------------------------------------------
    avogadro <- 6.023e23                    # molecules/mol, as written in Appendix S2
    kon <- kon_car * 60 * 60 * 24 / avogadro  # 1/M/s -> 1/(number/L)/day
    koff <- koff_car * 60 * 60 * 24           # 1/s -> 1/day

    # Allometric volumes (Park 2025 Table 1)
    vb <- vb_ref * (WT / 70)                # peripheral blood volume, L
    vbm <- vbm_ref * (WT / 70)              # bone marrow volume, L
    vu <- WT * 24 / 100                     # daily urine volume, dL/day (1 cc/kg/h)

    # -------------------------------------------------------------------
    # Individual parameters
    # -------------------------------------------------------------------
    kexp_max_cd8 <- exp(lkexp_max + etalkexp_max)
    kexp_max_cd4 <- kexp_max_cd8 / ratio_kexp
    kkill_max <- exp(lkkill_max + etalkkill_max)
    ec50_exp <- exp(lec50_exp + etalec50_exp)
    rm1 <- exp(lrm1)
    rm2 <- exp(lrm2)
    rm3 <- exp(lrm3)
    kel_e <- exp(lkel_e)
    kel_m <- exp(lkel_m)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kc50_cart <- exp(lkc50_cart)
    cart_rev_min <- exp(lcart_rev_min + etalcart_rev_min)
    k0_rev_cd19b <- exp(lk0_rev_cd19b + etalk0_rev_cd19b)
    cd19b_pb_max <- exp(lcd19b_pb_max + etalcd19b_pb_max)
    rbase_cd19b <- exp(lrbase_cd19b)

    kel_ldc <- exp(lkel_ldc)
    k0_cd4pre <- exp(lk0_cd4pre + etalk0_cd4pre)
    kout_cd4 <- exp(lkout_cd4 + etalkout_cd4)
    emax_cd4 <- exp(lemax_cd4)
    ec50_cd4 <- exp(lec50_cd4)
    k0_cd8pre <- exp(lk0_cd8pre)
    kout_cd8 <- exp(lkout_cd8 + etalkout_cd8)
    emax_cd8 <- exp(lemax_cd8)
    ec50_cd8 <- exp(lec50_cd8)
    rbase_nt <- exp(lrbase_nt + etalrbase_nt)
    mtt_nt <- exp(lmtt_nt + etalmtt_nt)
    gamma_nt <- exp(lgamma_nt)
    emax_nt <- exp(lemax_nt + etalemax_nt)
    ec50_nt <- exp(lec50_nt)

    ksyn_aab <- exp(lksyn_aab + etalksyn_aab)
    kdeg_ifna <- exp(lkdeg_ifna + etalkdeg_ifna)
    gamma_c3 <- exp(lgamma_c3)
    gamma_ptu <- exp(lgamma_ptu)
    gamma_ifna_cd4 <- exp(lgamma_ifna_cd4)
    gamma_ifna_cd8 <- exp(lgamma_ifna_cd8)

    rbase_aab <- exp(lrbase_aab)
    rbase_ifna <- exp(lrbase_ifna) * 1000   # pg/mL -> pg/L, the state units
    rbase_c3 <- exp(lrbase_c3)
    rbase_ptu <- exp(lrbase_ptu)

    # Host T cell baselines. Appendix S1 equations 29 and 30 define the
    # precursor pool at K0/Kprol and the circulating pool at Kprol * precursor
    # / Kout, which simplifies to K0/Kout.
    cd4t_pre_base <- k0_cd4pre / kprol_cd4pre
    cd4t_base <- k0_cd4pre / kout_cd4
    cd8t_pre_base <- k0_cd8pre / kprol_cd8pre
    cd8t_base <- k0_cd8pre / kout_cd8

    # -------------------------------------------------------------------
    # Initial conditions
    # -------------------------------------------------------------------
    # CAR-T product distributed across subset x immunophenotype. The 10000
    # divisor converts the two percentage factors to fractions (Appendix S2).
    cd8_tn_pb(0) <- dose_cart * WT * pct_cd8 * pct_tn / (vb * 10000)
    cd8_cm_pb(0) <- dose_cart * WT * pct_cd8 * pct_tcm / (vb * 10000)
    cd8_em_pb(0) <- dose_cart * WT * pct_cd8 * pct_tem / (vb * 10000)
    cd8_ef_pb(0) <- dose_cart * WT * pct_cd8 * pct_tef / (vb * 10000)
    cd4_tn_pb(0) <- dose_cart * WT * pct_cd4 * pct_tn / (vb * 10000)
    cd4_cm_pb(0) <- dose_cart * WT * pct_cd4 * pct_tcm / (vb * 10000)
    cd4_em_pb(0) <- dose_cart * WT * pct_cd4 * pct_tem / (vb * 10000)
    cd4_ef_pb(0) <- dose_cart * WT * pct_cd4 * pct_tef / (vb * 10000)

    # Bone marrow CAR-T pools and all CAR-target complexes start at the 1e-12
    # numerical floor used in Appendix S2 (Appendix S1 states an initial
    # condition of zero; the floor keeps the complexes-per-cell ratios finite).
    cd8_tn_bm(0) <- 1e-12
    cd8_cm_bm(0) <- 1e-12
    cd8_em_bm(0) <- 1e-12
    cd8_ef_bm(0) <- 1e-12
    cd4_tn_bm(0) <- 1e-12
    cd4_cm_bm(0) <- 1e-12
    cd4_em_bm(0) <- 1e-12
    cd4_ef_bm(0) <- 1e-12
    cplx_cd8_pb(0) <- 1e-12
    cplx_cd4_pb(0) <- 1e-12
    cplx_cd8_bm(0) <- 1e-12
    cplx_cd4_bm(0) <- 1e-12

    cd19b_pb(0) <- rbase_cd19b
    cd19b_rev_pb(0) <- 0
    cd19b_bm(0) <- cd19b_bm_0

    ldc(0) <- dose_ldc
    cd4t_pre(0) <- cd4t_pre_base
    cd4t(0) <- cd4t_base
    cd8t_pre(0) <- cd8t_pre_base
    cd8t(0) <- cd8t_base

    nt_prol(0) <- rbase_nt
    nt_tr1(0) <- rbase_nt
    nt_tr2(0) <- rbase_nt
    nt_tr3(0) <- rbase_nt
    nt(0) <- rbase_nt

    ifna(0) <- rbase_ifna
    aab(0) <- rbase_aab
    c3(0) <- rbase_c3
    ptu(0) <- rbase_ptu

    # -------------------------------------------------------------------
    # Numerical safeguards.
    #
    # Every state below is a non-negative cell count, receptor count or
    # concentration, so none of these guards changes the published solution in
    # exact arithmetic. They are needed because CAR-T mediated eradication
    # drives the CD19+ B cell pools down by more than thirty orders of
    # magnitude, far below any absolute solver tolerance, where round-off can
    # make a state or a free-species balance mildly negative. A negative value
    # entering a complexes-per-cell ratio, a Michaelis-Menten denominator or a
    # fractional power turns the right-hand side into NaN and aborts the
    # integration. The floor is 1e-12 cells/L, i.e. 1e-18 cells per microlitre.
    # -------------------------------------------------------------------
    tiny <- 1e-12
    cd19b_pb_pos <- max(cd19b_pb, tiny)
    cd19b_bm_pos <- max(cd19b_bm, tiny)
    nt_pos <- max(nt, tiny)
    ifna_pos <- max(ifna, tiny)
    aab_pos <- max(aab, tiny)
    c3_pos <- max(c3, tiny)

    # -------------------------------------------------------------------
    # Peripheral blood: target engagement and expansion drivers
    # (Appendix S1 equations 9-15)
    # -------------------------------------------------------------------
    cart_cd8_pb <- cd8_tn_pb + cd8_cm_pb + cd8_em_pb + cd8_ef_pb
    cart_cd4_pb <- cd4_tn_pb + cd4_cm_pb + cd4_em_pb + cd4_ef_pb
    cplx_per_cart_cd8_pb <- max(cplx_cd8_pb, 0) / max(cart_cd8_pb, tiny)
    cplx_per_cart_cd4_pb <- max(cplx_cd4_pb, 0) / max(cart_cd4_pb, tiny)
    kexp_cd8_pb <- kexp_max_cd8 * cplx_per_cart_cd8_pb / (ec50_exp + cplx_per_cart_cd8_pb)
    kexp_cd4_pb <- kexp_max_cd4 * cplx_per_cart_cd4_pb / (ec50_exp + cplx_per_cart_cd4_pb)
    car_free_cd8_pb <- max(cart_cd8_pb * density_car - cplx_cd8_pb, 0)
    car_free_cd4_pb <- max(cart_cd4_pb * density_car - cplx_cd4_pb, 0)
    ag_free_pb <- max(cd19b_pb * density_cd19b - (cplx_cd8_pb + cplx_cd4_pb), 0)

    # CD8+ CAR-T immunophenotypes in peripheral blood (Appendix S1 eq 1-4)
    d/dt(cd8_tn_pb) <- (-k12 * vb * cd8_tn_pb + k21 * vbm * cd8_tn_bm) / vb -
      kel_m * cd8_tn_pb + kexp_cd8_pb * cd8_tn_pb - rm1 * cd8_tn_pb
    d/dt(cd8_cm_pb) <- (-k12 * vb * cd8_cm_pb + k21 * vbm * cd8_cm_bm) / vb -
      kel_m * cd8_cm_pb + kexp_cd8_pb * cd8_cm_pb + rm1 * cd8_tn_pb - rm2 * cd8_cm_pb
    d/dt(cd8_em_pb) <- (-k12 * vb * cd8_em_pb + k21 * vbm * cd8_em_bm) / vb -
      kel_m * cd8_em_pb + kexp_cd8_pb * cd8_em_pb + rm2 * cd8_cm_pb - rm3 * cd8_em_pb
    d/dt(cd8_ef_pb) <- (-k12 * vb * cd8_ef_pb + k21 * vbm * cd8_ef_bm) / vb -
      kel_e * cd8_ef_pb + rm3 * cd8_em_pb

    # CD4+ CAR-T immunophenotypes in peripheral blood (Appendix S1 eq 1-4)
    d/dt(cd4_tn_pb) <- (-k12 * vb * cd4_tn_pb + k21 * vbm * cd4_tn_bm) / vb -
      kel_m * cd4_tn_pb + kexp_cd4_pb * cd4_tn_pb - rm1 * cd4_tn_pb
    d/dt(cd4_cm_pb) <- (-k12 * vb * cd4_cm_pb + k21 * vbm * cd4_cm_bm) / vb -
      kel_m * cd4_cm_pb + kexp_cd4_pb * cd4_cm_pb + rm1 * cd4_tn_pb - rm2 * cd4_cm_pb
    d/dt(cd4_em_pb) <- (-k12 * vb * cd4_em_pb + k21 * vbm * cd4_em_bm) / vb -
      kel_m * cd4_em_pb + kexp_cd4_pb * cd4_em_pb + rm2 * cd4_cm_pb - rm3 * cd4_em_pb
    d/dt(cd4_ef_pb) <- (-k12 * vb * cd4_ef_pb + k21 * vbm * cd4_ef_bm) / vb -
      kel_e * cd4_ef_pb + rm3 * cd4_em_pb

    # CAR-target complexes in peripheral blood (Appendix S1 eq 11)
    d/dt(cplx_cd8_pb) <- kon * car_free_cd8_pb * ag_free_pb - koff * cplx_cd8_pb
    d/dt(cplx_cd4_pb) <- kon * car_free_cd4_pb * ag_free_pb - koff * cplx_cd4_pb

    # Reactive CD19+ B cells in peripheral blood (Appendix S1 eq 13, 15, 16)
    cplx_per_b_pb <- max(cplx_cd8_pb + cplx_cd4_pb, 0) / cd19b_pb_pos
    kkill_pb <- kkill_max * cplx_per_b_pb / (kc50_cart + cplx_per_b_pb)
    d/dt(cd19b_pb) <- kg_cd19b * cd19b_pb - kkill_pb * cd19b_pb

    # Naive (non-autoreactive) CD19+ B cell reconstitution. Appendix S2 gates
    # the zero-order influx on the total CAR-T count expressed in cells/uL.
    total_cart <- (cart_cd8_pb + cart_cd4_pb) * 1e-6
    if (total_cart >= cart_rev_min) {
      k_rev <- 0
    } else {
      k_rev <- k0_rev_cd19b
    }
    d/dt(cd19b_rev_pb) <- k_rev * (1 - (cd19b_pb + cd19b_rev_pb) / cd19b_pb_max)

    # -------------------------------------------------------------------
    # Bone marrow: target engagement and expansion drivers
    # (Appendix S1 equations 19-26)
    # -------------------------------------------------------------------
    cart_cd8_bm <- cd8_tn_bm + cd8_cm_bm + cd8_em_bm + cd8_ef_bm
    cart_cd4_bm <- cd4_tn_bm + cd4_cm_bm + cd4_em_bm + cd4_ef_bm
    cplx_per_cart_cd8_bm <- max(cplx_cd8_bm, 0) / max(cart_cd8_bm, tiny)
    cplx_per_cart_cd4_bm <- max(cplx_cd4_bm, 0) / max(cart_cd4_bm, tiny)
    kexp_cd8_bm <- kexp_max_cd8 * cplx_per_cart_cd8_bm / (ec50_exp + cplx_per_cart_cd8_bm)
    kexp_cd4_bm <- kexp_max_cd4 * cplx_per_cart_cd4_bm / (ec50_exp + cplx_per_cart_cd4_bm)
    car_free_cd8_bm <- max(cart_cd8_bm * density_car - cplx_cd8_bm, 0)
    car_free_cd4_bm <- max(cart_cd4_bm * density_car - cplx_cd4_bm, 0)
    ag_free_bm <- max(cd19b_bm * density_cd19b - (cplx_cd8_bm + cplx_cd4_bm), 0)

    d/dt(cd8_tn_bm) <- (k12 * vb * cd8_tn_pb - k21 * vbm * cd8_tn_bm) / vbm +
      kexp_cd8_bm * cd8_tn_bm - rm1 * cd8_tn_bm
    d/dt(cd8_cm_bm) <- (k12 * vb * cd8_cm_pb - k21 * vbm * cd8_cm_bm) / vbm +
      kexp_cd8_bm * cd8_cm_bm + rm1 * cd8_tn_bm - rm2 * cd8_cm_bm
    d/dt(cd8_em_bm) <- (k12 * vb * cd8_em_pb - k21 * vbm * cd8_em_bm) / vbm +
      kexp_cd8_bm * cd8_em_bm + rm2 * cd8_cm_bm - rm3 * cd8_em_bm
    d/dt(cd8_ef_bm) <- (k12 * vb * cd8_ef_pb - k21 * vbm * cd8_ef_bm) / vbm +
      rm3 * cd8_em_bm

    d/dt(cd4_tn_bm) <- (k12 * vb * cd4_tn_pb - k21 * vbm * cd4_tn_bm) / vbm +
      kexp_cd4_bm * cd4_tn_bm - rm1 * cd4_tn_bm
    d/dt(cd4_cm_bm) <- (k12 * vb * cd4_cm_pb - k21 * vbm * cd4_cm_bm) / vbm +
      kexp_cd4_bm * cd4_cm_bm + rm1 * cd4_tn_bm - rm2 * cd4_cm_bm
    d/dt(cd4_em_bm) <- (k12 * vb * cd4_em_pb - k21 * vbm * cd4_em_bm) / vbm +
      kexp_cd4_bm * cd4_em_bm + rm2 * cd4_cm_bm - rm3 * cd4_em_bm
    d/dt(cd4_ef_bm) <- (k12 * vb * cd4_ef_pb - k21 * vbm * cd4_ef_bm) / vbm +
      rm3 * cd4_em_bm

    d/dt(cplx_cd8_bm) <- kon * car_free_cd8_bm * ag_free_bm - koff * cplx_cd8_bm
    d/dt(cplx_cd4_bm) <- kon * car_free_cd4_bm * ag_free_bm - koff * cplx_cd4_bm

    cplx_per_b_bm <- max(cplx_cd8_bm + cplx_cd4_bm, 0) / cd19b_bm_pos
    kkill_bm <- kkill_max * cplx_per_b_bm / (kc50_cart + cplx_per_b_bm)
    d/dt(cd19b_bm) <- kg_cd19b * cd19b_bm - kkill_bm * cd19b_bm

    # -------------------------------------------------------------------
    # Lymphodepletion chemotherapy K-PD (Appendix S1 equations 27-28)
    # -------------------------------------------------------------------
    d/dt(ldc) <- -kel_ldc * ldc
    c_ldc <- ldc / vb
    eff_cd4 <- emax_cd4 * c_ldc / (ec50_cd4 + c_ldc)
    eff_cd8 <- emax_cd8 * c_ldc / (ec50_cd8 + c_ldc)
    eff_nt <- emax_nt * c_ldc / (ec50_nt + c_ldc)

    # Host CD4+ / CD8+ T cells, precursor-dependent indirect response
    # (Appendix S1 equations 31-34)
    d/dt(cd4t_pre) <- k0_cd4pre - cd4t_pre * kprol_cd4pre * (1 - eff_cd4)
    d/dt(cd4t) <- cd4t_pre * kprol_cd4pre * (1 - eff_cd4) - cd4t * kout_cd4
    d/dt(cd8t_pre) <- k0_cd8pre - cd8t_pre * kprol_cd8pre * (1 - eff_cd8)
    d/dt(cd8t) <- cd8t_pre * kprol_cd8pre * (1 - eff_cd8) - cd8t * kout_cd8

    # Neutrophils, Friberg myelosuppression chain with three transit
    # compartments (Appendix S1 equations 35-41)
    ktr_nt <- (3 + 1) / mtt_nt
    d/dt(nt_prol) <- ktr_nt * nt_prol * (1 - eff_nt) * (rbase_nt / nt_pos)^gamma_nt -
      ktr_nt * nt_prol
    d/dt(nt_tr1) <- ktr_nt * nt_prol - ktr_nt * nt_tr1
    d/dt(nt_tr2) <- ktr_nt * nt_tr1 - ktr_nt * nt_tr2
    d/dt(nt_tr3) <- ktr_nt * nt_tr2 - ktr_nt * nt_tr3
    d/dt(nt) <- ktr_nt * nt_tr3 - ktr_nt * nt

    # -------------------------------------------------------------------
    # SLE disease biomarkers (Appendix S1 equations 42-53)
    # -------------------------------------------------------------------
    # Interferon-alpha. The pDC release rate is back-calculated so that the
    # baseline is a steady state (Appendix S1 eq 49).
    ifn_rel <- rbase_ifna * kdeg_ifna / n_pdc
    d/dt(ifna) <- n_pdc * ifn_rel * (aab_pos / rbase_aab) * (nt_pos / rbase_nt) -
      ifna * kdeg_ifna

    # Anti-dsDNA autoantibodies. Release switches off once the reactive B cell
    # pool falls below 1 cell/uL, and the degradation rate constant is derived
    # from the release rate so the untreated baseline is a steady state
    # (Appendix S1 eq 45-47).
    if (cd19b_pb > 1e6) {
      kgen_aab <- ksyn_aab
    } else {
      kgen_aab <- 0
    }
    kdeg_aab <- rbase_cd19b * kgen_aab / rbase_aab
    d/dt(aab) <- cd19b_pb * kgen_aab * (cd4t / cd4t_base) *
      (ifna_pos / rbase_ifna)^gamma_ifna_cd4 - aab * kdeg_aab

    # Complement C3 (Appendix S1 eq 50-51)
    ksyn_c3_day <- ksyn_c3 * 24 * WT              # mg/kg/h -> mg/day
    kdeg_c3 <- ksyn_c3_day / (rbase_c3 * vb)
    d/dt(c3) <- (ksyn_c3_day - kdeg_c3 * c3 * vb * (aab_pos / rbase_aab)^gamma_c3) / vb

    # Proteinuria (Appendix S1 eq 52-53). Appendix S2 applies the host CD8+ T
    # cell and neutrophil signals through exp(x/x0) rather than the bare ratio
    # printed in Appendix S1 eq 44; the control stream form is used here.
    kdeg_ptu <- ksyn_ptu / (vu * rbase_ptu)
    d/dt(ptu) <- (ksyn_ptu * (rbase_c3 / c3_pos)^gamma_ptu *
      (ifna_pos / rbase_ifna)^gamma_ifna_cd8 *
      exp(cd8t / cd8t_base) * exp(nt_pos / rbase_nt) -
      ptu * vu * kdeg_ptu) / vu

    # -------------------------------------------------------------------
    # Reported outputs (Appendix S2 OUTPUT block)
    # -------------------------------------------------------------------
    totalCart <- total_cart                         # total CAR-T cells, cells/uL
    tnCart <- cd8_tn_pb + cd4_tn_pb                 # naive CAR-T cells, cells/L
    cmCart <- cd8_cm_pb + cd4_cm_pb                 # central memory CAR-T cells, cells/L
    emCart <- cd8_em_pb + cd4_em_pb                 # effector memory CAR-T cells, cells/L
    efCart <- cd8_ef_pb + cd4_ef_pb                 # effector CAR-T cells, cells/L
    bCell <- (cd19b_pb + cd19b_rev_pb) * 1e-6       # total CD19+ B cells, cells/uL
    ifnaConc <- ifna / 1000                         # interferon-alpha, pg/mL

    # SLE-DAI translation from each surrogate biomarker (Park 2025 Figure 4A)
    sledaiC3 <- sledai_c3_int + sledai_c3_slope * c3
    sledaiIfna <- sledai_ifna_int + sledai_ifna_slope * ifnaConc
    sledaiAab <- sledai_aab_int + sledai_aab_slope * (aab / 1000)
  })
}
