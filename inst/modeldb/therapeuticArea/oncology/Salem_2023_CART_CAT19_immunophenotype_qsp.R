# Multiscale mechanistic CK-PD model of low-affinity CAT19-binder anti-CD19
# CAR T cells with CD4/CD8 subsets AND four immunophenotypes (step 3, the final
# model of Salem 2023). Each subset is split into stem-cell memory (TSCM),
# central memory (TCM), effector memory (TEM) and terminally differentiated
# effector (TEFF) cells linked by a catenary differentiation chain
# (Rm1, Rm2, Rm3) in blood and bone marrow; TSCM / TCM / TEM expand on
# CAR-target engagement and are eliminated at the memory rate, TEFF does not
# expand and is eliminated at the effector rate. A high-polyclonality infusion
# product multiplies Rm1, Rm2 and Rm3 by 0.08, 0.21 and 0.19.
#
# SOURCES. Equations: Appendix S1 Eqs. 43-66, written once and applying "for
# both CD4 and CD8 subsets" (so each subset carries its own complexes, as in
# step 2 and in the Appendix S2 code). Parameter values: main-text Table 2,
# "Parameters associated with clinical CK-PD model accounting for CD4/CD8 and
# immunophenotypes subsets (CAT19)". The Appendix S2 Monolix code for Figure 3
# is an earlier iteration that disagrees with Table 2 and Appendix S1 in
# several places (a single rm for every blood transition, KC50 4.05, initial
# bone-marrow tumour 7.8E10, transgene intercept / slope 6.42 / 0.823, no
# polyclonality effect); Table 2 and the Appendix S1 equations are used here.
# Table 2 does not list the qPCR link for step 3; the step-2 CAT19 estimates
# (7.26 / 0.74) are used, as the paper states that the CK parameters of this
# step were carried over from step 2.
#
# STATE UNITS. CAR T-cell pools are cell NUMBERS; complexes and tumour cells
# are concentrations (#/L, cells/L).
#
# LONG HORIZONS. After the tumour is eradicated its states sit at solver noise,
# which the exponential growth term amplifies; beyond roughly 1.5-2 years after
# eradication the trajectories depend on the solver tolerances (see the
# vignette, "Long-horizon simulations").

Salem_2023_CART_CAT19_immunophenotype_qsp <- function() {
  description <- "QSP (multiscale mechanistic cellular-kinetic / pharmacodynamic model). Low-affinity CAT19-binder anti-CD19 CAR T cells in children with relapsed / refractory B-cell acute lymphoblastic leukaemia with CD4+ / CD8+ subsets and four immunophenotypes (stem-cell memory, central memory, effector memory, terminal effector), Salem 2023 step 3 (final model). Each subset differentiates TSCM -> TCM -> TEM -> TEFF through first-order rate constants Rm1-Rm3 in blood and bone marrow; all phenotypes distribute between blood and bone marrow (K12 / K21) and are eliminated from blood (memory rate for TSCM / TCM / TEM, effector rate for TEFF). CAR-CD19 complexes per CAR T cell drive Emax expansion of TSCM, TCM and TEM (CD4+ rate = CD8+ rate / 1.12) and complexes per tumour cell drive Emax killing of CD19+ leukaemic cells in both compartments. High polyclonality of the infusion product lowers Rm1, Rm2 and Rm3 by 92%, 79% and 81%. Outputs: percentage of each immunophenotype among blood CAR T cells, blood CAR T cells/uL and transgene copies/ug genomic DNA."
  reference <- paste(
    "Salem AM, Mugundu GM, Singh AP (2023).",
    "Development of a multiscale mechanistic modeling framework integrating",
    "differential cellular kinetics of CAR T-cell subsets and immunophenotypes",
    "in cancer patients.",
    "CPT Pharmacometrics Syst Pharmacol 12(9):1285-1296.",
    "doi:10.1002/psp4.13009.",
    "Equations in Appendix S1 and Monolix code in Appendix S2.",
    "Fit to data digitised from Ghorashian et al. 2019",
    "(Nat Med 25:1408-1414; doi:10.1038/s41591-019-0549-5) and",
    "Biasco et al. 2021 (Nat Cancer 2:629-642; doi:10.1038/s43018-021-00207-7).",
    sep = " "
  )
  vignette <- "Salem_2023_CART_subsets_immunophenotypes"

  # Per-subset, per-immunophenotype CAR T-cell pools in blood and bone marrow,
  # the per-subset complexes and the two tumour pools are paper-mechanistic
  # states with no canonical analogue.
  paper_specific_compartments <- c(
    "tscm_cd8_pb",
    "tcm_cd8_pb",
    "tem_cd8_pb",
    "teff_cd8_pb",
    "tscm_cd4_pb",
    "tcm_cd4_pb",
    "tem_cd4_pb",
    "teff_cd4_pb",
    "tscm_cd8_bm",
    "tcm_cd8_bm",
    "tem_cd8_bm",
    "teff_cd8_bm",
    "tscm_cd4_bm",
    "tcm_cd4_bm",
    "tem_cd4_bm",
    "teff_cd4_bm",
    "complex_cd8_pb",
    "complex_cd4_pb",
    "complex_cd8_bm",
    "complex_cd4_bm",
    "tumor_pb",
    "tumor_bm"
  )

  units <- list(
    time = "day",
    dosing = "cells (CAR+ T cells of each immunophenotype in the infused product, i.e. total dose x the patient's product percentage of that phenotype; each phenotype is given as two records with the same amt, one on its _cd8_pb and one on its _cd4_pb compartment, and f() applies the 38% CD8 / 54% CD4 product fractions)",
    concentration = "percent of blood CAR T cells (pct_tscm_pb, pct_tcm_pb, pct_tem_pb, pct_teff_pb); cells/uL (cart_pb); transgene copies/ug genomic DNA (transgene)"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Scales the blood (5 L) and bone-marrow (1.6 L) volumes linearly from a 70 kg adult (Appendix S2 code: Vb = 5 * (WT/70)). The authors assumed 30 kg for every patient.",
      source_name = "WT"
    ),
    POLYCLONALITY_HIGH = list(
      description = "High polyclonality (clonal diversity) of the CAR T-cell infusion product",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = low polyclonality",
      notes = "Categorical covariate on the differentiation rate constants Rm1, Rm2 and Rm3. In the source the four patients with immunophenotype data split into two with high product polyclonality and long CAR T-cell persistence (patients 4 and 6 of Biasco et al. 2021) and two with low polyclonality and limited persistence (patients 10 and 17).",
      source_name = "polyclonality"
    )
  )

  compartmentData <- list(
    tscm_cd8_pb = list(
      analyte = "CD8+ stem-cell memory CAT19 CAR T cells",
      units = "cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    tcm_cd8_pb = list(
      analyte = "CD8+ central memory CAT19 CAR T cells",
      units = "cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    tem_cd8_pb = list(
      analyte = "CD8+ effector memory CAT19 CAR T cells",
      units = "cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    teff_cd8_pb = list(
      analyte = "CD8+ terminal effector CAT19 CAR T cells",
      units = "cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    tscm_cd4_pb = list(
      analyte = "CD4+ stem-cell memory CAT19 CAR T cells",
      units = "cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    tcm_cd4_pb = list(
      analyte = "CD4+ central memory CAT19 CAR T cells",
      units = "cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    tem_cd4_pb = list(
      analyte = "CD4+ effector memory CAT19 CAR T cells",
      units = "cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    teff_cd4_pb = list(
      analyte = "CD4+ terminal effector CAT19 CAR T cells",
      units = "cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    tscm_cd8_bm = list(
      analyte = "CD8+ stem-cell memory CAT19 CAR T cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    tcm_cd8_bm = list(
      analyte = "CD8+ central memory CAT19 CAR T cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    tem_cd8_bm = list(
      analyte = "CD8+ effector memory CAT19 CAR T cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    teff_cd8_bm = list(
      analyte = "CD8+ terminal effector CAT19 CAR T cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    tscm_cd4_bm = list(
      analyte = "CD4+ stem-cell memory CAT19 CAR T cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    tcm_cd4_bm = list(
      analyte = "CD4+ central memory CAT19 CAR T cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    tem_cd4_bm = list(
      analyte = "CD4+ effector memory CAT19 CAR T cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    teff_cd4_bm = list(
      analyte = "CD4+ terminal effector CAT19 CAR T cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_cd8_pb = list(
      analyte = "CD8+ CAR-CD19 complexes",
      units = "complexes/L",
      specimen = "whole blood",
      verified = TRUE
    ),
    complex_cd4_pb = list(
      analyte = "CD4+ CAR-CD19 complexes",
      units = "complexes/L",
      specimen = "whole blood",
      verified = TRUE
    ),
    complex_cd8_bm = list(
      analyte = "CD8+ CAR-CD19 complexes",
      units = "complexes/L",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_cd4_bm = list(
      analyte = "CD4+ CAR-CD19 complexes",
      units = "complexes/L",
      specimen = "tissue",
      verified = TRUE
    ),
    tumor_pb = list(analyte = "CD19+ leukaemic B cells", units = "cells/L", specimen = "whole blood", verified = TRUE),
    tumor_bm = list(analyte = "CD19+ leukaemic B cells", units = "cells/L", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 14L,
    n_studies = 1L,
    age_range = "children (median age about 9 years)",
    weight_range = "not reported; 30 kg assumed for every patient",
    disease_state = "Relapsed / refractory paediatric B-cell acute lymphoblastic leukaemia (CARPALL study)",
    dose_range = "Single infusion of a median 1 x 10^6 CAR+ T cells/kg",
    regions = "United Kingdom",
    notes = paste(
      "Cellular-kinetic data from all 14 patients (Ghorashian et al. 2019) and",
      "longitudinal immunophenotype percentages (pre-infusion product and days",
      "14 to 1080) from 4 of them (2 with long and 2 with short persistence;",
      "Biasco et al. 2021). Random effects other than those on Rm1-Rm3 were",
      "fixed to the step-2 CAT19 estimates; 20% variability on Rm1-Rm3 was",
      "assumed."
    )
  )

  ini({
    # ---- Estimated (Table 2, step 3) ----
    lkexp_max_cd8 <- log(3.44)
    label("Maximum first-order expansion rate constant of CD8+ CAR T cells, KmaxExpCD8 (1/day)") # Table 2 step 3: KmaxExpCD8 = 3.44 (23.3% RSE)
    lrm1 <- log(0.26)
    label("First-order TSCM-to-TCM differentiation rate constant, Rm1 (1/day)") # Table 2 step 3: Rm1 = 0.26 (21% RSE)
    lrm2 <- log(0.4)
    label("First-order TCM-to-TEM differentiation rate constant, Rm2 (1/day)") # Table 2 step 3: Rm2 = 0.4 (17.3% RSE)
    lrm3 <- log(0.29)
    label("First-order TEM-to-TEFF differentiation rate constant, Rm3 (1/day)") # Table 2 step 3: Rm3 = 0.29 (19% RSE)
    e_polyclonality_high_rm1 <- 0.08
    label("Multiplicative effect of high product polyclonality on Rm1 (unitless)") # Table 2 step 3: beta_polyclonality1 = 0.08 (10% RSE); Discussion: Rm1 reduced ~92%
    e_polyclonality_high_rm2 <- 0.21
    label("Multiplicative effect of high product polyclonality on Rm2 (unitless)") # Table 2 step 3: beta_polyclonality2 = 0.21 (14.5% RSE); Discussion: Rm2 reduced ~79%
    e_polyclonality_high_rm3 <- 0.19
    label("Multiplicative effect of high product polyclonality on Rm3 (unitless)") # Table 2 step 3: beta_polyclonality3 = 0.19 (13.7% RSE); Discussion: Rm3 reduced ~81%

    # ---- Fixed (Table 2, step 3) ----
    kexp_ratio_cd8_cd4 <- fixed(1.12)
    label("Ratio of the CD8+ to CD4+ maximum expansion rate constants (unitless)") # Table 2 step 3: KmaxExpCD4 = KmaxExpCD8/1.12 (fixed from Liso-cel estimates)
    pct_cd8 <- fixed(38)
    label("Percentage of CD8+ CAR T cells in the infused product (%)") # Table 2 step 3: PCTCD8 = 38 (fixed, ref 20)
    pct_cd4 <- fixed(54)
    label("Percentage of CD4+ CAR T cells in the infused product (%)") # Table 2 step 3: PCTCD4 = 54 (fixed, ref 20)
    lec50_exp <- fixed(log(10))
    label("CAR-target complexes per CAR T cell giving half-maximal expansion, EC50Exp (complexes/cell)") # Table 2 step 3: EC50Exp = 10 (fixed, ref 19 Singh 2021)
    lkel_e <- fixed(log(113))
    label("First-order elimination rate constant of TEFF CAR T cells from blood, Kele (1/day)") # Table 2 step 3: Kele = 113 (fixed, ref 19 Singh 2021)
    lkel_m <- fixed(log(0.013))
    label("First-order elimination rate constant of TSCM / TCM / TEM CAR T cells from blood, Kelm (1/day)") # Table 2 step 3: Kelm = 0.013 (fixed from CAT19 estimates)
    lk12 <- fixed(log(0.11))
    label("Blood-to-bone-marrow distribution rate constant, K12 (1/day)") # Table 2 step 3: K12 = 0.11 (fixed from CAT19 estimates)
    lk21 <- fixed(log(0.176))
    label("Bone-marrow-to-blood redistribution rate constant, K21 (1/day)") # Table 2 step 3: K21 = 0.176 (fixed from CAT19 estimates)
    lkkill_max <- fixed(log(1.8))
    label("Maximum first-order tumour-cell killing rate constant, KKillMax (1/day)") # Table 2 step 3: KKillMax = 1.8 (fixed from Liso-cel estimates)
    lkc50_kill <- fixed(log(0.28))
    label("CAR-target complexes per tumour cell giving half-maximal killing, KC50CAR-T (complexes/cell)") # Table 2 step 3: KC50CAR-T = 0.28 (fixed, CAT19 in vitro estimate)
    kon <- fixed(2.15e5)
    label("CAR-CD19 association rate constant, Kon (1/M/s)") # Table 2 step 3: Kon = 2.15E5 (fixed, ref 20)
    koff <- fixed(0.003)
    label("CAR-CD19 dissociation rate constant, Koff (1/s)") # Table 2 step 3: Koff = 0.003 (fixed, ref 20)
    ag_car <- fixed(15000)
    label("CAR density on CAR T cells, AgCAR (CARs/cell)") # Table 2 step 3: AgCAR = 15,000 (fixed, refs 19, 41)
    ag_tumor <- fixed(4000)
    label("CD19 density on leukaemic cells, AgTumor (antigens/cell)") # Table 2 step 3: AgTumor = 4000 (fixed, ref 39)
    dt_tumor <- fixed(10.7)
    label("Paediatric ALL tumour doubling time, DTTumor (day)") # Table 2 step 3: DTTumor = 10.7 (fixed, ref 49)
    kmig <- fixed(0.001)
    label("Bone-marrow-to-blood migration rate constant of tumour cells, Kmigration (1/day)") # Table 2 step 3: Kmigration = 0.001 (fixed, ref 44)
    tumor_bm0 <- fixed(1.06e11)
    label("Initial tumour burden in bone marrow, ICTumorBM (cells/L)") # Table 2 step 3: ICTumorBM = 1.06E11 (fixed, refs 45, 47)
    tumor_pb0 <- fixed(1.8e9)
    label("Initial tumour burden in blood, ICTumorPB (cells/L)") # Table 2 step 3: ICTumorPB = 1.8E9 (fixed, ref 50)
    v_blood <- fixed(5)
    label("Blood compartment volume for a 70 kg adult, Vb (L)") # Table 2 step 3: Vb = 5 x 30/70 (fixed, ref 40)
    v_bonemarrow <- fixed(1.6)
    label("Bone marrow compartment volume for a 70 kg adult, Vbm (L)") # Table 2 step 3: Vbm = 1.6 x 30/70 (fixed, ref 40)
    int_transgene <- fixed(7.26)
    label("Intercept of ln(transgene copies/ug gDNA) on ln(CAR T cells/uL), alpha") # Table 2 step 2 CAT19: alpha = 7.26 (carried to step 3)
    slope_transgene <- fixed(0.74)
    label("Slope of ln(transgene copies/ug gDNA) on ln(CAR T cells/uL), beta (unitless)") # Table 2 step 2 CAT19: beta = 0.74 (carried to step 3)

    # ---- Between-subject variability (all fixed in step 3) ----
    # Monolix omega is the SD of a log-normal random effect; variance = omega^2.
    etalkexp_max_cd8 ~ fixed(0.16) # Table 2 step 3: omega KmaxExpCD8 = 0.4 (not estimated); 0.4^2
    etalrm1 ~ fixed(0.04) # Table 2 step 3: omega Rm1 = 0.2 (assumed 20%, not estimated); 0.2^2
    etalrm2 ~ fixed(0.04) # Table 2 step 3: omega Rm2 = 0.2 (assumed 20%, not estimated); 0.2^2
    etalrm3 ~ fixed(0.04) # Table 2 step 3: omega Rm3 = 0.2 (assumed 20%, not estimated); 0.2^2
    etalkel_m ~ fixed(1.8225) # Table 2 step 3: omega Kelm = 1.35 (not estimated); 1.35^2
    etalk12 ~ fixed(4.3264) # Table 2 step 3: omega K12 = 2.08 (not estimated); 2.08^2
    etalk21 ~ fixed(1.1664) # Table 2 step 3: omega K21 = 1.08 (not estimated); 1.08^2
  })

  model({
    # ---- Individual parameters ----
    kexp_max_cd8 <- exp(lkexp_max_cd8 + etalkexp_max_cd8)
    kexp_max_cd4 <- kexp_max_cd8 / kexp_ratio_cd8_cd4
    rm1 <- exp(lrm1 + etalrm1) * e_polyclonality_high_rm1^POLYCLONALITY_HIGH
    rm2 <- exp(lrm2 + etalrm2) * e_polyclonality_high_rm2^POLYCLONALITY_HIGH
    rm3 <- exp(lrm3 + etalrm3) * e_polyclonality_high_rm3^POLYCLONALITY_HIGH
    kel_m <- exp(lkel_m + etalkel_m)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    ec50_exp <- exp(lec50_exp)
    kel_e <- exp(lkel_e)
    kkill_max <- exp(lkkill_max)
    kc50_kill <- exp(lkc50_kill)
    kg_tumor <- log(2) / dt_tumor

    # ---- Weight-scaled volumes (Appendix S2 code) ----
    vb <- v_blood * (WT / 70)
    vbm <- v_bonemarrow * (WT / 70)

    # ---- Unit conversion of the binding constants (Appendix S2 code) ----
    kon_d <- kon * (60 * 60 * 24) / 6.023e23
    koff_d <- koff * (60 * 60 * 24)

    # ---- Total CAR T cells per subset and compartment (cells/L) ----
    cart_cd8_pb_l <- (tscm_cd8_pb + tcm_cd8_pb + tem_cd8_pb + teff_cd8_pb) / vb
    cart_cd4_pb_l <- (tscm_cd4_pb + tcm_cd4_pb + tem_cd4_pb + teff_cd4_pb) / vb
    cart_cd8_bm_l <- (tscm_cd8_bm + tcm_cd8_bm + tem_cd8_bm + teff_cd8_bm) / vbm
    cart_cd4_bm_l <- (tscm_cd4_bm + tcm_cd4_bm + tem_cd4_bm + teff_cd4_bm) / vbm

    # ---- Blood binding and drivers (Appendix S1 Eqs. 51-57, per subset) ----
    car_free_cd8_pb <- cart_cd8_pb_l * ag_car - complex_cd8_pb
    car_free_cd4_pb <- cart_cd4_pb_l * ag_car - complex_cd4_pb
    ag_free_pb <- tumor_pb * ag_tumor - (complex_cd8_pb + complex_cd4_pb)
    cplx_per_cart_cd8_pb <- complex_cd8_pb / cart_cd8_pb_l
    cplx_per_cart_cd4_pb <- complex_cd4_pb / cart_cd4_pb_l
    kexp_cd8_pb <- kexp_max_cd8 * cplx_per_cart_cd8_pb / (ec50_exp + cplx_per_cart_cd8_pb)
    kexp_cd4_pb <- kexp_max_cd4 * cplx_per_cart_cd4_pb / (ec50_exp + cplx_per_cart_cd4_pb)
    # Eq. 57 x Tumor^PB, written so that an eradicated tumour gives zero kill
    # flux instead of 0/0.
    cplx_pb <- complex_cd8_pb + complex_cd4_pb
    kill_flux_pb <- kkill_max * cplx_pb * tumor_pb / (kc50_kill * tumor_pb + cplx_pb)

    # ---- Bone-marrow binding and drivers (Appendix S1 Eqs. 59-65, per subset) ----
    car_free_cd8_bm <- cart_cd8_bm_l * ag_car - complex_cd8_bm
    car_free_cd4_bm <- cart_cd4_bm_l * ag_car - complex_cd4_bm
    ag_free_bm <- tumor_bm * ag_tumor - (complex_cd8_bm + complex_cd4_bm)
    cplx_per_cart_cd8_bm <- complex_cd8_bm / cart_cd8_bm_l
    cplx_per_cart_cd4_bm <- complex_cd4_bm / cart_cd4_bm_l
    kexp_cd8_bm <- kexp_max_cd8 * cplx_per_cart_cd8_bm / (ec50_exp + cplx_per_cart_cd8_bm)
    kexp_cd4_bm <- kexp_max_cd4 * cplx_per_cart_cd4_bm / (ec50_exp + cplx_per_cart_cd4_bm)
    cplx_bm <- complex_cd8_bm + complex_cd4_bm
    kill_flux_bm <- kkill_max * cplx_bm * tumor_bm / (kc50_kill * tumor_bm + cplx_bm)

    # ---- CD8+ CAR T cells in blood (Appendix S1 Eqs. 43-46, cell numbers) ----
    d/dt(tscm_cd8_pb) <- -k12 * tscm_cd8_pb + k21 * tscm_cd8_bm - kel_m * tscm_cd8_pb +
      kexp_cd8_pb * tscm_cd8_pb - rm1 * tscm_cd8_pb
    d/dt(tcm_cd8_pb) <- -k12 * tcm_cd8_pb + k21 * tcm_cd8_bm - kel_m * tcm_cd8_pb +
      kexp_cd8_pb * tcm_cd8_pb + rm1 * tscm_cd8_pb - rm2 * tcm_cd8_pb
    d/dt(tem_cd8_pb) <- -k12 * tem_cd8_pb + k21 * tem_cd8_bm - kel_m * tem_cd8_pb +
      kexp_cd8_pb * tem_cd8_pb + rm2 * tcm_cd8_pb - rm3 * tem_cd8_pb
    d/dt(teff_cd8_pb) <- -k12 * teff_cd8_pb + k21 * teff_cd8_bm - kel_e * teff_cd8_pb +
      rm3 * tem_cd8_pb

    # ---- CD4+ CAR T cells in blood ----
    d/dt(tscm_cd4_pb) <- -k12 * tscm_cd4_pb + k21 * tscm_cd4_bm - kel_m * tscm_cd4_pb +
      kexp_cd4_pb * tscm_cd4_pb - rm1 * tscm_cd4_pb
    d/dt(tcm_cd4_pb) <- -k12 * tcm_cd4_pb + k21 * tcm_cd4_bm - kel_m * tcm_cd4_pb +
      kexp_cd4_pb * tcm_cd4_pb + rm1 * tscm_cd4_pb - rm2 * tcm_cd4_pb
    d/dt(tem_cd4_pb) <- -k12 * tem_cd4_pb + k21 * tem_cd4_bm - kel_m * tem_cd4_pb +
      kexp_cd4_pb * tem_cd4_pb + rm2 * tcm_cd4_pb - rm3 * tem_cd4_pb
    d/dt(teff_cd4_pb) <- -k12 * teff_cd4_pb + k21 * teff_cd4_bm - kel_e * teff_cd4_pb +
      rm3 * tem_cd4_pb

    # ---- CD8+ CAR T cells in bone marrow (Appendix S1 Eqs. 47-50) ----
    d/dt(tscm_cd8_bm) <- k12 * tscm_cd8_pb - k21 * tscm_cd8_bm + kexp_cd8_bm * tscm_cd8_bm -
      rm1 * tscm_cd8_bm
    d/dt(tcm_cd8_bm) <- k12 * tcm_cd8_pb - k21 * tcm_cd8_bm + kexp_cd8_bm * tcm_cd8_bm +
      rm1 * tscm_cd8_bm - rm2 * tcm_cd8_bm
    d/dt(tem_cd8_bm) <- k12 * tem_cd8_pb - k21 * tem_cd8_bm + kexp_cd8_bm * tem_cd8_bm +
      rm2 * tcm_cd8_bm - rm3 * tem_cd8_bm
    d/dt(teff_cd8_bm) <- k12 * teff_cd8_pb - k21 * teff_cd8_bm + rm3 * tem_cd8_bm

    # ---- CD4+ CAR T cells in bone marrow ----
    d/dt(tscm_cd4_bm) <- k12 * tscm_cd4_pb - k21 * tscm_cd4_bm + kexp_cd4_bm * tscm_cd4_bm -
      rm1 * tscm_cd4_bm
    d/dt(tcm_cd4_bm) <- k12 * tcm_cd4_pb - k21 * tcm_cd4_bm + kexp_cd4_bm * tcm_cd4_bm +
      rm1 * tscm_cd4_bm - rm2 * tcm_cd4_bm
    d/dt(tem_cd4_bm) <- k12 * tem_cd4_pb - k21 * tem_cd4_bm + kexp_cd4_bm * tem_cd4_bm +
      rm2 * tcm_cd4_bm - rm3 * tem_cd4_bm
    d/dt(teff_cd4_bm) <- k12 * teff_cd4_pb - k21 * teff_cd4_bm + rm3 * tem_cd4_bm

    # ---- CAR-CD19 complexes (Appendix S1 Eqs. 53 and 61, per subset) ----
    d/dt(complex_cd8_pb) <- kon_d * car_free_cd8_pb * ag_free_pb - koff_d * complex_cd8_pb
    d/dt(complex_cd4_pb) <- kon_d * car_free_cd4_pb * ag_free_pb - koff_d * complex_cd4_pb
    d/dt(complex_cd8_bm) <- kon_d * car_free_cd8_bm * ag_free_bm - koff_d * complex_cd8_bm
    d/dt(complex_cd4_bm) <- kon_d * car_free_cd4_bm * ag_free_bm - koff_d * complex_cd4_bm

    # ---- CD19+ tumour cells (Appendix S1 Eqs. 58 and 66) ----
    d/dt(tumor_pb) <- kg_tumor * tumor_pb - kill_flux_pb + kmig * tumor_bm
    d/dt(tumor_bm) <- kg_tumor * tumor_bm - kill_flux_bm - kmig * tumor_bm

    # ---- CD4 / CD8 split of each immunophenotype dose (Table 2 PCTCD8 / PCTCD4) ----
    f(tscm_cd8_pb) <- pct_cd8 / 100
    f(tcm_cd8_pb) <- pct_cd8 / 100
    f(tem_cd8_pb) <- pct_cd8 / 100
    f(teff_cd8_pb) <- pct_cd8 / 100
    f(tscm_cd4_pb) <- pct_cd4 / 100
    f(tcm_cd4_pb) <- pct_cd4 / 100
    f(tem_cd4_pb) <- pct_cd4 / 100
    f(teff_cd4_pb) <- pct_cd4 / 100

    # ---- Initial conditions (Appendix S2 code) ----
    # 1E-12 cells/L in the bone-marrow pools keeps the complexes-per-CAR-T
    # ratios finite before any cell has reached the marrow.
    tscm_cd8_bm(0) <- 1e-12 * vbm
    tcm_cd8_bm(0) <- 1e-12 * vbm
    tem_cd8_bm(0) <- 1e-12 * vbm
    teff_cd8_bm(0) <- 1e-12 * vbm
    tscm_cd4_bm(0) <- 1e-12 * vbm
    tcm_cd4_bm(0) <- 1e-12 * vbm
    tem_cd4_bm(0) <- 1e-12 * vbm
    teff_cd4_bm(0) <- 1e-12 * vbm
    complex_cd8_pb(0) <- 1e-12
    complex_cd4_pb(0) <- 1e-12
    complex_cd8_bm(0) <- 1e-12
    complex_cd4_bm(0) <- 1e-12
    tumor_pb(0) <- tumor_pb0
    tumor_bm(0) <- tumor_bm0

    # ---- Outputs (Appendix S1 text after Eq. 46; 1E-6 converts cells/L to cells/uL) ----
    cart_all_pb <- tscm_cd8_pb + tcm_cd8_pb + tem_cd8_pb + teff_cd8_pb +
      tscm_cd4_pb + tcm_cd4_pb + tem_cd4_pb + teff_cd4_pb
    pct_tscm_pb <- 100 * (tscm_cd8_pb + tscm_cd4_pb) / cart_all_pb
    pct_tcm_pb <- 100 * (tcm_cd8_pb + tcm_cd4_pb) / cart_all_pb
    pct_tem_pb <- 100 * (tem_cd8_pb + tem_cd4_pb) / cart_all_pb
    pct_teff_pb <- 100 * (teff_cd8_pb + teff_cd4_pb) / cart_all_pb
    cart_pb <- cart_all_pb / vb * 1e-6
    transgene <- exp(int_transgene + slope_transgene * log(cart_pb))
  })
}
