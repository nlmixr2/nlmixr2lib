# Multiscale mechanistic cellular-kinetic / pharmacodynamic (CK-PD) model of
# lisocabtagene maraleucel (Liso-cel, FMC63 binder, 1:1 CD4:CD8 product) with
# separate CD4+ and CD8+ CAR T-cell subsets (step 2 of Salem 2023). It extends
# the Singh 2021 anti-BCMA framework (Singh_2021_idecabtageneVicleucel_human.R)
# with a CD4/CD8 split, blood and bone-marrow CAR-target binding, and a CD19+
# tumour pool in both compartments. Fit to summary-level (mean / median) flow
# cytometry and qPCR profiles of adults with DLBCL and ALL.
#
# SOURCES. Equations: Appendix S1 Eqs. 6-42 and the Monolix code for Figure 2
# in Appendix S2. Parameter values: main-text Table 2, "Parameters associated
# with clinical CK-PD model accounting for differential CD4/CD8 expansion
# (Liso-cel)". Appendix S1 Eqs. 7 and 10 omit the memory-cell elimination term
# that the text ("first-order elimination rate constant (Kelm) occurs in PB
# compartment only") and the Monolix code both carry; it is included here.
#
# STATE UNITS. The Monolix code carries the CAR T-cell pools as concentrations
# (cells/L of blood or bone marrow) and doses cells/L. They are carried here
# as cell NUMBERS (concentration x compartment volume), so a dose record is a
# number of CAR+ T cells; every CAR T-cell flux in the code is a rate constant
# times a concentration times a fixed volume, so the transformation is exact.
# Complexes and tumour cells stay as the code's concentrations (#/L, cells/L).
#
# LONG HORIZONS. After the tumour is eradicated its states sit at solver noise,
# which the exponential growth term amplifies; beyond roughly 1.5-2 years after
# eradication the trajectories depend on the solver tolerances (see the
# vignette, "Long-horizon simulations").

Salem_2023_lisocabtageneMaraleucel_qsp <- function() {
  description <- "QSP (multiscale mechanistic cellular-kinetic / pharmacodynamic model). Anti-CD19 CAR T-cell therapy lisocabtagene maraleucel (Liso-cel, 1:1 CD4:CD8 product) in adults with diffuse large B-cell lymphoma (DLBCL) or acute lymphoblastic leukaemia (ALL), Salem 2023 step 2. Effector and memory CAR T cells of the CD4+ and CD8+ subsets distribute between blood and bone marrow (first-order K12 / K21), effector cells convert to memory at a net first-order rate, and both phenotypes are eliminated from blood only. In each compartment the CARs bind CD19 on tumour cells (second-order kon, first-order koff); complexes per CAR T cell drive subset-specific Emax expansion of effector cells and complexes per tumour cell drive Emax killing of exponentially growing CD19+ tumour cells, which migrate from bone marrow to blood. Disease-specific tumour growth, antigen density and initial burden switch on DIS_ALL. Outputs: blood CAR+ CD4, CD8 and CD3 cells/uL (flow cytometry) and transgene copies/ug genomic DNA (qPCR, log-linear link)."
  reference <- paste(
    "Salem AM, Mugundu GM, Singh AP (2023).",
    "Development of a multiscale mechanistic modeling framework integrating",
    "differential cellular kinetics of CAR T-cell subsets and immunophenotypes",
    "in cancer patients.",
    "CPT Pharmacometrics Syst Pharmacol 12(9):1285-1296.",
    "doi:10.1002/psp4.13009.",
    "Equations in Appendix S1 and Monolix code in Appendix S2.",
    "Fit to data digitised from Turtle et al. 2016 (J Clin Invest 126:2123-2138)",
    "and the Breyanzi (lisocabtagene maraleucel) FDA BLA review (2021).",
    sep = " "
  )
  vignette <- "Salem_2023_CART_subsets_immunophenotypes"

  # CAR T-cell effector / memory pools per subset and compartment, the
  # per-subset complexes and the two tumour pools are paper-mechanistic
  # states with no canonical analogue (Singh_2021_idecabtageneVicleucel_human.R
  # uses carte_pb / cartm_pb / carte_t / cartm_t for the subset-free form).
  paper_specific_compartments <- c(
    "carte_cd8_pb",
    "cartm_cd8_pb",
    "carte_cd4_pb",
    "cartm_cd4_pb",
    "carte_cd8_bm",
    "cartm_cd8_bm",
    "carte_cd4_bm",
    "cartm_cd4_bm",
    "complex_cd8_pb",
    "complex_cd4_pb",
    "complex_cd8_bm",
    "complex_cd4_bm",
    "tumor_pb",
    "tumor_bm"
  )

  units <- list(
    time = "day",
    dosing = "cells (total CAR+ T cells infused, given as two records with the same amt: one on carte_cd8_pb and one on carte_cd4_pb; f() applies the 50% / 50% subset split)",
    concentration = "cells/uL (cart_cd8_pb, cart_cd4_pb, cart_cd3_pb); transgene copies/ug genomic DNA (transgene)"
  )

  covariateData <- list(
    DIS_ALL = list(
      description = "Acute lymphoblastic leukaemia indicator; reference category is diffuse large B-cell lymphoma",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = DLBCL (Breyanzi BLA cohort)",
      notes = "Switches the disease-specific system parameters of Table 2: CD19 density (4000 ALL vs 20,000 DLBCL antigens/cell), tumour doubling time (5.8 vs 29 days) and initial tumour burden in bone marrow (1.06E11 vs 5E10 cells/L) and blood (5.5E9 vs 0 cells/L). Supplied to the Monolix fit as per-dataset regressors (Appendix S2: Tumor_BS_PB, DT, TAA_density, Tumor_BS_T).",
      source_name = "disease indication (ALL vs DLBCL)"
    )
  )

  compartmentData <- list(
    carte_cd8_pb = list(
      analyte = "effector CD8+ lisocabtagene maraleucel CAR T cells",
      units = "cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    cartm_cd8_pb = list(
      analyte = "memory CD8+ lisocabtagene maraleucel CAR T cells",
      units = "cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    carte_cd4_pb = list(
      analyte = "effector CD4+ lisocabtagene maraleucel CAR T cells",
      units = "cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    cartm_cd4_pb = list(
      analyte = "memory CD4+ lisocabtagene maraleucel CAR T cells",
      units = "cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    carte_cd8_bm = list(
      analyte = "effector CD8+ lisocabtagene maraleucel CAR T cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    cartm_cd8_bm = list(
      analyte = "memory CD8+ lisocabtagene maraleucel CAR T cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    carte_cd4_bm = list(
      analyte = "effector CD4+ lisocabtagene maraleucel CAR T cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    cartm_cd4_bm = list(
      analyte = "memory CD4+ lisocabtagene maraleucel CAR T cells",
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
    tumor_pb = list(analyte = "CD19+ malignant B cells", units = "cells/L", specimen = "whole blood", verified = TRUE),
    tumor_bm = list(analyte = "CD19+ malignant B cells", units = "cells/L", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 2L,
    age_range = "adults",
    weight_range = "not reported; 70 kg assumed to convert the ALL cells/kg doses",
    disease_state = "Relapsed / refractory diffuse large B-cell lymphoma (Breyanzi BLA, TRANSCEND NHL 001) and adult B-cell acute lymphoblastic leukaemia (Turtle et al. 2016)",
    dose_range = "DLBCL: single infusion of 50, 100 or 150 x 10^6 CAR+ T cells; ALL: 2 x 10^5 or 2 x 10^6 CAR+ T cells/kg",
    regions = "United States",
    notes = paste(
      "Summary-level data only (Table 1): mean +/- SEM blood CAR+CD4+ and",
      "CAR+CD8+ flow cytometry from Turtle 2016 (ALL, n = 13 and 15 per dose",
      "level), and median blood CAR+CD3+/CD4+/CD8+ flow cytometry (50 x 10^6",
      "cells, n = 45) and median qPCR (n = 45, 176 and 41 at 50, 100 and",
      "150 x 10^6 cells) from the Breyanzi BLA review (DLBCL). The product",
      "contains CD4+ and CD8+ CAR T cells in a 1:1 ratio."
    )
  )

  ini({
    # ---- Estimated (Table 2, Liso-cel) ----
    lkexp_max_cd8 <- log(2.14)
    label("Maximum first-order expansion rate constant of CD8+ CAR T cells, KmaxExpCD8 (1/day)") # Table 2 Liso-cel: KmaxExpCD8 = 2.14 (5.4% RSE)
    lkexp_max_cd4 <- log(1.92)
    label("Maximum first-order expansion rate constant of CD4+ CAR T cells, KmaxExpCD4 (1/day)") # Table 2 Liso-cel: KmaxExpCD4 = 1.92 (5.8% RSE)
    lrm <- log(0.000007)
    label("Net first-order effector-to-memory conversion rate constant, Rm (1/day)") # Table 2 Liso-cel: Rm = 0.000007 (52% RSE)
    lkel_m <- log(0.024)
    label("First-order elimination rate constant of memory CAR T cells from blood, Kelm (1/day)") # Table 2 Liso-cel: Kelm = 0.024 (44% RSE)
    lk12 <- log(0.4)
    label("Blood-to-bone-marrow distribution rate constant, K12 (1/day)") # Table 2 Liso-cel: K12 = 0.4 (51.8% RSE)
    lk21 <- log(0.15)
    label("Bone-marrow-to-blood redistribution rate constant, K21 (1/day)") # Table 2 Liso-cel: K21 = 0.15 (14.1% RSE)
    lkkill_max <- log(1.8)
    label("Maximum first-order tumour-cell killing rate constant, KKillMax (1/day)") # Table 2 Liso-cel: KKillMax = 1.8 (21.3% RSE)

    # ---- Fixed (Table 2, Liso-cel) ----
    lec50_exp <- fixed(log(10))
    label("CAR-target complexes per CAR T cell giving half-maximal expansion, EC50Exp (complexes/cell)") # Table 2 Liso-cel: EC50Exp = 10 (fixed, ref 19 Singh 2021)
    lkel_e <- fixed(log(113))
    label("First-order elimination rate constant of effector CAR T cells from blood, Kele (1/day)") # Table 2 Liso-cel: Kele = 113 (fixed, ref 19 Singh 2021)
    lkc50_kill <- fixed(log(9.16))
    label("CAR-target complexes per tumour cell giving half-maximal killing, KC50CAR-T (complexes/cell)") # Table 2 Liso-cel: KC50CAR-T = 9.16 (fixed, FMC63 in vitro estimate)
    kon <- fixed(2.07e5)
    label("CAR-CD19 association rate constant, Kon (1/M/s)") # Table 2 Liso-cel: Kon = 2.07E5 (fixed, ref 20)
    koff <- fixed(6.81e-5)
    label("CAR-CD19 dissociation rate constant, Koff (1/s)") # Table 2 Liso-cel: Koff = 6.81E-5 (fixed, ref 20)
    ag_car <- fixed(15000)
    label("CAR density on CAR T cells, AgCAR (CARs/cell)") # Table 2 Liso-cel: AgCAR = 15,000 (fixed, refs 19, 41)
    ag_tumor <- fixed(20000)
    label("CD19 density on DLBCL tumour cells, AgTumor (antigens/cell)") # Table 2 Liso-cel: AgTumor = 20,000 (DLBCL) (fixed, ref 39)
    ag_tumor_all <- fixed(4000)
    label("CD19 density on ALL tumour cells, AgTumor (antigens/cell)") # Table 2 Liso-cel: AgTumor = 4000 (ALL) (fixed, ref 39)
    dt_tumor <- fixed(29)
    label("DLBCL tumour doubling time, DTTumor (day)") # Table 2 Liso-cel: DTTumor = 29 (DLBCL) (fixed, ref 42)
    dt_tumor_all <- fixed(5.8)
    label("ALL tumour doubling time, DTTumor (day)") # Table 2 Liso-cel: DTTumor = 5.8 (ALL) (fixed, ref 43)
    kmig <- fixed(0.001)
    label("Bone-marrow-to-blood migration rate constant of tumour cells, Kmigration (1/day)") # Table 2 Liso-cel: Kmigration = 0.001 (fixed, ref 44)
    tumor_bm0 <- fixed(5e10)
    label("Initial DLBCL tumour burden in bone marrow, ICTumorBM (cells/L)") # Table 2 Liso-cel: ICTumorBM = 5E10 (DLBCL) (fixed, refs 45-47)
    tumor_bm0_all <- fixed(1.06e11)
    label("Initial ALL tumour burden in bone marrow, ICTumorBM (cells/L)") # Table 2 Liso-cel: ICTumorBM = 1.06E11 (ALL) (fixed, refs 45-47)
    tumor_pb0 <- fixed(0)
    label("Initial DLBCL tumour burden in blood, ICTumorPB (cells/L)") # Table 2 Liso-cel: ICTumorPB = 0 (DLBCL) (fixed; binding restricted to bone marrow)
    tumor_pb0_all <- fixed(5.5e9)
    label("Initial ALL tumour burden in blood, ICTumorPB (cells/L)") # Table 2 Liso-cel: ICTumorPB = 5.5E9 (ALL) (fixed, ref 48)
    v_blood <- fixed(5)
    label("Blood compartment volume, Vb (L)") # Table 2 Liso-cel: Vb = 5 (fixed, ref 40)
    v_bonemarrow <- fixed(1.6)
    label("Bone marrow compartment volume, Vbm (L)") # Table 2 Liso-cel: Vbm = 1.6 (fixed, ref 40)
    pct_cd8 <- fixed(50)
    label("Percentage of CD8+ CAR T cells in the infused product (%)") # Table 1: Liso-cel CD4/CD8 percentage 50%:50%
    pct_cd4 <- fixed(50)
    label("Percentage of CD4+ CAR T cells in the infused product (%)") # Table 1: Liso-cel CD4/CD8 percentage 50%:50%
    int_transgene <- fixed(7.6204)
    label("Intercept of ln(transgene copies/ug gDNA) on ln(CAR+ CD3 cells/uL), alpha") # Table 2 Liso-cel: alpha = 7.62 (fixed, ref 24); Appendix S2 code 7.6204
    slope_transgene <- fixed(0.8537)
    label("Slope of ln(transgene copies/ug gDNA) on ln(CAR+ CD3 cells/uL), beta (unitless)") # Table 2 Liso-cel: beta = 0.85 (fixed, ref 24); Appendix S2 code 0.8537

    # ---- Between-subject variability ----
    # Monolix omega is the SD of a log-normal random effect; variance = omega^2.
    etalkkill_max ~ fixed(0.25) # Table 2 Liso-cel: omega KKillMax = 0.5 (ALL vs DLBCL heterogeneity, not estimated); 0.5^2
  })

  model({
    # ---- Individual parameters ----
    kexp_max_cd8 <- exp(lkexp_max_cd8)
    kexp_max_cd4 <- exp(lkexp_max_cd4)
    rm <- exp(lrm)
    kel_m <- exp(lkel_m)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kkill_max <- exp(lkkill_max + etalkkill_max)
    ec50_exp <- exp(lec50_exp)
    kel_e <- exp(lkel_e)
    kc50_kill <- exp(lkc50_kill)

    # ---- Disease-specific system parameters (Table 2; DLBCL reference) ----
    ag_tum <- ag_tumor * (1 - DIS_ALL) + ag_tumor_all * DIS_ALL
    dt_tum <- dt_tumor * (1 - DIS_ALL) + dt_tumor_all * DIS_ALL
    kg_tumor <- log(2) / dt_tum

    # ---- Unit conversion of the binding constants (Appendix S2 code) ----
    # Kon (1/M/s) -> 1/(#/L)/day; Koff (1/s) -> 1/day.
    kon_d <- kon * (60 * 60 * 24) / 6.023e23
    koff_d <- koff * (60 * 60 * 24)

    # ---- CAR T-cell concentrations (cells/L) ----
    cart_cd8_pb_l <- (carte_cd8_pb + cartm_cd8_pb) / v_blood
    cart_cd4_pb_l <- (carte_cd4_pb + cartm_cd4_pb) / v_blood
    cart_cd8_bm_l <- (carte_cd8_bm + cartm_cd8_bm) / v_bonemarrow
    cart_cd4_bm_l <- (carte_cd4_bm + cartm_cd4_bm) / v_bonemarrow

    # ---- Blood binding and drivers (Appendix S1 Eqs. 19-29) ----
    car_free_cd8_pb <- cart_cd8_pb_l * ag_car - complex_cd8_pb
    car_free_cd4_pb <- cart_cd4_pb_l * ag_car - complex_cd4_pb
    ag_free_pb <- tumor_pb * ag_tum - (complex_cd8_pb + complex_cd4_pb)
    cplx_per_cart_cd8_pb <- complex_cd8_pb / cart_cd8_pb_l
    cplx_per_cart_cd4_pb <- complex_cd4_pb / cart_cd4_pb_l
    kexp_cd8_pb <- kexp_max_cd8 * cplx_per_cart_cd8_pb / (ec50_exp + cplx_per_cart_cd8_pb)
    kexp_cd4_pb <- kexp_max_cd4 * cplx_per_cart_cd4_pb / (ec50_exp + cplx_per_cart_cd4_pb)
    # Eq. 29 x Tumor^PB, written so that a zero DLBCL blood tumour burden
    # gives zero kill flux instead of 0/0.
    cplx_pb <- complex_cd8_pb + complex_cd4_pb
    kill_flux_pb <- kkill_max * cplx_pb * tumor_pb / (kc50_kill * tumor_pb + cplx_pb)

    # ---- Bone-marrow binding and drivers (Appendix S1 Eqs. 31-41) ----
    car_free_cd8_bm <- cart_cd8_bm_l * ag_car - complex_cd8_bm
    car_free_cd4_bm <- cart_cd4_bm_l * ag_car - complex_cd4_bm
    ag_free_bm <- tumor_bm * ag_tum - (complex_cd8_bm + complex_cd4_bm)
    cplx_per_cart_cd8_bm <- complex_cd8_bm / cart_cd8_bm_l
    cplx_per_cart_cd4_bm <- complex_cd4_bm / cart_cd4_bm_l
    kexp_cd8_bm <- kexp_max_cd8 * cplx_per_cart_cd8_bm / (ec50_exp + cplx_per_cart_cd8_bm)
    kexp_cd4_bm <- kexp_max_cd4 * cplx_per_cart_cd4_bm / (ec50_exp + cplx_per_cart_cd4_bm)
    cplx_bm <- complex_cd8_bm + complex_cd4_bm
    kill_flux_bm <- kkill_max * cplx_bm * tumor_bm / (kc50_kill * tumor_bm + cplx_bm)

    # ---- CAR T-cell kinetics in cell numbers (Appendix S1 Eqs. 6-17) ----
    d/dt(carte_cd8_pb) <- -k12 * carte_cd8_pb + k21 * carte_cd8_bm - kel_e * carte_cd8_pb +
      kexp_cd8_pb * carte_cd8_pb - rm * carte_cd8_pb
    d/dt(cartm_cd8_pb) <- -k12 * cartm_cd8_pb + k21 * cartm_cd8_bm - kel_m * cartm_cd8_pb +
      rm * carte_cd8_pb
    d/dt(carte_cd4_pb) <- -k12 * carte_cd4_pb + k21 * carte_cd4_bm - kel_e * carte_cd4_pb +
      kexp_cd4_pb * carte_cd4_pb - rm * carte_cd4_pb
    d/dt(cartm_cd4_pb) <- -k12 * cartm_cd4_pb + k21 * cartm_cd4_bm - kel_m * cartm_cd4_pb +
      rm * carte_cd4_pb
    d/dt(carte_cd8_bm) <- k12 * carte_cd8_pb - k21 * carte_cd8_bm + kexp_cd8_bm * carte_cd8_bm -
      rm * carte_cd8_bm
    d/dt(cartm_cd8_bm) <- k12 * cartm_cd8_pb - k21 * cartm_cd8_bm + rm * carte_cd8_bm
    d/dt(carte_cd4_bm) <- k12 * carte_cd4_pb - k21 * carte_cd4_bm + kexp_cd4_bm * carte_cd4_bm -
      rm * carte_cd4_bm
    d/dt(cartm_cd4_bm) <- k12 * cartm_cd4_pb - k21 * cartm_cd4_bm + rm * carte_cd4_bm

    # ---- CAR-CD19 complexes (Appendix S1 Eqs. 22-23, 34-35) ----
    d/dt(complex_cd8_pb) <- kon_d * car_free_cd8_pb * ag_free_pb - koff_d * complex_cd8_pb
    d/dt(complex_cd4_pb) <- kon_d * car_free_cd4_pb * ag_free_pb - koff_d * complex_cd4_pb
    d/dt(complex_cd8_bm) <- kon_d * car_free_cd8_bm * ag_free_bm - koff_d * complex_cd8_bm
    d/dt(complex_cd4_bm) <- kon_d * car_free_cd4_bm * ag_free_bm - koff_d * complex_cd4_bm

    # ---- CD19+ tumour cells (Appendix S1 Eqs. 30 and 42) ----
    d/dt(tumor_pb) <- kg_tumor * tumor_pb - kill_flux_pb + kmig * tumor_bm
    d/dt(tumor_bm) <- kg_tumor * tumor_bm - kill_flux_bm - kmig * tumor_bm

    # ---- Dose split (Table 1: 1:1 CD4:CD8 product) ----
    f(carte_cd8_pb) <- pct_cd8 / 100
    f(carte_cd4_pb) <- pct_cd4 / 100

    # ---- Initial conditions (Appendix S2 code) ----
    # 1E-12 cells/L in the empty CAR T-cell pools keeps the complexes-per-CAR-T
    # ratios finite before any cell has reached them; complexes start at 1E-12 #/L.
    cartm_cd8_pb(0) <- 1e-12 * v_blood
    cartm_cd4_pb(0) <- 1e-12 * v_blood
    carte_cd8_bm(0) <- 1e-12 * v_bonemarrow
    cartm_cd8_bm(0) <- 1e-12 * v_bonemarrow
    carte_cd4_bm(0) <- 1e-12 * v_bonemarrow
    cartm_cd4_bm(0) <- 1e-12 * v_bonemarrow
    complex_cd8_pb(0) <- 1e-12
    complex_cd4_pb(0) <- 1e-12
    complex_cd8_bm(0) <- 1e-12
    complex_cd4_bm(0) <- 1e-12
    tumor_pb(0) <- tumor_pb0 * (1 - DIS_ALL) + tumor_pb0_all * DIS_ALL
    tumor_bm(0) <- tumor_bm0 * (1 - DIS_ALL) + tumor_bm0_all * DIS_ALL

    # ---- Outputs (Appendix S1 Eqs. 8, 11-13; 1E-6 converts cells/L to cells/uL) ----
    cart_cd8_pb <- cart_cd8_pb_l * 1e-6
    cart_cd4_pb <- cart_cd4_pb_l * 1e-6
    cart_cd3_pb <- cart_cd8_pb + cart_cd4_pb
    transgene <- exp(int_transgene + slope_transgene * log(cart_cd3_pb))
  })
}
