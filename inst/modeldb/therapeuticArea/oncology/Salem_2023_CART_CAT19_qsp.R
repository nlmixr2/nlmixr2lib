# Multiscale mechanistic cellular-kinetic / pharmacodynamic (CK-PD) model of
# the low-affinity CAT19-binder anti-CD19 CAR T cell in children with acute
# lymphoblastic leukaemia (step 2 of Salem 2023, CAT19 arm). Same structure as
# Salem_2023_lisocabtageneMaraleucel_qsp.R (effector and memory CD4+ / CD8+
# CAR T cells in blood and bone marrow, CAR-CD19 binding and CD19+ tumour in
# both compartments), refit to the individual flow cytometry and qPCR data of
# Ghorashian et al. 2019 with between-subject variability. The CD4/CD8
# expansion ratio (1.12) is carried over from the Liso-cel fit.
#
# SOURCES. Equations: Appendix S1 Eqs. 6-42 and the Monolix code for Figure 2
# in Appendix S2. Parameter values: main-text Table 2, "Parameters associated
# with clinical CK-PD model accounting for differential CD4/CD8 expansion
# (CAT19)". That table reports the estimated rate constants per month; the
# per-day values used here are the ones the Results / Discussion text and the
# step-3 table print (KmaxExpCD8 5.3, Kelm 0.013, Rm 0.000037, K12 0.11 per
# day; i.e. 1 month = 30 days). Table 2's "5 * 3070" and "1.6 * 3070" volumes
# are 5 x 30/70 and 1.6 x 30/70 L, the adult volumes scaled to the 30 kg
# weight assumed for every patient; the step-3 Monolix code writes the same
# scaling as Vb = 5 * (WT/70), which is the form used here.
#
# STATE UNITS. CAR T-cell pools are cell NUMBERS (see the Liso-cel file);
# complexes and tumour cells are the code's concentrations (#/L, cells/L).
#
# LONG HORIZONS. After the tumour is eradicated its states sit at solver noise,
# which the exponential growth term amplifies; beyond roughly 1.5-2 years after
# eradication the trajectories depend on the solver tolerances (see the
# vignette, "Long-horizon simulations").

Salem_2023_CART_CAT19_qsp <- function() {
  description <- "QSP (multiscale mechanistic cellular-kinetic / pharmacodynamic model). Low-affinity CAT19-binder anti-CD19 CAR T cells in children with relapsed / refractory B-cell acute lymphoblastic leukaemia, Salem 2023 step 2. Effector and memory CAR T cells of the CD4+ and CD8+ subsets distribute between blood and bone marrow (first-order K12 / K21), effector cells convert to memory at a net first-order rate, and both phenotypes are eliminated from blood only. In each compartment the CARs bind CD19 on leukaemic cells; complexes per CAR T cell drive Emax expansion of effector cells (CD4+ rate = CD8+ rate / 1.12) and complexes per tumour cell drive Emax killing of exponentially growing CD19+ cells, which migrate from bone marrow to blood. Blood and bone-marrow volumes scale linearly with body weight. Between-subject variability on CD8+ expansion, CD4 / CD8 product fractions, conversion, memory elimination and distribution. Outputs: CAR T cells/uL and transgene copies/ug genomic DNA in blood and bone marrow."
  reference <- paste(
    "Salem AM, Mugundu GM, Singh AP (2023).",
    "Development of a multiscale mechanistic modeling framework integrating",
    "differential cellular kinetics of CAR T-cell subsets and immunophenotypes",
    "in cancer patients.",
    "CPT Pharmacometrics Syst Pharmacol 12(9):1285-1296.",
    "doi:10.1002/psp4.13009.",
    "Equations in Appendix S1 and Monolix code in Appendix S2.",
    "Fit to individual data digitised from Ghorashian et al. 2019",
    "(Nat Med 25:1408-1414; doi:10.1038/s41591-019-0549-5).",
    sep = " "
  )
  vignette <- "Salem_2023_CART_subsets_immunophenotypes"

  # Same paper-mechanistic states as Salem_2023_lisocabtageneMaraleucel_qsp.R.
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
    dosing = "cells (total CAR+ T cells infused, 1e6 cells/kg in the source, given as two records with the same amt: one on carte_cd8_pb and one on carte_cd4_pb; f() applies the 38% CD8 / 54% CD4 product fractions)",
    concentration = "cells/uL (cart_pb, cart_bm); transgene copies/ug genomic DNA (transgene, transgene_bm)"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Scales the blood (5 L) and bone-marrow (1.6 L) volumes linearly from a 70 kg adult. The individual weights were not available to the authors, who assumed 30 kg for every patient (median age about 9 years) to convert the 1e6 cells/kg dose and the volumes.",
      source_name = "WT"
    )
  )

  compartmentData <- list(
    carte_cd8_pb = list(
      analyte = "effector CD8+ CAT19 CAR T cells",
      units = "cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    cartm_cd8_pb = list(
      analyte = "memory CD8+ CAT19 CAR T cells",
      units = "cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    carte_cd4_pb = list(
      analyte = "effector CD4+ CAT19 CAR T cells",
      units = "cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    cartm_cd4_pb = list(
      analyte = "memory CD4+ CAT19 CAR T cells",
      units = "cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    carte_cd8_bm = list(
      analyte = "effector CD8+ CAT19 CAR T cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    cartm_cd8_bm = list(
      analyte = "memory CD8+ CAT19 CAR T cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    carte_cd4_bm = list(
      analyte = "effector CD4+ CAT19 CAR T cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    cartm_cd4_bm = list(
      analyte = "memory CD4+ CAT19 CAR T cells",
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
    disease_state = "Relapsed / refractory paediatric B-cell acute lymphoblastic leukaemia (CARPALL study, Ghorashian et al. 2019)",
    dose_range = "Single infusion of a median 1 x 10^6 CAR+ T cells/kg",
    regions = "United Kingdom",
    notes = paste(
      "Individual flow cytometry (n = 13; LLOQ 1 cell/uL) and qPCR (n = 14;",
      "LLOQ 100 copies/ug gDNA) CAR T-cell profiles in blood and bone marrow,",
      "digitised from Ghorashian 2019; censored values handled with the M3",
      "method. The product contained a median 54.65% CD4+ and 37.65% CD8+ CAR",
      "T cells (Table 1)."
    )
  )

  ini({
    # ---- Estimated (Table 2, CAT19; per-day values from the text) ----
    lkexp_max_cd8 <- log(5.3)
    label("Maximum first-order expansion rate constant of CD8+ CAR T cells, KmaxExpCD8 (1/day)") # Table 2 CAT19: KmaxExpCD8 = 159 /month (17.6% RSE); Results: 5.3 day-1
    lrm <- log(0.000037)
    label("Net first-order effector-to-memory conversion rate constant, Rm (1/day)") # Table 2 CAT19: Rm = 0.001 /month (81% RSE); Results: 0.000037 day-1
    lkel_m <- log(0.013)
    label("First-order elimination rate constant of memory CAR T cells from blood, Kelm (1/day)") # Table 2 CAT19: Kelm = 0.39 /month (29.5% RSE); Results: 0.013 day-1
    lk12 <- log(0.11)
    label("Blood-to-bone-marrow distribution rate constant, K12 (1/day)") # Table 2 CAT19: K12 = 3.3 /month (44.6% RSE); Table 2 step 3: K12 = 0.11 /day (from CAT19 estimates)
    int_transgene <- 7.26
    label("Intercept of ln(transgene copies/ug gDNA) on ln(CAR T cells/uL), alpha") # Table 2 CAT19: alpha = 7.26 (10.7% RSE)
    slope_transgene <- 0.74
    label("Slope of ln(transgene copies/ug gDNA) on ln(CAR T cells/uL), beta (unitless)") # Table 2 CAT19: beta = 0.74 (6.9% RSE)

    # ---- Fixed (Table 2, CAT19) ----
    lk21 <- fixed(log(0.176))
    label("Bone-marrow-to-blood redistribution rate constant, K21 (1/day)") # Table 2 CAT19: K21 = 0.176 /day (fixed, ref 19 Singh 2021)
    kexp_ratio_cd8_cd4 <- fixed(1.12)
    label("Ratio of the CD8+ to CD4+ maximum expansion rate constants (unitless)") # Table 2 CAT19: KmaxExpCD4 = KmaxExpCD8/1.12 (fixed from Liso-cel estimates)
    lpct_cd8 <- fixed(log(38))
    label("Percentage of CD8+ CAR T cells in the infused product (%)") # Table 2 CAT19: PCTCD8 = 38 (fixed, ref 20); Table 1: 37.65% median
    lpct_cd4 <- fixed(log(54))
    label("Percentage of CD4+ CAR T cells in the infused product (%)") # Table 2 CAT19: PCTCD4 = 54 (fixed, ref 20); Table 1: 54.65% median
    lec50_exp <- fixed(log(10))
    label("CAR-target complexes per CAR T cell giving half-maximal expansion, EC50Exp (complexes/cell)") # Table 2 CAT19: EC50Exp = 10 (fixed, ref 19 Singh 2021)
    lkel_e <- fixed(log(113))
    label("First-order elimination rate constant of effector CAR T cells from blood, Kele (1/day)") # Table 2 CAT19: Kele = 113 /day (fixed, ref 19 Singh 2021)
    lkkill_max <- fixed(log(1.8))
    label("Maximum first-order tumour-cell killing rate constant, KKillMax (1/day)") # Table 2 CAT19: KKillMax = 1.8 (fixed from Liso-cel estimates)
    lkc50_kill <- fixed(log(0.28))
    label("CAR-target complexes per tumour cell giving half-maximal killing, KC50CAR-T (complexes/cell)") # Table 2 CAT19: KC50CAR-T = 0.28 (fixed, CAT19 in vitro estimate)
    kon <- fixed(2.15e5)
    label("CAR-CD19 association rate constant, Kon (1/M/s)") # Table 2 CAT19: Kon = 2.15E5 (fixed, ref 20)
    koff <- fixed(0.003)
    label("CAR-CD19 dissociation rate constant, Koff (1/s)") # Table 2 CAT19: Koff = 0.003 (fixed, ref 20)
    ag_car <- fixed(15000)
    label("CAR density on CAR T cells, AgCAR (CARs/cell)") # Table 2 CAT19: AgCAR = 15,000 (fixed, refs 19, 41)
    ag_tumor <- fixed(4000)
    label("CD19 density on leukaemic cells, AgTumor (antigens/cell)") # Table 2 CAT19: AgTumor = 4000 (fixed, ref 39)
    dt_tumor <- fixed(10.7)
    label("Paediatric ALL tumour doubling time, DTTumor (day)") # Table 2 CAT19: DTTumor = 10.7 (fixed, ref 49)
    kmig <- fixed(0.001)
    label("Bone-marrow-to-blood migration rate constant of tumour cells, Kmigration (1/day)") # Table 2 CAT19: Kmigration = 0.001 (fixed, ref 44)
    tumor_bm0 <- fixed(1.06e11)
    label("Initial tumour burden in bone marrow, ICTumorBM (cells/L)") # Table 2 CAT19: ICTumorBM = 1.06E11 (fixed, refs 45, 47)
    tumor_pb0 <- fixed(1.8e9)
    label("Initial tumour burden in blood, ICTumorPB (cells/L)") # Table 2 CAT19: ICTumorPB = 1.8E9 (fixed, ref 50)
    v_blood <- fixed(5)
    label("Blood compartment volume for a 70 kg adult, Vb (L)") # Table 2 CAT19: Vb = 5 * 30/70 (fixed, ref 40)
    v_bonemarrow <- fixed(1.6)
    label("Bone marrow compartment volume for a 70 kg adult, Vbm (L)") # Table 2 CAT19: Vbm = 1.6 * 30/70 (fixed, ref 40)

    # ---- Between-subject variability ----
    # Monolix omega is the SD of a log-normal random effect; variance = omega^2.
    etalkexp_max_cd8 ~ 0.16 # Table 2 CAT19: omega KmaxExpCD8 = 0.4 (18.3% RSE); 0.4^2
    etalpct_cd8 ~ fixed(0.25) # Table 2 CAT19: omega PCTCD8 = 0.5 (not estimated); 0.5^2
    etalpct_cd4 ~ fixed(0.25) # Table 2 CAT19: omega PCTCD4 = 0.5 (not estimated); 0.5^2
    etalrm ~ 25.4016 # Table 2 CAT19: omega Rm = 5.04 (16.6% RSE); 5.04^2
    etalkel_m ~ 1.8225 # Table 2 CAT19: omega Kelm = 1.35 (16.4% RSE); 1.35^2
    etalk12 ~ 4.3264 # Table 2 CAT19: omega K12 = 2.08 (15.5% RSE); 2.08^2
    etalk21 ~ 1.1664 # Table 2 CAT19: omega K21 = 1.08 (17.8% RSE); 1.08^2
  })

  model({
    # ---- Individual parameters ----
    kexp_max_cd8 <- exp(lkexp_max_cd8 + etalkexp_max_cd8)
    kexp_max_cd4 <- kexp_max_cd8 / kexp_ratio_cd8_cd4
    rm <- exp(lrm + etalrm)
    kel_m <- exp(lkel_m + etalkel_m)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    pct_cd8 <- exp(lpct_cd8 + etalpct_cd8)
    pct_cd4 <- exp(lpct_cd4 + etalpct_cd4)
    ec50_exp <- exp(lec50_exp)
    kel_e <- exp(lkel_e)
    kkill_max <- exp(lkkill_max)
    kc50_kill <- exp(lkc50_kill)
    kg_tumor <- log(2) / dt_tumor

    # ---- Weight-scaled volumes (Appendix S2 step-3 code) ----
    vb <- v_blood * (WT / 70)
    vbm <- v_bonemarrow * (WT / 70)

    # ---- Unit conversion of the binding constants (Appendix S2 code) ----
    kon_d <- kon * (60 * 60 * 24) / 6.023e23
    koff_d <- koff * (60 * 60 * 24)

    # ---- CAR T-cell concentrations (cells/L) ----
    cart_cd8_pb_l <- (carte_cd8_pb + cartm_cd8_pb) / vb
    cart_cd4_pb_l <- (carte_cd4_pb + cartm_cd4_pb) / vb
    cart_cd8_bm_l <- (carte_cd8_bm + cartm_cd8_bm) / vbm
    cart_cd4_bm_l <- (carte_cd4_bm + cartm_cd4_bm) / vbm

    # ---- Blood binding and drivers (Appendix S1 Eqs. 19-29) ----
    car_free_cd8_pb <- cart_cd8_pb_l * ag_car - complex_cd8_pb
    car_free_cd4_pb <- cart_cd4_pb_l * ag_car - complex_cd4_pb
    ag_free_pb <- tumor_pb * ag_tumor - (complex_cd8_pb + complex_cd4_pb)
    cplx_per_cart_cd8_pb <- complex_cd8_pb / cart_cd8_pb_l
    cplx_per_cart_cd4_pb <- complex_cd4_pb / cart_cd4_pb_l
    kexp_cd8_pb <- kexp_max_cd8 * cplx_per_cart_cd8_pb / (ec50_exp + cplx_per_cart_cd8_pb)
    kexp_cd4_pb <- kexp_max_cd4 * cplx_per_cart_cd4_pb / (ec50_exp + cplx_per_cart_cd4_pb)
    # Eq. 29 x Tumor^PB, written so that an eradicated tumour gives zero kill
    # flux instead of 0/0.
    cplx_pb <- complex_cd8_pb + complex_cd4_pb
    kill_flux_pb <- kkill_max * cplx_pb * tumor_pb / (kc50_kill * tumor_pb + cplx_pb)

    # ---- Bone-marrow binding and drivers (Appendix S1 Eqs. 31-41) ----
    car_free_cd8_bm <- cart_cd8_bm_l * ag_car - complex_cd8_bm
    car_free_cd4_bm <- cart_cd4_bm_l * ag_car - complex_cd4_bm
    ag_free_bm <- tumor_bm * ag_tumor - (complex_cd8_bm + complex_cd4_bm)
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

    # ---- Dose split with 50% variability (Appendix S1, step 2 CAT19) ----
    f(carte_cd8_pb) <- pct_cd8 / 100
    f(carte_cd4_pb) <- pct_cd4 / 100

    # ---- Initial conditions (Appendix S2 code) ----
    cartm_cd8_pb(0) <- 1e-12 * vb
    cartm_cd4_pb(0) <- 1e-12 * vb
    carte_cd8_bm(0) <- 1e-12 * vbm
    cartm_cd8_bm(0) <- 1e-12 * vbm
    carte_cd4_bm(0) <- 1e-12 * vbm
    cartm_cd4_bm(0) <- 1e-12 * vbm
    complex_cd8_pb(0) <- 1e-12
    complex_cd4_pb(0) <- 1e-12
    complex_cd8_bm(0) <- 1e-12
    complex_cd4_bm(0) <- 1e-12
    tumor_pb(0) <- tumor_pb0
    tumor_bm(0) <- tumor_bm0

    # ---- Outputs (Eq. 1; 1E-6 converts cells/L to cells/uL) ----
    cart_pb <- (cart_cd8_pb_l + cart_cd4_pb_l) * 1e-6
    cart_bm <- (cart_cd8_bm_l + cart_cd4_bm_l) * 1e-6
    transgene <- exp(int_transgene + slope_transgene * log(cart_pb))
    transgene_bm <- exp(int_transgene + slope_transgene * log(cart_bm))
  })
}
