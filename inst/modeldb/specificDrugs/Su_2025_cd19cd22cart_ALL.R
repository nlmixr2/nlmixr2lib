# Multiscale mechanistic cellular-kinetic / pharmacodynamic (CK-PD) model for
# the dual-targeted bicistronic CD19/CD22 CAR T-cell product AUTO3 in pediatric
# and young-adult relapsed/refractory B-cell acute lymphoblastic leukemia
# (B-ALL), published by Su et al. 2025 (CPT Pharmacometrics Syst Pharmacol
# 14:229-245). This is the paper's clinical centrepiece (Figure 1b), fit to the
# phase 1 dose-escalation cellular-kinetic data of Cordoba et al. 2021 by
# simultaneously describing qPCR transgene copies and flow-cytometry CAR-T cell
# counts. The paper's three in-vitro cell-level killing sub-models (Figure 1a)
# are shipped as the companion files Su_2025_cd19cd22cart_invitro.R,
# Su_2025_gprc5dbcmacart_41BB_invitro.R and Su_2025_gprc5dbcmacart_CD28_invitro.R.
#
# SOURCES. Equations: Supplement 1 (Supplementary Text) Eqs. 15-52, and the
# as-run Monolix model code in Supplement 2 ("Clinical CKPD Model Accounting for
# Phenotypic Transformation (B-cell ALL)"). Parameter values: main-text Table 1,
# "Parameters associated with clinical CK-pd model accounting for phenotypic
# transformation (B-cell ALL)". The two residual-error magnitudes are not
# reported in the paper or its supplement; they were recovered by the
# maintainers from the Figure S2 goodness-of-fit panels (see the vignette
# Errata).
#
# STATE UNITS. The four CAR-T pools and the two tumour pools are carried as
# concentrations (cells/L of blood or of bone marrow); complexes are numbers/L.
# A CAR-T dose record is the number of infused CAR+ T cells per kg, scaled to a
# blood concentration by the bioavailability factor wt/Vb (which equals 70/5 = 14
# exactly, because Vb = 5*wt/70). The infused product is split between the memory
# and effector pools by the product's phenotypic composition (default 0.87:0.13,
# the cohort mean); see the vignette for per-patient dosing.

Su_2025_cd19cd22cart_ALL <- function() {
  description <- "QSP (multiscale mechanistic cellular-kinetic / pharmacodynamic model). Dual-targeted bicistronic CD19/CD22 CAR T-cell therapy (AUTO3) in pediatric and young-adult relapsed/refractory B-cell acute lymphoblastic leukemia. Effector and memory CAR T cells distribute between peripheral blood and bone marrow (first-order K12 / K21) and are eliminated from blood at phenotype-specific rates; effector cells convert to memory at a net first-order rate Rm. In each compartment the two CARs bind CD19 and CD22 on B-ALL cells (second-order kon / first-order koff), forming single-target and double-target complexes; complexes per CAR T cell drive Emax expansion of the effector pool and complexes per tumour cell drive Emax killing of a logistically growing tumour that migrates from bone marrow to blood. Outputs: blood CAR-T cells per microlitre (flow cytometry) and qPCR transgene copies per ug genomic DNA, linked by a fixed log-log regression. Baseline CD19/CD22 antigen densities are patient-level inputs encoded here at the cohort averages and overridable per patient. Between-subject variability was estimated on the expansion, elimination, distribution and killing rate constants."
  reference <- paste(
    "Su MC, Dey A, Maddah E, Mugundu GM, Singh AP (2025).",
    "Quantitative pharmacology of dual-targeted bicistronic CAR-T-cell therapy",
    "using multiscale mechanistic modeling.",
    "CPT Pharmacometrics Syst Pharmacol 14(2):229-245.",
    "doi:10.1002/psp4.13259.",
    "Equations and Monolix model code in the Supporting Information",
    "(PMC11812944 open-access package).",
    "Fit to the phase 1 CD19/CD22 CAR-T (AUTO3) cellular-kinetic data of",
    "Cordoba et al. 2021 (Nat Med 27:1797-1805; doi:10.1038/s41591-021-01497-1).",
    sep = " "
  )
  vignette <- "Su_2025_dual_targeted_cart"

  # Effector / memory CAR T-cell pools in blood and bone marrow, the six
  # antigen-specific CAR-target complexes, and the two tumour pools are
  # paper-mechanistic states with no canonical analogue; carte_pb / cartm_pb
  # follow the same pools in Singh_2021_idecabtageneVicleucel_human.R and the
  # Hardiansyah_2019 CAR-T models, and the complex_<antigen>_<site> / tumor_<site>
  # forms extend the registered complex / tumor families by antigen and site.
  paper_specific_compartments <- c(
    "carte_pb",
    "cartm_pb",
    "carte_bm",
    "cartm_bm",
    "complex_cd19_pb",
    "complex_cd22_pb",
    "complex_cd19_cd22_pb",
    "complex_cd19_bm",
    "complex_cd22_bm",
    "complex_cd19_cd22_bm",
    "tumor_pb",
    "tumor_bm"
  )

  units <- list(
    time = "day",
    dosing = "cells (CAR+ T cells per kg; input as amt on cartm_pb and carte_pb)",
    concentration = "cells/uL (cells output); transgene copies/ug genomic DNA (transgene output)"
  )
  dosing <- c("cartm_pb", "carte_pb")

  covariateData <- list(
    AGE = list(
      description = "Age; sets the age-based body weight (wt = 2*AGE + 8 for AGE <= 5, wt = 3*AGE + 7 for AGE >= 6) that scales the blood and bone-marrow volumes",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Age-based weight from Carasco et al. 2016; the dosing scale wt/Vb cancels to 70/5 since Vb = 5*wt/70.",
      source_name = "Age"
    )
  )

  compartmentData <- list(
    carte_pb = list(
      analyte = "effector CD19/CD22 CAR T cells (AUTO3)",
      units = "cells/L",
      specimen = "whole blood",
      verified = TRUE
    ),
    cartm_pb = list(
      analyte = "memory CD19/CD22 CAR T cells (AUTO3)",
      units = "cells/L",
      specimen = "whole blood",
      verified = TRUE
    ),
    carte_bm = list(
      analyte = "effector CD19/CD22 CAR T cells (AUTO3)",
      units = "cells/L",
      specimen = "tissue",
      verified = TRUE
    ),
    cartm_bm = list(
      analyte = "memory CD19/CD22 CAR T cells (AUTO3)",
      units = "cells/L",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_cd19_pb = list(
      analyte = "CAR-CD19 single-target complexes",
      units = "complexes/L",
      specimen = "whole blood",
      verified = TRUE
    ),
    complex_cd22_pb = list(
      analyte = "CAR-CD22 single-target complexes",
      units = "complexes/L",
      specimen = "whole blood",
      verified = TRUE
    ),
    complex_cd19_cd22_pb = list(
      analyte = "CAR-CD19-CD22 double-target complexes",
      units = "complexes/L",
      specimen = "whole blood",
      verified = TRUE
    ),
    complex_cd19_bm = list(
      analyte = "CAR-CD19 single-target complexes",
      units = "complexes/L",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_cd22_bm = list(
      analyte = "CAR-CD22 single-target complexes",
      units = "complexes/L",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_cd19_cd22_bm = list(
      analyte = "CAR-CD19-CD22 double-target complexes",
      units = "complexes/L",
      specimen = "tissue",
      verified = TRUE
    ),
    tumor_pb = list(analyte = "CD19+/CD22+ B-ALL cells", units = "cells/L", specimen = "whole blood", verified = TRUE),
    tumor_bm = list(analyte = "CD19+/CD22+ B-ALL cells", units = "cells/L", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 14L,
    n_studies = 1L,
    age_range = "4-16 years",
    disease_state = "Relapsed or refractory pediatric and young-adult B-cell acute lymphoblastic leukemia (phase 1 dose-escalation study of CD19/CD22 CAR-T AUTO3; Cordoba et al. 2021)",
    dose_range = "Single intravenous infusion of 0.3-2 x 10^6, 3 x 10^6 or 5 x 10^6 CAR+ T cells/kg",
    sex_female_pct = 29,
    regions = "United Kingdom (AUTO3 phase 1; Cordoba et al. 2021)",
    notes = paste(
      "14 patients contributed 157 qPCR samples; 5 mapped patients contributed",
      "43 flow-cytometry samples. The model accounts for the memory/effector",
      "phenotypic composition of each infused product (cohort mean 0.86:0.13,",
      "range 0.17-0.99), individual baseline CD19 and CD22 tumour-antigen",
      "densities, and age-based body weight. A 1:1 ratio of CD4+ to CD8+ CAR T",
      "cells was assumed. Missing individual inputs were set to the cohort mean."
    )
  )

  ini({
    # ---- Estimated parameters (Table 1, clinical CK-PD model) ----
    lkmax_exp <- log(2.04)
    label("Maximum first-order effector CAR-T expansion rate constant, KMaxExp (1/day)") # Table 1 clinical: KMaxExp = 2.04 (5.21% RSE)
    lec50_exp <- log(42.1)
    label("CAR-target complexes per effector cell giving half-maximal expansion, EC50Exp (complexes/cell)") # Table 1 clinical: EC50Exp = 42.1 (34.0% RSE)
    lkel_m <- log(0.015)
    label("First-order elimination rate constant of memory CAR T cells from blood, Kelm (1/day)") # Table 1 clinical: Kelm = 0.015 (31.7% RSE)
    lk12 <- log(2.91)
    label("Blood-to-bone-marrow distribution rate constant, K12 (1/day)") # Table 1 clinical: K12 = 2.91 (47.2% RSE)
    lkmax_kill <- log(0.58)
    label("Maximum first-order tumour-cell killing rate constant, KMaxKill (1/day)") # Table 1 clinical: KMaxKill = 0.58 (9.90% RSE)
    lkcarry <- log(0.39)
    label("Tumour carrying capacity, KcarryTumor (10^10 cells/L)") # Table 1 clinical: KcarryTumor = 0.39 (5.74% RSE)

    # ---- Fixed structural / system parameters (Table 1, clinical) ----
    lgamma_exp <- fixed(log(1))
    label("Sigmoidicity factor for CAR-T expansion, gamma_exp (unitless)") # Table 1 clinical: gammaExpCAR-T = 1 (fixed)
    lgamma_kill <- fixed(log(1))
    label("Sigmoidicity factor for tumour killing, gamma_kill (unitless)") # in-vitro composite gammaKillTumor = 1 (fixed); carried into the clinical model
    lrm <- fixed(log(0.00002))
    label("Net first-order effector-to-memory conversion rate constant, Rm (1/day)") # Table 1 clinical: Rm = 0.00002 (fixed)
    lkel_e <- fixed(log(113))
    label("First-order elimination rate constant of effector CAR T cells from blood, Kele (1/day)") # Table 1 clinical: Kele = 113 (fixed)
    lk21 <- fixed(log(0.18))
    label("Bone-marrow-to-blood redistribution rate constant, K21 (1/day)") # Table 1 clinical: K21 = 0.18 (fixed)
    lkc50_kill <- fixed(log(3.57))
    label("CAR-target complexes per tumour cell giving half-maximal killing, KC50Kill (complexes/cell)") # Table 1 clinical: KC50Kill = 3.57 (fixed, in vitro estimate)
    ldt_tumor <- fixed(log(10.7))
    label("Tumour-cell doubling time, DTTumor (day)") # Table 1 clinical: DTTumor = 10.7
    lkmigration <- fixed(log(0.001))
    label("Net first-order tumour-cell migration rate constant, bone marrow to blood, Kmigration (1/day)") # Table 1 clinical: Kmigration = 0.001 (fixed)
    dens_cd19_car <- fixed(15000)
    label("Anti-CD19 CAR density on CAR T cells (CARs/cell)") # Table 1 clinical: DensityCD19CAR = 15,000 (fixed)
    dens_cd22_car <- fixed(15000)
    label("Anti-CD22 CAR density on CAR T cells (CARs/cell)") # Table 1 clinical: DensityCD22CAR = 15,000 (fixed)
    kon_orig_cd19 <- fixed(1e5)
    label("Anti-CD19 CAR association rate constant, konCD19CAR (1/M/s)") # Table 1 clinical: konCD19CAR = 1e5 (fixed)
    kon_orig_cd22 <- fixed(1e5)
    label("Anti-CD22 CAR association rate constant, konCD22CAR (1/M/s)") # Table 1 clinical: konCD22CAR = 1e5 (fixed)
    kd_cd19 <- fixed(1.09)
    label("Anti-CD19 CAR equilibrium dissociation constant, KdCD19CAR (nM)") # Table 1 clinical: KdCD19CAR = 1.09 (fixed)
    kd_cd22 <- fixed(32)
    label("Anti-CD22 CAR equilibrium dissociation constant, KdCD22CAR (nM)") # Table 1 clinical: KdCD22CAR = 32 (fixed)
    alpha <- fixed(4.45)
    label("Intercept of the log qPCR vs log flow regression, alpha (log copies/ug gDNA)") # Table 1 clinical: alpha = 4.45 (fixed)
    beta <- fixed(1.41)
    label("Slope of the log qPCR vs log flow regression, beta (unitless)") # Table 1 clinical: beta = 1.41 (fixed)
    dens_cd19_tumor <- fixed(16666)
    label("Baseline CD19 tumour-antigen density (antigens/cell); patient-level input, set to the cohort average") # Table 1 clinical: DensityCD19Tumor average = 16,666 (1324-39,098)
    dens_cd22_tumor <- fixed(4269)
    label("Baseline CD22 tumour-antigen density (antigens/cell); patient-level input, set to the cohort average") # Table 1 clinical: DensityCD22Tumor average = 4269 (181-10,137)
    tumor0_pb <- fixed(1.8e9)
    label("Initial B-ALL cell concentration in peripheral blood (cells/L)") # Table 1 clinical: ICTumorPB = 1.8e9
    tumor0_bm <- fixed(1.6e11)
    label("Initial B-ALL cell concentration in bone marrow (cells/L)") # Table 1 clinical: ICTumorBM = 1.6e11

    # ---- Between-subject variability ----
    # Table 1's omega column is the Monolix standard deviation of the log-normal
    # random effect (SAEM, Monolix 2021R); variances below are omega^2.
    etalkmax_exp ~ 0.0169 # Table 1 clinical: omega KMaxExp = 0.13 (32.6% RSE); 0.13^2
    etalkel_m ~ 0.4489 # Table 1 clinical: omega Kelm = 0.67 (31.4% RSE); 0.67^2
    etalk12 ~ 2.4964 # Table 1 clinical: omega K12 = 1.58 (21.8% RSE); 1.58^2
    etalkmax_kill ~ 0.0961 # Table 1 clinical: omega KMaxKill = 0.31 (22.0% RSE); 0.31^2

    # ---- Residual error (recovered from Figure S2 goodness-of-fit panels) ----
    # Not reported in the paper or its supplement; the maintainers read the 90%
    # prediction-interval offsets off the Monolix observations-vs-predictions
    # panels (S2a/S2d for qPCR, S2b/S2f for flow). See the vignette Errata.
    propSd_transgene <- 0.726
    label("Proportional residual error on qPCR transgene copies (fraction)") # figure-derived from Figure S2a/S2d censored-point offsets
    addSd_cells <- 2.09
    label("Additive residual error on flow CAR-T cells (cells/uL)") # figure-derived from Figure S2b/S2f (combined1 additive term)
    propSd_cells <- 0.536
    label("Proportional residual error on flow CAR-T cells (fraction)") # figure-derived from Figure S2b/S2f (combined1 proportional term)
  })

  model({
    # ---- Individual parameters ----
    kmax_exp <- exp(lkmax_exp + etalkmax_exp)
    ec50_exp <- exp(lec50_exp)
    kel_m <- exp(lkel_m + etalkel_m)
    k12 <- exp(lk12 + etalk12)
    kmax_kill <- exp(lkmax_kill + etalkmax_kill)
    kcarry <- exp(lkcarry)
    gamma_exp <- exp(lgamma_exp)
    gamma_kill <- exp(lgamma_kill)
    rm <- exp(lrm)
    kel_e <- exp(lkel_e)
    k21 <- exp(lk21)
    kc50_kill <- exp(lkc50_kill)
    dt_tumor <- exp(ldt_tumor)
    kmigration <- exp(lkmigration)

    # ---- Age-based weight and compartment volumes (Monolix code) ----
    wt <- ifelse(AGE <= 5, 2 * AGE + 8, 3 * AGE + 7)
    vb <- 5 * wt / 70 # blood volume (L), Table 1: Vb = 5*(age-based wt/70)
    vbm <- 1.6 * wt / 70 # bone-marrow volume (L), Table 1: Vbm = 1.6*(age-based wt/70)

    # ---- Unit conversion of the binding constants (Supplement eq. notes) ----
    # kon: 1/M/s -> 1/(#/L)/day (divide by Avogadro, x seconds/day);
    # koff = kon*Kd: 1/s -> 1/day.
    kon_cd19 <- kon_orig_cd19 * (60 * 60 * 24) / 6.023e23
    koff_cd19 <- kon_orig_cd19 * (kd_cd19 * 1e-9) * (60 * 60 * 24)
    kon_cd22 <- kon_orig_cd22 * (60 * 60 * 24) / 6.023e23
    koff_cd22 <- kon_orig_cd22 * (kd_cd22 * 1e-9) * (60 * 60 * 24)

    # ---- Target engagement in peripheral blood (Supplement Eqs. 21-31) ----
    ag_car19_pb <- dens_cd19_car * (carte_pb + cartm_pb)
    ag_car22_pb <- dens_cd22_car * (carte_pb + cartm_pb)
    ag_tum19_pb <- dens_cd19_tumor * tumor_pb
    ag_tum22_pb <- dens_cd22_tumor * tumor_pb
    fcar19_pb <- ag_car19_pb - complex_cd19_pb - complex_cd19_cd22_pb
    fcar22_pb <- ag_car22_pb - complex_cd22_pb - complex_cd19_cd22_pb
    ftum19_pb <- ag_tum19_pb - complex_cd19_pb - complex_cd19_cd22_pb
    ftum22_pb <- ag_tum22_pb - complex_cd22_pb - complex_cd19_cd22_pb
    d/dt(complex_cd19_pb) <- kon_cd19 * fcar19_pb * ftum19_pb - koff_cd19 * complex_cd19_pb +
      koff_cd22 * complex_cd19_cd22_pb - kon_cd22 * complex_cd19_pb * ftum22_pb
    d/dt(complex_cd22_pb) <- kon_cd22 * fcar22_pb * ftum22_pb - koff_cd22 * complex_cd22_pb +
      koff_cd19 * complex_cd19_cd22_pb - kon_cd19 * complex_cd22_pb * ftum19_pb
    d/dt(complex_cd19_cd22_pb) <- kon_cd22 * complex_cd19_pb * ftum22_pb + kon_cd19 * complex_cd22_pb * ftum19_pb -
      koff_cd22 * complex_cd19_cd22_pb - koff_cd19 * complex_cd19_cd22_pb

    # ---- Target engagement in bone marrow (Supplement Eqs. 37-47) ----
    ag_car19_bm <- dens_cd19_car * (carte_bm + cartm_bm)
    ag_car22_bm <- dens_cd22_car * (carte_bm + cartm_bm)
    ag_tum19_bm <- dens_cd19_tumor * tumor_bm
    ag_tum22_bm <- dens_cd22_tumor * tumor_bm
    fcar19_bm <- ag_car19_bm - complex_cd19_bm - complex_cd19_cd22_bm
    fcar22_bm <- ag_car22_bm - complex_cd22_bm - complex_cd19_cd22_bm
    ftum19_bm <- ag_tum19_bm - complex_cd19_bm - complex_cd19_cd22_bm
    ftum22_bm <- ag_tum22_bm - complex_cd22_bm - complex_cd19_cd22_bm
    d/dt(complex_cd19_bm) <- kon_cd19 * fcar19_bm * ftum19_bm - koff_cd19 * complex_cd19_bm +
      koff_cd22 * complex_cd19_cd22_bm - kon_cd22 * complex_cd19_bm * ftum22_bm
    d/dt(complex_cd22_bm) <- kon_cd22 * fcar22_bm * ftum22_bm - koff_cd22 * complex_cd22_bm +
      koff_cd19 * complex_cd19_cd22_bm - kon_cd19 * complex_cd22_bm * ftum19_bm
    d/dt(complex_cd19_cd22_bm) <- kon_cd22 * complex_cd19_bm * ftum22_bm + kon_cd19 * complex_cd22_bm * ftum19_bm -
      koff_cd22 * complex_cd19_cd22_bm - koff_cd19 * complex_cd19_cd22_bm

    # ---- CAR-T cellular kinetics (Supplement Eqs. 15, 16, 19, 20) ----
    totcart_pb <- carte_pb + cartm_pb
    cplxcart_pb <- (complex_cd19_pb + complex_cd22_pb + complex_cd19_cd22_pb) / totcart_pb
    kexp_pb <- kmax_exp * cplxcart_pb^gamma_exp / (ec50_exp^gamma_exp + cplxcart_pb^gamma_exp)
    d/dt(carte_pb) <- (-k12 * carte_pb * vb + k21 * carte_bm * vbm) / vb - kel_e * carte_pb +
      kexp_pb * carte_pb - rm * carte_pb
    d/dt(cartm_pb) <- (-k12 * cartm_pb * vb + k21 * cartm_bm * vbm) / vb - kel_m * cartm_pb +
      rm * carte_pb

    totcart_bm <- carte_bm + cartm_bm
    cplxcart_bm <- (complex_cd19_bm + complex_cd22_bm + complex_cd19_cd22_bm) / totcart_bm
    kexp_bm <- kmax_exp * cplxcart_bm^gamma_exp / (ec50_exp^gamma_exp + cplxcart_bm^gamma_exp)
    d/dt(carte_bm) <- (k12 * carte_pb * vb - k21 * carte_bm * vbm) / vbm + kexp_bm * carte_bm -
      rm * carte_bm
    d/dt(cartm_bm) <- (k12 * cartm_pb * vb - k21 * cartm_bm * vbm) / vbm + rm * carte_bm

    # ---- Tumour dynamics (Supplement Eqs. 36, 52) ----
    cplxtum_pb <- (complex_cd19_pb + complex_cd22_pb + complex_cd19_cd22_pb) / tumor_pb
    kkill_pb <- kmax_kill * cplxtum_pb^gamma_kill / (kc50_kill^gamma_kill + cplxtum_pb^gamma_kill)
    d/dt(tumor_pb) <- log(2) / dt_tumor * tumor_pb * (1 - tumor_pb / (kcarry * 1e10)) -
      kkill_pb * tumor_pb + kmigration * tumor_bm * vbm / vb
    cplxtum_bm <- (complex_cd19_bm + complex_cd22_bm + complex_cd19_cd22_bm) / tumor_bm
    kkill_bm <- kmax_kill * cplxtum_bm^gamma_kill / (kc50_kill^gamma_kill + cplxtum_bm^gamma_kill)
    d/dt(tumor_bm) <- log(2) / dt_tumor * tumor_bm * (1 - tumor_bm / (kcarry * 1e10)) -
      kkill_bm * tumor_bm - kmigration * tumor_bm * vbm / vbm

    # ---- Initial conditions (Monolix code; small non-zero CAR-T / complex
    # pools keep the complexes-per-cell ratios finite before infusion) ----
    carte_pb(0) <- 1e-6
    cartm_pb(0) <- 1e-6
    carte_bm(0) <- 1e-6
    cartm_bm(0) <- 1e-6
    complex_cd19_pb(0) <- 1e-6
    complex_cd22_pb(0) <- 1e-6
    complex_cd19_cd22_pb(0) <- 1e-6
    complex_cd19_bm(0) <- 1e-6
    complex_cd22_bm(0) <- 1e-6
    complex_cd19_cd22_bm(0) <- 1e-6
    tumor_pb(0) <- tumor0_pb
    tumor_bm(0) <- tumor0_bm

    # ---- Dosing scale: cells/kg infused -> cells/L of blood (p = wt/Vb) ----
    f(cartm_pb) <- wt / vb
    f(carte_pb) <- wt / vb

    # ---- Outputs (Supplement Eq. 18; Monolix code TotCART_PB_uL / transgene) ----
    cells <- totcart_pb * 1e-6 # flow-cytometry CAR-T cells per uL of blood
    transgene <- exp(alpha + beta * log(cells)) # qPCR transgene copies per ug gDNA
    cells ~ add(addSd_cells) + prop(propSd_cells) + combined1()
    transgene ~ prop(propSd_transgene)
  })
}
