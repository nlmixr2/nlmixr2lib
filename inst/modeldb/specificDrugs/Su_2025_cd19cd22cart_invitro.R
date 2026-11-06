# In-vitro cell-level pharmacodynamic (tumour-killing) model for single- and
# dual-targeted CD19 and/or CD22 CAR T cells co-cultured with B-cell ALL cell
# lines, published by Su et al. 2025 (CPT Pharmacometrics Syst Pharmacol
# 14:229-245, Figure 1a / Figure 2). CAR T cells and tumour cells coexist in a
# well-mixed culture well; the two CARs bind CD19 and CD22 to form single- and
# double-target complexes, whose number per tumour cell drives Emax killing.
# This is the composite dual-targeting killing model (Figure 2) whose estimated
# potency (KC50 = 3.57 complexes/tumour cell) is carried forward as the fixed
# in-vivo killing potency of the clinical model Su_2025_cd19cd22cart_ALL.R.
#
# SOURCES. Equations: Supplement 1 (Supplementary Text) Eqs. 1-13, and the
# as-run Monolix model code in Supplement 2 ("Cell-Level Pharmacodynamic Model
# (B-cell ALL)"). Estimated values (KMaxKill, KC50, with IIV) and the fixed
# system constants: main-text Table 1, "Parameters associated with cell-level
# pharmacodynamic model (B-cell ALL)". No residual-error magnitude is reported;
# it is encoded as fixed(0) (see the vignette Errata). The cell-line antigen
# densities, CAR affinities and doubling time are per-scenario regressors; the
# defaults here are AUTO3 (dual CD19/CD22) co-cultured with RajiWT, and the
# vignette overrides them to reproduce every panel of Figure 2.
#
# STATE UNITS. n_tumor and n_cart are cell NUMBERS in a well of MV = 1e-4 L
# (100 uL); the three complexes are numbers/L. Viability is the percent of the
# initial tumour cells remaining, normalised to the non-transduced control.

Su_2025_cd19cd22cart_invitro <- function() {
  description <- "In vitro (B-cell ALL cell lines). Cell-level tumour-killing model for single- and dual-targeted CD19/CD22 CAR T cells (AUTO3, FMC63, CD19 CAR, CD22 CAR) co-cultured with leukemia cell lines. The two CARs bind CD19 and CD22 (second-order kon / first-order koff) to form single-target and double-target complexes; complexes per tumour cell drive Emax killing of the tumour-cell pool while CAR T cells grow exponentially. Output: percent tumour-cell viability relative to the non-transduced control. Estimated composite potency KC50 = 3.57 complexes/tumour cell is reused as the fixed in-vivo killing potency of the clinical CD19/CD22 CAR-T model. Cell-line antigen densities, CAR affinities and doubling time are regressors; defaults are AUTO3 on RajiWT."
  reference <- paste(
    "Su MC, Dey A, Maddah E, Mugundu GM, Singh AP (2025).",
    "Quantitative pharmacology of dual-targeted bicistronic CAR-T-cell therapy",
    "using multiscale mechanistic modeling.",
    "CPT Pharmacometrics Syst Pharmacol 14(2):229-245.",
    "doi:10.1002/psp4.13259.",
    "Equations and Monolix model code in the Supporting Information",
    "(PMC11812944 open-access package).",
    "In-vitro cytotoxicity data digitized from Cordoba et al. 2021",
    "(Nat Med 27:1797-1805).",
    sep = " "
  )
  vignette <- "Su_2025_dual_targeted_cart"

  # Tumour-cell and CAR-T-cell well populations and the three antigen-specific
  # CAR-target complexes are paper-mechanistic states with no canonical
  # analogue; the complex_<antigen> forms extend the registered complex family.
  paper_specific_compartments <- c(
    "n_tumor",
    "n_cart",
    "complex_cd19",
    "complex_cd22",
    "complex_cd19_cd22"
  )

  units <- list(
    time = "day",
    dosing = "cells (tumour and CAR T cells seeded into the well as amt)",
    concentration = "percent (tumour-cell viability relative to control)"
  )

  compartmentData <- list(
    n_tumor = list(analyte = "B-cell ALL tumour cells", units = "cells", specimen = "tissue", verified = TRUE),
    n_cart = list(analyte = "CD19/CD22 CAR T cells", units = "cells", specimen = "tissue", verified = TRUE),
    complex_cd19 = list(
      analyte = "CAR-CD19 single-target complexes",
      units = "complexes/L",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_cd22 = list(
      analyte = "CAR-CD22 single-target complexes",
      units = "complexes/L",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_cd19_cd22 = list(
      analyte = "CAR-CD19-CD22 double-target complexes",
      units = "complexes/L",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (B-cell ALL cell lines: RajiWT, RajiCD19KO, SupT1WT, CD19+/CD22+ SupT1)",
    n_subjects = NA_integer_,
    disease_state = "In vitro co-culture cytotoxicity assay (B-cell acute lymphoblastic leukemia)",
    notes = paste(
      "Cell viability at 24-72 h after co-culturing leukemia cell lines with",
      "non-transduced T cells, CD19 CAR-T, CD22 CAR-T or CD19/CD22 CAR-T (AUTO3)",
      "at a 1:1 effector:target ratio. KC50 and KMaxKill were estimated jointly",
      "across all cell lines and CAR-T treatments; non-targeting CARs were",
      "assigned a near-zero association rate (1e-10) with koff fixed to 1."
    )
  )

  ini({
    # ---- Estimated parameters (Table 1, cell-level B-ALL model) ----
    lkmax_kill <- log(0.8)
    label("Maximum first-order tumour-cell killing rate constant, KMaxKill (1/day)") # Table 1 B-ALL: KMaxKill = 0.8 (14.5% RSE)
    lkc50_kill <- log(3.57)
    label("CAR-target complexes per tumour cell giving half-maximal killing, KC50 (complexes/cell)") # Table 1 B-ALL: KC50CAR-T = 3.57 (53% RSE)

    # ---- Fixed structural / system parameters (Table 1, cell-level B-ALL) ----
    lgamma_kill <- fixed(log(1))
    label("Sigmoidicity factor for tumour killing, gamma (unitless)") # Table 1 B-ALL: gammaKillTumor = 1 (fixed)
    ldt_cart <- fixed(log(0.71))
    label("CAR-T-cell doubling time, DTCAR-T (day)") # dataset doubling_time_day for Raji lines = 0.71 (Table 1 B-ALL DTCAR-T = 2.3 for SupT1)
    dens_cd19_car <- fixed(15000)
    label("Anti-CD19 CAR density on CAR T cells (CARs/cell)") # Table 1 B-ALL: DensityCD19CAR = 15,000 (fixed)
    dens_cd22_car <- fixed(15000)
    label("Anti-CD22 CAR density on CAR T cells (CARs/cell)") # Table 1 B-ALL: DensityCD22CAR = 15,000 (fixed)
    dens_cd19_tumor <- fixed(50347)
    label("CD19 tumour-antigen density (antigens/cell); cell-line regressor, default RajiWT") # Table 1 B-ALL: DensityCD19Tumor RajiWT = 50,347
    dens_cd22_tumor <- fixed(55116)
    label("CD22 tumour-antigen density (antigens/cell); cell-line regressor, default RajiWT") # Table 1 B-ALL: DensityCD22Tumor RajiWT = 55,116
    kon_orig_cd19 <- fixed(1e5)
    label("Anti-CD19 CAR association rate constant, konCD19CAR (1/M/s); default AUTO3") # Table 1 B-ALL: konCD19CAR AUTO3 = 1e5
    kon_orig_cd22 <- fixed(1e5)
    label("Anti-CD22 CAR association rate constant, konCD22CAR (1/M/s); default AUTO3") # Table 1 B-ALL: konCD22CAR AUTO3 = 1e5
    kd_cd19 <- fixed(1.09)
    label("Anti-CD19 CAR equilibrium dissociation constant, KdCD19CAR (nM); default AUTO3") # Table 1 B-ALL: KdCD19CAR AUTO3 = 1.09
    kd_cd22 <- fixed(32)
    label("Anti-CD22 CAR equilibrium dissociation constant, KdCD22CAR (nM); default AUTO3") # Table 1 B-ALL: KdCD22CAR AUTO3 = 32
    tumor0 <- fixed(1e5)
    label("Initial tumour cells seeded in the well (cells)") # Table 1 B-ALL: Tumor0 = 100,000
    mv <- fixed(1e-4)
    label("Reaction-well volume (L)") # Table 1 B-ALL: Vwell = 100 uL = 1e-4 L

    # ---- Between-subject (between-replicate) variability ----
    # Monolix log-normal random-effect SDs; variances are omega^2.
    etalkmax_kill ~ 1.3456 # Table 1 B-ALL: omega KMaxKill = 1.16 (10.5% RSE); 1.16^2
    etalkc50_kill ~ 8.0656 # Table 1 B-ALL: omega KC50CAR-T = 2.84 (11.8% RSE); 2.84^2

    # ---- Residual error ----
    # No residual-error magnitude is reported for the in-vitro model (Table 1
    # carries no sigma row); encoded as fixed(0). See the vignette Errata.
    addSd_viability <- fixed(0)
    label("Additive residual error on tumour-cell viability (percent)") # not reported in the paper or its supplement
  })

  model({
    # ---- Individual parameters ----
    kmax_kill <- exp(lkmax_kill + etalkmax_kill)
    kc50_kill <- exp(lkc50_kill + etalkc50_kill)
    gamma <- exp(lgamma_kill)
    dt_cart <- exp(ldt_cart)

    # ---- Unit conversion of the binding constants (Supplement eq. notes) ----
    kon_cd19 <- kon_orig_cd19 * (60 * 60 * 24) / 6.023e23
    koff_cd19 <- kon_orig_cd19 * (kd_cd19 * 1e-9) * (60 * 60 * 24)
    kon_cd22 <- kon_orig_cd22 * (60 * 60 * 24) / 6.023e23
    koff_cd22 <- kon_orig_cd22 * (kd_cd22 * 1e-9) * (60 * 60 * 24)

    # ---- Overall and free antigen / CAR (Supplement Eqs. 1-6) ----
    ag_car19 <- dens_cd19_car * n_cart / mv
    ag_car22 <- dens_cd22_car * n_cart / mv
    ag_tum19 <- dens_cd19_tumor * n_tumor / mv
    ag_tum22 <- dens_cd22_tumor * n_tumor / mv
    fcar19 <- ag_car19 - complex_cd19 - complex_cd19_cd22
    fcar22 <- ag_car22 - complex_cd22 - complex_cd19_cd22
    ftum19 <- ag_tum19 - complex_cd19 - complex_cd19_cd22
    ftum22 <- ag_tum22 - complex_cd22 - complex_cd19_cd22

    # ---- Complex formation (Supplement Eqs. 7-9) ----
    d/dt(complex_cd19) <- kon_cd19 * ftum19 * fcar19 - koff_cd19 * complex_cd19 -
      kon_cd22 * ftum22 * complex_cd19 + koff_cd22 * complex_cd19_cd22
    d/dt(complex_cd22) <- kon_cd22 * ftum22 * fcar22 - koff_cd22 * complex_cd22 -
      kon_cd19 * ftum19 * complex_cd22 + koff_cd19 * complex_cd19_cd22
    d/dt(complex_cd19_cd22) <- kon_cd19 * fcar19 * complex_cd22 - koff_cd19 * complex_cd19_cd22 +
      kon_cd22 * fcar22 * complex_cd19 - koff_cd22 * complex_cd19_cd22

    # ---- Tumour killing and CAR-T growth (Supplement Eqs. 10-13) ----
    cplxcell <- (complex_cd19 + complex_cd22 + complex_cd19_cd22) * mv / n_tumor
    kkill <- kmax_kill * cplxcell^gamma / (kc50_kill^gamma + cplxcell^gamma)
    d/dt(n_tumor) <- -kkill * n_tumor # tumour growth rate Kg = 0 (viability is relative to control)
    d/dt(n_cart) <- log(2) / dt_cart * n_cart

    # ---- Initial conditions ----
    complex_cd19(0) <- 0
    complex_cd22(0) <- 0
    complex_cd19_cd22(0) <- 0

    # ---- Output: percent tumour-cell viability relative to control ----
    viability <- n_tumor / tumor0 * 100
    viability ~ add(addSd_viability)
  })
}
