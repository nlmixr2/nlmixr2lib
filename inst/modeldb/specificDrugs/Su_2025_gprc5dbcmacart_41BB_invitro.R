# In-vitro cell-level pharmacodynamic (tumour-killing) model for the dual-targeted
# bicistronic GPRC5D/BCMA CAR T-cell construct with matched 4-1BB co-stimulatory
# domains (GPRC5D-41BB-BCMA-41BB, "G4B4"), published by Su et al. 2025 (CPT
# Pharmacometrics Syst Pharmacol 14:229-245, Figure 3). This is the independent
# cytotoxicity dataset (Fernandez de Larrea et al. 2020) used to validate the
# dual-targeting killing-model structure of Su_2025_cd19cd22cart_invitro.R and to
# compare the potency of two co-stimulatory-domain designs; the CD28/41BB
# construct is the companion file Su_2025_gprc5dbcmacart_CD28_invitro.R.
#
# SOURCES. Equations: Supplement 1 Eqs. 1-14 and the as-run Monolix model code
# in Supplement 2 ("Cell-Level Pharmacodynamic Model (MM), GPRC5D-41BB-BCMA-41BB
# CAR-T"). Estimated values and fixed system constants: main-text Table 1,
# "Parameters associated with cell-level pharmacodynamic model (MM)", G4B4 rows.
# No residual-error magnitude is reported; encoded as fixed(0). Antigen densities
# are per-cell-line regressors; defaults here are the RPMI8226 line, and the
# vignette overrides them to reproduce Figure 3a-d.
#
# STATE UNITS. n_tumor and n_cart are cell NUMBERS in a well of MV = 1e-4 L; the
# three complexes are numbers/L. Output is percent tumour-cell viability (the
# paper plots its complement, cytotoxicity = 100 - viability).

Su_2025_gprc5dbcmacart_41BB_invitro <- function() {
  description <- "In vitro (multiple myeloma cell lines). Cell-level tumour-killing model for the dual-targeted bicistronic GPRC5D/BCMA CAR T-cell construct with matched 4-1BB co-stimulatory domains (GPRC5D-41BB-BCMA-41BB). The two CARs bind GPRC5D and BCMA (second-order kon derived from koff/KD) to form single-target and double-target complexes; complexes per tumour cell drive Emax killing while CAR T cells grow exponentially. Output: percent tumour-cell viability relative to control (paper reports cytotoxicity = 100 - viability). Estimated potency KC50 = 11.8 complexes/tumour cell, compared against the CD28/41BB construct to show the effect of co-stimulatory-domain design. Antigen densities are cell-line regressors; defaults are the RPMI8226 line."
  reference <- paste(
    "Su MC, Dey A, Maddah E, Mugundu GM, Singh AP (2025).",
    "Quantitative pharmacology of dual-targeted bicistronic CAR-T-cell therapy",
    "using multiscale mechanistic modeling.",
    "CPT Pharmacometrics Syst Pharmacol 14(2):229-245.",
    "doi:10.1002/psp4.13259.",
    "Equations and Monolix model code in the Supporting Information",
    "(PMC11812944 open-access package).",
    "In-vitro cytotoxicity data digitized from Fernandez de Larrea et al. 2020",
    "(Blood Cancer Discov 1:146-154).",
    sep = " "
  )
  vignette <- "Su_2025_dual_targeted_cart"

  paper_specific_compartments <- c(
    "n_tumor",
    "n_cart",
    "complex_bcma",
    "complex_gprc5d",
    "complex_bcma_gprc5d"
  )

  units <- list(
    time = "day",
    dosing = "cells (tumour and CAR T cells seeded into the well as amt)",
    concentration = "percent (tumour-cell viability relative to control)"
  )

  compartmentData <- list(
    n_tumor = list(analyte = "multiple myeloma tumour cells", units = "cells", specimen = "tissue", verified = TRUE),
    n_cart = list(
      analyte = "GPRC5D/BCMA CAR T cells (41BB/41BB)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_bcma = list(
      analyte = "CAR-BCMA single-target complexes",
      units = "complexes/L",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_gprc5d = list(
      analyte = "CAR-GPRC5D single-target complexes",
      units = "complexes/L",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_bcma_gprc5d = list(
      analyte = "CAR-BCMA-GPRC5D double-target complexes",
      units = "complexes/L",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (multiple myeloma cell lines: RPMI8226, OPM2WT, OPM2BCMA-KO, MM1S)",
    n_subjects = NA_integer_,
    disease_state = "In vitro co-culture cytotoxicity assay (multiple myeloma)",
    notes = paste(
      "Cytotoxicity at 24 h after co-culturing MM cell lines with the",
      "GPRC5D-41BB-BCMA-41BB dual-targeted CAR-T at effector:target ratios",
      "from 0:1 to 1:1. KC50, KMaxKill and the sigmoidicity factor were",
      "estimated; CAR affinities and antigen densities were fixed per cell line."
    )
  )

  ini({
    # ---- Estimated parameters (Table 1, cell-level MM model, G4B4 rows) ----
    lkmax_kill <- log(4.39)
    label("Maximum first-order tumour-cell killing rate constant, KMaxKill (1/day)") # Table 1 MM G4B4: KMaxKill = 4.39 (4.78% RSE)
    lkc50_kill <- log(11.8)
    label("CAR-target complexes per tumour cell giving half-maximal killing, KC50 (complexes/cell)") # Table 1 MM G4B4: KC50CAR-T = 11.8 (30.5% RSE)
    lgamma_kill <- log(0.69)
    label("Sigmoidicity factor for tumour killing, gamma (unitless)") # Table 1 MM G4B4: gammaKillTumor = 0.69 (2.04% RSE)

    # ---- Fixed structural / system parameters (Table 1, cell-level MM) ----
    ldt_cart <- fixed(log(1.4))
    label("CAR-T-cell doubling time, DTCAR-T (day)") # dataset DT = 1.4 (Table 1 MM: DTCAR-T = 1)
    dens_bcma_car <- fixed(15000)
    label("Anti-BCMA CAR density on CAR T cells (CARs/cell)") # Table 1 MM: DensityBCMACAR = 15,000 (assumed)
    dens_gprc5d_car <- fixed(15000)
    label("Anti-GPRC5D CAR density on CAR T cells (CARs/cell)") # Table 1 MM: DensityGPRC5DCAR = 15,000 (assumed)
    dens_bcma_tumor <- fixed(12590)
    label("BCMA tumour-antigen density (antigens/cell); cell-line regressor, default RPMI8226") # Table 1 MM: AgBCMATumor RPMI8226 = 12,590
    dens_gprc5d_tumor <- fixed(1184)
    label("GPRC5D tumour-antigen density (antigens/cell); cell-line regressor, default RPMI8226") # Table 1 MM: AgGPRC5DTumor RPMI8226 = 1184
    koff_orig_bcma <- fixed(2.39e-3)
    label("Anti-BCMA CAR dissociation rate constant, koffBCMACAR (1/s)") # Table 1 MM: koffBCMACAR = 2.39e-3
    kd_bcma <- fixed(8.1)
    label("Anti-BCMA CAR equilibrium dissociation constant, KdBCMACAR (nM)") # Table 1 MM: KdBCMACAR = 8.1
    koff_orig_gprc5d <- fixed(6.15e-4)
    label("Anti-GPRC5D CAR dissociation rate constant, koffGPRC5DCAR (1/s)") # Table 1 MM: koffGPRC5DCAR = 6.15e-4
    kd_gprc5d <- fixed(6.15)
    label("Anti-GPRC5D CAR equilibrium dissociation constant, KdGPRC5DCAR (nM)") # Table 1 MM: KdGPRC5DCAR = 6.15
    tumor0 <- fixed(2e4)
    label("Initial tumour cells seeded in the well (cells)") # Table 1 MM: Tumor0 = 20,000
    mv <- fixed(1e-4)
    label("Reaction-well volume (L)") # Table 1 MM: Vwell = 100 uL = 1e-4 L

    # ---- Between-replicate variability (Monolix log-normal SD; variance omega^2) ----
    etalkc50_kill ~ 1.2544 # Table 1 MM G4B4: omega KC50CAR-T = 1.12 (19.2% RSE); 1.12^2

    # ---- Residual error ----
    # No residual-error magnitude is reported (Table 1 carries no sigma row);
    # encoded as fixed(0). See the vignette Errata.
    addSd_viability <- fixed(0)
    label("Additive residual error on tumour-cell viability (percent)") # not reported in any on-disk source
  })

  model({
    # ---- Individual parameters ----
    kmax_kill <- exp(lkmax_kill)
    kc50_kill <- exp(lkc50_kill + etalkc50_kill)
    gamma <- exp(lgamma_kill)
    dt_cart <- exp(ldt_cart)

    # ---- Unit conversion; kon derived from koff and KD (Monolix MM code) ----
    kon_orig_bcma <- koff_orig_bcma / (kd_bcma * 1e-9)
    kon_bcma <- kon_orig_bcma * (60 * 60 * 24) / 6.023e23
    koff_bcma <- koff_orig_bcma * (60 * 60 * 24)
    kon_orig_gprc5d <- koff_orig_gprc5d / (kd_gprc5d * 1e-9)
    kon_gprc5d <- kon_orig_gprc5d * (60 * 60 * 24) / 6.023e23
    koff_gprc5d <- koff_orig_gprc5d * (60 * 60 * 24)

    # ---- Overall and free antigen / CAR (Supplement Eqs. 1-6) ----
    ag_car_bcma <- dens_bcma_car * n_cart / mv
    ag_car_gprc5d <- dens_gprc5d_car * n_cart / mv
    ag_tum_bcma <- dens_bcma_tumor * n_tumor / mv
    ag_tum_gprc5d <- dens_gprc5d_tumor * n_tumor / mv
    fcar_bcma <- ag_car_bcma - complex_bcma - complex_bcma_gprc5d
    fcar_gprc5d <- ag_car_gprc5d - complex_gprc5d - complex_bcma_gprc5d
    ftum_bcma <- ag_tum_bcma - complex_bcma - complex_bcma_gprc5d
    ftum_gprc5d <- ag_tum_gprc5d - complex_gprc5d - complex_bcma_gprc5d

    # ---- Complex formation (Supplement Eqs. 7-9, MM code) ----
    d/dt(complex_bcma) <- kon_bcma * fcar_bcma * ftum_bcma - koff_bcma * complex_bcma +
      koff_gprc5d * complex_bcma_gprc5d - kon_gprc5d * complex_bcma * ftum_gprc5d
    d/dt(complex_gprc5d) <- kon_gprc5d * fcar_gprc5d * ftum_gprc5d - koff_gprc5d * complex_gprc5d +
      koff_bcma * complex_bcma_gprc5d - kon_bcma * complex_gprc5d * ftum_bcma
    d/dt(complex_bcma_gprc5d) <- kon_gprc5d * complex_bcma * ftum_gprc5d - koff_gprc5d * complex_bcma_gprc5d +
      kon_bcma * complex_gprc5d * ftum_bcma - koff_bcma * complex_bcma_gprc5d

    # ---- Tumour killing and CAR-T growth (Supplement Eqs. 10-13) ----
    cplxcell <- (complex_bcma + complex_gprc5d + complex_bcma_gprc5d) * mv / n_tumor
    kkill <- kmax_kill * cplxcell^gamma / (kc50_kill^gamma + cplxcell^gamma)
    d/dt(n_tumor) <- -kkill * n_tumor # tumour growth rate Kg = 0 (viability is relative to control)
    d/dt(n_cart) <- log(2) / dt_cart * n_cart

    # ---- Initial conditions ----
    complex_bcma(0) <- 0
    complex_gprc5d(0) <- 0
    complex_bcma_gprc5d(0) <- 0

    # ---- Output: percent tumour-cell viability relative to control ----
    viability <- n_tumor / tumor0 * 100
    viability ~ add(addSd_viability)
  })
}
