# Cell-level in vitro pharmacodynamic model of anti-CD19 CAR T-cell killing
# (step 1 of the Salem 2023 multiscale CK-PD framework), CAT19-binder arm.
# The paper fits one model structure to the 4-h chromium-release cytotoxicity
# data of two anti-CD19 CARs (the low-affinity CAT19 and FMC63) against SupT1
# CD19+ cells; every estimated parameter (KKillMax, KC50CAR-T) and the binding
# constants are construct-specific and no random effects were estimated, so the
# joint fit separates exactly into one file per construct. The FMC63 arm is
# Salem_2023_CART_FMC63_invitro_qsp.R.
#
# SOURCES. Equations: Appendix S1 Eqs. 1-5 and the Monolix code for Figure 1 in
# Appendix S2. Parameter values: main-text Table 2, "Parameters associated with
# preclinical (in vitro) cell-level PD model".
#
# STATE UNITS. The Monolix code doses tumour and CAR T cells as concentrations
# in the 175 uL well (cells/L). They are carried here as cell NUMBERS, so a
# dose record is the number of cells seeded; the complex is the code's #/L.

Salem_2023_CART_CAT19_invitro_qsp <- function() {
  description <- "QSP. In vitro (SupT1 CD19+ cell line). Cell-level pharmacodynamic model of a low-affinity CAT19-binder anti-CD19 CAR T cell (Salem 2023, step 1). CAR T cells and tumour cells grow exponentially in a 175 uL coculture well; CARs bind CD19 (second-order kon, first-order koff) to form CAR-target complexes, and the number of complexes per tumour cell drives Emax tumour-cell killing. Output: percent cytolysis of the seeded tumour cells. Fit to 4-h chromium-release data at effector:target ratios 0.4:1 to 6.4:1."
  reference <- paste(
    "Salem AM, Mugundu GM, Singh AP (2023).",
    "Development of a multiscale mechanistic modeling framework integrating",
    "differential cellular kinetics of CAR T-cell subsets and immunophenotypes",
    "in cancer patients.",
    "CPT Pharmacometrics Syst Pharmacol 12(9):1285-1296.",
    "doi:10.1002/psp4.13009.",
    "Equations in Appendix S1 and Monolix code in Appendix S2.",
    "In vitro data digitised by the authors from Ghorashian et al. 2019",
    "(Nat Med 25:1408-1414; doi:10.1038/s41591-019-0549-5).",
    sep = " "
  )
  vignette <- "Salem_2023_CART_subsets_immunophenotypes"

  # The CAR T-cell pool of a coculture well has no canonical analogue.
  # tumor and complex are canonical.
  paper_specific_compartments <- c("cart")

  units <- list(
    time = "h",
    dosing = "cells (tumour cells as amt on tumor; CAR T cells as amt on cart, both at time 0)",
    concentration = "percent cytolysis of the seeded tumour cells (lysis)"
  )

  compartmentData <- list(
    cart = list(
      analyte = "low-affinity CAT19-binder anti-CD19 CAR T cells",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    tumor = list(
      analyte = "SupT1 CD19+ tumour cells",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    complex = list(
      analyte = "CAR-CD19 complexes",
      units = "complexes/L",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (SupT1 CD19+ cell line)",
    n_subjects = NA_integer_,
    disease_state = "SupT1 cells engineered to express CD19, cocultured with low-affinity CAT19-binder anti-CD19 CAR T cells",
    dose_range = "5000 tumour cells per well; effector:target ratios 0.4:1, 0.8:1, 1.6:1, 3.2:1 and 6.4:1",
    notes = paste(
      "Mean (+/- SEM) percent lysis at 4 h, digitised from Ghorashian et al. 2019",
      "(Table 2 references 20 and 21 of Salem 2023). No between-sample",
      "variability was estimated."
    )
  )

  ini({
    lkkill_max <- log(0.22)
    label("Maximum first-order tumour-cell killing rate constant, KKillMax (1/h)") # Table 2 in vitro: KKillMax CAT19 = 0.22 (4.5% RSE)
    lkc50_kill <- log(0.28)
    label("CAR-target complexes per tumour cell giving half-maximal killing, KC50CAR-T (complexes/cell)") # Table 2 in vitro: KC50CAR-T CAT19 = 0.28 (10.9% RSE)

    dt_tumor <- fixed(30)
    label("Doubling time of SupT1 CD19+ tumour cells, DTTumor (h)") # Table 2 in vitro: DTTumor = 30 (fixed, ref 37); code DT = 30
    dt_cart <- fixed(24)
    label("Doubling time of CAR T cells, DTCAR-T (h)") # Table 2 in vitro: DTCAR-T = 24 (fixed, ref 38); code DTCAR = 24
    ag_tumor <- fixed(4000)
    label("CD19 density on tumour cells, AgTumor (antigens/cell)") # Table 2 in vitro: AgTumor = 4000 (fixed, ref 39)
    ag_car <- fixed(15000)
    label("CAR density on CAR T cells, AgCAR (CARs/cell)") # Table 2 in vitro: AgCAR = 15,000 (fixed, ref 19)
    kon <- fixed(2.15e5)
    label("CAR-CD19 association rate constant, Kon (1/M/s)") # Table 2 in vitro: Kon CAT19 = 2.15E5 (fixed, ref 20)
    koff <- fixed(0.003)
    label("CAR-CD19 dissociation rate constant, Koff (1/s)") # Table 2 in vitro: Koff CAT19 = 0.003 (fixed, ref 20); Appendix S2 Figure 2 code comment gives 0.003096
    tumor0 <- fixed(5000)
    label("Tumour cells seeded per well, Tumor0 (cells)") # Table 2 in vitro: Tumor0 = 5000 (fixed, ref 21)
    v_well <- fixed(175e-6)
    label("Volume of the coculture well, Vwell (L)") # Table 2 in vitro: Vwell = 175 uL (fixed, ref 21); code p = 1/175E-6
  })

  model({
    kkill_max <- exp(lkkill_max)
    kc50_kill <- exp(lkc50_kill)

    # Appendix S2 code: Kon (1/M/s) -> 1/(#/L)/h; Koff (1/s) -> 1/h.
    kon_h <- kon * (60 * 60) / 6.023e23
    koff_h <- koff * (60 * 60)

    # Appendix S1 Eqs. 1 and 4-5 with cells converted to concentrations
    car_free <- cart / v_well * ag_car - complex
    ag_free <- tumor / v_well * ag_tumor - complex
    cplx_per_tumor <- complex / (tumor / v_well)
    kkill <- kkill_max * cplx_per_tumor / (kc50_kill + cplx_per_tumor)

    d/dt(cart) <- log(2) / dt_cart * cart # Appendix S1 Eq. 3
    d/dt(tumor) <- log(2) / dt_tumor * tumor - kkill * tumor # Appendix S1 Eq. 2
    d/dt(complex) <- kon_h * car_free * ag_free - koff_h * complex # Appendix S1 Eq. 1

    # Appendix S2 code: viability relative to the seeded count
    lysis <- 100 - tumor / tumor0 * 100
  })
}
