# Multiscale mechanistic cellular-kinetic / pharmacodynamic (PK-PD) model for
# the anti-BCMA CAR T-cell product bb2121 (idecabtagene vicleucel) in
# relapsed/refractory multiple myeloma, published by Singh et al. 2021
# (CPT Pharmacometrics Syst Pharmacol 10:362-376). This file carries the
# CLINICAL model (Figure 1b fitted to the phase 1 CRB-401 data), which is the
# paper's centrepiece and the only one of its three sub-models whose exact
# equations are on disk. The paper also reports a preclinical RPMI-8226
# xenograft fit (Figure 2b) and a cell-level in vitro cytotoxicity fit
# (Figure 2a); the authors deposited the Monolix code for the CLINICAL model
# only, and a faithful reconstruction of the preclinical model from Table 1
# does not reproduce Figure 2b (see the vignette Errata and the extraction
# report), so those two sub-models are documented but not shipped as validated
# model files.
#
# SOURCES. Equations: Supporting Information (PSP4-10-362-s001.pdf) Eqs. 6-20.
# Parameter values: main-text Table 1, "Parameters associated with the clinical
# PK-PD model". The authors' Monolix model code (PSP4-10-362-s002.docx, "Model
# Code in Monolix v2") is the as-run form and settles three details the printed
# equations leave open (see the trailing comments below and the vignette
# Errata): the tumour-kill term multiplies Tumor (not CARTe_T as SI Eq. 18
# prints), the initial tumour burden (2.5e9 cells/L of bone marrow), and the
# tiny non-zero initial tissue effector pool that keeps the complexes-per-CAR-T
# ratio finite before any cell has distributed.
#
# STATE UNITS. The Monolix code carries the four CAR-T pools as concentrations
# (cells/L of blood or of bone marrow) and doses "cells/L". They are carried
# here as cell NUMBERS (concentration x compartment volume), so a dose record
# is simply the number of CAR+ T cells infused; the transformation is exact
# because every CAR-T flux in the code is a rate constant times a
# concentration times a fixed volume. Complex and tumour stay as the code's
# bone-marrow concentrations (#/L, cells/L).

Singh_2021_idecabtageneVicleucel_human <- function() {
  description <- "QSP (multiscale mechanistic cellular-kinetic / pharmacodynamic model). Anti-BCMA CAR T-cell therapy bb2121 (idecabtagene vicleucel) in adults with relapsed/refractory multiple myeloma, clinical model. Effector and memory CAR T cells distribute between blood and bone marrow (first-order K12 / K21) and are eliminated from blood at phenotype-specific rates; in bone marrow the CARs bind BCMA on myeloma cells (second-order kon / first-order koff) to form CAR-target complexes, whose number per CAR T cell drives Emax expansion of the effector pool and whose number per tumour cell drives Emax killing of an exponentially growing tumour. Effector cells convert to memory at a net first-order rate. Serum M-protein and soluble BCMA are turnover biomarkers whose production scales with (tumour / baseline tumour)^gamma. Outputs: blood transgene copies per ug genomic DNA and percent change from baseline of soluble BCMA and serum M-protein. Mean parameters were fit to the mean phase 1 data; between-subject variability on expansion, killing and the biomarker exponents was fit to the individual M-protein response categories, and 50% variability on the four disposition rate constants is the authors' assumption."
  reference <- paste(
    "Singh AP, Chen W, Zheng X, Mody H, Carpenter TJ, Zong A, Heald DL (2021).",
    "Bench-to-bedside translation of chimeric antigen receptor (CAR) T cells",
    "using a multiscale systems pharmacokinetic-pharmacodynamic model:",
    "A case study with anti-BCMA CAR-T.",
    "CPT Pharmacometrics Syst Pharmacol 10(4):362-376.",
    "doi:10.1002/psp4.12598.",
    "Equations in the Supporting Information (PSP4-10-362-s001.pdf) and",
    "Monolix model code in PSP4-10-362-s002.docx (PMC8099446 open-access package).",
    "Fit to the phase 1 CRB-401 data of Raje et al. 2019 (N Engl J Med 380:1726-1737;",
    "doi:10.1056/NEJMoa1817226).",
    sep = " "
  )
  vignette <- "Singh_2021_idecabtageneVicleucel"

  # Effector / memory CAR T-cell pools in blood and tissue are
  # paper-mechanistic states with no canonical analogue; the same four names
  # are used by Hardiansyah_2019_CART_upn1_qsp.R for the same four pools.
  # sbcma (serum soluble BCMA as a tumour-burden turnover biomarker) is not
  # registered either. complex, tumor and mprotein are canonical.
  paper_specific_compartments <- c(
    "carte_pb",
    "carte_t",
    "cartm_pb",
    "cartm_t",
    "sbcma"
  )

  units <- list(
    time = "day",
    dosing = "cells (CAR+ T cells; input as amt on the carte_pb compartment)",
    concentration = "transgene copies/ug genomic DNA (transgene); percent change from baseline (sbcma_pctchg, mprotein_pctchg)"
  )

  compartmentData <- list(
    carte_pb = list(
      analyte = "effector anti-BCMA CAR T cells (idecabtagene vicleucel)",
      units = "cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    carte_t = list(
      analyte = "effector anti-BCMA CAR T cells (idecabtagene vicleucel)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    cartm_pb = list(
      analyte = "memory anti-BCMA CAR T cells (idecabtagene vicleucel)",
      units = "cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    cartm_t = list(
      analyte = "memory anti-BCMA CAR T cells (idecabtagene vicleucel)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    complex = list(
      analyte = "CAR-BCMA complexes on myeloma cells",
      units = "complexes/L",
      specimen = "tissue",
      verified = TRUE
    ),
    tumor = list(
      analyte = "BCMA-expressing multiple myeloma cells",
      units = "cells/L",
      specimen = "tissue",
      verified = TRUE
    ),
    mprotein = list(
      analyte = "serum M-protein",
      units = "g/L",
      specimen = "serum",
      verified = TRUE
    ),
    sbcma = list(
      analyte = "soluble B-cell maturation antigen",
      units = "ng/mL",
      specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 33L,
    n_studies = 1L,
    disease_state = "Relapsed or refractory multiple myeloma (phase 1 CRB-401 dose-escalation study of bb2121)",
    dose_range = "Single intravenous infusion of 50, 150, 450 or 800 x 10^6 CAR+ T cells (flat dose)",
    regions = "United States (CRB-401; Raje et al. 2019)",
    notes = paste(
      "The mean parameters were estimated from mean (+/- SD) blood transgene,",
      "soluble BCMA and serum M-protein profiles per dose cohort digitized from",
      "Raje et al. 2019; no individual cellular-kinetic data were available.",
      "Between-subject variability was then estimated from the individual",
      "IMWG response categories of all 33 patients (transformed into",
      "continuous M-protein changes), with the mean parameters fixed.",
      "Demographics are not reproduced in Singh 2021."
    )
  )

  ini({
    # ---- Estimated parameters (Table 1, clinical PK-PD model) ----
    lkexp_max <- log(1.73)
    label("Maximum first-order CAR T-cell expansion rate constant, Kexp_max (1/day)") # Table 1 clinical: K_Exp_max = 1.73 (10% RSE)
    lec50_exp <- log(10)
    label("CAR-target complexes per CAR T cell giving half-maximal expansion, EC50_Exp (complexes/cell)") # Table 1 clinical: EC_Exp_50 = 10 (18% RSE)
    lrm <- log(0.00002)
    label("Net first-order effector-to-memory conversion rate constant, Rm (1/day)") # Table 1 clinical: Rm = 0.00002 (66% RSE)
    lkel_e <- log(113)
    label("First-order elimination rate constant of effector CAR T cells from blood, Kel_e (1/day)") # Table 1 clinical: Kel_e = 113 (19% RSE)
    lkel_m <- log(0.219)
    label("First-order elimination rate constant of memory CAR T cells from blood, Kel_m (1/day)") # Table 1 clinical: Kel_m = 0.219 (13% RSE)
    lk12 <- log(1.71)
    label("Blood-to-bone-marrow distribution rate constant, K12 (1/day)") # Table 1 clinical: K12 = 1.71 (11% RSE)
    lk21 <- log(0.176)
    label("Bone-marrow-to-blood redistribution rate constant, K21 (1/day)") # Table 1 clinical: K21 = 0.176 (14% RSE)
    lkkill_max <- log(0.343)
    label("Maximum first-order tumour-cell killing rate constant, Kkill_max (1/day)") # Table 1 clinical: K_Kill_max = 0.343 (21% RSE)
    lgam_mprotein <- log(0.215)
    label("Exponent of relative tumour burden on M-protein production, gamma_m (unitless)") # Table 1 clinical: gamma_m = 0.215 (5% RSE)
    lgam_sbcma <- fixed(log(1))
    label("Exponent of relative tumour burden on soluble BCMA production, gamma_b (unitless)") # Table 1 clinical: gamma_b = 1 (Fixed)

    # ---- Fixed system- and drug-specific parameters (Table 1, clinical) ----
    lkc50_kill <- fixed(log(2.24))
    label("CAR-target complexes per tumour cell giving half-maximal killing, KC50_CAR-T (complexes/cell)") # Table 1 clinical: KC_CAR-T_50 = 2.24 (Fixed, in vitro estimate)
    lkg <- fixed(log(0.008))
    label("First-order tumour growth rate constant, Kg_Tumor (1/day)") # Table 1 clinical: Kg_Tumor = 0.008 (Fixed); code Kg_tumor0 = 0.008
    kon <- fixed(7.103e4)
    label("CAR-BCMA association rate constant, Kon (1/M/s)") # Table 1 clinical: Kon = 7.1E4 (Fixed); code Kon_orig = 7.103E+4 (as-run value)
    koff <- fixed(2.385e-3)
    label("CAR-BCMA dissociation rate constant, Koff (1/s)") # Table 1 clinical: Koff = 2.39E-3 (Fixed); code Koff_orig = 2.385E-3 (as-run value)
    ag_car <- fixed(15000)
    label("CAR density on CAR T cells, Ag_CAR (CARs/cell)") # Table 1 clinical: Ag_CAR = 15,000 (Fixed); code Density_CAR = 15000
    ag_tumor <- fixed(12590)
    label("BCMA density on myeloma cells, Ag_Tumor (antigens/cell)") # Table 1 clinical: Ag_Tumor = 12,590 (Fixed); code Density_TAA = 12590
    transc <- fixed(0.002)
    label("Conversion factor from blood CAR T cells (cells/L) to transgene copies/ug genomic DNA, TransC") # Table 1 clinical: TransC = 0.002 (Fixed); code TransC = 0.002
    v_blood <- fixed(5)
    label("Blood compartment volume, Vb (L)") # Table 1 clinical: Vb = 5 (Fixed); code Vb = 5
    v_bonemarrow <- fixed(3.65)
    label("Bone marrow (tissue) compartment volume, Vbm (L)") # Table 1 clinical: Vbm = 3.65 (Fixed); code Vt = 3.65
    tumor0 <- fixed(2.5e9)
    label("Baseline myeloma-cell concentration in bone marrow, Tumor0 (cells/L)") # Monolix code (s002): Tumor_T_0 = 2.5E9 'tumor in bone marrow (cells/L)'
    kdeg_mprotein <- fixed(0.117)
    label("Serum M-protein degradation rate constant, Km (1/day)") # Table 1 clinical: Km = 0.117 (Fixed); code 0.117
    ksyn_mprotein <- fixed(12.1)
    label("Serum M-protein production rate per tumour cell, Pm (pg/cell/day)") # Table 1 clinical: Pm = 12.1 (Fixed); code p = 12.1*2.5E9
    kdeg_sbcma <- fixed(0.7)
    label("Soluble BCMA degradation rate constant, Kb (1/day)") # Table 1 clinical: Kb = 0.7 (Fixed); code k = 0.7
    ksyn_sbcma <- fixed(0.175)
    label("Soluble BCMA production rate per tumour cell, Pb (pg/cell/day)") # Table 1 clinical: Pb = 0.175 (Fixed, calculated)

    # ---- Between-subject variability ----
    # Table 1's omega column is the Monolix standard deviation of the
    # log-normal random effect; variances below are omega^2. The paper itself
    # calls omega = 0.62 "~60% IIV", so the SI's "50% IIV" is omega = 0.5.
    etalkexp_max ~ 0.0484 # Table 1 clinical: omega K_Exp_max = 0.22 (38% RSE); 0.22^2
    etalrm ~ fixed(1.6e-9) # Table 1 clinical: omega Rm = 0.00004, not estimated; variance 0.00004^2
    etalkkill_max ~ 0.25 # Table 1 clinical: omega K_Kill_max = 0.50 (37% RSE); 0.50^2
    etalgam_mprotein ~ 0.0441 # Table 1 clinical: omega gamma_m = 0.21 (44% RSE); 0.21^2
    etalgam_sbcma ~ 1 # Table 1 clinical: omega gamma_b = 1.0 (31% RSE); 1.0^2
    etalk12 ~ fixed(0.25) # SI p.4 lines 169-171 and 181-182: assumed 50% IIV on K12; 0.5^2
    etalk21 ~ fixed(0.25) # SI p.4 lines 169-171 and 181-182: assumed 50% IIV on K21; 0.5^2
    etalkel_e ~ fixed(0.25) # SI p.4 lines 169-171 and 181-182: assumed 50% IIV on Kel_e; 0.5^2
    etalkel_m ~ fixed(0.25) # SI p.4 lines 169-171 and 181-182: assumed 50% IIV on Kel_m; 0.5^2
  })

  model({
    # ---- Individual parameters ----
    kexp_max <- exp(lkexp_max + etalkexp_max)
    ec50_exp <- exp(lec50_exp)
    rm <- exp(lrm + etalrm)
    kel_e <- exp(lkel_e + etalkel_e)
    kel_m <- exp(lkel_m + etalkel_m)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    kkill_max <- exp(lkkill_max + etalkkill_max)
    gam_mprotein <- exp(lgam_mprotein + etalgam_mprotein)
    gam_sbcma <- exp(lgam_sbcma + etalgam_sbcma)
    kc50_kill <- exp(lkc50_kill)
    kg <- exp(lkg)

    # ---- Unit conversion of the binding constants (Monolix code) ----
    # Kon (1/M/s) -> 1/(#/L)/day; Koff (1/s) -> 1/day.
    kon_d <- kon * (60 * 60 * 24) / 6.023e23
    koff_d <- koff * (60 * 60 * 24)

    # ---- Bone-marrow concentrations (SI Eqs. 11-15) ----
    cart_t_conc <- (carte_t + cartm_t) / v_bonemarrow # CAR T cells / L of marrow
    car_free <- cart_t_conc * ag_car - complex # SI Eq. 11, free CARs (#/L)
    ag_free <- tumor * ag_tumor - complex # SI Eq. 12, free BCMA (#/L)
    cplx_per_cart <- complex / cart_t_conc # SI Eq. 14
    cplx_per_tumor <- complex / tumor # SI Eq. 15

    # SI Eq. 16: antigen-driven effector expansion
    kexp <- kexp_max * cplx_per_cart / (ec50_exp + cplx_per_cart)
    # SI Eq. 17: CAR-T mediated tumour killing
    kkill <- kkill_max * cplx_per_tumor / (kc50_kill + cplx_per_tumor)

    # ---- CAR T-cell kinetics (SI Eqs. 6, 7, 9, 10 in cell numbers) ----
    d/dt(carte_pb) <- -k12 * carte_pb + k21 * carte_t - kel_e * carte_pb
    d/dt(cartm_pb) <- -k12 * cartm_pb + k21 * cartm_t - kel_m * cartm_pb
    d/dt(carte_t) <- k12 * carte_pb - k21 * carte_t + kexp * carte_t - rm * carte_t
    d/dt(cartm_t) <- k12 * cartm_pb - k21 * cartm_t + rm * carte_t

    # ---- CAR-BCMA binding in bone marrow (SI Eq. 13) ----
    d/dt(complex) <- kon_d * car_free * ag_free - koff_d * complex

    # ---- Tumour (Monolix code ddt_Tumor_T; see file header) ----
    d/dt(tumor) <- kg * tumor - kkill * tumor

    # ---- Biomarker turnover (SI Eqs. 19-20) ----
    # Production in pg/L/day (pg/cell/day x cells/L); 1e-12 converts M-protein
    # to g/L and 1e-6 converts sBCMA to ng/mL (1 ng/mL = 1e6 pg/L).
    d/dt(mprotein) <- ksyn_mprotein * 1e-12 * tumor0 * (tumor / tumor0)^gam_mprotein -
      kdeg_mprotein * mprotein
    d/dt(sbcma) <- ksyn_sbcma * 1e-6 * tumor0 * (tumor / tumor0)^gam_sbcma -
      kdeg_sbcma * sbcma

    # ---- Initial conditions ----
    # Monolix code: CARTe_T_0 = 1E-6 cells/L (keeps cplx_per_cart finite
    # before any CAR T cell has reached the marrow); Tumor_T_0 = tumor0.
    carte_t(0) <- 1e-6 * v_bonemarrow
    tumor(0) <- tumor0
    mprotein0 <- ksyn_mprotein * 1e-12 * tumor0 / kdeg_mprotein
    sbcma0 <- ksyn_sbcma * 1e-6 * tumor0 / kdeg_sbcma
    mprotein(0) <- mprotein0
    sbcma(0) <- sbcma0

    # ---- Outputs (SI Eq. 8; Monolix code BChange / MChange) ----
    transgene <- transc * (carte_pb + cartm_pb) / v_blood
    sbcma_pctchg <- (sbcma - sbcma0) * 100 / sbcma0
    mprotein_pctchg <- (mprotein - mprotein0) * 100 / mprotein0
  })
}
