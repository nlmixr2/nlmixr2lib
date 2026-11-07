Choi_2024_cd123nkce_invitro <- function() {
  description <- "In vitro (MOLM-13 AML cell line with resting human NK cells). Mechanistic synapse-driven PK/PD model of the trifunctional NKp46/CD16a x CD123 natural killer cell engager (CD123-NKCE): sequential mass-action binding of the engager to NKp46 on NK cells and CD123 on MOLM-13 cells forms drug-NKp46 and drug-CD123 binary complexes and the NKp46-drug-CD123 synapse (trimer), which drives first-order Hill-type killing of exponentially growing CD123+ target cells. NK cells are held static; CD123 receptors follow the surviving target cells; an untreated control arm grows in parallel so that percent cytotoxicity can be computed (Choi 2024 Appendix S1 Eqs 1-10)."
  reference <- paste(
    "Choi SL, Valente D, Virone-Oddos A, Mauriac C.",
    "Developing a mechanistic translational PK/PD model for a trifunctional NK cell",
    "engager to predict the first-in-human dose for acute myeloid leukemia.",
    "Clin Transl Sci. 2024;17(1):e13689. doi:10.1111/cts.13689",
    sep = " "
  )
  vignette <- "Choi_2024_cd123nkce"
  # Time is hours. The drug state holds the well concentration (nM): an
  # evid = 1 record with amt = B0 (nM) at time 0 sets the incubation
  # concentration, the same role as the initial condition B(0) = B0 of
  # Appendix S1 Eq 1. Cell states hold cells/uL.
  units <- list(time = "h", dosing = "nmol/L", concentration = "nmol/L")

  # The drug-receptor binary complexes, the target-cell pool and its untreated
  # control arm have no canonical compartment name; they are named after the
  # drug-<antigen>-<site> family of drug_cd3_central and the tcell_central
  # cell-pool family.
  paper_specific_compartments <- c(
    "drug_nkp46_central",
    "drug_cd123_central",
    "cd123cell_central",
    "cd123cell_control"
  )

  compartmentData <- list(
    central = list(analyte = "CD123-NKCE", units = "nmol/L", specimen = "not applicable", verified = TRUE),
    drug_nkp46_central = list(
      analyte = "CD123-NKCE-NKp46 binary complex",
      units = "nmol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    drug_cd123_central = list(
      analyte = "CD123-NKCE-CD123 binary complex",
      units = "nmol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    trimer = list(
      analyte = "NKp46-CD123-NKCE-CD123 synapse",
      units = "nmol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    cd123cell_central = list(
      analyte = "CD123+ MOLM-13 cells, drug-treated well",
      units = "cells/uL",
      specimen = "not applicable",
      verified = TRUE
    ),
    cd123cell_control = list(
      analyte = "CD123+ MOLM-13 cells, untreated control well",
      units = "cells/uL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (MOLM-13 AML cell line co-cultured with resting human NK cells)",
    n_subjects = NA_integer_,
    n_studies = 1,
    disease_state = "Acute myeloid leukemia cell line (CD123+ MOLM-13) as target cells",
    dose_range = "Eight CD123-NKCE concentrations, 1.5e-3 to 1.5e4 ng/mL (1.02e-5 to 102 nmol/L)",
    notes = paste(
      "51Cr-release cytotoxicity of resting NK cells from two donors against MOLM-13 at",
      "effector-to-target ratios 10:1, 5:1, 1:1 and 1:5 (EF0 = 200, 100, 20 and 4 NK cells/uL",
      "against TC0 = 20 MOLM-13 cells/uL) after 1, 4 and 20 h of incubation (Choi 2024 Methods,",
      "'In vitro cytotoxicity experiments', and Table 1). Fit in Monolix 2020 (SAEM)."
    )
  )

  ini({
    # Synapse-driven killing (Table 1, in vitro column; estimated)
    lkmax <- log(0.332)
    label("Maximal CD123+ cell killing rate constant of the synapse (1/h)") # Table 1 in vitro: kmax = 0.332 (RSE 8.42%)
    lec50 <- log(1.46e-9)
    label("Synapse concentration at half-maximal killing rate (nmol/L)") # Table 1 in vitro: EC50 = 1.46E-9 (RSE 21.7%)
    lhill <- log(1.04)
    label("Hill coefficient of synapse-driven killing (unitless)") # Table 1 in vitro: gamma = 1.04 (RSE 0.153%)

    # Target-cell growth (Table 1, in vitro column)
    kg <- fixed(0.0139)
    label("First-order growth rate constant of CD123+ MOLM-13 cells (1/h)") # Table 1 in vitro: kg = 0.0139 (footnote h: doubling time ~50 h)

    # Receptor densities and binding (Table 1, in vitro column)
    nkp46_receptors <- fixed(4600)
    label("NKp46 receptors per NK cell (receptors/cell)") # Table 1 in vitro: RTdE = 4600 (internal data)
    cd123_receptors <- fixed(4360)
    label("CD123 receptors per MOLM-13 cell (receptors/cell)") # Table 1 in vitro: RTdT = 4360 (internal data)
    kd_nkp46 <- fixed(16.6)
    label("Equilibrium dissociation constant of CD123-NKCE for NKp46 (nmol/L)") # Table 1 in vitro: KD1 = 16.6 (footnote l: internal measurement for humans)
    kd_cd123 <- fixed(0.404)
    label("Equilibrium dissociation constant of CD123-NKCE for CD123 (nmol/L)") # Table 1 in vitro: KD2 = 0.404 (footnote l: internal measurement for humans)
    koff_nkp46 <- fixed(0.001)
    label("Dissociation rate constant of CD123-NKCE from NKp46 (1/s)") # Table 1 in vitro: koff1 = 0.001 1/s (footnote o: ref 15)
    koff_cd123 <- fixed(0.001)
    label("Dissociation rate constant of CD123-NKCE from CD123 (1/s)") # Table 1 in vitro: koff2 = 0.001 1/s (footnote o: ref 15)

    # Experimental cell numbers (Table 1 footnote b: actual values used as regressors).
    # rbase_nk defaults to the 10:1 E:T arm; set it per well (200, 100, 20 or 4).
    rbase_nk <- fixed(200)
    label("NK (effector) cell concentration in the well (cells/uL)") # Table 1 in vitro: EF0 = 4-200 (regressor; 200 = E:T 10:1)
    rbase_cd123 <- fixed(20)
    label("Initial CD123+ MOLM-13 (target) cell concentration (cells/uL)") # Table 1 in vitro: TC0 = 20 (regressor)

    # Interindividual (inter-donor / inter-experiment) variability
    etalec50 ~ 0.996004 # Table 1 in vitro: IIV on EC50 = 0.998 (SD of the random effect, footnote a); 0.998^2
  })

  model({
    # Avogadro's number per nmol (Appendix S1 Table S2: NA = 6.022e14 #/nmol)
    navo <- 6.022e14

    kmax <- exp(lkmax)
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill)

    # koff in 1/s -> 1/h; kon = koff / KD (Appendix S1 'kon1 = koff1/KD1')
    koff1 <- koff_nkp46 * 3600
    koff2 <- koff_cd123 * 3600
    kon1 <- koff1 / kd_nkp46
    kon2 <- koff2 / kd_cd123

    # Free receptors (nmol/L). NK cells are static (Appendix S1 Eq 9). Total
    # CD123 follows the surviving treated target cells, as in the Table S2
    # script 'TC = TA * RTdT': Appendix S1 Eq 10 prints the baseline TC0 here,
    # but only the surviving-cell form reproduces the typical fits of Figure 3
    # (see the vignette).
    efree <- rbase_nk * nkp46_receptors * 1e6 / navo - drug_nkp46_central - trimer
    tfree <- cd123cell_central * cd123_receptors * 1e6 / navo - drug_cd123_central - trimer

    # Sequential binding (Appendix S1 Eqs 1-4). Eq 3 as printed has kon2*T*E as
    # its first term and Eq 4 a +koff1*TBE term; both are corrected to the
    # mass-balanced forms kon2*T*B and -koff1*TBE used in Table S2.
    d/dt(central) <- -kon1 * central * efree + koff1 * drug_nkp46_central -
      kon2 * central * tfree + koff2 * drug_cd123_central
    d/dt(drug_nkp46_central) <- kon1 * central * efree - koff1 * drug_nkp46_central -
      kon2 * tfree * drug_nkp46_central + koff2 * trimer
    d/dt(drug_cd123_central) <- kon2 * tfree * central - koff2 * drug_cd123_central -
      kon1 * drug_cd123_central * efree + koff1 * trimer
    d/dt(trimer) <- kon2 * tfree * drug_nkp46_central - koff2 * trimer +
      kon1 * drug_cd123_central * efree - koff1 * trimer

    # Synapse-driven killing rate (Appendix S1 Eq 8); synapse floored at 0 as
    # in Table S2 'TBE2 = max(0, TBE)'
    synapse <- max(0, trimer)
    kkill <- kmax * synapse^hill / (ec50^hill + synapse^hill)

    # Treated and untreated CD123+ cells (Appendix S1 Eqs 6-7)
    d/dt(cd123cell_central) <- kg * cd123cell_central - kkill * cd123cell_central
    d/dt(cd123cell_control) <- kg * cd123cell_control
    cd123cell_central(0) <- rbase_cd123
    cd123cell_control(0) <- rbase_cd123

    # Percent cytotoxicity (Appendix S1 Eq 5)
    cytotox <- (1 - cd123cell_central / cd123cell_control) * 100
    Cc <- central
  })
}
