SanchezDengra_2021_caffeine_rat_pbpk_mdck <- function() {
  description <- paste(
    "PBPK (semi-physiological neuro-PBPK). Preclinical (rat). Caffeine",
    "disposition in plasma, brain tissue and cerebrospinal fluid (CSF)",
    "after a constant-rate intravenous infusion into plasma (Equation 9).",
    "One plasma compartment with first-order elimination exchanges",
    "unbound drug with the brain across the blood-brain barrier (BBB) and",
    "with the CSF across the blood-CSF barrier (BCSFB); unbound brain",
    "drug also moves to the CSF with the interstitial bulk flow, and CSF",
    "drains back to plasma with the CSF sink flow. The barrier clearances",
    "are the in vitro MDCK apparent permeabilities (Table 2) times the",
    "rat barrier surface areas, rescaled by two fitted scaling factors",
    "(SC1 influx, SC2 efflux); a third fitted factor (SC3) rescales the",
    "in vitro pig-brain-homogenate unbound fraction in brain. Vd, kel and",
    "SC1-SC3 were fitted in Berkeley Madonna to literature mean plasma",
    "and brain profiles; the model is deterministic (no between-animal",
    "variability or residual error was estimated). One of three cell-line",
    "parameterisations (MDCK, MDCK-MDR1, hCMEC/D3) fitted for this drug;",
    "this file is the MDCK one."
  )
  reference <- paste(
    "Sanchez-Dengra B, Gonzalez-Alvarez I, Bermejo M, Gonzalez-Alvarez M.",
    "Physiologically Based Pharmacokinetic (PBPK) Modeling for Predicting",
    "Brain Levels of Drug in Rat. Pharmaceutics. 2021;13(9):1402.",
    "doi:10.3390/pharmaceutics13091402. Equations 4-14 and Tables 1-3 of",
    "the main text; logP for the QSPR scaling-factor sub-model is",
    "Supplementary Table S2."
  )
  vignette <- "SanchezDengra_2021_neuro_pbpk_rat"

  units <- list(time = "s", dosing = "ng", concentration = "ng/mL")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "caffeine", units = "ng", specimen = "plasma", verified = TRUE),
    brain = list(analyte = "caffeine", units = "ng", specimen = "tissue", verified = TRUE),
    csf = list(analyte = "caffeine", units = "ng", specimen = "CSF", verified = TRUE)
  )

  population <- list(
    species = "rat",
    n_subjects = NA_integer_,
    n_studies = 1L,
    weight_range = "300 g (Table 1)",
    disease_state = "healthy",
    dose_range = "constant-rate intravenous infusion at k0 = 833.333 ng/s (Table 1 k0; Equation 9)",
    in_vitro_system = "MDCK (Madin-Darby canine kidney) transwell monolayers; fu,brain from pig brain homogenate",
    notes = paste(
      "The caffeine plasma and brain concentration-time profiles were",
      "taken from the published literature (Methods 2.4, refs 16-20) and",
      "fitted as mean profiles, so no number of animals is reported. The",
      "brain data are the unbound brain concentration (Figure 2B,",
      "Cb,u_exp) and the plasma data are the total plasma concentration",
      "(Figure 2B, Cp_exp). Rat weight in the source study: 300 g (Table",
      "1). In vitro inputs are from the MDCK (Madin-Darby canine kidney)",
      "transwell assays with pig brain homogenate (Methods 2.2)."
    )
  )

  ini({
    # ---- Caffeine-specific in vivo parameters, fitted (Table 3) ----
    # Table 1 lists the initial estimates (marked '*', later adjusted); Table 3
    # holds the final fitted values, shared by all three cell-line fits.
    lvc <- log(273.6); label("Apparent volume of distribution of the plasma compartment Vd (mL)") # Table 3 Caffeine row: Vd = 273.6 cm^3
    lkel <- log(3.55e-5); label("First-order elimination rate constant from plasma kel (1/s)") # Table 3 Caffeine row: kel = 3.55 x 10^-5 s^-1

    # ---- Scaling factors for the MDCK cell line, fitted (Table 3) ----
    # All three were initialised at 1 (Results, paragraph after Figure 2).
    sc1 <- 3.85; label("Scaling factor SC1 on the in vitro influx permeability (unitless)") # Table 3 Caffeine row, MDCK cell line, SC1 = 3.85
    sc2 <- 1.00; label("Scaling factor SC2 on the in vitro efflux permeability (unitless)") # Table 3 Caffeine row, MDCK cell line, SC2 = 1.00
    sc3 <- 0.22; label("Scaling factor SC3 on the in vitro unbound fraction in brain (unitless)") # Table 3 Caffeine row, MDCK cell line, SC3 = 0.22

    # ---- MDCK in vitro inputs, held fixed in the fit (Table 2) ----
    papp_ab <- fixed(26.10); label("In vitro apparent permeability apical-to-basolateral, Papp A-B (x 1e-6 cm/s)") # Table 2 Caffeine row, MDCK cell line: Papp A-B = 26.10 x 10^-6 cm/s (measured in this work)
    papp_ba <- fixed(35.31); label("In vitro apparent permeability basolateral-to-apical, Papp B-A (x 1e-6 cm/s)") # Table 2 Caffeine row, MDCK cell line: Papp B-A = 35.31 x 10^-6 cm/s (measured in this work)
    fu_brain <- fixed(0.857); label("In vitro unbound fraction in brain, fu,brain (unitless)") # Table 2 Caffeine row, MDCK cell line: fu,brain = 0.857 (measured in this work)
    fu <- fixed(0.917); label("Unbound fraction in plasma, fu,plasma (unitless)") # Table 1 fu,plasma 0.917 (footnote b, ref 24)

    # ---- Rat CNS physiology, fixed and identical for all drugs ----
    # Methods 2.4, from Ball et al. (ref 14) and Engelhard et al. (ref 15).
    # The flows are printed in cm^3/s; they are used as printed because the
    # fitted profiles of Figure 2 are reproduced with them (see vignette).
    lvbrain <- fixed(log(1.28)); label("Brain volume Vb (mL)") # Methods 2.4: Vb = 1.28 cm^3 (ref 14)
    lvcsf <- fixed(log(0.25)); label("CSF volume VCSF (mL)") # Methods 2.4: VCSF = 0.25 cm^3 (ref 14)
    qbulk <- fixed(0.012); label("Brain interstitial bulk flow to CSF Qbulk (mL/s)") # Methods 2.4: Qbulk = 0.012 cm^3/s (ref 14)
    qsink <- fixed(0.132); label("CSF sink (drainage) flow back to plasma Qsink (mL/s)") # Methods 2.4: Qsink = 0.132 cm^3/s (ref 14)
    s_bbb <- fixed(187.5); label("Blood-brain barrier surface area SBBB (cm^2)") # Methods 2.4: SBBB = 187.5 cm^2 (ref 14)
    s_bcsfb <- fixed(0.0375); label("Blood-CSF barrier surface area SBCSFB (cm^2)") # Methods 2.4: SBCSFB = 0.0375 cm^2 (ref 15)
  })
  model({
    vc <- exp(lvc)
    kel <- exp(lkel)
    vbrain <- exp(lvbrain)
    vcsf <- exp(lvcsf)

    # Barrier permeability-surface-area products (mL/s), Equations 10-13.
    # Papp is carried in units of 1e-6 cm/s (Table 2 header).
    ps_bbb_in <- sc1 * papp_ab * 1e-6 * s_bbb
    ps_bcsfb_in <- sc1 * papp_ab * 1e-6 * s_bcsfb
    ps_bbb_out <- sc2 * papp_ba * 1e-6 * s_bbb
    ps_bcsfb_out <- sc2 * papp_ba * 1e-6 * s_bcsfb

    # Concentrations (ng/mL). Figure 1 and Equation 14: unbound plasma is
    # fu,plasma * Cp, unbound brain is SC3 * fu,brain * Cb, and all CSF drug
    # is unbound.
    Cc <- central / vc
    Cu <- fu * Cc
    Cbrain <- brain / vbrain
    Cu_brain <- sc3 * fu_brain * Cbrain
    Ccsf <- csf / vcsf

    # Equations 4-8 written for amounts (state = volume x concentration).
    # Equation 9 (the zero-order input k0 is supplied as an infusion event
    # into central rather than written into the ODE)
    d/dt(central) <- - ps_bbb_in * Cu + ps_bbb_out * Cu_brain -
      ps_bcsfb_in * Cu + ps_bcsfb_out * Ccsf + qsink * Ccsf - kel * central
    # Equation 5
    d/dt(brain) <- ps_bbb_in * Cu - ps_bbb_out * Cu_brain - qbulk * Cu_brain
    # Equation 6
    d/dt(csf) <- ps_bcsfb_in * Cu - ps_bcsfb_out * Ccsf - qsink * Ccsf + qbulk * Cu_brain
  })
}
