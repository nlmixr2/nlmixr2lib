Gall_2023_doxorubicin_mouse <- function() {
  description <- paste(
    "QSP. Preclinical (mouse, C57BL/6J small-intestinal epithelium). Gall 2023",
    "dynamic model of doxorubicin-induced intestinal toxicity fitted to mouse",
    "data. A concentration-based two-compartment doxorubicin plasma PK model",
    "(fitted to Simcyp mouse predictions) drives cell-cycle arrest and",
    "apoptosis of proliferative crypt cells through a linear function of the",
    "unbound doxorubicin plasma concentration (the Hill saturation used for",
    "5-FU is not identifiable at the mild doxorubicin exposures studied). The",
    "four-compartment epithelial model (stem, proliferative progenitors,",
    "non-proliferative crypt and villus cells) with villus-to-crypt feedback",
    "and arrested/apoptotic crypt and villus compartments also drives plasma",
    "citrulline and the probability of diarrhea. Deterministic typical-value",
    "model (no IIV, no residual error). Time is in hours; the plasma states",
    "are concentrations, so an intravenous dose in mg/kg is converted to its",
    "instantaneous plasma-concentration increment through f(central).",
    "Equations S8-S23 of Gall 2023. This linear doxorubicin-effect model",
    "should not be extrapolated to doses higher than those studied."
  )
  reference <- paste(
    "Gall L, Jardi F, Lammens L, Pinero J, Souza TM, Rodrigues D, Jennen DGJ,",
    "de Kok TM, Coyle L, Chung SW, Ferreira S, Jo H, Beattie KA, Kelly C,",
    "Duckworth CA, Pritchard DM, Pin C. (2023).",
    "A dynamic model of the intestinal epithelium integrates multiple sources",
    "of preclinical data and enables clinical translation of drug-induced",
    "toxicity. CPT Pharmacometrics Syst Pharmacol 12(10):1511-1528.",
    "doi:10.1002/psp4.13029.",
    "Numeric values are from the authors' MATLAB code deposited at BioModels",
    "(MODEL2212120003), which generated the published figures and agrees with",
    "Gall 2023 Tables 1 and 2 to the precision printed there."
  )
  vignette <- "Gall_2023_intestinal_epithelium_toxicity"

  # Mechanistic epithelial cell compartments of the Gall 2023 model with no
  # canonical analogue in inst/references/compartment-names.md.
  paper_specific_compartments <- c(
    "stem",
    "prolif",
    "nonprolif",
    "villus",
    "carr",
    "capo",
    "varr",
    "vapo",
    "citr"
  )

  units <- list(
    time = "h",
    dosing = "mg/kg (intravenous; converted to a plasma-concentration increment through f(central))",
    concentration = "ng/mL for doxorubicin plasma and citrulline; cell counts (two-dimensional crypt-villus strip) for epithelial states"
  )

  compartmentData <- list(
    central = list(analyte = "doxorubicin", units = "ng/mL", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "doxorubicin", units = "ng/mL", specimen = "tissue", verified = TRUE),
    stem = list(analyte = "stem cells", units = "cells", specimen = "tissue", verified = TRUE),
    prolif = list(
      analyte = "proliferative crypt progenitor cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    nonprolif = list(
      analyte = "non-proliferative (differentiated) crypt cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    villus = list(analyte = "healthy villus epithelial cells", units = "cells", specimen = "tissue", verified = TRUE),
    carr = list(analyte = "drug-arrested crypt cells", units = "cells", specimen = "tissue", verified = TRUE),
    capo = list(analyte = "apoptotic crypt cells", units = "cells", specimen = "tissue", verified = TRUE),
    varr = list(analyte = "drug-arrested villus cells", units = "cells", specimen = "tissue", verified = TRUE),
    vapo = list(analyte = "apoptotic villus cells", units = "cells", specimen = "tissue", verified = TRUE),
    citr = list(analyte = "citrulline", units = "ng/mL", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "mouse (C57BL/6J; small-intestinal epithelium)",
    n_subjects = "Literature group-mean data; no individual animals were modelled (deterministic typical-value fit).",
    n_studies = 1L,
    sex_female_pct = 0,
    disease_state = "Healthy small-intestinal epithelium challenged with doxorubicin",
    dose_range = "5 or 10 mg/kg intravenous bolus every 24 h for 2 days, then a 2-day recovery period (Gall 2023, newly generated mouse data)",
    regions = "Janssen Pharmaceutica, Belgium (preclinical)",
    notes = paste(
      "C57BL/6J male mice received doxorubicin 5 or 10 mg/kg by intravenous",
      "bolus on days 0 and 1. The epithelial cell counts and plasma citrulline",
      "were newly generated in Gall 2023; doxorubicin plasma PK was fitted to",
      "Simcyp mouse predictions. The epithelial, feedback and shape parameters",
      "are shared with the Gall 2023 5-FU mouse model; only the two",
      "doxorubicin potency constants, the doxorubicin-specific steady-state",
      "cell numbers and citrulline parameters, and the doxorubicin PK differ.",
      "The citrulline and diarrhea models were fitted to 5-FU data and used",
      "unchanged to predict doxorubicin responses. No inter-individual",
      "variability and no residual error: the model is deterministic."
    )
  )

  ini({
    # ==================================================================
    # Doxorubicin plasma PK -- Gall 2023 Table 2, 'Doxorubicin PK model'
    # (posterior means). The plasma states are concentrations (ng/mL); a
    # mg/kg intravenous dose is scaled to its ng/mL plasma increment through
    # f(central) = kadj.
    # ==================================================================
    kadj <- 1045.473
    label("Dose-to-concentration adjustment constant Ka (ng x kg/mg/mL)") # Table 2, Ka = 1045 +/- 23.49
    kedoxo <- 0.9746
    label("Doxorubicin plasma elimination rate KeDoxo (1/h)") # Table 2, KeDoxo = 0.9746 +/- 0.0196
    kper <- 2.138
    label("Doxorubicin plasma-to-peripheral transfer rate Kper (1/h)") # Table 2, Kper = 2.138 +/- 0.0453
    kfromper <- 0.2347
    label("Doxorubicin peripheral-to-plasma transfer rate Kfromper (1/h)") # Table 2, Kfromper = 0.2347 +/- 0.0040
    fu <- fixed(0.22)
    label("Fraction of unbound doxorubicin in plasma Fu (unitless)") # Table 1: 0.22 (Simcyp)

    # ==================================================================
    # Epithelial cell dynamics -- Gall 2023 Table 1, 'Mouse Doxorubicin'.
    # Proliferation and apoptotic-detachment rates are literature values;
    # steady-state cell numbers are direct data estimates. The feedback
    # strength and shape coefficient are shared with the 5-FU mouse fit.
    # ==================================================================
    lambda_s <- fixed(log(2) / 21)
    label("Stem-cell proliferation rate lambda_s (1/h)") # Table 1: ln(2)/21, literature value [54]
    lambda_p <- fixed(log(2) / 12)
    label("Proliferative-cell proliferation rate lambda_P (1/h)") # Table 1: ln(2)/12, literature value [9]
    tau_apo <- fixed(1 / 12)
    label("Apoptotic-cell detachment rate tau_Apo (1/h)") # Table 1: 1/12, literature value [55]
    gam <- 3.41394
    label("Villus-to-crypt proliferation feedback strength gamma (unitless)") # Table 1: 3.41 +/- 0.0793 (shared with 5-FU fit)
    mm <- 1.4685
    label("Crypt-to-villus migration shape coefficient m (unitless)") # Table 1: 1.469 +/- 0.0418 (shared with 5-FU fit)
    sss <- fixed(2)
    label("Steady-state stem-cell number s_ss (cells)") # Table 1: 2, literature value [8]
    pss <- fixed(14.6978)
    label("Steady-state proliferative-progenitor number p_ss (cells)") # Table 1: 14.70, direct data estimate (Ki67+ minus s_ss)
    dss <- fixed(7.1063)
    label("Steady-state non-proliferative crypt-cell number d_ss (cells)") # Table 1: 7.11, direct data estimate (crypt minus p_ss minus s_ss)
    vss <- 96.35
    label("Steady-state villus-cell number v_ss (cells)") # Table 1: 96.35 +/- 5.395

    # ==================================================================
    # Doxorubicin drug effect -- Gall 2023 Table 1, 'Mouse Doxorubicin'.
    # Linear in the unbound doxorubicin plasma concentration.
    # ==================================================================
    kinhib <- 0.002966
    label("Doxorubicin proliferation-inhibition (arrest) strength k_inhib (mL/ng)") # Table 1: 0.002966 +/- 0.000613
    kkill <- 0.00214
    label("Doxorubicin cell-death (apoptosis) strength k_kill (mL/ng)") # Table 1: 0.00214 +/- 0.000478

    # ==================================================================
    # Plasma citrulline -- Gall 2023 Table 1. kcitru is fitted; the clearance
    # Cl_Citru is derived at steady state in model() as kcitru v_ss^2/citru_ss.
    # ==================================================================
    citss <- 6600.072
    label("Steady-state plasma citrulline concentration citru_ss (ng/mL)") # Table 1: 6600 +/- 1061
    kcitru <- 0.022521
    label("Citrulline production rate k_Citru (ng/mL/h/cells^2)") # Table 1: 0.0232 +/- 0.00222 (code 0.022521; see vignette Errata)

    # ==================================================================
    # Probability of diarrhea -- Gall 2023 Table 1 logistic regression on
    # fractional villus-cell loss (fitted to 5-FU data, shared here).
    # ==================================================================
    b0 <- fixed(3.5)
    label("Diarrhea logistic intercept beta_0 (unitless)") # Table 1: 3.5 (fixed to achieve placebo probability ~0)
    b1 <- -5.1298
    label("Diarrhea logistic slope beta_1 (unitless)") # Table 1: -5.1298 +/- 0.8558
    loss50 <- 0.069511
    label("Villus-cell-loss value for half-maximal diarrhea odds Loss50 (unitless)") # Table 1: 0.06951 +/- 0.0421
  })
  model({
    # --- Steady-state-derived epithelial transfer and shedding rates. ---
    delta_sp <- lambda_s
    delta_pd <- lambda_p * (1 + gam / 2) + delta_sp * sss / pss
    delta_dv <- delta_pd * pss / dss
    tau_v <- delta_dv * dss / vss
    clcitru <- kcitru * vss^2 / citss

    # --- Doxorubicin plasma PK (Eqs S8-S9). -----------------------------
    d/dt(central) <- -(kedoxo + kper) * central + kfromper * peripheral1
    d/dt(peripheral1) <- kper * central - kfromper * peripheral1
    f(central) <- kadj

    # --- Drug effects on proliferative cells: linear in unbound
    # doxorubicin concentration (Eqs S19-S20, doxorubicin branch). --------
    doxo_unbound <- fu * central
    arrest <- kinhib * doxo_unbound
    kill <- kkill * doxo_unbound

    # --- Epithelial cell dynamics (Eqs S10-S18). ------------------------
    feedback <- gam * vss / (vss + villus)
    migration <- ((stem + prolif) / (sss + pss))^mm
    d/dt(stem) <- lambda_s * stem * (1 - arrest - kill) - delta_sp * stem * stem / sss
    d/dt(prolif) <- delta_sp * stem * stem / sss +
      lambda_p * prolif * (1 + feedback - arrest - kill) - delta_pd * prolif
    d/dt(nonprolif) <- delta_pd * prolif - delta_dv * nonprolif * migration
    d/dt(villus) <- delta_dv * nonprolif * migration - tau_v * villus
    d/dt(carr) <- arrest * (lambda_s * stem + lambda_p * prolif) - delta_dv * carr * migration
    d/dt(capo) <- kill * (lambda_s * stem + lambda_p * prolif) -
      delta_dv * capo * migration - tau_apo * capo
    d/dt(varr) <- delta_dv * carr * migration - tau_v * varr
    d/dt(vapo) <- delta_dv * capo * migration - tau_v * vapo - tau_apo * vapo

    # --- Plasma citrulline (Eq S21). ------------------------------------
    d/dt(citr) <- kcitru * villus^2 - clcitru * citr

    stem(0) <- sss
    prolif(0) <- pss
    nonprolif(0) <- dss
    villus(0) <- vss
    citr(0) <- citss

    # --- Clinical-readout observables. ----------------------------------
    Cc <- central
    vtot <- villus + varr + vapo
    cellloss <- max(0, 1 - vtot / vss)
    prob_diarrhea <- 1 - expit(b0 + b1 * cellloss / (cellloss + loss50))
  })
}
