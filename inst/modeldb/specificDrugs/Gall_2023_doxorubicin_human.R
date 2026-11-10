Gall_2023_doxorubicin_human <- function() {
  description <- paste(
    "QSP. Human small-intestinal epithelium. Gall 2023 translational model of",
    "doxorubicin-induced clinical intestinal toxicity. A human epithelial",
    "model with the same structure as the Gall 2023 human 5-FU model",
    "translates the doxorubicin cytotoxic activity quantified in mice into",
    "human epithelial injury and clinical diarrhea risk. Cell-cycle arrest",
    "and apoptosis of proliferative crypt cells are driven by a linear",
    "function of the unbound doxorubicin plasma concentration, computed from",
    "the time-varying total plasma doxorubicin concentration",
    "CP_DOXORUBICIN_NGML (supplied externally; the authors simulated it with",
    "Simcyp for a 75 mg/m2 clinical regimen) and the human unbound fraction.",
    "The four-compartment epithelial model (stem, proliferative progenitors,",
    "non-proliferative crypt and villus cells) with villus-to-crypt feedback",
    "and arrested/apoptotic compartments drives plasma citrulline and the",
    "probability of diarrhea. Deterministic typical-value model (no IIV, no",
    "residual error). Time is in hours. Equations S8-S23 of Gall 2023. The",
    "drug-effect parameters are assumed from the mouse fit, so this linear",
    "model should not be extrapolated beyond the exposures studied."
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
    dosing = "n/a (the plasma doxorubicin exposure is supplied externally via the time-varying covariate CP_DOXORUBICIN_NGML; the model carries no doxorubicin dosing compartment)",
    concentration = "ng/mL for citrulline; cell counts (two-dimensional crypt-villus strip) for epithelial states"
  )

  covariateData <- list(
    CP_DOXORUBICIN_NGML = list(
      description = "Time-varying total plasma doxorubicin concentration; multiplied by the unbound fraction to drive the cytotoxic effect.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Supplied by the user at each event record. Gall 2023 simulated it",
        "with Simcyp for 75 mg/m2 intravenous doxorubicin every 3 weeks; the",
        "deposited profile is in the BioModels supplement (Pk_Doxo_Hu.mat).",
        "Interpolated linearly between records. Set to 0 outside the exposure",
        "window. The model multiplies it by the human unbound fraction fu to",
        "obtain the unbound concentration driving the effect."
      ),
      source_name = "Doxongml"
    )
  )

  compartmentData <- list(
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

  population <- list(
    species = "human (ileum epithelium)",
    n_subjects = "Typical-value translational model; parametrised from published human ileal cell-composition and proliferation data, not fitted to individual patients.",
    n_studies = NA_integer_,
    disease_state = "Human small-intestinal epithelium under clinical doxorubicin chemotherapy",
    dose_range = "Simulated clinical regimen: 75 mg/m2 intravenous doxorubicin every 3 weeks",
    regions = "Translational (human-ileum parameters from the literature)",
    notes = paste(
      "The human epithelial model shares the structure and human-ileum",
      "parameters of the Gall 2023 human 5-FU model. The doxorubicin drug",
      "effect is linear in the unbound plasma concentration; its potency",
      "constants are assumed from the mouse doxorubicin fit. The citrulline",
      "and diarrhea models are shared with the human 5-FU model. No",
      "inter-individual variability and no residual error: the model is",
      "deterministic."
    )
  )

  ini({
    # ==================================================================
    # Human epithelial cell dynamics -- Gall 2023 Table 1, 'Human' column
    # (shared with the human 5-FU model).
    # ==================================================================
    lambda_s <- fixed(log(2) / (21 / 12 * 29.7))
    label("Stem-cell proliferation rate lambda_s (1/h)") # Table 1: ln(2)/(21/12 x 29.7), scaled using mouse stem/proliferative cycle ratio
    lambda_p <- fixed(log(2) / 29.7)
    label("Proliferative-cell proliferation rate lambda_P (1/h)") # Table 1: ln(2)/29.7, human value [30]
    tau_apo <- fixed(1 / 12)
    label("Apoptotic-cell detachment rate tau_Apo (1/h)") # Table 1: 1/12, literature value [55]
    gam <- fixed(3.41394)
    label("Villus-to-crypt proliferation feedback strength gamma (unitless)") # Table 1: 3.41, assumed from mouse
    mm <- fixed(1.468566)
    label("Crypt-to-villus migration shape coefficient m (unitless)") # Table 1: 1.469, assumed from mouse
    sss <- fixed(2)
    label("Steady-state stem-cell number s_ss (cells)") # Table 1: 2, human value [25]
    pss <- fixed(19.446)
    label("Steady-state proliferative-progenitor number p_ss (cells)") # Table 1: 19.446, human value [30]
    dss <- fixed(11.154)
    label("Steady-state non-proliferative crypt-cell number d_ss (cells)") # Table 1: 11.154, human value [30]
    vss <- fixed(4 * 32.6)
    label("Steady-state villus-cell number v_ss (cells)") # Table 1: 4 x 32.6, scaled using the mouse villus/crypt ratio

    # ==================================================================
    # Doxorubicin drug effect -- assumed from the Gall 2023 mouse
    # doxorubicin fit (Table 1). Linear in the unbound concentration.
    # ==================================================================
    fu <- fixed(0.28)
    label("Fraction of unbound doxorubicin in plasma Fu (unitless)") # Table 1: 0.28 (Simcyp, human)
    kinhib <- fixed(0.002966)
    label("Doxorubicin proliferation-inhibition (arrest) strength k_inhib (mL/ng)") # Table 1: 0.002966, assumed from mouse
    kkill <- fixed(0.00214)
    label("Doxorubicin cell-death (apoptosis) strength k_kill (mL/ng)") # Table 1: 0.00214, assumed from mouse

    # ==================================================================
    # Plasma citrulline -- shared with the human 5-FU model.
    # ==================================================================
    citss <- fixed(7175.064)
    label("Steady-state plasma citrulline concentration citru_ss (ng/mL)") # Table 1: 7175.064, human value [18]
    clcitru <- fixed(1.066380278)
    label("Plasma citrulline elimination rate Cl_Citru (1/h)") # Table 1: 1.066380278, human value [56]

    # ==================================================================
    # Probability of diarrhea -- shared with the human 5-FU model.
    # ==================================================================
    b0 <- fixed(3.5)
    label("Diarrhea logistic intercept beta_0 (unitless)") # Table 1: 3.5, assumed from mouse
    b1 <- fixed(-5.1298)
    label("Diarrhea logistic slope beta_1 (unitless)") # Table 1: -5.1298, assumed from mouse
    loss50 <- fixed(0.069511)
    label("Villus-cell-loss value for half-maximal diarrhea odds Loss50 (unitless)") # Table 1: 0.06951, assumed from mouse
  })
  model({
    # --- Steady-state-derived epithelial transfer and shedding rates. ---
    delta_sp <- lambda_s
    delta_pd <- lambda_p * (1 + gam / 2) + delta_sp * sss / pss
    delta_dv <- delta_pd * pss / dss
    tau_v <- delta_dv * dss / vss
    kcitru <- clcitru * citss / vss^2

    # --- Drug effect: linear in unbound doxorubicin concentration
    # (Eqs S19-S20, doxorubicin branch). --------------------------------
    doxo_unbound <- fu * CP_DOXORUBICIN_NGML
    arrest <- kinhib * doxo_unbound
    kill <- kkill * doxo_unbound

    # --- Human epithelial cell dynamics (Eqs S10-S18). ------------------
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
    vtot <- villus + varr + vapo
    cellloss <- max(0, 1 - vtot / vss)
    prob_diarrhea <- 1 - expit(b0 + b1 * cellloss / (cellloss + loss50))
  })
}
