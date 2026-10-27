Gall_2023_fluorouracil_human <- function() {
  description <- paste(
    "QSP. Human small-intestinal epithelium. Gall 2023 translational model of",
    "5-fluorouracil (5-FU)-induced clinical intestinal toxicity. A human",
    "epithelial model with the same structure as the Gall 2023 mouse 5-FU",
    "model but parametrised for the human ileum (cell composition and",
    "proliferation rates from published human data) translates the cytotoxic",
    "activity of 5-FU quantified in mice into human epithelial injury and",
    "clinical diarrhea risk. The intracellular metabolic cascade",
    "(FUMP -> FUDP -> FUTP) is driven by the time-varying plasma 5-FU",
    "concentration CP_5FU_NGML (supplied externally; the authors simulated it",
    "with Simcyp for relevant clinical regimens), using metabolite rate",
    "constants assumed equal to the mouse values. FUTP drives, through a Hill",
    "function, cell-cycle arrest and apoptosis of proliferative crypt cells in",
    "the four-compartment epithelial model (stem, proliferative progenitors,",
    "non-proliferative crypt and villus cells) with villus-to-crypt feedback",
    "and arrested/apoptotic compartments; the epithelial state drives plasma",
    "citrulline and the probability of diarrhea. Deterministic typical-value",
    "model (no IIV, no residual error). Time is in hours. Equations S1-S23 of",
    "Gall 2023. The drug-effect parameters are assumed from the mouse fit."
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
    "fump",
    "fudp",
    "futp",
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
    dosing = "n/a (the plasma 5-FU exposure is supplied externally via the time-varying covariate CP_5FU_NGML; the model carries no 5-FU dosing compartment)",
    concentration = "ng/mL for metabolites and citrulline; cell counts (two-dimensional crypt-villus strip) for epithelial states"
  )

  covariateData <- list(
    CP_5FU_NGML = list(
      description = "Time-varying plasma 5-fluorouracil concentration driving the intracellular metabolite cascade.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Supplied by the user at each event record. Gall 2023 simulated it",
        "with Simcyp for clinical regimens (a fortnightly 46-hour infusion",
        "giving AUC 18 or 28 mg.h/L, and 450 mg daily for 5 days every 28",
        "days); the deposited profiles are in the BioModels supplement",
        "(Pk_Mets_Hu.mat). Interpolated linearly between records. Set to 0",
        "outside the exposure window."
      ),
      source_name = "FUplasma"
    )
  )

  compartmentData <- list(
    fump = list(analyte = "fluorouridine monophosphate (FUMP)", units = "ng/mL", specimen = "tissue", verified = TRUE),
    fudp = list(analyte = "fluorouridine diphosphate (FUDP)", units = "ng/mL", specimen = "tissue", verified = TRUE),
    futp = list(analyte = "fluorouridine triphosphate (FUTP)", units = "ng/mL", specimen = "tissue", verified = TRUE),
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
    disease_state = "Human small-intestinal epithelium under clinical 5-FU chemotherapy",
    dose_range = "Simulated clinical regimens: 46-hour intravenous 5-FU infusion every 2 weeks (AUC 18 or 28 mg.h/L); 450 mg intravenous bolus daily for 5 days every 28 days",
    regions = "Translational (human-ileum parameters from the literature)",
    notes = paste(
      "The human epithelial model keeps the structure of the mouse model but",
      "is parametrised for the human ileum: stem, proliferative and",
      "non-proliferative crypt cell numbers and proliferation rates from",
      "published human data, villus cell number scaled using the mouse",
      "villus/crypt ratio, and the stem-cell cycle scaled by the mouse",
      "stem/proliferative cycle ratio. The feedback strength, migration shape",
      "coefficient, drug-effect parameters and diarrhea logistic parameters",
      "are assumed from the mouse fit; metabolite rate constants are assumed",
      "equal to the mouse values. Citrulline clearance uses a human half-life",
      "(0.65 h); the production rate is derived at the human citrulline",
      "steady state. No inter-individual variability and no residual error:",
      "the model is deterministic."
    )
  )

  ini({
    # ==================================================================
    # 5-FU intracellular metabolism -- metabolite rate constants assumed
    # equal to the Gall 2023 mouse values (Table 2), applied to the
    # externally supplied human plasma 5-FU concentration.
    # ==================================================================
    kfump <- fixed(1.555)
    label("5-FU conversion rate into FUMP KFUMP (1/h)") # Table 2, mouse value assumed for human
    kfudp <- fixed(0.2973)
    label("FUMP conversion rate into FUDP KFUDP (1/h)") # Table 2, mouse value assumed for human
    kfutp <- fixed(0.8478)
    label("FUDP conversion rate into FUTP KFUTP (1/h)") # Table 2, mouse value assumed for human
    kefutp <- fixed(1.010)
    label("FUTP clearance rate KeFUTP (1/h)") # Table 2, mouse value assumed for human

    # ==================================================================
    # Human epithelial cell dynamics -- Gall 2023 Table 1, 'Human' column.
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
    # Drug-effect parameters -- assumed from the Gall 2023 mouse 5-FU fit
    # (Table 1).
    # ==================================================================
    kinhib <- fixed(0.793992)
    label("5-FU proliferation-inhibition (arrest) strength k_inhib (unitless)") # Table 1: 0.7940, assumed from mouse
    kkill <- fixed(0.55509)
    label("5-FU cell-death (apoptosis) strength k_kill (unitless)") # Table 1: 0.5551, assumed from mouse
    c50 <- fixed(1317.88)
    label("FUTP concentration for half-maximal effect C50 (ng/mL)") # Table 1: 1318, assumed from mouse
    hill <- fixed(5)
    label("Hill coefficient of the FUTP effect (unitless)") # Table 1: 5, assumed from mouse

    # ==================================================================
    # Plasma citrulline -- Gall 2023 Table 1, 'Human'. Clearance uses a human
    # citrulline half-life of 0.65 h [56]; the production rate is derived at
    # the human steady state in model().
    # ==================================================================
    citss <- fixed(7175.064)
    label("Steady-state plasma citrulline concentration citru_ss (ng/mL)") # Table 1: 7175.064, human value [18]
    clcitru <- fixed(1.066380278)
    label("Plasma citrulline elimination rate Cl_Citru (1/h)") # Table 1: 1.066380278, human value [56]

    # ==================================================================
    # Probability of diarrhea -- assumed from the Gall 2023 mouse fit.
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

    # --- Intracellular 5-FU metabolism driven by external plasma 5-FU
    # (Eqs S3-S5). ------------------------------------------------------
    d/dt(fump) <- kfump * CP_5FU_NGML - kfudp * fump
    d/dt(fudp) <- kfudp * fump - kfutp * fudp
    d/dt(futp) <- kfutp * fudp - kefutp * futp

    # --- Drug effects on proliferative cells (Eqs S19-S20). -------------
    emax_futp <- futp^hill / (futp^hill + c50^hill)
    arrest <- kinhib * emax_futp
    kill <- kkill * emax_futp

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
