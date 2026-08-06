Mody_2023_doxorubicin_dexrazoxane_static_mdamb468 <- function() {
  description <- "In vitro (MDA-MB-468 human breast cancer cell line; triple-negative breast cancer). Static 72-h concentration-response model for doxorubicin (DOX) and dexrazoxane (DEX) alone and in combination (Mody 2023 Eqs 1 and 2; Table 1, MDA-MB-468 columns). This is the ALGEBRAIC endpoint model, not a time-course model: it predicts percent cell viability at 72 h as a function of the two nominal applied concentrations, which enter as the covariates CONC_DOXORUBICIN and CONC_DEXRAZOXANE rather than as ODE states. The combination is described by the Chakraborty-Jusko / Pawaskar competitive interaction model (Eq 2), in which the interaction parameter psi scales the half-maximal inhibitory concentration of BOTH drugs; psi > 1 is antagonistic, psi < 1 synergistic and psi = 1 additive. Estimated psi is 0.84 (modestly synergistic) in MDA-MB-468. The single-agent Imax and IC50 values were FIXED from the Eq-1 single-agent concentration-response fits and only psi was estimated (Monolix 2016R1). Setting either concentration to zero and psi to 1 recovers the single-agent inhibitory Hill model of Eq 1. Sibling models: Mody_2023_doxorubicin_dexrazoxane_static_jimt1 (the other cell line) and Mody_2023_doxorubicin_dexrazoxane_mdamb468 (the time-course PD model fitted to the same cells)."
  reference <- paste(
    "Mody H, Vaidya TR, Lezeau J, Taha K, Ait-Oudhia S (2023).",
    "In vitro to clinical translation of combinatorial effects of doxorubicin",
    "and dexrazoxane in breast cancer: a mechanism-based pharmacokinetic/",
    "pharmacodynamic modeling approach.",
    "Frontiers in Pharmacology 14:1239141.",
    "doi:10.3389/fphar.2023.1239141.",
    sep = " "
  )
  vignette <- "Mody_2023_doxorubicin_dexrazoxane_breast_cancer"

  units <- list(
    time          = "h (the model is evaluated at the single 72 h assay endpoint)",
    dosing        = "not applicable (drug exposure enters as the concentration covariates)",
    concentration = "uM (both drug concentration covariates); % (cell viability, the PD readout)"
  )

  # CONC_DOXORUBICIN / CONC_DEXRAZOXANE are the experimentally applied,
  # nominal well concentrations of the 72-h static assay. They are
  # supplied as covariates from the event data, exactly as the in-vitro
  # precedent Landersdorfer_2013_nisin_amikacin.R supplies Cnis / Cami.
  depends <- c("CONC_DOXORUBICIN", "CONC_DEXRAZOXANE")

  compartmentData <- list()

  covariateData <- list(
    CONC_DOXORUBICIN = list(
      description        = "Nominal doxorubicin concentration applied to the well",
      units              = "uM",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Experimentally controlled assay input, constant over the 72 h exposure. Levels used: 0.005, 0.01, 0.05, 0.1, 0.5 and 1 uM (Figure 2A), crossed with the six DEX levels to give 36 combinations. In-vitro experimental input -- not in inst/references/covariate-columns.md, whose canonical register covers human population-PK covariates and does not apply to this in-vitro assay model (same treatment as Cnis / Cami in Landersdorfer_2013_nisin_amikacin.R).",
      source_name        = "DOX concentration (Methods, CCK-8 cell viability assay; Figure 2A)"
    ),
    CONC_DEXRAZOXANE = list(
      description        = "Nominal dexrazoxane concentration applied to the well",
      units              = "uM",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Experimentally controlled assay input, constant over the 72 h exposure. Levels used: 6.25, 12.5, 25, 50, 100 and 200 uM (Figure 2A). In-vitro experimental input -- not in inst/references/covariate-columns.md (see the CONC_DOXORUBICIN note).",
      source_name        = "DEX concentration (Methods, CCK-8 cell viability assay; Figure 2A)"
    )
  )

  population <- list(
    species        = "in vitro (MDA-MB-468 human breast cancer cell line)",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    disease_state  = "triple-negative breast cancer. MDA-MB-468 cells seeded at 5 x 10^3 cells per well (100 uL) of a 96-well plate, incubated overnight, then exposed to DOX, DEX or their combination and assayed for viability with CCK-8 (absorbance 450 nm).",
    dose_range     = "DOX 0.005-1 uM, DEX 6.25-200 uM, and 36 DOX x DEX combinations of those levels.",
    cell_line      = "MDA-MB-468 (ATCC, Manassas, VA); triple-negative breast cancer.",
    culture        = "DMEM with 10% fetal bovine serum and 1% penicillin/streptomycin; 37 C, humidified 5% CO2; passaged at confluency with 0.25% trypsin / 2.21 nM EDTA.",
    notes          = "Static 72-h endpoint assay, at least triplicate wells per condition against matched vehicle controls. Only psi was estimated in Eq 2; every other parameter was fixed from the single-agent Eq-1 fits."
  )

  ini({
    # ================================================================
    # SINGLE-AGENT CONCENTRATION-RESPONSE -- Mody 2023 Table 1
    # ================================================================
    # Eq 2's baseline is stated in the text as "R0 is % cell viability
    # at baseline (i.e., 100%)", so the competitive interaction model
    # uses 100 rather than the per-drug Eq-1 baselines (Table 1 reports
    # 100 (Fixed) for the DOX curve and 110 for the DEX curve).
    lrbase <- fixed(log(100)) ; label("Baseline percent cell viability R0 (%)")   # Eq 2 text: "R0 is % cell viability at baseline (i.e., 100%)"

    # Both drugs' Imax and IC50 were fixed into Eq 2 from the Eq-1
    # single-agent fits ("All individual drug-related parameters were
    # fixed from the DOX and DEX concentration-response curve fittings").
    limax_dox <- fixed(log(0.959))  ; label("Maximal fractional inhibition of viability by DOX, Imax,DOX (unitless)")  # Table 1: 0.959 (1.27% RSE)
    lec50_dox <- fixed(log(0.0211)) ; label("DOX concentration for half-maximal inhibition IC50,DOX (uM)")            # Table 1: 21.1 nM = 0.0211 uM (19.9% RSE)

    limax_dexrazoxane <- fixed(log(0.825))  ; label("Maximal fractional inhibition of viability by DEX, Imax,DEX (unitless)")  # Table 1: 0.825 (2.17% RSE)
    lec50_dexrazoxane <- fixed(log(36))  ; label("DEX concentration for half-maximal inhibition IC50,DEX (uM)")             # Table 1: 36 uM (16.2% RSE)

    # Eq 2 carries Hill coefficients gamma_A / gamma_B, but Table 1
    # reports none, and Eq 1 -- the source of the fixed single-agent
    # parameters -- has no Hill exponent. Both are therefore 1.
    lhill_dox         <- fixed(log(1)) ; label("Hill coefficient gamma_A for DOX (unitless)")   # not reported; Eq 1 has no Hill exponent, so 1
    lhill_dexrazoxane <- fixed(log(1)) ; label("Hill coefficient gamma_B for DEX (unitless)")   # not reported; Eq 1 has no Hill exponent, so 1

    # The only parameter estimated in Eq 2.
    lpsi <- log(0.84) ; label("DOX-DEX interaction parameter psi (unitless; >1 antagonistic, <1 synergistic, =1 additive)")  # Results / Figure 2 caption: psi = 0.84

    # Table 1 reports %RSE on point estimates only; no residual-error
    # magnitude is reported. Encoded as fixed(0) rather than invented.
    addSd_viability <- fixed(0) ; label("Additive residual SD on percent cell viability (%; ZERO - not reported in source)")
  })

  model({
    rbase             <- exp(lrbase)
    imax_dox          <- exp(limax_dox)
    ec50_dox          <- exp(lec50_dox)
    imax_dexrazoxane  <- exp(limax_dexrazoxane)
    ec50_dexrazoxane  <- exp(lec50_dexrazoxane)
    hill_dox          <- exp(lhill_dox)
    hill_dexrazoxane  <- exp(lhill_dexrazoxane)
    psi               <- exp(lpsi)

    # Competitive interaction model, Eq 2. Each drug contributes a
    # psi-scaled, Hill-transformed concentration ratio; the two ratios
    # compete in a shared denominator that also carries the +1 for the
    # drug-free state.
    ratio_dox         <- (CONC_DOXORUBICIN / (psi * ec50_dox))^hill_dox
    ratio_dexrazoxane <- (CONC_DEXRAZOXANE / (psi * ec50_dexrazoxane))^hill_dexrazoxane

    viability <- rbase *
      (1 - (imax_dox * ratio_dox + imax_dexrazoxane * ratio_dexrazoxane) /
             (ratio_dox + ratio_dexrazoxane + 1))

    viability ~ add(addSd_viability)
  })
}
