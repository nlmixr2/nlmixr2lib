Sun_2026_tilmicosin_pkpd_pif <- function() {
  description <- paste0(
    "Veterinary (pig). Hollow-fiber infection model (HFIM) sigmoid Imax PK/PD-integration model ",
    "for the antibacterial effect of tilmicosin against Pasteurella multocida (challenge isolate ",
    "ZJWZ-A) under simulated PULMONARY INTERSTITIAL FLUID (PIF) exposure in swine. PIF is the ",
    "site of action for porcine respiratory pathogens, and Sun 2026 reports the first collection ",
    "of PIF from conscious pigs by microdialysis. Sun 2026 eq 8 parameterises the drug effect ",
    "over a 24 h window as E = E0 - (Imax * C^gamma) / (C^gamma + IC50^gamma), where E is the ",
    "SIGNED log10(CFU/mL) change from the initial inoculum (positive = net growth, negative = net ",
    "kill), E0 is that same change in the drug-free control, and C is the PK/PD index ",
    "AUC0-24/MIC. Parameters from Sun 2026 Table 3, column 'PIF': Imax = 6.82 log10 CFU/mL, ",
    "IC50 = 16.06 h, E0 = 2.36 log10 CFU/mL, gamma = 2.27. The parameterisation was confirmed ",
    "numerically against the paper's own independently printed target rows: solving E = 0, E = -3 ",
    "and E = -4 returns AUC0-24/MIC = 12.133, 28.481 and 51.081 h versus the published 12.13, ",
    "28.48 and 51.08 h, i.e. agreement to better than 0.03%. The bacterial density bact (linear ",
    "CFU/mL) is integrated as d/dt(bact) = ln(10) * (E / 24) * bact so that log10(bact) changes ",
    "by exactly E across each 24 h window, reproducing the paper's model at the times bacteria ",
    "were counted. There is NO PK component in this file: exposure enters as the externally ",
    "supplied covariate AUCMIC_TILM, which can be generated from the companion whole-body PBPK ",
    "model Sun_2026_tilmicosin_pbpk. That covariate carries the AUC0-24/MIC RATIO directly rather ",
    "than an absolute AUC divided by a model MIC parameter, because the reported isolate MIC of ",
    "8 ug/mL does not reconstruct the paper's own index column; see the vignette Errata. Neither ",
    "between-subject variability nor a residual error magnitude was reported (the +/- values in ",
    "Table 3 are standard errors of the point estimates, not estimated variance components), so ",
    "there are no eta parameters and addSd is FIXED at 0; the model is intended for typical-value ",
    "simulation. Compared with the plasma fit in the sibling model ",
    "Sun_2026_tilmicosin_pkpd_plasma, the PIF fit is far steeper (gamma 2.27 vs 0.46) and reaches ",
    "bactericidal and eradication effects at substantially lower exposure, which is the paper's ",
    "basis for concluding that tilmicosin performs better at the infection site than in plasma."
  )
  reference <- paste(
    "Sun L, Zhang C, Mi K, Wang H, Pan Y, Tao Y, Huang L.",
    "Dose optimization of tilmicosin against Pasteurella multocida in swine",
    "by physiologically based pharmacokinetic-pharmacodynamic model.",
    "J Agric Food Chem. 2026;74(8):4754-4766.",
    "doi:10.1021/acs.jafc.5c11368.",
    "Model equation from eq 8; parameter values from Table 3, column 'PIF';",
    "PK/PD index data from Supporting Information Table S10.",
    sep = " "
  )
  vignette <- "Sun_2026_tilmicosin"

  units <- list(
    time = "h",
    dosing = "h (tilmicosin AUC0-24/MIC index, supplied as a covariate)",
    concentration = "log10 CFU/mL (observation)"
  )

  depends <- c("AUCMIC_TILM")
  paper_specific_compartments <- c("bact")

  compartmentData <- list(
    bact = list(analyte = "Pasteurella multocida", units = "CFU/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    AUCMIC_TILM = list(
      description        = "Tilmicosin PK/PD index: pulmonary interstitial fluid area under the concentration-time curve over a 24 h dosing interval divided by the MIC of the challenge isolate (AUC0-24/MIC)",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Sun 2026 eq 8 defines the sigmoid Imax driver C as the AUC0-24/MIC ratio, and Table 3 ",
        "reports the model's outputs (IC50 and the bacteriostatic / bactericidal / eradication ",
        "thresholds) in those same ratio units (h). Supporting Information Table S10 tabulates the ",
        "index at 24, 48 and 72 h for simulated 20, 40, 50 and 60 mg/kg regimens. In this file the ",
        "AUC is the PIF (unbound, site-of-action) exposure predicted by the companion PBPK model, ",
        "not the plasma exposure. The ratio is carried directly rather than as an absolute AUC ",
        "divided by a model MIC parameter even though Sun 2026 DOES report the challenge isolate's ",
        "MIC (ZJWZ-A, 8 ug/mL, Supporting Information 'MIC and MBC determination'), because that ",
        "MIC does not reconstruct the paper's own index: back-solving MIC = AUC/index against the ",
        "companion PBPK returns 1.082 ug/mL from the PIF column of Table S10 and 0.743 ug/mL from ",
        "the plasma column, mutually inconsistent and roughly 8-10x away from 8. Adopting the ",
        "absolute-AUC form with mic = 8 would put the PIF index near 23 h at the recommended ",
        "40 mg/kg, below the 28.48 h bactericidal target, where Table S10 reports 172.68 h and the ",
        "paper concludes eradication. Set to 0 for drug-free control records so the sigmoid term ",
        "vanishes and the predicted 24 h change reduces to E0. See the vignette Assumptions and ",
        "deviations section."
      ),
      source_name        = "AUC/MIC (h) -- Sun 2026 eq 8 driver C; Table 3 rows 'AUC24/MIC for bacteriostatic effect', '... bactericidal effect', '... eradication effect'; Supporting Information Table S10 column 'Pulmonary interstitial fluid / AUC/MIC (h)'"
    )
  )

  population <- list(
    species             = "in vitro (hollow-fiber infection model simulating swine pulmonary interstitial fluid exposure)",
    n_subjects          = NA_integer_,
    n_studies           = 1L,
    organism            = "Pasteurella multocida isolate ZJWZ-A, isolated in Wuhan, Hubei Province, China in 2023 and deposited in the National Reference Laboratory for Veterinary Drug Residues at Huazhong Agricultural University; tilmicosin MIC 8 ug/mL by CLSI broth microdilution over a 0.25-64 ug/mL series",
    system              = paste0(
      "Hollow-fiber infection model (Cartridge C2011, FiberCell Systems Inc., Frederick, MD) with ",
      "a central chamber, dilution chamber and elimination chamber, each driven by its own ",
      "peristaltic pump; pump settings were derived from the Cmax, Tmax and t1/2 predicted by the ",
      "companion PBPK model in PIF at each simulated dose, entered into the public HFIM calculator ",
      "at https://pkpdia.shinyapps.io/hfs_app/ (Supporting Information Table S5). Growth medium ",
      "was TSB supplemented with 5% fetal bovine serum. Fifteen millilitres of ZJWZ-A suspension ",
      "at approximately 1 x 10^6 CFU/mL was inoculated into the cartridge and equilibrated for 1 h ",
      "at 37 C with 5% CO2. Bacterial counts were taken over 72 h; PK samples were drawn at 0, ",
      "0.25, 0.5, 0.75, 1, 2, 3, 4, 6, 8, 12, 24, 36, 48 and 72 h"
    ),
    disease_state       = "in vitro infection model",
    dose_range          = "simulated oral swine regimens of 20, 40, 50 and 60 mg/kg (Supporting Information Tables S5 and S10)",
    regions             = "China (Huazhong Agricultural University, Wuhan)",
    notes               = paste0(
      "Sun 2026 fitted the sigmoid Imax model separately to the plasma and pulmonary interstitial ",
      "fluid arms of the same HFIM experiment, so the two fits are packaged as ",
      "Sun_2026_tilmicosin_pkpd_plasma and Sun_2026_tilmicosin_pkpd_pif, sharing the vignette ",
      "Sun_2026_tilmicosin with the whole-body PBPK model Sun_2026_tilmicosin_pbpk. The PK/PD ",
      "modelling was performed in Phoenix 64 version 8.1.0. Sun 2026 also evaluated f%T > MIC and ",
      "reports that AUC0-24/MIC showed the strongest correlation with bactericidal efficacy in ",
      "both matrices, which is why only the AUC0-24/MIC fits are parameterised. In vitro ",
      "confirmation of the optimised 40 mg/kg once-daily 3-day regimen gave a 5.01 log10 CFU/mL ",
      "reduction at 72 h with no regrowth through 120 h (Supporting Information Figure S8)."
    )
  )

  ini({
    # =============================================================
    # Sun 2026 Table 3 -- sigmoid Imax PK/PD integration
    # column: PIF (pulmonary interstitial fluid)
    # =============================================================
    # Sun 2026 eq 8:
    #   E = E0 - (Imax * C^gamma) / (C^gamma + IC50^gamma)
    # E is the SIGNED log10(CFU/mL) change from the initial inoculum
    # over the 24 h window: it equals E0 (net growth of the drug-free
    # control) at zero exposure and falls to E0 - Imax at saturating
    # exposure. The minus sign is dropped by the symbol font in the
    # published PDF ("E = E0 (Imax x C^Gamma)/..."); the sign convention
    # is fixed independently by Table 3's three target rows, which this
    # parameterisation reproduces to better than 0.03%.
    le0 <- log(2.36)
    label("Change in bacterial count in the drug-free control over 24 h E0 (log10 CFU/mL)")  # Sun 2026 Table 3, PIF, E0 = 2.36 +/- 0.28

    limax <- log(6.82)
    label("Maximum antibacterial effect Imax, the span from control growth to maximum kill (log10 CFU/mL)")  # Sun 2026 Table 3, PIF, Imax = 6.82 +/- 0.58

    lic50 <- log(16.06)
    label("AUC0-24/MIC producing half the maximum effect IC50 (h)")  # Sun 2026 Table 3, PIF, IC50 = 16.06 +/- 0.16

    lhill <- log(2.27)
    label("Hill coefficient gamma, steepness of the AUC0-24/MIC effect curve (unitless)")  # Sun 2026 Table 3, PIF, gamma = 2.27 +/- 0.03

    # =============================================================
    # Starting bacterial density
    # =============================================================
    # FIXED experimental design input, not an estimated parameter.
    log10_cfu0 <- fixed(6)
    label("log10 starting bacterial density in the hollow-fiber cartridge (log10 CFU/mL)")  # Sun 2026 Supporting Information 'Establishment of hollow fiber infection model': "Fifteen milliliters of P. multocida ZJWZ-A bacterial solution (density: ~ 10^6 CFU/mL) was inoculated into the hollow fiber cartridge"

    # =============================================================
    # Residual error
    # =============================================================
    # Sun 2026 reports only the point estimates and their standard
    # errors in Table 3 and gives no residual standard deviation for the
    # PK/PD fit, so the density-scale residual SD is held at zero for
    # deterministic typical-value simulation. See the vignette
    # Assumptions and deviations section.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (0; not reported in Sun 2026)")  # Sun 2026 Table 3 reports standard errors of the estimates only, no residual SD
  })

  model({
    e0 <- exp(le0)
    imax <- exp(limax)
    ic50 <- exp(lic50)
    hill <- exp(lhill)

    # Sun 2026 eq 8 sigmoid Imax equation. `effect` is the SIGNED change
    # in log10(CFU/mL) accrued over a 24 h window (positive = net
    # growth, negative = net kill); it equals e0 at zero exposure and
    # approaches e0 - imax at saturating exposure.
    effect <- e0 - imax * AUCMIC_TILM^hill / (AUCMIC_TILM^hill + ic50^hill)

    # Sun 2026 fitted the change in log10 CFU/mL accrued over a 24 h
    # window against the AUC0-24/MIC of that window. Spreading that
    # change uniformly across the window makes log10(bact) move by
    # exactly `effect` over 24 h, so the trajectory matches the paper's
    # model at every counted time point:
    #   d(log10 N)/dt = effect / 24
    #   => d(N)/dt    = ln(10) * (effect / 24) * N
    d/dt(bact) <- log(10) * (effect / 24) * bact
    bact(0) <- 10^log10_cfu0

    # log10 CFU/mL observation with a 1-CFU/mL floor (matches the
    # Lee 2023 / Chen 2023 in-vitro PD convention so the log10 stays
    # finite if bact is driven below 1 CFU/mL).
    Cc <- log10(bact + 1)
    Cc ~ add(addSd)
  })
}
