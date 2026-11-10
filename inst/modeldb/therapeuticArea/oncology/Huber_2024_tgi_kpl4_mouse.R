Huber_2024_tgi_kpl4_mouse <- function() {
  description <- paste(
    "Preclinical (mouse KPL-4 breast-cancer cell-line xenograft).",
    "Mayneord-like semi-mechanistic tumour growth inhibition model: the",
    "tumour radius grows linearly at rate g and is reduced at a maximal rate",
    "d scaled by the in-vitro Hill function of the unbound plasma",
    "concentration, dR/dt = g - d * C^hill / (C^hill + IC50^hill).",
    "Xenograft-specific g and d from Huber and Mistry's re-analysis of the",
    "Diegmiller 2022 PI3K-pathway inhibitor xenograft studies; IC50 defaults",
    "to taselisib (GDC-0032) in KPL-4. PD-only: the unbound concentration",
    "is supplied as the time-varying covariate CEFFECT.",
    sep = " "
  )
  reference <- paste(
    "Huber HJ, Mistry HB. Explaining in-vitro to in-vivo efficacy",
    "correlations in oncology pre-clinical development via a",
    "semi-mechanistic mathematical model.",
    "J Pharmacokinet Pharmacodyn. 2024;51(2):169-185.",
    "doi:10.1007/s10928-023-09891-7.",
    "Xenograft data re-analysed from Diegmiller R et al. Growth-rate model",
    "predicts in vivo tumor response from in vitro data.",
    "CPT Pharmacometrics Syst Pharmacol. 2022;11(9):1183-1193.",
    "doi:10.1002/psp4.12836.",
    sep = " "
  )
  vignette <- "Huber_2024_ivivc_tumor_stasis"
  units <- list(
    time = "day",
    dosing = "n/a (no dose events; the unbound plasma concentration enters as the covariate CEFFECT)",
    concentration = "umol/L"
  )

  covariateData <- list(
    CEFFECT = list(
      description = "Time-varying unbound (free) plasma concentration of the anti-proliferative compound",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Huber 2024 assumes the free-plasma-drug hypothesis with no delay",
        "between the unbound plasma concentration and the effect on the",
        "tumour cells (Methods, 'Semi-mechanistic model'), so the unbound",
        "plasma concentration is the effect-site driver. Must be in the same",
        "unit as IC50 (only the ratio CEFFECT / IC50 enters the model). The",
        "Appendix A re-analysis held it constant at the published average",
        "unbound concentration of each dose group (column 'Cavg' of the",
        "supplementary CSV); the main-text analysis uses a repeated",
        "mono-exponential steady-state profile from cmax to ctrough. Set to 0",
        "for vehicle controls. Supply a time-varying profile on a dense grid:",
        "rxode2 does not restart the integrator where a covariate changes."
      ),
      source_name = "c_plasma,free(t) (Eq. 7, 19, 30); 'Cavg' in the Appendix A supplementary CSV"
    )
  )

  compartmentData <- list(
    tumor_size = list(
      analyte = "tumour radius",
      units = "mm",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (KPL-4 cell-line xenograft)",
    n_subjects = NA_integer_,
    n_studies = 9L,
    disease_state = paste(
      "KPL-4 human breast-cancer cell-line",
      "xenografts treated with the PI3K-pathway inhibitors alpelisib, GDC-0032",
      "(taselisib), GDC-0077 (inavolisib) and",
      "GDC-0941 (pictilisib)"
    ),
    dose_range = "0.2-150 mg/kg (compound-dependent; 29 dose groups in the Appendix A supplementary CSV)",
    notes = paste(
      "Not an animal-level fit. Huber and Mistry (Appendix A) took the",
      "published summary data of Diegmiller 2022 (Genentech) and converted",
      "each dose group's initial and final tumour volumes to a spherical",
      "radius, giving an observed net radial growth rate per group. g was",
      "taken from the vehicle controls and d was optimised so that",
      "g - d * Cavg / (Cavg + IC50) reproduced the observed net rates across",
      "all compounds and doses. n_studies counts the distinct Diegmiller",
      "study IDs in the CSV (1, 2, 6, 11, 12, 33, 36, 43, 52); the",
      "number of animals is not reported. No between-animal variability and",
      "no residual-error model were reported."
    )
  )

  ini({
    # ---- Xenograft-specific in-vivo parameters (Appendix A) ---------------
    lp <- log(0.20)
    label("Linear tumour-radius growth rate g (mm/day)") # Appendix A text: '0.20 and 0.30 mm/day for KPL-4'; CSV column 'Growth.Rate..mm.day.' = 0.2
    lemax <- log(0.30)
    label("Maximal drug-induced tumour-radius decay rate d (mm/day)") # Appendix A text (KPL-4 d = 0.30 mm/day); CSV column 'd..mm.day.' = 0.3

    # ---- Compound-specific in-vitro pharmacology --------------------------
    # Default compound: GDC-0032 (taselisib), the compound tested in both
    # xenografts. Other KPL-4 values from the same CSV column: alpelisib
    # 0.921081556, GDC-0077 0.037580179, GDC-0941 0.216294289 umol/L.
    # The CSV does not label the unit; umol/L is the unit Diegmiller 2022
    # reports for these in-vitro potencies.
    lic50 <- log(0.062188168)
    label("In-vitro IC50 (GI50) of the compound in the cell line (umol/L)") # Appendix A supplementary CSV, column 'GI50', rows 'GDC-0032,KPL-4'
    hill <- fixed(1)
    label("Hill coefficient of the in-vitro dose-response (unitless)") # Eq. 19-20 (hill = 1 in the Appendix A re-analysis); general form Eq. 7 and 30

    # ---- Initial condition ----------------------------------------------
    lrbase <- fixed(log(3))
    label("Initial tumour radius R0 (mm)") # Supplement 1 (PDF) section 1: 'R0 ~ 3 cm' for the authors' studies, read as 3 mm (see vignette)
  })

  model({
    p <- exp(lp)
    emax <- exp(lemax)
    ic50 <- exp(lic50)
    rbase <- exp(lrbase)

    # In-vitro Hill function of the unbound plasma concentration (Eq. 7)
    drugEffect <- CEFFECT^hill / (CEFFECT^hill + ic50^hill)

    # Mayneord-like linear radial growth with drug-induced decay (Eq. 8, 19, 30)
    tumor_size(0) <- rbase
    d/dt(tumor_size) <- p - emax * drugEffect

    # Spherical tumour volume (Eq. 1; used for TGI in Eq. 9)
    tumor_vol <- 4 / 3 * pi * tumor_size^3
  })
}
