McKenna_2017_doxorubicin_mdamb453_pk <- function() {
  description <- "In vitro (MDA-MB-453 human triple-negative breast cancer cell line). Three-compartment cellular pharmacokinetic model of doxorubicin uptake and DNA binding (McKenna 2017 Eqs 1-3): drug in the well medium (conc_dox, the extracellular compartment) exchanges with a free intracellular pool (intracellular_free) at rate constants kEF and kFE scaled by the medium-to-cell volume ratio, and free drug binds irreversibly to DNA (intracellular_bound) at kFB. All states are concentrations (nM). A treatment is a bolus of the applied concentration into conc_dox at t = 0, and media replacement is an evid = 5 replacement of conc_dox by 0 at the end of exposure. The maximum bound concentration CB,max (the bound concentration after the free pool has washed out) is the summary statistic that drives the companion treatment-response model McKenna_2017_doxorubicin_mdamb453_pd. Deterministic: rate constants are per-cell-line point estimates (Fig 2f-h, digitised); no between-well variability or residual error was modelled. The intracellular volume vI is not printed in the source and was back-solved from the CB,max values the paper prints in Figs 3-4. Sibling cell lines: McKenna_2017_doxorubicin_sum149pt_pk, McKenna_2017_doxorubicin_mdamb231_pk, McKenna_2017_doxorubicin_mdamb468_pk."
  reference <- "McKenna MT, Weis JA, Barnes SL, Tyson DR, Miga MI, Quaranta V, Yankeelov TE. A Predictive Mathematical Modeling Approach for the Study of Doxorubicin Treatment in Triple Negative Breast Cancer. Sci Rep. 2017;7:5725. doi:10.1038/s41598-017-05902-z"
  vignette <- "McKenna_2017_doxorubicin"

  units <- list(
    time = "h",
    dosing = "nM (doxorubicin concentration applied to the well medium, as a bolus into conc_dox)",
    concentration = "nM"
  )

  # The three states are the paper's own compartments (Fig 2a): the well
  # medium, free intracellular drug and DNA-bound intracellular drug. The
  # medium state follows the conc_dox precedent of
  # Mody_2023_doxorubicin_dexrazoxane_mdamb468.R; none of the three is a
  # registered canonical compartment.
  paper_specific_compartments <- c("conc_dox", "intracellular_free", "intracellular_bound")
  dosing <- c("conc_dox")

  compartmentData <- list(
    conc_dox = list(analyte = "doxorubicin", units = "nM", specimen = "administration site", verified = TRUE),
    intracellular_free = list(analyte = "doxorubicin (free)", units = "nM", specimen = "tumor", verified = TRUE),
    intracellular_bound = list(analyte = "doxorubicin (DNA-bound)", units = "nM", specimen = "tumor", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "in vitro (MDA-MB-453 human triple-negative breast cancer cell line, parental)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    disease_state = "triple-negative breast cancer (MDA-MB-453, Lehmann luminal androgen receptor subtype). About 10,000 cells per well of a 96-well plate, imaged every ~15 min by brightfield and fluorescence microscopy from 1 h before to ~24 h after drug application; doxorubicin's intrinsic fluorescence gives the intracellular and extracellular concentration time courses.",
    dose_range = "312-2500 nM (two-fold series) applied for 6 or 12 h then removed by media replacement; each condition in duplicate.",
    cell_line = "MDA-MB-453 (ATCC).",
    culture = "96-well microtiter plates, 250 uL medium per well, maintained per ATCC recommendations.",
    notes = "The three rate constants were fitted jointly to all conditions of the cell line by nonlinear least squares with an AUC-normalised objective (Supplement); the mean percent error of the fit was 23.5%. Values are digitised from the bar charts of Fig 2f-h (the paper prints no table)."
  )

  ini({
    lkef  <- log(3.1e-07) ; label("Rate constant from the medium into the free intracellular compartment kEF (1/h)")  # Fig 2f bar, MDA-MB-453 (digitised; axis x10^-5 1/h)
    lkfe  <- log(0.092) ; label("Rate constant from the free intracellular compartment back to the medium kFE (1/h)")  # Fig 2g bar, MDA-MB-453 (digitised)
    lkfb  <- log(0.0807) ; label("Binding rate constant from free to DNA-bound intracellular doxorubicin kFB (1/h)")  # Fig 2h bar, MDA-MB-453 (digitised; axis x10^-1 1/h)
    lvext <- fixed(log(250)) ; label("Extracellular (medium) volume vE (uL)")  # Methods, Doxorubicin compartment modeling: vE set to 250 uL, the media volume per well
    # vI = 10,000 seeded cells x an ellipsoid-fit cell volume (Methods); the cell volume is not printed. Back-solved by the maintainers so that the bound concentration after washout per unit AUC (kFB/(kFE+kFB) * kEF * vE/vI = 0.01935 nM per nM*h) matches the geometric mean of the CB,max/AUC ratios printed in Figs 3-4: 156.2 nM x 6 h -> 18.1 (Fig4); 312.5 nM x 6 h -> 36.3 (Fig4); 625 nM x 6 h -> 72.6 (Fig4).
    lvint <- fixed(log(0.001872)) ; label("Intracellular volume vI (uL)")  # back-solved from the Figs 3-4 printed CB,max values; see comment above
  })

  model({
    kef  <- exp(lkef)
    kfe  <- exp(lkfe)
    kfb  <- exp(lkfb)
    vext <- exp(lvext)
    vint <- exp(lvint)

    # Eqs 1-3. The free and bound pools share the intracellular volume vI.
    d/dt(conc_dox)            <- kfe * vint / vext * intracellular_free - kef * conc_dox
    d/dt(intracellular_free)  <- kef * vext / vint * conc_dox - kfe * intracellular_free - kfb * intracellular_free
    d/dt(intracellular_bound) <- kfb * intracellular_free

    # Total intracellular concentration, the quantity imaged in Fig 2e
    cint <- intracellular_free + intracellular_bound
  })
}
