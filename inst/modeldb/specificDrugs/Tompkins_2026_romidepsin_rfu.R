Tompkins_2026_romidepsin_rfu <- function() {
  description <- paste(
    "In vitro (N6 cells, a Jurkat-derived T-cell clone latently infected with",
    "an HIV-1 luciferase reporter virus). Sigmoid inhibitory (Imax)",
    "concentration-response model for romidepsin cytotoxicity, measured as",
    "CellTiter-Fluor viable-cell fluorescence (relative fluorescence units,",
    "RFU) after a 24-hour static exposure in the two-step ONE-Glo + Tox",
    "assay (Figure 1f): fluorescence = e0 - imax * C^hill / (ic50^hill + C^hill),",
    "where ic50 is the paper's CC50 (5.1 ng/mL) and the CC90 is 8.4 ng/mL.",
    "Fitted to 0.78-200 ng/mL; at 300-600 ng/mL RFU falls below the fitted",
    "plateau and the curve does not apply. Static model: no exposure-time",
    "dependence. No inter-individual variability and no residual error are",
    "reported. The luminescence readout of the same wells is",
    "Tompkins_2026_romidepsin_rlu_oneglo.R; the paper's therapeutic window",
    "for latency reversal is the gap between its EC50 and this CC90.",
    sep = " "
  )

  reference <- paste(
    "Tompkins LAR, Poliseno A, Ward A, Browne EP, Gilbert C, Sykes C,",
    "Kashuba ADM, Cottrell ML.",
    "Empirical PK/PD Model for Differentiating Intrinsic Dosing Rules by",
    "Therapeutic Indication With Relevance for HIV-1 Latency Reversal.",
    "Clin Pharmacol Ther. 2026;120(5):1414-1426. doi:10.1002/cpt.70488.",
    "The model equation is in Methods 'Statistics'; the estimates are in",
    "Results 'RMD-mediated latency-reversing activity displays a narrow",
    "effective concentration range' and in the Figure 1f panel.",
    sep = " "
  )

  vignette <- "Tompkins_2026_romidepsin"

  units <- list(
    time = "h",
    dosing = "(none; static in vitro concentration-response model driven by an external romidepsin concentration covariate)",
    concentration = "(observation fluorescence is CellTiter-Fluor viable-cell fluorescence in relative fluorescence units per culture well; driving covariate CONC_ROMIDEPSIN_NGML is the romidepsin concentration in the culture medium in ng/mL)"
  )

  covariateData <- list(
    CONC_ROMIDEPSIN_NGML = list(
      description = "Romidepsin concentration applied to the culture medium, held constant over the 24-hour incubation.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Same twofold dilution scheme as Figure 1a (0.78125-600 ng/mL); this Imax fit used 0.78-200 ng/mL (Results; Figure 1f caption and panel).",
        "Concentrations are nominal medium concentrations in RPMI with 10% FBS; the Discussion estimates an unbound fraction of up to about 15% in this medium, so the potency is NOT an unbound or a human-plasma potency."
      ),
      source_name = "RMD concentration (C)"
    )
  )

  population <- list(
    species = "in vitro (N6 cells: Jurkat E6-1-derived T-cell clone carrying a single integrated HIV-1 NL4-3 provirus expressing firefly luciferase)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = NA_character_,
    weight_range = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state = "Not applicable -- clonal cell-line model of HIV-1 latency.",
    dose_range = "Romidepsin 0.78-200 ng/mL for 24 hours (fitted range); 0.78-600 ng/mL tested.",
    regions = NA_character_,
    notes = paste(
      "One experiment (Figure 1c-f), about 100,000 viable N6 cells per well, each concentration in triplicate wells with one reading per well. The Results state the fit used 0.78-200 ng/mL (51 wells, 17 concentrations x 3); the Supporting Data sheet 'Figure 1F' lists the 48 wells up to 150 ng/mL. Refitting the Methods equation to either set reproduces the printed estimates. Twenty-four untreated control wells read a mean of 1.42e7 RFU (SD 5.50e5) at 24 hours, about 78 RFU per Trypan-excluding live cell (Results).",
      "RFU is proportional to the number of viable, membrane-intact cells (Figure S2A). The loss of RFU at high concentration reflects an increase in dead cells rather than a fall in live-cell number: after 24 hours at 600 ng/mL viability fell 24% by cell count while live cells were 0.95-fold of baseline (Results).",
      "Fitted in Phoenix PK/PD (Certara)."
    )
  )

  ini({
    # Tompkins 2026 Methods 'Statistics':
    #   RFU = E0 - Imax * C^gamma / (CC50^gamma + C^gamma)
    # Results: 'The CC50 ... was estimated at 5.1 ng/mL (CV%: 2.4%), with
    # mean (%CV) parameter estimates of gamma = -4.4 ng/mL (9.9%),
    # E0 = 1.4e7 RFU (0.8%), and Imax = 6e6 RFU (2.0%), and a CC90 ...
    # calculated at 8.4 ng/mL.' The %CV values are relative standard errors.

    lic50 <- log(5.1)
    label("Log romidepsin concentration giving half-maximal RFU loss, CC50 (ng/mL)") # Results; Figure 1f panel 'CC50=5.1ng/mL'

    # Printed as gamma = -4.4 (Results) and 'Hill slope=-4.4' (Figure 1f
    # panel). The minus sign is the plotting software's sign convention for
    # a falling curve, not part of the printed equation: substituting -4.4
    # into the Methods equation with a positive Imax makes RFU RISE from
    # E0 - Imax to E0, the opposite of the data. The magnitude is encoded.
    # Refitting the Methods equation to the 51 well readings at
    # 0.78-200 ng/mL in the Supporting Data returns gamma = +4.37
    # (RSE 10.0%), CC50 5.08 (RSE 2.4%), E0 1.364e7 (RSE 0.7%) and
    # Imax 6.04e6 (RSE 2.0%), and the printed CC90 = CC50 * 9^(1/4.4)
    # = 8.4 ng/mL holds only for the positive exponent.
    lhill <- log(4.4)
    label("Log Hill coefficient of the RFU inhibitory sigmoid (unitless)")

    le0 <- log(1.4e7)
    label("Log viable-cell fluorescence at zero romidepsin (RFU)") # Results 'E0 = 1.4e7 RFU (0.8%)'

    limax <- log(6e6)
    label("Log maximal RFU loss below baseline (RFU)") # Results 'Imax = 6e6 RFU (2.0%)'

    # Not reported: single curve through the well readings, no
    # residual-error estimate in the source.
    addSd <- fixed(0)
    label("Additive residual error on RFU (zero; not reported in the source)")
  })

  model({
    ic50 <- exp(lic50)
    hill <- exp(lhill)
    e0 <- exp(le0)
    imax <- exp(limax)

    fluorescence <- e0 -
      imax *
        CONC_ROMIDEPSIN_NGML^hill /
        (ic50^hill + CONC_ROMIDEPSIN_NGML^hill)

    fluorescence ~ add(addSd)
  })
}
