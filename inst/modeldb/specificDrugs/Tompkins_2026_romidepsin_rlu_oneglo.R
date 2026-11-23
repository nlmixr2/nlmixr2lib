Tompkins_2026_romidepsin_rlu_oneglo <- function() {
  description <- paste(
    "In vitro (N6 cells, a Jurkat-derived T-cell clone latently infected with",
    "an HIV-1 luciferase reporter virus). Sigmoid Emax concentration-response",
    "model for HIV-1 latency-reversing activity of romidepsin, measured as",
    "luciferase luminescence (relative light units, RLU) after a 24-hour",
    "static exposure in the two-step ONE-Glo + Tox assay (Figure 1d), which",
    "reads viability (RFU) and RLU from the same wells:",
    "luminescence = e0 + emax * C^hill / (ec50^hill + C^hill).",
    "This is the paper's replicate of the Figure 1b fit",
    "(Tompkins_2026_romidepsin_rlu.R) in a different assay; the potency agrees",
    "(EC50 3.7 vs 3.2 ng/mL) but the RLU scale does not, because the assay",
    "chemistry and luminometer differ. Fitted to 0.78-6.25 ng/mL only; above",
    "about 9 ng/mL RLU falls with cytotoxicity and this curve does not apply.",
    "Static model: no exposure-time dependence. No inter-individual",
    "variability and no residual error are reported. The viability readout",
    "of the same wells is Tompkins_2026_romidepsin_rfu.R.",
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
    "effective concentration range' and in the Figure 1d panel.",
    sep = " "
  )

  vignette <- "Tompkins_2026_romidepsin"

  units <- list(
    time = "h",
    dosing = "(none; static in vitro concentration-response model driven by an external romidepsin concentration covariate)",
    concentration = "(observation luminescence is ONE-Glo luciferase luminescence in relative light units per culture well; driving covariate CONC_ROMIDEPSIN_NGML is the romidepsin concentration in the culture medium in ng/mL)"
  )

  covariateData <- list(
    CONC_ROMIDEPSIN_NGML = list(
      description = "Romidepsin concentration applied to the culture medium, held constant over the 24-hour incubation.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Same twofold dilution scheme as Figure 1a (0.78125-600 ng/mL); this Emax fit used only 0.78-6.25 ng/mL (Results; Figure 1d caption).",
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
    dose_range = "Romidepsin 0.78-6.25 ng/mL for 24 hours (fitted range); 0.78-600 ng/mL tested.",
    regions = NA_character_,
    notes = paste(
      "One experiment (Figure 1c-f), about 100,000 viable N6 cells per well, each concentration in triplicate wells with one reading per well; 21 wells (7 concentrations x 3) inform this curve. Twenty-four untreated control wells read a mean of 85 RLU (SD 17) at 24 hours (Results).",
      "Readout: two-step ONE-Glo + Tox assay (Promega) read on a CLARIOstar Plus plate reader -- CellTiter-Fluor viability fluorescence first, then ONE-Glo luminescence in the same well (Supplemental Methods).",
      "Fitted in Phoenix PK/PD (Certara)."
    )
  )

  ini({
    # Tompkins 2026 Methods 'Statistics':
    #   RLU = E0 + Emax * C^gamma / (EC50^gamma + C^gamma)
    # Results: 'Like our previous estimates, the mean (%CV) parameter
    # estimates were: EC50 = 3.7 ng/mL (0.9%), gamma = 7.6 ng/mL (3.9%),
    # E0 = 177 RLU (20.9%), and Emax = 7,626 RLU (1.3%), with an EC90
    # calculated at 5.0 ng/mL.' The %CV values are relative standard errors
    # of the estimates. A refit of the same equation to the 21 well readings
    # in the Supporting Data (sheet 'Figure 1D') returns EC50 3.737
    # (SE 0.033), gamma 7.60 (SE 0.30), E0 176 (SE 37), Emax 7,627 (SE 102).

    lec50 <- log(3.7)
    label("Log romidepsin concentration giving half-maximal RLU increase (ng/mL)") # Results; Figure 1d panel 'EC50=3.7ng/mL'

    lhill <- log(7.6)
    label("Log Hill coefficient of the RLU sigmoid (unitless)") # Results 'gamma = 7.6' (printed with ng/mL units in error; Figure 1d panel 'Hill slope=7.6')

    le0 <- log(177)
    label("Log baseline luminescence at zero romidepsin (RLU)") # Results 'E0 = 177 RLU (20.9%)'

    lemax <- log(7626)
    label("Log maximal RLU increase above baseline (RLU)") # Results 'Emax = 7,626 RLU (1.3%)'

    # Not reported: single curve through the well readings, no
    # residual-error estimate in the source.
    addSd <- fixed(0)
    label("Additive residual error on RLU (zero; not reported in the source)")
  })

  model({
    ec50 <- exp(lec50)
    hill <- exp(lhill)
    e0 <- exp(le0)
    emax <- exp(lemax)

    luminescence <- e0 +
      emax *
        CONC_ROMIDEPSIN_NGML^hill /
        (ec50^hill + CONC_ROMIDEPSIN_NGML^hill)

    luminescence ~ add(addSd)
  })
}
