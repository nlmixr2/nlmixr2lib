Tompkins_2026_romidepsin_rlu <- function() {
  description <- paste(
    "In vitro (N6 cells, a Jurkat-derived T-cell clone latently infected with",
    "an HIV-1 luciferase reporter virus). Sigmoid Emax concentration-response",
    "model for HIV-1 latency-reversing activity of romidepsin, measured as",
    "luciferase luminescence (relative light units, RLU) after a 24-hour",
    "static exposure (Luciferase Assay System readout, Figure 1b):",
    "luminescence = e0 + emax * C^hill / (ec50^hill + C^hill).",
    "The fit was restricted to 0.78-6.25 ng/mL so that it describes latency",
    "reversal before cytotoxicity sets in; above about 9 ng/mL the measured",
    "RLU falls again and this curve does not apply. The model is static: it",
    "maps one 24-hour constant concentration to one 24-hour endpoint and does",
    "not represent the exposure-time dependence the paper reports from its",
    "dose-fractionation experiments. No inter-individual variability and no",
    "residual error are reported. Companion files",
    "Tompkins_2026_romidepsin_rlu_oneglo.R (replicate RLU fit in the",
    "two-step assay) and Tompkins_2026_romidepsin_rfu.R (cytotoxicity fit)",
    "come from the same paper.",
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
    "effective concentration range' and in the Figure 1b panel.",
    sep = " "
  )

  vignette <- "Tompkins_2026_romidepsin"

  units <- list(
    time = "h",
    dosing = "(none; static in vitro concentration-response model driven by an external romidepsin concentration covariate)",
    concentration = "(observation luminescence is luciferase luminescence in relative light units per culture well; driving covariate CONC_ROMIDEPSIN_NGML is the romidepsin concentration in the culture medium in ng/mL)"
  )

  covariateData <- list(
    CONC_ROMIDEPSIN_NGML = list(
      description = "Romidepsin concentration applied to the culture medium, held constant over the 24-hour incubation.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Twofold serial dilutions from 0.78125 to 600 ng/mL were tested (Methods 'Exposure-response studies'); this Emax fit used only 0.78-6.25 ng/mL (Results; Figure 1b caption).",
        "The paper also expresses each concentration as AUC0-24h = 24 * C ng*h/mL (19-150 ng*h/mL for the fitted range); that is a relabelling of the same static concentration, not a second driver.",
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
      "About 100,000 viable N6 cells per well in 96-well plates, incubated 24 hours with romidepsin; each concentration in triplicate wells with one RLU reading per well (Methods 'Exposure-response studies').",
      "Two independent experiments are plotted in Figure 1a/1b (X and O symbols) and all individual well readings were pooled in the fit; 42 wells (2 experiments x 7 concentrations x 3 wells) inform this curve.",
      "Readout: Luciferase Assay System (Promega) after overnight lysis, read on a Veritas microplate luminometer (Supplemental Methods). The two-step ONE-Glo + Tox readout of Figure 1d uses a different luminometer and assay chemistry, so its RLU scale is not interchangeable with this one.",
      "Fitted in Phoenix PK/PD (Certara)."
    )
  )

  ini({
    # Tompkins 2026 Methods 'Statistics':
    #   RLU = E0 + Emax * C^gamma / (EC50^gamma + C^gamma)
    # Results: 'An Emax model (Figure 1b) was fit for RLU vs. RMD
    # concentrations constrained to 0.78-6.25 ng/mL ... The mean (%CV)
    # parameter estimates were EC50 = 3.2 ng/mL (1.4%), gamma = 8.3 ng/mL
    # (14.6%), E0 = 746 RLU (86.7%), and Emax = 49,035 RLU (2.7%) with an
    # EC90 ... calculated at 4.2 ng/mL.' The %CV values are the relative
    # standard errors of the estimates, not between-well variability: a
    # refit of the same equation to the 42 well readings in the paper's
    # Supporting Data (sheet 'Figure 1B') returns EC50 3.204 (SE 0.045),
    # gamma 8.27 (SE 1.20), E0 745 (SE 646), Emax 49,035 (SE 1,322).

    lec50 <- log(3.2)
    label("Log romidepsin concentration giving half-maximal RLU increase (ng/mL)") # Results; Figure 1b panel 'EC50=3.2ng/mL'

    lhill <- log(8.3)
    label("Log Hill coefficient of the RLU sigmoid (unitless)") # Results 'gamma = 8.3' (printed with ng/mL units in error; Figure 1b panel 'Hill slope=8.3')

    le0 <- log(746)
    label("Log baseline luminescence at zero romidepsin (RLU)") # Results 'E0 = 746 RLU (86.7%)'

    lemax <- log(49035)
    label("Log maximal RLU increase above baseline (RLU)") # Results 'Emax = 49,035 RLU (2.7%)'

    # Not reported: the fit is a single curve through pooled well readings,
    # and the paper gives no residual-error estimate.
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
