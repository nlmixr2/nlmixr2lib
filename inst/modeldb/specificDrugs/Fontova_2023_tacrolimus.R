Fontova_2023_tacrolimus <- function() {
  description <- paste(
    "Direct inhibitory Imax model (no baseline, Hill coefficient 1) relating",
    "the intracellular tacrolimus concentration in peripheral blood",
    "mononuclear cells (PBMCs) to calcineurin (CN) activity in 23 stable",
    "adult kidney transplant recipients studied at steady state on",
    "twice-daily immediate-release tacrolimus (TAC-IR) and, four weeks after",
    "conversion, on once-daily extended-release tacrolimus (TAC-LCP); one",
    "pooled model for both formulations. The response is the percentage by",
    "which a patient's CN activity exceeds that patient's own nadir CN",
    "activity (the paper's '% inhibition Inadir'), and it falls from E0 at",
    "zero concentration towards the nadir as the concentration rises:",
    "E = E0 * (1 - C / (IC50 + C)). Naive-pooled NONMEM fit with no",
    "inter-individual variability. PD-only model: the measured intracellular",
    "concentration is supplied per observation row as the covariate",
    "CONC_TACROLIMUS_PGMC (pg/million cells). The paper analysed",
    "pharmacokinetics non-compartmentally and fitted no PK model.",
    sep = " "
  )

  reference <- paste(
    "Fontova P, Colom H, Rigo-Bonnin R, Bestard O, Vidal-Alabro A,",
    "van Merendonk LN, Cerezo G, Polo C, Montero N, Melilli E, Coloma A,",
    "Manonelles A, Cruzado JM, Grinyo JM, Lloberas N. The Effect of",
    "Intracellular Tacrolimus Exposure on Calcineurin Inhibition in",
    "Immediate- and Extended-Release Tacrolimus Formulations.",
    "Pharmaceutics 2023;15(5):1481. doi:10.3390/pharmaceutics15051481.",
    sep = " "
  )

  vignette <- "Fontova_2023_tacrolimus"

  units <- list(
    time = "h",
    dosing = "(none; PD-only model driven by the measured intracellular tacrolimus concentration covariate)",
    concentration = "pg/million cells (driving covariate CONC_TACROLIMUS_PGMC); the observation calcineurin_pct_above_nadir is in percent"
  )

  covariateData <- list(
    CONC_TACROLIMUS_PGMC = list(
      description = "Intracellular tacrolimus concentration in peripheral blood mononuclear cells (PBMCs) at the time of each calcineurin-activity observation, supplied per row.",
      units = "pg/million cells",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying per observation row; the C of Fontova 2023 Equation 4.",
        "PBMCs were isolated by Ficoll density gradient, lysed in hypotonic buffer, and tacrolimus was quantified by UHPLC-MS/MS (Methods 2.2).",
        "Total intracellular (not unbound) concentration.",
        "Steady-state geometric means (Table 2): TAC-IR C0 26.78, Cmax 93.06; TAC-LCP C0 25.63, Cmax 70.65 pg/million cells; the concentrations used in the fit spanned roughly 5 to 420 pg/million cells (Figure 2D).",
        "Set to 0 for a drug-free state (E then equals E0)."
      ),
      source_name = "C"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 23L,
    n_studies = 1L,
    age_range = "median 58.80 years (IQR 48.44-71.12), enrolled cohort of 25 (Table 1)",
    weight_range = "median 73.00 kg (IQR 62.75-79.90), enrolled cohort of 25 (Table 1)",
    sex_female_pct = 28,
    race_ethnicity = NA_character_,
    disease_state = "Stable adult kidney transplant recipients at least six months after transplantation (median 1.84 years, IQR 0.97-3.88) on tacrolimus with whole-blood C0 between 5 and 10 ng/mL at steady state; one of 25 had a prior kidney transplant.",
    dose_range = paste(
      "TAC-IR (Prograf or Adoport) twice daily, total daily dose median 3.0 mg (IQR 2.25-5.00);",
      "then TAC-LCP (Envarsus) once daily after a 1:0.7 dose conversion, total daily dose median 2.0 mg (IQR 1.63-3.50);",
      "doses taken fasting, at least 1 h before meals (Methods 2.1, Table 1)."
    ),
    regions = "Spain (Bellvitge University Hospital, single centre)",
    cyp3a5_distribution = "*1/*3 n = 5 (20%), *3/*3 n = 20 (80%) (Table 1)",
    abcb1_distribution = "T/T n = 6 (24%), C-carriers n = 19 (76%) (Table 1)",
    sampling_window = "Paired PBMC tacrolimus and PBMC calcineurin activity sampled pre-dose and at 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 12.5, 13, 14, 15, 20 and 24 h post-dose, once on TAC-IR and once four weeks after conversion to TAC-LCP (Methods 2.1).",
    notes = paste(
      "Post-hoc analysis of a prospective, open-label, single-centre crossover trial (NCT02961608, 2016-2018).",
      "25 enrolled; 23 completed both intracellular PK/PD profiles (Results 3.1).",
      "Calcineurin activity was measured in PBMC lysate as dephosphorylation of the RIIp peptide (pmol RII/min/mg protein; Methods 2.3) and converted per patient to the percentage above that patient's nadir activity before fitting (Equation 2).",
      "Demographics in Table 1 are for the 25 enrolled patients."
    )
  )

  ini({
    # Fontova 2023 Methods 2.4.2 Equation 4 with gamma = 1, selected in
    # Results 3.4 as the best model ('simple inhibitory Imax model without
    # baseline'). Naive-pooled NONMEM 7.4 FOCE-I fit; the paper reports
    # the two estimates in the text only, without their RSEs.
    le0 <- log(30.4)
    label("Calcineurin activity above the patient's nadir at zero intracellular concentration, E0 (percent)")
    # Results 3.4 and Discussion: 'The maximum reduction of % inhibition
    # Inadir (E0 ...) was 30.4 %'. Figure 2E's predicted values (about 3
    # to 26 over the observed concentrations) agree with this value.

    lic50 <- log(43.9)
    label("Intracellular tacrolimus concentration giving half of E0, IC50 (pg/million cells)")
    # Abstract, Results 3.4 and Discussion: IC50 = 43.9 pg/million cells.

    # The residual-error model and its magnitude are not reported (Methods
    # 2.4.2 names only the FOCE-I method); held at zero.
    addSd <- fixed(0)
    label("Additive residual SD on the percentage above nadir (percentage points); not reported")
  })

  model({
    e0 <- exp(le0)
    ic50 <- exp(lic50)

    # Fontova 2023 Equation 4, gamma = 1: E = E0 * (1 - C / (IC50 + C))
    calcineurin_pct_above_nadir <- e0 *
      (1 - CONC_TACROLIMUS_PGMC / (ic50 + CONC_TACROLIMUS_PGMC))
    calcineurin_pct_above_nadir ~ add(addSd)
  })
}
