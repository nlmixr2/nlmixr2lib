Zhu_2026_denosumab_reference <- function() {
  description <- paste(
    "One-compartment population PK model with first-order subcutaneous absorption and parallel linear",
    "and Michaelis-Menten elimination from the central compartment for US-sourced reference denosumab",
    "(Prolia), fitted as the active comparator arm of a Phase III biosimilarity study in Chinese",
    "postmenopausal women with osteoporosis at high risk of fracture (Zhu 2026, CTR20201555). Body",
    "weight enters linear clearance as a power term normalised to 55.5 kg. Inter-individual variability",
    "was estimated on clearance and on the maximum Michaelis-Menten elimination rate; residual error is",
    "combined proportional and additive. Zhu 2026 fitted the reference product and the KN012 biosimilar",
    "as two separate models on the two treatment arms; the companion biosimilar model is",
    "Zhu_2026_denosumab_kn012."
  )
  reference <- paste(
    "Zhu X, Liu J, Mao Y, Li J, Shao F, Lin H.",
    "Comparison of KN012, a denosumab biosimilar, versus reference denosumab in Chinese postmenopausal",
    "women with osteoporosis: efficacy, safety, and population pharmacokinetics in a 12-month phase III",
    "study. Bone Rep. 2026;101916. doi:10.1016/j.bonr.2026.101916"
  )
  vignette <- "Zhu_2026_denosumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters linear clearance as a normalised power term, CL = CL_TV * (WT/55.5)^1.16, printed in",
        "the footnote to Zhu 2026 Supplementary Table 5. The centering constant is 55.5 kg, which is",
        "the mean body weight of the KN012 arm rather than this arm's own mean of 56.3 kg",
        "(Supplementary Table 3), so 55.5 kg is a shared model centering value used for both products",
        "and is reproduced here as printed. Body weight was the only covariate retained by the stepwise",
        "forward-addition / backward-elimination covariate search."
      ),
      source_name        = "WT"
    )
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Zhu 2026 Supplementary Fig. 1 (model
  # schema: Dose -> 'SC Site' -KA-> 'Central', with CL and Vmax/(Km+Conc)
  # leaving 'Central').
  compartmentData <- list(
    depot = list(
      analyte = "denosumab (US reference product)", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "denosumab (US reference product)", units = "mg",
      specimen = "serum", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 138L,
    n_observations = 925L,
    n_studies      = 1L,
    age_range      = "50.0-80.0 years",
    age_mean       = "65.3 years",
    weight_range   = "38.0-77.0 kg",
    weight_mean    = "56.3 kg",
    sex_female_pct = 100,
    race_ethnicity = c(`Han Chinese` = 98.6, Other = 1.4),
    disease_state  = "Postmenopausal osteoporosis at high risk of fracture",
    dose_range     = "60 mg subcutaneously at month 0 and month 6",
    regions        = "China (25 centres)",
    renal_function = "Creatinine clearance mean 75.1 mL/min (range 36.6-131); 23.2% normal (CrCL >= 90), 53.6% mild (60 <= CrCL < 90), 23.2% moderate (30 <= CrCL < 60) impairment",
    co_medication  = "Daily calcium carbonate D3 600 mg and vitamin D 800 IU from screening to end of study",
    notes          = paste(
      "US-denosumab arm of the population-PK analysis set (Zhu 2026 Results 3.5 and Supplementary",
      "Table 3). Twelve of the 138 participants were intensively sampled and the remainder sparsely",
      "sampled (Supplementary Table 1). Serum concentrations were measured by electrochemiluminescence",
      "on the Meso Scale Discovery platform. Estimation was by FOCE with interaction in NONMEM 7.4;",
      "the minimum success rate of 1000 bootstrap runs was 97.9%. Anti-drug antibodies were positive",
      "in 9 of 138 (6.5%) participants and ADA was screened but not retained as a covariate. Age,",
      "weight, BMI, BSA, baseline ALT, AST, albumin, total protein, total bilirubin, creatinine",
      "clearance and ADA status were all tested; only body weight on clearance was retained."
    )
  )

  ini({
    # ---- Structural PK: Zhu 2026 Supplementary Table 5, 'Final model' column.
    # Bioavailability was not estimated, so Vc and CL are apparent values
    # (Vc/F, CL/F) for the 60 mg subcutaneous route; the source table prints
    # them without the /F qualifier.
    lvc <- log(0.485)
    label("Apparent central volume of distribution (L)")  # Supplementary Table 5, row 'Vc, L' = 0.485 [RSE 7.0%, 95% CI 0.418-0.552]
    lcl <- log(0.116)
    label("Apparent linear clearance at 55.5 kg (L/day)")  # Supplementary Table 5, row 'CL, L/day' = 0.116 [RSE 3.4%, 95% CI 0.108-0.124]
    lka <- log(0.0149)
    label("First-order subcutaneous absorption rate constant (1/day)")  # Supplementary Table 5, row 'KA, 1/day' = 0.0149 [RSE 2.4%, 95% CI 0.0142-0.0156]

    # ---- Parallel Michaelis-Menten elimination from the central compartment.
    # Vmax is carried verbatim from the source table, but its printed unit is
    # wrong: it is a CONCENTRATION elimination rate, not an amount rate. See
    # the note in model() for the arithmetic that establishes this.
    lvmax <- log(0.304)
    label("Maximum Michaelis-Menten elimination rate Vmax (ug/mL/day)")  # Supplementary Table 5, row 'Vmax, mg/day' = 0.304 [RSE 10.6%, 95% CI 0.241-0.367]; printed unit corrected to ug/mL/day, see model()
    lkm <- log(0.058)
    label("Michaelis-Menten constant Km (ug/mL)")  # Supplementary Table 5, row 'Km, ug/mL' = 0.058 [RSE 5.9%, 95% CI 0.0513-0.0647]

    # ---- Covariate effect.
    e_wt_cl <- 1.16
    label("Power exponent on (WT/55.5) for linear clearance (unitless)")  # Supplementary Table 5, row 'Weight on CL' = 1.16 [RSE 14.1%, 95% CI 0.839-1.48]; footnote 'CL=CL_TV(WT/55.5)^1.16'

    # ---- Inter-individual variability.
    # Zhu 2026 Supplementary Table 5 reports IIV as a percentage under the
    # symbol omega, i.e. as a coefficient of variation for log-normally
    # distributed parameters. Variances below are omega^2 = log(1 + CV^2).
    etalcl ~ 0.0343715  # Supplementary Table 5, row 'omega(CL),%' = 18.7% CV [RSE 10.0%]; shrinkage 21.6%
    etalvmax ~ 0.2003587  # Supplementary Table 5, row 'omega(Vmax), %' = 47.1% CV [RSE 8.8%]; shrinkage 2.8%

    # ---- Residual error. Combined proportional and additive; the companion
    # KN012 model carries the proportional term only (Supplementary Table 4).
    propSd <- 0.202
    label("Proportional residual error (fraction)")  # Supplementary Table 5, row 'sigma (Prop), %' = 20.2% [RSE 2.9%, 95% CI 19.1-21.3]; shrinkage 11.4%
    addSd <- 0.014
    label("Additive residual error (ug/mL)")  # Supplementary Table 5, row 'sigma (Add), ug/mL' = 0.014 [RSE 8.8%, 95% CI 0.0116-0.0164]
  })

  model({
    # Individual parameters. Body weight enters linear clearance only, as the
    # normalised power term printed in the Supplementary Table 5 footnote.
    # Inter-individual variability was estimated on CL and Vmax only.
    cl <- exp(lcl + etalcl) * (WT / 55.5)^e_wt_cl
    vc <- exp(lvc)
    ka <- exp(lka)
    vmax <- exp(lvmax + etalvmax)
    km <- exp(lkm)

    # One-compartment disposition with first-order absorption from the
    # subcutaneous site and parallel linear plus Michaelis-Menten elimination
    # (Zhu 2026 Supplementary Fig. 1, which labels the saturable arm
    # 'Vmax/(Km+Conc)'). Amounts are in mg and volume in L, so Cc is in mg/L,
    # which equals ug/mL, and Km is a concentration in ug/mL.
    #
    # UNIT CORRECTION on Vmax. Supplementary Table 5 prints 'Vmax, mg/day',
    # which would make the saturable term vmax * Cc / (km + Cc) an amount rate
    # directly. Taken that way the model does NOT reproduce the paper's own
    # simulated exposures for a typical 60 kg subject given 60 mg
    # (Supplementary Table 6): it returns Cmax 3.71 vs 4.84 (-23%), AUC0-6mon
    # 153.0 vs 266.00 (-43%) and AUCinf 153.8 vs 267.49 (-43%). Treating Vmax
    # as a CONCENTRATION rate in ug/mL/day, so that the amount rate is
    # vmax * vc * Cc / (km + Cc), reproduces all three to within 0.16%
    # (4.84 / 266.1 / 267.9). The same correction reproduces all three
    # published values for the companion KN012 model to within 0.22% using
    # that model's own different Vc of 0.448 L, so the factor is Vc and not a
    # single fitted fudge constant. A change in Km cannot explain the gap: no
    # single Km value reproduces Cmax and AUCinf jointly. The printed 'mg/day'
    # is therefore a unit-tag slip and Vmax is the saturable elimination rate
    # on the concentration scale. See the vignette Errata.
    Cc <- central / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl * Cc - vmax * vc * Cc / (km + Cc)

    Cc ~ add(addSd) + prop(propSd)
  })
}
