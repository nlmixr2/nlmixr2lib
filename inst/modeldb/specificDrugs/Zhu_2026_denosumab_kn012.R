Zhu_2026_denosumab_kn012 <- function() {
  description <- paste(
    "One-compartment population PK model with first-order subcutaneous absorption and parallel linear",
    "and Michaelis-Menten elimination from the central compartment for KN012, a candidate denosumab",
    "biosimilar, in Chinese postmenopausal women with osteoporosis at high risk of fracture (Zhu 2026,",
    "Phase III study CTR20201555). Body weight enters linear clearance as a power term normalised to",
    "55.5 kg. Inter-individual variability was estimated on clearance and on the maximum",
    "Michaelis-Menten elimination rate; residual error is proportional only. Zhu 2026 fitted KN012 and",
    "the US reference denosumab as two separate models on the two treatment arms; the companion",
    "reference-product model is Zhu_2026_denosumab_reference."
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
        "Enters linear clearance as a normalised power term, CL = CL_TV * (WT/55.5)^1.34, printed in",
        "the footnote to Zhu 2026 Supplementary Table 4. The centering constant 55.5 kg is the mean",
        "body weight of the KN012 population-PK analysis set (Supplementary Table 3); the footnote to",
        "the US reference-product model (Supplementary Table 5) uses the same 55.5 kg constant even",
        "though that arm's own mean weight was 56.3 kg, so 55.5 kg is a shared model centering value",
        "and not a per-arm re-derived mean. Body weight was the only covariate retained by the stepwise",
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
      analyte = "KN012 (denosumab biosimilar)", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "KN012 (denosumab biosimilar)", units = "mg",
      specimen = "serum", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 139L,
    n_observations = 922L,
    n_studies      = 1L,
    age_range      = "51.0-78.0 years",
    age_mean       = "64.8 years",
    weight_range   = "35.0-76.0 kg",
    weight_mean    = "55.5 kg",
    sex_female_pct = 100,
    race_ethnicity = c(`Han Chinese` = 98.6, Other = 1.4),
    disease_state  = "Postmenopausal osteoporosis at high risk of fracture",
    dose_range     = "60 mg subcutaneously at month 0 and month 6",
    regions        = "China (25 centres)",
    renal_function = "Creatinine clearance mean 73.3 mL/min (range 33.1-134); 20.1% normal (CrCL >= 90), 54.0% mild (60 <= CrCL < 90), 25.9% moderate (30 <= CrCL < 60) impairment",
    co_medication  = "Daily calcium carbonate D3 600 mg and vitamin D 800 IU from screening to end of study",
    notes          = paste(
      "KN012 arm of the population-PK analysis set (Zhu 2026 Results 3.5 and Supplementary Table 3).",
      "Twelve of the 139 KN012 participants were intensively sampled and the remainder sparsely",
      "sampled (Supplementary Table 1). Serum concentrations were measured by electrochemiluminescence",
      "on the Meso Scale Discovery platform. Estimation was by FOCE with interaction in NONMEM 7.4;",
      "the minimum success rate of 1000 bootstrap runs was 94.1%. Anti-drug antibodies were positive",
      "in 3 of 139 (2.2%) KN012 participants and ADA was screened but not retained as a covariate.",
      "Age, weight, BMI, BSA, baseline ALT, AST, albumin, total protein, total bilirubin, creatinine",
      "clearance and ADA status were all tested; only body weight on clearance was retained."
    )
  )

  ini({
    # ---- Structural PK: Zhu 2026 Supplementary Table 4, 'Final model' column.
    # Bioavailability was not estimated, so Vc and CL are apparent values
    # (Vc/F, CL/F) for the 60 mg subcutaneous route; the source table prints
    # them without the /F qualifier.
    lvc <- log(0.448)
    label("Apparent central volume of distribution (L)")  # Supplementary Table 4, row 'Vc, L' = 0.448 [RSE 7.8%, 95% CI 0.38-0.516]
    lcl <- log(0.0983)
    label("Apparent linear clearance at 55.5 kg (L/day)")  # Supplementary Table 4, row 'CL, L/day' = 0.0983 [RSE 3.6%, 95% CI 0.0914-0.105]
    lka <- log(0.0136)
    label("First-order subcutaneous absorption rate constant (1/day)")  # Supplementary Table 4, row 'KA, 1/day' = 0.0136 [RSE 2.9%, 95% CI 0.0128-0.0144]

    # ---- Parallel Michaelis-Menten elimination from the central compartment.
    # Vmax is carried verbatim from the source table, but its printed unit is
    # wrong: it is a CONCENTRATION elimination rate, not an amount rate. See
    # the note in model() for the arithmetic that establishes this.
    lvmax <- log(0.321)
    label("Maximum Michaelis-Menten elimination rate Vmax (ug/mL/day)")  # Supplementary Table 4, row 'Vmax, mg/day' = 0.321 [RSE 12.1%, 95% CI 0.245-0.397]; printed unit corrected to ug/mL/day, see model()
    lkm <- log(0.0549)
    label("Michaelis-Menten constant Km (ug/mL)")  # Supplementary Table 4, row 'Km, ug/mL' = 0.0549 [RSE 5.8%, 95% CI 0.0487-0.0611]

    # ---- Covariate effect.
    e_wt_cl <- 1.34
    label("Power exponent on (WT/55.5) for linear clearance (unitless)")  # Supplementary Table 4, row 'Weight on CL' = 1.34 [RSE 13.0%, 95% CI 0.999-1.68]; footnote 'CL=CL_TV(WT/55.5)^1.34'

    # ---- Inter-individual variability.
    # Zhu 2026 Supplementary Table 4 reports IIV as a percentage under the
    # symbol omega, i.e. as a coefficient of variation for log-normally
    # distributed parameters. Variances below are omega^2 = log(1 + CV^2).
    etalcl ~ 0.0573713  # Supplementary Table 4, row 'omega(CL),%' = 24.3% CV [RSE 8.3%]; shrinkage 14.9%
    etalvmax ~ 0.1567849  # Supplementary Table 4, row 'omega(Vmax), %' = 41.2% CV [RSE 7.7%]; shrinkage 4.3%

    # ---- Residual error. The KN012 final model carries a proportional term
    # only; the US reference-product model additionally carries an additive
    # term (Supplementary Table 5).
    propSd <- 0.248
    label("Proportional residual error (fraction)")  # Supplementary Table 4, row 'sigma (Prop), %' = 24.8% [RSE 2.0%, 95% CI 23.8-25.8]; shrinkage 11.6%
  })

  model({
    # Individual parameters. Body weight enters linear clearance only, as the
    # normalised power term printed in the Supplementary Table 4 footnote.
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
    # UNIT CORRECTION on Vmax. Supplementary Table 4 prints 'Vmax, mg/day',
    # which would make the saturable term vmax * Cc / (km + Cc) an amount rate
    # directly. Taken that way the model does NOT reproduce the paper's own
    # simulated exposures for a typical 60 kg subject given 60 mg
    # (Supplementary Table 6): it returns Cmax 3.59 vs 5.07 (-29%), AUC0-6mon
    # 147.1 vs 297.05 (-50%) and AUCinf 148.1 vs 298.91 (-50%). Treating Vmax
    # as a CONCENTRATION rate in ug/mL/day, so that the amount rate is
    # vmax * vc * Cc / (km + Cc), reproduces all three to within 0.22%
    # (5.07 / 297.2 / 299.6). The same correction reproduces all three
    # published values for the companion reference-product model to within
    # 0.16% using that model's own different Vc of 0.485 L, so the factor is
    # Vc and not a single fitted fudge constant. A change in Km cannot explain
    # the gap: no single Km value reproduces Cmax and AUCinf jointly. The
    # printed 'mg/day' is therefore a unit-tag slip and Vmax is the saturable
    # elimination rate on the concentration scale. See the vignette Errata.
    Cc <- central / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl * Cc - vmax * vc * Cc / (km + Cc)

    Cc ~ prop(propSd)
  })
}
