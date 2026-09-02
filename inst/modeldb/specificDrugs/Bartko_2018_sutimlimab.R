Bartko_2018_sutimlimab <- function() {
  description <- paste(
    "Population pharmacodynamic sigmoidal inhibitory Emax (Imax) model for sutimlimab (BIVV009 / TNT009), a humanized monoclonal IgG4 antibody against complement factor C1s, describing knockdown of classical complement pathway (CP) activity in healthy volunteers.",
    "The direct-effect model has no delay component: Bartko 2018 first checked for hysteresis, found none, and therefore time-matched each individual sutimlimab serum concentration with the CP activity measured in the same sample.",
    "PD-only model: the sutimlimab serum concentration is supplied as a time-varying covariate CP_SUTIMLIMAB_UGML (ug/mL).",
    "The source publication characterised sutimlimab PK by model-independent NCA only (Cmax / tmax / AUC / half-life in Table 2) and did not develop a structural population PK model, so the PD model has no coupled PK component. The authors note the PK is nonlinear below about 100 ug/mL, consistent with target-mediated elimination.",
    "Population: 48 healthy volunteers who received sutimlimab in a phase I, first-in-human, double-blind, randomized, placebo-controlled single- (part A) and multiple- (part B) ascending-dose trial (NCT02502903)."
  )
  reference <- paste(
    "Bartko J, Schoergenhofer C, Schwameis M, Firbas C, Beliveau M, Chang C,",
    "Marier JF, Nix D, Gilbert JC, Panicker S, Jilma B (2018).",
    "A Randomized, First-in-Human, Healthy Volunteer Trial of sutimlimab, a Humanized Antibody",
    "for the Specific Inhibition of the Classical Complement Pathway.",
    "Clin Pharmacol Ther 104(4):655-663.",
    "doi:10.1002/cpt.1111.",
    sep = " "
  )
  vignette <- "Bartko_2018_sutimlimab"
  units <- list(
    time = "h",
    dosing = "(none; PD-only model fed by an external sutimlimab serum-concentration covariate)",
    concentration = "(observation CPactivity is classical complement pathway activity expressed as a percentage of the assay's normal reference, %; the driving covariate CP_SUTIMLIMAB_UGML is in ug/mL)"
  )

  covariateData <- list(
    CP_SUTIMLIMAB_UGML = list(
      description        = "Instantaneous sutimlimab (BIVV009 / TNT009) serum concentration at the time of each CP-activity observation, supplied as a time-varying covariate from observed serum samples or an upstream PK source.",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying per event row. Drives the sigmoidal Imax expression",
        "CPactivity = e0 - imax * CP_SUTIMLIMAB_UGML^hill / (ic50^hill + CP_SUTIMLIMAB_UGML^hill).",
        "In Bartko 2018 this was the individual serum sutimlimab concentration measured by a validated immunoassay at a GLP-certified laboratory (Vela Laboratories, Vienna, Austria; Bartko 2018 Methods, 'Pharmacokinetics'), time-matched to the CP activity measured in the same sample (Methods, 'PK/PD').",
        "Reference values observed (Bartko 2018 Table 2): mean Cmax rose from 40 ug/mL at 3 mg/kg to 2036 ug/mL at 100 mg/kg after a single 60-minute i.v. infusion; concentrations were below the limit of quantification at the 0.3 and 1 mg/kg dose levels, so those two cohorts contribute no concentration-effect pairs.",
        "Set to 0 outside the drug-exposure window (the inhibition term then collapses to 0 and CPactivity returns to e0)."
      ),
      source_name        = NA_character_
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 48L,
    n_studies      = 1L,
    age_range      = "19-59 years (part A sutimlimab arms, median 32); 22-41 years (part B sutimlimab arms, median 27)",
    age_median     = "32 years (part A sutimlimab); 27 years (part B sutimlimab)",
    weight_range   = "not tabulated as a range; cohort means (SD) ran 55 (3) to 79 (4) kg. Protocol excluded body weight > 98 kg, and > 58 kg for the 100 mg/kg part A cohort",
    weight_median  = NA_character_,
    sex_female_pct = 45.8,
    race_ethnicity = c(Caucasian = 95.8, African = 4.2),
    disease_state  = "healthy volunteers (no complement-mediated disorder); all vaccinated against encapsulated bacterial pathogens before dosing",
    dose_range     = paste(
      "Part A (single ascending dose): 0.3, 1, 3, 10, 30, 60, or 100 mg/kg sutimlimab or placebo as a single ~60-minute i.v. infusion, 3:1 active:placebo (n = 3 active for 0.3 and 1 mg/kg, n = 6 active for the remainder).",
      "Part B (multiple ascending dose): 30 or 60 mg/kg sutimlimab or placebo as four once-weekly ~60-minute i.v. infusions (n = 6 active per dose level), with a 2-week follow-up observation period."
    ),
    regions        = "Single centre: Department of Clinical Pharmacology, Medical University of Vienna, Austria",
    notes          = paste(
      "Baseline demographics are Bartko 2018 Table 1; 48 of the 64 enrolled volunteers received sutimlimab (36 in part A, 12 in part B) and 16 received placebo.",
      "Part C of the integrated protocol (patients with a complement-mediated disorder) was ongoing and is not reported in this publication, so this model is fit to healthy volunteers only.",
      "Baseline CP activity was normal in every subject (part A: placebo 95% +/- 10%, sutimlimab 97% +/- 14%; part B: placebo 99% +/- 7%, sutimlimab 94% +/- 18%).",
      "CP activity was measured semiquantitatively in serum with the commercial Complement System Classical Pathway WIESLAB enzyme immunoassay (Euro Diagnostica, Malmo, Sweden).",
      "PK/PD modelling was performed with Phoenix NLME v7 (Certara, Princeton, NJ); Bartko 2018 Methods, 'PK/PD'.",
      "The trial is registered as NCT02502903."
    )
  )

  ini({
    # Sigmoidal inhibitory Emax (Imax) model fit by Bartko 2018 to
    # time-matched pairs of individual serum sutimlimab concentration and
    # classical complement pathway (CP) activity pooled across parts A and B.
    # All four typical values are Supplementary Table S3 ("PK/PD parameters of
    # BIVV009 and CP activity - parts A and B"); the parenthesised figures in
    # that table are relative standard errors of the estimate (RSE%), i.e.
    # parameter precision, and are NOT interindividual variability.

    le0 <- log(94.8)
    label("Baseline classical complement pathway activity in the absence of sutimlimab (E0, % of the assay's normal reference)")
    # Bartko 2018 Supplementary Table S3: E0 (%) = 94.8 (RSE 1.1%).
    # Consistent with the observed baselines quoted in Results,
    # 'Pharmacodynamics' (part A sutimlimab 97 +/- 14%; part B sutimlimab
    # 94 +/- 18%).

    limax <- log(90.2)
    label("Maximum reduction in CP activity attributable to sutimlimab (Imax, percentage points of CP activity)")
    # Bartko 2018 Supplementary Table S3: Imax (%) = 90.2 (RSE 1.1%).
    # NOTE ON PARAMETERISATION: Imax enters ADDITIVELY, in the same units as
    # E0, so the saturating (high-concentration) asymptote is E0 - Imax =
    # 94.8 - 90.2 = 4.6%. It is NOT a fraction multiplying E0 (which would
    # asymptote at 94.8 * (1 - 0.902) = 9.29%). The published fitted curve in
    # Figure 5 settles this: it spans 94.5% at its left end down to a
    # high-concentration plateau of 4.4%, a span of 90.1 percentage points
    # that reproduces Imax = 90.2 directly. See the validation vignette,
    # section 'Which Imax parameterisation?', for the digitisation evidence.

    lic50 <- log(6.2)
    label("Sutimlimab serum concentration producing half of the maximum CP-activity reduction (IC50, ug/mL)")
    # Bartko 2018 Supplementary Table S3: IC50 (ug/mL) = 6.2 (RSE 27.5%);
    # also quoted in Results, 'Pharmacokinetic/pharmacodynamic (PD)
    # correlations': "a 50% knockdown of CP activity (IC50) predicted at a
    # sutimlimab concentration of 6.2 ug/mL".

    lhill <- log(2.4)
    label("Hill coefficient of the sutimlimab concentration / CP-activity inhibition curve (unitless)")
    # Bartko 2018 Supplementary Table S3: H = 2.4 (RSE 19.9%); also quoted in
    # the Abstract ("a Hill coefficient of 2.4") and in Results.
    # Cross-check: the paper's separately reported IC90 of 15.5 ug/mL is the
    # concentration giving 90% of Imax, i.e. IC50 * 9^(1/H) =
    # 6.2 * 9^(1/2.4) = 15.49 ug/mL. That identity over-determines IC50 and H
    # jointly and confirms both transcriptions.

    # Residual unexplained variability. Bartko 2018 Supplementary Table S3
    # reports only the four typical values above, each with its RSE%; no
    # residual error model, no sigma, and no interindividual variability are
    # tabulated anywhere in the paper or its five supplementary files.
    # Encoded as fixed(0) rather than invented: the model therefore predicts
    # typical values only. See the validation vignette's 'Assumptions and
    # deviations' section.
    addSd <- fixed(0)
    label("Additive residual error on CP activity (percentage points; not reported by Bartko 2018)")
  })

  model({
    e0 <- exp(le0)
    imax <- exp(limax)
    ic50 <- exp(lic50)
    hill <- exp(lhill)

    # Direct-effect sigmoidal Imax inhibition of classical complement pathway
    # activity. CP_SUTIMLIMAB_UGML is supplied per event row (ug/mL); the
    # response is CP activity as a percentage of the assay's normal reference,
    # falling from e0 at zero concentration toward the plateau e0 - imax.
    CPactivity <- e0 - imax * CP_SUTIMLIMAB_UGML^hill / (ic50^hill + CP_SUTIMLIMAB_UGML^hill)
    CPactivity ~ add(addSd)
  })
}
