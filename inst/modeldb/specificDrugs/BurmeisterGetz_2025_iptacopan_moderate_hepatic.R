BurmeisterGetz_2025_iptacopan_moderate_hepatic <- function() {
  description <- paste(
    "Two-site competitive equilibrium plasma-protein-binding model for iptacopan in adults with moderate hepatic impairment (Child-Pugh class B, score 7-9; CLNP023A2105 hepatic-impairment study).",
    "Saturable high-affinity binding to the pharmacological target complement factor B (capacity bmax, dissociation constant kd) plus linear non-saturable binding to other plasma proteins (maximum unbound fraction fumax).",
    "Closed-form, time-independent model: the total plasma iptacopan concentration is supplied as the time-varying covariate CP_IPTACOPAN_NGML and the model returns the unbound concentration Cu (the single meaningful root of the binding quadratic) and the fractional occupancy of available factor B, occ.",
    "The binding model was fitted separately to each hepatic-function group; this file carries the moderate hepatic impairment parameter set."
  )
  reference <- paste(
    "Burmeister Getz E, Stein RR, Fink M, Kulmatycki K, Baltcheva I, Weis W, Shah B, Lawitz E, Schmouder R.",
    "Effect of Target-Mediated Disposition on Iptacopan Clinical Pharmacokinetics in Participants with Normal or Impaired Hepatic Function.",
    "Clin Pharmacol Ther. 2025;117(5):1358-1368. doi:10.1002/cpt.3559.",
    sep = " "
  )
  vignette <- "BurmeisterGetz_2025_iptacopan"
  units <- list(
    time = "h",
    dosing = "(none; binding model fed by an external total iptacopan plasma-concentration covariate)",
    concentration = "ng/mL (output Cu, unbound iptacopan; the binding parameters bmax and kd are in nM, and the model converts with the iptacopan molecular weight 422.52 g/mol)"
  )

  covariateData <- list(
    CP_IPTACOPAN_NGML = list(
      description = "Total (bound plus unbound) iptacopan plasma concentration at the time of each record, supplied as a time-varying covariate from observed plasma samples or an upstream PK source.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying per event row; drives the closed-form unbound-concentration expression (Burmeister Getz 2025 eq. 2 of the Supplementary Appendix, reproduced in the main text).",
        "Converted inside the model to nM with the iptacopan molecular weight 422.52 g/mol (Table 3 footnote b).",
        "Observed range in the source study: about 10 to 5,000 ng/mL after a single oral 200 mg dose (Figure 3); moderate hepatic impairment geometric-mean Cmax 3,230 ng/mL (Table 2).",
        "Set to 0 for drug-free records (Cu and occ are then 0)."
      ),
      source_name = "Ctot"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 8L,
    n_studies = 1L,
    age_range = "53-70 years",
    age_median = "61.5 years",
    weight_range = "65.5-112.8 kg",
    weight_median = "90.0 kg",
    sex_female_pct = 50.0,
    race_ethnicity = c(White = 100),
    disease_state = "Adults with moderate hepatic impairment (Child-Pugh class B; median (range) Child-Pugh score 7.5 (7-8))",
    dose_range = "Single oral 200 mg iptacopan, fasted",
    regions = "Single clinical site (CLNP023A2105, NCT05078580)",
    hepatic_function = "Moderate impairment, Child-Pugh B; median (range) albumin 38.0 (30-46) g/L, bilirubin 16.3 (6.9-29.1) umol/L, prothrombin INR 1.10 (1.0-1.2); hepatic encephalopathy grade 1 in 7/8; slight or moderate ascites in 7/8",
    notes = paste(
      "Baseline demographics from Burmeister Getz 2025 Table 1 (moderate-HI column).",
      "The binding model was fitted with minpack.lm::nlsLM and a proportional error model to paired total and unbound iptacopan concentrations from 22 post-dose samples per participant (15 min to 240 h; Supplementary Appendix, Statistical analyses and PK assessments).",
      "No inter-individual variability was estimated: the fit is a pooled per-group nonlinear regression."
    )
  )

  ini({
    # Burmeister Getz 2025 Table 3, moderate-HI column. All three were
    # estimated (unconstrained nlsLM fit); RSEs are parameter precision,
    # not variability.
    lbmax <- log(1403.16)
    label("Factor B plasma concentration available for iptacopan binding, fB,avail (nM)") # Table 3 'fB,avail, nM' = 1403.16 (126.28 ug/mL at 90 kDa), RSE 6.47%
    lkd <- log(8.00)
    label("Equilibrium dissociation constant between iptacopan and factor B, KD (nM)") # Table 3 'KD, nM' = 8.00, RSE 10.31%
    fumax <- 0.34
    label("Maximum unbound fraction of plasma iptacopan, fu,max (unitless)") # Table 3 'fu,max' = 0.34, RSE 4.69%

    # Supplementary Appendix, Statistical analyses: the fit used a
    # proportional error model whose magnitude is not reported.
    propSd <- fixed(0)
    label("Proportional residual error on unbound iptacopan (fraction; magnitude not reported)")
  })

  model({
    bmax <- exp(lbmax)
    kd <- exp(lkd)

    # Total iptacopan in nM; molecular weight 422.52 g/mol (Table 3
    # footnote b), so 1 ng/mL = 1000 / 422.52 nM.
    ctot_nm <- CP_IPTACOPAN_NGML * 1000 / 422.52

    # Binding model (Figure 2a; Supplementary Appendix eq. 1):
    #   Ctot = Cu + bmax * Cu / (kd + Cu) + Cu * (1 - fumax) / fumax
    # solved for Cu as the single meaningful (non-negative) root of the
    # quadratic (Supplementary Appendix eq. 2; main-text equation image):
    #   Cu = fumax * b / 2 + sqrt(fumax^2 / 4 * b^2 + fumax * Ctot * kd),
    #   b  = Ctot - kd / fumax - bmax
    bq <- ctot_nm - kd / fumax - bmax
    cu_nm <- fumax * bq / 2 + sqrt(fumax^2 / 4 * bq^2 + fumax * ctot_nm * kd)

    # Fractional occupancy of available factor B (Supplementary Appendix
    # eq. 3, as a fraction rather than a percentage); 1 - occ is the
    # fractional factor B activity plotted in Figure 4.
    occ <- cu_nm / (kd + cu_nm)

    Cu <- cu_nm * 422.52 / 1000
    Cu ~ prop(propSd)
  })
}
