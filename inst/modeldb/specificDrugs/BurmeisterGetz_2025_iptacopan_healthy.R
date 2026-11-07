BurmeisterGetz_2025_iptacopan_healthy <- function() {
  description <- paste(
    "Two-site competitive equilibrium plasma-protein-binding model for iptacopan in adults with normal hepatic function (healthy controls of the CLNP023A2105 hepatic-impairment study).",
    "Saturable high-affinity binding to the pharmacological target complement factor B (capacity bmax, dissociation constant kd) plus linear non-saturable binding to other plasma proteins (maximum unbound fraction fumax).",
    "Closed-form, time-independent model: the total plasma iptacopan concentration is supplied as the time-varying covariate CP_IPTACOPAN_NGML and the model returns the unbound concentration Cu (the single meaningful root of the binding quadratic) and the fractional occupancy of available factor B, occ.",
    "The binding model was fitted separately to each hepatic-function group; this file carries the healthy-control parameter set."
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
        "Observed range in the source study: about 10 to 5,000 ng/mL after a single oral 200 mg dose (Figure 3); healthy-control geometric-mean Cmax 3,490 ng/mL (Table 2).",
        "Set to 0 for drug-free records (Cu and occ are then 0)."
      ),
      source_name = "Ctot"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 16L,
    n_studies = 1L,
    age_range = "31-73 years",
    age_median = "60.5 years",
    weight_range = "59.5-112.6 kg",
    weight_median = "86.8 kg",
    sex_female_pct = 43.8,
    race_ethnicity = c(White = 93.8, Other = 6.3),
    disease_state = "Healthy adults with normal hepatic function, each demographically matched (sex, age +/- 10 years, body weight +/- 15%, smoking status) to a participant with hepatic impairment",
    dose_range = "Single oral 200 mg iptacopan, fasted",
    regions = "Single clinical site (CLNP023A2105, NCT05078580)",
    hepatic_function = "Normal; median (range) albumin 44.5 (38-47) g/L, bilirubin 6.8 (3.4-22.3) umol/L, prothrombin INR 1.00 (0.9-1.0)",
    notes = paste(
      "Baseline demographics from Burmeister Getz 2025 Table 1 (healthy-control column).",
      "The binding model was fitted with minpack.lm::nlsLM and a proportional error model to paired total and unbound iptacopan concentrations from 22 post-dose samples per participant (15 min to 240 h; Supplementary Appendix, Statistical analyses and PK assessments).",
      "No inter-individual variability was estimated: the fit is a pooled per-group nonlinear regression."
    )
  )

  ini({
    # Burmeister Getz 2025 Table 3, healthy-control column. All three were
    # estimated (unconstrained nlsLM fit); RSEs are parameter precision,
    # not variability.
    lbmax <- log(1874.30)
    label("Factor B plasma concentration available for iptacopan binding, fB,avail (nM)") # Table 3 'fB,avail, nM' = 1874.30 (168.69 ug/mL at 90 kDa), RSE 3.85%
    lkd <- log(9.63)
    label("Equilibrium dissociation constant between iptacopan and factor B, KD (nM)") # Table 3 'KD, nM' = 9.63, RSE 6.29%
    fumax <- 0.25
    label("Maximum unbound fraction of plasma iptacopan, fu,max (unitless)") # Table 3 'fu,max' = 0.25, RSE 3.36%

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
