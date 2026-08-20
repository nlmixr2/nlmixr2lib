Lee_2023_tylosin_Pmultocida_infected <- function() {
  description <- paste0(
    "Preclinical (pig, crossbred Duroc x (Landrace x Yorkshire)). ",
    "Ex vivo sigmoidal Emax PK/PD-integration model for the antibacterial effect of tylosin ",
    "against Pasteurella multocida (challenge isolate BA1700127) in plasma drawn from co-infected pigs after a ",
    "single 20 mg/kg intramuscular dose. Lee 2023 Section 2.8 parameterises the drug effect over a ",
    "24 h ex vivo incubation as E = E0 - (Emax * C^gamma) / (C^gamma + EC50^gamma), where E is the ",
    "SIGNED log10(CFU/mL) difference between 0 and 24 h of incubation (positive = net growth, ",
    "negative = net kill), E0 is that same 24 h change in the drug-free control, and C is the ",
    "PK/PD index AUC24h/MIC. Parameters from Lee 2023 Table 3, column 'P. multocida / Infected': ",
    "Emax = 7.04 log10 CFU/mL, EC50 = 1.06 h, E0 = 3.41 log10 CFU/mL, gamma = 2.90. ",
    "The parameterisation was confirmed numerically against the paper's own independently printed ",
    "rows: Emax - E0 reproduces exactly (3.63), and solving E = 0 and E = -3 returns ",
    "AUC24h/MIC = 1.037 and 2.359 h versus the published 1.12 and 2.36 h. ",
    "The bacterial density bact (linear CFU/mL) is integrated as ",
    "d/dt(bact) = ln(10) * (E / 24) * bact so that log10(bact) changes by exactly E across each 24 h ",
    "window, reproducing the paper's model at the only times bacteria were counted. There is NO PK ",
    "component: Lee 2023 analysed the plasma concentrations non-compartmentally in WinNonlin 8.3 ",
    "(Section 2.7; Table 2 reports only Cmax, Tmax, T1/2, AUC, Vz/F, Cl/F and MRT) and published no ",
    "structural PK model, so exposure enters as the externally supplied covariate AUCMIC_TYLO. That ",
    "covariate carries the AUC24h/MIC RATIO directly because Lee 2023 never reports the challenge ",
    "isolates' own MIC -- only that strains 'with MIC values similar to the MIC90' were chosen -- so ",
    "no MIC can be sourced to split the ratio into numerator and denominator. Neither between-subject ",
    "variability nor a residual error magnitude was reported (the +/- values in Table 3 are standard ",
    "deviations of point estimates across the n = 3 animals whose plasma was used, not estimated ",
    "variance components), so there are no eta parameters and addSd is FIXED at 0; the model is ",
    "intended for typical-value simulation."
  )
  reference <- paste(
    "Lee E-B, Abbas MA, Park J, Tassew DD, Park S-C. (2023).",
    "Optimizing tylosin dosage for co-infection of Actinobacillus pleuropneumoniae",
    "and Pasteurella multocida in pigs using pharmacokinetic/pharmacodynamic modeling.",
    "Frontiers in Pharmacology 14:1258403.",
    "doi:10.3389/fphar.2023.1258403.",
    sep = " "
  )
  vignette <- "Lee_2023_tylosin"
  units <- list(
    time = "h",
    dosing = "h (tylosin AUC24h/MIC index, supplied as a covariate)",
    concentration = "log10 CFU/mL (observation)"
  )

  depends <- c("AUCMIC_TYLO")
  paper_specific_compartments <- c("bact")

  compartmentData <- list(
    bact = list(analyte = "Pasteurella multocida", units = "CFU/mL", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    AUCMIC_TYLO = list(
      description        = "Tylosin PK/PD index: plasma area under the concentration-time curve over 24 h divided by the MIC of the challenge isolate (AUC24h/MIC)",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Lee 2023 Section 2.8 defines the sigmoid Emax driver C as the AUC24h/MIC ratio, and Table 3 ",
        "reports the model's outputs (EC50 and the bacteriostatic / bactericidal thresholds) in those ",
        "same ratio units (h). The ratio is carried directly rather than as an absolute AUC divided by ",
        "a model MIC parameter, because Lee 2023 does not report the MIC of either challenge isolate: ",
        "Sections 2.4 and 3.1 state only that strains 'with MIC values similar to the MIC90' were ",
        "selected (BA2000013 for A. pleuropneumoniae, BA1700127 for P. multocida). Substituting the ",
        "collection MIC90s (16 and 32 ug/mL) is contradicted by the paper's own Table 3, whose ",
        "T > MIC, AUC/MIC and Cmax/MIC rows are mutually inconsistent with those values and with each ",
        "other. Set to 0 for drug-free control records so the sigmoid term vanishes and the predicted ",
        "24 h change reduces to E0. See the vignette Assumptions and deviations section."
      ),
      source_name        = "AUC24 h/MIC (Lee 2023 Section 2.8 equation; Table 3 rows 'AUC 24 h /MIC for bacteriostatic activity' and 'AUC 24 h /MIC for bactericidal activity')"
    )
  )

  population <- list(
    species             = "pig (crossbred Duroc x (Landrace x Yorkshire), male)",
    n_subjects          = 3L,
    n_studies           = 1L,
    age_range           = "5-6 weeks",
    weight_median       = "9.5 kg",
    weight_range        = "9.5 +/- 1.1 kg (mean +/- SD)",
    sex_female_pct      = 0,
    organism            = "Pasteurella multocida challenge isolate BA1700127 (Animal and Plant Quarantine Agency, Kimchen, Korea), selected in Section 3.1 as a strain 'with MIC values similar to the MIC90'; the isolate's own MIC is NOT reported. Collection-wide tylosin MIC50/MIC90 were 16/16 ug/mL for A. pleuropneumoniae (89 isolates) and 16/32 ug/mL for P. multocida (363 isolates); MBC was 32 ug/mL for both species; the EUCAST ECOFFinder ECOFF was 64 ug/mL for both",
    system              = "Ex vivo time-kill assay. Plasma collected from dosed pigs at 0, 0.25, 0.5, 0.75, 1, 2, 4, 6, 8, 12 and 24 h post-dose, pre-filtered through a 0.22 um membrane, then inoculated to a final density of 1 x 10^6 CFU/mL and incubated at 37 C with plate counts at 1, 2, 4, 8, 12 and 24 h (Section 2.6)",
    disease_state       = "Experimental intranasal co-infection 12 h before dosing with a 1 mL mixed suspension, 0.5 mL per naris, containing 2.0 x 10^9 CFU/mL A. pleuropneumoniae BA2000013 and 2.0 x 10^9 CFU/mL P. multocida BA1700127; confirmed by clinical signs of depression, coughing and mild dyspnoea, by fever of 40.1 +/- 1.1 C versus 37.2 +/- 1.5 C in healthy pigs with p less than 0.001, and by PCR amplification of apxIVA and kmt1 (Sections 2.3 and 3.4)",
    dose_range          = "20 mg/kg tylosin (50 mg/mL injectable solution) as a single intramuscular dose (Section 2.3)",
    regions             = "Republic of Korea (Kyungpook National University; Gyeongsangbuk-do Veterinary Service Laboratory, Daegu)",
    notes               = paste0(
      "The parent PK study enrolled 14 pigs, randomised 7 to the healthy (non-infected) group and 7 to ",
      "the co-infected group; the ex vivo PD assay used plasma from n = 3 pigs per group (Section 2.6). ",
      "Ethics approval PTB-2022-IACUC013-A (Petobio Clinical Institute). Group-level ",
      "non-compartmental PK for co-infected pigs (Table 2): Cmax 3.59 ug/mL, Tmax 0.25 h, ",
      "T1/2 1.83 h, AUC 10.46 h*ug/mL, Cl/F 1905.66 mL/h/kg. The paper reports four independent ",
      "pathogen x health-state fits, packaged here as Lee_2023_tylosin_Apleuropneumoniae_healthy, ",
      "Lee_2023_tylosin_Apleuropneumoniae_infected, Lee_2023_tylosin_Pmultocida_healthy and ",
      "Lee_2023_tylosin_Pmultocida_infected, all sharing the vignette Lee_2023_tylosin."
    )
  )

  ini({
    # =============================================================
    # Lee 2023 Table 3 -- sigmoid Emax PK/PD integration
    # column: P. multocida / Infected
    # =============================================================
    # Section 2.8 equation:
    #   E = E0 - (Emax * C^gamma) / (C^gamma + EC50^gamma)
    # E is the SIGNED log10(CFU/mL) difference between 0 and 24 h of
    # ex vivo incubation: it equals E0 (net growth of the drug-free
    # control) at zero exposure and falls to E0 - Emax at saturating
    # exposure. Table 3's separately printed 'E max - E o' row equals
    # Emax minus E0 exactly for all four fits, which is what fixes this
    # sign convention rather than the (garbled) Table 3 footnote.
    le0 <- log(3.41)
    label("Change in bacterial count in the drug-free control over 24 h E0 (log10 CFU/mL)")  # Lee 2023 Table 3, P. multocida / Infected, E 0 = 3.41 +/- SD
    lemax <- log(7.04)
    label("Maximum drug effect Emax, the span from control growth to maximum kill (log10 CFU/mL)")  # Lee 2023 Table 3, P. multocida / Infected, E max = 7.04 +/- SD
    lec50 <- log(1.06)
    label("AUC24h/MIC producing half the maximum effect EC50 (h)")  # Lee 2023 Table 3, P. multocida / Infected, EC 50 = 1.06 +/- SD
    lhill <- log(2.90)
    label("Hill coefficient gamma, steepness of the AUC24h/MIC effect curve (unitless)")  # Lee 2023 Table 3, P. multocida / Infected, gamma = 2.90 +/- SD

    # =============================================================
    # Starting bacterial density
    # =============================================================
    # FIXED experimental design input, not an estimated parameter.
    log10_cfu0 <- fixed(6)
    label("log10 starting bacterial density in the ex vivo plasma tube (log10 CFU/mL)")  # Lee 2023 Section 2.6: "5 uL of bacterial culture in the stationary phase was introduced into 0.5 mL of plasma, resulting in a final inoculum of 1 x 10^6 CFU/mL"

    # =============================================================
    # Residual error
    # =============================================================
    # Lee 2023 reported only the coefficient of determination of the
    # sigmoid fit (R^2 = 0.99; Discussion) and gave no residual standard
    # deviation, so the density-scale residual SD is held at zero for
    # deterministic typical-value simulation. See the vignette
    # Assumptions and deviations section.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (0; not reported in Lee 2023)")  # Lee 2023 reported R^2 only, no residual SD
  })

  model({
    e0   <- exp(le0)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Lee 2023 Section 2.8 sigmoid Emax equation. `effect` is the SIGNED
    # change in log10(CFU/mL) accrued over a 24 h ex vivo incubation
    # (positive = net growth, negative = net kill); it equals e0 at zero
    # exposure and approaches e0 - emax at saturating exposure.
    effect <- e0 - emax * AUCMIC_TYLO^hill / (AUCMIC_TYLO^hill + ec50^hill)

    # Lee 2023 counted bacteria only at the 0 and 24 h boundaries of the
    # ex vivo incubation and fitted the 24 h difference. Spreading that
    # difference uniformly across the window makes log10(bact) change by
    # exactly `effect` over 24 h, so the trajectory matches the paper's
    # model at every counted time point:
    #   d(log10 N)/dt = effect / 24
    #   => d(N)/dt    = ln(10) * (effect / 24) * N
    d/dt(bact) <- log(10) * (effect / 24) * bact
    bact(0) <- 10^log10_cfu0

    # log10 CFU/mL observation with a 1-CFU/mL floor (matches the
    # Wen 2016 / Chen 2023 in-vitro PD convention so the log10 stays
    # finite if bact is driven below 1 CFU/mL).
    Cc <- log10(bact + 1)
    Cc ~ add(addSd)
  })
}
