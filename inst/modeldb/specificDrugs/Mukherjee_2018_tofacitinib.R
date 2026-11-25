Mukherjee_2018_tofacitinib <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and a",
    "fixed absorption lag time for oral tofacitinib in adults with moderately",
    "to severely active ulcerative colitis from the 8-week phase 2",
    "dose-ranging induction trial A3921063 (NCT00787202; 0.5, 3, 10 or 15 mg",
    "twice daily). Parameterized in apparent oral clearance (CL/F) and",
    "apparent volume (V/F); no covariates were retained (dose, visit week,",
    "faecal calprotectin, C-reactive protein, albumin and total Mayo score",
    "were tested on CL/F and none was significant). Exponential IIV on CL/F",
    "and Ka, no IIV on V/F, single residual-error term. The model's purpose",
    "in the source was to supply individual Cav,ss = Dose / (12 * CL/F) and",
    "Ctrough,ss for the exposure-response models of the same paper",
    "(Mukherjee_2018_tofacitinib_<endpoint>_<metric>).",
    sep = " "
  )
  reference <- paste(
    "Mukherjee A, Hazra A, Smith MK, Martin SW, Mould DR, Su C, Niezychowski",
    "W. Exposure-response characterization of tofacitinib efficacy in",
    "moderate to severe ulcerative colitis: Results from a dose-ranging phase",
    "2 trial. Br J Clin Pharmacol. 2018;84(6):1136-1145.",
    "doi:10.1111/bcp.13523. PMID: 29377257. PMCID: PMC5980553.",
    sep = " "
  )
  vignette <- "Mukherjee_2018_tofacitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(
      analyte = "tofacitinib",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "tofacitinib",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    DOSE = list(
      description = "Assigned tofacitinib dose level (mg twice daily)",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F to assess dose proportionality; not significant (Results, 'Pharmacokinetics'). No estimate printed.",
      source_name = "Dose"
    ),
    CRP = list(
      description = "Baseline C-reactive protein",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F as a disease-severity marker; not significant (Results, 'Pharmacokinetics'). No estimate printed.",
      source_name = "CRP"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F; not significant (Results, 'Pharmacokinetics'). No estimate printed. Table 1 prints the unit as 'mg dl-1' but the values (mean 4.2) are g/dL.",
      source_name = "albumin"
    ),
    SCORE_MAYO = list(
      description = "Baseline total Mayo score (0-12)",
      units = "(score, 0-12)",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F; not significant (Results, 'Pharmacokinetics'). No estimate printed.",
      source_name = "total Mayo score"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 141L,
    n_studies = 1L,
    age_range = "19-77 years (all 188 analysed patients; Table 1)",
    age_mean = "42.1 years (SD 13.6), all 188 analysed patients (Table 1)",
    weight_range = "41.4-126 kg (Table 1)",
    weight_mean = "74.8 kg (SD 15.5), all 188 analysed patients (Table 1)",
    sex_female_pct = 45.2,
    disease_state = paste(
      "Moderately to severely active ulcerative colitis (total Mayo score",
      "6-12, endoscopic subscore >= 2, diagnosis >= 3 months); stable oral",
      "mesalamine and/or oral corticosteroids (<= 30 mg/day prednisone",
      "equivalent) allowed."
    ),
    dose_range = "Tofacitinib 0.5, 3, 10 or 15 mg orally twice daily for 8 weeks (randomized 2:2:2:3:3 with placebo).",
    regions = "Multinational phase 2 trial A3921063 (NCT00787202).",
    notes = paste(
      "141 tofacitinib-treated patients with plasma concentrations. Samples",
      "pre-dose and 0.25, 0.5, 1 and 2-3 h post-dose at the baseline and week",
      "8 visits were modelled; week 2 and week 4 samples were excluded for",
      "missing dosing information. LLOQ 0.100 ng/mL. Demographics in Table 1",
      "are for all 188 analysed patients (placebo included); 103/188 (54.8%)",
      "male. Fitted with NONMEM 7.3, FOCE."
    )
  )

  ini({
    lka <- log(2.83); label("Absorption rate constant Ka (1/h)")                   # Results 'Pharmacokinetics': Ka = 2.83 (SE 0.46) h-1
    lcl <- log(22.4); label("Apparent oral clearance CL/F (L/h)")                 # Results 'Pharmacokinetics': CL/F = 22.4 (SE 0.95) L/h
    lvc <- log(94.2); label("Apparent volume of distribution V/F (L)")            # Results 'Pharmacokinetics': V/F = 94.2 (SE 2.35) L
    ltlag <- fixed(log(0.16)); label("Absorption lag time (h)")                    # Results 'Pharmacokinetics': lag 'estimated to be 0.16 h ... fixed to this value to allow successful convergence'

    # IIV printed as %CV; converted with omega^2 = log(CV^2 + 1).
    etalcl ~ 0.094034                                                             # Results 'Pharmacokinetics': BSV CL/F 31.4% CV -> log(1 + 0.314^2)
    etalka ~ 0.568505                                                             # Results 'Pharmacokinetics': BSV Ka 87.5% CV -> log(1 + 0.875^2)

    propSd <- 0.475; label("Proportional residual error (fraction)")              # Results 'Pharmacokinetics': residual variability 47.5%
  })

  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    tlag <- exp(ltlag)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
