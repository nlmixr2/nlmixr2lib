Rolsma_2025_tazobactam <- function() {
  description <- paste(
    "One-compartment population PK model with linear elimination for the",
    "tazobactam component of intravenously infused piperacillin-tazobactam in",
    "children and adults with cystic fibrosis (Rolsma 2025; 31 participants /",
    "34 enrollments, 107 plasma concentrations shared with piperacillin, ages",
    "5 to 54 years). No covariate reached significance: clearance is 7.67 L/h",
    "and the central volume 13.0 L for every subject. Strongly correlated",
    "interindividual variability on clearance and volume (correlation 0.93)",
    "plus interoccasion variability on the volume across 3 occasions, with a",
    "combined proportional plus additive residual error. The companion",
    "piperacillin model is Rolsma_2025_piperacillin.",
    sep = " "
  )
  reference <- paste(
    "Rolsma SL, Sokolow A, Patel PC, Sokolow K, Jimenez-Truque N, Fissell WH,",
    "Ryan V, Kirkpatrick CM, Nation RL, Gu K, Teresi M, Fishbane N, Kontos M,",
    "An G, Winokur P, Landersdorfer CB, Creech CB (2025). Population",
    "Pharmacokinetic Modeling of Cefepime, Meropenem, and",
    "Piperacillin-Tazobactam in Patients With Cystic Fibrosis. The Journal of",
    "Infectious Diseases 231(2):e364-e374. doi:10.1093/infdis/jiae451.",
    "Final parameter estimates from Supplemental Table 6D.",
    sep = " "
  )
  vignette <- "Rolsma_2025_betalactams_cysticfibrosis"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    OCC = list(
      description        = "Integer-valued occasion indicator for interoccasion variability on the volume of distribution",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Three occasions. Supplemental Table 6D caption: 'Parameter estimates",
        "of the final tazobactam model with IOV on V (n=3 occasions)'; the",
        "Supplemental Methods state 'The inclusion of inter-occasion",
        "variability on CL and V was evaluated for each compound' and the",
        "Results confirm the retained term: 'the interoccasion variability on",
        "volume of distribution was significant for the tazobactam model'.",
        "Note that this differs from the companion piperacillin model, which",
        "carries its IOV on clearance instead. The source reports a single IOV",
        "variance for V rather than one per occasion, which is the NONMEM",
        "'$OMEGA BLOCK(1) SAME' idiom; nlmixr2 has no SAME shortcut, so",
        "occasions 2 and 3 get their own etas with the variance fixed equal to",
        "the occasion-1 estimate, following the Jonsson_2011_ethambutol.R /",
        "Smythe_2013_gatifloxacin.R / Chen_2023_nemonoxacin.R precedents.",
        "Decomposed inside model() into mutually exclusive binary indicators",
        "oc1, oc2 and oc3. The paper does not define what physically separates",
        "the occasions; participants could enroll more than once (34",
        "enrollments from 31 participants per Supplemental Table 1) and",
        "sampling was opportunistic across a hospitalization, so an occasion",
        "is most plausibly a sampling occasion. For single-occasion simulation",
        "pass OCC = 1."
      ),
      source_name        = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened and NOT retained -- the headline negative finding for this",
        "model. Abstract: 'No covariates were identified for piperacillin and",
        "tazobactam.' Results 'PopPK Modeling': 'The elimination clearance and",
        "volume of distribution of piperacillin and tazobactam did not appear",
        "to correlate with biologically plausible covariates or any of the",
        "collected demographic characteristics.' Cohort 55.74 +/- 16.37 kg",
        "(median 54.05, range 17.6 to 93.3) per Supplemental Table 3."
      )
    ),
    FFM = list(
      description = "Fat-free mass by the Janmahasatian formula (reported by the source as lean body weight)",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened as part of the pre-defined base covariate set and not",
        "retained. Cohort 41.01 +/- 11.72 kg (median 38.25, range 13.6 to",
        "59.3) per Supplemental Table 3."
      )
    ),
    CRCL = list(
      description = "Estimated creatinine clearance computed on lean body weight (raw, not BSA-normalized)",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Screened as part of the pre-defined base covariate set and not",
        "retained, unlike in the companion cefepime and meropenem models where",
        "CLCR,LBW is the strongest covariate on clearance. Results 'PopPK",
        "Modeling': 'The elimination clearance and volume of distribution of",
        "piperacillin and tazobactam did not appear to correlate with",
        "biologically plausible covariates.'"
      )
    ),
    CFTRMOD = list(
      description = "CFTR modulator use",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened and explicitly rejected. Results 'PopPK Modeling': 'Other",
        "covariates, such as sample acquisition site, CFTR mutation, use of",
        "CFTR modulators ..., or complications of CF (eg, diabetes) did not",
        "have substantial impact on overall PK.' 35% of",
        "piperacillin-tazobactam enrollments reported any modulator use",
        "(Supplemental Table 2). No point estimate is reported, so no effect",
        "can be encoded."
      )
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened (Methods 'PopPK Model Development' covariate list) but not",
        "retained. Cohort 24.8 +/- 11.9 years (median 22.0, range 5 to 54) per",
        "Supplemental Table 3."
      )
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Screened (Methods covariate list) but not retained. Cohort 21.04 +/-",
        "4.36 kg/m^2 (median 20.15, range 15.4 to 33.1) per Supplemental",
        "Table 3."
      )
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = paste(
        "Screened (Methods covariate list) but not retained. Cohort 1.57 +/-",
        "0.31 m^2 (median 1.60, range 0.7 to 2.1) per Supplemental Table 3."
      )
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = paste(
        "Screened (Methods covariate list) but not retained. Cohort 161.20",
        "+/- 16.06 cm (median 164.00, range 99.1 to 186.0) per Supplemental",
        "Table 3."
      )
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened (Methods covariate list) but not retained. 53% female among",
        "piperacillin-tazobactam enrollments (18 of 34) per Supplemental",
        "Table 2."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 31L,
    n_studies      = 1L,
    n_enrollments  = 34L,
    age_range      = "5 to 54 years",
    age_median     = "22.0 years (mean 24.8 +/- 11.9)",
    weight_range   = "17.6 to 93.3 kg total body weight",
    weight_median  = "54.05 kg total body weight (mean 55.74 +/- 16.37); lean body weight median 38.25 kg (mean 41.01 +/- 11.72)",
    sex_female_pct = 100 * 18 / 34,
    race_ethnicity = c(White = 91, `Multi-Racial` = 9, `Hispanic or Latino (ethnicity)` = 6),
    disease_state  = paste(
      "Cystic fibrosis, admitted for a pulmonary exacerbation or for",
      "microbial eradication therapy. 97% of the piperacillin-tazobactam",
      "group carried at least one copy of the DF508 mutation. 35% reported",
      "any CFTR modulator use. Participants receiving piperacillin-tazobactam",
      "had more CF-related complications than the cefepime group (CF-related",
      "diabetes, hepatic disease, pancreatic insufficiency, and 10 or more",
      "pulmonary exacerbations in the previous 2 years; Supplemental Table",
      "4). Staphylococcus aureus and Pseudomonas aeruginosa were the most",
      "common sputum isolates. Participants with a history of solid-organ or",
      "hematologic transplantation were excluded."
    ),
    dose_range     = paste(
      "Intravenous infusion of the piperacillin-tazobactam combination",
      "product, standard of care as ordered by the treating team. The paper",
      "reports doses in terms of the piperacillin component: the majority",
      "(n = 454) were 4000 mg piperacillin and the remaining 190 doses 1867",
      "to 3556 mg piperacillin, with most doses in participants under 17",
      "years being 86 to 106 mg/kg piperacillin. The commercial product is",
      "formulated at an 8:1 piperacillin:tazobactam ratio, so a 4000 mg",
      "piperacillin dose corresponds to 500 mg tazobactam; the paper does not",
      "state the tazobactam amounts explicitly, so dose amounts must be",
      "supplied by the user when simulating this model. Infusion durations in",
      "hours were 0.5 (108 doses), 3 (6 doses), 4 (524 doses) and 6 (6",
      "doses)."
    ),
    regions        = "United States (Vanderbilt University Medical Center, Nashville TN; University of Iowa Hospital, Iowa City IA).",
    notes          = paste(
      "Opportunistic sampling during hospitalization, January 2018 to March",
      "2020. 107 of the 667 total plasma samples in the study were",
      "piperacillin/tazobactam. Assay LC-MS/MS, linear 0.25 to 150 mg/L for",
      "tazobactam. NONMEM 7.4.3 with FOCE+I; outliers removed per the FDA",
      "popPK guidance (February 2022) at a base-model weighted residual",
      "greater than 5. One- and two-compartment models were evaluated and the",
      "one-compartment model with linear elimination was selected. Drug",
      "administration was modeled as a zero-order (infusion) process. Model",
      "evaluation by prediction-corrected VPC (Supplemental Figure 2).",
      "Baseline demographics from Supplemental Tables 2 and 3; final",
      "parameter estimates from Supplemental Table 6D. Unlike the other three",
      "drugs in the paper, tazobactam has no PK/PD target-attainment analysis",
      "-- it is a beta-lactamase inhibitor rather than an antibacterial, so",
      "the published fT>MIC targets and breakpoints apply to piperacillin",
      "only. Supplemental Table 6D reports a 'Standard error of population",
      "estimate' column header where Tables 6A-6C report '%RSE of population",
      "estimate'; the values (10.6%, 19.3%) are in percent and are read here",
      "as relative standard errors on the same basis as the other three",
      "tables. Piperacillin and tazobactam were fitted as two separate models",
      "with different interoccasion-variability structures (IOV on CL for",
      "piperacillin, IOV on V for tazobactam) and are packaged as two model",
      "files accordingly."
    )
  )

  ini({
    # Structural parameters -- Rolsma 2025 Supplemental Table 6D, 'Population
    # estimate' column. No covariate was retained for tazobactam, so these are
    # the typical values for every subject.
    #
    # Supplemental Table 6D verbatim (Parameter / Unit / Population estimate):
    #   CL   / L/h  / 7.67  (SE 10.6%, IIV 55%, correlation (CL,V) 93%)
    #   V    / L    / 13.0  (SE 19.3%, IIV 79%, IOV 54%)
    #   CVCP / -    / 43%
    #   SDCP / mg/L / 0.548
    lcl <- log(7.67); label("Clearance (L/h)")                     # Suppl Table 6D row CL = 7.67 L/h (SE 10.6%)
    lvc <- log(13.0); label("Central volume of distribution (L)")  # Suppl Table 6D row V  = 13.0 L   (SE 19.3%)

    # Interindividual variability. The Supplemental Methods state the IIV model
    # explicitly ('The inter-individual variability (IIV) of the PK parameters
    # for all antibiotics in this study was estimated using an exponential
    # model'), so the tabulated percentages are CV% on the log-normal scale and
    # convert as omega^2 = log(CV^2 + 1).
    # cov = 0.93 * sqrt(0.264285 * 0.484954) = 0.332943.
    etalcl + etalvc ~ c(0.264285,
                        0.332943, 0.484954)  # Suppl Table 6D: IIV CL 55% CV -> log(0.55^2+1) = 0.264285 (RSE 39%, shrinkage 13%); IIV V 79% CV -> 0.484954 (RSE 55%, shrinkage 19%); correlation 0.93

    # Interoccasion variability on the log volume of distribution across the 3
    # occasions of Supplemental Table 6D's caption. Note the contrast with the
    # companion piperacillin model, whose IOV sits on clearance. The source
    # reports a single IOV variance for V (54% CV), i.e. the NONMEM
    # '$OMEGA BLOCK(1) SAME' idiom; nlmixr2 has no SAME shortcut, so occasions
    # 2 and 3 carry their own etas with the variance fixed equal to the
    # occasion-1 estimate.
    etaiov_vc_1 ~ 0.255882       # Suppl Table 6D: IOV V 54% CV -> log(0.54^2+1) = 0.255882, occasion 1 (estimated)
    etaiov_vc_2 ~ fix(0.255882)  # SAME-equivalent: occasion 2 shares the occasion-1 IOV variance
    etaiov_vc_3 ~ fix(0.255882)  # SAME-equivalent: occasion 3 shares the occasion-1 IOV variance

    # Residual error. Supplemental Methods: 'To describe the residual
    # variability (RV) of each antibiotic, a combined additive and proportional
    # error model was utilized.'
    propSd <- 0.43;  label("Proportional residual error (fraction)")  # Suppl Table 6D row CVCP = 43%     (RSE 37%)
    addSd  <- 0.548; label("Additive residual error (mg/L)")          # Suppl Table 6D row SDCP = 0.548 mg/L (RSE 127%)
  })

  model({
    # 1. Decompose the integer occasion column into mutually exclusive binary
    # indicators multiplexing the IOV etas on the log volume of distribution.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2 + oc3 * etaiov_vc_3

    # 2. Individual PK parameters. No covariate was retained for tazobactam
    # (Abstract: 'No covariates were identified for piperacillin and
    # tazobactam'), so the typical values apply to every subject and all
    # between-subject differences come from the correlated IIV plus the
    # per-occasion IOV on the volume of distribution.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc + iov_vc)

    # 3. Micro-constant.
    kel <- cl / vc

    # 4. One-compartment disposition with linear elimination. Tazobactam is
    # given by intravenous infusion as part of the piperacillin-tazobactam
    # combination product, so dose records enter the central compartment
    # directly as a zero-order input (Methods: 'Since all antibiotics were
    # administered by intravenous infusions at a constant rate, drug
    # administration was modeled as a zero-order process'); there is no depot.
    # Dose amounts are of tazobactam, NOT of the combination product and not of
    # piperacillin -- at the commercial 8:1 ratio a 4000 mg piperacillin dose
    # is 500 mg tazobactam.
    d/dt(central) <- -kel * central

    # 5. Observation. Total (not unbound) plasma tazobactam in mg/L, matching
    # the LC-MS/MS assay the model was fitted to. Combined proportional plus
    # additive residual error.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
