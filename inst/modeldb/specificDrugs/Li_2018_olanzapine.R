Li_2018_olanzapine <- function() {
  description <- "Two-compartment population PK model for oral olanzapine with first-order absorption and an absorption lag time, developed in Han Chinese healthy male volunteers (single 10 mg dose, rich sampling) and adults with schizophrenia on routine therapeutic drug monitoring (Li 2018). Apparent central volume scales with body weight as a power function around a 60.59 kg reference; three olanzapine products from different manufacturers ('formulation #0', '#1' and '#2') carry separate effects on ka, Vc/F and CL/F, with formulation #2 as the reference. Between-subject variability is correlated between CL/F and Vc/F."
  reference <- paste(
    "Li A, Ji S, Yue W, Yan H, Dong F, Ruan C, Li W, Lu W, Zhang D, Wang C.",
    "Development of a population pharmacokinetic model of olanzapine for Chinese",
    "health volunteers and patients with schizophrenia.",
    "BMJ Open. 2018;8(8):e020070. doi:10.1136/bmjopen-2017-020070.",
    "Structural and covariate model: Results Equations (5)-(9);",
    "parameter estimates and non-parametric bootstrap: Table 2.",
    sep = " "
  )
  vignette <- "Li_2018_olanzapine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters Vc/F as an estimated power function normalised to 60.59 kg",
        "(Results Equations 5-7: Vc/F = 2390 * (WT/60.59)^0.579 * ...). The 60.59 kg",
        "reference is printed in the equations and equals the total-cohort mean body",
        "weight of Table 1 (60.6 +/- 11.1 kg); cohort medians are 62.5 kg (cohort A)",
        "and 60.0 kg (cohort B), overall range 36.0-98.0 kg. The exponent 0.579",
        "(RSE 25.7%) was estimated, not fixed to an allometric 1 or 0.75. Body weight",
        "was not retained on CL/F. Treated as a baseline, time-fixed covariate: Table 1",
        "reports a single weight per subject."
      ),
      source_name = "WT"
    ),
    FORM_OLZ_F0 = list(
      description = "Olanzapine product 'formulation #0' indicator (1 = the dose was given as formulation #0, 0 = any other product)",
      units = "(binary)",
      type = "binary",
      reference_category = "0; FORM_OLZ_F0 = FORM_OLZ_F1 = 0 selects formulation #2, the model reference",
      notes = paste(
        "Per-dose-record indicator. The paper's NONMEM covariate DRUG0 (Table 2 rows",
        "'DRUG0 ON Ka', 'DRUG0 ON Vc/F', 'DRUG0 ON CL/F'). Olanzapine from three",
        "manufacturers was used; the products are labelled 'formulation #0', '#1' and",
        "'#2' and all are oral tablets (Discussion). Cohort A (22 healthy male",
        "volunteers) received formulations #0 and #1 in a two-period crossover with a",
        "3-week washout, so this indicator varies within subject across periods;",
        "cohort B (234 patients) received only formulation #2.",
        "IMPORTANT CONFOUNDING: because formulation #2 was given exclusively to the",
        "patient cohort and formulations #0/#1 exclusively to the healthy-volunteer",
        "cohort, the formulation effects on Vc/F and CL/F are completely confounded with",
        "cohort, study design and sampling richness. The authors state this explicitly",
        "('the influences of formulation and population on PK parameters were mixed and",
        "could not be distinguished due to data limitations', Discussion). Read this",
        "column as a product-and-population stratum indicator, not as a bioequivalence",
        "statement about the tablets."
      ),
      source_name = "DRUG0"
    ),
    FORM_OLZ_F1 = list(
      description = "Olanzapine product 'formulation #1' indicator (1 = the dose was given as formulation #1, 0 = any other product)",
      units = "(binary)",
      type = "binary",
      reference_category = "0; FORM_OLZ_F0 = FORM_OLZ_F1 = 0 selects formulation #2, the model reference",
      notes = paste(
        "Per-dose-record indicator. The paper's NONMEM covariate DRUG1 (Table 2 rows",
        "'DRUG1 ON Vc/F', 'DRUG1 ON CL/F'). There is deliberately no 'DRUG1 ON Ka' term:",
        "the final model gives formulation #0 its own ka and lets formulations #1 and #2",
        "share the typical value, which the Discussion justifies on an objective-function",
        "comparison (OFV 4126.82 vs 4133.16). There is likewise no FORM_OLZ_F2 column --",
        "formulation #2 is the reference and is selected when both indicators are 0.",
        "The same cohort/formulation confounding described under FORM_OLZ_F0 applies."
      ),
      source_name = "DRUG1"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened (Methods, covariate model) but not retained. Total cohort median 29.5 years (range 18.0-48.0; mean 30.4 +/- 7.9), Table 1."
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      notes = "Screened but not retained. Total cohort median 165.0 cm (range 145.0-190.0; mean 164.8 +/- 8.1), Table 1."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened but not retained. Total cohort median 21.9 kg/m^2 (range 15.0-36.2; mean 22.3 +/- 3.4), Table 1."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened but not retained. 134 men / 122 women overall; cohort A was entirely",
        "male (22/0) and cohort B was 112 men / 122 women (Table 1). The paper reports",
        "'Male/female', so a SEXM-oriented source column needs value inversion.",
        "The Discussion contrasts this null result with two earlier analyses that found",
        "male elimination 30-38% higher than female, and attributes the difference to the",
        "all-male cohort A, to the Han Chinese population, and to data limitations."
      )
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained. Total cohort median 17.0 U/L (range 0.7-182.0; mean 23.1 +/- 20.2), Table 1. Subjects with severe liver-function abnormality were excluded; the Discussion states the observed transaminase range did not affect the PK parameters."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained. Total cohort median 20.9 U/L (range 2.0-147.0; mean 23.8 +/- 13.5), Table 1."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened but not retained. Total cohort median 3.9 mmol/L (range 1.1-9.8; mean 4.2 +/- 1.5), Table 1. No subject had obviously abnormal kidney function (Discussion)."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "olanzapine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "olanzapine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "olanzapine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 256,
    n_studies = 2,
    age_range = "18.0-48.0 years",
    age_median = "29.5 years",
    weight_range = "36.0-98.0 kg",
    weight_median = "60.0 kg",
    height_range = "145.0-190.0 cm",
    sex_female_pct = 47.7,
    race_ethnicity = c(Asian = 100),
    disease_state = "Cohort A: healthy male volunteers. Cohort B: adults with schizophrenia (DSM-IV-TR structured clinical interview diagnosis) on routine therapeutic drug monitoring.",
    dose_range = "Cohort A: a single 10 mg oral dose in each of two crossover periods (formulation #0 then formulation #1, or the reverse, separated by a 3-week washout). Cohort B: therapeutic oral doses titrated over the first 2 weeks and then held constant; the paper does not tabulate the administered doses.",
    regions = "China (Han Chinese throughout). Cohort A: Beijing Anding Hospital phase I clinical trial unit. Cohort B: a multicentre therapeutic-drug-monitoring study run from five main centres; the abstract cites 12 hospitals as contributing the olanzapine subset.",
    notes = paste(
      "Registry ChiCTR-TRC-10000934. Baseline demographics are Table 1; the values",
      "recorded above are the Total column. Cohort A contributed 616 concentrations",
      "from 22 healthy male volunteers sampled at 1, 2, 3, 4, 6, 8, 12, 24, 36, 48, 72,",
      "96, 120 and 168 h after each single dose (February-June 2001). Cohort B",
      "contributed 458 concentrations from 234 patients (May 2010-December 2011) at the",
      "end of treatment weeks 4 and 6; all but 19 were steady-state samples, and most",
      "were drawn at 06:00 with the remainder at 22:00. Patients were not allowed liver",
      "enzyme inducers or inhibitors (rifampin, warfarin, carbamazepine, phenobarbital,",
      "phenytoin, St John's Wort) for 2 weeks before or during the study.",
      "Assay: validated HPLC-MS/MS, calibrated 2-400 ng/mL.",
      "Estimation used NONMEM 7 level 2 with FOCE-I and the ADVAN4 subroutine.",
      "Qualification: goodness-of-fit plots (Figure 1), a 500-replicate non-parametric",
      "bootstrap with a 98.8% success rate (Table 2), normalised prediction distribution",
      "errors from 1000 simulations (Figure 2) and a prediction-corrected VPC (Figure 3).",
      "Smoking status and genotype were NOT collected and so could not be screened,",
      "which the authors list as the study's main limitation given that smoking is a",
      "known inducer of olanzapine metabolism."
    )
  )

  ini({
    # Structural parameters -- Table 2 'Estimate' column. The typical values are
    # those of formulation #2 (the reference product); Results Equations (7) and
    # (9) print the reference forms Vc/F = 2390 * (WT/60.59)^0.579 * exp(eta) and
    # CL/F = 25.4 * exp(eta) with no formulation multiplier.
    lka <- log(2.85)
    label("Absorption rate constant ka for formulations #1 and #2 (1/h)") # Table 2: Ka = 2.85 1/h (RSE 8.7%; bootstrap mean 2.24, 95% CI 1.21-3.87)
    lcl <- log(25.4)
    label("Apparent oral clearance CL/F for formulation #2 (L/h)") # Table 2: CL/F = 25.4 L/h (RSE 3.6%; bootstrap mean 25.4, 95% CI 23.7-27.4); Eq. 9
    lvc <- log(2390)
    label("Apparent central volume Vc/F for formulation #2 at 60.59 kg (L)") # Table 2: Vc/F = 2390 L (RSE 12.6%; bootstrap mean 2377.8, 95% CI 1479-3690); Eq. 7
    lq <- log(8.41)
    label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2: Q/F = 8.41 L/h (RSE 13.3%; bootstrap mean 9.80, 95% CI 6.54-16.20)
    lvp <- log(168)
    label("Apparent peripheral volume Vp/F (L)") # Table 2: Vp/F = 168 L (RSE 11.7%; bootstrap mean 175.9, 95% CI 135-239)
    ltlag <- log(0.877)
    label("Absorption lag time (h)") # Table 2: ALAG1 = 0.877 h (RSE 4.7%; bootstrap mean 0.800, 95% CI 0.672-0.928); Abstract and Results: 'The typical absorption time delay was 0.877 hour'

    # Body-weight effect on Vc/F. Results Equations (5)-(7) all carry the same
    # (WT/60.59)^0.579 term, so the exponent is shared across the three products
    # and the 60.59 kg reference is printed in the equations themselves.
    e_wt_vc <- 0.579
    label("Power exponent on (WT/60.59) for Vc/F (unitless)") # Table 2: 'WT ON Vc' = 0.579 (RSE 25.7%; bootstrap mean 0.590, 95% CI 0.191-1.06); Eqs. 5-7

    # Formulation effects. Table 2 reports these as multipliers on the typical
    # value, and Results Equations (5), (6) and (8) show them entering
    # multiplicatively; they are stored here as log ratios so they add on the log
    # scale inside model(). Formulation #2 is the reference (multiplier 1) and has
    # no row in Table 2 and no term in Equations (7) and (9).
    e_form_olz_f0_ka <- log(1.89)
    label("Log ratio of the formulation #0 absorption rate constant to the reference (unitless; exp = 1.89)") # Table 2: 'DRUG0 ON Ka' = 1.89 (RSE 22.0%; bootstrap mean 2.00, 95% CI 1.23-2.88). 2.85 * 1.89 = 5.39 1/h, the upper end of the Abstract's '2.85 h-1 to 5.39 h-1' range. Formulation #1 shares the typical ka; Table 2 has no 'DRUG1 ON Ka' row.
    e_form_olz_f0_vc <- log(0.216)
    label("Log ratio of the formulation #0 apparent central volume to the reference (unitless; exp = 0.216)") # Table 2: 'DRUG0 ON Vc/F' = 0.216 (RSE 13.1%; bootstrap mean 0.220, 95% CI 0.138-0.359); Eq. 5 prints the 2-significant-figure rounding 0.22
    e_form_olz_f1_vc <- log(0.207)
    label("Log ratio of the formulation #1 apparent central volume to the reference (unitless; exp = 0.207)") # Table 2: 'DRUG1 ON Vc/F' = 0.207 (RSE 12.8%; bootstrap mean 0.210, 95% CI 0.128-0.328); Eq. 6 prints the 2-significant-figure rounding 0.21
    e_form_olz_f0_cl <- log(0.610)
    label("Log ratio of the formulation #0 apparent oral clearance to the reference (unitless; exp = 0.610)") # Table 2: 'DRUG0 ON CL/F' = 0.610 (RSE 6.6%; bootstrap mean 0.610, 95% CI 0.537-0.696); Eq. 8
    e_form_olz_f1_cl <- log(0.610)
    label("Log ratio of the formulation #1 apparent oral clearance to the reference (unitless; exp = 0.610)") # Table 2: 'DRUG1 ON CL/F' = 0.610 (RSE 5.9%; bootstrap mean 0.610, 95% CI 0.545-0.695); Eq. 8. Estimated separately from the DRUG0 effect (different RSE) but numerically identical to 3 significant figures.

    # Between-subject variability. Methods Equation (1) is P_i = P_pop * EXP(eta_i)
    # with eta of 'mean zero and variance omega^2', i.e. exponential IIV. Table 2
    # prints the diagonal elements as percentages and the CL/F-Vc/F element as the
    # raw covariance 0.174, so the percentages are read as omega (the log-scale SD)
    # and the variances below are their squares. See the vignette Errata for the
    # arithmetic that supports this reading over an exact-lognormal-CV reading.
    etalka ~ 0.777924 # Table 2, row 'IIV - Ka' = 88.2% (RSE 24.0%); variance = 0.882^2
    etalvp ~ 0.274576 # Table 2, row 'IIV - Vp/F' = 52.4% (RSE 29.3%); variance = 0.524^2

    # CL/F and Vc/F share an OMEGA BLOCK(2). Diagonals: 'IIV - CL/F' = 49.1%
    # (0.491^2 = 0.241081) and 'IIV - Vc/F' = 40.8% (0.408^2 = 0.166464).
    # Off-diagonal: 'omega cov CL/F-Vc/F' = 0.174, which implies a correlation of
    # 0.174 / sqrt(0.241081 * 0.166464) = 0.869. The Discussion devotes a paragraph
    # to why CL/F and Vc/F correlate after oral dosing and states that the
    # correlation was estimated 'via the $OMEGA BLOCK syntax, as presented in this
    # study'.
    etalcl + etalvc ~ c(
      0.241081,
      0.174, 0.166464
    )

    # Table 2, row 'IIV - Q/F' = '0 FIXED' -- the paper held the Q/F
    # between-subject variance at zero, with no RSE and no bootstrap row.
    etalq ~ fixed(0)

    # Residual error. Methods Equation (4), the mixed model that the Results say
    # best described the data: Obs = Pred * (1 + eps1) + eps2.
    propSd <- 0.216
    label("Proportional residual error (fraction)") # Table 2: 'Prop-Error' = 0.216 (RSE 5.0%; bootstrap mean 0.200, 95% CI 0.189-0.232); Results: 'The proportional error was 21.6%'
    addSd <- 0.303
    label("Additive residual error (ng/mL)") # Table 2: 'Add-Error' = 0.303 (RSE 14.5%; bootstrap mean 0.300, 95% CI 0.215-0.382); Results: 'the additive error was 0.303 ng/mL'
  })

  model({
    # 1. Individual parameters. Results Equations (5)-(9), with the formulation
    #    multipliers carried as log ratios so they add on the log scale:
    #      ka   = 2.85 * 1.89^FORM_OLZ_F0
    #      Vc/F = 2390 * (WT/60.59)^0.579 * 0.216^FORM_OLZ_F0 * 0.207^FORM_OLZ_F1
    #      CL/F = 25.4  * 0.610^FORM_OLZ_F0 * 0.610^FORM_OLZ_F1
    #    FORM_OLZ_F0 and FORM_OLZ_F1 are mutually exclusive; both zero selects
    #    formulation #2, whose equations (7) and (9) carry no multiplier.
    ka <- exp(lka + e_form_olz_f0_ka * FORM_OLZ_F0 + etalka)
    cl <- exp(lcl + e_form_olz_f0_cl * FORM_OLZ_F0 + e_form_olz_f1_cl * FORM_OLZ_F1 + etalcl)
    vc <- exp(lvc + e_form_olz_f0_vc * FORM_OLZ_F0 + e_form_olz_f1_vc * FORM_OLZ_F1 + etalvc) *
      (WT / 60.59)^e_wt_vc
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    tlag <- exp(ltlag)

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODE system -- two-compartment disposition with first-order oral
    #    absorption (NONMEM ADVAN4). Bioavailability is not identifiable from
    #    oral-only data, so cl, vc, q and vp are all apparent (/F) quantities and
    #    no separate F term is applied; the paper names them CL/F, Vc/F, Q/F and
    #    Vp/F throughout.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Absorption lag time, NONMEM ALAG1 on the depot.
    alag(depot) <- tlag

    # 5. Observation. Doses are in mg and vc is in L, giving mg/L; the factor of
    #    1000 converts to the ng/mL in which the paper reports olanzapine
    #    concentrations (HPLC-MS/MS assay range 2-400 ng/mL).
    Cc <- 1000 * central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
