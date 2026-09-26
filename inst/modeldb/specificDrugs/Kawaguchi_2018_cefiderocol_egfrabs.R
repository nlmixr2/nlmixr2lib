Kawaguchi_2018_cefiderocol_egfrabs <- function() {
  description <- "Three-compartment population PK model for intravenous cefiderocol in healthy subjects, subjects spanning normal renal function to end-stage renal disease, and patients with complicated urinary tract infection or acute uncomplicated pyelonephritis, with absolute (non-BSA-normalized) estimated glomerular filtration rate on CL, body weight on V1 and V2, and an infected-versus-uninfected disease-status factor on CL and V1"
  reference <- paste(
    "Kawaguchi N, Katsube T, Echols R, Wajima T.",
    "Population pharmacokinetic analysis of cefiderocol, a parenteral",
    "siderophore cephalosporin, in healthy subjects, subjects with various",
    "degrees of renal function, and patients with complicated urinary tract",
    "infection or acute uncomplicated pyelonephritis.",
    "Antimicrob Agents Chemother. 2018;62(1):e01391-17.",
    "doi:10.1128/AAC.01391-17",
    sep = " "
  )
  vignette <- "Kawaguchi_2018_cefiderocol"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Absolute estimated glomerular filtration rate (the paper's",
        "'eGFRabs'): the BSA-adjusted eGFR de-normalized back to whole-body",
        "mL/min by multiplying by the individual's body surface area and",
        "dividing by 1.73 m^2. RAW mL/min and NOT BSA-normalized."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The underlying BSA-adjusted eGFR was calculated by the Modification",
        "of Diet in Renal Disease (MDRD) equation, or by the equation of",
        "Matsuo et al. for Japanese subjects (Introduction), then converted",
        "to this absolute form as eGFRabs = eGFRadj * BSA / 1.73. The",
        "two-equation split matters: the phase 1 ascending-dose study was",
        "run in Japan and the renal-impairment study in the United States.",
        "Reference 83.0 mL/min is the OVERALL cohort median from Table 1",
        "('eGFRabs (ml/min)', Overall, Median (range) 83.0 (5-148));",
        "subgroup medians are 99.0 without infection and 78.0 with",
        "infection, so the centring constant is the pooled value.",
        "Kawaguchi 2018 fitted THREE parallel final models that differ only",
        "in which renal-function marker enters CL. This model is the",
        "eGFRabs arm; the siblings are Kawaguchi_2018_cefiderocol_clcr",
        "(Cockcroft-Gault creatinine clearance, raw mL/min, reference 90.0)",
        "and Kawaguchi_2018_cefiderocol_egfradj (the same MDRD/Matsuo",
        "estimate left BSA-normalized, mL/min/1.73 m^2, reference 77.0).",
        "The three are NOT interchangeable inputs to one another's models:",
        "each carries its own centring constant and exponent, and the",
        "eGFRadj arm additionally retains a different covariate set. In",
        "particular, eGFRabs and eGFRadj are the SAME underlying estimate",
        "differing only by the BSA factor, so supplying one where the other",
        "is expected silently rescales CL.",
        "This model was not the one the authors used for their post hoc",
        "analysis -- the CLCR model had the lowest objective function value",
        "(9,363.552 against 9,377.486 here) -- but the paper states that the",
        "difference in predictive performance 'would not be clinically",
        "significant' and that any of the three markers can be used for dose",
        "adjustment.",
        "Time-fixed: baseline subject characteristics were used (Materials",
        "and Methods, Data). Observed range 5-148 mL/min."
      ),
      source_name = "eGFRabs"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference 74.1 kg is the OVERALL cohort median from Table 1",
        "('Body wt (kg)', Overall, Median (range) 74.1 (45.1-138.0));",
        "subgroup medians are 68.4 without infection and 76.4 with",
        "infection.",
        "Enters V1 and V2 only in this model. Body weight was tested on CL",
        "as well but was NOT retained here, because -- as the Discussion",
        "explains -- an absolute (non-BSA-normalized) filtration rate",
        "already carries body scale, so a separate weight term on CL is",
        "redundant. The sibling Kawaguchi_2018_cefiderocol_egfradj model",
        "DOES carry a weight term on CL, precisely because its marker has",
        "had body size divided out. Do not add a weight effect on CL to",
        "this model.",
        "Time-fixed (baseline). The authors judged the weight effect on V1",
        "not clinically significant: post hoc median V1 relative to the",
        "typical infected-patient value ran 0.85, 0.87, 0.96 and 1.27 across",
        "the <55, 55-<70, 70-<90 and >=90 kg groups (Discussion, Fig. 4)."
      ),
      source_name = "body weight"
    ),
    DIS_INFECT_ACTIVE = list(
      description = paste(
        "Active clinical infection indicator: 1 = patient enrolled in the",
        "phase 2 study with complicated urinary tract infection (cUTI) or",
        "acute uncomplicated pyelonephritis (AUP) caused by a Gram-negative",
        "pathogen; 0 = subject without infection."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (subject without infection: healthy volunteer or renal-impairment-study subject)",
      notes = paste(
        "The paper calls this column 'disease status' and states the coding",
        "explicitly in the Table 2 footnote: 'disease status = 0 for",
        "subjects without infection and disease status = 1 for patients with",
        "infection'. Same orientation as the canonical, so no value",
        "transformation is needed.",
        "Time-fixed per subject in this analysis, unlike the time-varying",
        "per-record use in Kloos_2021_pegasparaginase.R: cohort membership",
        "is study membership. The 91 uninfected subjects came from two",
        "phase 1 studies and the 238 infected patients from one multinational",
        "phase 2 cUTI study (Table S1).",
        "Enters BOTH CL (x1.15) and V1 (x1.36) in this model. Note the",
        "clearance factor is 1.15 here against 1.26 in the CLCR sibling: the",
        "same biological contrast is partitioned differently depending on",
        "which renal marker absorbs it, which is a good illustration of why",
        "the three models must not be mixed and matched."
      ),
      source_name = "disease status"
    )
  )

  # Covariates Kawaguchi 2018 screened but did not retain in this final model.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on both CL and V1 (Materials and Methods, Population",
        "pharmacokinetic analyses) and not retained. Overall median 59.0",
        "years, range 18-93 (Table 1); the uninfected arm was much younger",
        "(median 36.0) than the infected arm (median 65.0), so age is",
        "strongly confounded with disease status in this cohort."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Screened on both CL and V1 and not retained. 44.4% female overall",
        "(Table 1). Sex does not enter the MDRD/Matsuo estimate the way it",
        "enters Cockcroft-Gault in the CLCR sibling -- MDRD carries its own",
        "sex factor, so the confounding is structurally similar."
      )
    ),
    ALB = list(
      description = "Serum albumin concentration",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on both CL and V1 and not retained. Overall median 4.2",
        "g/dL, range 2.5-5.3 (Table 1). The later, larger analysis by the",
        "same group (Katsube/Wajima, AAC 2021, doi:10.1128/AAC.01437-20)",
        "DID retain albumin on V1 with an exponent of -0.617 in a cohort",
        "including pneumonia and bloodstream-infection patients whose",
        "albumin ran much lower (mean 2.8 g/dL); this cohort simply did not",
        "span enough hypoalbuminaemia to detect it."
      )
    ),
    AST = list(
      description = "Aspartate aminotransferase concentration",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL and not retained. Overall median 18.0 U/L, range 6-101 (Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase concentration",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL and not retained. Overall median 16.0 U/L, range 4-111 (Table 1)."
    ),
    BILI = list(
      description = "Total bilirubin concentration",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL and not retained. Overall median 0.57 mg/dL, range 0.19-2.88 (Table 1)."
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-White)",
      notes = paste(
        "Race was screened on both CL and V1 and not retained. 76.9% White",
        "overall (Table 1), but the split is almost perfectly confounded",
        "with disease status -- 25.3% White in the uninfected arm against",
        "96.6% in the infected arm -- so this cohort cannot separate a race",
        "effect from the disease-status effect that was retained."
      )
    )
  )

  compartmentData <- list(
    central = list(analyte = "cefiderocol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefiderocol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "cefiderocol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 329L,
    n_studies = 3L,
    n_observations = 2571L,
    age_range = "18-93 years (overall median 59.0; 20-74, median 36.0 without infection; 18-93, median 65.0 with infection)",
    weight_range = "45.1-138.0 kg (overall median 74.1)",
    sex_female_pct = 44.4,
    race_ethnicity = c(White = 76.9, Asian = 17.0, Black = 5.2, `Native American or Alaska Native` = 0.3, Other = 0.6),
    disease_state = paste(
      "Pooled: 91 subjects without infection (healthy volunteers plus",
      "subjects spanning normal renal function, mild, moderate and severe",
      "renal impairment, and end-stage renal disease requiring hemodialysis)",
      "and 238 patients with complicated urinary tract infection (n = 175) or",
      "acute uncomplicated pyelonephritis (n = 63) caused by Gram-negative",
      "pathogens."
    ),
    renal_function = paste(
      "Absolute eGFR 5-148 mL/min overall (median 83.0); BSA-adjusted eGFR",
      "4-146 mL/min/1.73 m^2 (median 77.0); Cockcroft-Gault CLCR 7-186",
      "mL/min (median 90.0). Renal function was an enrolment axis of one of",
      "the three contributing studies, so the cohort deliberately spans",
      "end-stage renal disease through augmented renal function."
    ),
    dose_range = paste(
      "Phase 1 part 1: single cefiderocol 0.1, 0.25, 0.5, 1 and 2 g doses",
      "infused over 1 h. Phase 1 part 2: 1 or 2 g every 8 h over 1 h for 10",
      "days. Phase 1 renal-impairment study: single 1 g dose over 1 h.",
      "Phase 2 cUTI study: 2 g every 8 h over 1 h for 7 or 14 days, with the",
      "dose reduced for renal function and body size per Table S2 (regimens",
      "actually given ran 0.75 g q6h through 2 g q8h). All doses",
      "intravenous. Cefiderocol PK are linear over 100-2,000 mg."
    ),
    regions = "Japan (phase 1 ascending-dose study), United States (phase 1 renal-impairment study), multinational (phase 2 cUTI study).",
    bioanalytic_methods = paste(
      "Plasma treated 1:1 by volume with 0.2 mol/L ammonium acetate buffer,",
      "pH 5, then assayed by validated LC-MS/MS. Linear 0.1-100 ug/mL;",
      "lower limit of quantification 0.1 ug/mL; precision 1.2-6.2% and",
      "accuracy -5.3% to 2.1%."
    ),
    protein_binding = "In vitro plasma protein binding 57.8%; the paper used a fixed unbound fraction of 0.422 for its fT>MIC calculations because the unbound fraction was not measured in the phase 2 study.",
    notes = paste(
      "Baseline subject characteristics are in Table 1 and the study designs",
      "in Table S1. 2,571 plasma concentrations were analysed after excluding",
      "264 below-limit-of-quantification samples, 156 samples stored at -20 C",
      "for more than 7 days before buffer stabilisation, 8 anomalous",
      "concentrations, and 3 samples with unidentified sampling times."
    )
  )

  ini({
    # Structural parameters. Table 2, column 'Final model with eGFRabs'. CL,
    # V1 and V2 are the leading constants of the Table 2 footnote b
    # equations, i.e. the typical values for an UNINFECTED subject at the
    # reference eGFRabs of 83.0 mL/min and reference body weight of 74.1 kg.
    lcl <- log(4.56); label("Typical total clearance in an uninfected subject at absolute eGFR 83.0 mL/min (L/h)") # Table 2, 'CL', final model with eGFRabs = 4.56 (RSE 1.8%); footnote b leading constant
    lvc <- log(7.92); label("Typical central volume of distribution V1 in an uninfected subject at 74.1 kg (L)") # Table 2, 'V 1', final model with eGFRabs = 7.92 (RSE 3.2%); footnote b leading constant
    lq <- log(5.78); label("Typical intercompartmental clearance Q2 between central and peripheral1 (L/h)") # Table 2, 'Q 2', final model with eGFRabs = 5.78 (RSE 5.6%)
    lvp <- log(5.41); label("Typical first peripheral volume of distribution V2 at 74.1 kg (L)") # Table 2, 'V 2', final model with eGFRabs = 5.41 (RSE 3.4%); footnote b leading constant
    lq2 <- log(0.109); label("Typical intercompartmental clearance Q3 between central and peripheral2 (L/h)") # Table 2, 'Q 3', final model with eGFRabs = 0.109 (RSE 17.2%)
    lvp2 <- log(0.735); label("Typical second peripheral volume of distribution V3 (L)") # Table 2, 'V 3', final model with eGFRabs = 0.735 (RSE 8.6%)

    # Covariate effects. Table 2 footnote b prints the full equations:
    #   CL = 4.56 * (eGFRabs/83.0)^0.621 * 1.15^disease_status
    #   V1 = 7.92 * (body weight/74.1)^0.789 * 1.36^disease_status
    #   V2 = 5.41 * (body weight/74.1)^0.673
    # Each exponent / factor also appears as its own Table 2 row, so every
    # value below is printed twice in the paper and the two agree.
    e_crcl_cl <- 0.621; label("Power exponent on (CRCL / 83.0) for CL (unitless)") # Table 2, 'Effect of renal function marker on CL', eGFRabs column = 0.621 (RSE 3.5%); footnote b
    e_wt_vc <- 0.789; label("Power exponent on (WT / 74.1) for V1 (unitless)") # Table 2, 'Effect of body wt on V 1', eGFRabs column = 0.789 (RSE 12.8%); footnote b
    e_wt_vp <- 0.673; label("Power exponent on (WT / 74.1) for V2 (unitless)") # Table 2, 'Effect of body wt on V 2', eGFRabs column = 0.673 (RSE 10.8%); footnote b
    e_infect_cl <- 1.15; label("Multiplicative effect of active infection on CL (unitless)") # Table 2, 'Effect of disease status on CL', eGFRabs column = 1.15 (RSE 3.2%); footnote b exponentiated by disease status
    e_infect_vc <- 1.36; label("Multiplicative effect of active infection on V1 (unitless)") # Table 2, 'Effect of disease status on V 1', eGFRabs column = 1.36 (RSE 5.1%); footnote b exponentiated by disease status

    # Interindividual variability. Table 2 reports IIV as a percent CV. For
    # THIS research group the percent CV is the omega standard deviation
    # times 100 -- NOT the log-normal sqrt(exp(omega^2)-1) form that the
    # skill's verification checklist warns about as the usual default. The
    # companion analysis by the same authors (Katsube/Wajima, AAC 2021,
    # doi:10.1128/AAC.01437-20) settles this by printing CV%, the omega
    # covariances AND the implied correlations: all three covariance pairs
    # agree with omega = CV/100 to within 0.1% and reject the log-normal
    # reading by 6-10%. See Kawaguchi_2018_cefiderocol_clcr.R for the
    # arithmetic. The variances below are therefore (CV/100)^2.
    #
    # Kawaguchi 2018 reports no covariances between the etas, so the omega
    # matrix is encoded as diagonal.
    etalcl ~ 0.106276 # Table 2, '% CV for IIV for CL', eGFRabs column = 32.6% (shrinkage 2.9%, RSE 14.4%); 0.326^2
    etalvc ~ 0.214369 # Table 2, '% CV for IIV for V1', eGFRabs column = 46.3% (shrinkage 11.0%, RSE 27.7%); 0.463^2
    etalvp ~ 0.146689 # Table 2, '% CV for IIV for V2', eGFRabs column = 38.3% (shrinkage 34.2%, RSE 35.6%); 0.383^2

    # Residual error. The proportional model was selected over the combined
    # additive-plus-proportional model because the latter gave PK parameter
    # RSEs of 4.1% to 1,721%, which the authors read as non-robust (Results).
    propSd <- 0.151; label("Proportional residual error (fraction)") # Table 2, '% CV for proportional residual error', eGFRabs column = 15.1% (shrinkage 14.2%, RSE 12.6%)
  })

  model({
    # Individual PK parameters. Covariate forms are Table 2 footnote b
    # verbatim. As in the CLCR sibling, the renal marker is an absolute
    # (whole-body mL/min) rate that already carries body scale, so CL takes
    # no separate body-weight term. Q2, Q3 and V3 carry neither covariates
    # nor IIV.
    cl <- exp(lcl + etalcl) * (CRCL / 83.0)^e_crcl_cl * e_infect_cl^DIS_INFECT_ACTIVE
    vc <- exp(lvc + etalvc) * (WT / 74.1)^e_wt_vc * e_infect_vc^DIS_INFECT_ACTIVE
    q <- exp(lq)
    vp <- exp(lvp + etalvp) * (WT / 74.1)^e_wt_vp
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)

    # Three-compartment intravenous PK with first-order elimination from the
    # central compartment (NONMEM ADVAN11 / TRANS4 mass balance). Doses land
    # in `central`; the 1 h infusion is encoded on the dose record as a rate
    # or duration.
    d/dt(central) <- q / vp * peripheral1 + q2 / vp2 * peripheral2 -
      (cl + q + q2) / vc * central
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1
    d/dt(peripheral2) <- q2 / vc * central - q2 / vp2 * peripheral2

    # Total plasma cefiderocol concentration. Dose in mg and vc in L give
    # mg/L, which is ug/mL. The paper's fT>MIC analysis multiplies this by a
    # fixed unbound fraction of 0.422 to obtain the free concentration; that
    # scaling is a post-processing step and is deliberately not baked into
    # the observation here.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
