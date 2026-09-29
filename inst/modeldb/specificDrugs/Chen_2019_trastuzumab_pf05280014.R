Chen_2019_trastuzumab_pf05280014 <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order linear elimination from the central",
    "compartment for intravenous PF-05280014, a trastuzumab biosimilar (Trazimera), in patients with",
    "HER2-positive metastatic breast cancer treated with PF-05280014 plus paclitaxel (Chen 2019,",
    "NCT01989676). Baseline body weight enters clearance and central volume as power terms normalised",
    "to the 68.2 kg arm median. Inter-individual variability was estimated on CL, V1, V2 and Q with a",
    "diagonal omega matrix, and the residual error is additive on the natural-log scale, i.e.",
    "log-normal. Chen 2019 fitted the biosimilar and the EU-sourced reference product as two separate",
    "models on the two treatment arms; the companion reference-product model is",
    "Chen_2019_trastuzumab_reference."
  )
  reference <- paste(
    "Chen X, Li C, Ewesuedo R, Yin D.",
    "Population pharmacokinetics of PF-05280014 (a trastuzumab biosimilar) and reference trastuzumab",
    "(Herceptin) in patients with HER2-positive metastatic breast cancer.",
    "Cancer Chemother Pharmacol. 2019;84(1):83-92. doi:10.1007/s00280-019-03850-1.",
    "Correction: Cancer Chemother Pharmacol. 2019;84(3):667. doi:10.1007/s00280-019-03890-7",
    "(open-access licence change only; no model value was revised).",
    "Baseline demographics, including the arm median body weight used here as the covariate centering",
    "value, are from the companion trial report cited as reference 7 by Chen 2019:",
    "Pegram MD, Bondarenko I, Zorzetto MMC, et al. Br J Cancer. 2019;120(2):172-182.",
    "doi:10.1038/s41416-018-0340-2."
  )
  vignette <- "Chen_2019_trastuzumab"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed baseline value. Enters clearance and central volume as the normalised power model",
        "of Chen 2019 Eq. 1, TVP = Ppop * (COV/COVmedian)^theta, with theta = 0.637 on CL and 0.507 on",
        "V1 (Chen 2019 Table 2). Baseline body weight was the only covariate retained by the stepwise",
        "forward-inclusion (alpha = 0.05) / backward-elimination (alpha = 0.001) search.",
        "IMPORTANT: Chen 2019 never prints the median body weight that Eq. 1 normalises to, and the",
        "paper contains no baseline demographics table at all. The 68.2 kg used here is the median",
        "body weight of the PF-05280014 arm reported in Table 1 of the companion trial report for the",
        "same study (Pegram 2019, Br J Cancer 120:172-182, reference 7 of Chen 2019; mean 69.1 kg,",
        "SD 17.1 kg, range 29-147 kg, n = 352 ITT). The value is therefore companion-paper derived",
        "rather than printed in Chen 2019 itself; see the vignette Errata. It is corroborated by",
        "Chen 2019 Fig. 3, whose simulated cycle-1 day-1 peak median of 83.8 mg/L matches the",
        "4 mg/kg * 68.2 kg / 3.15 L loading-dose prediction."
      ),
      source_name = "BWT"
    )
  )

  # Issue #482: what each compartment holds, in what amount units, in what
  # biological matrix. Chen 2019 measured serum drug concentrations by a
  # validated ELISA (Methods, 'Pharmacokinetic evaluations'), so the matrix is
  # serum rather than plasma. There is no depot: every dose is an intravenous
  # infusion into the central compartment. The model is solved with linCmt(),
  # so `central` is the only named state; the peripheral compartment of the
  # two-compartment system (V2 = 5.55 L, reached through Q = 0.0194 L/h) is
  # internal to the analytic kernel and is not separately addressable.
  compartmentData <- list(
    central = list(
      analyte = "PF-05280014 (trastuzumab biosimilar)",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    )
  )

  # Screened by the stepwise covariate search (Chen 2019 Table 1) on both CL
  # and V1, and retained by neither. Chen 2019 Discussion attributes the
  # difference from Bruno 2005 -- which did retain the number of metastatic
  # sites and baseline HER2 ECD on CL -- to the more stringent backward
  # elimination criterion used here (p < 0.001 versus p < 0.005). None of
  # these is referenced in model().
  covariatesDataExcluded <- list(
    HER2_ECD = list(
      description = "Baseline circulating (shed) HER2 extracellular domain concentration",
      units = "ng/mL",
      type = "continuous",
      notes = paste(
        "Source column HER2STAT. Screened on CL and V1 with the Eq. 1 power model and not retained.",
        "Measured at pre-dose of cycles 1, 3, 5 and 8 and at end of treatment with a HER2 ELISA",
        "development kit (Nuclea Diagnostics). Chen 2019 reports no point estimate for its effect."
      ),
      source_name = "HER2STAT"
    ),
    RACE_JAPANESE = list(
      description = "Japanese heritage indicator",
      units = "(binary)",
      type = "categorical",
      notes = paste(
        "Source column JAPA (Japanese versus non-Japanese). Screened with the Eq. 2 fractional-change",
        "model and not retained. Only 18 of the 349 PF-05280014 patients were Japanese; Chen 2019",
        "Fig. 4a instead overlays their observed concentrations on the all-patient VPC and concludes",
        "the PK is consistent."
      ),
      source_name = "JAPA"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "categorical",
      notes = paste(
        "Source column RACE_STAT (Asian versus non-Asian). Screened with the Eq. 2 fractional-change",
        "model and not retained. 104 of 352 (29.5%) of the ITT arm were Asian (Pegram 2019 Table 1)."
      ),
      source_name = "RACE_STAT"
    ),
    ADA_POS = list(
      description = "Baseline antidrug antibody positive indicator",
      units = "(binary)",
      type = "categorical",
      notes = paste(
        "Source column ADA_BL. Screened with the Eq. 2 fractional-change model and not retained.",
        "Anti-PF-05280014 antibodies were assayed by a validated electrochemiluminescent immunoassay",
        "with a screen / confirm / titre tiered approach."
      ),
      source_name = "ADA_BL"
    ),
    WHO_PS = list(
      description = "Baseline ECOG performance-status score",
      units = "(integer score)",
      type = "continuous",
      notes = paste(
        "Source column ECOG. Screened with the Eq. 2 fractional-change model and not retained.",
        "Pegram 2019 Table 1 reports 52.8% at 0, 42.6% at 1 and 4.5% at 2 for this arm."
      ),
      source_name = "ECOG"
    ),
    MET_GE4 = list(
      description = "Indicator of baseline number of metastatic sites >= 4",
      units = "(binary)",
      type = "categorical",
      notes = paste(
        "Source column N_META. Screened with the Eq. 2 fractional-change model and not retained.",
        "Chen 2019 does not restate the dichotomisation; the >= 4 cut is the definition used by",
        "Bruno 2005 (reference 8 of Chen 2019), the trastuzumab popPK analysis whose covariate list",
        "Chen 2019 adopted and which did retain this effect on CL. See modellib('Bruno_2005_trastuzumab')."
      ),
      source_name = "N_META"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 349L,
    n_studies = 1L,
    age_median = "55.0 years",
    age_range = "19-80 years",
    weight_median = "68.2 kg",
    weight_range = "29-147 kg",
    weight_mean = "69.1 kg (SD 17.1)",
    sex_female_pct = 100,
    race_ethnicity = c(White = 65.9, Asian = 29.5, Black = 1.4, Other = 3.1),
    disease_state = "HER2-positive metastatic breast cancer, first-line treatment setting",
    dose_range = paste(
      "4 mg/kg IV loading dose infused over 90 min on cycle 1 day 1, then 2 mg/kg infused over",
      "30-90 min on days 8, 15 and 22 of cycle 1 and on days 1, 8, 15 and 22 of each subsequent",
      "28-day cycle, until at least week 33. From no earlier than week 33 the regimen could be",
      "changed at investigator discretion to 6 mg/kg infused over 30-90 min every 3 weeks."
    ),
    regions = "Multinational (global phase III trial NCT01989676)",
    co_medication = "Paclitaxel 80 mg/m2 IV over 60 min on days 1, 8 and 15 of each 28-day cycle for at least six cycles",
    ecog_status = "ECOG 0 in 52.8%, 1 in 42.6%, 2 in 4.5% (Pegram 2019 Table 1)",
    reference_subject = "Baseline body weight 68.2 kg (the arm median; see covariateData$WT)",
    notes = paste(
      "PF-05280014 arm of the population PK analysis (Chen 2019 Results). 349 of the 352 randomised",
      "PF-05280014 patients received study treatment. All serum concentration-time data up to and",
      "including cycle 17 day 1 as of the 24 August 2016 primary completion date were used, excluding",
      "the end-of-treatment visit and unplanned records. Sampling was peak-and-trough: pre-dose within",
      "4 h of the start of infusion on day 1 of cycles 1, 3, 4, 5, 7 and 8 and on day 8 of cycles 1 and",
      "5, end-of-infusion samples 1 h after the end of infusion on day 1 of cycles 1 and 5, and",
      "pre-dose on day 1 every 3 cycles thereafter. Concentrations were measured by validated ELISA",
      "over a 0.5-100 mg/L calibration range; 43 of 7098 post-dose observations were below the LLOQ",
      "and were excluded (M1 method). Concentrations were log-transformed before fitting. Estimation",
      "was FOCE in NONMEM 7.2 with PsN 4.2.0 for stepwise covariate modelling, VPC and bootstrap; 916",
      "of 1000 bootstrap replicates minimised successfully. Baseline demographics are from Table 1 of",
      "the companion trial report (Pegram 2019) and are reported there for the 352-patient ITT arm;",
      "Chen 2019 itself contains no demographics table. Sex is not tabulated by Pegram 2019, whose",
      "eligibility criteria describe a metastatic breast cancer population; 100% female is an",
      "assumption recorded in the vignette Errata rather than a reported figure."
    )
  )

  ini({
    # ---- Structural PK: Chen 2019 Table 2, 'PF-05280014 / NONMEM results'
    # column. Values are for the typical patient at the arm median baseline
    # body weight. Doses are in mg and volumes in L, so Cc is in mg/L, which is
    # the unit Chen 2019 uses on every concentration axis and equals ug/mL.
    lcl <- log(0.0104)
    label("Linear clearance at 68.2 kg baseline body weight (L/h)") # Table 2, row 'CL (L/h)' = 0.0104 [95% CI 0.0098-0.0110]; bootstrap median 0.0103
    lvc <- log(3.15)
    label("Central volume of distribution V1 at 68.2 kg baseline body weight (L)") # Table 2, row 'V 1 (L)' = 3.15 [95% CI 2.99-3.31]; bootstrap median 3.15
    lq <- log(0.0194)
    label("Intercompartmental clearance Q (L/h)") # Table 2, row 'Q (L/h)' = 0.0194 [95% CI 0.0160-0.0228]; bootstrap median 0.0192; no covariate retained
    lvp <- log(5.55)
    label("Peripheral volume of distribution V2 (L)") # Table 2, row 'V 2 (L)' = 5.55 [95% CI 5.24-5.86]; bootstrap median 5.59; no covariate retained

    # ---- Covariate effects. Chen 2019 Eq. 1 is the normalised power model
    # TVP_j = Ppop * (COV_j / COVmedian)^theta, applied to CL and V1 only.
    # COVmedian is not printed in Chen 2019; see covariateData$WT for the
    # 68.2 kg sourcing.
    e_wt_cl <- 0.637
    label("Power exponent of baseline body weight on linear clearance (unitless; reference 68.2 kg)") # Table 2, row 'BWT effect on CL' = 0.637 [95% CI 0.450-0.824]; bootstrap median 0.638
    e_wt_vc <- 0.507
    label("Power exponent of baseline body weight on central volume V1 (unitless; reference 68.2 kg)") # Table 2, row 'BWT effect on V 1' = 0.507 [95% CI 0.316-0.698]; bootstrap median 0.513

    # ---- Inter-individual variability. Chen 2019 Methods, 'Structural PK
    # model and variability models': IIV on CL, Q, V1 and V2 was log-normal
    # (exponential model) and the omega matrix was DIAGONAL, because no
    # significant correlation between CL and V1 was found during model
    # development. Table 2 prints each row as 'X omega 2 (%CV)', so the leading
    # number is the VARIANCE on the log scale and the parenthesised figure is
    # its square root as a percentage; sqrt(0.0934) = 0.306 -> 31%, and the
    # same identity holds for all four rows. The variances are therefore used
    # verbatim, with no CV-to-variance conversion.
    etalcl ~ 0.0934 # Table 2, row 'CL omega 2 (%CV)' = 0.0934 (31) [95% CI 0.072-0.115]; bootstrap median 0.091
    etalvc ~ 0.0405 # Table 2, row 'V 1 omega 2 (%CV)' = 0.0405 (20) [95% CI 0.011-0.070]; bootstrap median 0.040
    etalq ~ 0.504 # Table 2, row 'Q omega 2 (%CV)' = 0.504 (71) [95% CI 0.332-0.676]; bootstrap median 0.491
    etalvp ~ 1.06 # Table 2, row 'V 2 omega 2 (%CV)' = 1.06 (103) [95% CI 0.770-1.350]; bootstrap median 1.057

    # ---- Residual error. Chen 2019 Methods: 'The residual error was described
    # using an additive error model, after log-transformation of the PK data',
    # which is a log-normal residual, encoded here as lnorm(expSd).
    #
    # SCALE. Table 2 prints 'Res Add Err' = 0.272 without marking it as a
    # variance or an SD. It is the NONMEM $SIGMA VARIANCE, so the SD entered
    # here is sqrt(0.272) = 0.5215. Three independent lines of evidence:
    # (1) Chen 2019 Fig. 1a, the log-observed versus log-individual-predicted
    #     panel, is a direct picture of this residual. Digitising its blue
    #     scatter against the printed axes gives a 2.5th-to-97.5th percentile
    #     residual span of about +/- 1.2 natural-log units. An SD of 0.5215
    #     predicts +/- 1.02; an SD of 0.272 predicts only +/- 0.53, less than
    #     half what the figure shows.
    # (2) Every other variance component in Table 2 is a NONMEM variance, and
    #     the table column is headed 'NONMEM results Estimate'. NONMEM reports
    #     $SIGMA on the variance scale.
    # (3) The relative standard error implied by the printed 95% CI (the table
    #     footnote states the CI is estimate +/- 1.96 x SE) is 5.4%. For a
    #     variance that implies about 680 effective observations, consistent
    #     with 7098 correlated peak-and-trough samples in 349 patients; read as
    #     an SD it would imply only about 170, fewer than one independent
    #     observation per patient.
    # See the vignette Errata for the full argument and the digitisation.
    expSd <- 0.5215
    label("Log-normal residual error SD on the natural-log scale (unitless)") # Table 2, row 'Res Add Err' = 0.272 [95% CI 0.243-0.301] read as a variance; sqrt(0.272) = 0.5215
  })

  model({
    # Individual PK parameters. Baseline body weight enters CL and V1 only, as
    # the Eq. 1 normalised power term centred on the 68.2 kg arm median; Q and
    # V2 carry no covariate (Chen 2019 Table 2 lists no BWT effect for either).
    # Inter-individual variability is exponential on all four parameters.
    cl <- exp(lcl + etalcl) * (WT / 68.2)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 68.2)^e_wt_vc
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    # Two-compartment disposition with first-order elimination from the central
    # compartment and zero-order intravenous input (Chen 2019 Results,
    # 'Determination of the structural PK model'). Every dose is an IV infusion
    # directly into the central compartment, so there is no depot and no
    # absorption process; encode the infusion with rate or dur on the dose
    # record. The four-parameter cl / vc / q / vp set is exactly rxode2's
    # two-compartment analytic kernel, so linCmt() is used explicitly rather
    # than writing d/dt() bodies that rxode2 would silently replace with the
    # same closed form.
    Cc <- linCmt()

    Cc ~ lnorm(expSd)
  })
}
