Chen_2019_trastuzumab_reference <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order linear elimination from the central",
    "compartment for intravenous EU-sourced reference trastuzumab (Herceptin), fitted as the active",
    "comparator arm of a phase III biosimilarity study in patients with HER2-positive metastatic breast",
    "cancer treated with trastuzumab plus paclitaxel (Chen 2019, NCT01989676). Baseline body weight",
    "enters clearance and central volume as power terms normalised to the 66.0 kg arm median.",
    "Inter-individual variability was estimated on CL, V1, V2 and Q with a diagonal omega matrix, and",
    "the residual error is additive on the natural-log scale, i.e. log-normal. Chen 2019 fitted the",
    "reference product and the PF-05280014 biosimilar as two separate models on the two treatment",
    "arms; the companion biosimilar model is Chen_2019_trastuzumab_pf05280014."
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
        "of Chen 2019 Eq. 1, TVP = Ppop * (COV/COVmedian)^theta, with theta = 0.673 on CL and 0.512 on",
        "V1 (Chen 2019 Table 2). Baseline body weight was the only covariate retained by the stepwise",
        "forward-inclusion (alpha = 0.05) / backward-elimination (alpha = 0.001) search.",
        "IMPORTANT: Chen 2019 never prints the median body weight that Eq. 1 normalises to, and the",
        "paper contains no baseline demographics table at all. The 66.0 kg used here is the median",
        "body weight of the trastuzumab-EU arm reported in Table 1 of the companion trial report for",
        "the same study (Pegram 2019, Br J Cancer 120:172-182, reference 7 of Chen 2019; mean 68.1 kg,",
        "SD 16.1 kg, range 36-139 kg, n = 355 ITT). The value is therefore companion-paper derived",
        "rather than printed in Chen 2019 itself; see the vignette Errata. Because Chen 2019 fitted the",
        "two arms as separate models, each arm is normalised to its own median, which is why this model",
        "uses 66.0 kg where the companion biosimilar model uses 68.2 kg. It is corroborated by",
        "Chen 2019 Fig. 3, whose simulated cycle-1 day-1 peak median of 83.8 mg/L matches the",
        "4 mg/kg * 66.0 kg / 3.10 L loading-dose prediction."
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
  # two-compartment system (V2 = 5.66 L, reached through Q = 0.0186 L/h) is
  # internal to the analytic kernel and is not separately addressable.
  compartmentData <- list(
    central = list(
      analyte = "trastuzumab (EU-sourced reference product, Herceptin)",
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
        "model and not retained. Only 14 of the 353 trastuzumab-EU patients were Japanese; Chen 2019",
        "Fig. 4b instead overlays their observed concentrations on the all-patient VPC and concludes",
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
        "model and not retained. 84 of 355 (23.7%) of the ITT arm were Asian (Pegram 2019 Table 1)."
      ),
      source_name = "RACE_STAT"
    ),
    ADA_POS = list(
      description = "Baseline antidrug antibody positive indicator",
      units = "(binary)",
      type = "categorical",
      notes = paste(
        "Source column ADA_BL. Screened with the Eq. 2 fractional-change model and not retained.",
        "Anti-trastuzumab-EU antibodies were assayed by a validated electrochemiluminescent immunoassay",
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
        "Pegram 2019 Table 1 reports 54.6% at 0, 41.1% at 1 and 4.2% at 2 for this arm."
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
    n_subjects = 353L,
    n_studies = 1L,
    age_median = "54.0 years",
    age_range = "25-85 years",
    weight_median = "66.0 kg",
    weight_range = "36-139 kg",
    weight_mean = "68.1 kg (SD 16.1)",
    sex_female_pct = 100,
    race_ethnicity = c(White = 68.7, Asian = 23.7, Black = 2.3, Other = 5.4),
    disease_state = "HER2-positive metastatic breast cancer, first-line treatment setting",
    dose_range = paste(
      "4 mg/kg IV loading dose infused over 90 min on cycle 1 day 1, then 2 mg/kg infused over",
      "30-90 min on days 8, 15 and 22 of cycle 1 and on days 1, 8, 15 and 22 of each subsequent",
      "28-day cycle, until at least week 33. From no earlier than week 33 the regimen could be",
      "changed at investigator discretion to 6 mg/kg infused over 30-90 min every 3 weeks."
    ),
    regions = "Multinational (global phase III trial NCT01989676)",
    co_medication = "Paclitaxel 80 mg/m2 IV over 60 min on days 1, 8 and 15 of each 28-day cycle for at least six cycles",
    ecog_status = "ECOG 0 in 54.6%, 1 in 41.1%, 2 in 4.2% (Pegram 2019 Table 1)",
    reference_subject = "Baseline body weight 66.0 kg (the arm median; see covariateData$WT)",
    notes = paste(
      "Trastuzumab-EU arm of the population PK analysis (Chen 2019 Results). 353 of the 355 randomised",
      "trastuzumab-EU patients received study treatment. All serum concentration-time data up to and",
      "including cycle 17 day 1 as of the 24 August 2016 primary completion date were used, excluding",
      "the end-of-treatment visit and unplanned records. Sampling was peak-and-trough: pre-dose within",
      "4 h of the start of infusion on day 1 of cycles 1, 3, 4, 5, 7 and 8 and on day 8 of cycles 1 and",
      "5, end-of-infusion samples 1 h after the end of infusion on day 1 of cycles 1 and 5, and",
      "pre-dose on day 1 every 3 cycles thereafter. Concentrations were measured by validated ELISA",
      "over a 0.5-100 mg/L calibration range; 43 of 7098 post-dose observations across both arms were",
      "below the LLOQ and were excluded (M1 method). Concentrations were log-transformed before",
      "fitting. Estimation was FOCE in NONMEM 7.2 with PsN 4.2.0 for stepwise covariate modelling, VPC",
      "and bootstrap; 932 of 1000 bootstrap replicates minimised successfully. Baseline demographics",
      "are from Table 1 of the companion trial report (Pegram 2019) and are reported there for the",
      "355-patient ITT arm; Chen 2019 itself contains no demographics table. Sex is not tabulated by",
      "Pegram 2019, whose eligibility criteria describe a metastatic breast cancer population; 100%",
      "female is an assumption recorded in the vignette Errata rather than a reported figure."
    )
  )

  ini({
    # ---- Structural PK: Chen 2019 Table 2, 'Trastuzumab-EU / NONMEM results'
    # column. Values are for the typical patient at the arm median baseline
    # body weight. Doses are in mg and volumes in L, so Cc is in mg/L, which is
    # the unit Chen 2019 uses on every concentration axis and equals ug/mL.
    lcl <- log(0.00948)
    label("Linear clearance at 66.0 kg baseline body weight (L/h)") # Table 2, row 'CL (L/h)' = 0.00948 [95% CI 0.009-0.010]; bootstrap median 0.00948
    lvc <- log(3.10)
    label("Central volume of distribution V1 at 66.0 kg baseline body weight (L)") # Table 2, row 'V 1 (L)' = 3.10 [95% CI 2.91-3.29]; bootstrap median 3.10
    lq <- log(0.0186)
    label("Intercompartmental clearance Q (L/h)") # Table 2, row 'Q (L/h)' = 0.0186 [95% CI 0.017-0.020]; bootstrap median 0.0186; no covariate retained
    lvp <- log(5.66)
    label("Peripheral volume of distribution V2 (L)") # Table 2, row 'V 2 (L)' = 5.66 [95% CI 5.12-6.20]; bootstrap median 5.65; no covariate retained

    # ---- Covariate effects. Chen 2019 Eq. 1 is the normalised power model
    # TVP_j = Ppop * (COV_j / COVmedian)^theta, applied to CL and V1 only.
    # COVmedian is not printed in Chen 2019; see covariateData$WT for the
    # 66.0 kg sourcing.
    e_wt_cl <- 0.673
    label("Power exponent of baseline body weight on linear clearance (unitless; reference 66.0 kg)") # Table 2, row 'BWT effect on CL' = 0.673 [95% CI 0.430-0.916]; bootstrap median 0.672
    e_wt_vc <- 0.512
    label("Power exponent of baseline body weight on central volume V1 (unitless; reference 66.0 kg)") # Table 2, row 'BWT effect on V 1' = 0.512 [95% CI 0.026-0.998]; bootstrap median 0.518

    # ---- Inter-individual variability. Chen 2019 Methods, 'Structural PK
    # model and variability models': IIV on CL, Q, V1 and V2 was log-normal
    # (exponential model) and the omega matrix was DIAGONAL, because no
    # significant correlation between CL and V1 was found during model
    # development. Table 2 prints each row as 'X omega 2 (%CV)', so the leading
    # number is the VARIANCE on the log scale and the parenthesised figure is
    # its square root as a percentage; sqrt(0.0687) = 0.262 -> 26%, and the
    # same identity holds for all four rows. The variances are therefore used
    # verbatim, with no CV-to-variance conversion.
    etalcl ~ 0.0687 # Table 2, row 'CL omega 2 (%CV)' = 0.0687 (26) [95% CI 0.054-0.084]; bootstrap median 0.067
    etalvc ~ 0.123 # Table 2, row 'V 1 omega 2 (%CV)' = 0.123 (35) [95% CI 0.061-0.185]; bootstrap median 0.122
    etalq ~ 0.528 # Table 2, row 'Q omega 2 (%CV)' = 0.528 (73) [95% CI 0.339-0.717]; bootstrap median 0.519
    etalvp ~ 1.08 # Table 2, row 'V 2 omega 2 (%CV)' = 1.08 (104) [95% CI 0.772-1.388]; bootstrap median 1.056

    # ---- Residual error. Chen 2019 Methods: 'The residual error was described
    # using an additive error model, after log-transformation of the PK data',
    # which is a log-normal residual, encoded here as lnorm(expSd).
    #
    # SCALE. Table 2 prints 'Res Add Err' = 0.292 without marking it as a
    # variance or an SD. It is the NONMEM $SIGMA VARIANCE, so the SD entered
    # here is sqrt(0.292) = 0.5404. Three independent lines of evidence:
    # (1) Chen 2019 Fig. 1b, the log-observed versus log-individual-predicted
    #     panel, is a direct picture of this residual. Digitising its blue
    #     scatter against the printed axes gives a 2.5th-to-97.5th percentile
    #     residual span of about +/- 1.2 natural-log units. An SD of 0.5404
    #     predicts +/- 1.06; an SD of 0.292 predicts only +/- 0.57, about half
    #     what the figure shows.
    # (2) Every other variance component in Table 2 is a NONMEM variance, and
    #     the table column is headed 'NONMEM results Estimate'. NONMEM reports
    #     $SIGMA on the variance scale.
    # (3) The relative standard error implied by the printed 95% CI (the table
    #     footnote states the CI is estimate +/- 1.96 x SE) is 7.5%. For a
    #     variance that implies a few hundred effective observations,
    #     consistent with correlated peak-and-trough sampling in 353 patients;
    #     read as an SD it would imply fewer than one independent observation
    #     per patient.
    # See the vignette Errata for the full argument and the digitisation.
    expSd <- 0.5404
    label("Log-normal residual error SD on the natural-log scale (unitless)") # Table 2, row 'Res Add Err' = 0.292 [95% CI 0.249-0.335] read as a variance; sqrt(0.292) = 0.5404
  })

  model({
    # Individual PK parameters. Baseline body weight enters CL and V1 only, as
    # the Eq. 1 normalised power term centred on the 66.0 kg arm median; Q and
    # V2 carry no covariate (Chen 2019 Table 2 lists no BWT effect for either).
    # Inter-individual variability is exponential on all four parameters.
    cl <- exp(lcl + etalcl) * (WT / 66.0)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 66.0)^e_wt_vc
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
