Takechi_2025_nemolizumab_ppnrs <- function() {
  description <- paste0(
    "Population PD model of the weekly average Peak Pruritus Numerical Rating ",
    "Scale (PP-NRS, 0-10) under nemolizumab in Japanese adolescent and adult ",
    "patients with prurigo nodularis. The score falls from each subject's own ",
    "observed baseline by two independent mono-exponential reductions that ",
    "approach constant maxima: a placebo response Pmax * (1 - exp(-Kp * t)) ",
    "present in every arm, and a nemolizumab response Edrug * (1 - exp(-Kd * ",
    "t)) present only in the active arms. Both maxima are CONSTANTS: the ",
    "source tested an exposure-response relationship and found none, because ",
    "the 30 mg and 60 mg Q4W arms of the pivotal phase II/III study responded ",
    "alike, so this model has no PK layer and no dose-response - the ",
    "DOSE_NEMOLIZUMAB_MG column acts only as an on/off treatment switch. ",
    "Nemolizumab is roughly ten times faster than placebo (Kd 28.7e-3 vs Kp ",
    "2.74e-3 per day) as well as larger. IMPORTANT READING OF TABLE 2: the ",
    "printed Pmax = 0.618 and Edrug = 1.41 are LOG-SCALE thetas of the ",
    "source's exponential IIV parameterisation, so the typical maxima are ",
    "exp(0.618) = 1.86 and exp(1.41) = 4.10 PP-NRS points; taking the printed ",
    "numbers as point reductions misses the observed active-arm median in the ",
    "source's own Fig. 2 VPC by 3.7 points. Log-normal interindividual ",
    "variability is estimated on all four parameters and is very large, ",
    "particularly on the placebo onset rate; no covariate reached ",
    "significance. Sister models from the same paper: ",
    "modellib('Takechi_2025_nemolizumab') (popPK) and ",
    "modellib('Takechi_2025_nemolizumab_mbma_iga') (MBMA of IGA success rates)."
  )

  reference <- paste(
    "Takechi T, Shimizu J, Kabashima K, Ieiri I.",
    "Quantitative evaluation of nemolizumab pharmacokinetics and efficacy in",
    "prurigo nodularis: a population pharmacokinetics and model-based",
    "meta-analysis approach.",
    "Dermatol Ther (Heidelb). 2025;15(12):3615-3632.",
    "doi:10.1007/s13555-025-01554-4.",
    "Structural model: Methods 'Pharmacodynamic Model', the two displayed",
    "'Placebo response' / 'Drug response' equations.",
    "Parameter values: Table 2, 'PopPD' block, 'Original data / Estimate'",
    "column.",
    sep = " "
  )

  vignette <- "Takechi_2025_nemolizumab"

  units <- list(
    time          = "day",
    dosing        = paste0(
      "mg/administration (nemolizumab 30 mg Q4W after a 60 mg loading dose, ",
      "or 60 mg Q4W, supplied through the DOSE_NEMOLIZUMAB_MG covariate ",
      "column and NOT as rxode2 dose events; the source found no ",
      "dose-response, so any non-zero value produces the same drug effect)"
    ),
    concentration = paste0(
      "points/subject (weekly average Peak Pruritus Numerical Rating Scale ",
      "score on the 0-10 patient-reported scale, for one subject; this output ",
      "is NOT a drug concentration - the slash satisfies ",
      "checkModelConventions parsing)"
    )
  )

  covariateData <- list(
    SCORE_PPNRS = list(
      description        = paste0(
        "Subject's own observed weekly average Peak Pruritus Numerical Rating ",
        "Scale score at baseline, used as the initial value of the response."
      ),
      units              = "(score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Takechi 2025 Results 'PopPD Analysis': 'Model initialization was ",
        "based on the observed weekly average PP-NRS scores at baseline.' The ",
        "baseline is therefore DATA, not an estimated parameter, and carries ",
        "no typical value and no eta of its own - which is why Table 2's ",
        "PopPD block has four structural rows and no baseline row. Median 8.6 ",
        "with a range of 6.4-10 in the 229-patient PopPD dataset (Table 1, ",
        "'Patients with PN' column, footnote b). The daily PP-NRS is averaged ",
        "over each week before use (Methods 'Pharmacokinetic and ",
        "Pharmacodynamic Assessments')."
      ),
      source_name        = "weekly average PP-NRS at baseline (Takechi 2025 Results 'PopPD Analysis'; Table 1 'PP-NRS score' row)"
    ),
    DOSE_NEMOLIZUMAB_MG = list(
      description        = paste0(
        "Per-subject assigned subcutaneous nemolizumab maintenance dose; 0 ",
        "identifies a placebo subject."
      ),
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Acts ONLY as an on/off switch in this model: the drug effect enters ",
        "as Edrug * (DOSE_NEMOLIZUMAB_MG > 0), with no dose-response. Takechi ",
        "2025 Discussion: 'an initial attempt was made to evaluate the ",
        "exposure-response relationship; however, no clear association ",
        "between exposure and response was identified ... two dosage groups, ",
        "30 mg and 60 mg, both of which exhibited similar responses.' The ",
        "values used in the source study are 0 (placebo, 76 subjects), 30 ",
        "(with a 60 mg loading dose) and 60. Because there is no ",
        "dose-response, entering 30 or 60 gives an identical prediction; the ",
        "column is kept numeric rather than binary so the arm is recorded ",
        "faithfully and so a future dose-response extension has somewhere to ",
        "read the dose from."
      ),
      source_name        = "treatment arm (Takechi 2025 Results 'Data Summary'; Table 3 Yokozeki et al. rows)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description        = "Body weight. Screened as a PD covariate but NOT retained in the final model.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Takechi 2025 Methods 'Pharmacodynamic Model' lists body weight among ",
        "the covariates evaluated on the PD parameters. Results 'PopPD ",
        "Analysis': 'No covariates were included in the model, as no ",
        "statistically significant effects were identified.' No point ",
        "estimate is reported. Median 61.2 kg, range 32.8-109.6 (Table 1, ",
        "footnote b)."
      ),
      source_name        = "body weight (Takechi 2025 Methods 'Pharmacodynamic Model' covariate list)"
    ),
    ALB = list(
      description        = "Serum albumin. Screened as a PD covariate but NOT retained in the final model.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Listed as 'ALB' in the Takechi 2025 Methods 'Pharmacodynamic Model' ",
        "covariate screen; no effect retained and no point estimate reported. ",
        "ALB IS retained in the companion popPK model on CL/F - see ",
        "modellib('Takechi_2025_nemolizumab')."
      ),
      source_name        = "ALB (Takechi 2025 Methods 'Pharmacodynamic Model' covariate list)"
    ),
    AGE = list(
      description        = "Patient age. Screened as a PD covariate but NOT retained in the final model.",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Listed in the Takechi 2025 Methods 'Pharmacodynamic Model' covariate ",
        "screen; no effect retained and no point estimate reported. Median 51 ",
        "years, range 13-84 (Table 1, footnote b)."
      ),
      source_name        = "age (Takechi 2025 Methods 'Pharmacodynamic Model' covariate list)"
    ),
    SEXF = list(
      description        = "Sex. Screened as a PD covariate but NOT retained in the final model.",
      units              = "(binary)",
      type               = "binary",
      reference_category = NULL,
      notes              = paste0(
        "Listed as 'sex' in the Takechi 2025 Methods 'Pharmacodynamic Model' ",
        "covariate screen; no effect retained and no point estimate reported. ",
        "The PopPD cohort is 107 male / 122 female (Table 1, footnote b)."
      ),
      source_name        = "sex (Takechi 2025 Methods 'Pharmacodynamic Model' covariate list)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 229L,
    n_studies      = 1L,
    n_observations = 3833L,
    age_range      = "13-84 years",
    age_median     = "51 years",
    weight_range   = "32.8-109.6 kg",
    weight_median  = "61.2 kg",
    sex_female_pct = 53.3,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste0(
      "Prurigo nodularis diagnosed more than 6 months earlier, with limb ",
      "lesions, at least 20 bilateral prurigo nodules, and inadequate ",
      "response to (or inability to receive) high-potency topical ",
      "corticosteroids for at least 4 weeks and oral antihistamines for at ",
      "least 2 weeks. Baseline severity: median weekly average PP-NRS 8.6 ",
      "(range 6.4-10) and Investigator's Global Assessment 3 in 112 patients ",
      "and 4 in 117 patients."
    ),
    dose_range     = paste0(
      "Subcutaneous nemolizumab 30 mg Q4W after a 60 mg loading dose, 60 mg ",
      "Q4W, or placebo, for 16 weeks (initial treatment phase only). Arm ",
      "sizes 77 / 76 / 76 (Table 3, Yokozeki et al. rows)."
    ),
    regions        = "Japan",
    notes          = paste0(
      "Single study: the phase II/III randomized, placebo-controlled, ",
      "double-blind, multicentre trial M525101-11 (jRCT2011200017; Yokozeki ",
      "et al. 2024, Br J Dermatol 191:200-208), in patients aged 13 years and ",
      "older. Demographics are Takechi 2025 Table 1, 'Patients with PN' ",
      "column, footnote-b rows (the PopPD analysis dataset). The 229 subjects ",
      "are the whole randomized population; the 153 non-placebo subjects of ",
      "the same study are the ones that entered the companion popPK external ",
      "validation. PP-NRS was recorded daily to day 112 and averaged by week ",
      "before modelling."
    )
  )

  ini({
    # ========================================================================
    # All values are Takechi 2025 Table 2, 'PopPD' block, 'Original data /
    # Estimate' column. Bootstrap medians and 95% CIs in the adjacent columns
    # are used only as a transcription check.
    #
    # SIGN CONVENTION. The source defines Pmax and Edrug as "the maximum
    # effects, defined as the reduction in PP-NRS scores attributable to the
    # placebo and nemolizumab" (Methods 'Pharmacodynamic Model'). Both are
    # therefore POSITIVE reductions and are SUBTRACTED from baseline in
    # model().
    #
    # ========================================================================
    # Pmax AND Edrug ARE LOG-SCALE THETAS, NOT POINT REDUCTIONS.
    #
    # This is the single most consequential reading in this file, so the
    # evidence is set out in full. Table 2 prints Pmax = 0.618 and Edrug =
    # 1.41 with no units. Taking those at face value as PP-NRS point
    # reductions is WRONG; they are the log-scale thetas of the source's
    # stated exponential (log-normal) IIV parameterisation, so the typical
    # values are exp(0.618) = 1.855 and exp(1.41) = 4.096 points.
    #
    # 1. The units column. In the same table Kp and Kd are labelled
    #    "(x 10^-3/day)" while Pmax and Edrug carry no unit at all. A
    #    log-scale theta is unitless; a PP-NRS point reduction is not.
    #
    # 2. Fig. 2 falsifies the linear reading outright. Simulating the source's
    #    own 229-patient design and comparing to the observed percentiles that
    #    Fig. 2 plots, at days 0/28/56/84/112:
    #      nemolizumab arm, observed median   8.5 / 6.0 / 4.6 / 3.6 / 3.1
    #      linear reading   simulated median  8.6 / 7.3 / 7.0 / 6.9 / 6.8
    #      log-scale        simulated median  8.6 / 5.6 / 4.8 / 4.3 / 3.9
    #    The linear reading misses the day-112 median by 3.7 points, which is
    #    more than the entire observed placebo-arm change. It cannot be
    #    rescued by the IIV: with Edrug = 1.41 as a point reduction the median
    #    subject's drug effect is CAPPED at 1.41 points however large the etas
    #    are, and the observed median falls by about 5.4.
    #
    # 3. The log-scale reading also reproduces the 95th percentile
    #    (simulated 10.4/9.0/8.7/8.5/8.4 vs observed 10/8.9/8.5/8.9/8.2) and
    #    the placebo arm's median (7.6 vs 7.7 at day 112) and its descending
    #    lower prediction band, none of which the linear reading matches.
    #
    # 4. It reproduces the source's own stated MISFIT. Results 'PopPD
    #    Analysis': "The observed 50th and 95th percentiles were almost within
    #    the 95% prediction interval ... The simulated 5th percentile slightly
    #    underestimated the severity of pruritus during early treatment,
    #    suggesting a modest overprediction of early response in some
    #    patients." That is exactly the log-scale signature - its 50th and
    #    95th percentiles land on the observed values while its 5th percentile
    #    runs far too low (below zero in the active arm, which is why Fig. 2's
    #    lower band sits clipped against the axis). The linear reading has the
    #    opposite deficiency: its 5th percentile fits and its median does not.
    #
    # Kp and Kd are NOT log-scale. They are printed with units and with an
    # explicit 10^-3 multiplier, and reading them as logs would give
    # exp(0.00274) = 1.003 /day, i.e. a placebo effect fully expressed within
    # a few days - flatly contradicted by Fig. 2A, where the placebo arm
    # declines slowly across the whole 112-day observation period. They are
    # therefore wrapped in log() below while Pmax and Edrug are not.
    #
    # The vignette runs this whole comparison as a reproducible check.
    # ========================================================================
    #
    # OMEGA SCALE. Table 2 prints the IIV rows as bare numbers with no CV%
    # label; they are VARIANCES. The proof is in the same table's PopPK block,
    # where reading the IIV rows as variances is the only reading under which
    # the printed CL/F-V/F covariance yields a legal correlation. The PopPD
    # variances are extremely large (omega for Kp is sqrt(10.4) = 3.2 on the
    # log scale); this is a real feature of the fit, not a transcription slip,
    # and it is what lets the model's population mean response be much larger
    # than its typical-subject response. The vignette checks the simulated
    # cohort against the source's own Fig. 2 VPC.
    # ========================================================================

    lpmax <- 0.618
    label("Log maximum placebo effect Pmax, as a reduction in PP-NRS points (log points); back-transform Pmax = exp(0.618) = 1.855 points")
    # Table 2 Pmax = 0.618 (%RSE 17.2; bootstrap median 0.610, 95% CI
    # 0.431-0.894). Already on the log scale - NOT wrapped in log(). See the
    # log-scale block above.

    lkp <- log(2.74e-3)
    label("Log rate constant for the onset of the placebo effect (log 1/day); back-transform Kp = 0.00274 1/day, onset half-life 253 days")
    # Table 2 Kp = 2.74 x 10^-3 /day (%RSE 47.4; bootstrap median 2.75, 95% CI
    # 1.04-6.70). The column header carries the 10^-3 multiplier.

    lemax <- 1.41
    label("Log maximum nemolizumab effect Edrug, as a reduction in PP-NRS points (log points); back-transform Edrug = exp(1.41) = 4.096 points")
    # Table 2 Edrug = 1.41 (%RSE 9.9; bootstrap median 1.41, 95% CI
    # 1.14-1.73). Already on the log scale - NOT wrapped in log(). See the
    # log-scale block above. exp(1.41) = 4.096 points is 2.2 times the typical
    # placebo maximum of exp(0.618) = 1.855 points, and because nemolizumab is
    # also about ten times faster (Kd / Kp = 10.5) it expresses essentially
    # all of that maximum within the 112-day observation window while the
    # placebo expresses only about a quarter of its own.

    lkdrug <- log(28.7e-3)
    label("Log rate constant for the onset of the nemolizumab effect (log 1/day); back-transform Kd = 0.0287 1/day, onset half-life 24.1 days")
    # Table 2 Kd = 28.7 x 10^-3 /day (%RSE 13.8; bootstrap median 28.7, 95% CI
    # 21.8-37.6). The column header carries the 10^-3 multiplier. Kd / Kp =
    # 10.5, the source's basis for describing nemolizumab as having a rapid
    # onset of itch relief relative to placebo.

    # ---- Interindividual variability ----------------------------------------
    # "Interindividual variability (IIV) of each parameter assumed a
    # log-normal distribution and was estimated with an exponential error
    # model" (Methods 'Pharmacodynamic Model'), so each eta is additive on the
    # log scale. Diagonal: Table 2 reports no PopPD covariance row.
    etalpmax  ~ 1.55
    # Table 2 IIV Pmax = 1.55 (%RSE 39.9; bootstrap median 1.56, 95% CI
    # 0.641-3.18).
    etalkp    ~ 10.4
    # Table 2 IIV Kp = 10.4 (%RSE 26.2; bootstrap median 10.2, 95% CI
    # 5.80-16.6).
    etalemax  ~ 1.19
    # Table 2 IIV Edrug = 1.19 (%RSE 15.3; bootstrap median 1.18, 95% CI
    # 0.812-1.59).
    etalkdrug ~ 2.10
    # Table 2 IIV Kd = 2.10 (%RSE 19.6; bootstrap median 2.09, 95% CI
    # 1.37-3.06).

    # ---- Residual error -----------------------------------------------------
    addSd <- 0.762889
    label("Additive residual standard deviation on the weekly average PP-NRS score (points)")
    # Table 2 'Additive error' = 0.582 (%RSE 5.0; bootstrap median 0.580, 95%
    # CI 0.525-0.635), read as a VARIANCE, so the SD is sqrt(0.582) =
    # 0.762889.
    #
    # WHY VARIANCE, NOT SD. Table 2 is a NONMEM output table and every other
    # variability row in it is on the raw NONMEM scale: the PopPK IIV rows are
    # variances (proved by the covariance row - see the companion popPK model
    # file), and the ONLY row the authors back-transformed for reporting is
    # explicitly relabelled 'log normal error (CV%)'. The 'Additive error' row
    # carries no such relabelling, so it is the raw $SIGMA, i.e. a variance.
    # Reading it instead as an SD of 0.582 changes the residual by 24% and is
    # documented as the alternative in the vignette's Assumptions and
    # deviations section; it does not change any structural conclusion.
  })

  model({
    # ---- 1. Individual parameters -------------------------------------------
    # Log-normal IIV throughout (Methods: exponential error model).
    pmaxPlacebo <- exp(lpmax  + etalpmax)
    kp          <- exp(lkp    + etalkp)
    edrug       <- exp(lemax  + etalemax)
    kdrug       <- exp(lkdrug + etalkdrug)

    # ---- 2. Treatment switch -------------------------------------------------
    # Constant drug effect above zero dose; there is no exposure-response.
    onDrug <- (DOSE_NEMOLIZUMAB_MG > 0)

    # ---- 3. Response components (Methods 'Pharmacodynamic Model') -----------
    #   Placebo response_i(t) = Pmax * (1 - exp(-Kp * time))
    #   Drug response_i(t)    = Edrug * (1 - exp(-Kd * time))
    # Both are reductions in the PP-NRS score. `time` is days since the first
    # dose, matching the 1/day units of Kp and Kd and the 0-112 day axis of
    # the source's Fig. 2.
    placeboResponse <- pmaxPlacebo * (1 - exp(-kp    * time))
    drugResponse    <- edrug       * (1 - exp(-kdrug * time)) * onDrug

    # ---- 4. Observation and error -------------------------------------------
    # The score starts at the subject's own observed baseline and is reduced by
    # the placebo response in every arm plus the drug response in the active
    # arms. The source prints the two response components but never prints the
    # line that combines them with the baseline; that this is the ADDITIVE
    # combination below follows from the source's own definition of Pmax and
    # Edrug as reductions in PP-NRS points, and it is checked against the Fig.
    # 2 VPC in the vignette. See the vignette's Assumptions and deviations
    # section.
    #
    # No 0-10 clamp is applied: the source describes none, and clamping would
    # change the residual distribution near the floor. Simulated values below
    # 0 are possible for strong responders and are an artefact of the additive
    # residual, exactly as in the source's own simulations (its Fig. 2 lower
    # prediction interval reaches 0).
    #
    # NAMING. The single output is named `Cc` because that is the package's
    # canonical single-output observation variable (R/conventions.R
    # `observationVar`), NOT because it is a concentration - it is a weekly
    # average PP-NRS score in points. Same convention as
    # `Goteti_2024_SLE_mbma.R`, whose `Cc` is a responder probability. The
    # descriptively-named `ppnrs` below is exposed as a derived model variable
    # so it appears in the rxSolve output alongside `Cc`.
    ppnrs <- SCORE_PPNRS - placeboResponse - drugResponse
    Cc <- ppnrs
    Cc ~ add(addSd)
  })
}
