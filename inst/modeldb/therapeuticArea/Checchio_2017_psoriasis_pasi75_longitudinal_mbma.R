Checchio_2017_psoriasis_pasi75_longitudinal_mbma <- function() {
  description <- paste0(
    "MBMA. LONGITUDINAL model-based meta-analysis of the PASI75 responder ",
    "rate (proportion of patients achieving a 75% reduction from baseline in ",
    "the Psoriasis Area and Severity Index) over time in moderate-to-severe ",
    "plaque psoriasis, fitted to study-arm summary data from 57 trials of ",
    "systemic agents published 1998-2015 plus an internal Pfizer tofacitinib ",
    "database. The logit of the PASI75 rate is the sum of a placebo component ",
    "that rises exponentially to a plateau, f0 = A + B*(1 - exp(-kpbo*t)), and ",
    "a drug component that is an Emax function of dose multiplied by its own ",
    "exponential onset, fdrug = Emax_d*(1 - exp(-kdrug_d*t))*Dose/(Dose + ",
    "ED50_d), with Emax, ED50 and kdrug estimated separately for each of 13 ",
    "drugs (adalimumab, certolizumab, etanercept, infliximab, briakinumab, ",
    "ustekinumab, brodalumab, ixekizumab, secukinumab, alefacept, apremilast, ",
    "methotrexate, tofacitinib). Body weight acts as a power function ",
    "centred at the 90 kg dataset median on BOTH the placebo plateau ",
    "(exponent -0.198) and the drug onset rate (exponent -0.834): heavier ",
    "patients respond less and more slowly. Dose enters through one ",
    "CONMED_<drug>_DOSE covariate column per drug rather than through rxode2 ",
    "dose events; there is no PK layer. The model reproduces every PASI75 ",
    "value in the source Table 1 to within 1.1 percentage points and every ",
    "ET50 to within 0.3 weeks. A secondary output prob_pasi90 reproduces the ",
    "source's external validation (Figure 5), which predicts the PASI90 time ",
    "course by importing two scaling factors from the companion landmark ",
    "model. Simulation scope is STUDY-ARM-MEAN responder trajectories, NOT ",
    "individual patients. The companion dose-response model from the same ",
    "paper is modellib('Checchio_2017_psoriasis_pasi_landmark_mbma')."
  )

  reference <- paste(
    "Checchio T, Ahadieh S, Gupta P, Mandema J, Puig L, Wolk R, Valdez H,",
    "Tan H, Krishnaswami S, Tallman A, Kaur M, Ito K.",
    "Quantitative Evaluations of Time-Course and Treatment Effects of Systemic",
    "Agents for Psoriasis: A Model-Based Meta-Analysis.",
    "Clin Pharmacol Ther. 2017;102(6):1006-1016.",
    "doi:10.1002/cpt.732. PMC5697570.",
    "Structural model: Equations 1-4 of Methods 'Longitudinal model'.",
    "Covariate model: Equation 9 of Methods 'Covariate model'.",
    "Residual model: Equations 11-13 of Methods 'Residual error model'.",
    "All display equations are rasterised in the published PDF and are LOST",
    "from the preprocessed markdown; they were recovered with",
    "pdftotext -layout. Drug-specific parameter values are Supplementary",
    "Table S1.1 and common parameter values are Supplementary Table S1.2, both",
    "in the Supplementary Appendix (CPT-102-1006-s001.docx), obtained from the",
    "EuropePMC supplementaryFiles endpoint for PMC5697570.",
    sep = " "
  )

  vignette <- "Checchio_2017_psoriasis_systemic_agents"

  # `eta_study_corr` is the arm-level correlated RESIDUAL of Checchio 2017
  # Equation 13, not inter-study variability on a structural parameter: it has
  # no typical value to pair with, because it is one of the two multipliers on
  # the binomial standard error W and the source tabulates only its variance
  # ('sigma corr'). It is realised once per study arm and is therefore an eta
  # rather than an epsilon in nlmixr2 -- the same device and the same naming
  # problem as the 'correlation coefficient' row of
  # Chen_2025_methotrexate_acr50_mbma. Declared here so the convention checker
  # does not look for a fixed effect named `corr`.
  paper_specific_etas <- c("eta_study_corr")

  # This model has no PK layer, no concentration and no rxode2 dose events:
  # dose is a per-arm covariate column. The placeholder `units` entries follow
  # the arm-level responder-rate MBMA convention already used by
  # Serrano_2026_atopicDermatitis_placebo_mbma and
  # Dodds_2013_psoriasis_biologics_mbma so that checkModelConventions() sees a
  # parseable dosing / concentration pair.
  units <- list(
    time = "week (weeks since first dose; kpbo and kdrug are reported in week^-1 and the source's primary read-outs are Weeks 4 and 12)",
    dosing = "mg/administration (dose per administration of the named agent, supplied in the CONMED_<drug>_DOSE covariate columns; infliximab is mg/kg per administration. This model consumes NO rxode2 dose events.)",
    concentration = "probability/arm (prob_pasi75 is the STUDY-ARM probability that a patient achieves a 75% reduction from baseline in the Psoriasis Area and Severity Index, on a 0-1 scale; it is NOT a drug concentration. The slash satisfies checkModelConventions unit parsing.)"
  )

  covariateData <- list(
    WT = list(
      description = "Study-arm MEAN body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "TRIAL-ARM-LEVEL, not subject-level: this is the arm's reported mean body weight, the grain at which an aggregate-data meta-analysis operates. Enters as the power ratio (WT/90)^theta of Checchio 2017 Equation 9, where 90 kg is stated in the same paragraph to be the approximate median of the dataset; the same 90 kg anchor is restated in the Table 1, Table 2 and Figure 4 footnotes ('generated assuming a typical weight of 90 kg'). The effect acts on BOTH the placebo plateau B (exponent -0.198) and the drug onset rate kdrug (exponent -0.834); both exponents are negative, i.e. heavier arms reach a lower plateau and get there more slowly, which is the direction the Discussion reports. Checchio 2017 Methods 'Covariate model' states that when body weight was not reported (10-15% of studies) WEIG was set to 90 so the covariate term collapses to 1; downstream simulation code should do the same for an arm with unknown weight. The source cautions explicitly (Discussion) that this is an AGGREGATE-level covariate effect whose magnitude is likely attenuated relative to the individual-level truth, and that the attribution to the placebo component 'should not be interpreted in a mechanistic context'.",
      source_name = "WEIG (Checchio 2017 Equations 9 and 10); 'body weight' (Supplementary Table S1.2 rows 'Body weight effect on placebo (B)' and 'Body weight effect on drug onset (kdrug)')"
    ),
    N_ARM = list(
      description = "Number of patients contributing to the study-arm PASI75 proportion at a given timepoint.",
      units = "participants",
      type = "count",
      reference_category = NULL,
      notes = "Study-design quantity supplied per observation row, not estimated. It is the meta-analytic weight of Checchio 2017 Equation 12, W = sqrt(Pr*(1 - Pr)/N), the binomial standard error of an arm proportion. Both residual terms of Equation 13 are multiplied by W, so the arm size scales BOTH the independent residual and the arm-level correlated residual. N_ARM is carried as a data column rather than applied downstream because the correlated term eta_study_corr is a RANDOM EFFECT whose SD therefore varies from arm to arm, and a per-arm random-effect SD cannot be reproduced by rescaling after a single rxSolve() -- this is the case the N_ARM register entry reserves the column for. Must be strictly positive. Checchio 2017 Methods: 'Because studies were weighted by inverse of standard errors, larger studies had more influence on estimating the parameter means.'",
      source_name = "N (Checchio 2017 Equations 11 and 12, 'the number of subjects for each arm within each study')"
    ),
    CONMED_ADALIMUMAB_DOSE = list(
      description = "Adalimumab dose per administration in the study arm; 0 if the arm did not receive adalimumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose PER ADMINISTRATION, not total daily or weekly dose. This is settled arithmetically: Table 1's clinical-dose column reads '40mg Q2W' and only Dose = 40 reproduces the published Week 4 and Week 12 PASI75 values. Zero for every arm not randomised to adalimumab, which collapses the adalimumab Emax term to exactly zero. Checchio 2017 Methods notes that doses were normalised to the approved (or, for unapproved agents, most frequently reported) regimen to handle titration schemes, so an arm on a non-standard titration is outside the studied range.",
      source_name = "Dose (Checchio 2017 Equation 4); 'Clinical dose' column of Table 1"
    ),
    CONMED_CERTOLIZUMAB_DOSE = list(
      description = "Certolizumab pegol dose per administration in the study arm; 0 if the arm did not receive certolizumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration; Table 1's clinical dose is 200 mg Q2W. Certolizumab is the one drug for which Supplementary Table S1.1 reports NO onset rate ('NAb -- Not estimated as only one study for certolizumab was published'), so its drug effect has no time-course term; see the ini() note.",
      source_name = "Dose (Checchio 2017 Equation 4); 'Clinical dose' column of Table 1"
    ),
    CONMED_ETANERCEPT_DOSE = list(
      description = "Etanercept dose per administration in the study arm; 0 if the arm did not receive etanercept.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration. Table 1 reports two etanercept arms, 25 mg and 50 mg twice weekly, that share one Emax / ED50 / kdrug triplet and differ only in this column; both are reproduced to within 1.2 percentage points, which is the cleanest internal check that the column carries dose per administration rather than a weekly total.",
      source_name = "Dose (Checchio 2017 Equation 4); 'Clinical dose' column of Table 1"
    ),
    CONMED_INFLIXIMAB_DOSE = list(
      description = "Infliximab dose per administration in the study arm, per kilogram of body weight; 0 if the arm did not receive infliximab.",
      units = "mg/kg",
      type = "continuous",
      reference_category = NULL,
      notes = "UNIQUE UNITS within this model: infliximab is dosed on a mg/kg basis and its ED50 of 0.315 is therefore in mg/kg, per the Supplementary Table S1.1 footnote 'ED50, dose (mg or mg/kg) to achieve half of the maximal response'. Table 1's clinical dose is 5 mg/kg Q8W, and Dose = 5 reproduces the published Week 4 (32.1% vs 32.7%) and Week 12 (75.6% vs 75.3%) values. Supplying a milligram amount here would overstate the dose roughly 90-fold and saturate the Emax term.",
      source_name = "Dose (Checchio 2017 Equation 4); 'Clinical dose' column of Table 1; Supplementary Table S1.1 footnote"
    ),
    CONMED_BRIAKINUMAB_DOSE = list(
      description = "Briakinumab dose per administration in the study arm; 0 if the arm did not receive briakinumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration; Table 1's clinical dose is 100 mg Q4W. Briakinumab was discontinued in development (Table 1 footnote b) and is retained here because the source models and predicts it.",
      source_name = "Dose (Checchio 2017 Equation 4); 'Clinical dose' column of Table 1"
    ),
    CONMED_USTEKINUMAB_DOSE = list(
      description = "Ustekinumab dose per administration in the study arm; 0 if the arm did not receive ustekinumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration; Table 1's clinical dose is 45 mg Q12W. Note that the approved ustekinumab regimen is itself weight-banded (45 mg at or below 100 kg, 90 mg above 100 kg), which the source handles in its Figure 4 simulation by splitting the cohort at 100 kg (Figure 4 caption); that banding is a DOSING rule applied to this column by downstream code, not a model term.",
      source_name = "Dose (Checchio 2017 Equation 4); 'Clinical dose' column of Table 1"
    ),
    CONMED_BRODALUMAB_DOSE = list(
      description = "Brodalumab dose per administration in the study arm; 0 if the arm did not receive brodalumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration; Table 1's clinical dose is 210 mg Q2W. Investigational at the time of the analysis (Table 1 footnote a).",
      source_name = "Dose (Checchio 2017 Equation 4); 'Clinical dose' column of Table 1"
    ),
    CONMED_IXEKIZUMAB_DOSE = list(
      description = "Ixekizumab dose per administration in the study arm; 0 if the arm did not receive ixekizumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration; Table 1's clinical dose is 80 mg Q4W. Ixekizumab has the highest predicted Week 12 PASI75 of any agent in the longitudinal analysis (81.7%) together with the shortest ET50 (3.7 weeks).",
      source_name = "Dose (Checchio 2017 Equation 4); 'Clinical dose' column of Table 1"
    ),
    CONMED_SECUKINUMAB_DOSE = list(
      description = "Secukinumab dose per administration in the study arm; 0 if the arm did not receive secukinumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration; Table 1's clinical dose is 150 mg QM.",
      source_name = "Dose (Checchio 2017 Equation 4); 'Clinical dose' column of Table 1"
    ),
    CONMED_ALEFACEPT_DOSE = list(
      description = "Alefacept dose per administration in the study arm; 0 if the arm did not receive alefacept. Any positive value selects the full alefacept effect, because no ED50 was estimable.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED: Supplementary Table S1.1 reports 'NAa -- Not estimated as limited dose-response data were available' for the alefacept ED50, so this model has no alefacept dose-response curve and the drug term is a single-step offset triggered by any positive dose, exactly as the source describes for agents without estimable dose-response (Methods 'Landmark model', which applies the same device). The magnitude therefore corresponds to Table 1's 10 mg QW clinical dose and must NOT be read as a prediction for any other alefacept dose. Alefacept has since been withdrawn from the market.",
      source_name = "Dose (Checchio 2017 Equation 4); 'Clinical dose' column of Table 1; Supplementary Table S1.1 footnote a"
    ),
    CONMED_APREMILAST_DOSE = list(
      description = "Apremilast dose per administration in the study arm; 0 if the arm did not receive apremilast.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration, NOT total daily dose: Table 1's clinical dose is 30 mg b.i.d. and Dose = 30 reproduces the published Week 4 (5.58% vs 5.66%) and Week 12 (29.3% vs 28.5%) values, whereas Dose = 60 does not.",
      source_name = "Dose (Checchio 2017 Equation 4); 'Clinical dose' column of Table 1"
    ),
    CONMED_MTX_DOSE = list(
      description = "Methotrexate dose per administration in the study arm; 0 if the arm did not receive methotrexate. Any positive value selects the full methotrexate effect, because no ED50 was estimable.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, for the same reason as alefacept: Supplementary Table S1.1 reports 'NAa' for the methotrexate ED50 because the traditional oral agents were titrated per patient on safety and efficacy grounds, leaving insufficient literature dose-response information (Methods 'Landmark model'). The magnitude corresponds to Table 1's 18 mg QW clinical dose. Register canonical CONMED_MTX_DOSE already exists for randomised methotrexate in immune-mediated-inflammatory-disease dose-response meta-analyses and is reused unchanged.",
      source_name = "Dose (Checchio 2017 Equation 4); 'Clinical dose' column of Table 1; Supplementary Table S1.1 footnote a"
    ),
    CONMED_TOFACITINIB_DOSE = list(
      description = "Tofacitinib dose per administration in the study arm; 0 if the arm did not receive tofacitinib.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration, NOT total daily dose: Table 1 reports 5 mg b.i.d. and 10 mg b.i.d. arms that share one Emax / ED50 / kdrug triplet, and Dose = 5 and Dose = 10 reproduce both published Week 12 values (41.7% vs 40.8% and 61.8% vs 60.7%) while 10 and 20 do not. Tofacitinib is the only agent whose data came from an internal Pfizer database rather than the literature (Methods 'Database development'), which the source states was done to maintain numerical accuracy; all contributing tofacitinib studies are published.",
      source_name = "Dose (Checchio 2017 Equation 4); 'Clinical dose' column of Table 1"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 57L,
    age_range = "not reported at arm level",
    weight_range = "not reported at arm level; the dataset median body weight is stated to be approximately 90 kg (Methods 'Covariate model'; Equation 9 denominator) and 90 kg is the typical weight used for every published prediction",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported at arm level",
    disease_state = "adults with moderate to severe plaque psoriasis enrolled in randomised placebo- or active-controlled trials",
    dose_range = "clinical doses per Table 1: adalimumab 40 mg Q2W, certolizumab 200 mg Q2W, etanercept 25 and 50 mg BIW, infliximab 5 mg/kg Q8W, brodalumab 210 mg Q2W, ixekizumab 80 mg Q4W, secukinumab 150 mg QM, briakinumab 100 mg Q4W, ustekinumab 45 mg Q12W, methotrexate 18 mg QW, tofacitinib 5 and 10 mg BID, alefacept 10 mg QW, apremilast 30 mg BID",
    timepoints = "all repeated measures of the arm PASI75 responder rate over each trial's double-blind period; the source's read-outs are Weeks 4 and 12, and alefacept data extend to 24 weeks",
    regions = "international; literature 1998-2015 identified by an Ovid Medline / Summary Basis of Approval / European Public Assessment Report search following the Cochrane approach, plus an internal Pfizer tofacitinib database",
    notes = "MBMA at the STUDY-ARM level: each modelled observation is one trial arm's PASI75 responder rate at one timepoint, so the random effect is BETWEEN-STUDY and this model must not be used to simulate individual patients. n_subjects is NA because the source reports patient counts only graphically (Figure 2 circle areas are proportional to arm N) and never tabulates a pooled total. The literature search yielded 912 abstracts, of which 151 studies were screened in; 57 studies survived into the longitudinal analysis and 71 into the companion landmark analysis (Figure 1 and Results 'Available data'). Cyclosporine, acitretin and baricitinib are in the landmark model but NOT here -- Results 'Available data' states they were excluded from the longitudinal model for insufficient longitudinal data -- which is why this file carries 13 drugs and its landmark companion carries 16 drug arms. Placebo-arm data pool oral and injectable/infusible trials (Figure 2 caption). Fitted in NONMEM 7.3 with the first-order conditional expectation method; the bootstrap used PsN 3.5.4."
  )

  ini({
    # ========================================================================
    # PROVENANCE FOR THE WHOLE ini() BLOCK
    #
    # Structural model, Checchio 2017 Methods 'Longitudinal model',
    # Equations 1-4. Every display equation in this paper is RASTERISED and is
    # dropped by the markdown preprocessor; the forms below were recovered
    # with `pdftotext -layout`:
    #
    #   (1)  Pr(PASI75) = g{f0 + fdrug}
    #   (2)  g          = 1 / (1 + exp(-(f0 + fdrug)))          [inverse logit]
    #   (3)  f0         = A + B * (1 - exp(-kpbo * time)) * exp(eta)
    #   (4)  fdrug      = Emax_d * (1 - exp(-kdrug_d * time)) *
    #                       Dose^gamma / (Dose^gamma + ED50_d^gamma)
    #
    # Covariate model, Equation 9:  WTeffect = (WEIG / 90)^theta,
    # applied to B and to kdrug per Supplementary Table S1.2 row labels.
    #
    # Common parameter values are Supplementary Table S1.2 and drug-specific
    # values are Supplementary Table S1.1, both from the Supplementary
    # Appendix (CPT-102-1006-s001.docx). The parenthesised number beside each
    # tabulated value is that table's CV(%) column, which is the RELATIVE
    # STANDARD ERROR of the estimate (it is reported for fixed effects such as
    # A and B, where a between-subject CV would be meaningless).
    #
    # NO HILL COEFFICIENT. Equation 4 carries an exponent gamma, and Methods
    # says 'The Hill coefficient (c) was also tested in the model.' No gamma
    # appears anywhere in Supplementary Table S1.1 or S1.2, so it was tested
    # and not retained; this model uses gamma = 1. That reading is confirmed
    # numerically rather than assumed: with gamma = 1 the model reproduces all
    # 15 Week-4 PASI75 values of Table 1 to within 0.7 percentage points, all
    # 15 Week-12 values to within 1.1 points, and all 15 ET50 values to within
    # 0.3 weeks. The companion LANDMARK model does estimate a Hill
    # coefficient, so the absence here is a real difference between the two
    # analyses, not a reporting gap.
    #
    # SOURCE-TRACE CONFIRMATION. Nothing below is tuned. The parameters are
    # transcribed from the two supplementary tables and then checked against
    # 45 published numbers that are NOT parameter estimates (Table 1's Week 4,
    # Week 12 and ET50 columns for all 15 drug / dose arms); the vignette
    # tabulates every one.
    # ========================================================================

    # ---- Placebo component (Equation 3; Supplementary Table S1.2) ----------
    a_pbo <- -7.61
    label("Intercept of the placebo effect on the PASI75 logit scale (paper: A; unitless log-odds)")
    # Supplementary Table S1.2, 'Intercept of placebo effect (A)' = -7.61
    # (CV 5.11%). expit(-7.61) = 0.049%, i.e. essentially nobody is a PASI75
    # responder at time zero, as must be the case since PASI75 is a 75%
    # reduction from each patient's own baseline.

    b_pbo <- 4.84
    label("Asymptote of the placebo effect, i.e. the plateau increment on the PASI75 logit scale for a 90 kg arm (paper: B; unitless log-odds)")
    # Supplementary Table S1.2, 'Asymptote of placebo effect (B)' = 4.84
    # (CV 7.56%). At 90 kg the weight term is 1, so the placebo logit
    # plateaus at -7.61 + 4.84 = -2.77, i.e. a 5.9% PASI75 rate; at Week 12 it
    # has reached 4.68%, which is the placebo response every drug prediction
    # in Table 1 sits on top of.

    lkpbo <- log(0.249)
    label("Log of the placebo-effect onset rate constant (paper: kpbo; back-transform 0.249 /week)")
    # Supplementary Table S1.2, 'Rate of onset of placebo effect (kpbo)' =
    # 0.249 (CV 7.31%). Held in log space so the rate stays positive; the
    # source reports it on the linear scale. Half-time to the placebo plateau
    # is log(2)/0.249 = 2.78 weeks.

    # ---- Body weight effects (Equation 9; Supplementary Table S1.2) --------
    # Equation 9 is a POWER model, WTeffect = (WEIG/90)^theta, so both values
    # below are EXPONENTS, not fractional changes. Both are negative: a
    # heavier arm reaches a lower placebo plateau AND has a slower drug onset,
    # which is the direction the Discussion reports ('heavier patients tend to
    # achieve lower efficacy and may also experience slower onset of effect
    # compared with lighter patients'). This sign agreement is the only
    # available check on the covariate orientation, because every published
    # prediction is made at exactly 90 kg where both terms collapse to 1.
    e_wt_b_pbo <- -0.198
    label("Power exponent on (WT/90) for the placebo asymptote B (unitless)")
    # Supplementary Table S1.2, 'Body weight effect on placebo (B)' = -0.198
    # (CV 17.9%). The Results text quotes a slightly different value, -0.193
    # with 95% CI [-0.277, -0.119], which is the MEDIAN OF A NONPARAMETRIC
    # BOOTSTRAP (N = 1,000), not the NONMEM point estimate. The supplementary
    # table is the final-estimates table and is used here; the bootstrap
    # median is recorded in the vignette. Both exclude zero.

    e_wt_kdrug <- -0.834
    label("Power exponent on (WT/90) for every drug onset rate constant kdrug (unitless)")
    # Supplementary Table S1.2, 'Body weight effect on drug onset (kdrug)' =
    # -0.834 (CV 19.3%). The Results text quotes the bootstrap median -0.867
    # with 95% CI [-1.32, -0.496]. ONE exponent is shared by all 13 drugs:
    # Supplementary Table S1.2 is headed 'Common Parameter Estimates' and
    # lists a single row, and Supplementary Table S1.1 carries no per-drug
    # weight column. Applied to kdrug only, NOT to Emax or ED50 -- Results
    # 'Landmark model' states that the weight effect gave no significant
    # improvement on any drug-effect magnitude term.

    # ========================================================================
    # DRUG-SPECIFIC PARAMETERS -- Supplementary Table S1.1, 'Drug-Specific
    # Parameter Estimates from NONMEM Output'. Each row gives Emax, ED50 and
    # kdrug with its CV(%) in parentheses. ED50 values are in mg per
    # administration, EXCEPT infliximab which is mg/kg (Table S1.1 footnote:
    # 'ED50, dose (mg or mg/kg) to achieve half of the maximal response').
    # ED50 is held in log space so it stays positive; the source reports it on
    # the linear scale. kdrug is likewise held in log space.
    #
    # TWO STRUCTURAL GAPS, both flagged by the source's own footnotes:
    #   * methotrexate and alefacept have NO ED50 ('NAa -- Not estimated as
    #     limited dose-response data were available'), so their drug term is a
    #     single-step offset at any positive dose.
    #   * certolizumab has NO kdrug ('NAb -- Not estimated as only one study
    #     for certolizumab was published'), so its drug term has no onset
    #     time course and steps to full effect immediately after time zero.
    #     That is the kdrug -> infinity limit of Equation 4, and it is the
    #     reading the published numbers force: Table 1's certolizumab Week 4
    #     PASI75 of 38.8% requires a drug term of 4.102 logit units, which
    #     ALREADY EXCEEDS the term's own ceiling Emax * Dose/(Dose + ED50) =
    #     4.097, so no finite onset rate can reproduce it. The step reading
    #     gives 38.7% at Week 4, 74.7% at Week 12 (published 74.5%), ET50 4.1
    #     weeks (published 4.2) and ET90 9.8 weeks (published 9.9).
    # ========================================================================

    # ---- TNF-alpha inhibitors ----
    emax_adalimumab <- 5.84
    label("Maximum adalimumab effect on the PASI75 logit scale (unitless log-odds)")  # Supplementary Table S1.1, adalimumab Emax = 5.84 (CV 4.28%)
    led50_adalimumab <- log(22.2)
    label("Log adalimumab ED50 (log mg per administration); back-transform 22.2 mg per administration")  # Supplementary Table S1.1, adalimumab ED50 = 22.2 mg (CV 13.1%)
    lkdrug_adalimumab <- log(0.373)
    label("Log adalimumab drug-effect onset rate (log week^-1); back-transform 0.373 /week")  # Supplementary Table S1.1, adalimumab kdrug = 0.373 (CV 19.6%)

    emax_certolizumab <- 5.30
    label("Maximum certolizumab effect on the PASI75 logit scale (unitless log-odds)")  # Supplementary Table S1.1, certolizumab Emax = 5.30 (CV 1.60%)
    led50_certolizumab <- log(58.7)
    label("Log certolizumab ED50 (log mg per administration); back-transform 58.7 mg per administration")  # Supplementary Table S1.1, certolizumab ED50 = 58.7 mg (CV 2.37%)

    emax_etanercept <- 4.53
    label("Maximum etanercept effect on the PASI75 logit scale (unitless log-odds)")  # Supplementary Table S1.1, etanercept Emax = 4.53 (CV 9.62%)
    led50_etanercept <- log(22.0)
    label("Log etanercept ED50 (log mg per administration); back-transform 22.0 mg per administration")  # Supplementary Table S1.1, etanercept ED50 = 22.0 mg (CV 22.1%)
    lkdrug_etanercept <- log(0.244)
    label("Log etanercept drug-effect onset rate (log week^-1); back-transform 0.244 /week")  # Supplementary Table S1.1, etanercept kdrug = 0.244 (CV 17.7%)

    emax_infliximab <- 4.41
    label("Maximum infliximab effect on the PASI75 logit scale (unitless log-odds)")  # Supplementary Table S1.1, infliximab Emax = 4.41 (CV 7.72%)
    led50_infliximab <- log(0.315)
    label("Log infliximab ED50; back-transform 0.315 mg/kg per administration (NOT mg)")  # Supplementary Table S1.1, infliximab ED50 = 0.315 mg/kg (CV 79.0%)
    lkdrug_infliximab <- log(0.626)
    label("Log infliximab drug-effect onset rate (log week^-1); back-transform 0.626 /week")  # Supplementary Table S1.1, infliximab kdrug = 0.626 (CV 13.2%)

    # ---- IL-12/23 inhibitors ----
    emax_briakinumab <- 5.49
    label("Maximum briakinumab effect on the PASI75 logit scale (unitless log-odds)")  # Supplementary Table S1.1, briakinumab Emax = 5.49 (CV 3.83%)
    led50_briakinumab <- log(22.6)
    label("Log briakinumab ED50 (log mg per administration); back-transform 22.6 mg per administration")  # Supplementary Table S1.1, briakinumab ED50 = 22.6 mg (CV 19.5%)
    lkdrug_briakinumab <- log(0.292)
    label("Log briakinumab drug-effect onset rate (log week^-1); back-transform 0.292 /week")  # Supplementary Table S1.1, briakinumab kdrug = 0.292 (CV 17.1%)

    emax_ustekinumab <- 4.35
    label("Maximum ustekinumab effect on the PASI75 logit scale (unitless log-odds)")  # Supplementary Table S1.1, ustekinumab Emax = 4.35 (CV 4.78%)
    led50_ustekinumab <- log(5.05)
    label("Log ustekinumab ED50 (log mg per administration); back-transform 5.05 mg per administration")  # Supplementary Table S1.1, ustekinumab ED50 = 5.05 mg (CV 67.8%)
    lkdrug_ustekinumab <- log(0.235)
    label("Log ustekinumab drug-effect onset rate (log week^-1); back-transform 0.235 /week")  # Supplementary Table S1.1, ustekinumab kdrug = 0.235 (CV 11.5%)

    # ---- IL-17 inhibitors ----
    emax_brodalumab <- 7.18
    label("Maximum brodalumab effect on the PASI75 logit scale (unitless log-odds)")  # Supplementary Table S1.1, brodalumab Emax = 7.18 (CV 3.60%)
    led50_brodalumab <- log(129)
    label("Log brodalumab ED50 (log mg per administration); back-transform 129 mg per administration")  # Supplementary Table S1.1, brodalumab ED50 = 129 mg (CV 7.02%)
    lkdrug_brodalumab <- log(1.16)
    label("Log brodalumab drug-effect onset rate (log week^-1); back-transform 1.16 /week")  # Supplementary Table S1.1, brodalumab kdrug = 1.16 (CV 13.1%)

    emax_ixekizumab <- 5.12
    label("Maximum ixekizumab effect on the PASI75 logit scale (unitless log-odds)")  # Supplementary Table S1.1, ixekizumab Emax = 5.12 (CV 3.58%)
    led50_ixekizumab <- log(10.5)
    label("Log ixekizumab ED50 (log mg per administration); back-transform 10.5 mg per administration")  # Supplementary Table S1.1, ixekizumab ED50 = 10.5 mg (CV 15.6%)
    lkdrug_ixekizumab <- log(1.17)
    label("Log ixekizumab drug-effect onset rate (log week^-1); back-transform 1.17 /week")  # Supplementary Table S1.1, ixekizumab kdrug = 1.17 (CV 15.8%)

    emax_secukinumab <- 5.75
    label("Maximum secukinumab effect on the PASI75 logit scale (unitless log-odds)")  # Supplementary Table S1.1, secukinumab Emax = 5.75 (CV 5.01%)
    led50_secukinumab <- log(73.3)
    label("Log secukinumab ED50 (log mg per administration); back-transform 73.3 mg per administration")  # Supplementary Table S1.1, secukinumab ED50 = 73.3 mg (CV 18.5%)
    lkdrug_secukinumab <- log(0.467)
    label("Log secukinumab drug-effect onset rate (log week^-1); back-transform 0.467 /week")  # Supplementary Table S1.1, secukinumab kdrug = 0.467 (CV 14.0%)

    # ---- CD2 antagonist (no ED50 estimable) ----
    emax_alefacept <- 1.87
    label("Maximum alefacept effect on the PASI75 logit scale, realised as a single-step offset at any positive dose (unitless log-odds)")  # Supplementary Table S1.1, alefacept Emax = 1.87 (CV 16.3%); ED50 'NAa'
    lkdrug_alefacept <- log(0.109)
    label("Log alefacept drug-effect onset rate (log week^-1); back-transform 0.109 /week")  # Supplementary Table S1.1, alefacept kdrug = 0.109 (CV 72.3%)

    # ---- PDE4 inhibitor ----
    emax_apremilast <- 9.03
    label("Maximum apremilast effect on the PASI75 logit scale (unitless log-odds)")  # Supplementary Table S1.1, apremilast Emax = 9.03 (CV 55.5%)
    led50_apremilast <- log(96.2)
    label("Log apremilast ED50 (log mg per administration); back-transform 96.2 mg per administration")  # Supplementary Table S1.1, apremilast ED50 = 96.2 mg (CV 85.1%)
    lkdrug_apremilast <- log(0.409)
    label("Log apremilast drug-effect onset rate (log week^-1); back-transform 0.409 /week")  # Supplementary Table S1.1, apremilast kdrug = 0.409 (CV 29.0%)

    # ---- Dihydrofolate reductase inhibitor (no ED50 estimable) ----
    emax_methotrexate <- 2.42
    label("Maximum methotrexate effect on the PASI75 logit scale, realised as a single-step offset at any positive dose (unitless log-odds)")  # Supplementary Table S1.1, methotrexate Emax = 2.42 (CV 4.73%); ED50 'NAa'
    lkdrug_methotrexate <- log(0.180)
    label("Log methotrexate drug-effect onset rate (log week^-1); back-transform 0.180 /week")  # Supplementary Table S1.1, methotrexate kdrug = 0.180 (CV 22.2%)

    # ---- JAK inhibitor ----
    emax_tofacitinib <- 5.04
    label("Maximum tofacitinib effect on the PASI75 logit scale (unitless log-odds)")  # Supplementary Table S1.1, tofacitinib Emax = 5.04 (CV 7.62%)
    led50_tofacitinib <- log(4.37)
    label("Log tofacitinib ED50 (log mg per administration); back-transform 4.37 mg per administration")  # Supplementary Table S1.1, tofacitinib ED50 = 4.37 mg (CV 19.6%)
    lkdrug_tofacitinib <- log(0.466)
    label("Log tofacitinib drug-effect onset rate (log week^-1); back-transform 0.466 /week")  # Supplementary Table S1.1, tofacitinib kdrug = 0.466 (CV 22.2%)

    # ========================================================================
    # PASI90 SCALING FACTORS -- IMPORTED FROM THE COMPANION LANDMARK MODEL.
    #
    # These two values are NOT estimated by the longitudinal analysis. They
    # are Supplementary Table 2 estimates of the LANDMARK model, and Checchio
    # 2017 imports them into this model to perform its external validation:
    # 'the final longitudinal model for PASI75 was used to simulate (N = 1,000)
    # the PASI90 time-course by incorporating the scaling factor between
    # PASI75 and PASI90 obtained from the landmark model' (Methods 'External
    # validation'), and the Figure 5 caption specifies the scaling acts 'on
    # the placebo and drug term'. They are fixed() here because they are
    # inherited from a different fit, exactly as an upstream-model parameter
    # would be.
    # ========================================================================
    i_pbo_pasi90 <- fixed(-1.211)
    label("Additive shift of the placebo component from the PASI75 to the PASI90 scale, imported from the landmark model (paper: I2; unitless log-odds)")
    # Supplementary Table 2, e_o block, row I(Endpoint = 'PASI90') = -1.211
    # (90% CI -1.349, -1.074). Negative: a PASI90 response is a stricter
    # threshold than PASI75, so fewer placebo patients meet it.

    e_drug_pasi90 <- fixed(0.026)
    label("Multiplicative scaling of the drug component from the PASI75 to the PASI90 scale, imported from the landmark model (paper: I6 / em90; unitless)")
    # Supplementary Table 2, scaling block, row em90 = 0.026 (90% CI -0.012,
    # 0.063), entering Equation 8 as Edrug = em * (1 + PASI90 * I6). Its 90%
    # CI includes zero, so the drug term is close to unchanged between the two
    # endpoints and essentially all of the PASI75-to-PASI90 shift is carried
    # by the placebo term above.

    # ========================================================================
    # BETWEEN-STUDY RANDOM EFFECT (Equation 3 eta; Supplementary Table S1.2).
    # This is an aggregate-data BETWEEN-STUDY variance, NOT popPK
    # between-subject variability: one draw describes one whole trial.
    # Methods 'Longitudinal model': 'The between-study variability was
    # described by [eta], having a normal probability distribution with mean 0
    # and variance [omega]^2.' Discussion: 'a placebo component was modeled as
    # a population mean response with between-study variability as a random
    # effect. The estimated random effects in the longitudinal model were
    # nearly normally distributed.'
    #
    # SCALE. The table row is labelled 'omega' while the Methods text names
    # the variance 'omega^2', so the printed 0.0594 could be read either as a
    # variance or as an SD. It is read as the VARIANCE here, for two reasons.
    # (i) The table is headed 'Common Parameter Estimates from NONMEM Output',
    # and a NONMEM $OMEGA estimate IS a variance; relabelling it 'omega'
    # without taking a square root is the ordinary way this drift happens.
    # (ii) The SD reading is physically too tight. Under it, omega = 0.0594
    # moves the placebo plateau B only between 4.56 and 5.13 across a 1-SD
    # band, giving a Week-12 placebo PASI75 range of just 4.0-5.5%, whereas
    # published psoriasis placebo arms genuinely span roughly 1-8% and the
    # source's own companion landmark analysis measures an I-squared of 82%
    # for between-study placebo heterogeneity. The variance reading
    # (omega = 0.244) gives 2.3-9.4%, which matches. Note that ini() takes the
    # VARIANCE, so under this reading the tabulated number is entered
    # unchanged. The alternative is recorded in the vignette Errata; it
    # changes no typical-value prediction, only the simulated between-study
    # spread.
    #
    # PLACEMENT. Equation 3 puts exp(eta) on the whole time-varying placebo
    # term, B * (1 - exp(-kpbo*t)) * exp(eta). Because exp(eta) does not
    # depend on t, that is algebraically identical to a log-normal random
    # effect on B alone, which is how it is written in model() and how the
    # Supplementary Table S1.2 footnote names it ('omega, between-study
    # variability on placebo').
    # ========================================================================
    eta_study_b_pbo ~ 0.0594
    # Supplementary Table S1.2, 'omega' = 0.0594 (CV 31.1%); read as the
    # NONMEM variance, so SD = 0.244 on the log scale of B.

    # ========================================================================
    # RESIDUAL ERROR (Equations 11-13; Supplementary Table S1.2).
    #
    #   (11)  PASI response(%) = Pr(PASI{50,75,90,100}) + W * epsilon
    #   (12)  W                = sqrt(Pr * (1 - Pr) / N)
    #   (13)  PASI75           = Pr(PASI75) + W * (epsilon + epsilon_corr)
    #
    # W is the binomial standard error of an arm proportion, so BOTH residual
    # terms are scaled by the arm size N. Methods: 'Because studies were
    # weighted by inverse of standard errors, larger studies had more
    # influence on estimating the parameter means.'
    #
    # epsilon_corr is realised ONCE PER STUDY ARM and is constant across that
    # arm's timepoints -- Methods: 'due to the longitudinal nature of the data
    # and correlation within timepoint measurements of the same group of
    # patients, a third level of variability ... was also incorporated'. A
    # term constant within an arm is an eta, not an epsilon, in nlmixr2, so it
    # is encoded as eta_study_corr below (same device as
    # Chen_2025_methotrexate_acr50_mbma, whose 'correlation coefficient' row
    # plays the identical role).
    #
    # Unlike most of the library's MBMA extractions, which store the residual
    # UNWEIGHTED and leave the 1/sqrt(N) to downstream code, this model
    # carries N_ARM as a covariate column and applies W inside model(). That
    # is required rather than stylistic: epsilon_corr is a RANDOM EFFECT whose
    # SD varies from arm to arm, and a per-arm random-effect SD cannot be
    # reproduced by rescaling the output of a single rxSolve(). It is exactly
    # the case the N_ARM register entry reserves the column for.
    #
    # Both values are read as NONMEM VARIANCES for the same reasons given for
    # omega above, and both are converted to SDs here because they multiply W
    # rather than entering an ini() variance slot. Under this reading the two
    # residual terms together give a total multiplier on the binomial standard
    # error of sqrt(0.563 + 1.07) = 1.28, i.e. the arm-to-arm scatter is 28%
    # wider than pure binomial sampling -- the mild extra-binomial
    # heterogeneity an aggregate-data meta-analysis should show. (The SD
    # reading gives 1.21, so this choice is nearly immaterial for the TOTAL
    # residual and matters only for how it splits between the two terms.)
    # ========================================================================
    eta_study_corr ~ 1.07
    # Supplementary Table S1.2, 'sigma corr' = 1.07 (CV 9.63%), footnoted
    # 'residual error to account for correlation between time points with
    # longitudinal data'. Read as the NONMEM variance; SD = 1.034 in units of
    # the binomial standard error W, i.e. the arm-level offset is 1.034 * W on
    # the probability scale.

    addSd_prob_pasi75 <- 0.75033326
    label("Multiplier on the binomial standard error W = sqrt(Pr*(1-Pr)/N_ARM) giving the independent per-observation residual SD of the arm PASI75 proportion (unitless)")
    # Supplementary Table S1.2, 'sigma' = 0.563 (CV 23.1%), footnoted
    # 'residual error'. Read as the NONMEM variance, so the SD is
    # sqrt(0.563) = 0.75033326. It is NOT an SD on the probability scale: the
    # per-observation residual SD is this number times W (Equation 11), which
    # model() forms explicitly.
  })

  model({
    # ======================================================================
    # Body weight, Equation 9: WTeffect = (WEIG / 90)^theta. 90 kg is the
    # approximate dataset median (Methods 'Covariate model') and is the weight
    # at which every published prediction is made, so both factors are exactly
    # 1 when reproducing Table 1. An arm with unreported weight should be
    # given WT = 90, which is what the source itself did for the 10-15% of
    # studies that did not report it.
    # ======================================================================
    wtPbo <- (WT / 90)^e_wt_b_pbo
    wtKdrug <- (WT / 90)^e_wt_kdrug

    # ======================================================================
    # Placebo component, Equation 3. The between-study random effect is
    # log-normal on the placebo asymptote B (see the ini() PLACEMENT note).
    # ======================================================================
    bStudy <- b_pbo * wtPbo * exp(eta_study_b_pbo)
    f0 <- a_pbo + bStudy * (1 - exp(-exp(lkpbo) * time))

    # ======================================================================
    # Drug component, Equation 4, summed over the 13 drugs. Each arm supplies
    # a positive dose in exactly one CONMED_<drug>_DOSE column and zero in the
    # rest; a zero dose drives that drug's saturating fraction to exactly
    # zero, so a placebo arm (all columns zero) reduces to f0 alone. Setting
    # two columns positive at once is outside the source's design and would
    # make the model ADD the two drug effects.
    #
    # The per-drug onset rate carries the shared body weight exponent:
    # kdrug_d(WT) = kdrug_d * (WT/90)^-0.834.
    # ======================================================================
    fdAdalimumab <- emax_adalimumab *
      (1 - exp(-exp(lkdrug_adalimumab) * wtKdrug * time)) *
      CONMED_ADALIMUMAB_DOSE / (CONMED_ADALIMUMAB_DOSE + exp(led50_adalimumab))

    # Certolizumab has no estimable onset rate, so its time course is the
    # kdrug -> infinity limit of Equation 4: zero at time 0, full effect
    # thereafter. See the ini() note for the arithmetic that forces this.
    fdCertolizumab <- emax_certolizumab *
      (time > 0) *
      CONMED_CERTOLIZUMAB_DOSE / (CONMED_CERTOLIZUMAB_DOSE + exp(led50_certolizumab))

    fdEtanercept <- emax_etanercept *
      (1 - exp(-exp(lkdrug_etanercept) * wtKdrug * time)) *
      CONMED_ETANERCEPT_DOSE / (CONMED_ETANERCEPT_DOSE + exp(led50_etanercept))

    fdInfliximab <- emax_infliximab *
      (1 - exp(-exp(lkdrug_infliximab) * wtKdrug * time)) *
      CONMED_INFLIXIMAB_DOSE / (CONMED_INFLIXIMAB_DOSE + exp(led50_infliximab))

    fdBriakinumab <- emax_briakinumab *
      (1 - exp(-exp(lkdrug_briakinumab) * wtKdrug * time)) *
      CONMED_BRIAKINUMAB_DOSE / (CONMED_BRIAKINUMAB_DOSE + exp(led50_briakinumab))

    fdUstekinumab <- emax_ustekinumab *
      (1 - exp(-exp(lkdrug_ustekinumab) * wtKdrug * time)) *
      CONMED_USTEKINUMAB_DOSE / (CONMED_USTEKINUMAB_DOSE + exp(led50_ustekinumab))

    fdBrodalumab <- emax_brodalumab *
      (1 - exp(-exp(lkdrug_brodalumab) * wtKdrug * time)) *
      CONMED_BRODALUMAB_DOSE / (CONMED_BRODALUMAB_DOSE + exp(led50_brodalumab))

    fdIxekizumab <- emax_ixekizumab *
      (1 - exp(-exp(lkdrug_ixekizumab) * wtKdrug * time)) *
      CONMED_IXEKIZUMAB_DOSE / (CONMED_IXEKIZUMAB_DOSE + exp(led50_ixekizumab))

    fdSecukinumab <- emax_secukinumab *
      (1 - exp(-exp(lkdrug_secukinumab) * wtKdrug * time)) *
      CONMED_SECUKINUMAB_DOSE / (CONMED_SECUKINUMAB_DOSE + exp(led50_secukinumab))

    fdApremilast <- emax_apremilast *
      (1 - exp(-exp(lkdrug_apremilast) * wtKdrug * time)) *
      CONMED_APREMILAST_DOSE / (CONMED_APREMILAST_DOSE + exp(led50_apremilast))

    # Alefacept and methotrexate have no estimable ED50, so their dose term is
    # a single-step offset: any positive dose gives the full effect. The
    # magnitude corresponds to the Table 1 clinical dose and must not be read
    # as a prediction for any other dose of these two agents.
    fdAlefacept <- emax_alefacept *
      (1 - exp(-exp(lkdrug_alefacept) * wtKdrug * time)) *
      (CONMED_ALEFACEPT_DOSE > 0)

    fdMethotrexate <- emax_methotrexate *
      (1 - exp(-exp(lkdrug_methotrexate) * wtKdrug * time)) *
      (CONMED_MTX_DOSE > 0)

    fdTofacitinib <- emax_tofacitinib *
      (1 - exp(-exp(lkdrug_tofacitinib) * wtKdrug * time)) *
      CONMED_TOFACITINIB_DOSE / (CONMED_TOFACITINIB_DOSE + exp(led50_tofacitinib))

    fdrug <- fdAdalimumab + fdCertolizumab + fdEtanercept + fdInfliximab +
      fdBriakinumab + fdUstekinumab + fdBrodalumab + fdIxekizumab +
      fdSecukinumab + fdApremilast + fdAlefacept + fdMethotrexate +
      fdTofacitinib

    # ======================================================================
    # Equations 1-2: the PASI75 responder probability is the inverse logit of
    # the summed placebo and drug components.
    # ======================================================================
    lp_pasi75 <- f0 + fdrug
    pPasi75 <- expit(lp_pasi75)

    # ======================================================================
    # PASI90 time course -- the source's EXTERNAL VALIDATION (Figure 5). The
    # PASI75 model is carried over unchanged and the two landmark scaling
    # factors are applied, additively on the placebo term and multiplicatively
    # on the drug term. This is a secondary output with no residual of its
    # own: the longitudinal model was FITTED to PASI75 only, and the PASI90
    # trajectory is a prediction that the source then checked against held-out
    # data.
    # ======================================================================
    lp_pasi90 <- (f0 + i_pbo_pasi90) + fdrug * (1 + e_drug_pasi90)
    prob_pasi90 <- expit(lp_pasi90)

    # ======================================================================
    # Residual model, Equations 11-13. W is the binomial standard error of
    # the arm proportion; the arm-level correlated term is added on the
    # probability scale, and the independent term becomes the record's
    # residual SD. Setting the random effects to zero (rxode2::zeroRe())
    # returns the typical-value trajectory pPasi75 exactly, which is what the
    # vignette compares against Table 1.
    # ======================================================================
    wArm <- sqrt(pPasi75 * (1 - pPasi75) / N_ARM)
    sdArm <- wArm * addSd_prob_pasi75

    prob_pasi75 <- pPasi75 + wArm * eta_study_corr
    prob_pasi75 ~ add(sdArm)
  })
}
