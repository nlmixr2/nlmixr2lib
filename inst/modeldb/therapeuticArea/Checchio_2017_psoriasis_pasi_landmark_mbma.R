Checchio_2017_psoriasis_pasi_landmark_mbma <- function() {
  description <- paste0(
    "MBMA. LANDMARK (Week 12) dose-response model-based meta-analysis of the ",
    "PASI50, PASI75, PASI90 and PASI100 responder rates in moderate-to-severe ",
    "plaque psoriasis, fitted jointly across all four endpoints to study-arm ",
    "summary data from 71 trials of systemic agents published 1998-2015 plus ",
    "an internal Pfizer tofacitinib database. The logit of the responder rate ",
    "is a study placebo effect plus a drug effect: the placebo effect is a ",
    "per-trial intercept shifted by a fixed per-endpoint offset and by body ",
    "weight, and the drug effect is a SIGMOIDAL Emax function of dose with ",
    "one Emax per drug CLASS and one ED50 per DRUG, rescaled per endpoint. A ",
    "single Hill coefficient (0.785) is shared by every class except the ",
    "IL-17 inhibitors, whose dose-response is steeper (2.93). Sixteen drug ",
    "arms are covered: adalimumab, certolizumab, etanercept, infliximab, ",
    "brodalumab, ixekizumab, secukinumab, briakinumab, ustekinumab, ",
    "tofacitinib, baricitinib, alefacept, apremilast, methotrexate, ",
    "ciclosporin and acitretin. The traditional oral agents (methotrexate, ",
    "ciclosporin, acitretin) have no dose-response and enter as single-step ",
    "offsets, because per-patient titration left insufficient literature ",
    "dose-response information. Dose enters through one CONMED_<drug>_DOSE ",
    "covariate column per drug rather than through rxode2 dose events; there ",
    "is no PK layer and no time axis, since this is a single Week-12 ",
    "landmark. The model reproduces all 36 published Week-12 responder rates ",
    "of the source Table 2 with a median absolute error of 0.4 percentage ",
    "points. TWO values are NOT printed by the source and are BACK-SOLVED ",
    "from its own Table 2 predictions, flagged inline and in the vignette: ",
    "the typical study placebo intercept and the briakinumab ED50. ",
    "Simulation scope is STUDY-ARM-MEAN responder rates, NOT individual ",
    "patients. The companion time-course model from the same paper is ",
    "modellib('Checchio_2017_psoriasis_pasi75_longitudinal_mbma')."
  )

  reference <- paste(
    "Checchio T, Ahadieh S, Gupta P, Mandema J, Puig L, Wolk R, Valdez H,",
    "Tan H, Krishnaswami S, Tallman A, Kaur M, Ito K.",
    "Quantitative Evaluations of Time-Course and Treatment Effects of Systemic",
    "Agents for Psoriasis: A Model-Based Meta-Analysis.",
    "Clin Pharmacol Ther. 2017;102(6):1006-1016.",
    "doi:10.1002/cpt.732. PMC5697570.",
    "Structural model: Equations 5-8 of Methods 'Landmark model'.",
    "Covariate model: Equation 10 of Methods 'Covariate model'.",
    "Residual model: Equations 11-12 of Methods 'Residual error model'.",
    "All display equations are rasterised in the published PDF and are LOST",
    "from the preprocessed markdown; they were recovered with",
    "pdftotext -layout. Parameter values are Supplementary Table 2 of the",
    "Supplementary Appendix (CPT-102-1006-s001.docx), obtained from the",
    "EuropePMC supplementaryFiles endpoint for PMC5697570.",
    sep = " "
  )

  vignette <- "Checchio_2017_psoriasis_systemic_agents"

  # No PK layer, no concentration, no rxode2 dose events and no time axis:
  # dose is a per-arm covariate column and the read-out is a single Week-12
  # landmark. The placeholder `units` entries follow the arm-level
  # responder-rate MBMA convention already used by
  # Dodds_2013_psoriasis_biologics_mbma and
  # Serrano_2026_atopicDermatitis_placebo_mbma.
  units <- list(
    time = "week (this is a LANDMARK model with no time axis; every prediction is the Week-12 read-out, the primary efficacy timepoint in the majority of contributing studies, range 10-24 weeks)",
    dosing = "mg/administration (dose per administration of the named agent, supplied in the CONMED_<drug>_DOSE covariate columns; infliximab and ciclosporin are per kilogram. This model consumes NO rxode2 dose events.)",
    concentration = "probability/arm (prob_pasi50, prob_pasi75, prob_pasi90 and prob_pasi100 are STUDY-ARM probabilities of achieving the named reduction from baseline in the Psoriasis Area and Severity Index, on a 0-1 scale; none is a drug concentration. The slash satisfies checkModelConventions unit parsing.)"
  )

  covariateData <- list(
    WT = list(
      description = "Study-arm MEAN body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "TRIAL-ARM-LEVEL, not subject-level. Enters ADDITIVELY and CENTERED at 90 kg via Checchio 2017 Equation 10, WTeffect = (WEIG - 90) * I4 -- a different functional form from the companion longitudinal model, which uses the power form (WEIG/90)^theta of Equation 9. In the landmark model the weight effect was retained on the PLACEBO term ONLY: Results 'Landmark model' states 'There was no significant improvement in the model fit when the weight effect was applied to any of the drug effect terms (Emax, ED50), while it was significant for the placebo effect term. Therefore, the presented predictions are based on a model in which body weight was retained for the placebo effect term alone.' The coefficient is negative, so heavier arms show a lower placebo response. Methods 'Covariate model' states that when body weight was not reported (10-15% of studies) WEIG was set to 90 so the term vanishes; downstream code should do the same for an arm with unknown weight. The source cautions (Discussion) that this is an AGGREGATE-level covariate whose magnitude is likely attenuated relative to the individual-level truth.",
      source_name = "WEIG / weightI (Checchio 2017 Equation 10; Supplementary Table 2 row 'I(weightI - 90)')"
    ),
    ROUTE_IV = list(
      description = "1 = the alefacept arm was dosed intravenously, 0 = intramuscularly.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (intramuscular alefacept, the non-IV comparator in this analysis). Note the register's founding contrast is IV-vs-SC; as in Clegg_2024_tixagevimab_cilgavimab.R the non-IV route here is intramuscular, not subcutaneous.",
      notes = "Applies to ALEFACEPT ONLY: it is the sole agent for which Supplementary Table 2 carries a route term, the ED50 row 'alefacetp.I(bio.route = 'IV')' = -1.502 (the source's spelling of the drug name in that row is a typo). It lowers the alefacept log ED50 by 1.502, i.e. an IV arm is about 4.5-fold more potent per milligram than an IM arm, which is the expected direction for a fusion protein whose IM bioavailability is incomplete. ROUTE_IV has NO effect for any other drug in this model and should be left at 0 for them. IMPORTANT: the source's own Table 2 prediction for 'Alefacept 10 mg QW' is reproduced by the IV branch, not the IM one -- ROUTE_IV = 1 gives 24.1% PASI75 against the published 22.7%, whereas ROUTE_IV = 0 gives 11.8%. Set ROUTE_IV = 1 to reproduce any published alefacept number.",
      source_name = "bio.route = 'IV' (Checchio 2017 Supplementary Table 2, ed (ED50) block)"
    ),
    REGI_Q4W = list(
      description = "1 = the arm was dosed once every 4 weeks, 0 = the per-drug comparator interval.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the comparator interval; every 2 weeks for brodalumab, which is the only drug this indicator affects in this model).",
      notes = "Applies to BRODALUMAB ONLY: it is the sole agent for which Supplementary Table 2 carries a dosing-interval term, the ED50 row 'brodalumab.I(bio.freq = 'Q4W')' = 0.151, which raises the brodalumab log ED50 by 0.151 (an ED50 of 104.1 mg on Q4W against 89.5 mg on Q2W, i.e. a Q4W arm needs about 16% more drug per administration for the same effect). Its 90% CI (-0.033, 0.334) includes zero, so the interval effect is weak. REGI_Q4W has NO effect for any other drug in this model and should be left at 0 for them -- in particular, ixekizumab's Table 1 clinical regimen is also Q4W but ixekizumab has a single ED50 with no interval term, so setting REGI_Q4W = 1 for an ixekizumab arm would be wrong. Table 2's published brodalumab prediction is the Q2W arm, so REGI_Q4W = 0 reproduces it. Member of the REGI_* regimen-indicator family; the REGI_QM register entry explicitly anticipates registering a REGI_Q4W sibling.",
      source_name = "bio.freq = 'Q4W' (Checchio 2017 Supplementary Table 2, ed (ED50) block)"
    ),
    CONMED_ADALIMUMAB_DOSE = list(
      description = "Adalimumab dose per administration in the study arm; 0 if the arm did not receive adalimumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose PER ADMINISTRATION, not total daily or weekly dose; Table 2's clinical dose is 40 mg Q2W and Dose = 40 reproduces its published Week-12 PASI75 (65.4% vs 64.9%) and PASI90 (38.1% vs 37.9%). Zero for every arm not randomised to adalimumab, which collapses the adalimumab term to exactly zero. Methods notes doses were normalised to the approved (or most frequently reported) regimen to handle titration schemes.",
      source_name = "DOSE_d (Checchio 2017 Equation 7); 'Clinical dose' column of Table 2"
    ),
    CONMED_CERTOLIZUMAB_DOSE = list(
      description = "Certolizumab pegol dose per administration in the study arm; 0 if the arm did not receive certolizumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration; Table 2's clinical dose is 200 mg Q2W. Investigational for psoriasis at the time of the analysis (Table 2 footnote a).",
      source_name = "DOSE_d (Checchio 2017 Equation 7); 'Clinical dose' column of Table 2"
    ),
    CONMED_ETANERCEPT_DOSE = list(
      description = "Etanercept dose per administration in the study arm; 0 if the arm did not receive etanercept.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration. Table 2 reports 25 mg and 50 mg twice-weekly arms that share one ED50 and differ only in this column; both are reproduced to within 0.7 percentage points, which is the cleanest internal check that the column carries dose per administration rather than a weekly total.",
      source_name = "DOSE_d (Checchio 2017 Equation 7); 'Clinical dose' column of Table 2"
    ),
    CONMED_INFLIXIMAB_DOSE = list(
      description = "Infliximab dose per administration in the study arm, per kilogram of body weight; 0 if the arm did not receive infliximab.",
      units = "mg/kg",
      type = "continuous",
      reference_category = NULL,
      notes = "UNIQUE UNITS: infliximab is dosed on a mg/kg basis and its landmark ED50 of exp(-0.247) = 0.781 is therefore in mg/kg. Table 2's clinical dose is 5 mg/kg Q8W and Dose = 5 reproduces its published Week-12 PASI75 (77.9% vs 77.0%). Supplying a milligram amount would overstate the dose roughly 90-fold and saturate the Emax term.",
      source_name = "DOSE_d (Checchio 2017 Equation 7); 'Clinical dose' column of Table 2"
    ),
    CONMED_BRIAKINUMAB_DOSE = list(
      description = "Briakinumab dose per administration in the study arm; 0 if the arm did not receive briakinumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration; Table 2's clinical dose is 100 mg Q4W. Briakinumab was discontinued in development (Table 2 footnote b). NOTE that briakinumab is the one drug whose landmark ED50 is NOT printed in Supplementary Table 2; the value used here is back-solved from Table 2 and flagged in ini().",
      source_name = "DOSE_d (Checchio 2017 Equation 7); 'Clinical dose' column of Table 2"
    ),
    CONMED_USTEKINUMAB_DOSE = list(
      description = "Ustekinumab dose per administration in the study arm; 0 if the arm did not receive ustekinumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration; Table 2's clinical dose is 45 mg Q12W. The approved ustekinumab regimen is weight-banded (45 mg at or below 100 kg, 90 mg above 100 kg), which the source handles in its Figure 4 forest plot by simulating 500 arms per weight stratum and averaging (Figure 4 caption); that banding is a DOSING rule applied to this column by downstream code, not a model term.",
      source_name = "DOSE_d (Checchio 2017 Equation 7); 'Clinical dose' column of Table 2"
    ),
    CONMED_BRODALUMAB_DOSE = list(
      description = "Brodalumab dose per administration in the study arm; 0 if the arm did not receive brodalumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration; Table 2's clinical dose is 210 mg Q2W (set REGI_Q4W = 0 to match it). Brodalumab and the other IL-17 inhibitors carry the steeper Hill coefficient of 2.93 rather than the 0.785 shared by every other class.",
      source_name = "DOSE_d (Checchio 2017 Equation 7); 'Clinical dose' column of Table 2"
    ),
    CONMED_IXEKIZUMAB_DOSE = list(
      description = "Ixekizumab dose per administration in the study arm; 0 if the arm did not receive ixekizumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration; Table 2's clinical dose is 80 mg Q4W. Ixekizumab has the highest predicted Week-12 placebo-adjusted response of any agent in the analysis (82.7 percentage points for PASI75). Leave REGI_Q4W = 0 even though the clinical regimen is Q4W: the interval term in Supplementary Table 2 is brodalumab-specific.",
      source_name = "DOSE_d (Checchio 2017 Equation 7); 'Clinical dose' column of Table 2"
    ),
    CONMED_SECUKINUMAB_DOSE = list(
      description = "Secukinumab dose per administration in the study arm; 0 if the arm did not receive secukinumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration; Table 2's clinical dose is 150 mg QM.",
      source_name = "DOSE_d (Checchio 2017 Equation 7); 'Clinical dose' column of Table 2"
    ),
    CONMED_ALEFACEPT_DOSE = list(
      description = "Alefacept dose per administration in the study arm; 0 if the arm did not receive alefacept.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration; Table 2's clinical dose is 10 mg QW. Unlike the companion longitudinal model, which could not estimate an alefacept ED50 and used a single-step offset, the landmark model DOES estimate one (two, in fact -- one per route; see ROUTE_IV). Alefacept has since been withdrawn from the market.",
      source_name = "DOSE_d (Checchio 2017 Equation 7); 'Clinical dose' column of Table 2"
    ),
    CONMED_APREMILAST_DOSE = list(
      description = "Apremilast dose per administration in the study arm; 0 if the arm did not receive apremilast.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration, NOT total daily dose: Table 2's clinical dose is 30 mg b.i.d. and Dose = 30 reproduces its published Week-12 PASI75 (25.9% vs 26.8%) whereas Dose = 60 does not.",
      source_name = "DOSE_d (Checchio 2017 Equation 7); 'Clinical dose' column of Table 2"
    ),
    CONMED_TOFACITINIB_DOSE = list(
      description = "Tofacitinib dose per administration in the study arm; 0 if the arm did not receive tofacitinib.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration, NOT total daily dose: Table 2 reports 5 mg and 10 mg b.i.d. arms sharing one ED50, and Dose = 5 and Dose = 10 reproduce both published PASI75 values (34.8% vs 35.2% and 54.2% vs 53.8%) while 10 and 20 do not. Tofacitinib is the only agent whose data came from an internal Pfizer database rather than the literature (Methods 'Database development').",
      source_name = "DOSE_d (Checchio 2017 Equation 7); 'Clinical dose' column of Table 2"
    ),
    CONMED_BARICITINIB_DOSE = list(
      description = "Baricitinib dose per administration in the study arm; 0 if the arm did not receive baricitinib.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration; Table 2's clinical dose is 10 mg QD. Baricitinib appears in the LANDMARK model only: Results 'Available data' states the longitudinal model excluded baricitinib (with ciclosporin and acitretin) for insufficient longitudinal data, which is why the companion longitudinal file carries no baricitinib column. Investigational for psoriasis (Table 2 footnote a); its 90% CI on PASI75 (21.4, 49.0) is the widest of any drug in the table.",
      source_name = "DOSE_d (Checchio 2017 Equation 7); 'Clinical dose' column of Table 2"
    ),
    CONMED_MTX_DOSE = list(
      description = "Methotrexate dose per administration in the study arm; 0 if the arm did not receive methotrexate. Any positive value selects the full methotrexate effect, because the traditional oral agents have no dose-response in this model.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED. Methods 'Landmark model': 'For the traditional oral systemic treatments (methotrexate, cyclosporine, and acitretin), the optimal PASI responses were achieved in practice by flexible dose adjustment for each patient, based on safety and efficacy. Therefore, there were insufficient data in the literature to fully describe the dose-response relationship for these drugs. As such, the maximum responses for traditional oral systemic agents were modeled as single-step offsets (Equation 7).' The magnitude corresponds to Table 2's 18 mg QW clinical dose and must not be read as a prediction for any other methotrexate dose. Register canonical CONMED_MTX_DOSE already exists for randomised methotrexate in immune-mediated-inflammatory-disease dose-response meta-analyses and is reused unchanged.",
      source_name = "DOSE_d (Checchio 2017 Equation 7, traditional-agent branch); 'Clinical dose' column of Table 2"
    ),
    CONMED_CSA_DOSE = list(
      description = "Ciclosporin dose in the study arm; 0 if the arm did not receive ciclosporin. Any positive value selects the full ciclosporin effect, because the traditional oral agents have no dose-response in this model.",
      units = "mg/kg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, for the same reason as methotrexate (see that entry). The magnitude corresponds to Table 2's clinical dose of 2.5-5 mg/kg/day, which the source reports as a RANGE rather than a single value precisely because ciclosporin was titrated per patient. Note the units: the register's CONMED_CSA_DOSE entry describes a total daily dose in mg, whereas this arm-level MBMA records ciclosporin per kilogram per day as the source does; because the column acts only as a step here, the discrepancy changes no prediction, but downstream code must not compare this column numerically against a milligram CONMED_CSA_DOSE from another model. Ciclosporin appears in the LANDMARK model only (Results 'Available data').",
      source_name = "DOSE_d (Checchio 2017 Equation 7, traditional-agent branch); 'Clinical dose' column of Table 2"
    ),
    CONMED_ACITRETIN_DOSE = list(
      description = "Acitretin dose per administration in the study arm; 0 if the arm did not receive acitretin. Any positive value selects the full acitretin effect, because the traditional oral agents have no dose-response in this model.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, for the same reason as methotrexate (see that entry). The magnitude corresponds to Table 2's 30 mg QD clinical dose. Acitretin carries its OWN maximum-effect intercept in Supplementary Table 2 ('vitamin A analog' = 1.783), separate from the 'dmard.(intercept)' shared by methotrexate and ciclosporin, so it is not a dmard offset; see the ini() note. Acitretin appears in the LANDMARK model only (Results 'Available data').",
      source_name = "DOSE_d (Checchio 2017 Equation 7, traditional-agent branch); 'Clinical dose' column of Table 2"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 71L,
    age_range = "not reported at arm level",
    weight_range = "not reported at arm level; 90 kg is stated to be the approximate median of the dataset (Methods 'Landmark model', Equation 6 commentary) and is the typical weight used for every published prediction",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported at arm level",
    disease_state = "adults with moderate to severe plaque psoriasis enrolled in randomised placebo- or active-controlled trials",
    dose_range = "clinical doses per Table 2: adalimumab 40 mg Q2W, certolizumab 200 mg Q2W, etanercept 25 and 50 mg BIW, infliximab 5 mg/kg Q8W, brodalumab 210 mg Q2W, ixekizumab 80 mg Q4W, secukinumab 150 mg QM, briakinumab 100 mg Q4W, ustekinumab 45 mg Q12W, methotrexate 18 mg QW, tofacitinib 5 and 10 mg BID, baricitinib 10 mg QD, acitretin 30 mg QD, alefacept 10 mg QW, apremilast 30 mg BID, ciclosporin 2.5-5 mg/kg/day",
    timepoints = "a SINGLE landmark read-out per arm, at each trial's primary efficacy measurement time; this was 12 weeks in the majority of studies, with a range of 10-24 weeks across all studies (Methods 'Database development')",
    regions = "international; literature 1998-2015 identified by an Ovid Medline / Summary Basis of Approval / European Public Assessment Report search following the Cochrane approach, plus an internal Pfizer tofacitinib database",
    notes = "MBMA at the STUDY-ARM level: each modelled observation is one trial arm's responder rate for one PASI endpoint at the landmark time, so the random effects are BETWEEN-STUDY and this model must not be used to simulate individual patients. n_subjects is NA because the source reports patient counts only graphically (Figure 3 circle areas are proportional to arm N) and never tabulates a pooled total. The literature search yielded 912 abstracts, of which 151 studies were screened in; 71 studies survived into this landmark analysis and 57 into the companion longitudinal analysis (Figure 1 and Results 'Available data'). Between-study heterogeneity in the placebo response is substantial: the source reports an I-squared statistic of 82% for the study-specific placebo effects and states that the 71 fitted placebo effects 'were shown to be normally distributed around a single value (data not shown)' -- that single value is not printed, and the typical intercept used in this file is back-solved from Table 2 (see ini()). Fitted with a nonlinear mixed-effects routine in S-PLUS 8.0.4, NOT in NONMEM; the companion longitudinal model was fitted in NONMEM 7.3."
  )

  ini({
    # ========================================================================
    # PROVENANCE FOR THE WHOLE ini() BLOCK
    #
    # Structural model, Checchio 2017 Methods 'Landmark model', Equations 5-8.
    # Every display equation in this paper is RASTERISED and is dropped by the
    # markdown preprocessor; the forms below were recovered with
    # `pdftotext -layout`:
    #
    #   (5)  P(event)_ijk = g{E0_i + Edrug + eta_i,k}
    #   (6)  E0_i         = trial_i + PASI50*I1 + PASI90*I2 + PASI100*I3
    #   (7)  traditional oral agent:  em = em_d
    #        all other drugs:         em = em_d * DOSE_d^gamma_d /
    #                                        (exp(ED50_d)^gamma_d + DOSE_d^gamma_d)
    #   (8)  Edrug        = em * (1 + PASI50*I5 + PASI90*I6 + PASI100*I7)
    #
    # Covariate model, Equation 10:  WTeffect = (WEIG - 90) * I4, ADDITIVE on
    # the placebo term (contrast the companion longitudinal model, which uses
    # the POWER form of Equation 9).
    #
    # PASI50 / PASI90 / PASI100 in Equations 6 and 8 are ENDPOINT INDICATORS,
    # not covariates: exactly one is 1 for a given observation and all three
    # are 0 for the PASI75 reference endpoint. Rather than carry three
    # indicator columns, this file emits all four endpoint probabilities as
    # separate outputs, which is the same information in a form a downstream
    # user can read directly.
    #
    # All values are Supplementary Table 2, 'Parameter Estimate' column, with
    # the 90% confidence interval from the flanking columns. Per that table's
    # footnote a, 'gamma and ED50 were parameterized as exponential in the
    # model', so the tabulated ed and nt rows are on the LOG scale and are
    # entered here unchanged (NOT re-logged).
    #
    # SOURCE-TRACE CONFIRMATION. Nothing below is tuned except the two
    # explicitly flagged back-solved values. The parameters are transcribed
    # from Supplementary Table 2 and then checked against the 36 published
    # Week-12 responder rates of Table 2 (PASI75 and PASI90 for 18 drug / dose
    # arms), which they reproduce with a median absolute error of 0.4
    # percentage points and a maximum of 4.4; the vignette tabulates every
    # one. The two largest residuals (brodalumab and ustekinumab) are in the
    # direction expected from the fact that Table 2's entries are MEANS over
    # 1,000 simulations that include the between-study random effect, whereas
    # the values here are typical-value predictions.
    # ========================================================================

    # ---- Placebo component (Equation 6; Supplementary Table 2, e_o block) --
    #
    # DERIVED VALUE -- NOT PRINTED BY THE SOURCE.
    #
    # Equation 6 makes the placebo intercept trial-specific: 'A different
    # placebo response was estimated for every trial' (Methods) -- 71 of them.
    # None is printed, and no typical value is printed either; the source only
    # states that the 71 estimates 'were shown to be normally distributed
    # around a single value (data not shown)' (Discussion). A library model
    # needs that single value.
    #
    # It is recovered by INVERTING the source's own published predictions
    # through the source's own published equation. Every other parameter in
    # Equations 5-8 is printed, so Table 2's 36 Week-12 responder rates are a
    # one-parameter family in trial_i; a least-squares fit over all 36 gives
    # -2.8986. The fit is heavily over-determined (36 observations, 1 unknown)
    # and the residuals are small and unstructured, which is the evidence that
    # the inversion is sound rather than a curve-fit: median absolute error
    # 0.4 percentage points across the 36 values. Equivalently and much more
    # simply, Table 2's own 'difference from placebo' columns imply a typical
    # placebo PASI75 of 64.9 - 59.5 = 5.4% and a typical placebo PASI90 of
    # 37.9 - 36.3 = 1.6%, against this value's 5.2% and 1.6%.
    #
    # The vignette re-runs the fit from scratch so a reviewer can audit it.
    # This is a value the source COULD have printed and did not; it is not a
    # value invented from outside the paper.
    e0_trial <- -2.8986
    label("Typical study placebo effect on the PASI75 logit scale at 90 kg (paper: trial_i of Equation 6; unitless log-odds). BACK-SOLVED from the source's Table 2 -- see the note above; the source prints only per-trial values, none of which is tabulated.")

    i_pbo_pasi50 <- 1.133
    label("Additive shift of the placebo effect from the PASI75 to the PASI50 scale (paper: I1; unitless log-odds)")
    # Supplementary Table 2, e_o block, I(Endpoint = 'PASI50') = 1.133
    # (90% CI 1.028, 1.238). POSITIVE: PASI50 is a laxer threshold than
    # PASI75, so more patients meet it. Methods: 'A different placebo response
    # was estimated for every trial with a constant fixed shift between PASI
    # levels' -- these three offsets are that fixed shift.

    i_pbo_pasi90 <- -1.211
    label("Additive shift of the placebo effect from the PASI75 to the PASI90 scale (paper: I2; unitless log-odds)")
    # Supplementary Table 2, e_o block, I(Endpoint = 'PASI90') = -1.211
    # (90% CI -1.349, -1.074). This is the value the companion longitudinal
    # model imports for its external validation of the PASI90 time course.

    i_pbo_pasi100 <- -2.646
    label("Additive shift of the placebo effect from the PASI75 to the PASI100 scale (paper: I3; unitless log-odds)")
    # Supplementary Table 2, e_o block, I(Endpoint = 'PASI100') = -2.646
    # (90% CI -2.907, -2.385). Most negative of the three, as it must be:
    # PASI100 is complete clearance.

    e_wt_e0 <- -0.014
    label("Additive shift in the placebo effect per kilogram of arm mean body weight above 90 kg (paper: I4; log-odds per kg)")
    # Supplementary Table 2, e_o block, I(weightI - 90) = -0.014 (90% CI
    # -0.016, -0.011), entering Equation 10 as WTeffect = (WEIG - 90) * I4.
    # NEGATIVE: heavier arms show a lower placebo response. Small but
    # precisely estimated -- the CI excludes zero and is only 0.005 wide.
    # Retained on the PLACEBO term only (Results 'Landmark model').

    # ========================================================================
    # MAXIMUM DRUG EFFECT (Equation 7 em_d; Supplementary Table 2, em block).
    # One Emax per drug CLASS, expressed as a reference intercept plus a
    # per-class offset, exactly as the source's regression-style row labels
    # do. The reference class carries no offset row.
    #
    # WHICH CLASSES SHARE THE REFERENCE. Supplementary Table 2 prints offsets
    # for only four of the six biologic-side classes: IL-12/23, JAK, CD2 and
    # IL-17. The TNF-alpha inhibitors and the PDE4 inhibitor (apremilast) have
    # no offset row and therefore take bio.(intercept) unchanged. That reading
    # is confirmed numerically rather than assumed: apremilast's Table 2
    # Week-12 PASI75 of 26.8% is reproduced as 25.9% by the bare intercept.
    # ========================================================================
    em_bio <- 5.126
    label("Maximum drug effect for the reference biologic classes -- TNF-alpha inhibitors and the PDE4 inhibitor -- on the PASI75 logit scale (paper: bio.(intercept); unitless log-odds)")
    # Supplementary Table 2, em block, bio.(intercept) = 5.126 (90% CI 4.120,
    # 6.132).

    e_class_il1223_em <- -0.370
    label("Offset on the maximum drug effect for IL-12/23 inhibitors relative to the reference classes (unitless log-odds)")
    # Supplementary Table 2, bio.I(bio.class = 'IL-12/23 inhibitor') = -0.370
    # (90% CI -1.172, 0.432). CI includes zero. Applies to briakinumab and
    # ustekinumab; class Emax = 5.126 - 0.370 = 4.756.

    e_class_jak_em <- 0.830
    label("Offset on the maximum drug effect for JAK inhibitors relative to the reference classes (unitless log-odds)")
    # Supplementary Table 2, bio.I(bio.class = 'JAK inhibitor') = 0.830
    # (90% CI -1.124, 2.784). Very wide CI, spanning zero. Applies to
    # tofacitinib and baricitinib; class Emax = 5.126 + 0.830 = 5.956.

    e_class_cd2_em <- -2.622
    label("Offset on the maximum drug effect for the CD2 antagonist relative to the reference classes (unitless log-odds)")
    # Supplementary Table 2, bio.I(bio.class = 'CD2 antagonist') = -2.622
    # (90% CI -4.096, -1.148). The only class offset whose CI excludes zero,
    # and much the largest: alefacept is far less efficacious than the other
    # injectables. Class Emax = 5.126 - 2.622 = 2.504.

    e_class_il17_em <- -0.201
    label("Offset on the maximum drug effect for IL-17 inhibitors relative to the reference classes (unitless log-odds)")
    # Supplementary Table 2, bio.I(bio.class = 'IL-17 inhibitor') = -0.201
    # (90% CI -1.219, 0.817). Applies to brodalumab, ixekizumab and
    # secukinumab; class Emax = 5.126 - 0.201 = 4.925. Note that the IL-17
    # inhibitors are the MOST efficacious class in the analysis despite the
    # slightly lower Emax: their advantage is carried by the much steeper
    # Hill coefficient below, which puts the clinical doses further up the
    # dose-response curve.

    em_dmard <- 2.326
    label("Maximum drug effect for the reference traditional oral agent, methotrexate, on the PASI75 logit scale (paper: dmard.(intercept); unitless log-odds)")
    # Supplementary Table 2, em block, dmard.(intercept) = 2.326 (90% CI
    # 2.071, 2.580). A single-step offset with no dose-response (Equation 7,
    # traditional-agent branch).

    e_dmard_ciclosporin_em <- 0.453
    label("Offset on the maximum drug effect for ciclosporin relative to methotrexate (unitless log-odds)")
    # Supplementary Table 2, dmard.I(dmard.name = 'ciclosporin') = 0.453
    # (90% CI -0.037, 0.944); CI includes zero. Ciclosporin step offset =
    # 2.326 + 0.453 = 2.779.

    em_vitamina <- 1.783
    label("Maximum drug effect for the vitamin A analog acitretin on the PASI75 logit scale (unitless log-odds)")
    # Supplementary Table 2, em block, 'vitamin A analog' = 1.783 (90% CI
    # 1.069, 2.497). This row is a SEPARATE INTERCEPT, not an offset on
    # dmard.(intercept): it is printed at the same indentation level as
    # dmard.(intercept) and bio.(intercept), and the arithmetic confirms it --
    # 1.783 alone reproduces Table 2's acitretin PASI75 of 25.0% as 24.7%,
    # whereas 2.326 + 1.783 would give 77%.

    # ========================================================================
    # HILL COEFFICIENT (Equation 7 gamma_d; Supplementary Table 2, nt block).
    # Per that table's footnote a, gamma was 'parameterized as exponential in
    # the model', so these rows are LOG gamma and are entered unchanged.
    # Methods: 'a single Hill coefficient was applied to all drug classes
    # except for IL-17 inhibitors. The dose-response curve for IL-17
    # inhibitors was sufficiently different from other drugs of interest that
    # an individual offset was required'.
    # ========================================================================
    lgamma <- -0.242
    label("Log Hill coefficient shared by every drug class except the IL-17 inhibitors (paper: nt.(intercept); back-transform gamma = 0.785)")
    # Supplementary Table 2, nt.(intercept) = -0.242 (90% CI -0.550, 0.065).
    # gamma = exp(-0.242) = 0.785, i.e. SHALLOWER than a plain Emax curve.

    e_class_il17_lgamma <- 1.318
    label("Offset on the log Hill coefficient for IL-17 inhibitors (unitless); gives gamma = 2.93 for that class")
    # Supplementary Table 2, nt.I(bio.class = 'IL-17 inhibitor') = 1.318
    # (90% CI 0.932, 1.705); CI excludes zero by a wide margin, which is the
    # statistical basis for the Methods statement that the IL-17 curve is
    # 'sufficiently different' to need its own coefficient. IL-17 gamma =
    # exp(-0.242 + 1.318) = exp(1.076) = 2.93, a markedly steeper curve.

    # ========================================================================
    # ED50 (Equation 7; Supplementary Table 2, ed block). Per footnote a these
    # are LOG-SCALE values and Equation 7 applies exp() to them, so they are
    # entered here exactly as printed. Units are mg per administration, EXCEPT
    # infliximab (mg/kg). The traditional oral agents have no ED50 row because
    # they have no dose-response.
    # ========================================================================
    led50_adalimumab <- 2.673
    label("Log adalimumab ED50 (log mg per administration); back-transform 14.5 mg per administration")  # Supplementary Table 2, adalimumab = 2.673 (90% CI 2.033, 3.313)
    led50_certolizumab <- 3.773
    label("Log certolizumab ED50 (log mg per administration); back-transform 43.5 mg per administration")  # Supplementary Table 2, certolizumab = 3.773 (90% CI 2.364, 5.181)
    led50_etanercept <- 3.428
    label("Log etanercept ED50 (log mg per administration); back-transform 30.8 mg per administration")  # Supplementary Table 2, etanercept = 3.428 (90% CI 2.907, 3.949)
    led50_infliximab <- -0.247
    label("Log infliximab ED50; back-transform 0.781 mg/kg per administration (NOT mg)")  # Supplementary Table 2, infliximab = -0.247 (90% CI -1.150, 0.656)
    led50_ustekinumab <- 2.292
    label("Log ustekinumab ED50 (log mg per administration); back-transform 9.90 mg per administration")  # Supplementary Table 2, ustekinumab = 2.292 (90% CI 1.958, 2.626)
    led50_ixekizumab <- 2.334
    label("Log ixekizumab ED50 (log mg per administration); back-transform 10.3 mg per administration")  # Supplementary Table 2, ixekizumab = 2.334 (90% CI 2.087, 2.581)
    led50_secukinumab <- 4.497
    label("Log secukinumab ED50 (log mg per administration); back-transform 89.7 mg per administration")  # Supplementary Table 2, secukinumab = 4.497 (90% CI 4.378, 4.615)
    led50_apremilast <- 4.131
    label("Log apremilast ED50 (log mg per administration); back-transform 62.2 mg per administration")  # Supplementary Table 2, apremilast = 4.131 (90% CI 3.486, 4.775)
    led50_tofacitinib <- 2.227
    label("Log tofacitinib ED50 (log mg per administration); back-transform 9.27 mg per administration")  # Supplementary Table 2, tofacitinib = 2.227 (90% CI 1.135, 3.319)
    led50_baricitinib <- 3.017
    label("Log baricitinib ED50 (log mg per administration); back-transform 20.4 mg per administration")  # Supplementary Table 2, baricitinib = 3.017 (90% CI 1.791, 4.242)

    led50_alefacept <- 2.732
    label("Log alefacept ED50 for the reference INTRAMUSCULAR route (log mg per administration); back-transform 15.4 mg per administration")  # Supplementary Table 2, alefacept.(intercept) = 2.732 (90% CI 1.100, 4.364)
    e_route_iv_led50_alefacept <- -1.502
    label("Offset on the log alefacept ED50 for intravenous administration (unitless); gives 3.42 mg IV")
    # Supplementary Table 2, alefacetp.I(bio.route = 'IV') = -1.502 (90% CI
    # -2.896, -0.108); the drug name is misspelled in the source's row label.
    # NEGATIVE, i.e. IV is about 4.5-fold more potent per milligram than IM,
    # the expected direction for incomplete IM bioavailability. Table 2's own
    # alefacept prediction is reproduced by the IV branch (ROUTE_IV = 1), not
    # the IM one -- see the ROUTE_IV covariateData note.

    led50_brodalumab <- 4.494
    label("Log brodalumab ED50 for the reference every-2-weeks interval (log mg per administration); back-transform 89.5 mg per administration")  # Supplementary Table 2, brodalumab.(intercept) = 4.494 (90% CI 4.407, 4.581)
    e_regi_q4w_led50_brodalumab <- 0.151
    label("Offset on the log brodalumab ED50 for every-4-weeks dosing (unitless); gives 104.1 mg on Q4W")
    # Supplementary Table 2, brodalumab.I(bio.freq = 'Q4W') = 0.151 (90% CI
    # -0.033, 0.334); CI includes zero, so the interval effect is weak.
    # Applies to BRODALUMAB ONLY -- see the REGI_Q4W covariateData note.

    # ------------------------------------------------------------------------
    # DERIVED VALUE -- NOT PRINTED BY THE SOURCE.
    #
    # Briakinumab is predicted in Table 2 (100 mg Q4W: PASI75 80.8%, PASI90
    # 58.4%) and is described by Equation 7's sigmoidal branch like every other
    # biologic, but Supplementary Table 2's ed block has NO briakinumab row.
    # The omission is in the source, not in the conversion: the raw
    # WordprocessingML of the supplement was enumerated row by row and the ed
    # block contains 14 rows, none of them briakinumab.
    #
    # Its ED50 is recovered the same way as e0_trial above, by inverting the
    # source's own Table 2 predictions through its own Equation 7 with every
    # other parameter fixed at its printed value. Briakinumab's two published
    # rates are then reproduced exactly (0.0 percentage points on both), which
    # is unsurprising since two observations determine one unknown, and is why
    # the supporting evidence for this value is the SOUNDNESS OF THE JOINT FIT
    # rather than its own residual: e0_trial and this ED50 were fitted
    # together over all 36 Table 2 values, and the other 34 residuals -- which
    # this parameter cannot influence -- have a median absolute error of 0.4
    # percentage points. The vignette re-runs the fit so a reviewer can audit
    # it.
    #
    # Sanity check on the value: 5.11 mg is close to ustekinumab's longitudinal
    # ED50 of 5.05 mg and to its landmark ED50 of 9.90 mg, i.e. it lands where
    # the other IL-12/23 inhibitor does, and the resulting briakinumab clinical
    # dose of 100 mg sits high on the saturating curve, consistent with its
    # being one of the most efficacious agents in Table 2.
    # ------------------------------------------------------------------------
    led50_briakinumab <- 1.6318
    label("Log briakinumab ED50 (log mg per administration); back-transform 5.11 mg per administration. BACK-SOLVED from the source's Table 2 -- Supplementary Table 2 prints no briakinumab ED50 row; see the note above.")

    # ========================================================================
    # ENDPOINT SCALING OF THE DRUG EFFECT (Equation 8; Supplementary Table 2,
    # scaling block). Equation 8 rescales the drug term per endpoint:
    # Edrug = em * (1 + PASI50*I5 + PASI90*I6 + PASI100*I7). All three are
    # small and two have CIs including zero, i.e. the drug term is nearly
    # unchanged across endpoints and essentially the whole PASI75-to-PASI90
    # shift is carried by the placebo offsets above. Methods: 'As with E0, a
    # scaling factor was applied to the estimated Emax to account for
    # differences in PASI75, PASI50, PASI90, and PASI100 scores.'
    # ========================================================================
    e_drug_pasi50 <- -0.040
    label("Multiplicative scaling of the drug effect from the PASI75 to the PASI50 endpoint (paper: I5 / em50; unitless)")  # Supplementary Table 2, em50 = -0.040 (90% CI -0.076, -0.004)
    e_drug_pasi90 <- 0.026
    label("Multiplicative scaling of the drug effect from the PASI75 to the PASI90 endpoint (paper: I6 / em90; unitless)")  # Supplementary Table 2, em90 = 0.026 (90% CI -0.012, 0.063)
    e_drug_pasi100 <- 0.075
    label("Multiplicative scaling of the drug effect from the PASI75 to the PASI100 endpoint (paper: I7 / em100; unitless)")  # Supplementary Table 2, em100 = 0.075 (90% CI 0.009, 0.140)

    # ========================================================================
    # RANDOM EFFECTS AND RESIDUAL ERROR -- NOT REPORTED BY THE SOURCE.
    #
    # Equation 5 carries a random effect eta_i,k on the overall drug response
    # with an endpoint-specific variance omega^2_k, and Methods adds that 'The
    # model included terms to account for the correlation between different
    # PASI scores ... within the same study and to account for the correlation
    # between PASI scores among different arms within the same study.'
    # Equations 11-12 give the same arm-size-weighted residual as the
    # companion longitudinal model.
    #
    # NONE of these variances is printed. Supplementary Table 2 tabulates only
    # the fixed effects; unlike Supplementary Table S1.2 for the longitudinal
    # model, it has no omega or sigma rows, and no correlation matrix appears
    # anywhere in the paper or supplement. No variance is therefore encoded
    # here -- inventing one would be unauditable. The consequence is that this
    # model produces TYPICAL-VALUE arm predictions only; the companion
    # longitudinal model is the one to use when between-study spread is
    # needed. This gap is recorded in the vignette Errata.
    #
    # The placeholder below exists only so the nlmixr2 likelihood machinery
    # accepts the model for forward simulation, the same device used by
    # Serrano_2026_atopicDermatitis_placebo_mbma and
    # Bhatnagar_2024_upadacitinib_asas20_as.
    # ========================================================================
    addSd_prob_pasi75 <- fixed(0.001)
    label("Placeholder additive residual SD on the arm-level probability output prob_pasi75 (unitless); NOT from the source, which reports no residual variance for the landmark model")
  })

  model({
    # ======================================================================
    # Placebo component, Equations 6 and 10. The trial intercept is the
    # typical value across the 71 fitted trials (back-solved; see ini()), and
    # body weight shifts it additively, centered at 90 kg. The three endpoint
    # offsets are applied below, one per output.
    # ======================================================================
    e0Pasi75 <- e0_trial + e_wt_e0 * (WT - 90)

    # ======================================================================
    # Class-level maximum drug effect, Equation 7. The reference intercept
    # em_bio covers the TNF-alpha inhibitors and the PDE4 inhibitor; the four
    # printed class offsets cover the rest. Traditional oral agents have their
    # own intercepts and no dose-response.
    # ======================================================================
    emTnf <- em_bio
    emPde4 <- em_bio
    emIl1223 <- em_bio + e_class_il1223_em
    emJak <- em_bio + e_class_jak_em
    emCd2 <- em_bio + e_class_cd2_em
    emIl17 <- em_bio + e_class_il17_em

    # ======================================================================
    # Hill coefficients, Equation 7 with Supplementary Table 2 footnote a
    # (gamma parameterised as exponential). One shared value, plus a steeper
    # one for the IL-17 inhibitors.
    # ======================================================================
    gam <- exp(lgamma)
    gamIl17 <- exp(lgamma + e_class_il17_lgamma)

    # ======================================================================
    # Per-drug saturating dose fractions, Equation 7:
    #   DOSE^gamma / (ED50^gamma + DOSE^gamma)
    # Each arm supplies a positive dose in exactly one CONMED_<drug>_DOSE
    # column and zero in the rest; a zero dose drives that drug's fraction to
    # exactly zero, so a placebo arm (all columns zero) reduces to the placebo
    # term alone. Setting two columns positive at once is outside the source's
    # design and would make the model ADD the two drug effects.
    #
    # Alefacept's ED50 depends on route and brodalumab's on dosing interval;
    # both indicators are zero for every other drug and have no effect there.
    # ======================================================================
    fAdalimumab <- CONMED_ADALIMUMAB_DOSE^gam /
      (exp(led50_adalimumab)^gam + CONMED_ADALIMUMAB_DOSE^gam)
    fCertolizumab <- CONMED_CERTOLIZUMAB_DOSE^gam /
      (exp(led50_certolizumab)^gam + CONMED_CERTOLIZUMAB_DOSE^gam)
    fEtanercept <- CONMED_ETANERCEPT_DOSE^gam /
      (exp(led50_etanercept)^gam + CONMED_ETANERCEPT_DOSE^gam)
    fInfliximab <- CONMED_INFLIXIMAB_DOSE^gam /
      (exp(led50_infliximab)^gam + CONMED_INFLIXIMAB_DOSE^gam)

    fBriakinumab <- CONMED_BRIAKINUMAB_DOSE^gam /
      (exp(led50_briakinumab)^gam + CONMED_BRIAKINUMAB_DOSE^gam)
    fUstekinumab <- CONMED_USTEKINUMAB_DOSE^gam /
      (exp(led50_ustekinumab)^gam + CONMED_USTEKINUMAB_DOSE^gam)

    ed50Brodalumab <- exp(led50_brodalumab + e_regi_q4w_led50_brodalumab * REGI_Q4W)
    fBrodalumab <- CONMED_BRODALUMAB_DOSE^gamIl17 /
      (ed50Brodalumab^gamIl17 + CONMED_BRODALUMAB_DOSE^gamIl17)
    fIxekizumab <- CONMED_IXEKIZUMAB_DOSE^gamIl17 /
      (exp(led50_ixekizumab)^gamIl17 + CONMED_IXEKIZUMAB_DOSE^gamIl17)
    fSecukinumab <- CONMED_SECUKINUMAB_DOSE^gamIl17 /
      (exp(led50_secukinumab)^gamIl17 + CONMED_SECUKINUMAB_DOSE^gamIl17)

    fTofacitinib <- CONMED_TOFACITINIB_DOSE^gam /
      (exp(led50_tofacitinib)^gam + CONMED_TOFACITINIB_DOSE^gam)
    fBaricitinib <- CONMED_BARICITINIB_DOSE^gam /
      (exp(led50_baricitinib)^gam + CONMED_BARICITINIB_DOSE^gam)

    ed50Alefacept <- exp(led50_alefacept + e_route_iv_led50_alefacept * ROUTE_IV)
    fAlefacept <- CONMED_ALEFACEPT_DOSE^gam /
      (ed50Alefacept^gam + CONMED_ALEFACEPT_DOSE^gam)

    fApremilast <- CONMED_APREMILAST_DOSE^gam /
      (exp(led50_apremilast)^gam + CONMED_APREMILAST_DOSE^gam)

    # ======================================================================
    # Drug effect on the PASI75 reference endpoint, Equation 7. The
    # traditional oral agents enter as single-step offsets triggered by any
    # positive dose, per the traditional-agent branch of Equation 7.
    # ======================================================================
    emArm <-
      emTnf * (fAdalimumab + fCertolizumab + fEtanercept + fInfliximab) +
      emIl1223 * (fBriakinumab + fUstekinumab) +
      emIl17 * (fBrodalumab + fIxekizumab + fSecukinumab) +
      emJak * (fTofacitinib + fBaricitinib) +
      emCd2 * fAlefacept +
      emPde4 * fApremilast +
      em_dmard * (CONMED_MTX_DOSE > 0) +
      (em_dmard + e_dmard_ciclosporin_em) * (CONMED_CSA_DOSE > 0) +
      em_vitamina * (CONMED_ACITRETIN_DOSE > 0)

    # ======================================================================
    # The four endpoint probabilities, Equations 5, 6 and 8. The placebo term
    # takes an ADDITIVE per-endpoint offset and the drug term a
    # MULTIPLICATIVE per-endpoint scaling; PASI75 is the reference, where both
    # corrections vanish.
    # ======================================================================
    prob_pasi50 <- expit(e0Pasi75 + i_pbo_pasi50 + emArm * (1 + e_drug_pasi50))
    prob_pasi90 <- expit(e0Pasi75 + i_pbo_pasi90 + emArm * (1 + e_drug_pasi90))
    prob_pasi100 <- expit(e0Pasi75 + i_pbo_pasi100 + emArm * (1 + e_drug_pasi100))

    prob_pasi75 <- expit(e0Pasi75 + emArm)
    prob_pasi75 ~ add(addSd_prob_pasi75)
  })
}
