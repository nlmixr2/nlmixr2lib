Volkova_2023_lipidLowering_mace_mbma <- function() {
  description <- "MBMA. Random-effects meta-regression of the treatment effect (log risk ratio vs comparator) on ten individual major-adverse-cardiac-event (MACE) components for statin and anti-PCSK9 lipid-lowering therapy, plus the companion trial sample-size calculator used to optimise MACE endpoint composition. Fit with metafor to 54 randomised controlled trials (270,471 patients; 47 statin trials, 7 anti-PCSK9 trials) identified by a PRISMA systematic search. Predicts log(RR) for coronary revascularization (CR), cardiovascular mortality (CVM), heart failure (HF), ischemic stroke (iST), myocardial infarction (MI), stroke (ST), total mortality (TM), coronary mortality (CM), fatal stroke (fST) and nonfatal myocardial infarction (nfMI) from therapy class, comparator-adjusted LDL-C lowering, baseline HDL-C and remnant cholesterol, hypertension prevalence, renal-disease enrolment, prevention category and mean age. The n_required output implements the paper's Supplementary File 2 sample-size formula for a single endpoint or, when supplied the Supplementary File 3 aggregated inputs, for a composite MACE endpoint. This is a STUDY-LEVEL model: predictions are trial-arm mean treatment effects across a trial population, never individual-patient risks or concentrations. Only the ten endpoints for which the paper retained a covariate are implemented; the five covariate-free endpoints (fMI, hST, nfST, TIA, UA) are reported graphically only and are out of scope (see vignette Errata)."

  reference <- paste(
    "Volkova A, Shulgin B, Helmlinger G, Peskov K, Sokolov V.",
    "Optimization of the MACE endpoint composition to increase power in studies",
    "of lipid-lowering therapies - a model-based meta-analysis.",
    "Front Cardiovasc Med. 2023;10:1242845 (published online 08 January 2024).",
    "doi:10.3389/fcvm.2023.1242845.",
    sep = " "
  )

  vignette <- "Volkova_2023_lipidLowering_mace_mbma"

  # Algebraic study-level MBMA: no rxode2 dose events are consumed (every input
  # arrives as a trial-level covariate column) and the outputs are dimensionless
  # risk ratios plus a subject count, not drug concentrations. The `dosing` /
  # `concentration` strings below are placeholders chosen so that
  # checkModelConventions() sees a dimensionally consistent pair; a more
  # descriptive label trips the dimensional-compatibility check. Same device as
  # Yoshioka_2018_FXa_inhibitors_mbma.
  units <- list(
    time          = "year",
    dosing        = "probability",
    concentration = "probability/probability"
  )

  covariateData <- list(
    TRT_PCSK9I = list(
      description        = "Indicator that the trial's active arm was an anti-PCSK9 therapy (monoclonal antibody or siRNA) rather than a statin.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (statin trial -- the reference therapy class for every beta_1 in Volkova 2023 Table 1)",
      notes              = "MBMA trial-level covariate. 1 = anti-PCSK9 trial (evolocumab, alirocumab, bococizumab or inclisiran), 0 = statin trial. Volkova 2023 Supplementary Table 3 gives 7 anti-PCSK9 and 47 statin trials. The Table 1 footnote defines beta_0 as 'the impact of statins therapy' and beta_1 as 'the impact of PCSK9 inhibitors relative to statins', so no sign flip is needed: setting TRT_PCSK9I = 1 adds beta_1 to the statin intercept.",
      source_name        = "PCSK9 inhibitors (Volkova 2023 Table 1 beta_1 column); 'Therapy' in the Supplementary Table 4 covariate search"
    ),
    LDLC_DELTA = list(
      description        = "Trial-level comparator-adjusted change in LDL-C from baseline to the latest available time point. Negative values indicate LDL-C lowering relative to the comparator arm.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "MBMA trial-level covariate, defined in Volkova 2023 Section 2.2 as delta-LDLc = (LDLc - LDLc_baseline)_treatment - (LDLc - LDLc_baseline)_comparator, so it is a placebo/comparator-adjusted absolute change rather than a raw change from baseline. Volkova 2023 Supplementary Table 3 reports mean -36.91 mg/dL (SD 16.57; range -75 to +1.4) across the 54 trials. The unit is load-bearing: the paper's headline statement 'a 1 mmol/L (38.67 mg/dL) decrease in LDLc was associated with RR changes of -24% for CR, -14% for CM and -27% for nfMI' reproduces exactly from the mg/dL coefficients (exp(0.007 * -38.67) = 0.763; exp(0.004 * -38.67) = 0.857; exp(0.008 * -38.67) = 0.734). Do NOT supply mmol/L. Enters the CR, CM and nfMI models only. Distinct from the baseline LDLC canonical, which carries the pre-treatment concentration.",
      source_name        = "delta-LDLc (Volkova 2023 Table 1 / Section 2.2)"
    ),
    HDLC = list(
      description        = "Trial-level mean baseline (pre-treatment) HDL-C concentration.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "MBMA trial-level covariate. Volkova 2023 Supplementary Table 3 reports mean 45.82 mg/dL (SD 6.05; range 35.96-59.55) across 53 trials with the value available. NOTE that the Supplementary Table 3 subsection headings are transposed: the block headed 'Treatment-related lipid measurements' holds the BASELINE values (its LDLc mean of 136.57 mg/dL is the same number the main-text Results calls 'mean baseline LDLc level of 136.6 mg/dl'), while the block headed 'Baseline lipid measurements' holds the comparator-adjusted delta values. The value used here is the baseline one. Enters the CVM, MI and TM models, uncentred, with a linear coefficient on log(RR). Volkova 2023 Results notes the model implies no significant TM/CVM benefit above a baseline HDLc of 48 mg/dL for statins.",
      source_name        = "HDLc (Volkova 2023 Table 1 beta_2 column)"
    ),
    REMC = list(
      description        = "Trial-level mean baseline (pre-treatment) remnant cholesterol concentration.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "MBMA trial-level covariate. Remnant cholesterol is the cholesterol carried by triglyceride-rich lipoproteins; Volkova 2023 Section 2.2 derives it, when not reported, from the Friedewald relations remC = TC - HDLc - LDLc with TC = LDLc + HDLc + TG/5. Supplementary Table 3 reports baseline mean 30.72 mg/dL (SD 7.06; range 13.05-60) across 53 trials (see the transposed-heading note on HDLC). Enters the HF model only, uncentred, with a linear coefficient on log(RR); higher baseline remnant cholesterol predicts greater heart-failure benefit.",
      source_name        = "remC (Volkova 2023 Table 1 beta_2 column)"
    ),
    DIS_HYPERT_PERCENT = list(
      description        = "Percentage of the trial population with hypertension at baseline; range 0-100.",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = "MBMA trial-level cohort-prevalence covariate, scaled in percent (NOT a fraction) -- the percent scale is confirmed by the paper's own Results: 'RR reduction of iST in a population with 75% subjects experiencing high blood pressure was 1.7 [95% CI 1.3-2.3] fold less effective compared to a population with 25% hypertensive patients', which reproduces as exp(0.011 * 50) = 1.73 with CI bounds exp(50 * (0.011 -/+ 1.96 * 0.003)) = 1.29 and 2.33. Volkova 2023 Supplementary Table 3 reports mean 54.71% (SD 22.61; range 15.7-100) across 52 trials. Enters the iST model only. Family precedent: DIS_CHD_PERCENT (Vargo 2014).",
      source_name        = "Hypertension (Volkova 2023 Table 1 beta_2 column); 'Hypertension (%)' in Supplementary Table 3"
    ),
    RENALIMP_SEV = list(
      description        = "Indicator that the trial enrolled patients with severe (late-stage) renal disease.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (trial enrolled patients with normal renal function only)",
      notes              = "MBMA trial-level ELIGIBILITY flag, not an individual-patient renal stratum: 1 = the trial's inclusion criteria admitted patients with severe renal disease, 0 = they did not. Volkova 2023 Supplementary Table 3 records 6 such trials out of 54. The paper does not state the eGFR / CrCL threshold used to classify a trial. This is the same arm-level reading of a cohort-composition flag that Vargo_2014_statins_ezetimibe_mbma applies to DIS_ACS and DIS_HEFH. Enters the fST model as a main effect (exp(0.654) = 1.92-fold less treatment benefit, matching the paper's '1.9 [95% CI 1.1-3.3] times' statement) and the nfMI model ONLY as an interaction with LDLC_DELTA, with no nfMI main effect.",
      source_name        = "RD (Volkova 2023 Table 1 beta_2 / beta_3 columns and footnote); 'Inclusion of patients with renal disease' in Supplementary Table 3"
    ),
    PREVENT_PRIMARY = list(
      description        = "Indicator that the trial enrolled a primary-prevention population (no previous history of a major adverse cardiac event).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not a primary-prevention trial)",
      notes              = "MBMA trial-level covariate; one of the mutually exclusive PREVENT_PRIMARY / PREVENT_SECONDARY / PREVENT_MIXED indicator trio that decomposes the paper's three-level 'Prevention' category. Exactly one of the three must be 1. Volkova 2023 Supplementary Table 3: 15 primary, 21 secondary, 18 mixed trials. Volkova 2023 Section 2.2 assigns the mixed category when the publication does not state the prevention type explicitly. Used ONLY in the ST (unspecified stroke) model, and only through its interaction with AGE -- the paper retained no main effect for either prevention category or age.",
      source_name        = "Prevention = Primary (Volkova 2023 Table 1 beta_2 column)"
    ),
    PREVENT_SECONDARY = list(
      description        = "Indicator that the trial enrolled a secondary-prevention population (previous history of a major adverse cardiac event).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not a secondary-prevention trial)",
      notes              = "MBMA trial-level covariate; see PREVENT_PRIMARY for the indicator-trio convention. Volkova 2023 Supplementary Table 3 records 21 secondary-prevention trials of 54. Used ONLY in the ST model through the interaction with AGE.",
      source_name        = "Prevention = Secondary (Volkova 2023 Table 1 beta_2 column)"
    ),
    PREVENT_MIXED = list(
      description        = "Indicator that the trial enrolled a mixed primary-and-secondary prevention population, or did not state the prevention type.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not a mixed-prevention trial)",
      notes              = "MBMA trial-level covariate; see PREVENT_PRIMARY for the indicator-trio convention. Volkova 2023 Table 1 labels this level 'Both'; Section 2.2 states that a trial was assigned to it when the prevention type was not explicitly indicated in the published source, so it pools genuinely mixed cohorts with unreported ones. Volkova 2023 Supplementary Table 3 records 18 such trials of 54. Used ONLY in the ST model through the interaction with AGE.",
      source_name        = "Prevention = Both (Volkova 2023 Table 1 beta_2 column)"
    ),
    AGE = list(
      description        = "Trial-level mean age of the enrolled population.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "MBMA trial-level covariate, entered UNCENTRED (no reference age is subtracted). Volkova 2023 Supplementary Table 3 reports mean 61.2 years (SD 4.9; range 49.75-75.35) across all 54 trials. Used ONLY in the ST model, and only multiplied by the prevention-category indicator. The uncentred form is confirmed by the paper's Results, which reproduce as RR = 0.70 (statin, primary, age 60) and RR = 0.84 (statin, secondary, age 60) from int_st + e_age_st_* * 60 with no centring term.",
      source_name        = "Age (Volkova 2023 Table 1 beta_2 column); 'Age, years' in Supplementary Table 3"
    ),
    LOGRR_TARGET = list(
      description        = "The log risk ratio the planned trial must show to be statistically significant; input to the sample-size calculator only.",
      units              = "(unitless log ratio)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Sample-size calculator input (symbol theta_i, equivalently u_i, in Volkova 2023 Supplementary File 2). Supply EITHER one of this model's own logrr_* outputs, for sizing a single endpoint, OR the Supplementary File 3 composite aggregate theta_sum = log(sum(k_i * exp(theta_i)) / sum(k_i)) when sizing a composite MACE endpoint. Must be strictly negative for a finite, meaningful n_required: the formula divides by LOGRR_TARGET^2 and the paper states that as the effect size approaches zero the required sample size tends to infinity. Because the effect size is itself a model output, reproducing the paper's sample-size figures is a deliberate two-pass exercise (solve for logrr_*, then re-solve supplying it here) -- the same two-stage structure the paper uses.",
      source_name        = "theta_i / u_i (Volkova 2023 Section 2.4 and Supplementary File 2)"
    ),
    K_CTRL = list(
      description        = "Proportion of control-arm patients who experience the endpoint being sized over the planned trial follow-up; a fraction in (0, 1).",
      units              = "(fraction)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Sample-size calculator input (symbol k_i(t) in Volkova 2023 Section 2.4 and Supplementary File 2). A FRACTION, not a percentage. For a composite endpoint supply the Supplementary File 3 aggregate k_sum = sum(k_i), which assumes the number of patients experiencing several component events within one trial is negligible. Volkova 2023 derived k_i(t) from a separate meta-regression of control-arm event rates against follow-up duration, whose fitted lines are reported ONLY in Supplementary Figure 1 -- no coefficients are printed anywhere in the paper, so this model takes k as an external input rather than computing it from trial duration. A value can be recovered for individual endpoints by inverting the paper's own published sample sizes; the vignette does this for nfMI (k ~ 0.043 at a 4-year follow-up, recovered consistently from two different Figure 3 points). Must be > 0.",
      source_name        = "k_i(t) (Volkova 2023 Section 2.4, Supplementary File 2 and Supplementary Figure 1)"
    )
  )

  population <- list(
    species             = "human",
    n_subjects          = 270471L,
    n_studies           = 54L,
    age_range           = "Trial-level mean age 49.75-75.35 years (mean of means 61.2, SD 4.9) across the 54 trials (Volkova 2023 Supplementary Table 3).",
    weight_range        = "Body weight was not collected; body mass index was, with trial-level means 23.5-30.15 kg/m^2 (mean 26.89, SD 1.75) across 48 trials (Volkova 2023 Supplementary Table 3). BMI was screened but not retained in any final model.",
    sex_female_pct      = 30.3,
    disease_state       = "Dyslipidemia / hypercholesterolemia, treated with a statin or an anti-PCSK9 agent against a placebo or standard-of-care comparator. 15 trials enrolled primary-prevention populations (no previous MACE history), 21 were secondary-prevention studies, and 18 enrolled a mixed or unstated population; 6 of 54 trials admitted patients with severe renal disease.",
    dose_range          = "Trial-level; dose and dose schedule were digitised per cohort but dose itself was NOT a covariate in any final model -- the treatment effect enters through therapy class (statin vs anti-PCSK9) and through the achieved comparator-adjusted LDL-C lowering.",
    regions             = "Multinational; 54 randomised controlled trials published in English and indexed in PubMed or ClinicalTrials.gov up to 28 June 2022.",
    baseline_lipids     = "Trial-level baseline means (Volkova 2023 Supplementary Table 3): LDL-C 136.57 mg/dL (SD 27.71; range 87.9-195.56), HDL-C 45.82 mg/dL (SD 6.05; range 35.96-59.55), remnant cholesterol 30.72 mg/dL (SD 7.06; range 13.05-60), triglycerides 151 mg/dL (SD 25.32; range 113.27-264).",
    treatment_response  = "Trial-level comparator-adjusted changes (Volkova 2023 Supplementary Table 3): delta-LDL-C -36.91 mg/dL (SD 16.57; range -75 to +1.4), delta-remnant-C -4.37 mg/dL, delta-triglycerides -16.27 mg/dL, delta-HDL-C +1.53 mg/dL.",
    follow_up           = "1.5-10 years (mean 3.87, SD 1.83). Trials with mean cardiovascular follow-up shorter than 1.5 years, or fewer than 100 participants per treatment group, were excluded.",
    therapy_split       = "47 statin trials and 7 anti-PCSK9 trials (evolocumab, alirocumab, bococizumab, inclisiran).",
    other_demographics  = "Trial-level means (Volkova 2023 Supplementary Table 3): males 69.73%, hypertension 54.71%, diabetes 22.72%, smokers 21.49%, systolic blood pressure 135.27 mmHg. Missing covariate values (at most 13% of any measurement) were multiply imputed by predictive mean matching (100 datasets, mice 3.14.0) and all inferences pooled across imputations.",
    notes               = "Study-level (aggregate) meta-analytic data only -- no individual patient data. Each observation is a trial-level risk ratio computed from the per-arm event and subject counts, and the model predicts trial-level mean treatment effects. The composite MACE endpoint itself was DROPPED from the modelling after the Egger test showed significant publication bias (p < 0.001) for it alone; all downstream work is on the stand-alone components. Per-trial characteristics are tabulated in Volkova 2023 Supplementary Table 2."
  )

  ini({
    # =========================================================================
    # Volkova 2023 Table 1 -- final random-effects meta-regression per endpoint.
    # Common structure: log(RR) = beta_0 + beta_1 * PCSK9 + beta_2 * covariate
    # (+ beta_3 * interaction for nfMI). beta_0 is the statin-therapy intercept
    # at a zero covariate value; the covariates are entered UNCENTRED, so
    # several intercepts sit far from any observed log(RR) (HF at +0.554 is the
    # clearest example -- remnant cholesterol never approaches 0 mg/dL).
    #
    # These are ESTIMATED point estimates (metafor 3.0.2, pooled over 100
    # multiple imputations), not fixed constants, so they are NOT wrapped in
    # fixed(). The standard error printed beside each value is load-bearing:
    # Volkova 2023 Section 2.4 samples each effect size from N(mu, SE^2) when
    # propagating uncertainty into the sample-size calculation, and the
    # vignette reproduces that step.
    #
    # The p-values printed in Table 1 are NOT consistent with the ratio of the
    # estimate to its standard error under a normal reference (iST beta_0 is
    # the widest gap: -0.791/0.17 = -4.65 gives p ~ 3e-6, against a printed
    # p = 0.01). The estimates and standard errors are nevertheless confirmed
    # exactly by six independent back-calculations of the paper's own Results
    # narrative, so the p-value column alone is treated as unreliable; see
    # vignette Errata.
    # =========================================================================

    # -------- Coronary revascularization (CR) --------------------------------
    # log(RR) = beta_0 + beta_1 * PCSK9 + beta_2 * dLDLc
    int_cr        <-  -0.035 ; label("CR: log(RR) intercept, statin therapy at delta-LDLc = 0 (unitless)")                       # Volkova 2023 Table 1 CR beta_0 (SE 0.057; p = 0.55)
    e_pcsk9_cr    <-   0.22  ; label("CR: additive shift in log(RR) for anti-PCSK9 vs statin therapy (unitless)")                # Volkova 2023 Table 1 CR beta_1 (SE 0.067; p < 0.001). exp(0.22) = 1.246 reproduces the Results statement that treatment benefit is decreased by 24.6% [95% CI 9.3-42.1] in anti-PCSK9 populations
    e_dldlc_cr    <-   0.007 ; label("CR: change in log(RR) per 1 mg/dL of comparator-adjusted LDL-C change (1/(mg/dL))")        # Volkova 2023 Table 1 CR beta_2 (SE 0.002; p < 0.001). exp(0.007 * -38.67) = 0.763 reproduces the Results statement of -24% RR change per 1 mmol/L LDL-C decrease

    # -------- Cardiovascular mortality (CVM) ---------------------------------
    # log(RR) = beta_0 + beta_1 * PCSK9 + beta_2 * HDLc
    int_cvm       <-  -0.625 ; label("CVM: log(RR) intercept, statin therapy at baseline HDL-C = 0 mg/dL (unitless)")            # Volkova 2023 Table 1 CVM beta_0 (SE 0.253; p = 0.02)
    e_pcsk9_cvm   <-   0.082 ; label("CVM: additive shift in log(RR) for anti-PCSK9 vs statin therapy (unitless)")               # Volkova 2023 Table 1 CVM beta_1 (SE 0.067; p = 0.23; not significant)
    e_hdlc_cvm    <-   0.012 ; label("CVM: change in log(RR) per 1 mg/dL of baseline HDL-C (1/(mg/dL))")                         # Volkova 2023 Table 1 CVM beta_2 (SE 0.006; p = 0.05). Crosses log(RR) = 0 at HDL-C = 52.1 mg/dL for statins

    # -------- Heart failure (HF) ---------------------------------------------
    # log(RR) = beta_0 + beta_1 * PCSK9 + beta_2 * remC
    int_hf        <-   0.554 ; label("HF: log(RR) intercept, statin therapy at baseline remnant cholesterol = 0 mg/dL (unitless)")  # Volkova 2023 Table 1 HF beta_0 (SE 0.253; p = 0.05)
    e_pcsk9_hf    <-   0.097 ; label("HF: additive shift in log(RR) for anti-PCSK9 vs statin therapy (unitless)")                   # Volkova 2023 Table 1 HF beta_1 (SE 0.111; p = 0.4; not significant)
    e_remc_hf     <-  -0.023 ; label("HF: change in log(RR) per 1 mg/dL of baseline remnant cholesterol (1/(mg/dL))")               # Volkova 2023 Table 1 HF beta_2 (SE 0.009; p = 0.03). See vignette Errata: the Results gloss of "an increase of 13 mg/dl in baseline remC resulted in an approximately 2-fold decrease in the RR of HF" is a loose paraphrase, not an exact restatement of this coefficient

    # -------- Ischemic stroke (iST) ------------------------------------------
    # log(RR) = beta_0 + beta_1 * PCSK9 + beta_2 * Hypertension(%)
    int_ist       <-  -0.791 ; label("iST: log(RR) intercept, statin therapy at 0% hypertension prevalence (unitless)")          # Volkova 2023 Table 1 iST beta_0 (SE 0.17; p = 0.01)
    e_pcsk9_ist   <-  -0.309 ; label("iST: additive shift in log(RR) for anti-PCSK9 vs statin therapy (unitless)")               # Volkova 2023 Table 1 iST beta_1 (SE 0.122; p = 0.06; not significant)
    e_hypert_ist  <-   0.011 ; label("iST: change in log(RR) per 1 percentage point of hypertension prevalence (1/%)")           # Volkova 2023 Table 1 iST beta_2 (SE 0.003; p = 0.02). exp(0.011 * 50) = 1.73 with 95% CI 1.29-2.33 reproduces the Results statement that a 75%-hypertensive population is "1.7 [95% CI 1.3-2.3] fold less effective" than a 25%-hypertensive one

    # -------- Myocardial infarction, all (MI) --------------------------------
    # log(RR) = beta_0 + beta_1 * PCSK9 + beta_2 * HDLc
    int_mi        <-   0.407 ; label("MI: log(RR) intercept, statin therapy at baseline HDL-C = 0 mg/dL (unitless)")             # Volkova 2023 Table 1 MI beta_0 (SE 0.323; p = 0.22)
    e_pcsk9_mi    <-   0.033 ; label("MI: additive shift in log(RR) for anti-PCSK9 vs statin therapy (unitless)")                # Volkova 2023 Table 1 MI beta_1 (SE 0.104; p = 0.75; not significant). NOTE the Results sentence "For CR and MI, mean treatment benefit was decreased ... by 24.6% ... and 27.4% ..." refers to nfMI, not MI: exp(0.242) = 1.274 comes from the nfMI beta_1, whereas exp(0.033) = 1.034. See vignette Errata
    e_hdlc_mi     <-  -0.017 ; label("MI: change in log(RR) per 1 mg/dL of baseline HDL-C (1/(mg/dL))")                          # Volkova 2023 Table 1 MI beta_2 (SE 0.008; p = 0.04). Sign is opposite to CVM/TM: higher baseline HDL-C predicts GREATER MI benefit

    # -------- Stroke, unspecified (ST) ---------------------------------------
    # log(RR) = beta_0 + beta_1 * PCSK9 + beta_2(Prevention) * Age
    # Age has no main effect and Prevention has no main effect; the retained
    # term is the interaction alone, i.e. a prevention-category-specific slope
    # on uncentred mean age. Confirmed against the Results: at age 60 the model
    # gives RR 0.70 (statin/primary), 0.84 (statin/secondary), 0.66
    # (PCSK9/primary) and 0.80 (PCSK9/secondary) against published
    # 0.71 / 0.82 / 0.67 / 0.78.
    int_st              <-  -0.773 ; label("ST: log(RR) intercept, statin therapy at age 0 years (unitless)")                              # Volkova 2023 Table 1 ST beta_0 (SE 0.403; p = 0.06)
    e_pcsk9_st          <-  -0.056 ; label("ST: additive shift in log(RR) for anti-PCSK9 vs statin therapy (unitless)")                    # Volkova 2023 Table 1 ST beta_1 (SE 0.096; p = 0.56; not significant)
    e_age_st_primary    <-   0.007 ; label("ST: change in log(RR) per year of mean age, primary-prevention trials (1/year)")               # Volkova 2023 Table 1 ST beta_2 'Primary' (SE 0.007; p = 0.29)
    e_age_st_secondary  <-   0.010 ; label("ST: change in log(RR) per year of mean age, secondary-prevention trials (1/year)")             # Volkova 2023 Table 1 ST beta_2 'Secondary' (SE 0.006; p = 0.15)
    e_age_st_mixed      <-   0.011 ; label("ST: change in log(RR) per year of mean age, mixed/unstated-prevention trials (1/year)")        # Volkova 2023 Table 1 ST beta_2 'Both' (SE 0.006; p = 0.09)

    # -------- Total (all-cause) mortality (TM) -------------------------------
    # log(RR) = beta_0 + beta_1 * PCSK9 + beta_2 * HDLc
    int_tm        <-  -0.448 ; label("TM: log(RR) intercept, statin therapy at baseline HDL-C = 0 mg/dL (unitless)")             # Volkova 2023 Table 1 TM beta_0 (SE 0.17; p = 0.01)
    e_pcsk9_tm    <-   0.029 ; label("TM: additive shift in log(RR) for anti-PCSK9 vs statin therapy (unitless)")                # Volkova 2023 Table 1 TM beta_1 (SE 0.063; p = 0.65; not significant)
    e_hdlc_tm     <-   0.008 ; label("TM: change in log(RR) per 1 mg/dL of baseline HDL-C (1/(mg/dL))")                          # Volkova 2023 Table 1 TM beta_2 (SE 0.004; p = 0.03). Crosses log(RR) = 0 at HDL-C = 56 mg/dL for statins

    # -------- Coronary mortality (CM) ----------------------------------------
    # log(RR) = beta_0 + beta_1 * PCSK9 + beta_2 * dLDLc
    int_cm        <-  -0.019 ; label("CM: log(RR) intercept, statin therapy at delta-LDLc = 0 (unitless)")                       # Volkova 2023 Table 1 CM beta_0 (SE 0.064; p = 0.77)
    e_pcsk9_cm    <-   0.147 ; label("CM: additive shift in log(RR) for anti-PCSK9 vs statin therapy (unitless)")                # Volkova 2023 Table 1 CM beta_1 (SE 0.109; p = 0.19; not significant)
    e_dldlc_cm    <-   0.004 ; label("CM: change in log(RR) per 1 mg/dL of comparator-adjusted LDL-C change (1/(mg/dL))")        # Volkova 2023 Table 1 CM beta_2 (SE 0.002; p = 0.02). exp(0.004 * -38.67) = 0.857 reproduces the Results statement of -14% RR change per 1 mmol/L LDL-C decrease

    # -------- Fatal stroke (fST) ---------------------------------------------
    # log(RR) = beta_0 + beta_1 * PCSK9 + beta_2 * RD
    int_fst       <-  -0.113 ; label("fST: log(RR) intercept, statin therapy in trials without severe renal disease (unitless)")  # Volkova 2023 Table 1 fST beta_0 (SE 0.089; p = 0.23)
    e_pcsk9_fst   <-   0.05  ; label("fST: additive shift in log(RR) for anti-PCSK9 vs statin therapy (unitless)")                # Volkova 2023 Table 1 fST beta_1 (SE 0.289; p = 0.87; not significant)
    e_rd_fst      <-   0.654 ; label("fST: additive shift in log(RR) when the trial enrolled severe renal disease (unitless)")    # Volkova 2023 Table 1 fST beta_2 (SE 0.283; p = 0.04). exp(0.654) = 1.92 with 95% CI 1.10-3.35 reproduces the Results statement that fST benefit is "notably less explicit [by 1.9 (95% CI 1.1-3.3) times]" with late-stage renal disease

    # -------- Nonfatal myocardial infarction (nfMI) --------------------------
    # log(RR) = beta_0 + beta_1 * PCSK9 + beta_2 * dLDLc + beta_3 * RD * dLDLc
    # RD enters ONLY through the interaction -- there is no beta for RD alone.
    int_nfmi          <-  -0.057 ; label("nfMI: log(RR) intercept, statin therapy at delta-LDLc = 0 (unitless)")                                     # Volkova 2023 Table 1 nfMI beta_0 (SE 0.091; p = 0.54)
    e_pcsk9_nfmi      <-   0.242 ; label("nfMI: additive shift in log(RR) for anti-PCSK9 vs statin therapy (unitless)")                              # Volkova 2023 Table 1 nfMI beta_1 (SE 0.1; p = 0.02). exp(0.242) = 1.274 with 95% CI 1.05-1.55 reproduces the Results value of 27.4% [95% CI 4.7-55] (printed against the label "MI"; it is nfMI)
    e_dldlc_nfmi      <-   0.008 ; label("nfMI: change in log(RR) per 1 mg/dL of comparator-adjusted LDL-C change (1/(mg/dL))")                      # Volkova 2023 Table 1 nfMI beta_2 (SE 0.003; p < 0.001). exp(0.008 * -38.67) = 0.734 reproduces the Results statement of -27% RR change per 1 mmol/L LDL-C decrease
    e_rd_dldlc_nfmi   <-  -0.005 ; label("nfMI: additional change in log(RR) per 1 mg/dL of LDL-C change when the trial enrolled severe renal disease (1/(mg/dL))")  # Volkova 2023 Table 1 nfMI beta_3 (SE 0.002; p = 0.02). The renal-disease slope is 0.008 - 0.005 = 0.003, giving exp(0.003 * -38.67) = 0.890, which reproduces the Results statement that renal disease diminishes the decrease in RR to -11% [95% CI -39%-30%]

    # -------- Sample-size calculator constant --------------------------------
    # The two-sided 95% normal quantile that defines "statistically significant"
    # in Volkova 2023 Section 2.4: the upper bound of the 95% CI for the
    # meta-analysis average must be below zero, i.e. 1.96 * SE(theta) + theta = 0
    # at the minimal N. A structural constant of the method, not a fitted value.
    z975 <- fixed(1.96) ; label("Two-sided 95% standard-normal quantile used in the sample-size formula (unitless)")  # Volkova 2023 Section 2.4 and Supplementary File 2
  })

  model({
    # =======================================================================
    # Endpoint-specific meta-regression predictions (Volkova 2023 Table 1).
    # Each logrr_* is a TRIAL-LEVEL mean treatment effect on the log risk-ratio
    # scale; rr_* is the same quantity as a risk ratio. rr_* < 1 favours the
    # active treatment over the comparator.
    # =======================================================================
    logrr_cr   <- int_cr  + e_pcsk9_cr  * TRT_PCSK9I + e_dldlc_cr  * LDLC_DELTA
    logrr_cvm  <- int_cvm + e_pcsk9_cvm * TRT_PCSK9I + e_hdlc_cvm  * HDLC
    logrr_hf   <- int_hf  + e_pcsk9_hf  * TRT_PCSK9I + e_remc_hf   * REMC
    logrr_ist  <- int_ist + e_pcsk9_ist * TRT_PCSK9I + e_hypert_ist * DIS_HYPERT_PERCENT
    logrr_mi   <- int_mi  + e_pcsk9_mi  * TRT_PCSK9I + e_hdlc_mi   * HDLC
    logrr_tm   <- int_tm  + e_pcsk9_tm  * TRT_PCSK9I + e_hdlc_tm   * HDLC
    logrr_cm   <- int_cm  + e_pcsk9_cm  * TRT_PCSK9I + e_dldlc_cm  * LDLC_DELTA
    logrr_fst  <- int_fst + e_pcsk9_fst * TRT_PCSK9I + e_rd_fst    * RENALIMP_SEV

    # ST carries a prevention-category-specific slope on uncentred mean age and
    # no main effect for either age or prevention category.
    slope_age_st <-
      e_age_st_primary   * PREVENT_PRIMARY +
      e_age_st_secondary * PREVENT_SECONDARY +
      e_age_st_mixed     * PREVENT_MIXED
    logrr_st <- int_st + e_pcsk9_st * TRT_PCSK9I + slope_age_st * AGE

    # nfMI carries a renal-disease interaction on the LDL-C slope, with no
    # renal-disease main effect.
    slope_dldlc_nfmi <- e_dldlc_nfmi + e_rd_dldlc_nfmi * RENALIMP_SEV
    logrr_nfmi <- int_nfmi + e_pcsk9_nfmi * TRT_PCSK9I + slope_dldlc_nfmi * LDLC_DELTA

    rr_cr   <- exp(logrr_cr)
    rr_cvm  <- exp(logrr_cvm)
    rr_hf   <- exp(logrr_hf)
    rr_ist  <- exp(logrr_ist)
    rr_mi   <- exp(logrr_mi)
    rr_st   <- exp(logrr_st)
    rr_tm   <- exp(logrr_tm)
    rr_cm   <- exp(logrr_cm)
    rr_fst  <- exp(logrr_fst)
    rr_nfmi <- exp(logrr_nfmi)

    # =======================================================================
    # Minimal trial size required for statistical significance
    # (Volkova 2023 Section 2.4, derived in Supplementary File 2).
    #
    #            z^2 * (2 - 4 * k * exp(theta) + 2 * exp(theta))
    #   N(t) =  ------------------------------------------------
    #                    k * exp(theta) * theta^2
    #
    # N is the TOTAL enrolment across both arms, assuming equal allocation.
    # theta = LOGRR_TARGET is the effect size being powered for and k = K_CTRL
    # is the control-arm event proportion over the follow-up period. For a
    # composite MACE endpoint, Supplementary File 3 shows the identical formula
    # applies once theta and k are replaced by
    #   theta_sum = log(sum(k_i * exp(theta_i)) / sum(k_i))  and  k_sum = sum(k_i),
    # which is why this calculator takes both as external inputs rather than
    # hard-coding one endpoint composition.
    # =======================================================================
    rr_target  <- exp(LOGRR_TARGET)
    n_required <-
      z975 * z975 * (2 - 4 * K_CTRL * rr_target + 2 * rr_target) /
      (K_CTRL * rr_target * LOGRR_TARGET * LOGRR_TARGET)
  })
}
