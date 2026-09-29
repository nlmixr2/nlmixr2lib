Zierhut_2016_hcc_antiangiogenic_os_mbma <- function() {
  description <- paste0(
    "MBMA. Mixed-effects meta-regression of MEDIAN OVERALL SURVIVAL (mOS) in ",
    "patients with advanced hepatocellular carcinoma (aHCC) receiving systemic ",
    "antiangiogenic therapy (AAT), fitted to 68 treatment arms from 59 ",
    "published studies (4813 patients) identified by an OVID Medline / Embase ",
    "literature search through late 2012. The response is ln(mOS) in months, ",
    "so every covariate is additive on the log scale and therefore ",
    "PROPORTIONAL on the mOS scale. The model carries two mutually exclusive ",
    "intercepts - one for AAT arms (8.49 months) and one for placebo arms ",
    "(7.06 months), both quoted at the analysis-population median of 34% ",
    "hepatitis-B-positive patients and no prior systemic therapy - plus five ",
    "additive covariate effects: sorafenib rather than another AAT (+21%), ",
    "concomitant locoregional therapy (+42%), median-centred percentage of ",
    "hepatitis-B-positive patients (-0.4% of mOS per percentage point), prior ",
    "systemic therapy (-6.7%) and concomitant chemotherapy (-4.0%). The last ",
    "two were FORCED into the model on clinical grounds and are not ",
    "statistically significant; their 95% confidence intervals span both ",
    "benefit and harm. This is a STUDY-LEVEL model: one record is one ",
    "published trial arm and the prediction is that arm's median overall ",
    "survival, never an individual patient's survival time and never a drug ",
    "concentration. It has no dose, no exposure and no time axis - mOS is a ",
    "single scalar per arm - so it cannot be used to simulate a survival ",
    "curve. Variability is BETWEEN-STUDY (xi = 0.216 on the ln(mOS) scale) ",
    "plus a residual that the source SCALES BY EACH ARM'S OWN REPORTED ",
    "STANDARD ERROR, with the unit-weight SD held at 1; see the residual ",
    "label and the vignette for why a bare stochastic rxSolve() overstates ",
    "the residual roughly five-fold. The paper's own clinical trial ",
    "simulations for the axitinib phase II study NCT01210495 are reproduced ",
    "exactly by the typical-value predictions: 6.16, 7.40 and 8.95 months for ",
    "placebo, non-sorafenib AAT and sorafenib in a 50%-hepatitis-B, ",
    "all-prior-systemic-therapy population."
  )

  reference <- paste(
    "Zierhut ML, Chen Y, Pithavala YK, Nickens DJ, Valota O, Amantea MA.",
    "Clinical trial simulations from a model-based meta-analysis of studies in",
    "patients with advanced hepatocellular carcinoma receiving antiangiogenic",
    "therapy.",
    "CPT Pharmacometrics Syst Pharmacol. 2016;5(5):274-282.",
    "doi:10.1002/psp4.12078.",
    "Structural model: the displayed ln(mOS_ij) equation in Results 'Final",
    "model' and the displayed fixed/random-effect equation in Methods 'Model",
    "building'.",
    "Parameter values: Table 3.",
    "Covariate medians and analysis-set summaries: Tables 1 and 2.",
    sep = " "
  )

  vignette <- "Zierhut_2016_hcc_antiangiogenic_os_mbma"

  # The single random effect is a BETWEEN-STUDY effect (gamma_j) added
  # directly to the ln(mOS) response, not an IIV on a named fixed-effect
  # parameter, so it cannot follow the popPK eta<transformed-param-name>
  # convention. Declared paper-specific per the Goteti_2024_SLE_mbma and
  # Takechi_2025_nemolizumab_mbma_iga precedent.
  paper_specific_etas <- c("eta_study_lnmos")

  # Study-level MBMA with no drug input: nothing is dosed, nothing is
  # absorbed and the output is a survival duration rather than a
  # concentration. The dosing / concentration strings are declared "n/a"
  # in the Zhang_2025_bietti_crystalline_dystrophy_mbma style so
  # checkModelConventions() is not asked to reconcile a dimensional pair
  # that does not exist in this model.
  units <- list(
    time = "month",
    dosing = "n/a (no dosing; treatment enters as trial-arm indicators)",
    concentration = "n/a (no drug concentration is modelled)",
    response = "month (mos, the study arm's median overall survival; the fitted response lnmos is its natural logarithm and is dimensionless)"
  )

  covariateData <- list(
    ON_TREATMENT = list(
      description = "Indicator that the trial arm received a systemic antiangiogenic therapy (AAT) rather than placebo. 1 = AAT arm, 0 = placebo arm.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo arm, which selects the separately estimated lmos_pbo intercept rather than a shift away from lmos_aat)",
      notes = paste0(
        "MBMA trial-arm indicator. Zierhut 2016 carries the paper's AAT and ",
        "PBO indicators as two SEPARATE intercepts with no shared intercept, ",
        "and states in Results 'Final model' that the two are 'mutually ",
        "exclusive and collectively exhaustive indicator variables'. This ",
        "model therefore consumes a SINGLE column: AAT = ON_TREATMENT and the ",
        "paper's PBO = 1 - ON_TREATMENT, computed inside model(). Carrying ",
        "one column rather than two makes the exhaustiveness structural, so ",
        "an input row cannot silently set both indicators to 1 (or to 0) and ",
        "get a meaningless sum or a zero intercept. Every patient in the ",
        "analysis received best supportive care in addition to the listed ",
        "treatment, so BSC is implicit in BOTH intercepts and is not a ",
        "separate covariate. The AAT class pools 15 agents across 62 arms ",
        "(Table 1: 26 sorafenib arms and 36 'other AAT agent' arms, the ",
        "latter spanning 14 agents at 1-9 arms each); sorafenib arms carry ",
        "SORAFENIB = 1 on top of this indicator. Setting ON_TREATMENT = 0 ",
        "together with SORAFENIB = 1 is outside the source's calibration - no ",
        "placebo arm was a sorafenib arm."
      ),
      source_name = "AAT / PBO (Zierhut 2016 Results 'Final model' equation and Table 1)"
    ),
    SORAFENIB = list(
      description = "Indicator that the antiangiogenic therapy administered in this trial arm was sorafenib rather than one of the other antiangiogenic agents.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (a non-sorafenib AAT arm, or a placebo arm)",
      notes = paste0(
        "MBMA per-drug trial-arm indicator, applied ON TOP OF ON_TREATMENT ",
        "rather than instead of it: Zierhut 2016 Table 3 reports ",
        "exp(theta_SOR) = 121%, i.e. a sorafenib arm's median overall ",
        "survival is 1.21 times that of a non-sorafenib AAT arm with the same ",
        "covariates, so a sorafenib arm needs BOTH ON_TREATMENT = 1 and ",
        "SORAFENIB = 1. Table 1 records 26 sorafenib arms of the 68 total. ",
        "NAMING: the register answers 'what do I call a per-drug MBMA ",
        "trial-arm indicator?' two incompatible ways - the bare-INN family ",
        "(NAPROXEN, TRAMADOL, TAPENTADOL, ACETAMINOPHEN, DICLOFENAC) with the ",
        "rule statement in the TRT_BENRALIZUMAB notes ('a bare INN only for ",
        "MBMA trial-arm indicators'), versus the TRT_<INN> form the Sandra ",
        "2024 anti-HBV siRNA batch adopted for its 13 MBMA arm indicators. ",
        "This model follows the register's explicit rule statement for this ",
        "exact case (bare INN); the conflict between the two forms is ",
        "already under review by the maintainers. If that review resolves ",
        "in favour of TRT_<INN>, rename ",
        "this column to TRT_SORAFENIB along with the rest of the bare-INN ",
        "family. Distinct from CP_SORAFENIB_NGML, which is an individual ",
        "patient's measured sorafenib plasma concentration driving a ",
        "concentration-response model; this column carries no exposure ",
        "information at all and no sorafenib dose was a covariate in this ",
        "meta-analysis."
      ),
      source_name = "SOR (Zierhut 2016 Results 'Final model' equation, Table 1 and Table 3)"
    ),
    CONMED_LOCOREGIONAL = list(
      description = "Indicator that the trial arm received a concomitant locoregional liver-directed therapy alongside its systemic treatment.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (systemic therapy alone; no concomitant locoregional procedure)",
      notes = paste0(
        "MBMA trial-arm study-design indicator. Zierhut 2016 Results 'Final ",
        "model' defines the locoregional class explicitly: 'LOC therapy ",
        "consisted of transarterial chemoembolization, SIRspheres, or ",
        "cryoablation', so the column pools intra-arterial, radioembolic and ",
        "ablative liver-directed procedures under one effect. Table 1 records ",
        "7 arms with concomitant locoregional therapy of the 68 total. The ",
        "effect is large and significant - exp(theta_LOC) = 142%, the biggest ",
        "single covariate in the model - but it rests on those 7 arms, so ",
        "treat it as the source does, as an arm-level study-design contrast ",
        "rather than a transportable per-patient treatment effect. CONCOMITANT ",
        "only: a patient who had locoregional therapy BEFORE enrolment is not ",
        "captured here (the paper carries no prior-locoregional covariate at ",
        "all, which matters in aHCC because prior transarterial ",
        "chemoembolization is common). Modality-class member of the CONMED_ ",
        "family alongside CONMED_CHEMO and CONMED_PLATIN, which likewise name ",
        "a class of concomitant therapy rather than a single INN."
      ),
      source_name = "LOC / 'Concomitant LOC' (Zierhut 2016 Results 'Final model' equation, Table 1 and Table 3)"
    ),
    DIS_CHB_PERCENT = list(
      description = "Percentage (0-100) of the trial arm's enrolled population who were hepatitis B virus positive at baseline.",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "MBMA trial-arm cohort-prevalence covariate, scaled in PERCENT and ",
        "not as a fraction: Zierhut 2016 Table 3 gives the coefficient in ",
        "'d.u./%' and the Results gloss it as 'mOS decreases by ~0.4% for ",
        "every 1% increase in HBV-positive patients', which reproduces as ",
        "exp(-0.00418) - 1 = -0.417%. MEDIAN-CENTRED at 34.0%, per Methods ",
        "'Model building' ('all continuous covariates were centered to the ",
        "median value') and Table 2, whose HBV row reads 38.2 (SD 26.1) mean ",
        "and 34.0 [4.8, 100] median, reported in 70% of arms. The centring ",
        "value is confirmed independently by the Results sentence defining ",
        "the typical population as '34% of patients with HBV' and by the ",
        "clinical trial simulation, whose published placebo median of 6.16 ",
        "months reproduces to three digits only when 34.0 is subtracted. ",
        "Arm-level HBV percentages were IMPUTED where unreported, as the ",
        "patient-number-weighted mean of the arms that did report a value ",
        "(Methods 'Data processing'), so 30% of the arms carry an imputed ",
        "value. Supply 0-100, not 0-1."
      ),
      source_name = "HBV / 'Percent with HBV' (Zierhut 2016 Results 'Final model' equation, Table 2 and Table 3)"
    ),
    PRIOR_SYSTEMIC = list(
      description = "Indicator that the ENTIRE enrolled population of the trial arm had received prior systemic anticancer therapy, i.e. the arm is wholly second-line or later.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the arm was not wholly pretreated: first-line arms and arms of mixed treatment-line composition both take 0)",
      notes = paste0(
        "MBMA trial-arm cohort-composition flag, NOT an individual patient's ",
        "treatment history. The threshold is all-or-nothing: Zierhut 2016 ",
        "Table 1 counts 7 arms of 'All second line patients' and 43 of 'All ",
        "first line patients' out of 68, and the Discussion names the same 7 ",
        "as the arms in which 'the entire population had received prior ",
        "therapy (n = 7)'. The remaining 18 arms are of mixed line and take 0 ",
        "here, so this column is NOT 1 - LINE_1L (which would make 25 arms ",
        "positive) and the underlying continuous 'percent with prior systemic ",
        "tx' of Table 2 (mean 16.6, median 0 [0, 100]) is NOT what the model ",
        "consumes - Methods 'Model building' classes this among the binary ",
        "covariates, whose coefficients are a shift rather than a per-percent ",
        "slope, and Table 3 accordingly prints exp(theta_PTx) as a percentage ",
        "multiplier (93.3%) rather than in the per-percent 'd.u./%' units it ",
        "uses for the HBV slope. SYSTEMIC only: prior locoregional therapy, ",
        "surgery and radiotherapy do not set this flag, which is why the ",
        "broader PRIOR_ANTICANCER (any modality, including radiotherapy and ",
        "debulking surgery) is the wrong column in aHCC, where prior ",
        "transarterial chemoembolization is near-universal. FORCED into the ",
        "model on clinical grounds and not statistically significant: the 95% ",
        "confidence interval on exp(theta_PTx), 74.0-118%, spans both benefit ",
        "and harm, and the Discussion cautions that with only 7 informative ",
        "arms 'their true impact may not have been fully accounted for'."
      ),
      source_name = "PTx / 'All second line patients' (Zierhut 2016 Results 'Final model' equation, Table 1 and Table 3)"
    ),
    CONMED_CHEMO = list(
      description = "Indicator that the trial arm received concomitant cytotoxic chemotherapy alongside its antiangiogenic or placebo treatment.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant chemotherapy)",
      notes = paste0(
        "MBMA trial-arm study-design indicator; Zierhut 2016 Table 1 records ",
        "13 arms with concomitant chemotherapy of the 68 total. Reuses the ",
        "registered CONMED_CHEMO column in its general sense - any ",
        "chemotherapy backbone collapsed into one binary study-design flag - ",
        "outside the anti-PD-(L)1 setting that founded the entry; no ",
        "chemotherapy agent, dose or schedule is resolved by this ",
        "meta-analysis. The functional form here is additive on the log ",
        "response, i.e. multiplicative on mOS: exp(theta_CTx) = 96.0%. FORCED ",
        "into the model on clinical grounds and not statistically significant ",
        "(95% CI 80.7-114%). The Discussion warns explicitly against ",
        "over-reading the apparent 4% detriment, and notes that whatever ",
        "effect exists is 'in addition to the benefit of AAT' - i.e. it may ",
        "reflect chemotherapy removing some of the antiangiogenic benefit ",
        "rather than harming survival outright. Distinct from CONMED_PLATIN ",
        "and the per-agent CONMED_<INN> columns, none of which this source ",
        "resolves."
      ),
      source_name = "CTx / 'Concomitant CTx' (Zierhut 2016 Results 'Final model' equation, Table 1 and Table 3)"
    )
  )

  # Covariates Zierhut 2016 prospectively tested in the forward-selection
  # step but did NOT retain in the final model. Recorded for provenance
  # only; none is referenced in model(), and the paper reports no usable
  # point estimate for any of them.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Trial-arm mean age.",
      units = "years",
      type = "continuous",
      notes = "Zierhut 2016 Table 2: 60.2 years (SD 6.8), median 60 [47, 75], reported in 83% of arms. Screened as 'median age' in the forward-selection covariate set (Methods 'Model building') and not retained."
    ),
    SEXF = list(
      description = "Percentage of the trial arm that was male; the paper screened the male rather than the female percentage.",
      units = "%",
      type = "continuous",
      notes = "Zierhut 2016 Table 2 reports 'Percent male' 78.6 (SD 9.7), median 80.1 [55, 94.7], reported in 79% of arms; the arm-level female percentage is its complement. Screened as 'percent of population ... male' and not retained."
    ),
    RACE_ASIAN_PCT = list(
      description = "Asian trial-site indicator; the source screened a whole-trial Asian-site flag rather than a continuous Asian percentage.",
      units = "(binary)",
      type = "binary",
      notes = "Zierhut 2016 Table 1 records 22 trials of 59 conducted at Asian study sites. Not significant in the forward selection. Re-examined post hoc as a placebo-arm-specific covariate after the axitinib trial reported an Asian-site benefit, where it showed a trend only (delta MOF = -2.63, short of the 3.84 entry threshold); the Discussion reports the resulting model-predicted mOS ratios of 0.821 (Asian) and 0.902 (non-Asian) but no fitted coefficient, so nothing is encodable."
    ),
    PS_ECOG_0_PCT = list(
      description = "Percentage of the trial arm with an ECOG performance status of 0.",
      units = "%",
      type = "continuous",
      notes = "Zierhut 2016 Table 2: 43.1% (SD 17.4), median 36.3 [0, 100], reported in 58% of arms. A companion 'Percent with ECOG < 2' row reads 90.8% (SD 9.5), median 94.4 [76.8, 100], reported in 73% of arms. Neither was retained."
    ),
    DIS_CHILDPUGH_A_PCT = list(
      description = "Percentage of the trial arm with Child-Pugh class A liver function.",
      units = "%",
      type = "continuous",
      notes = "Zierhut 2016 Table 2: 82.7% (SD 16.1), median 83.9 [7.4, 100], reported in 80% of arms; the class B row reads 13.7% (SD 15.5), median 8.5 [0, 92.6]. Child-Pugh B percentage is named in Methods 'Model building' as one of the screened covariates and was not retained. No canonical column is minted for either, because neither is used by any model."
    ),
    YEAR_PUB = list(
      description = "Publication year of the source study.",
      units = "year",
      type = "continuous",
      notes = "Zierhut 2016 Table 2: 2010.6 (SD 1.70), median 2011 [2005, 2012]. Tested prospectively as a proxy for improving best supportive care over time and not significant (delta MOF = -0.449). The Discussion argues this is the most likely explanation for the phase II trial's higher observed mOS in BOTH arms, while noting the dataset shows no evidence of the trend."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 4813L,
    n_studies = 59L,
    n_arms = 68L,
    age_range = "Trial-arm mean ages 47-75 years (mean of arm means 60.2, SD 6.8; median 60), reported in 83% of arms (Zierhut 2016 Table 2).",
    weight_range = "Body weight was not collected and was not among the screened covariates.",
    sex_female_pct = 21.4,
    disease_state = paste0(
      "Advanced (unresectable or metastatic) hepatocellular carcinoma; 98.1% ",
      "of the pooled population had advanced disease (median 100% per arm, ",
      "range 69-100%). Liver function was predominantly preserved: 82.7% ",
      "Child-Pugh A and 13.7% Child-Pugh B on average across arms. Hepatitis ",
      "B positivity averaged 38.2% (median 34.0%, range 4.8-100%), prior ",
      "vascular invasion 35.2%, ECOG performance status 0 in 43.1% and ECOG ",
      "below 2 in 90.8%. Prior chemotherapy 20.5% and prior systemic therapy ",
      "16.6% (median 0%), i.e. most arms were treatment-naive."
    ),
    dose_range = "Not modelled. Dose was not a covariate: treatment enters only as the AAT / placebo / sorafenib / concomitant-therapy arm indicators, and no exposure metric was available at the aggregate level.",
    regions = "Multinational. 22 of the 59 trials were conducted at Asian study sites; the remainder were non-Asian or mixed.",
    treatments = paste0(
      "68 arms: 6 placebo, 26 sorafenib and 36 other-AAT arms spanning 14 ",
      "further agents at 1-9 arms each. The agents tested in covariate ",
      "selection were brivanib, vandetanib, bevacizumab, erlotinib, ",
      "lapatinib, sunitinib and sorafenib, plus eight grouped as 'other' ",
      "(cabozantinib, cetuximab, cediranib, PTK787/ZK222584, imatinib, ",
      "linifanib, tivantinib and TSU-68). 24 arms used combination therapy ",
      "(one of them a triple), 13 had concomitant chemotherapy and 7 ",
      "concomitant locoregional therapy. All patients received best ",
      "supportive care in addition to the listed treatment."
    ),
    trial_design = "10 blinded trials, 9 randomised, 43 phase II or III; 49 of the 59 were single-arm studies. Arm sizes ranged from 10 to 544 patients (mean 70.8, median 42).",
    endpoint = paste0(
      "Median overall survival. Across the 68 arms the median of the arm mOS ",
      "values was 9.4 months, range 4.2-20.8; on the modelled ln scale the ",
      "mean was 2.18 (SD 0.36) with median 2.24 [1.44, 3.04]. The per-arm ",
      "standard error of ln(mOS) averaged 0.186 (SD 0.11), median 0.183 ",
      "[0.045, 0.61]. NOTE that Table 2's printed arithmetic mean mOS of 6.6 ",
      "months (SD 3.3) is not reproducible and is treated here as a ",
      "typographical error: with 68 arms, a minimum of 4.2 and a median of ",
      "9.4, the arithmetic mean cannot fall below 6.8, and the paper's own ",
      "ln(mOS) row implies about 9.4 months. Only the mean is affected; the ",
      "median and range are internally consistent with the ln(mOS) row via ",
      "exp(2.24) = 9.4, exp(1.44) = 4.2 and exp(3.04) = 20.9. See the ",
      "vignette Errata."
    ),
    notes = paste0(
      "Study-level (aggregate) literature data only - no individual patient ",
      "data. One record is one published treatment arm, and the model ",
      "predicts that arm's median overall survival. Sources were screened ",
      "from 350 publications via OVID Medline, Embase, Embase Alerts and ",
      "Medline in Process, covering data available through late 2012; ",
      "retrospective analyses, protocols, reviews, case studies, duplicate ",
      "trials and arms not reporting mOS were excluded, as were arms whose ",
      "population characteristics would have confounded simulation of the ",
      "axitinib phase II trial. Variables reported for fewer than 65% of the ",
      "pooled patients were dropped; those retained were imputed when missing ",
      "as the patient-number-weighted mean across the reporting arms, so ",
      "several covariate columns carry imputed values for a substantial ",
      "fraction of arms (HBV, the model's only continuous covariate, was ",
      "reported in 70% of arms). Standard errors of ln(mOS) were likewise ",
      "derived rather than reported: from the published 95% CI where ",
      "available ((p97.5 - p2.5)/3.92 on the ln scale), otherwise from a ",
      "reported range ((ln(max) - ln(min))/(X * sqrt(N)) with X = 3, 4, 5 or ",
      "6 by arm size), otherwise imputed as sqrt(median(SD^2)/N) over the ",
      "arms with a calculable SE. Estimated in NONMEM 7.2 with FOCE-I after ",
      "model building in R's nlme. No significant publication bias was found ",
      "by funnel plot or metabias() (p = 0.22 overall)."
    )
  )

  ini({
    # =======================================================================
    # Zierhut 2016 Table 3 -- final MBMA parameter estimates. The response is
    # ln(mOS) in months, so:
    #   * the two intercepts are reported as exp(theta) in MONTHS, and are
    #     entered here as log(<months>);
    #   * the four binary covariate effects are reported as exp(theta) in
    #     PERCENT of the uncovariated value, and are entered as
    #     log(<percent>/100);
    #   * the one continuous covariate effect is reported on the log domain
    #     directly and is entered unchanged.
    # All seven are ESTIMATED point estimates, not fixed constants, so none
    # is wrapped in fixed(). The printed 95% CIs are 'parameter estimate
    # +/- 1.96 SE' (Table 3 footnote a) on the scale each row is printed in.
    # =======================================================================

    # -------- Treatment-arm intercepts ------------------------------------
    # AAT and PBO are mutually exclusive and collectively exhaustive, so
    # these are two parallel intercepts rather than a reference plus an
    # offset -- hence the symmetric stratum suffixes and no bare `lmos`.
    # Both are quoted at the analysis-set median of 34% HBV positivity with
    # no prior systemic therapy, no concomitant locoregional therapy and no
    # concomitant chemotherapy.
    lmos_aat <- log(8.49)  ; label("Typical median overall survival on antiangiogenic therapy (mOS, months)")  # Zierhut 2016 Table 3, exp(theta_AAT) = 8.49 mo (95% CI 7.72-9.36)
    lmos_pbo <- log(7.06)  ; label("Typical median overall survival on placebo (mOS, months)")                  # Zierhut 2016 Table 3, exp(theta_PBO) = 7.06 mo (95% CI 6.28-7.86)

    # -------- Additive covariate effects on ln(mOS) ------------------------
    # Because the response is log-transformed, an additive shift here is a
    # PROPORTIONAL change in mOS (Methods 'Model building': 'the additive
    # covariates were interpreted as proportional to mOS').
    e_sorafenib_mos <- log(1.21)   ; label("Proportional change in mOS when the antiangiogenic agent is sorafenib (unitless)")            # Zierhut 2016 Table 3, exp(theta_SOR) = 121% (95% CI 114-128)
    e_loc_mos       <- log(1.42)   ; label("Proportional change in mOS with concomitant locoregional therapy (unitless)")                 # Zierhut 2016 Table 3, exp(theta_LOC) = 142% (95% CI 122-165)
    e_hbv_mos       <- -0.00418    ; label("Change in ln(mOS) per percentage point of hepatitis-B-positive patients above 34% (1/%)")     # Zierhut 2016 Table 3, theta_HBV = -0.00418 d.u./% (95% CI -7.02e-3 to -1.34e-3); exp(-0.00418) - 1 = -0.417% of mOS per percentage point
    e_ptx_mos       <- log(0.933)  ; label("Proportional change in mOS when the whole arm had prior systemic therapy (unitless)")         # Zierhut 2016 Table 3, exp(theta_PTx) = 93.3% (95% CI 74.0-118); forced into the model and not significant
    e_ctx_mos       <- log(0.960)  ; label("Proportional change in mOS with concomitant chemotherapy (unitless)")                         # Zierhut 2016 Table 3, exp(theta_CTx) = 96.0% (95% CI 80.7-114); forced into the model and not significant

    # -------- Between-study random effect ----------------------------------
    # gamma_j in the Results 'Final model' equation: one draw per TRIAL,
    # shared by every arm of that trial, additive on ln(mOS). Table 3
    # reports xi as the SD of gamma (footnote d), so the variance entered
    # here is 0.216^2. This is BETWEEN-STUDY variability, not between-subject
    # IIV -- two arms of the same trial share a draw, and no draw describes
    # a patient.
    eta_study_lnmos ~ 0.046656  # Zierhut 2016 Table 3, xi = 0.216 (95% CI 0.166-0.266), shrinkage 21.2%; variance = 0.216^2

    # -------- Residual ------------------------------------------------------
    # WARNING - this SD is the UNIT-WEIGHT residual and is roughly five times
    # the residual any real trial arm carries. Zierhut 2016 scales epsilon by
    # each arm's OWN reported standard error of ln(mOS) (the eps * SE term of
    # the Results 'Final model' equation), and then fixes the unit-weight SD
    # at 1 because a value below 1 would mean the model was less variable
    # than the published arms themselves (Methods 'Model building'). Across
    # the analysis set that per-arm SE had median 0.183 and range 0.045-0.61
    # (Table 2), so an unscaled draw from N(0, 1) is about 5-fold too wide at
    # the median arm. The scaling factor is a property of each published arm
    # rather than of the model, so it has no ini() home; supply it downstream,
    # or use rxode2::zeroRe(mod, 'sigma') to drop the residual entirely and
    # keep the between-study effect, as the validation vignette does.
    #
    # Declared as expSd / lnorm() rather than addSd / add(): the source's
    # residual is additive on ln(mOS), which is exactly a log-normal residual
    # on mOS itself, and observing mos rather than lnmos keeps the model's
    # output on the months scale the paper reports.
    expSd <- fixed(1)  ; label("Residual SD of ln(mOS) at unit study weight (dimensionless); the source scales it per record by that arm's own standard error of ln(mOS)")  # Zierhut 2016 Table 3, sigma = 1 (fixed), shrinkage 41.5%; fixing it changed no fixed effect by more than 10% and reduced xi by about 20% relative to estimating sigma = 0.568
  })

  model({
    # =======================================================================
    # Zierhut 2016, Results 'Final model':
    #
    #   ln(mOS_ij) = theta_AAT * AAT + theta_PBO * PBO + gamma_j
    #                + theta_SOR * SOR + theta_LOC * LOC + theta_HBV * HBV
    #                + theta_PTx * PTx + theta_CTx * CTx + eps_ij * SE_ij
    #
    # AAT and PBO are mutually exclusive and collectively exhaustive, so
    # PBO is carried as (1 - ON_TREATMENT) and only one input column is
    # needed. The eps * SE weighting is documented in the ini() residual
    # label and is applied downstream, not here.
    # =======================================================================

    # Median-centred continuous covariate. Methods 'Model building': 'All
    # continuous covariates were centered to the median value to simplify
    # interpretation of parameter estimates.' The median hepatitis-B
    # positivity of the analysis set is 34.0% (Table 2), which the Results
    # restate as the typical population's '34% of patients with HBV'.
    hbvCentred <- DIS_CHB_PERCENT - 34.0

    lnmos <-
      lmos_aat * ON_TREATMENT +
      lmos_pbo * (1 - ON_TREATMENT) +
      e_sorafenib_mos * SORAFENIB +
      e_loc_mos * CONMED_LOCOREGIONAL +
      e_hbv_mos * hbvCentred +
      e_ptx_mos * PRIOR_SYSTEMIC +
      e_ctx_mos * CONMED_CHEMO +
      eta_study_lnmos

    # The study arm's median overall survival on its natural scale, in
    # months. The paper models ln(mOS) because mOS was assumed log-normally
    # distributed (Methods 'Data processing'), so observing mos under
    # lnorm() is the same statement as observing lnmos under add().
    mos <- exp(lnmos)

    mos ~ lnorm(expSd)
  })
}
