Chan_2017_ipf_fvc_mbma <- function() {
  description <- paste(
    "MBMA. Longitudinal model-based meta-analysis of the placebo-corrected",
    "change from baseline in percent-predicted forced vital capacity",
    "(%predicted FVC) for 15 treatment regimens in idiopathic pulmonary",
    "fibrosis, fit by maximum likelihood (R 3.1.2, nlme::gnls) to arm-level",
    "summary data from 43 arms (148 arm-timepoints, 4,919 subjects) in 20",
    "published trials with treatment durations of 8-104 weeks (Chan 2017).",
    "Each treatment carries its own maximum effect Emax in %predicted FVC;",
    "a single shared empirical time course 1 - exp(-kel * time) governs how",
    "fast that maximum is approached (kel = exp(lambda) = 0.0417 /week, so",
    "50% of the maximum at 16.6 weeks and 90% at 55 weeks); and the arm's",
    "mean baseline %predicted FVC scales the whole treatment effect through",
    "the power term (FVC_PCTPRED / 74)^1.48. Pirfenidone is the only",
    "treatment with a dose-response: a step function selected by DOSE_HIGH",
    "gives 2.42 at 1,197 mg/day and 3.87 at 2,403 mg/day. Pirfenidone and",
    "nintedanib were the only regimens whose 95% CI excluded zero.",
    "CRITICAL SCOPE LIMIT: the output fvcppcfb is the PLACEBO-CORRECTED",
    "change from baseline, i.e. the active-minus-placebo difference only.",
    "Chan 2017 modelled the placebo arm NONPARAMETRICALLY -- one free",
    "estimate per trial per timepoint, deliberately so that no distributional",
    "assumption was imposed on the highly variable IPF placebo response --",
    "and those per-trial placebo estimates are not tabulated in the paper.",
    "This model therefore CANNOT produce an absolute %predicted FVC",
    "trajectory; to reproduce a figure like Chan 2017 Figure 3 the user must",
    "supply the trial's own observed placebo arm. Suitable simulation scope",
    "is the study-arm-mean placebo-corrected treatment effect, NOT",
    "individual-patient FVC.",
    sep = " "
  )

  reference <- paste(
    "Chan P, Bax L, Chen C, Zhang N, Huang SP, Soares H, Rosen G, AbuTarif M.",
    "Model-based Meta-Analysis on the Efficacy of Pharmacological Treatments",
    "for Idiopathic Pulmonary Fibrosis.",
    "CPT Pharmacometrics Syst Pharmacol. 2017;6(10):695-704.",
    "doi:10.1002/psp4.12227.",
    "Final parameter estimates are in Table 2; the treatment-effect ranking",
    "is duplicated in Figure 2; the model equation is Eq. 1 with the",
    "time-course and baseline-covariate forms given in the Table 2 note.",
    sep = " "
  )

  vignette <- "Chan_2017_ipf_fvc"

  # Chan 2017 Eq. 1 places the between-trial-arm random effect eta_ij
  # ADDITIVELY on the response, not on any structural parameter, so it has no
  # typical-value partner for the eta<x> <-> x pairing check to find. This is
  # the documented `paper_specific_etas` case: "IIV names whose typical-value
  # parameter is a paper-mechanistic structural equation rather than a 1-to-1
  # lX ini parameter".
  paper_specific_etas <- c("eta_arm_fvcppcfb")

  units <- list(
    time = "week (weeks since randomization; lambda / kel are reported per week)",
    dosing = "(no dose events; each treatment's regimen is fixed at the single dose level studied and is identified by its arm indicator, with the pirfenidone 1,197-vs-2,403 mg/day step selected by DOSE_HIGH)",
    concentration = "% predicted (placebo-corrected change from baseline in FVC percent predicted; the modelled quantity is a difference between two spirometric percentages, not a drug concentration, so the dosing-versus-concentration dimensional check is not applicable and the dosing string is parenthesised to skip it)"
  )

  covariateData <- list(
    # ------------------------------------------------------------------
    # Arm-mean baseline lung function. The ONLY covariate Chan 2017
    # retained out of the six prespecified ones (baseline %predicted FVC,
    # disease duration, age, and the proportions male, white, and
    # current/former smokers).
    # ------------------------------------------------------------------
    FVC_PCTPRED = list(
      description = "Arm-mean (or arm-median) baseline forced vital capacity expressed as a percent of the sex / age / height / ethnicity reference-equation predicted value, at randomization.",
      units = "% predicted (numeric percentage, e.g. 76.2 for the CAPACITY-2 arms; not a fraction 0.762)",
      type = "continuous",
      reference_category = "n/a -- enters as the power form (FVC_PCTPRED / 74)^e_fvc_pctpred_emax. 74 is the approximate median baseline %predicted FVC across the analysis dataset (Chan 2017 Results, 'Model development', and the Table 2 note), so an arm at 74% carries a covariate factor of exactly 1 and its treatment effect equals the tabulated Emax.",
      notes = paste(
        "MBMA study-arm-level covariate: the published arm's MEAN baseline,",
        "not an individual patient's. Chan 2017 imputed this value by",
        "multivariate linear regression for the trials that did not report it",
        "(0-36% of trials were missing a prespecified covariate; Results,",
        "'Exploratory analysis'), and derived %predicted FVC from observed",
        "absolute FVC via the Crapo 1981 reference equation for the three",
        "nintedanib trials (Methods, 'Database construction'). Observed",
        "arm-level values seen in the paper: 76.2 (CAPACITY-2 / PIPF-004) and",
        "74.2 (TOMORROW) from the Figure 3 panel annotations, and the",
        "stratified low / high baseline group medians used for the Figure 5",
        "simulation. NOTE a source inconsistency on those two medians: the",
        "Results text ('Model application') says 66% and 76%, while the",
        "Figure 5 legend says FVC=67 and FVC=76 -- see the vignette Errata.",
        "The covariate acts on the TREATMENT EFFECT, not on the placebo",
        "trajectory or on a baseline intercept: a higher baseline %predicted",
        "FVC predicts a LARGER treatment effect (Chan 2017 Results).",
        "Second member of the percent-predicted spirometry surface whose",
        "first member is the registered FEV1_PCTPRED; see the register entry",
        "Notes for the family reading under which it was added."
      ),
      source_name = "arm level baseline %predicted FVC (Chan 2017 Table 2 note)"
    ),

    # ------------------------------------------------------------------
    # Treatment-arm indicators, one per active treatment. Exactly one is 1
    # on an active arm; ALL FOURTEEN are 0 on a placebo arm, which makes
    # the treatment effect vanish -- correct, because the output is the
    # placebo-corrected difference and a placebo arm's value is zero by
    # construction.
    #
    # These are bare-INN names, following the register's explicit written
    # rule (recorded in the TRT_BENRALIZUMAB entry Notes): "use TRT_<INN>
    # when the drug is the study treatment under investigation,
    # CONMED_<INN> when it is background co-medication, and a bare INN only
    # for MBMA trial-arm indicators". They join NAPROXEN / TRAMADOL /
    # TAPENTADOL / ACETAMINOPHEN / DICLOFENAC in that family. The register
    # also carries a COMPETING precedent -- the Sandra 2024 anti-HBV siRNA
    # batch registered 13 MBMA arm indicators as TRT_<drug> -- and the two
    # have not been reconciled upstream; see each register entry's Notes.
    # ------------------------------------------------------------------
    NINTEDANIB = list(
      description = "Binary study-arm treatment indicator: 1 = the arm received nintedanib, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the arm received placebo or a different active treatment).",
      notes = "MBMA study-arm-level treatment indicator. Selects emax_nintedanib = 3.31 %predicted FVC. Chan 2017 Table 1: 3 arms, 21 timepoints, 723 subjects, 300 mg/day, from TOMORROW and INPULSIS 1/2 (refs 34-35). One of only two treatments whose 95% CI excluded zero. The Figure 4 forest plot labels the two INPULSIS arms 281.7 and 280.5 mg/day -- the exposure-weighted mean daily dose after interruptions and reductions -- while Table 1 and the model use the nominal 300 mg/day; the dose does not enter the model, only the indicator does.",
      source_name = "Nintedanib treatment (Chan 2017 Table 2)"
    ),
    PIRFENIDONE = list(
      description = "Binary study-arm treatment indicator: 1 = the arm received pirfenidone at either studied dose, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the arm received placebo or a different active treatment).",
      notes = "MBMA study-arm-level treatment indicator. The ONLY treatment in the analysis with an estimated dose-response, so it must be paired with DOSE_HIGH: this indicator selects the pirfenidone branch and DOSE_HIGH picks 2.42 (1,197 mg/day) or 3.87 (2,403 mg/day) within it. Chan 2017 Table 1: 4 arms, 22 timepoints, 710 subjects, from CAPACITY-1/2 and ASCEND (refs 20, 25). One of only two treatments whose 95% CI excluded zero, and the numerically largest effect among the approved drugs at its recommended 2,403 mg/day dose.",
      source_name = "Pirfenidone 1,197 / 2,403 mg/day treatment (Chan 2017 Table 2)"
    ),
    DOSE_HIGH = list(
      description = "Pirfenidone high-dose indicator: 1 = the arm received pirfenidone 2,403 mg/day, 0 = the arm received pirfenidone 1,197 mg/day (or is not a pirfenidone arm).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (pirfenidone 1,197 mg/day, the low dose; also 0 on every non-pirfenidone arm, where it is inert because the PIRFENIDONE indicator is 0).",
      notes = "Registered DOSE_HIGH semantics exactly: 'a step-function switch in a parameter at the top of the dose range'. Chan 2017 Results, 'Model development': 'the dose-response relationship of pirfenidone was characterized using a step-function, with the high dose (2,403 mg/day) having an estimated greater treatment effect than the low dose (1,197 mg/day)'. Only two pirfenidone dose levels exist in the analysis dataset, so no numerical mg threshold is needed or implied -- the indicator names the level, it does not bin a continuous dose. Only meaningful when PIRFENIDONE = 1. Chan 2017 attempted a dose-response for PRM-151 as well (1, 5 and 10 mg/kg) but the estimation failed to minimize, so PRM-151 carries a single Emax and does NOT use this column.",
      source_name = "Pirfenidone dose level (Chan 2017 Table 2 rows 'Pirfenidone 1,197 mg/day treatment' and 'Pirfenidone 2,403 mg/day treatment')"
    ),
    PRM151 = list(
      description = "Binary study-arm treatment indicator: 1 = the arm received PRM-151 (recombinant human pentraxin-2 / serum amyloid P; later zinpentraxin alfa), 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the arm received placebo or a different active treatment).",
      notes = "MBMA study-arm-level treatment indicator, named by the development code the source uses throughout. Selects emax_prm151 = 12.5 %predicted FVC -- numerically the largest effect in the analysis, but the least reliable: Chan 2017 Table 1 records only 3 arms, 3 timepoints and 15 subjects in total, at a single 8-week timepoint (trial PRM151F-12GL, ref 27), and the 95% CI is -5.13 to 30.1. Chan 2017 Discussion: 'the Emax estimate of this compound was not precise or reliable'. Three dose strengths (1, 5 and 10 mg/kg) were in the dataset but a dose-response could not be estimated -- 'attempts to model a dose-response relationship resulted in minimization failures' -- so all three dose arms share this one indicator and one Emax.",
      source_name = "PRM-151 treatment (Chan 2017 Table 2)"
    ),
    INTERFERON_GAMMA = list(
      description = "Binary study-arm treatment indicator: 1 = the arm received interferon gamma-1b, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the arm received placebo or a different active treatment).",
      notes = "MBMA study-arm-level treatment indicator. Selects emax_ifngamma = -0.618 %predicted FVC. Chan 2017 Table 1: 4 arms, 8 timepoints, 607 subjects, 200 ug (refs 18, 23, 32, 33) -- the most-studied treatment in the dataset by number of trials. Chan 2017 Discussion flags it as 'an obvious outlier': two trials reported a large positive effect and two reported none, with the large INSPIRE trial (ref 23) dominating the pooled estimate. Load-bearing beyond its own arms, because colchicine and prednisone had no placebo-controlled data and reach placebo only INDIRECTLY through small active-control trials against interferon gamma.",
      source_name = "Interferon-gamma treatment (Chan 2017 Table 2)"
    ),
    COLCHICINE = list(
      description = "Binary study-arm treatment indicator: 1 = the arm received colchicine, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the arm received placebo or a different active treatment).",
      notes = "MBMA study-arm-level treatment indicator. Selects emax_colchicine = -7.63 %predicted FVC, the second most negative effect in the analysis. Chan 2017 Table 1: 3 arms, 7 timepoints, 39 subjects, 1 mg (refs 18, 19, 33); the Figure 4 forest plot labels the Douglas 1998 arm 0.9 mg. INDIRECTLY estimated: Chan 2017 Discussion records that colchicine 'lacked placebo-controlled trial data in subjects with IPF' and was compared to placebo only through small active-control trials against interferon gamma, which itself showed no effect in INSPIRE -- hence the large negative point estimate and wide CI (-13 to -2.3).",
      source_name = "Colchicine treatment (Chan 2017 Table 2)"
    ),
    PREDNISONE = list(
      description = "Binary study-arm treatment indicator: 1 = the arm received prednisone, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the arm received placebo or a different active treatment).",
      notes = "MBMA study-arm-level treatment indicator. Selects emax_prednisone = -11.9 %predicted FVC, the most negative effect in the analysis. Chan 2017 Table 1: 1 arm, 1 timepoint, 12 subjects, 40 mg (ref 19). INDIRECTLY estimated through the same interferon-gamma active-control chain as colchicine (Chan 2017 Discussion), which is why the CI (-25.8 to 1.95) spans 28 %predicted FVC on 12 subjects. Distinct from the registered DRUG_PRED and PRED_DOSE columns, which concern prednisone/prednisolone as a comparator arm and as a concomitant daily dose respectively.",
      source_name = "Prednisone treatment (Chan 2017 Table 2)"
    ),
    SILDENAFIL = list(
      description = "Binary study-arm treatment indicator: 1 = the arm received sildenafil, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the arm received placebo or a different active treatment).",
      notes = "MBMA study-arm-level treatment indicator. Selects emax_sildenafil = 2.67 %predicted FVC (95% CI -7.15 to 12.5, not distinguishable from zero). Chan 2017 Table 1: 1 arm, 1 timepoint, 14 subjects, 20 mg (ref 21, trial CLIN-009-05F). Distinct from the registered CONMED_SILDENAFIL, which flags sildenafil as background co-medication; here it is the randomized study treatment.",
      source_name = "Sildenafil treatment (Chan 2017 Table 2)"
    ),
    BOSENTAN = list(
      description = "Binary study-arm treatment indicator: 1 = the arm received bosentan, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the arm received placebo or a different active treatment).",
      notes = "MBMA study-arm-level treatment indicator. Selects emax_bosentan = 1.75 %predicted FVC (95% CI -4.37 to 7.86). Chan 2017 Table 1: 1 arm, 1 timepoint, 74 subjects, 125 mg (ref 22, trial BUILD-1). An endothelin receptor antagonist: Chan 2017 Results records that this class -- ambrisentan, bosentan and macitentan -- had a HIGHER predicted residual SD than other classes in the separate SD-imputation model, which affects the study weights but is not a parameter of this model. Distinct from the registered CONMED_BOSENTAN (background co-medication).",
      source_name = "Bosentan treatment (Chan 2017 Table 2)"
    ),
    AMBRISENTAN = list(
      description = "Binary study-arm treatment indicator: 1 = the arm received ambrisentan, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the arm received placebo or a different active treatment).",
      notes = "MBMA study-arm-level treatment indicator. Selects emax_ambrisentan = -2.03 %predicted FVC (95% CI -4.18 to 0.126, the narrowest CI of any non-significant treatment). Chan 2017 Table 1: 1 arm, 7 timepoints, 330 subjects, 10 mg (ref 28, trial ARTEMIS-IPF) -- one of the eight treatments contributing longitudinal data. An endothelin receptor antagonist; see BOSENTAN on the class-level residual-SD finding.",
      source_name = "Ambrisentan treatment (Chan 2017 Table 2)"
    ),
    ACETYLCYSTEINE = list(
      description = "Binary study-arm treatment indicator: 1 = the arm received N-acetylcysteine, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the arm received placebo or a different active treatment).",
      notes = "MBMA study-arm-level treatment indicator, under the INN acetylcysteine (the source writes N-acetylcysteine). Selects emax_acetylcysteine = -0.0568 %predicted FVC, the effect closest to exactly zero in the analysis. Chan 2017 Table 1: 1 arm, 4 timepoints, 133 subjects, 1,800 mg/day (ref 24, trial PANTHER-IPF) -- one of the eight treatments contributing longitudinal data.",
      source_name = "N-acetylcysteine treatment (Chan 2017 Table 2)"
    ),
    WARFARIN = list(
      description = "Binary study-arm treatment indicator: 1 = the arm received warfarin, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the arm received placebo or a different active treatment).",
      notes = "MBMA study-arm-level treatment indicator. Selects emax_warfarin = -0.276 %predicted FVC (95% CI -4.89 to 4.33). Chan 2017 Table 1: 1 arm, 3 timepoints, 72 subjects, 1.8 mg (ref 26, trial ACE-IPF) -- one of the eight treatments contributing longitudinal data. The 1.8 mg figure is the trial's mean achieved daily dose under INR-guided titration, not a fixed assigned dose; the dose does not enter the model, only the indicator does.",
      source_name = "Warfarin treatment (Chan 2017 Table 2)"
    ),
    ETANERCEPT = list(
      description = "Binary study-arm treatment indicator: 1 = the arm received etanercept, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the arm received placebo or a different active treatment).",
      notes = "MBMA study-arm-level treatment indicator. Selects emax_etanercept = 1.47 %predicted FVC (95% CI -2.3 to 5.24). Chan 2017 Table 1: 1 arm, 5 timepoints, 46 subjects, 25 mg (ref 29, trial 0881A4-203) -- one of the eight treatments contributing longitudinal data, and the only tumor necrosis factor inhibitor in the dataset.",
      source_name = "Etanercept treatment (Chan 2017 Table 2)"
    ),
    AZATHIOPRINE = list(
      description = "Binary study-arm treatment indicator: 1 = the arm received azathioprine, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the arm received placebo or a different active treatment).",
      notes = "MBMA study-arm-level treatment indicator. Selects emax_azathioprine = 5.89 %predicted FVC (95% CI -4.14 to 15.9), the second largest point estimate but on 14 subjects. Chan 2017 Table 1: 1 arm, 1 timepoint, 14 subjects, '3 mg' (ref 30, Raghu 1991); the Figure 4 forest plot gives the same arm as 3 mg/kg, which is the clinically standard azathioprine dosing and is almost certainly what Table 1's unqualified '3 mg' abbreviates -- see the vignette Errata. The dose does not enter the model, only the indicator does.",
      source_name = "Azathioprine treatment (Chan 2017 Table 2)"
    ),
    COTRIMOXAZOLE = list(
      description = "Binary study-arm treatment indicator: 1 = the arm received co-trimoxazole (sulfamethoxazole / trimethoprim), 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the arm received placebo or a different active treatment).",
      notes = "MBMA study-arm-level treatment indicator. Selects emax_cotrimoxazole = 0.172 %predicted FVC (95% CI -3.72 to 4.07). Chan 2017 Table 1: 1 arm, 1 timepoint, 95 subjects, 960 mg (ref 31, Shulgina 2013) -- the 960 mg is the standard double-strength co-trimoxazole tablet (800 mg sulfamethoxazole + 160 mg trimethoprim), so the single number covers the fixed combination and no component split is needed.",
      source_name = "Co-trimoxazole treatment (Chan 2017 Table 2)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Arm-mean age at randomization.",
      units = "years",
      type = "continuous",
      notes = "One of the six covariates Chan 2017 prespecified and screened for an association with the treatment effect (Methods, 'Exploratory analysis'). Not retained: Results, 'Model development' states that 'only baseline %predicted FVC ... was included in the final model'. No coefficient is reported, so there is nothing to encode. Chan 2017 also used age as an INPUT to the multivariate regression that imputed missing baseline %predicted FVC and smoking status, but that imputation happened during dataset construction and is not part of the fitted model."
    ),
    DIS_DUR = list(
      description = "Arm-mean duration of idiopathic pulmonary fibrosis since diagnosis at randomization.",
      units = "years",
      type = "continuous",
      notes = "Prespecified and screened (Chan 2017 Methods, 'Exploratory analysis'); not retained in the final model and no coefficient is reported."
    ),
    SEXF = list(
      description = "Arm-level proportion of subjects who were female, derived as 1 minus the reported male proportion.",
      units = "fraction 0-1",
      type = "continuous",
      notes = "Chan 2017 prespecified 'proportions of subjects who were men' and screened it (Methods, 'Exploratory analysis'); not retained and no coefficient is reported. Recorded here under the canonical female-coded SEXF so the screen's provenance survives, with the value inversion noted; the direction is untestable because no estimate exists."
    ),
    RACE_WHITE_PCT = list(
      description = "Arm-level percentage of subjects who were white.",
      units = "% of subjects in the arm (0-100)",
      type = "continuous",
      notes = "Prespecified and screened (Chan 2017 Methods, 'Exploratory analysis'); not retained and no coefficient is reported."
    ),
    SMOKE = list(
      description = "Arm-level proportion of subjects who were current or former smokers.",
      units = "fraction 0-1",
      type = "continuous",
      notes = "Prespecified and screened (Chan 2017 Methods, 'Exploratory analysis'); not retained and no coefficient is reported. Chan 2017 imputed missing arm-level smoking status by multivariate regression during dataset construction, and used it alongside age and region as an input to the baseline-%predicted-FVC imputation, but neither step is part of the fitted model. Note the pooled current-OR-former coding, which is broader than the current-smoker-only coding the canonical SMOKE carries in the COPD models."
    ),
    N_ARM = list(
      description = "Number of subjects contributing to the study arm.",
      units = "(count of subjects)",
      type = "continuous",
      notes = "Not a covariate on any structural parameter, but load-bearing for the residual: Chan 2017 Eq. 1 defines the within-arm residual as normally distributed 'with a variance dependent on the sample size and observed SD of the trial arm', and the fit used the inverse of that estimated variance as the study weight. Documented here rather than in covariateData because the model carries the UNIT-WEIGHT residual scale (addSd) and leaves the per-arm reweighting to downstream code -- the same convention as every other MBMA in this library. See the ini() note on addSd and the vignette Assumptions and deviations section."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 4919L,
    n_studies = 20L,
    n_arms = 43L,
    n_observations = 148L,
    age_range = "not tabulated at arm level; IPF most commonly occurs over age 50 (Chan 2017 Introduction). Arm-mean age was prespecified as a covariate and screened but not retained.",
    weight_range = "not available -- Chan 2017 Results records that weight, height and body mass index were NOT imputed because more than 60% of trial arms had missing values",
    sex_female_pct = NA_real_,
    race_ethnicity = "not tabulated at arm level; the arm-level white proportion was screened as a covariate and not retained",
    disease_state = "idiopathic pulmonary fibrosis",
    baseline_endpoint = "baseline %predicted FVC, approximate median 74% across the analysis dataset (the model's normalizing constant). Observed arm-level values named in the paper: 76.2% (CAPACITY-2 / PIPF-004) and 74.2% (TOMORROW) from the Figure 3 panel annotations; the stratified low / high baseline group medians used for Figure 5 are given as 66% and 76% in the Results text and as 67 and 76 in the Figure 5 legend.",
    dose_range = "one fixed regimen per treatment except pirfenidone (1,197 and 2,403 mg/day) and PRM-151 (1, 5, 10 mg/kg, pooled into one effect). Chan 2017 Table 1: ambrisentan 10 mg, azathioprine 3 mg (3 mg/kg per Figure 4), bosentan 125 mg, colchicine 1 mg (0.9 mg per Figure 4), co-trimoxazole 960 mg, etanercept 25 mg, interferon gamma 200 ug, N-acetylcysteine 1,800 mg/day, nintedanib 300 mg/day, prednisone 40 mg, sildenafil 20 mg, warfarin 1.8 mg",
    regions = "Europe, North America and Australia; Chan 2017 Results notes that none of the 20 trials was conducted in Asia",
    treatments = "placebo plus 14 active treatments in 12 drug classes, giving 15 treatment regimens once the pirfenidone dose step is counted: ambrisentan, azathioprine, bosentan, colchicine, co-trimoxazole, etanercept, interferon gamma, N-acetylcysteine, nintedanib, pirfenidone (2 doses), prednisone, PRM-151, sildenafil, warfarin",
    timepoints = "148 arm-timepoints; treatment durations 8-104 weeks. Eight treatments contributed longitudinal data (etanercept, warfarin, interferon gamma, ambrisentan, nintedanib, N-acetylcysteine, pirfenidone, colchicine); the rest contributed a single timepoint each.",
    notes = paste(
      "MBMA at the study-arm level: each modelled observation is one",
      "published trial arm's MEAN change from baseline %predicted FVC at one",
      "timepoint, inverse-variance weighted by arm size and arm SD. The",
      "analysis dataset was selected from a larger 40-citation / 32-trial",
      "database built by a PRISMA-style systematic review of PubMed plus the",
      "FDA and ClinicalTrials.gov websites (search run September 2015);",
      "17 of the 20 analysis trials were double-blinded and placebo-controlled",
      "and none had stratified data within arms. Placebo contributed the",
      "largest single block: 17 arms, 63 timepoints, 2,035 subjects. 38% of",
      "records in the augmented database had a missing SD, imputed from a",
      "separate SD model with linear time and drug class as predictors",
      "(observed-vs-predicted r = 0.70); a sensitivity analysis using a",
      "mixture of observed and predicted SDs left the predictions 'virtually",
      "indistinguishable'. The model predicts study-arm-mean PLACEBO-CORRECTED",
      "responses and is NOT suitable for individual-patient simulation, nor",
      "for absolute %predicted FVC -- see the description and the vignette."
    )
  )

  ini({
    # ================================================================
    # Treatment maximum effects (Chan 2017 Table 2, 'Estimate (95% CI)'
    # column; every value is duplicated in the Figure 2 ranking plot and
    # the two agree exactly, which is the cross-check used for this
    # block). Each is the maximal change in %predicted FVC relative to
    # placebo, at the reference baseline of 74% predicted FVC and at
    # infinite time.
    #
    # These are on the LINEAR scale, not log-transformed: six of the
    # fifteen are NEGATIVE, so no log parameterization is possible. They
    # are NOT constrained positive and a re-fit must not constrain them.
    #
    # Listed in Table 2 order.
    # ================================================================
    emax_ifngamma <- -0.618
    label("Maximum interferon gamma-1b effect on placebo-corrected change from baseline %predicted FVC (% predicted)")  # Chan 2017 Table 2 'Interferon-gamma treatment' = -0.618 (95% CI -2.46 to 1.23); Figure 2 'ifn gamma' = -0.618 [-2.46;1.23]

    emax_colchicine <- -7.63
    label("Maximum colchicine effect on placebo-corrected change from baseline %predicted FVC (% predicted)")  # Chan 2017 Table 2 'Colchicine treatment' = -7.63 (95% CI -13 to -2.3); Figure 2 'colchicine' = -7.63 [-13;-2.3]. Estimated indirectly via active-control trials against interferon gamma

    emax_prednisone <- -11.9
    label("Maximum prednisone effect on placebo-corrected change from baseline %predicted FVC (% predicted)")  # Chan 2017 Table 2 'Prednisone treatment' = -11.9 (95% CI -25.8 to 1.95); Figure 2 'prednisone' = -11.9 [-25.8;1.95]. Estimated indirectly via active-control trials against interferon gamma, on 12 subjects

    emax_sildenafil <- 2.67
    label("Maximum sildenafil effect on placebo-corrected change from baseline %predicted FVC (% predicted)")  # Chan 2017 Table 2 'Sildenafil treatment' = 2.67 (95% CI -7.15 to 12.5); Figure 2 'sildenafil' = 2.67 [-7.15;12.5]

    emax_bosentan <- 1.75
    label("Maximum bosentan effect on placebo-corrected change from baseline %predicted FVC (% predicted)")  # Chan 2017 Table 2 'Bosentan treatment' = 1.75 (95% CI -4.37 to 7.86); Figure 2 'bosentan' = 1.75 [-4.37;7.86]

    emax_ambrisentan <- -2.03
    label("Maximum ambrisentan effect on placebo-corrected change from baseline %predicted FVC (% predicted)")  # Chan 2017 Table 2 'Ambrisentan treatment' = -2.03 (95% CI -4.18 to 0.126); Figure 2 'ambrisentan' = -2.03 [-4.18;0.126]

    emax_acetylcysteine <- -0.0568
    label("Maximum N-acetylcysteine effect on placebo-corrected change from baseline %predicted FVC (% predicted)")  # Chan 2017 Table 2 'N-acetylcysteine treatment' = -0.0568 (95% CI -2.01 to 1.9); Figure 2 'n-acetylcysteine' = -0.0568 [-2.01;1.9]

    emax_warfarin <- -0.276
    label("Maximum warfarin effect on placebo-corrected change from baseline %predicted FVC (% predicted)")  # Chan 2017 Table 2 'Warfarin treatment' = -0.276 (95% CI -4.89 to 4.33); Figure 2 'warfarin' = -0.276 [-4.89;4.33]

    emax_etanercept <- 1.47
    label("Maximum etanercept effect on placebo-corrected change from baseline %predicted FVC (% predicted)")  # Chan 2017 Table 2 'Etanercept treatment' = 1.47 (95% CI -2.3 to 5.24); Figure 2 'etanercept' = 1.47 [-2.3;5.24]

    emax_azathioprine <- 5.89
    label("Maximum azathioprine effect on placebo-corrected change from baseline %predicted FVC (% predicted)")  # Chan 2017 Table 2 'Azathioprine treatment' = 5.89 (95% CI -4.14 to 15.9); Figure 2 'azathioprine' = 5.89 [-4.14;15.9]

    emax_cotrimoxazole <- 0.172
    label("Maximum co-trimoxazole effect on placebo-corrected change from baseline %predicted FVC (% predicted)")  # Chan 2017 Table 2 'Co-trimoxazole treatment' = 0.172 (95% CI -3.72 to 4.07); Figure 2 'co-trimoxazole' = 0.172 [-3.72;4.07]

    emax_nintedanib <- 3.31
    label("Maximum nintedanib 300 mg/day effect on placebo-corrected change from baseline %predicted FVC (% predicted)")  # Chan 2017 Table 2 'Nintedanib treatment' = 3.31 (95% CI 2.15 to 4.47); Figure 2 'nintedanib' = 3.31 [2.15;4.47]. One of only two CIs excluding zero

    emax_pirfenidone_low <- 2.42
    label("Maximum pirfenidone 1,197 mg/day effect on placebo-corrected change from baseline %predicted FVC (% predicted)")  # Chan 2017 Table 2 'Pirfenidone 1,197 mg/day treatment' = 2.42 (95% CI 0.733 to 4.11); Figure 2 'pirfenidone 1197 mg/day' = 2.42 [0.733;4.11]

    emax_pirfenidone_high <- 3.87
    label("Maximum pirfenidone 2,403 mg/day effect on placebo-corrected change from baseline %predicted FVC (% predicted)")  # Chan 2017 Table 2 'Pirfenidone 2,403 mg/day treatment' = 3.87 (95% CI 2.68 to 5.06); Figure 2 'pirfenidone 2403 mg/day' = 3.87 [2.68;5.06]. Largest effect among the two approved drugs

    emax_prm151 <- 12.5
    label("Maximum PRM-151 effect on placebo-corrected change from baseline %predicted FVC (% predicted)")  # Chan 2017 Table 2 'PRM-151 treatment' = 12.5 (95% CI -5.13 to 30.1); Figure 2 'prm-151' = 12.5 [-5.13;30.1]. Largest point estimate but on 15 subjects at one 8-week timepoint; Discussion calls it 'not precise or reliable'

    # ================================================================
    # Shared treatment-effect time course (Chan 2017 Results, 'Model
    # development', and the Table 2 note): the fraction of the maximum
    # effect attained by time t is 1 - exp(-kel * time), with a SINGLE
    # shared kel across all fifteen regimens. Chan 2017 tested
    # treatment-specific and class-specific time courses and reports that
    # they 'did not result in successful minimization of the model run,
    # potentially due to limited longitudinal data'.
    #
    # lambda IS ALREADY ON THE LOG SCALE, so this ini() value is the
    # tabulated number itself and NOT log() of it. The paper never says
    # so, and the Table 2 note's bare '1-exp(-lambda*time)' reads as if
    # lambda were the rate. Three checks settle it:
    #   (a) a rate of -3.18 /week makes 1-exp(+3.18*t) diverge to -Inf,
    #       which is not a treatment-effect time course at all;
    #   (b) exp(-3.18) = 0.041674 /week gives ln(2)/kel = 16.6 weeks and
    #       ln(10)/kel = 55.3 weeks, matching the paper's own '50% maximum
    #       efficacy of 16.5 weeks' and '90% ... at 56 weeks' (Results);
    #   (c) a rate of 3.18 /week would put 50% of maximum at 0.22 weeks,
    #       contradicting the same sentence by a factor of 75.
    # ================================================================
    lkel <- -3.18
    label("Log first-order rate constant of the shared approach to the maximum treatment effect (log of 1/week; kel = 0.0417 /week, 50% of maximum at 16.6 weeks and 90% at 55.3 weeks)")  # Chan 2017 Table 2 'Time-varying effect (lambda)' = -3.18 (95% CI -3.69 to -2.68); Results quotes t50 = 16.5 weeks and t90 at about 56 weeks, which back-solve kel = exp(lambda)

    # ================================================================
    # Baseline-lung-function effect on the treatment effect (Chan 2017
    # Table 2 note): the whole treatment effect is multiplied by
    # (arm level baseline %predicted FVC / 74)^embase. A positive
    # exponent means a higher baseline %predicted FVC predicts a LARGER
    # treatment effect, which is the direction the Results paragraph
    # describes.
    #
    # VALUE CONFLICT IN THE SOURCE, resolved in favour of Table 2. Table
    # 2 row 'Baseline FVC effect (embase)' gives 1.48 (95% CI -2.7 to
    # 5.66); the Results covariate paragraph instead gives 3.86 (95% CI
    # -0.631 to 8.35). Both are internally well formed, so the CI cannot
    # adjudicate. Chan 2017 Figure 5 does, and decisively -- it plots the
    # simulated nintedanib and pirfenidone time courses at two baselines
    # (legend FVC=67 and FVC=76), so the ratio between its two curves is
    # (76/67)^embase with everything else cancelling:
    #   embase = 1.48 predicts a ratio of 1.21
    #   embase = 3.86 predicts a ratio of 1.63
    #   Figure 5 reads about 3.0 / 2.4 = 1.25 for nintedanib and
    #   3.5 / 2.85 = 1.23 for pirfenidone.
    # The absolute levels agree too: nintedanib at FVC=76 and 52 weeks is
    # 3.31 * 0.885 * (76/74)^1.48 = 3.05 against a read of about 3.0,
    # and the dotted 95% prediction interval reads about 2.0 to 4.05
    # against the Emax CI scaled the same way, 1.98 to 4.12. Chan 2017
    # Figure 3 confirms it independently on a different baseline: at
    # week 24 with the annotated base=76.2, the pirfenidone 2,403 mg/day
    # fitted curve sits 2.55 %predicted FVC above the placebo curve and
    # the model gives 3.87 * 0.632 * 1.044 = 2.555.
    # See the vignette Errata; the gates are reproduced there.
    # ================================================================
    e_fvc_pctpred_emax <- 1.48
    label("Power exponent on (FVC_PCTPRED / 74) scaling every treatment effect (unitless)")  # Chan 2017 Table 2 'Baseline FVC effect (embase)' = 1.48 (95% CI -2.7 to 5.66); the Results text's conflicting 3.86 is falsified by Figures 3 and 5 -- see the block comment above

    # ================================================================
    # Between-trial-arm variability (Chan 2017 Eq. 1 eta_ij, 'the random
    # residual due to between arm variability'; Table 2 last row). This
    # is an MBMA ARM-level random effect -- one draw per published trial
    # arm, not per patient -- and it enters Eq. 1 ADDITIVELY on the
    # response rather than on any structural parameter, which is why it
    # is named for the output and not for a parameter.
    #
    # The estimate is 1.07e-6, i.e. numerically indistinguishable from
    # zero: on the variance reading its SD is 0.00103 %predicted FVC,
    # and on the SD reading 1.07e-6 %predicted FVC. Chan 2017 does not
    # say which scale Table 2 reports and its SE 'was not estimated', so
    # the two readings cannot be separated -- but they differ by three
    # orders of magnitude on a quantity that is itself three to six
    # orders below the effects it perturbs, so the choice has no
    # practical consequence. Encoded as the VARIANCE, per the nlmixr2
    # '~' convention. Retained rather than dropped because the paper
    # reports it; see the vignette Errata.
    # ================================================================
    eta_arm_fvcppcfb ~ 1.07e-6  # Chan 2017 Table 2 'Between trial-arm variability' = 1.07 x 10^-6; the SE was not estimated and the scale is not stated, so this is read as a variance

    # ================================================================
    # Residual. Chan 2017 Eq. 1's within-arm residual epsilon_ijt is
    # 'assumed to be normally distributed with a variance dependent on
    # the sample size and observed SD of the trial arm', and 'the model
    # used the inverse of the estimated variance as weights'. That
    # variance is therefore a property of each published arm -- its size
    # and its observed or SD-model-predicted standard deviation -- and
    # not a parameter of the model. Chan 2017 tabulates no sigma at all.
    #
    # Encoded as the UNIT-WEIGHT residual scale, the convention used by
    # every other MBMA in this library (Boucher 2018 naproxen, Vargo 2014
    # statins, Sato 2024 SGLT2, Yao 2023 SGLT2).
    #
    # USABILITY TRAP, and it is a large one: a bare stochastic rxSolve()
    # applies a 1.0 %predicted-FVC additive residual to every arm, which
    # is far too small for a 12-subject arm and too large for a
    # 550-subject one. The operative per-observation SD is the arm's own
    # SD / sqrt(N_ARM). Anything reproducing the paper's fitted curves
    # should use rxode2::zeroRe(mod, 'sigma'), which keeps the arm-level
    # omega and drops the residual -- the vignette does exactly that.
    # ================================================================
    addSd <- fixed(1)
    label("Residual SD at unit study weight (% predicted); the source scales it per record by each arm's standard error, the arm SD divided by the square root of the arm size, which is supplied downstream rather than in ini()")  # Chan 2017 Eq. 1: epsilon_ijt has 'a variance dependent on the sample size and observed SD of the trial arm'; no sigma is tabulated
  })

  model({
    # ================================================================
    # Study-arm covariates supplied per row:
    #   FVC_PCTPRED       -- arm-mean baseline %predicted FVC
    #   NINTEDANIB, PIRFENIDONE, PRM151, INTERFERON_GAMMA, COLCHICINE,
    #   PREDNISONE, SILDENAFIL, BOSENTAN, AMBRISENTAN, ACETYLCYSTEINE,
    #   WARFARIN, ETANERCEPT, AZATHIOPRINE, COTRIMOXAZOLE
    #                     -- treatment indicators; exactly one is 1 on an
    #                        active arm, all are 0 on a placebo arm
    #   DOSE_HIGH         -- 1 for pirfenidone 2,403 mg/day, 0 for
    #                        1,197 mg/day; inert unless PIRFENIDONE is 1
    #
    # This model has NO ODE states and NO dose events: Chan 2017 Eq. 1 is
    # an algebraic regression in trial time, fitted with nlme::gnls.
    # ================================================================

    # Back-transform the shared time-course rate. lambda is tabulated on
    # the log scale; see the ini() block comment for the three checks.
    kel <- exp(lkel)

    # ----------------------------------------------------------------
    # Treatment selection. Chan 2017 estimated a separate Emax per
    # treatment rather than per drug class, noting that with 12 classes
    # and 14 drugs 'estimating the maximal effect of each treatment
    # instead of each treatment class did not drastically increase the
    # complexity of the model'.
    #
    # Pirfenidone is the only treatment with a dose step, so its two
    # estimates are selected by DOSE_HIGH within the PIRFENIDONE branch.
    # PRM-151 pools its 1, 5 and 10 mg/kg arms into one Emax because the
    # dose-response would not minimize.
    #
    # On a placebo arm every indicator is 0, so emax_arm is 0 and the
    # whole treatment effect vanishes -- which is correct, because the
    # modelled quantity is the placebo-corrected difference.
    # ----------------------------------------------------------------
    emax_arm <-
      emax_ifngamma * INTERFERON_GAMMA +
      emax_colchicine * COLCHICINE +
      emax_prednisone * PREDNISONE +
      emax_sildenafil * SILDENAFIL +
      emax_bosentan * BOSENTAN +
      emax_ambrisentan * AMBRISENTAN +
      emax_acetylcysteine * ACETYLCYSTEINE +
      emax_warfarin * WARFARIN +
      emax_etanercept * ETANERCEPT +
      emax_azathioprine * AZATHIOPRINE +
      emax_cotrimoxazole * COTRIMOXAZOLE +
      emax_nintedanib * NINTEDANIB +
      emax_prm151 * PRM151 +
      emax_pirfenidone_low * PIRFENIDONE * (1 - DOSE_HIGH) +
      emax_pirfenidone_high * PIRFENIDONE * DOSE_HIGH

    # ----------------------------------------------------------------
    # Shared empirical time course (Chan 2017 Table 2 note). Rises from
    # 0 at randomization to 1 at infinite time; 0.5 at 16.6 weeks.
    # ----------------------------------------------------------------
    tcourse <- 1 - exp(-kel * time)

    # ----------------------------------------------------------------
    # Baseline-lung-function scaling (Chan 2017 Table 2 note):
    # (arm level baseline %predicted FVC / 74)^embase. Equals 1 at the
    # 74% reference, so an arm at the dataset median realises exactly the
    # tabulated Emax.
    # ----------------------------------------------------------------
    fvcfactor <- (FVC_PCTPRED / 74)^e_fvc_pctpred_emax

    # ----------------------------------------------------------------
    # Chan 2017 Eq. 1 without its nonparametric placebo term E0_it,
    # which the paper leaves as one free estimate per trial per timepoint
    # and does not tabulate. What remains is exactly the paper's
    # placebo-corrected change from baseline -- the quantity plotted in
    # Figures 2, 4 and 5 and the difference between the fitted curves of
    # Figure 3. eta_arm_fvcppcfb is Eq. 1's additive between-arm term.
    # ----------------------------------------------------------------
    fvcppcfb <- emax_arm * tcourse * fvcfactor + eta_arm_fvcppcfb

    # Unit-weight additive residual; see the ini() note. Reproducing the
    # paper's fitted curves needs rxode2::zeroRe(mod, "sigma").
    fvcppcfb ~ add(addSd)
  })
}
