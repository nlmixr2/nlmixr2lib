Guo_2025_glp1ReceptorAgonists_mbma <- function() {
  description <- paste0(
    "MBMA. Model-based meta-analysis of the placebo-adjusted weight-reduction ",
    "time course of 12 glucagon-like peptide-1 receptor agonists (GLP-1RAs) in ",
    "adults with overweight or obesity, fit to arm-level summary data digitised ",
    "from 55 randomised, double-blind, placebo-controlled trials (16,269 ",
    "participants). The endpoint is the 'pure' drug effect: the change from ",
    "baseline in body weight in the drug arm minus the change from baseline in ",
    "the concurrent placebo arm (the paper's ddWeight), which removes the ",
    "trial-specific diet and exercise background. The time course is a single ",
    "exponential approach to an asymptote, E(t) = Emax * (1 - exp(-k * t)), with ",
    "a drug-specific Emax and a single shared onset rate k = 0.0578 /week (the ",
    "time to half of Emax is 0.693/k = 12.0 weeks). Six drugs (cotadutide, ",
    "danuglipron, JNJ-64565111, retatrutide, orforglipron and injectable ",
    "semaglutide) had an estimable dose-response and carry an Emax-in-dose term ",
    "Emax * Dose / (ED50 + Dose) with the ED50 fixed at the published value; the ",
    "other six (BI 456906, exenatide, liraglutide, oral semaglutide, tirzepatide ",
    "and mazdutide) were studied over too narrow a dose range for a ",
    "dose-response to be estimable and carry a flat per-arm Emax. Mean trial age ",
    "is the only retained covariate and acts exponentially on Emax, centred on ",
    "the 53.6-year median of the trial means; baseline weight, baseline BMI and ",
    "male ratio were screened and not retained. Between-STUDY (not ",
    "between-subject) variability is carried as study-level etas on Emax and on ",
    "log k; the residual is additive at unit study weight and the paper weights ",
    "it by 1/sqrt(N) for an arm of N participants. Suitable simulation scope is ",
    "the arm-mean placebo-adjusted weight-reduction time course; the model is ",
    "NOT suitable for individual-subject simulation. Parameter values are ",
    "Supplementary Table S6 (NONMEM 7.4, FOCEI)."
  )

  reference <- paste(
    "Guo H, Yang J, Huang J, Xu L, Lv Y, Wang Y, Ren J, Feng Y, Zheng Q, Li L.",
    "Comparative efficacy and safety of GLP-1 receptor agonists for weight",
    "reduction: A model-based meta-analysis of placebo-controlled trials.",
    "Obes Pillars. 2025 Feb 6;13:100162.",
    "doi:10.1016/j.obpill.2025.100162.",
    sep = " "
  )
  vignette <- "Guo_2025_glp1ReceptorAgonists_mbma"
  units <- list(
    time = paste0(
      "week (time since randomisation; the paper reports k in 1/week and ",
      "simulates the time course over 0-52 weeks)"
    ),
    dosing = paste0(
      "mg (nominal assigned dose level of the arm's GLP-1RA, supplied as the ",
      "DOSE_<DRUG>_MG covariate columns; the dosing interval differs by drug ",
      "-- once daily for the oral small molecules, once weekly or once daily ",
      "for the peptides -- and the model works from the assigned dose level, ",
      "not a dose rate. This MBMA does not consume rxode2 dose events)"
    ),
    concentration = paste0(
      "kg/kg (arm-mean placebo-adjusted change from baseline in body weight, ",
      "in kg; negative values are weight LOSS, so Cc = -7.57 means the drug ",
      "arm lost 7.57 kg more than its placebo arm. The observation is NOT a ",
      "drug concentration; the slash in the unit string satisfies ",
      "checkModelConventions parsing)"
    )
  )

  covariateData <- list(
    AGE = list(
      description = paste0(
        "Arm-mean (trial-mean) participant age. This is a study-arm-level ",
        "aggregate covariate, not an individual subject age."
      ),
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "The single covariate retained in the final model. Enters Emax ",
        "exponentially, Emax = Emax_typical * exp(-0.0304 * (AGE - 53.6)) ",
        "(Guo 2025 Equation 7), so older trial cohorts lose less weight. The ",
        "centring constant 53.6 years is the median of the 55 trial mean ages ",
        "(Guo 2025 Results 3.1); trial mean ages ranged 29.5-64.7 years, which ",
        "is the domain over which the coefficient is calibrated. The paper's ",
        "own check of this term: raising mean age from 40 to 50 years gives ",
        "exp(-0.0304 * 10) = 0.738, i.e. the 26.2 percent decrease in Emax ",
        "quoted in the Abstract and in Results 3.2."
      ),
      source_name = "Age (Guo 2025 Equation 7 / Supplementary Table S6 'thetaAge on Emax')"
    ),
    DOSE_LIRAGLUTIDE_MG = list(
      description = "Per-arm assigned liraglutide dose (mg per once-daily SC injection; 0 if the arm did not receive liraglutide).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "MBMA study-arm-level covariate. Liraglutide had NO estimable ",
        "dose-response over the studied range, so the model uses this column ",
        "only as a presence indicator (DOSE_LIRAGLUTIDE_MG > 0) selecting the ",
        "flat Emax of -4.25 kg; the numeric value does not change the ",
        "prediction. Studied range 1.2-3 mg (Guo 2025 Discussion); the paper's ",
        "52-week simulation tier quotes 1.8-3 mg (Results 3.3). 18 trials."
      ),
      source_name = "Liraglutide dose (Guo 2025 Results 3.1 / Supplementary Table S3)"
    ),
    DOSE_SEMAGLUTIDE_INJ_MG = list(
      description = "Per-arm assigned injectable (subcutaneous) semaglutide dose (mg per once-weekly injection; 0 if the arm did not receive injectable semaglutide).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "MBMA study-arm-level covariate. Injectable semaglutide HAS an ",
        "estimable dose-response: Emax_dose = -11.7 * Dose / (0.384 + Dose) ",
        "(Guo 2025 Equation 6). Studied doses 0.05-2.4 mg; the paper's worked ",
        "example uses 1.0 mg (Results 3.3). Route-qualified to keep it distinct ",
        "from DOSE_SEMAGLUTIDE_PO_MG, which has its own Emax (-5.36 kg): the ",
        "two formulations are separate drugs in this analysis. 11 trials."
      ),
      source_name = "Semaglutide (INJ) dose (Guo 2025 Equation 6 / Supplementary Table S6)"
    ),
    DOSE_SEMAGLUTIDE_PO_MG = list(
      description = "Per-arm assigned oral semaglutide dose (mg per once-daily tablet; 0 if the arm did not receive oral semaglutide).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "MBMA study-arm-level covariate. Oral semaglutide had NO estimable ",
        "dose-response, so the model uses this column only as a presence ",
        "indicator selecting the flat Emax of -5.36 kg. Studied range 1-40 mg ",
        "(Guo 2025 Results 3.3). 2 trials."
      ),
      source_name = "Semaglutide (P.O) dose (Guo 2025 Results 3.1 / Supplementary Table S3)"
    ),
    DOSE_EXENATIDE_MG = list(
      description = "Per-arm assigned exenatide dose (mg per administration; 0 if the arm did not receive exenatide).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "MBMA study-arm-level covariate. Exenatide had NO estimable ",
        "dose-response, so the model uses this column only as a presence ",
        "indicator selecting the flat Emax of -6.05 kg. Studied range ",
        "0.01-2 mg (Guo 2025 Results 3.3), spanning both the twice-daily ",
        "immediate-release and the once-weekly extended-release products, ",
        "which this analysis does not distinguish. 3 trials."
      ),
      source_name = "Exenatide dose (Guo 2025 Results 3.1 / Supplementary Table S3)"
    ),
    DOSE_DANUGLIPRON_MG = list(
      description = "Per-arm assigned danuglipron dose (mg per administration; 0 if the arm did not receive danuglipron).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "MBMA study-arm-level covariate. Danuglipron HAS an estimable ",
        "dose-response: Emax_dose = -9.29 * Dose / (80 + Dose) (Guo 2025 ",
        "Equation 2). Maximum administered dose across the included trials was ",
        "200 mg, which reaches 66.4 percent of Emax (Guo 2025 Discussion); the ",
        "52-week simulation tier quotes 100 mg (Results 3.3)."
      ),
      source_name = "Danuglipron dose (Guo 2025 Equation 2 / Supplementary Table S6)"
    ),
    DOSE_ORFORGLIPRON_MG = list(
      description = "Per-arm assigned orforglipron dose (mg per once-daily tablet; 0 if the arm did not receive orforglipron).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "MBMA study-arm-level covariate. Orforglipron HAS an estimable ",
        "dose-response: Emax_dose = -14.7 * Dose / (14.6 + Dose) (Guo 2025 ",
        "Equation 5). Maximum administered dose 45 mg (75.2 percent of Emax); ",
        "the 52-week simulation tier quotes 24 mg (Guo 2025 Results 3.3). ",
        "2 trials."
      ),
      source_name = "Orforglipron dose (Guo 2025 Equation 5 / Supplementary Table S6)"
    ),
    DOSE_TIRZEPATIDE_MG = list(
      description = "Per-arm assigned tirzepatide dose (mg per once-weekly SC injection; 0 if the arm did not receive tirzepatide).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "MBMA study-arm-level covariate. Tirzepatide had NO estimable ",
        "dose-response, so the model uses this column only as a presence ",
        "indicator selecting the flat Emax of -12.9 kg. Studied range 5-15 mg ",
        "(Guo 2025 Results 3.3). 4 trials."
      ),
      source_name = "Tirzepatide dose (Guo 2025 Results 3.1 / Supplementary Table S3)"
    ),
    DOSE_COTADUTIDE_MG = list(
      description = "Per-arm assigned cotadutide dose (mg per once-daily SC injection; 0 if the arm did not receive cotadutide).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "MBMA study-arm-level covariate. Cotadutide HAS an estimable ",
        "dose-response: Emax_dose = -10.5 * Dose / (0.219 + Dose) (Guo 2025 ",
        "Equation 1). Maximum administered dose 0.6 mg (70.8 percent of Emax); ",
        "the 52-week simulation tier quotes 0.2 mg (Guo 2025 Results 3.3). ",
        "4 trials."
      ),
      source_name = "Cotadutide dose (Guo 2025 Equation 1 / Supplementary Table S6)"
    ),
    DOSE_MAZDUTIDE_MG = list(
      description = "Per-arm assigned mazdutide dose (mg per once-weekly SC injection; 0 if the arm did not receive mazdutide).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "MBMA study-arm-level covariate. Mazdutide had NO estimable ",
        "dose-response, so the model uses this column only as a presence ",
        "indicator selecting the flat Emax of -7.75 kg. Studied range 3-10 mg ",
        "(Guo 2025 Results 3.3). 3 trials."
      ),
      source_name = "Mazdutide dose (Guo 2025 Results 3.1 / Supplementary Table S3)"
    ),
    DOSE_BI456906_MG = list(
      description = "Per-arm assigned BI 456906 (survodutide) dose (mg per once-weekly SC injection; 0 if the arm did not receive BI 456906).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "MBMA study-arm-level covariate. BI 456906 had NO estimable ",
        "dose-response, so the model uses this column only as a presence ",
        "indicator selecting the flat Emax of -13.5 kg. Studied range ",
        "1.8-4.8 mg (Guo 2025 Results 3.3). The development code BI 456906 is ",
        "retained in the column name because that is the name used throughout ",
        "Guo 2025; the INN survodutide was assigned later. 2 trials."
      ),
      source_name = "BI 456906 dose (Guo 2025 Results 3.1 / Supplementary Table S3)"
    ),
    DOSE_JNJ64565111_MG = list(
      description = "Per-arm assigned JNJ-64565111 dose (mg per once-weekly SC injection; 0 if the arm did not receive JNJ-64565111).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "MBMA study-arm-level covariate. JNJ-64565111 HAS an estimable ",
        "dose-response: Emax_dose = -18.6 * Dose / (6.73 + Dose) (Guo 2025 ",
        "Equation 3). Maximum administered dose 10 mg, reaching only 65.8 ",
        "percent of Emax -- the paper singles this drug out as the one whose ",
        "effect could be materially improved by dose escalation (Discussion). ",
        "The 52-week simulation tier quotes 7.4 mg (Results 3.3). Guo 2025 ",
        "Equations 3 and the Table S6 row spell the code 'JNJ-6456111' and ",
        "'JNJ-65465111' in places; the compound is JNJ-64565111, the spelling ",
        "used in Results 3.1 and 3.2. 2 trials."
      ),
      source_name = "JNJ-64565111 dose (Guo 2025 Equation 3 / Supplementary Table S6)"
    ),
    DOSE_RETATRUTIDE_MG = list(
      description = "Per-arm assigned retatrutide dose (mg per once-weekly SC injection; 0 if the arm did not receive retatrutide).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "MBMA study-arm-level covariate. Retatrutide HAS an estimable ",
        "dose-response: Emax_dose = -22.6 * Dose / (4 + Dose) (Guo 2025 ",
        "Equation 4). Retatrutide is the only GLP-1/GIP/GCG tri-agonist in the ",
        "analysis and carries the largest Emax of the 12 drugs. Maximum ",
        "administered dose 12 mg (60.7 percent of Emax); the 52-week simulation ",
        "tier quotes 6.5 mg (Guo 2025 Results 3.3). 2 trials."
      ),
      source_name = "Retatrutide dose (Guo 2025 Equation 4 / Supplementary Table S6)"
    )
  )

  # Covariates the paper screened in the forward-backward covariate search but
  # did NOT retain in the final model. Guo 2025 Results 3.1 and the Discussion
  # report that no significant effect of baseline weight, baseline BMI, or male
  # ratio on the weight-reduction effect was found over the studied ranges, and
  # Supplementary Figure S5 shows the corresponding correlation plots. No point
  # estimates are reported for them, so they cannot be implemented; they are
  # documented here so the provenance of the covariate screen is preserved.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Arm-mean baseline body weight.",
      units = "kg",
      type = "continuous",
      notes = paste0(
        "Screened as a candidate covariate on Emax and not retained (Guo 2025 ",
        "Results 3.1 / Discussion: 'no significant effects of baseline weight, ",
        "baseline BMI, or proportion of males ... were found'). Arm-mean ",
        "baseline weight ranged 72.2-121 kg (median 95.8 kg) across the 55 ",
        "trials; the paper states the null result holds over that range only. ",
        "No coefficient is reported, so the effect cannot be implemented."
      )
    ),
    BMI = list(
      description = "Arm-mean baseline body mass index.",
      units = "kg/m^2",
      type = "continuous",
      notes = paste0(
        "Screened as a candidate covariate on Emax and not retained (Guo 2025 ",
        "Results 3.1 / Discussion). Arm-mean baseline BMI ranged 24.1-45.1 ",
        "kg/m^2 (median 33.9) across the 55 trials. No coefficient is reported."
      )
    ),
    SEXF = list(
      description = "Arm-level proportion of participants who are female, expressed here on the canonical female-referenced scale.",
      units = "fraction",
      type = "continuous",
      notes = paste0(
        "Guo 2025 screened the arm-level MALE ratio and did not retain it (Guo ",
        "2025 Results 3.1 / Discussion). Recorded here on the canonical ",
        "female-referenced SEXF scale, so the source values invert as ",
        "SEXF = 1 - male_ratio. The source male ratio ranged 19.1-100.0 ",
        "percent with a median of 51.6 percent, i.e. a female fraction of ",
        "0-80.9 percent with a median of 48.4 percent. This is an ARM-LEVEL ",
        "aggregate proportion, not the individual-level binary SEXF of a popPK ",
        "model. No coefficient is reported, so the effect cannot be implemented."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 16269L,
    n_studies = 55L,
    n_drugs = 12L,
    age_range = paste0(
      "arm-mean age 29.5-64.7 years across the 55 trials; the median of the ",
      "trial means is 53.6 years, which is the centring constant used by the ",
      "age covariate on Emax"
    ),
    weight_range = "arm-mean baseline weight 72.2-121 kg; median 95.8 kg",
    bmi_range = "arm-mean baseline BMI 24.1-45.1 kg/m^2; median 33.9 kg/m^2",
    sex_female_pct = NA_real_,
    sex_male_pct = paste0(
      "arm-level male proportion 19.1-100.0 percent, median 51.6 percent ",
      "(equivalently a female proportion of 0-80.9 percent, median 48.4 ",
      "percent). Reported by the source as a male ratio; recorded as such ",
      "here because that is the direction the paper screened."
    ),
    race_ethnicity = paste0(
      "Not reported at the arm level and not screened as a covariate. The ",
      "inclusion criteria are BMI-threshold-adjusted for Japanese (BMI >= 23 ",
      "kg/m^2) and Chinese (BMI >= 24 kg/m^2) participants, so East Asian ",
      "cohorts are represented."
    ),
    disease_state = paste0(
      "Adults (>= 18 years) with overweight or obesity, defined as BMI >= 25 ",
      "kg/m^2 (>= 23 kg/m^2 for Japanese and >= 24 kg/m^2 for Chinese ",
      "participants). Trials of GLP-1RAs in participants with and without type ",
      "2 diabetes are pooled; the analysis does not stratify on diabetes status."
    ),
    dose_range = paste0(
      "By drug, over the 52-week simulation tiers of Guo 2025 Results 3.3: ",
      "liraglutide 1.2-3 mg (studied range per Discussion; 1.8-3 mg simulated), ",
      "injectable semaglutide 0.05-2.4 mg, oral semaglutide 1-40 mg, exenatide ",
      "0.01-2 mg, danuglipron up to 200 mg, orforglipron up to 45 mg, ",
      "tirzepatide 5-15 mg, cotadutide up to 0.6 mg, mazdutide 3-10 mg, ",
      "BI 456906 1.8-4.8 mg, JNJ-64565111 up to 10 mg, retatrutide up to 12 mg."
    ),
    treatment_duration = "6-104 weeks across trials; median 26 weeks",
    regions = paste0(
      "International. PubMed and Embase searched from inception to 2024-01-20; ",
      "English-language clinical trials only, which the paper notes as a ",
      "potential source of publication bias."
    ),
    notes = paste0(
      "MBMA at the study-arm level: each modelled data point is the arm-mean ",
      "placebo-adjusted change from baseline in body weight (ddWeight) in one ",
      "trial arm at one follow-up time, digitised with Engauge Digitizer 11.3 ",
      "when only published figures were available (extraction error held below ",
      "2 percent). Intention-to-treat results were preferred over per-protocol. ",
      "The etas are BETWEEN-TRIAL-GROUP, not between-subject: the model is ",
      "intended for simulating arm-mean weight-reduction time courses and is ",
      "NOT suitable for individual-subject simulation. The residual is weighted ",
      "by the inverse square root of the arm sample size (Supplementary Methods ",
      "1 Equation 1), so an arm of N participants has residual SD addSd / ",
      "sqrt(N); addSd in ini() is the unit-weight value and the N weighting is ",
      "applied downstream (same convention as Mercier_2014, Yao_2023 and ",
      "Asiimwe_2025 in this package). 22 of the 55 trials (40 percent) were ",
      "rated moderate risk of bias and 33 (60 percent) low risk, by Cochrane ",
      "RoB 2.0. Parameter robustness was confirmed by sampling importance ",
      "resampling (1000 resamples, Supplementary Table S6 SIR columns) and by a ",
      "prediction-corrected VPC (Figure 2). The dropout and adverse-event ",
      "analyses in Guo 2025 Table 1 are pairwise meta-analyses of endpoint ",
      "relative risks, not a time-course model, and are NOT part of this model ",
      "file -- see the vignette Errata."
    )
  )

  ini({
    # ==================================================================
    # Drug-specific maximum weight-reduction effect, Emax (kg).
    #
    # SIGN: negative values are weight LOSS, matching the sign convention
    # of Supplementary Table S6 and of Equations 1-6, where every printed
    # numerator is negative. Emax is the asymptote of the placebo-ADJUSTED
    # change from baseline, so Emax = -4.25 kg means liraglutide arms
    # eventually lose 4.25 kg more than their placebo arms.
    #
    # SCALE: kept on the natural (linear, signed) scale rather than log,
    # because the values are negative and cannot be log-transformed. This
    # follows the Vargo_2014_statins_ezetimibe_mbma precedent, where the
    # signed Emax intercept is also carried linearly.
    #
    # Six of the twelve values below are the AMPLITUDE of an Emax-in-dose
    # term (the effect approached at an infinite dose) and six are the flat
    # per-arm effect for a drug with no estimable dose-response. Which is
    # which is set in model(), not here; the tabulated Emax value has the
    # same meaning in Table S6 either way.
    # ==================================================================

    # --- Drugs WITH an estimable dose-response (Guo 2025 Equations 1-6) ---
    emax_cotadutide <- -10.5
    label("Cotadutide maximum placebo-adjusted weight change at infinite dose (kg; negative = loss)")  # Guo 2025 Supplementary Table S6 Emax_Cotadutide = -10.5 (RSE 14.8%; SIR median -10.6, 95% CI -13.3 to -7.86); also the numerator of Equation 1

    emax_danuglipron <- -9.29
    label("Danuglipron maximum placebo-adjusted weight change at infinite dose (kg; negative = loss)")  # Guo 2025 Supplementary Table S6 Emax_Danuglipron = -9.29 (RSE 24.7%; SIR median -9.37, 95% CI -13.9 to -5.30); also the numerator of Equation 2

    emax_jnj64565111 <- -18.6
    label("JNJ-64565111 maximum placebo-adjusted weight change at infinite dose (kg; negative = loss)")  # Guo 2025 Supplementary Table S6 Emax_JNJ-64565111 = -18.6 (RSE 6.10%; SIR median -18.7, 95% CI -20.7 to -16.3); also the numerator of Equation 3

    emax_retatrutide <- -22.6
    label("Retatrutide maximum placebo-adjusted weight change at infinite dose (kg; negative = loss). Largest Emax of the 12 drugs; the only GLP-1/GIP/GCG tri-agonist.")  # Guo 2025 Supplementary Table S6 Emax_Retatrutide = -22.6 (RSE 19.3%; SIR median -22.6, 95% CI -30.9 to -15.2); also the numerator of Equation 4; the 22.6 kg upper end of the Abstract's "4.25 kg to 22.6 kg" range

    emax_orforglipron <- -14.7
    label("Orforglipron maximum placebo-adjusted weight change at infinite dose (kg; negative = loss)")  # Guo 2025 Supplementary Table S6 Emax_Orforglipron = -14.7 (RSE 7.40%; SIR median -14.6, 95% CI -16.7 to -12.4); also the numerator of Equation 5

    emax_semaglutide_inj <- -11.7
    label("Injectable (SC) semaglutide maximum placebo-adjusted weight change at infinite dose (kg; negative = loss)")  # Guo 2025 Supplementary Table S6 Emax_Semaglutide(INJ) = -11.7 (RSE 12.0%; SIR median -11.7, 95% CI -14.3 to -8.81); also the numerator of Equation 6

    # --- Drugs with NO estimable dose-response (flat per-arm Emax) ---
    # Guo 2025 Discussion: "no dose-response relationship was observed
    # within the studied dosage range for the other six GLP-1RA drugs ...
    # indicating that the dosages of these drugs have reached their efficacy
    # plateaus." For these six the tabulated Emax IS the arm effect at any
    # studied dose; the DOSE_* column acts only as a presence indicator.
    emax_bi456906 <- -13.5
    label("BI 456906 (survodutide) flat placebo-adjusted maximum weight change (kg; negative = loss); no dose-response estimable over 1.8-4.8 mg")  # Guo 2025 Supplementary Table S6 Emax_BI 456906 = -13.5 (RSE 14.8%; SIR median -13.8, 95% CI -17.9 to -9.66)

    emax_exenatide <- -6.05
    label("Exenatide flat placebo-adjusted maximum weight change (kg; negative = loss); no dose-response estimable over 0.01-2 mg")  # Guo 2025 Supplementary Table S6 Emax_Exenatide = -6.05 (RSE 8.60%; SIR median -6.10, 95% CI -7.09 to -5.12)

    emax_liraglutide <- -4.25
    label("Liraglutide flat placebo-adjusted maximum weight change (kg; negative = loss); no dose-response estimable over 1.2-3 mg. Smallest Emax of the 12 drugs.")  # Guo 2025 Supplementary Table S6 Emax_Liraglutide = -4.25 (RSE 9.30%; SIR median -4.26, 95% CI -5.08 to -3.46); the 4.25 kg lower end of the Abstract's "4.25 kg to 22.6 kg" range

    emax_semaglutide_po <- -5.36
    label("Oral semaglutide flat placebo-adjusted maximum weight change (kg; negative = loss); no dose-response estimable over 1-40 mg")  # Guo 2025 Supplementary Table S6 Emax_Semaglutide(P.O) = -5.36 (RSE 22.8%; SIR median -5.41, 95% CI -7.88 to -3.16)

    emax_tirzepatide <- -12.9
    label("Tirzepatide flat placebo-adjusted maximum weight change (kg; negative = loss); no dose-response estimable over 5-15 mg")  # Guo 2025 Supplementary Table S6 Emax_Tirzepatide = -12.9 (RSE 9.10%; SIR median -12.9, 95% CI -15.4 to -10.6)

    emax_mazdutide <- -7.75
    label("Mazdutide flat placebo-adjusted maximum weight change (kg; negative = loss); no dose-response estimable over 3-10 mg")  # Guo 2025 Supplementary Table S6 Emax_Mazdutide = -7.75 (RSE 17.9%; SIR median -7.73, 95% CI -10.8 to -5.11)

    # ==================================================================
    # ED50 -- dose producing half of Emax, for the six drugs with an
    # estimable dose-response. Every one is reported as "Fixed" in the
    # Supplementary Table S6 "RSE(%)" column and has no SIR interval, so
    # all six are wrapped in fixed(). Carried on the log scale so the
    # back-transform in model() is exp(); the log() is INSIDE fixed()
    # per the package convention.
    # ==================================================================
    led50_cotadutide <- fixed(log(0.219))
    label("Log cotadutide ED50 (mg). Back-transformed 0.219 mg.")  # Guo 2025 Supplementary Table S6 "thetaDose on Emax_Cotadutide" = 0.219, Fixed; Results 3.2 lists it as the effective dose; denominator of Equation 1

    led50_danuglipron <- fixed(log(80))
    label("Log danuglipron ED50 (mg). Back-transformed 80 mg.")  # Guo 2025 Supplementary Table S6 "thetaDose on Emax_Danuglipron" = 80, Fixed; denominator of Equation 2

    led50_jnj64565111 <- fixed(log(6.73))
    label("Log JNJ-64565111 ED50 (mg). Back-transformed 6.73 mg.")  # Guo 2025 Supplementary Table S6 "thetaDose on Emax_JNJ-64565111" = 6.73, Fixed; denominator of Equation 3

    led50_retatrutide <- fixed(log(4))
    label("Log retatrutide ED50 (mg). Back-transformed 4 mg.")  # Guo 2025 Supplementary Table S6 "thetaDose on Emax_Retatrutide" = 4, Fixed; denominator of Equation 4

    led50_orforglipron <- fixed(log(14.6))
    label("Log orforglipron ED50 (mg). Back-transformed 14.6 mg.")  # Guo 2025 Supplementary Table S6 "thetaDose on Emax_Orforglipron" = 14.6, Fixed; denominator of Equation 5

    led50_semaglutide_inj <- fixed(log(0.384))
    label("Log injectable semaglutide ED50 (mg). Back-transformed 0.384 mg.")  # Guo 2025 Supplementary Table S6 "thetaDose on Emax_Semaglutide(INJ)" = 0.384, Fixed; denominator of Equation 6

    # ==================================================================
    # Onset of effect.
    # ==================================================================
    lkel <- log(0.0578)
    label("Log first-order onset rate constant k (1/week) for the exponential approach of the effect to Emax. Back-transformed 0.0578 /week, giving a time to half of Emax of ET50 = 0.693/k = 12.0 weeks. A SINGLE shared k is estimated in the final model: Guo 2025 Results 3.2 states that 'owing to the limited number of time points available for some drugs, it was not possible to estimate the k values individually for each drug', and drug-specific k values were instead recovered post hoc by Bayesian feedback combined with single-arm meta-analysis (Supplementary Figure S3, ET50 6.4 weeks for orforglipron to 19.5 weeks for tirzepatide). Those per-drug values are NOT part of the final model and are not implemented here; see the vignette Errata.")  # Guo 2025 Supplementary Table S6 k = 0.0578 /week (RSE 8.60%; SIR median 0.0572, 95% CI 0.0480 to 0.0666). Name follows the Mercier_2014_tramadol_tapentadol_mbma precedent, which uses lkel for the same approach-to-plateau rate constant.

    # ==================================================================
    # Covariate effect: arm-mean age on Emax (Guo 2025 Equation 7).
    # ==================================================================
    e_age_emax <- -0.0304
    label("Exponential coefficient of arm-mean age on Emax, per year, centred on the 53.6-year median trial age: Emax = Emax_typical * exp(e_age_emax * (AGE - 53.6)). Negative on a negative Emax means the magnitude of weight loss SHRINKS with increasing age. The paper's own arithmetic check: exp(-0.0304 * 10) = 0.738, the '26.2% decrease in the Emax value' quoted for a 40-to-50-year increase in mean age.")  # Guo 2025 Supplementary Table S6 "thetaAge on Emax" = -0.0304 (RSE 23.0%; SIR median -0.0306, 95% CI -0.0428 to -0.0160); Equation 7

    # ==================================================================
    # BETWEEN-STUDY (inter-trial-group) variability -- NOT between-subject.
    #
    # Supplementary Methods 1 Equations 2-3 define both as exponential on
    # the parameter: Emax_i = Emax_typical * exp(eta_Emax) and
    # k_i = k_typical * exp(eta_k), with "eta_Emax and eta_k ... both
    # follow normal distributions centered at 0 with variances of omega1^2
    # and omega2^2". The Supplementary Table S6 rows labelled "etaEmax" and
    # "etak" are therefore VARIANCES, and are used here as nlmixr2 eta
    # variances without transformation.
    #
    # The variance (rather than SD) reading is confirmed by the reported
    # RSEs. NONMEM's RSE on a variance is approximately sqrt(2/n_eff): the
    # 11.1% RSE on both omegas implies n_eff = 2/0.111^2 = 162 study arms,
    # which matches an analysis of 55 trials contributing multiple dose
    # arms each. Read as SDs the same RSEs would imply only ~41 arms,
    # fewer than the 55 trials.
    # ==================================================================
    eta_study_emax ~ 0.330
    label("Between-trial-group variance of Emax (exponential; Emax_i = Emax_typical * exp(eta))")  # Guo 2025 Supplementary Table S6 etaEmax = 0.330 (RSE 11.1%; SIR median 0.339, 95% CI 0.271 to 0.404); variance per Supplementary Methods 1 Eq 2

    eta_study_lkel ~ 0.627
    label("Between-trial-group variance of the onset rate constant k (exponential; k_i = k_typical * exp(eta))")  # Guo 2025 Supplementary Table S6 etak = 0.627 (RSE 11.1%; SIR median 0.626, 95% CI 0.479 to 0.744); variance per Supplementary Methods 1 Eq 3

    # ==================================================================
    # Residual error.
    #
    # Supplementary Methods 1 Equation 1 is
    #   E_ij = Emax_i * (1 - exp(-k_i * time_j)) + eps_ij / sqrt(N_ij)
    # with eps "normally distributed with a mean of 0 and a variance of
    # sigma^2 ... weighted by the inverse of the square root of the sample
    # size". Supplementary Table S6 reports eps = 0.391, which is that
    # sigma^2, so the SD at UNIT study weight is sqrt(0.391) = 0.625.
    # The expression is left as sqrt(0.391) so the tabulated variance stays
    # visible at the assignment site.
    #
    # nlmixr2's add() takes a constant SD, so the 1/sqrt(N) arm-size
    # weighting cannot be expressed inside the error model. It is applied
    # DOWNSTREAM: to simulate an arm of N participants, scale the residual
    # SD by 1/sqrt(N). This is the same convention already used by
    # Mercier_2014_tramadol_tapentadol_mbma, the three Yao_2023 SGLT2
    # files, and the two Asiimwe_2025 ADC files in this package.
    # ==================================================================
    addSd <- sqrt(0.391)
    label("Additive residual SD on arm-mean placebo-adjusted weight change (kg) at UNIT study weight; the operative SD for an arm of N participants is addSd / sqrt(N)")  # Guo 2025 Supplementary Table S6 eps = 0.391 (RSE 6.30%; SIR median 0.389, 95% CI 0.342 to 0.440), a variance per Supplementary Methods 1; sqrt(0.391) = 0.625
  })

  model({
    # ----- Centring constant -----
    # Median of the 55 trial mean ages (Guo 2025 Results 3.1: "a median age
    # of 53.6 years"), the COV_median of Supplementary Methods 1 Equation 7.
    ref_age_years <- 53.6

    # ----- ED50 back-transforms (mg) -----
    ed50_cotadutide       <- exp(led50_cotadutide)
    ed50_danuglipron      <- exp(led50_danuglipron)
    ed50_jnj64565111      <- exp(led50_jnj64565111)
    ed50_retatrutide      <- exp(led50_retatrutide)
    ed50_orforglipron     <- exp(led50_orforglipron)
    ed50_semaglutide_inj  <- exp(led50_semaglutide_inj)

    # ----- Per-drug Emax contribution -----
    #
    # Standard MBMA usage: exactly ONE of the twelve DOSE_* columns is
    # non-zero in any study arm (Guo 2025 analysed placebo-controlled
    # monotherapy trials; no arm combined two GLP-1RAs), so the sum below
    # collapses to the single active drug's term. If a downstream user
    # codes an arm with two agonists at once the model returns the additive
    # sum of their individual Emax contributions, which is outside the
    # paper's calibration.
    #
    # Six drugs use the Emax-in-dose form of Supplementary Methods 1
    # Equation 4, Emax * Dose / (ED50 + Dose), which printed as Equations
    # 1-6 in the main text. Each term is automatically zero when its dose
    # column is zero, so no separate indicator is needed.
    emax_dose_cotadutide      <- emax_cotadutide      * DOSE_COTADUTIDE_MG      / (ed50_cotadutide      + DOSE_COTADUTIDE_MG)
    emax_dose_danuglipron     <- emax_danuglipron     * DOSE_DANUGLIPRON_MG     / (ed50_danuglipron     + DOSE_DANUGLIPRON_MG)
    emax_dose_jnj64565111     <- emax_jnj64565111     * DOSE_JNJ64565111_MG     / (ed50_jnj64565111     + DOSE_JNJ64565111_MG)
    emax_dose_retatrutide     <- emax_retatrutide     * DOSE_RETATRUTIDE_MG     / (ed50_retatrutide     + DOSE_RETATRUTIDE_MG)
    emax_dose_orforglipron    <- emax_orforglipron    * DOSE_ORFORGLIPRON_MG    / (ed50_orforglipron    + DOSE_ORFORGLIPRON_MG)
    emax_dose_semaglutide_inj <- emax_semaglutide_inj * DOSE_SEMAGLUTIDE_INJ_MG / (ed50_semaglutide_inj + DOSE_SEMAGLUTIDE_INJ_MG)

    # The remaining six drugs had no estimable dose-response, so their
    # tabulated Emax applies flat at any dose actually studied. The
    # comparison evaluates to 1 when the arm received the drug and 0
    # otherwise, giving the same "zero unless present" behaviour as the
    # Emax-in-dose terms above.
    emax_flat_bi456906       <- emax_bi456906       * (DOSE_BI456906_MG > 0)
    emax_flat_exenatide      <- emax_exenatide      * (DOSE_EXENATIDE_MG > 0)
    emax_flat_liraglutide    <- emax_liraglutide    * (DOSE_LIRAGLUTIDE_MG > 0)
    emax_flat_semaglutide_po <- emax_semaglutide_po * (DOSE_SEMAGLUTIDE_PO_MG > 0)
    emax_flat_tirzepatide    <- emax_tirzepatide    * (DOSE_TIRZEPATIDE_MG > 0)
    emax_flat_mazdutide      <- emax_mazdutide      * (DOSE_MAZDUTIDE_MG > 0)

    emax_drug <-
      emax_dose_cotadutide + emax_dose_danuglipron + emax_dose_jnj64565111 +
      emax_dose_retatrutide + emax_dose_orforglipron + emax_dose_semaglutide_inj +
      emax_flat_bi456906 + emax_flat_exenatide + emax_flat_liraglutide +
      emax_flat_semaglutide_po + emax_flat_tirzepatide + emax_flat_mazdutide

    # ----- Age effect on Emax and between-trial-group variability -----
    # Guo 2025 Equation 7 (age) and Supplementary Methods 1 Equations 2-3
    # (etas). Both are multiplicative on Emax, so the order of the two
    # factors does not matter.
    emax_i <- emax_drug * exp(e_age_emax * (AGE - ref_age_years)) * exp(eta_study_emax)

    # ----- Onset rate constant for this trial group -----
    k_i <- exp(lkel + eta_study_lkel)

    # ----- Arm-mean placebo-adjusted weight change (kg) -----
    # Supplementary Methods 1 Equation 1. The effect is 0 at t = 0 and
    # approaches emax_i as t grows; ET50 = 0.693 / k_i.
    #
    # Cc is the package's canonical single-output observation name. It is
    # overloaded here and is NOT a drug concentration: it is the arm-mean
    # difference between the drug arm's and the placebo arm's change from
    # baseline in body weight, in kg, negative for weight loss.
    Cc <- emax_i * (1 - exp(-k_i * t))

    # Residual SD is the UNIT-WEIGHT value; divide by sqrt(arm N) downstream.
    Cc ~ add(addSd)
  })
}
