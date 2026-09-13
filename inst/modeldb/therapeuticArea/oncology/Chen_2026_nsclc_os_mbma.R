Chen_2026_nsclc_os_mbma <- function() {
  description <- paste(
    "MBMA.",
    "Mixed-effects non-parametric conditional-probability model-based",
    "meta-analysis (MBMA) of OVERALL SURVIVAL in treatment-naive,",
    "locally advanced or metastatic non-small cell lung cancer (NSCLC)",
    "patients who are ineligible for platinum-doublet chemotherapy",
    "(ECOG performance status >= 2, or age >= 70 years, or",
    "investigator-determined platinum ineligibility). Chen 2026 digitised",
    "Kaplan-Meier overall-survival curves from 26 published trials",
    "(41 arms, 3637 participants) of single-agent docetaxel, gemcitabine,",
    "paclitaxel, pemetrexed or vinorelbine, and compared them against the",
    "control arm of the phase III IPSOS trial (NCT03191786, Lee 2023),",
    "in which investigators chose gemcitabine or vinorelbine. This file",
    "encodes the paper's FINAL model, Model 010, whose complete parameter",
    "set is Chen 2026 Table S8.",
    "The survival model is semi-parametric proportional hazards with a",
    "NON-PARAMETRIC reference survival curve: the logit of the reference",
    "survival probability follows an estimated random walk over 30 knots",
    "(months 1-24, then 27, 30, 33, 36, 42 and 48), and every study arm is",
    "a proportional rescaling of it, S_arm(t) = S_ref(t)^exp(log(HR)).",
    "The reference is a historical single-agent-chemotherapy control arm",
    "composed entirely of ECOG performance status 0 patients. log(HR) is a",
    "linear combination of the arm's ECOG-PS-1 and ECOG-PS-2/3 percentages,",
    "an IPSOS-control-arm indicator, and a between-trial random effect",
    "(SD 0.257). ECOG performance status was the ONLY covariate that",
    "reached statistical significance; age, sex, disease stage, squamous",
    "and adenocarcinoma histology, Asian region and publication year were",
    "all screened and rejected (see covariatesDataExcluded).",
    "All parameter values are wrapped in fixed() because the model is a",
    "downstream user of the published fit, not a re-estimation of it.",
    "Simulation scope: per-ARM overall-survival curves over the 1-48 month",
    "window the digitised data support. The random effect is BETWEEN-TRIAL,",
    "not between-subject, so this model predicts trial-arm mean survival",
    "curves and is NOT suitable for individual-subject time-to-event",
    "simulation. No residual error is estimated: Chen 2026 fits a",
    "user-supplied binomial -2LL on the number of deaths per interval",
    "rather than a residual-error model, so there is no $SIGMA to encode",
    "and no observation endpoint is declared.",
    "Companion paper in the same issue: Franzese_2026_pdl1_nsclc_mbma."
  )

  reference <- paste(
    "Chen J, Wada R, Zhang N, Graupner V, Morris S, Hu Y, Zhang W,",
    "Kassir N, Wu B, Chan P.",
    "Model-Based Meta-Analysis of Overall Survival in Vulnerable",
    "Platinum-Ineligible NSCLC Populations.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15:e70197.",
    "doi:10.1002/psp4.70197.",
    "Parameter values are Table S8 ('Final Model Parameters (mod010)') of",
    "the Supporting Information; the model structure is the NONMEM $PRED",
    "block distributed as Supporting Information file",
    "psp470197-sup-0002-Supinfo2.ctl.",
    sep = " "
  )

  vignette <- "Chen_2026_nsclc_os_mbma"

  # `eta_study` is an MBMA BETWEEN-TRIAL random effect on log(HR). It does not
  # pair 1:1 with a fixed-effect parameter of the same name -- it is an
  # additive shift on the whole linear predictor, exactly as ETA(1)*OM1 enters
  # LOGRRI in the source control stream. The pbpk-qsp-mbma convention
  # prescribes eta_study_<name> for MBMA study-level etas and requires them to
  # be distinguished from between-subject variability. Same device as
  # Franzese_2026_pdl1_nsclc_mbma.
  paper_specific_etas <- c("eta_study")

  # Algebraic study-arm-level MBMA: there is no PK compartment, no dose event
  # and no drug concentration anywhere in the model -- the only output is a
  # dimensionless survival probability in [0, 1]. The `dosing` /
  # `concentration` strings below are placeholders chosen so that
  # checkModelConventions() sees a dimensionally consistent pair; a more
  # descriptive label trips the dimensional-compatibility check. Same device as
  # Volkova_2023_lipidLowering_mace_mbma and Yoshioka_2018_FXa_inhibitors_mbma.
  units <- list(
    time          = "month",
    dosing        = "probability",
    concentration = "probability/probability"
  )

  covariateData <- list(
    PS_ECOG_1_PCT = list(
      description        = "Study-arm-level percentage (0-100) of the enrolled cohort with an ECOG performance status of 1 at baseline.",
      units              = "%",
      type               = "continuous",
      reference_category = "0% (no ECOG-1 patients). The model's reference arm is PS_ECOG_1_PCT = PS_ECOG_2_PCT = PS_ECOG_3_PCT = 0, i.e. an all-ECOG-0 arm.",
      source_name        = "PS1 (Chen 2026 NONMEM dataset psp470197-sup-0001-Supinfo1.csv $INPUT); 'PS.1' (Chen 2026 Tables S5 and S6); 'ECOG PS 1 (%)' (Chen 2026 Table 2)",
      notes              = paste(
        "MBMA study-arm-level covariate, scaled in PERCENT (0-100), not a fraction -- the model divides by 100 internally, matching ECOGEFF = COVECOG1*PS1/100 in the source control stream.",
        "Enters log(HR) with coefficient e_ecog1 = 0.004, i.e. hazard ratio exp(0.004) = 1.004 for an arm that is 100% ECOG 1 relative to an all-ECOG-0 arm (Chen 2026 Table 3, Model 010, ECOG1HR column).",
        "The estimate is essentially null and very imprecisely determined (RSE 11922%, 95% CI 0.412-2.443), which the paper reports as-is; it is retained because Model 010 estimates the ECOG effects from the MBMA data rather than fixing them.",
        "Across the 41 analysis arms the median is 42.9% (range 0.0-69.6; Chen 2026 Table 2).",
        "PS_ECOG_0_PCT, PS_ECOG_1_PCT, PS_ECOG_2_PCT and PS_ECOG_3_PCT sum to 100 within an arm, so supplying all four is redundant; the model consumes only PS 1, 2 and 3 and treats ECOG 0 as the reference level.",
        "Distinct from the per-subject binaries ECOG_GE1 / ECOG_GE2, which are individual-level indicators. Family member of PS_ECOG_0_PCT (founding example Franzese 2026).",
        "Chen 2026 Methods 2.2 notes that 30 of the ECOG PS percentages in the database were IMPUTED -- first by mapping Karnofsky Performance Scores where available, otherwise by splitting composite ECOG categories with a logistic regression calibrated on studies with similar entry criteria.",
        sep = " "
      )
    ),
    PS_ECOG_2_PCT = list(
      description        = "Study-arm-level percentage (0-100) of the enrolled cohort with an ECOG performance status of 2 at baseline.",
      units              = "%",
      type               = "continuous",
      reference_category = "0% (no ECOG-2 patients); ECOG 0 is the model's reference performance-status level.",
      source_name        = "PS2 (Chen 2026 NONMEM dataset psp470197-sup-0001-Supinfo1.csv $INPUT); 'PS.2' (Chen 2026 Tables S5 and S6)",
      notes              = paste(
        "MBMA study-arm-level covariate, scaled in PERCENT (0-100).",
        "PS_ECOG_2_PCT and PS_ECOG_3_PCT SHARE a single coefficient in this model: the source control stream computes ECOGEFF using COVECOG2*(PS2+PS3)/100, so ECOG 2 and ECOG 3 are pooled into one composite effect (e_ecog23 = 0.769, hazard ratio 2.158; Chen 2026 Table 3 Model 010 ECOG2,3HR column).",
        "They are kept as two separate columns rather than a single pooled PS_ECOG_23_PCT because the source dataset carries them separately and because Chen 2026 Table 2 reports them pooled while Tables S5/S6 report them apart.",
        "Pooled across ECOG 2 and 3, the median across the 41 analysis arms is 28.8% (range 0.0-100.0; Chen 2026 Table 2).",
        "The pooled ECOG-2/3 percentage is the single most load-bearing covariate in the paper: the IPSOS control arm is 87.4% ECOG 2 or 3 (76.8% + 10.6%), against a database median of 28.8%, which is why adjusting for performance status moves the IPSOS-versus-historical hazard ratio from 0.812 (Model 007a, unadjusted) to 0.531 (Model 010, adjusted).",
        "See PS_ECOG_1_PCT for the imputation caveat, which applies equally here.",
        sep = " "
      )
    ),
    PS_ECOG_3_PCT = list(
      description        = "Study-arm-level percentage (0-100) of the enrolled cohort with an ECOG performance status of 3 at baseline.",
      units              = "%",
      type               = "continuous",
      reference_category = "0% (no ECOG-3 patients); ECOG 0 is the model's reference performance-status level.",
      source_name        = "PS3 (Chen 2026 NONMEM dataset psp470197-sup-0001-Supinfo1.csv $INPUT); 'PS.3' (Chen 2026 Tables S5 and S6)",
      notes              = paste(
        "MBMA study-arm-level covariate, scaled in PERCENT (0-100).",
        "Shares the coefficient e_ecog23 with PS_ECOG_2_PCT -- see the PS_ECOG_2_PCT notes for why the two are pooled and why they are nevertheless carried as separate columns.",
        "ECOG 3 is rare in this database: it is non-zero in only four of the 41 arms (Chen 2026 Table S6), the largest being the IPSOS trial arms (10.6% in the control arm).",
        sep = " "
      )
    ),
    TRT_IPSOS_CONTROL = list(
      description        = "Indicator that the study arm is the control arm of the phase III IPSOS trial (NCT03191786, reported as Lee 2023), in which investigators chose single-agent gemcitabine or vinorelbine.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (a historical single-agent-chemotherapy control arm from the published literature -- the model's reference treatment)",
      source_name        = "ID == 15 with the atezolizumab arm excluded (Chen 2026 NONMEM dataset psp470197-sup-0001-Supinfo1.csv; the source control stream selects the trial arm with IF (ID.EQ.15))",
      notes              = paste(
        "MBMA trial-arm indicator, NOT a patient characteristic. 1 = the IPSOS control arm (Lee 2023, ID 15, ARM 2, n = 151, drug code 3 = 'gemcitabine or vinorelbine'); 0 = any of the 39 historical control arms.",
        "Enters log(HR) with coefficient e_ipsosctl = -0.633, i.e. hazard ratio exp(-0.633) = 0.531 (95% CI 0.424-0.664), reproducing the IPSOSHR column of Chen 2026 Table 3 for Model 010.",
        "Model 010 is fitted with the IPSOS ATEZOLIZUMAB arm EXCLUDED from the dataset (Chen 2026 Table 5, Data column), so within this model the IPSOS trial identifier and the IPSOS control arm are the same thing. Do NOT set this flag on an atezolizumab arm: the atezolizumab effect belongs to Model 041, whose reference survival curve Chen 2026 does not publish (see vignette Errata).",
        "The paper's headline hazard ratios of 0.543 (IPSOS control) and 0.418 (IPSOS atezolizumab) versus historical trials come from Model 041, not Model 010; Model 010's 0.531 is the corresponding quantity in the final model.",
        "Family member of the TRT_<arm> treatment-arm-indicator canonicals (TRT_EPHEDRINE, TRT_PBT, TRT_PCSK9I, ...).",
        sep = " "
      )
    )
  )

  # Covariates Chen 2026 extracted per study arm and formally SCREENED as
  # candidate effects on log(HR), but which did not reduce the objective
  # function by the pre-specified 6.635 points (alpha = 0.01) and are therefore
  # absent from the final model. Documentation only -- they are deliberately
  # not referenced in model(). Chen 2026 Section 3.2.2 and Table S7.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Median age of the patients in the study arm.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a candidate covariate on log(HR) (tested as an age category) and not retained. Median across arms 74.0 years, range 61.0-79.0 (Chen 2026 Table 2). Source column 'Age' (Table S1)."
    ),
    SEXM_PCT = list(
      description = "Percentage (0-100) of the study arm who are male.",
      units       = "%",
      type        = "continuous",
      notes       = "Screened and not retained; Model 012 gave a covariate HR of 1.070 (95% CI 0.387-2.954) with dOBJ = +0.0 versus Model 008 (Chen 2026 Table S7). Median across arms 74.7%, range 38.9-94.2 (Table 2). Source column 'MaleP' (Table S1). Recorded as a male percentage because that is the direction the source reports; the individual-level canonical SEXF is coded 1 = female."
    ),
    DIS_STAGE4_PCT = list(
      description = "Percentage (0-100) of the study arm with Stage IV (metastatic) disease.",
      units       = "%",
      type        = "continuous",
      notes       = "Screened and not retained; Model 011 gave a covariate HR of 0.624 (95% CI 0.210-1.852) (Chen 2026 Table S7). Median across arms 73.1%, range 41.3-88.9 (Table 2). Source column 'Stage4' (Tables S5 and S6)."
    ),
    TUMTP_SQUAM_PCT = list(
      description = "Percentage (0-100) of the study arm with squamous-cell histology.",
      units       = "%",
      type        = "continuous",
      notes       = "Screened and not retained; Model 013 gave a covariate HR of 1.125 (95% CI 0.454-2.788) (Chen 2026 Table S7). Median across arms 39.3%, range 0.0-54.4 (Table 2). Three squamous percentages were imputed (Chen 2026 Section 3.1.2). Source column 'SquamousP' (Table S1)."
    ),
    TUMTP_ADENO_PCT = list(
      description = "Percentage (0-100) of the study arm with adenocarcinoma histology.",
      units       = "%",
      type        = "continuous",
      notes       = "Screened and not retained; Model 014 gave a covariate HR of 0.625 (95% CI 0.250-1.560) (Chen 2026 Table S7). Median across arms 39.3%, range 14.0-87.0 (Table 2). Five adenocarcinoma percentages were imputed (Chen 2026 Section 3.1.2). Source column 'AdenoP' (Table S1)."
    ),
    RACE_ASIAN_PCT = list(
      description = "Percentage (0-100) of the study arm enrolled in Asia.",
      units       = "%",
      type        = "continuous",
      notes       = "Screened and not retained; Model 015 gave a covariate HR of 0.817 (95% CI 0.658-1.015) (Chen 2026 Table S7). Bimodal across arms -- median 0.0% with range 0.0-100.0 (Table 2), because trials were either wholly Asian or wholly non-Asian. Tested because Asian patients may survive longer on immune-checkpoint-inhibitor therapy (Chen 2026 Section 2.2). Source column 'Asia' (Tables S5 and S6)."
    ),
    YEAR_PUBLICATION = list(
      description = "Calendar year in which the trial was published.",
      units       = "year",
      type        = "continuous",
      notes       = "Screened as a time-trend covariate and not retained: Model 049 reduced the objective function by only 3.3 points against the 6.635 required (Chen 2026 Table 5). The coefficient was negative (more recent publication associated with lower hazard), consistent with medical practice improving over time, and including it widened the IPSOS hazard ratio to 0.780 (95% CI 0.479-1.269), a CI that includes 1. Source column 'YEAR' (NONMEM dataset)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 3637L,
    n_studies      = 26L,
    n_arms         = 41L,
    age_range      = "Arm-level median age 74.0 years, range across arms 61.0-79.0 (Chen 2026 Table 2). Individual ages are not available: this is a summary-level meta-analysis.",
    weight_range   = "Not collected. Body weight is not a covariate in any Chen 2026 model and is not reported in the analysis database.",
    sex_female_pct = 25.3,
    disease_state  = paste(
      "Treatment-naive, locally advanced or metastatic non-small cell lung cancer in patients unsuitable for platinum-doublet chemotherapy.",
      "Eligibility for the literature database required ECOG performance status >= 2, or age >= 70 years, or platinum ineligibility as defined by the respective published study.",
      "Arm-level medians (range) across the 41 arms: ECOG PS 0 14.7% (0.0-57.4), ECOG PS 1 42.9% (0.0-69.6), ECOG PS 2 or 3 28.8% (0.0-100.0), Stage IV disease 73.1% (41.3-88.9), squamous histology 39.3% (0.0-54.4), adenocarcinoma histology 39.3% (14.0-87.0) (Chen 2026 Table 2).",
      sep = " "
    ),
    dose_range     = "Not modelled. Chen 2026 carries treatment identity, not dose: no dose, schedule or exposure metric enters any model.",
    regions        = "Multinational. Arm-level percentage enrolled in Asia has median 0.0% and range 0.0-100.0 (Chen 2026 Table 2) -- the constituent trials were each wholly Asian or wholly non-Asian.",
    treatments     = paste(
      "Single-agent chemotherapy by drug (Chen 2026 Table 1, studies / arms / n):",
      "docetaxel 7 / 11 / 711; gemcitabine 10 / 12 / 969; gemcitabine or vinorelbine 2 / 2 / 341;",
      "paclitaxel 3 / 3 / 311; pemetrexed 2 / 2 / 90; vinorelbine 9 / 10 / 913; atezolizumab 1 / 1 / 302.",
      "Studies sum to more than 26 because some trials contributed more than one drug.",
      sep = " "
    ),
    notes          = paste(
      "SUMMARY-LEVEL meta-analysis. Each 'subject' in nlmixr2 corresponds to one study-strata-arm, carrying that arm's mean covariate values and its digitised Kaplan-Meier overall-survival curve; the random effect eta_study is BETWEEN-TRIAL, not between-subject.",
      "sex_female_pct is DERIVED as 100 minus the arm-level MEDIAN male percentage of 74.7% (Chen 2026 Table 2); it is a median across arms rather than a participant-weighted pooled percentage, which the paper does not report.",
      "Data were obtained by digitising published Kaplan-Meier curves at monthly intervals with Engauge Digitizer 12.1 (Chen 2026 Section 2.2); censoring tick marks could not be recovered, so the analysis is of conditional death probabilities per interval rather than of individual event times.",
      "Model 010 as encoded here is fitted with the IPSOS ATEZOLIZUMAB arm (n = 302) excluded, leaving 3335 participants across 40 arms (Chen 2026 Table 5, Data column).",
      "Curve tails were trimmed before fitting wherever fewer than 10 patients remained at risk or the survival probability fell below 0.1 (the FLAG = 1 and FLAG = 2 exclusions in the distributed dataset), which is why the knot grid stops at 48 months and why the last interval is unidentifiable and fixed to zero.",
      "The literature search covered PubMed to 18 September 2023 and returned 275 references, of which 28 (26 unique studies) were retained; a funnel plot showed no material publication bias (Chen 2026 Figure S3).",
      sep = " "
    )
  )

  ini({
    # =========================================================================
    # NON-PARAMETRIC REFERENCE SURVIVAL CURVE (Chen 2026 Table S8, rows 1-29).
    #
    # The reference arm is a HISTORICAL single-agent-chemotherapy control arm
    # composed entirely of ECOG performance status 0 patients, i.e.
    # PS_ECOG_1_PCT = PS_ECOG_2_PCT = PS_ECOG_3_PCT = 0 and
    # TRT_IPSOS_CONTROL = 0, at the typical value of the between-trial random
    # effect.
    #
    # The curve is parameterised as a monotone random walk on the LOGIT of the
    # reference survival probability, exactly as in the source $PRED block:
    #
    #     X1 = THETA(1)                 -> logit S_ref(1 month)
    #     Xk = X(k-1) - THETA(k)        -> logit S_ref(k-th knot)
    #     S_ref = exp(X)/(1 + exp(X))
    #
    # Every decrement is constrained positive in the source control stream
    # ($THETA lower bound 0), which is what forces the curve to be monotone
    # non-increasing. The knots are months 1-24 and then 27, 30, 33, 36, 42
    # and 48, matching the TIME grid of the distributed analysis dataset.
    #
    # All values are wrapped in fixed(): this file is a downstream user of the
    # published fit, not a re-estimation of it. The standard error quoted in
    # each trailing comment is the Table S8 SE; the printed RSE column is
    # reproduced by SE / Value * 100 throughout, which is how Table S8 was
    # checked for transcription errors.
    # =========================================================================

    lgtsurv_01  <- fixed(2.802) ; label("Logit of reference overall-survival probability at 1 month (logit scale)")                            # Chen 2026 Table S8 'logit(Survival at 1-mo)' = 2.802 (SE 0.378, RSE 13.5%). inv-logit = 0.9428
    dlgtsurv_02 <- fixed(0.831) ; label("Decrement in logit reference survival from month 1 to month 2 (logit scale)")                          # Chen 2026 Table S8 'logit(Surv at 1-mo)-logit(Surv at 2-mo)' = 0.831 (SE 0.042, RSE 5.1%)
    dlgtsurv_03 <- fixed(0.520) ; label("Decrement in logit reference survival from month 2 to month 3 (logit scale)")                          # Chen 2026 Table S8 = 0.520 (SE 0.029, RSE 5.7%)
    dlgtsurv_04 <- fixed(0.401) ; label("Decrement in logit reference survival from month 3 to month 4 (logit scale)")                          # Chen 2026 Table S8 = 0.401 (SE 0.029, RSE 7.2%)
    dlgtsurv_05 <- fixed(0.301) ; label("Decrement in logit reference survival from month 4 to month 5 (logit scale)")                          # Chen 2026 Table S8 = 0.301 (SE 0.029, RSE 9.6%)
    dlgtsurv_06 <- fixed(0.258) ; label("Decrement in logit reference survival from month 5 to month 6 (logit scale)")                          # Chen 2026 Table S8 = 0.258 (SE 0.028, RSE 11.0%)
    dlgtsurv_07 <- fixed(0.189) ; label("Decrement in logit reference survival from month 6 to month 7 (logit scale)")                          # Chen 2026 Table S8 = 0.189 (SE 0.019, RSE 10.1%)
    dlgtsurv_08 <- fixed(0.196) ; label("Decrement in logit reference survival from month 7 to month 8 (logit scale)")                          # Chen 2026 Table S8 = 0.196 (SE 0.024, RSE 12.3%)
    dlgtsurv_09 <- fixed(0.159) ; label("Decrement in logit reference survival from month 8 to month 9 (logit scale)")                          # Chen 2026 Table S8 = 0.159 (SE 0.027, RSE 16.7%)
    dlgtsurv_10 <- fixed(0.152) ; label("Decrement in logit reference survival from month 9 to month 10 (logit scale)")                         # Chen 2026 Table S8 = 0.152 (SE 0.018, RSE 11.5%)
    dlgtsurv_11 <- fixed(0.143) ; label("Decrement in logit reference survival from month 10 to month 11 (logit scale)")                        # Chen 2026 Table S8 = 0.143 (SE 0.027, RSE 18.8%)
    dlgtsurv_12 <- fixed(0.114) ; label("Decrement in logit reference survival from month 11 to month 12 (logit scale)")                        # Chen 2026 Table S8 = 0.114 (SE 0.022, RSE 19.3%)
    dlgtsurv_13 <- fixed(0.099) ; label("Decrement in logit reference survival from month 12 to month 13 (logit scale)")                        # Chen 2026 Table S8 = 0.099 (SE 0.020, RSE 20.7%)
    dlgtsurv_14 <- fixed(0.168) ; label("Decrement in logit reference survival from month 13 to month 14 (logit scale)")                        # Chen 2026 Table S8 = 0.168 (SE 0.034, RSE 20.1%)
    dlgtsurv_15 <- fixed(0.102) ; label("Decrement in logit reference survival from month 14 to month 15 (logit scale)")                        # Chen 2026 Table S8 = 0.102 (SE 0.022, RSE 21.2%)
    dlgtsurv_16 <- fixed(0.098) ; label("Decrement in logit reference survival from month 15 to month 16 (logit scale)")                        # Chen 2026 Table S8 = 0.098 (SE 0.025, RSE 25.4%)
    dlgtsurv_17 <- fixed(0.094) ; label("Decrement in logit reference survival from month 16 to month 17 (logit scale)")                        # Chen 2026 Table S8 = 0.094 (SE 0.026, RSE 28.0%)
    dlgtsurv_18 <- fixed(0.083) ; label("Decrement in logit reference survival from month 17 to month 18 (logit scale)")                        # Chen 2026 Table S8 = 0.083 (SE 0.017, RSE 20.0%)
    dlgtsurv_19 <- fixed(0.115) ; label("Decrement in logit reference survival from month 18 to month 19 (logit scale)")                        # Chen 2026 Table S8 = 0.115 (SE 0.032, RSE 27.7%)
    dlgtsurv_20 <- fixed(0.187) ; label("Decrement in logit reference survival from month 19 to month 20 (logit scale)")                        # Chen 2026 Table S8 = 0.187 (SE 0.047, RSE 25.1%)
    dlgtsurv_21 <- fixed(0.051) ; label("Decrement in logit reference survival from month 20 to month 21 (logit scale)")                        # Chen 2026 Table S8 = 0.051 (SE 0.015, RSE 29.6%)
    dlgtsurv_22 <- fixed(0.093) ; label("Decrement in logit reference survival from month 21 to month 22 (logit scale)")                        # Chen 2026 Table S8 = 0.093 (SE 0.022, RSE 23.5%)
    dlgtsurv_23 <- fixed(0.105) ; label("Decrement in logit reference survival from month 22 to month 23 (logit scale)")                        # Chen 2026 Table S8 = 0.105 (SE 0.055, RSE 52.9%)
    dlgtsurv_24 <- fixed(0.135) ; label("Decrement in logit reference survival from month 23 to month 24 (logit scale)")                        # Chen 2026 Table S8 = 0.135 (SE 0.069, RSE 51.1%)
    dlgtsurv_27 <- fixed(0.137) ; label("Decrement in logit reference survival from month 24 to month 27 (logit scale)")                        # Chen 2026 Table S8 = 0.137 (SE 0.035, RSE 25.5%). Knot spacing widens to 3 months here
    dlgtsurv_30 <- fixed(0.225) ; label("Decrement in logit reference survival from month 27 to month 30 (logit scale)")                        # Chen 2026 Table S8 = 0.225 (SE 0.107, RSE 47.8%)
    dlgtsurv_33 <- fixed(0.260) ; label("Decrement in logit reference survival from month 30 to month 33 (logit scale)")                        # Chen 2026 Table S8 = 0.260 (SE 0.090, RSE 34.8%)
    dlgtsurv_36 <- fixed(0.442) ; label("Decrement in logit reference survival from month 33 to month 36 (logit scale)")                        # Chen 2026 Table S8 = 0.442 (SE 0.115, RSE 25.9%)
    dlgtsurv_42 <- fixed(0.517) ; label("Decrement in logit reference survival from month 36 to month 42 (logit scale)")                        # Chen 2026 Table S8 = 0.517 (SE 0.143, RSE 27.8%). Knot spacing widens to 6 months here
    dlgtsurv_48 <- fixed(0)     ; label("Decrement in logit reference survival from month 42 to month 48 (logit scale); held constant at zero")  # Source control stream psp470197-sup-0002-Supinfo2.ctl $THETA row 30: '0 FIX ; TIME48'. Absent from Table S8 because it was not estimated -- too few arms survive past 42 months once the FLAG = 1 / FLAG = 2 tail trimming is applied, so the reference curve is flat over 42-48 months

    # =========================================================================
    # COVARIATE EFFECTS ON log(HR) (Chen 2026 Table S8, rows 31-33).
    #
    # log(HR) is the log hazard ratio of an arm relative to the reference
    # curve, and enters as an exponent on the reference survival function:
    # S_arm(t) = S_ref(t)^exp(log(HR)). Because S = exp(-H), raising S to the
    # power exp(log(HR)) multiplies the cumulative hazard by the hazard ratio,
    # so this is an ordinary proportional-hazards model written on the
    # survival scale.
    #
    # Each effect below is on the LOG scale; exponentiating reproduces the
    # hazard ratios printed in Chen 2026 Table 3 for Model 010:
    #     exp(0.004)  = 1.004   ECOG1HR     (published 1.004)
    #     exp(0.769)  = 2.158   ECOG2,3HR   (published 2.158)
    #     exp(-0.633) = 0.531   IPSOSHR     (published 0.531)
    # =========================================================================

    e_ecog1    <- fixed( 0.004) ; label("Change in log(HR) for overall survival when an arm is 100% ECOG performance status 1 rather than 100% ECOG 0 (log hazard ratio)")          # Chen 2026 Table S8 'ECOG1 effect' = 0.004 (SE 0.454). exp(0.004) = 1.004, matching Table 3 Model 010 ECOG1HR 1.004 [0.412, 2.443]
    e_ecog23   <- fixed( 0.769) ; label("Change in log(HR) for overall survival when an arm is 100% ECOG performance status 2 or 3 rather than 100% ECOG 0 (log hazard ratio)")      # Chen 2026 Table S8 'ECOG2,3 effect' = 0.769 (SE 0.391). exp(0.769) = 2.158, matching Table 3 Model 010 ECOG2,3HR 2.158 [1.003, 4.639]
    e_ipsosctl <- fixed(-0.633) ; label("Change in log(HR) for overall survival in the IPSOS control arm relative to a historical single-agent-chemotherapy control arm (log hazard ratio)")  # Chen 2026 Table S8 'IPSOS control arm effect' = -0.633 (SE 0.114). exp(-0.633) = 0.531, matching Table 3 Model 010 IPSOSHR 0.531 [0.424, 0.664]

    # =========================================================================
    # BETWEEN-TRIAL RANDOM EFFECT (Chen 2026 Table S8, 'Between-trial SD').
    #
    # This is an MBMA BETWEEN-TRIAL variance component on log(HR), NOT
    # between-subject variability: it describes how far one trial arm's whole
    # survival curve sits from the model prediction, and it is the only
    # stochastic element in the model.
    #
    # Table S8 reports the STANDARD DEVIATION, 0.257. The source control
    # stream confirms the scale unambiguously: $OMEGA is '1 FIX', so ETA(1) is
    # standard normal, and the log(HR) contribution is written ETA(1)*OM1 with
    # OM1 = THETA(31). A THETA multiplying a unit-variance ETA is a standard
    # deviation, so no variance-versus-SD ambiguity arises here.
    #
    # nlmixr2 parameterises an eta by its VARIANCE, so the value below is
    # 0.257^2 = 0.066049.
    # =========================================================================

    eta_study ~ fixed(0.066049)
  })

  model({
    # =======================================================================
    # Reference survival curve S_ref(t), piecewise constant between knots.
    #
    # The random walk of the source $PRED block is evaluated here as a single
    # cumulative decrement, using the same indicator arithmetic the sibling
    # Franzese_2026_pdl1_nsclc_mbma uses for its monthly baseline hazards.
    # (t >= k) contributes the k-th decrement once the k-th knot has been
    # reached, so dlgtsurv_cum at any t equals the sum of every decrement up
    # to and including the most recent knot -- identical to the chained
    # Xk = X(k-1) - THETA(k) recursion.
    #
    # Between knots the curve is held FLAT. That is the model as fitted: the
    # source likelihood only ever evaluates the curve at the knot times, and a
    # step function is also what a Kaplan-Meier curve actually is. Times
    # between 24 and 27 months therefore report S_ref(24), and so on for the
    # wider late knots.
    # =======================================================================
    dlgtsurv_cum <-
      (t >=  2) * dlgtsurv_02 + (t >=  3) * dlgtsurv_03 +
      (t >=  4) * dlgtsurv_04 + (t >=  5) * dlgtsurv_05 +
      (t >=  6) * dlgtsurv_06 + (t >=  7) * dlgtsurv_07 +
      (t >=  8) * dlgtsurv_08 + (t >=  9) * dlgtsurv_09 +
      (t >= 10) * dlgtsurv_10 + (t >= 11) * dlgtsurv_11 +
      (t >= 12) * dlgtsurv_12 + (t >= 13) * dlgtsurv_13 +
      (t >= 14) * dlgtsurv_14 + (t >= 15) * dlgtsurv_15 +
      (t >= 16) * dlgtsurv_16 + (t >= 17) * dlgtsurv_17 +
      (t >= 18) * dlgtsurv_18 + (t >= 19) * dlgtsurv_19 +
      (t >= 20) * dlgtsurv_20 + (t >= 21) * dlgtsurv_21 +
      (t >= 22) * dlgtsurv_22 + (t >= 23) * dlgtsurv_23 +
      (t >= 24) * dlgtsurv_24 + (t >= 27) * dlgtsurv_27 +
      (t >= 30) * dlgtsurv_30 + (t >= 33) * dlgtsurv_33 +
      (t >= 36) * dlgtsurv_36 + (t >= 42) * dlgtsurv_42 +
      (t >= 48) * dlgtsurv_48

    lgtsurv0 <- lgtsurv_01 - dlgtsurv_cum

    # Before the first knot no interval has elapsed, so the reference curve is
    # at 1 by construction (S0LAST = 1 in the source $PRED block).
    surv0 <- (t < 1) * 1 + (t >= 1) * exp(lgtsurv0) / (1 + exp(lgtsurv0))

    # =======================================================================
    # log(HR) for this arm relative to the reference curve.
    #
    # ECOG performance status enters as the arm's COMPOSITION: the percentage
    # of the arm at each performance-status level, divided by 100. ECOG 2 and
    # ECOG 3 share one coefficient, matching
    #   ECOGEFF = COVECOG1*PS1/100 + COVECOG2*(PS2+PS3)/100
    # in the source $PRED block. An all-ECOG-0 historical control arm gives
    # log(HR) = 0 at the typical value of eta_study, which is the reference.
    # =======================================================================
    lhr <-
      e_ecog1    * (PS_ECOG_1_PCT / 100) +
      e_ecog23   * ((PS_ECOG_2_PCT + PS_ECOG_3_PCT) / 100) +
      e_ipsosctl * TRT_IPSOS_CONTROL +
      eta_study

    hr <- exp(lhr)

    # Proportional-hazards rescaling of the reference curve. surv is the
    # arm-level overall-survival probability at time t, a dimensionless number
    # in [0, 1]. No residual error is declared: Chen 2026 maximises a
    # user-supplied binomial -2LL over per-interval death counts and estimates
    # no $SIGMA (see description).
    surv <- surv0^hr
  })
}
attr(Chen_2026_nsclc_os_mbma, "message") <-
  paste(
    "MBMA of overall survival in platinum-ineligible advanced NSCLC (Chen 2026, Model 010;",
    "26 studies, 41 arms, 3637 participants). Semi-parametric proportional hazards with a",
    "non-parametric 30-knot reference survival curve on the logit scale.",
    "Inputs: PS_ECOG_1_PCT (%), PS_ECOG_2_PCT (%), PS_ECOG_3_PCT (%), TRT_IPSOS_CONTROL (binary).",
    "Outputs: surv0 = reference (historical control, all-ECOG-0) survival curve; hr = arm hazard",
    "ratio relative to that reference; surv = arm overall-survival curve.",
    "Time is in MONTHS and the supported window is 1-48 months. Simulation scope: per-ARM mean",
    "survival curves -- the random effect is between-trial, NOT between-subject, so this model",
    "cannot simulate individual-subject event times. All parameter values wrapped in fixed()",
    "because the model is a downstream user of the published fit, not a re-estimation of it.",
    sep = " "
  )
Chen_2026_nsclc_os_mbma
