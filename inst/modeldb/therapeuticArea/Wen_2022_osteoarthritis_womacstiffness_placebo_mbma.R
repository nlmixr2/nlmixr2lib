Wen_2022_osteoarthritis_womacstiffness_placebo_mbma <- function() {
  description <- paste0(
    "MBMA. Longitudinal model-based meta-analysis of the PLACEBO-arm ",
    "response on the WOMAC stiffness subscale (0-20 standardized scale: two ",
    "items scored 0-10) in randomized, double-blind, placebo-controlled ",
    "trials of ORALLY administered treatments for osteoarthritis, fitted to ",
    "the placebo arms of 96 trials (of 130 included trials, 12,673 ",
    "participants) published 1991-2022, using data up to week 36. There is ",
    "no drug, no dose and no PK layer: the arm-mean score falls from its ",
    "observed baseline as an exponential approach to a plateau, ",
    "womacstiffness = BASE - Emax * (1 - exp(-k * t)), with typical Emax 1.76 ",
    "points at the median arm baseline of 10.23 and k 0.327 /week. The arm ",
    "mean baseline score (SCORE_WOMAC_STIFFNESS) is both the starting value and ",
    "the only retained covariate, acting linearly on Emax: ",
    "Emax = 1.76 * (1 + 0.0836 * (SCORE_WOMAC_STIFFNESS - 10.23)). Between-study ",
    "variability is additive on Emax and exponential on k. The residual SD ",
    "of an arm mean is addSd / sqrt(N_ARM), so N_ARM (arm sample size) must ",
    "be supplied. Simulation scope is STUDY-ARM-MEAN placebo trajectories, ",
    "NOT individual patients. Companion models of the same paper: ",
    "Wen_2022_osteoarthritis_womacpain_placebo_mbma, ",
    "Wen_2022_osteoarthritis_womacfunction_placebo_mbma."
  )
  reference <- paste0(
    "Wen X, Luo J, Mai Y, Li Y, Cao Y, Li Z, Han S, Fu Q, Zheng Q, Ding L, ",
    "Zhang Z, Li L. Placebo Response to Oral Administration in ",
    "Osteoarthritis Clinical Trials and Its Associated Factors: A ",
    "Model-Based Meta-analysis. JAMA Netw Open. 2022;5(10):e2235060. ",
    "doi:10.1001/jamanetworkopen.2022.35060. PMC9552894. Structural, ",
    "random-effect and residual models: Supplement eMethods 3 Equations 1-4 ",
    "and the NONMEM code of eMethods 5 ('NONMEM codes of WOMAC stiffness ",
    "model'). Covariate equation: Supplement eResults Equation 10. Parameter ",
    "estimates: main-text Table 1, 'WOMAC stiffness' columns."
  )
  vignette <- "Wen_2022_osteoarthritis_placebo_mbma"

  # No drug, no dose events and no concentration. The units entries follow
  # the placeholder convention of the arm-level MBMA siblings
  # (Serrano_2026_atopicDermatitis_placebo_mbma, Yu_2022_crohns_dcdai_mbma)
  # so checkModelConventions() sees a parseable dosing / concentration pair.
  units <- list(
    time = "week (weeks since randomization; k is reported in 1/week in eResults, and only data up to week 36 were modelled)",
    dosing = "n/a (placebo arms only; this model consumes NO rxode2 dose events and has no exposure driver)",
    concentration = "points/arm (womacstiffness is the STUDY-ARM mean WOMAC stiffness subscale score on the 0-20 standardized scale; it is NOT a drug concentration)"
  )

  covariateData <- list(
    SCORE_WOMAC_STIFFNESS = list(
      description = "Study-arm mean BASELINE WOMAC stiffness subscale score on the 0-20 standardized scale (two items, each scored 0-10).",
      units = "(score, 0-20)",
      type = "continuous",
      reference_category = NULL,
      notes = "TRIAL-ARM-LEVEL, not subject-level. Plays two roles, exactly as BASE does in the eMethods 5 NONMEM code: it is the model's starting value (womacstiffness at time 0 equals SCORE_WOMAC_STIFFNESS) and it scales the typical Emax linearly, Emax = 1.76 * (1 + 0.0836 * (SCORE_WOMAC_STIFFNESS - 10.23)) (eResults Equation 10). 10.23 is the median arm baseline across the 96 trials reporting WOMAC stiffness (Results, 'Characteristics of the Included Studies': range 2.10-16.75, Q1-Q3 7.81-12.00). The source standardized every trial's answers to 0-10 per item (Methods, 'Data Extraction'), so a trial reported on the 0-4 Likert version is multiplied by 2.5 and one reported on a 0-100 mm VAS per item is divided by 10 (the paper's own conversion: Likert 2 / 4 / 6 = 5 / 10 / 15). The covariate must be held constant over an arm's records. The typical Emax stays positive over the whole 0-20 scale (it would turn negative only below 10.23 - 1/0.0836 = -1.7).",
      source_name = "BASE (eMethods 5 NONMEM $INPUT and $PRED); 'Baseline_stiffness' (eResults Equation 10)"
    ),
    N_ARM = list(
      description = "Number of patients in the placebo arm contributing to the arm-mean score.",
      units = "participants",
      type = "count",
      reference_category = NULL,
      notes = "Study-design quantity supplied per observation row; the meta-analytic weight of the residual. eMethods 3 Equation 4 and the eMethods 5 code (W = 1 / SQRT(SIZE); Y = EFT + W * ERR(1)) divide the residual SD by sqrt(N). It does not affect the typical-value or between-study prediction.",
      source_name = "SIZE (eMethods 5 NONMEM $INPUT); N_i,j (eMethods 3 Equation 4)"
    )
  )

  # Screened in the source's stepwise covariate search (eMethods 3,
  # 'Covariate Model Establishment') and not retained. Only covariates that
  # already have a register canonical are listed; the rest are named in
  # population$notes.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Study-arm mean patient age. Screened but NOT retained in the final model.",
      units = "year",
      type = "continuous",
      reference_category = NULL,
      notes = "Listed among the examined covariates in eMethods 3; the final model retains only the baseline score (Results, 'Model Establishment and Assessment'). No point estimate is reported. Overall mean age 59.9 years (Results).",
      source_name = "AGE (eMethods 5 NONMEM $INPUT)"
    ),
    SEXF_PCT = list(
      description = "Percentage of female patients in the study arm. Screened (as the male ratio) but NOT retained.",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "The source screened the male ratio (MALE column; SEXF_PCT = 100 - male %). Not retained; no point estimate. Overall 68.9% women (Results).",
      source_name = "MALE (eMethods 5 NONMEM $INPUT)"
    ),
    YEAR_PUB = list(
      description = "Publication year of the trial. Screened but NOT retained in the final model.",
      units = "year",
      type = "continuous",
      reference_category = NULL,
      notes = "Publication year. Not retained as a model covariate; the post hoc subgroup analysis (Figure 3, eFigure 6) shows a trend toward larger placebo responses in trials published 2011-2022, which the model does not carry.",
      source_name = "YEAR (eMethods 5 NONMEM $INPUT)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 12673L,
    n_studies = 96L,
    age_range = "overall mean 59.9 years (Results); arm means in eTable 4",
    sex_female_pct = 68.9,
    race_ethnicity = "Proportion of White patients reported by some trials only (eTable 4); it was the only racial or ethnic group with complete enough data for the subgroup analysis, which found no difference (Figure 3 / eFigure 6).",
    disease_state = "Osteoarthritis (knee and/or hip) in randomized, double-blind, placebo-controlled trials whose interventions and placebo were given ORALLY and that reported at least one WOMAC subscale. Arm-mean baseline WOMAC stiffness 2.10-16.75 (median 10.23; Q1-Q3 7.81-12.00) on the 0-20 standardized scale.",
    dose_range = "n/a (placebo arms only; placebo forms were caplet, tablet or pill vs powder, package or granule)",
    regions = "International; PubMed, EMBASE and the Cochrane Library searched from 1 January 1991 to 2 July 2022 (eTable 1).",
    notes = "MBMA at the STUDY-ARM level: 130 trials with 12,673 participants were included, of which 96 reported WOMAC stiffness. Only records up to week 36 were modelled because just 10 trials reported data beyond 36 weeks (Results). The random effects are BETWEEN-STUDY, not between-subject. Covariates screened and not retained (eMethods 3): age, male ratio, Kellgren-Lawrence grade, proportion with prior NSAID use, proportion with prior supplement use, proportion of White patients, placebo form, administration frequency, publication year, funding source and literature quality; covariates with more than 30% missing were not evaluated and those with up to 30% missing were median-imputed. Fitted in NONMEM 7.4 with FOCE-I. The placebo response includes regression to the mean and natural disease course (Discussion limitations)."
  )

  ini({
    # ========================================================================
    # eMethods 5, 'NONMEM codes of WOMAC stiffness model', $PRED:
    #   EMCOV = (1 + THETA(3) * (BASE - 10.23))
    #   TVEM  = EMCOV * THETA(1)
    #   EM    = TVEM + ETA(1)
    #   K     = THETA(2) * EXP(ETA(2))
    #   W     = 1 / SQRT(SIZE)
    #   EFT   = BASE - EM * (1 - EXP(-k * TIME))
    #   Y     = EFT + W * ERR(1)
    # Final estimates are main-text Table 1, 'WOMAC stiffness' 'Estimates (RSE %)'.
    # ========================================================================
    emax <- 1.76
    label("Typical maximum placebo response (decrease in WOMAC stiffness) at the median arm baseline of 10.23 (paper: Emax; points, 0-20 scale)") # Table 1 WOMAC stiffness 'Emax' = 1.76 (RSE 8.80%); eMethods 5 THETA(1)

    lkpbo <- log(0.327)
    label("Log onset rate constant of the placebo response (paper: K; back-transform 0.327 /week)") # Table 1 WOMAC stiffness 'K' = 0.327 (RSE 19.7%); eResults 'onset rate (k) ... 0.327 ... week-1'; eMethods 5 THETA(2)

    e_score_womac_stiffness_emax <- 0.0836
    label("Fractional change in Emax per point of arm baseline WOMAC stiffness above 10.23 (paper: theta baseline; 1/point)") # Table 1 WOMAC stiffness 'theta baseline' = 0.0836 (RSE 19.4%); eResults Equation 10; eMethods 5 THETA(3)

    # Between-STUDY random effects (eMethods 3 Equations 2-3): additive on
    # Emax, exponential on k. Table 1 reports the square roots of the OMEGA
    # diagonal: 'eta Emax' in score points and 'eta k, %' as 100 * omega.
    eta_study_emax ~ 1.9881 # Table 1 WOMAC stiffness 'eta Emax' = 1.41 (RSE 8.90%), an SD in points; variance = 1.41^2
    eta_study_lkpbo ~ 1.1881 # Table 1 WOMAC stiffness 'eta k, %' = 109 (RSE 13.4%), read as 100 * omega; variance = 1.09^2

    # Residual (eMethods 3 Equation 4; eMethods 5 Y = EFT + W * ERR(1)).
    addSd <- 1.74
    label("Residual SD of the arm-mean WOMAC stiffness for an arm of 1 patient; the per-record SD is addSd / sqrt(N_ARM) (points)") # Table 1 WOMAC stiffness 'epsilon' = 1.74 (RSE 10.9%), an SD in points
  })

  model({
    # eResults Equation 10 / eMethods 5 EMCOV and TVEM; ETA(1) additive.
    emax_i <- emax * (1 + e_score_womac_stiffness_emax * (SCORE_WOMAC_STIFFNESS - 10.23)) + eta_study_emax

    # eMethods 3 Equation 3; ETA(2) exponential.
    kpbo <- exp(lkpbo + eta_study_lkpbo)

    # eMethods 3 Equation 1: the arm-mean score falls from its baseline.
    womacstiffness <- SCORE_WOMAC_STIFFNESS - emax_i * (1 - exp(-kpbo * time))

    # eMethods 3 Equation 4: the residual SD of an arm mean is eps / sqrt(N).
    sdArm <- addSd / sqrt(N_ARM)
    womacstiffness ~ add(sdArm)
  })
}
