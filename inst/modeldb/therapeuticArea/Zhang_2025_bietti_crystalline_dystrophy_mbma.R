Zhang_2025_bietti_crystalline_dystrophy_mbma <- function() {
  description <- paste0(
    "MBMA. Linear natural-history (disease-progression) model for the change ",
    "from baseline in best-corrected visual acuity (BCVA, LogMAR scale) in ",
    "untreated Bietti crystalline corneoretinal dystrophy (BCD), an ",
    "autosomal-recessive retinal degeneration caused by biallelic CYP4V2 ",
    "mutations. Zhang 2025 pooled INDIVIDUAL eye-level BCVA time courses ",
    "digitised from 14 published studies (117 patients, 193 study eyes) and ",
    "fit a single-parameter structural model, deltaBCVA = slope * time ",
    "(supplement equation 1, where the supplement calls the slope K). ",
    "Higher LogMAR = worse vision, so a positive slope encodes vision LOSS. ",
    "Unusually for an MBMA the random effect is BETWEEN-EYE, not ",
    "between-study: because individual participant data were reconstructed ",
    "from each source publication, the paper estimated an exponential ",
    "between-individual random effect on the slope (supplement equation 2, ",
    "omega = 0.898 on the log scale) rather than a study-level effect, and ",
    "reported no study-level variance term. The model is therefore suitable ",
    "for simulating eye-level trajectories, but it carries no study effect ",
    "and cannot separate between-study from between-eye heterogeneity. ",
    "Residual error is additive on the LogMAR scale (supplement equation 4). ",
    "There is no drug input and no dosing: the model describes untreated ",
    "natural history and was built to support the design and efficacy ",
    "evaluation of CYP4V2 gene-therapy trials. Covariate screening (age, age ",
    "at onset, disease duration, baseline BCVA, sex, race, CYP4V2 genotype, ",
    "family history) retained NOTHING, so the base model is the final model ",
    "(supplement appendix 6). The paper additionally reports post-hoc ",
    "empirical-Bayes subgroup slopes (0.0766 for baseline BCVA < 0.5 LogMAR, ",
    "0.0979 for baseline BCVA >= 0.5 LogMAR, 0.0675 for onset age < 40 y, ",
    "0.0892 for onset age >= 40 y, and 0.0900 for baseline BCVA >= 0.5 ",
    "LogMAR combined with disease duration >= 10 y); these are summaries of ",
    "individual slopes, not fitted covariate effects, and are reproduced in ",
    "the validation vignette by overriding lslope rather than encoded here."
  )
  reference <- paste(
    "Zhang H, Yin S, Guan N, Wang J, Cheng Q, Zhang L, Zheng Q, Lv H, Wei W.",
    "Natural history of progressive vision loss in Bietti crystalline",
    "dystrophy: a model-based meta-analysis.",
    "BMJ Open Ophthalmology 2025;10:e001908.",
    "doi:10.1136/bmjophth-2024-001908.",
    sep = " "
  )
  vignette <- "Zhang_2025_bietti_crystalline_dystrophy_mbma"

  units <- list(
    time          = "year",
    dosing        = "n/a (no dosing; untreated natural history)",
    concentration = "n/a (no drug concentration)",
    response      = "LogMAR (change from baseline in best-corrected visual acuity; observation deltaBCVA)"
  )

  # Every covariate Zhang 2025 examined was screened OUT: no covariate
  # reduced the objective function by the forward-inclusion threshold of
  # 3.84 (supplement appendix 6), so the base model is the final model and
  # model() references no covariate at all. They are recorded here as
  # documentation of the paper's covariate screen. Age at onset, disease
  # duration, CYP4V2 genotype and family history have no entry in
  # inst/references/covariate-columns.md and no canonical name is minted
  # for them here (they are not used by any model); they are described in
  # population$notes and in the validation vignette instead.
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Age of the subject at the first BCVA recording.",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on the slope in Zhang 2025 supplement appendix 6 run 201 (dOFV -3.131, short of the -3.84 forward-inclusion threshold) and not retained. Cohort median 49.0 y [15.0, 76.0] (supplement appendix 2).",
      source_name        = "age"
    ),
    SCORE_BCVA = list(
      description        = "Baseline best-corrected visual acuity. NOTE: Zhang 2025 records BCVA on the LogMAR scale (higher = worse vision), NOT on the ETDRS-letter scale (0-100, higher = better) that the covariate register's SCORE_BCVA entry defines. Because the covariate was screened out and is never referenced in model(), no scale reconciliation was required and none was performed.",
      units              = "LogMAR",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The single closest covariate to significance: Zhang 2025 supplement appendix 6 run 203 gave dOFV -3.969, which passed the forward-inclusion threshold of -3.84, but backward elimination (run 209) returned dOFV +3.969, short of the 6.63 retention threshold, so it was dropped. Cohort median 0.150 LogMAR [-0.200, 2.60] (supplement appendix 2).",
      source_name        = "bcva"
    ),
    SEXF = list(
      description        = "Sex. 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male).",
      notes              = "Screened on the slope in Zhang 2025 supplement appendix 6 run 208 (dOFV -0.568) and not retained. Cohort 77/117 female (65.8%) (Zhang 2025 table 1).",
      source_name        = "sex"
    ),
    RACE_ASIAN_NORTHEAST = list(
      description        = "East Asian heritage indicator. 1 = East Asian, 0 = non-East Asian.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-East Asian).",
      notes              = "Screened on the slope in Zhang 2025 supplement appendix 6 run 205 ('ethnicity on K', dOFV -0.829) and not retained. Cohort 80/117 East Asian (68.4%) (Zhang 2025 table 1). Zhang 2025 dichotomises race as East Asian vs non-East Asian, which matches the RACE_ASIAN_NORTHEAST definition (worldwide Chinese, Japanese or Korean heritage).",
      source_name        = "ethnicity"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 117L,
    n_eyes         = 193L,
    n_studies      = 14L,
    age_range      = "15.0-76.0 years (Zhang 2025 supplement appendix 2, age at the first BCVA recording)",
    age_median     = "49.0 years (Zhang 2025 supplement appendix 2)",
    weight_range   = NA_character_,
    weight_median  = NA_character_,
    sex_female_pct = 65.8,
    race_ethnicity = c(`East Asian` = 68.4, `Non-East Asian` = 31.6),
    disease_state  = "Bietti crystalline corneoretinal dystrophy (BCD), an autosomal-recessive progressive retinal degeneration caused by biallelic CYP4V2 mutations. All eyes untreated (a study-inclusion criterion was that the population had not previously been treated); there is no approved therapy for BCD.",
    dose_range     = "n/a (untreated natural history; the model has no drug input)",
    regions        = "Six of the 14 included studies enrolled East Asian populations (China, Japan, Korea) and nine enrolled non-East Asian populations; one study contributed patients to both strata, which is why the two study counts sum to 15 rather than 14 (Zhang 2025 Results and table 1).",
    onset_age_range  = "11.0-76.0 years (Zhang 2025 supplement appendix 2; where the age of onset was missing the age at the first BCVA recording was substituted)",
    onset_age_median = "47.0 years (Zhang 2025 supplement appendix 2)",
    disease_duration_range  = "0.500-47.0 years (Zhang 2025 supplement appendix 2; duration = age - age of onset + longest visit)",
    disease_duration_median = "8.00 years (Zhang 2025 supplement appendix 2)",
    baseline_bcva  = "Median 0.150 LogMAR [-0.200, 2.60] (Zhang 2025 supplement appendix 2). Of the 193 study eyes, 151 had baseline BCVA < 0.5 LogMAR and 42 had >= 0.5 LogMAR (supplement appendix 10).",
    notes          = paste0(
      "Individual eye-level BCVA time courses were reconstructed from 14 ",
      "published studies retrieved by a PRISMA-style search of PubMed, ",
      "MEDLINE, EMBASE, CINAHL and CNKI updated 20 November 2023 (Zhang ",
      "2025 Methods). Values reported only in figures were digitised with ",
      "xyscan v4.1; Snellen and decimal acuities were converted to LogMAR; ",
      "off-chart acuities were assigned 1.9 (counting fingers), 2.3 (hand ",
      "motion), 2.7 (light perception) and 3.0 (no light perception) ",
      "LogMAR. The ANALYSIS UNIT is the study eye (193 eyes from 117 ",
      "patients; 84 left and 109 right), so the between-individual random ",
      "effect etalslope is between-EYE, and eyes from the same patient are ",
      "treated as independent. Screened-but-not-retained covariates with no ",
      "canonical register name: age at onset (supplement appendix 6 run ",
      "202, dOFV -2.412), disease duration (run 204, dOFV -0.139), CYP4V2 ",
      "genotype coded as c.802-8_810de117insGC [exon 7 deletion] homozygous ",
      "or compound heterozygous vs other vs unknown (run 207, dOFV -0.200; ",
      "53 / 36 / 28 patients), and family history of BCD (run 206, dOFV ",
      "-0.663; 11 of 117 patients, 9.4%). Study quality was graded with the ",
      "JBI critical appraisal checklist for case reports (supplement ",
      "appendix 12) because most included studies are retrospective case ",
      "series or case reports with no control group."
    )
  )

  ini({
    # Structural model (Zhang 2025 supplement, "Data analysis", equation 1):
    #   Effect = K * Time
    # a single-parameter linear disease-progression model in which "Effect"
    # is the change from baseline in BCVA (LogMAR) and K is the slope.
    # LogMAR is a log-scale acuity measure on which HIGHER = WORSE vision,
    # so a positive slope encodes progressive vision loss.
    #
    # The slope is log-transformed here because the source paper places an
    # EXPONENTIAL between-individual random effect on it (supplement
    # equation 2, P_i = P_TV * exp(eta_i); supplement appendix 4 model 102,
    # "ETA: Exponential type", selected as the base model). log() of the
    # published linear-scale estimate is the exact nlmixr2 equivalent.
    lslope <- log(0.0566); label("Log rate of BCVA loss on the LogMAR scale (LogMAR/year)")
    # Zhang 2025 supplement appendix 7, "Final model / Parameter of population / K": estimate 0.0566 (RSE 10.1%), 95% CI 0.0454 to 0.0678. Bootstrap (1000 successes) median 0.0564, 95% CI 0.0452 to 0.0694.

    # Between-EYE (not between-study) exponential random effect on the
    # slope. The supplement defines the reported omega as a STANDARD
    # DEVIATION, not a variance: "eta_i is a random effect between
    # individuals, conforming to a normal distribution with mean 0 and
    # variance omega^2" (supplement text following equations 2 and 3).
    # nlmixr2's ini() takes the VARIANCE, hence the square. Written as
    # 0.898^2 rather than 0.806404 so the published number stays legible in
    # the source trace.
    etalslope ~ 0.898^2
    # Zhang 2025 supplement appendix 7, "Between-individual random effect parameters / omega (K)": estimate 0.898 (RSE 8.7%), 95% CI 0.745 to 1.05. Bootstrap median 0.894, 95% CI 0.700 to 1.08. omega is an SD per the supplement's own definition (variance = omega^2), so the ini() variance is 0.898^2 = 0.806404.

    # Additive residual error on the LogMAR scale (supplement equation 4,
    # Y_obs = Y_pred + eps; supplement appendix 4 model 102, "SIGMA:
    # Additive type"). As with omega, the supplement defines the reported
    # sigma as a STANDARD DEVIATION: "eps_i,1 and eps_i,2 conform to a
    # normal distribution with mean 0 and variance sigma_i,1^2 and
    # sigma_i,2^2" (supplement text following equations 4 to 6). nlmixr2's
    # add() also takes an SD, so the published value is used unchanged.
    addSd <- 0.336; label("Additive residual error standard deviation on the LogMAR scale (LogMAR)")
    # Zhang 2025 supplement appendix 7, "Intra-individual random effect parameters / sigma (ADD)": estimate 0.336 (RSE 5.3%), 95% CI 0.301 to 0.371. Bootstrap median 0.336, 95% CI 0.205 to 0.446.
  })

  model({
    # Individual (per-eye) rate of BCVA loss. Exponential between-individual
    # random effect per Zhang 2025 supplement equation 2:
    #   K_i = K_TV * exp(eta_i)
    # Because exp() is strictly positive, every simulated eye deteriorates
    # monotonically; the observed BCVA improvements in the source data
    # (changes from baseline as low as -0.30 LogMAR, supplement appendix 3)
    # are absorbed entirely by the additive residual error.
    slope <- exp(lslope + etalslope)

    # Zhang 2025 supplement equation 1: Effect = K * Time. There is no
    # intercept term, so the predicted change from baseline is exactly 0 at
    # time 0 by construction (time 0 is the first visual-acuity evaluation,
    # i.e. the baseline visit). Dimensional check:
    # (LogMAR/year) * year = LogMAR.
    deltaBCVA <- slope * time

    deltaBCVA ~ add(addSd)
  })
}
