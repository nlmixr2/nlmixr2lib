Oniki_2018_bmi <- function() {
  description <- paste0(
    "Population prediction model for body mass index (BMI, kg/m^2) in ",
    "341 elderly Japanese health-screening participants (Oniki 2018). ",
    "BMI is parameterised as the typical value at age 70.8 years for a ",
    "male non-T/T reference subject, multiplied by a power-of-age scalar ",
    "and a female sex multiplier, with an additive shift for carriers of ",
    "the DsbA-L (GSTK1) rs1917760 -1308G>T T/T genotype (Eq. 1 of the ",
    "source paper). Log-normal between-subject variability acts on the ",
    "typical BMI, and log-normal (exponential) residual error is added ",
    "on the linear scale. The fit was performed with NONMEM 7.2.0 ",
    "$PRED METHOD=COND INTER (s010 control stream); MINIMIZATION ",
    "SUCCESSFUL, OBJ 1454.644. There is no drug input. Companion ",
    "BMI-driven NAFLD-risk model in Oniki_2018_nafld_risk."
  )
  reference <- paste(
    "Oniki K, Watanabe T, Kudo M, Izuka T, Ono T, Matsuda K, Sakamoto Y,",
    "Nagaoka K, Imafuku T, Ishima Y, Watanabe H, Maruyama T, Otake K,",
    "Ogata Y, Saruwatari J.",
    "Modeling of the Weight Status and Risk of Nonalcoholic Fatty Liver",
    "Disease in Elderly Individuals: The Potential Impact of the Disulfide",
    "Bond-Forming Oxidoreductase A-Like Protein (DsbA-L) Polymorphism on",
    "the Weight Status.",
    "CPT Pharmacometrics Syst Pharmacol. 2018 Jun;7(6):384-393.",
    "doi:10.1002/psp4.12292.",
    "Companion NAFLD-risk model in Oniki_2018_nafld_risk."
  )
  vignette <- "Oniki_2018_BMI_NAFLD"
  units <- list(
    time          = "year",
    dosing        = "n/a (population BMI-prediction model; no drug input)",
    concentration = "BMI (kg/m^2)"
  )

  covariateData <- list(
    AGE = list(
      description        = "Subject age at the time of the BMI observation (years). The Oniki 2018 longitudinal NONMEM dataset carries AGE per record; the BMI model uses the per-record value. Because the typical BMI changes by < 1 kg/m^2 across the 5.5-year follow-up window for this elderly cohort (the age power exponent is -0.0709 centred on 70.8 years), the practical distinction between baseline and per-record age is small.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centring value 70.8 years (Oniki 2018 Eq. 1). The power-form covariate effect is (AGE / 70.8)^e_age_bmi; a typical 80-year-old is predicted to weigh ~0.85% less than a 70.8-year-old at the same sex and genotype.",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Sex indicator; 1 = female, 0 = male. The Oniki 2018 dataset GENDER column already encodes 0 = male / 1 = female, matching the canonical SEXF orientation without inversion.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Used as a multiplicative scalar on the typical BMI (a = e_sexf_bmi when SEXF = 1, a = 1 when SEXF = 0); the female typical BMI is 96.8% of the male typical BMI at the same age and genotype (Oniki 2018 Eq. 1).",
      source_name        = "GENDER"
    ),
    DSBAL_TT = list(
      description        = "Indicator for the DsbA-L (GSTK1) rs1917760 -1308G>T T/T genotype; 1 = subject carries the T/T genotype, 0 = subject carries the G/G or G/T genotype (the pooled reference). Time-fixed per subject (germline genotype).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (G/G or G/T, pooled)",
      notes              = "Derive from the source DsbAL three-level column (0 = G/G, 1 = G/T, 2 = T/T) as DSBAL_TT = as.integer(DsbAL == 2). G/G and G/T are pooled in Oniki 2018 because the T/T allele is the functional minor-allele state most strongly associated with elevated BMI. Used as an additive +1.5 kg/m^2 shift on the typical BMI (Oniki 2018 Eq. 1).",
      source_name        = "DsbAL"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 341L,
    n_studies      = 1L,
    n_observations = "2015 BMI records across 342 subject-ids in the NONMEM dataset (Oniki 2018 s010 .lst NO. OF DATA RECS); 341 unique subjects after excluding those with habitual alcohol intake or hepatitis B/C virus positivity (Oniki 2018 Methods, Subjects and study protocol).",
    age_range      = "Elderly Japanese cohort; baseline age mean 67.7 years (SD ~5.9) pooled across DsbA-L genotypes (Oniki 2018 Table 1).",
    age_median     = NA_character_,
    weight_range   = NA_character_,
    weight_median  = NA_character_,
    sex_female_pct = 100 * (80 + 56 + 9) / (192 + 129 + 20),
    race_ethnicity = c(Asian = 100),
    disease_state  = paste0(
      "General elderly Japanese cohort participating in the Japanese ",
      "Red Cross Kumamoto Health Care Center elderly health-screening ",
      "program; baseline mean BMI 22.5-23.9 kg/m^2 across DsbA-L ",
      "genotype strata (Oniki 2018 Table 1)."
    ),
    dose_range     = "n/a (no drug input; population disease-risk model)",
    regions        = "Japan (Kumamoto)",
    notes          = paste0(
      "Retrospective longitudinal observation, 5.5 +/- 1.1 years of ",
      "follow-up. The Japanese Red Cross Kumamoto Health Care Center ",
      "screening program collected baseline demographics (Oniki 2018 ",
      "Table 1) and annual BMI and laboratory values. Genotyping was ",
      "performed on DsbA-L (GSTK1) rs1917760 and PNPLA3 rs738409 ",
      "via TaqMan allelic-discrimination assay (Oniki 2018 Methods, ",
      "Genotyping). Subjects with habitual alcohol intake (> 30 g/day ",
      "in men, > 20 g/day in women) or positive hepatitis B/C serology ",
      "were excluded per the Japanese NAFLD practical guidelines."
    )
  )

  ini({
    # ============================================================
    # All structural-parameter values are from the s010 NONMEM
    # control stream's FINAL PARAMETER ESTIMATE block (in agreement
    # with Eq. 1 of the published Oniki 2018 paper). MINIMIZATION
    # SUCCESSFUL; OBJ 1454.644. The covariance step ran successfully
    # (no covariance-step caveat for this sub-model).
    # ============================================================

    # ----- Structural parameters (Oniki 2018 Eq. 1; s010 .lst FINAL PARAMETER ESTIMATE) -----
    e0_bmi          <- 22.6      ; label("Typical BMI at age 70.8 for the male non-T/T reference (kg/m^2)")  # s010 .lst THETA(1) E0; SE 0.223
    e_age_bmi       <- -0.0709   ; label("Power exponent for (AGE / 70.8) on typical BMI (unitless)")         # s010 .lst THETA(2) Age; SE 0.015
    e_sexf_bmi      <- 0.968     ; label("Multiplicative factor on typical BMI for female (vs male reference; unitless)")  # s010 .lst THETA(3) Gender; SE 0.0133
    e_dsbal_tt_bmi  <- 1.50      ; label("Additive shift on typical BMI for DsbA-L T/T (vs G/G or G/T pooled reference; kg/m^2)")  # s010 .lst THETA(4) DsbA-L; SE 0.693

    # ----- Between-subject variability -----
    # NONMEM: E1 = E0 * EXP(ETA(1)); $OMEGA 1.50E-02 ; ETA(1)
    # Log-normal IIV on the typical BMI; variance on the natural-log
    # scale is 0.0150, corresponding to roughly 12.3% CV on the
    # back-transformed BMI scale (sqrt(exp(0.0150) - 1) ~= 0.123).
    etae0_bmi ~ 0.0150  # s010 .lst OMEGA(1,1); SE 9.87e-4

    # ----- Residual error -----
    # NONMEM: Y = E1 * EXP(EPS(1)); $SIGMA 6.59E-04 ; EPS(1).
    # Log-normal residual on BMI; the source-reported correlation-form
    # diagnostic (2.57%) equals sqrt(6.59e-4), the SD on the
    # log-normal scale. nlmixr2's `~ lnorm(expSd)` matches the source's
    # Y = E1 * EXP(EPS(1)) form exactly.
    expSd <- sqrt(6.59e-4)  ; label("Log-normal residual SD on BMI (unitless, on natural-log scale; ~ 2.57%)")  # s010 .lst SIGMA(1,1); SE 1.02e-5
  })

  model({
    # ----- Sex multiplier (Oniki 2018 Eq. 1; s010 .lst $PRED IF (GENDER...)) -----
    # In NONMEM: A = 1 if GENDER == 0 (male), A = THETA(3) if GENDER == 1 (female).
    # The (1 - SEXF) * 1 + SEXF * e_sexf_bmi pattern compiles cleanly in
    # rxode2 (no IF/ELSE) and gives a = 1 when SEXF = 0 and a = e_sexf_bmi
    # when SEXF = 1.
    sex_mult <- (1 - SEXF) * 1 + SEXF * e_sexf_bmi

    # ----- DsbA-L additive shift (Oniki 2018 Eq. 1; s010 .lst $PRED IF (DsbAL...)) -----
    # NONMEM: B = 0 if DsbAL <= 1 (G/G or G/T), B = THETA(4) if DsbAL == 2 (T/T).
    # The DSBAL_TT binary indicator is constructed upstream as
    # as.integer(DsbAL == 2), so multiplying by e_dsbal_tt_bmi reproduces
    # the IF/ELSE branch directly.
    dsbal_shift <- DSBAL_TT * e_dsbal_tt_bmi

    # ----- Typical BMI for this subject (Oniki 2018 Eq. 1) -----
    # bmi_typical = e0_bmi * sex_mult * (AGE / 70.8)^e_age_bmi + dsbal_shift
    bmi_typical <- e0_bmi * sex_mult * (AGE / 70.8)^e_age_bmi + dsbal_shift

    # ----- Individual BMI with log-normal between-subject variability -----
    # NONMEM: E1 = E0 * EXP(ETA(1)).
    bmi <- bmi_typical * exp(etae0_bmi)

    # ----- Residual error -----
    # NONMEM: Y = E1 * EXP(EPS(1)) -> rxode2 log-normal residual.
    bmi ~ lnorm(expSd)
  })
}
