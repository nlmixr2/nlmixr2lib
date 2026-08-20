Garcia_2025_garadacimab_hae_attack <- function() {
  description <- "Repeated-time-to-event (RTTE) exposure-response model for hereditary angioedema (HAE) attacks under the anti-activated-factor-XII (anti-FXIIa) monoclonal antibody garadacimab, from Garcia 2025. The log hazard is the sum of a constant baseline (run-in) hazard, a constant on-treatment (placebo/study) effect, and an Imax inhibitory effect driven by continuous-time garadacimab plasma concentration: h(t) = h0 * (1 - iplac * ON_TREATMENT) * (1 - imax * Cc^hill / (ec50^hill + Cc^hill)). Imax and the Hill coefficient are fixed. The two-compartment garadacimab PK model of the companion PopPK analysis is embedded verbatim so the hazard is driven by simulated concentration, exactly as in the source NONMEM ADVAN13 control stream. The cumulative hazard is carried as an ODE state and exposed together with the survivor function and the relative risk of attack versus the untreated run-in period. Sister model file from the same paper (PopPK plus FXIIa-mediated kallikrein activity PD): modellib('Garcia_2025_garadacimab')."
  reference <- paste(
    "Garcia R, Cheng S, Glassman F, Sharma A, De Miguel-Lillo B, Wiens M,",
    "Johnston C, Lawo JP, Pragst I, French J, Polhamus D, Nandy P.",
    "Population pharmacokinetic/pharmacodynamic and exposure-response modeling",
    "of garadacimab in healthy volunteers and patients with hereditary angioedema.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(5):954-963.",
    "doi:10.1002/psp4.70009.",
    "Open Access under CC BY-NC 4.0.",
    "The hazard equations and the NONMEM control stream are in Methods S2;",
    "final parameter estimates are in Tables S8-S9 of the supplement",
    "(Wiley file PSP4-14-954-s001.docx); the embedded PK parameters are from",
    "Tables S4-S5.",
    "Sister model file from the same paper:",
    "modellib('Garcia_2025_garadacimab').",
    sep = " "
  )
  vignette <- "Garcia_2025_garadacimab"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline total body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Acts on the EMBEDDED PK model only (power exponents on CL and Q, and on Vc and Vp, both centred at 70 kg), and thereby on the hazard indirectly through garadacimab concentration. Body weight was explicitly REMOVED as a direct predictor of the ER parameters: the sampling-importance-resampling CIs identified baseline body weight as the only covariate containing the null effect for both ER parameters, and because it was correlated with exposure it was dropped (Results Section 3.5). This matches the Methods S2 control stream, where THETA(6) (weight on baseline hazard) and THETA(7) (weight on EC50) are both 0 FIX. Observed range 43.3-153 kg, median 79.2 kg.",
      source_name        = "BLWT"
    ),
    RACE_JAPANESE = list(
      description        = "Japanese-heritage indicator: 1 = Japanese, 0 = non-Japanese.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Japanese)",
      notes              = "Time-fixed per subject. Acts on the EMBEDDED PK model only (multiplicative effect on CL). Its direct ER effects (THETA(12) on the baseline hazard and THETA(13) on EC50) are estimated in the Methods S2 control stream but are not reported -- see covariatesDataExcluded.",
      source_name        = "JPN"
    ),
    RACE_CHINESE = list(
      description        = "Chinese-heritage indicator: 1 = Chinese, 0 = non-Chinese.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Chinese)",
      notes              = "Time-fixed per subject. Acts on the EMBEDDED PK model only (multiplicative effect on CL). Its direct ER effects (THETA(14), THETA(15)) are estimated in the Methods S2 control stream but are not reported -- see covariatesDataExcluded.",
      source_name        = "CHN"
    ),
    DIS_HAE = list(
      description        = "Hereditary angioedema patient indicator: 1 = patient with HAE, 0 = healthy volunteer.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer)",
      notes              = "Time-fixed per subject. Acts on the EMBEDDED PK model only (multiplicative effect on CL). Every subject in the ER analysis population is by construction a patient with HAE (DIS_HAE = 1); the covariate is retained here so the embedded PK model reproduces the companion PopPK model exactly and so healthy-volunteer PK can still be simulated.",
      source_name        = "PAT"
    ),
    CREAT = list(
      description        = "Baseline serum creatinine.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Acts on the EMBEDDED PK model only (power effect on CL centred at 0.75 mg/dL). US-convention units, not SI umol/L (Section 2.4 reference subject).",
      source_name        = "BLCREAT"
    ),
    ALT = list(
      description        = "Baseline serum alanine aminotransferase activity.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Acts on the EMBEDDED PK model only (power effect on CL centred at 25 U/L).",
      source_name        = "BLALT"
    ),
    TBILI = list(
      description        = "Baseline total serum bilirubin.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Acts on the EMBEDDED PK model only (power effect on CL centred at 8 umol/L). SI units (Section 2.4 reference subject).",
      source_name        = "BLBILI"
    ),
    ON_TREATMENT = list(
      description        = "On-treatment indicator: 1 = the subject is on study after the first garadacimab or placebo administration, 0 = the untreated observational run-in period.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (untreated run-in period)",
      notes              = "TIME-VARYING within a subject: it switches from 0 to 1 at the first administration. This reproduces the source's `PLAEFF = 1; IF (TAFD .GT. 0) PLAEFF = EXP(THETA(3))` switch (Methods S2 $PK block), where TAFD is time after first dose. The effect applies to placebo recipients as well as garadacimab recipients, so it is the study/placebo effect rather than a drug effect; the drug effect enters separately through concentration. Encoding it as an explicit covariate rather than via rxode2's `tafd()` keeps the untreated run-in period simulatable (rxode2 returns NA for `tafd()` before the first dose), which matters because the run-in hazard is the denominator of the relative risk the paper reports.",
      source_name        = "TAFD > 0 (derived)"
    )
  )

  # Covariates the source screened on the ER parameters but for which no
  # point estimate is reported anywhere on disk. Documentation only --
  # checkModelConventions() does not require these to be referenced in
  # model(). See the vignette's Assumptions and deviations section.
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Baseline age (years), centred at 41 years.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "ESTIMATED in the source's final ER model on both the baseline hazard (THETA(8)) and EC50 (THETA(9)), entering as THETA*(BLAGE - 41) (Methods S2 $PK block). Neither coefficient is reported: Table S8 tabulates only the five structural ER parameters, and Figure 3 renders conditional predictions of the monthly attack rate at fixed covariate levels rather than the coefficients themselves. Reference value 41 years (Table S10)."
    ),
    SEXF = list(
      description        = "Female sex indicator: 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female; the source's reference subject per Table S10)",
      notes              = "ESTIMATED in the source's final ER model as a MALE effect on both the baseline hazard (THETA(16)) and EC50 (THETA(17)); the source column SEX codes 2 = male (Methods S2 $PK block). Neither coefficient is reported. Note the reference-category inversion relative to the canonical SEXF encoding: the source's reference subject is female (Table S10, Figure 3 caption), so a canonical SEXF term would be `theta * (1 - SEXF)`."
    ),
    HAERATE_BL = list(
      description        = "Baseline HAE attack rate during the study run-in period, centred at 0.61 attacks per week.",
      units              = "attacks per week",
      type               = "continuous",
      reference_category = NULL,
      notes              = "SCREENED but NOT retained: THETA(10) (on the baseline hazard) and THETA(11) (on EC50) are both 0 FIX in the Methods S2 control stream. Observed range 0.13-3.11 attacks per week, median 0.61 (Results Section 3.1). Note the window difference from the canonical register entry, which is normalised per 4 weeks; the source centres on a per-week rate."
    ),
    HAE_PRIOR_LTP = list(
      description        = "Indicator that the patient had used HAE long-term prophylaxis or on-demand treatment prior to study entry.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no prior HAE treatment; the source's reference subject per Table S10)",
      notes              = "ESTIMATED in the source's final ER model on both the baseline hazard (THETA(18)) and EC50 (THETA(19)) (Methods S2 $PK block, source column HAETRT). Neither coefficient is reported. No canonical covariate column is proposed because the coefficient cannot be populated."
    ),
    HAE_SUBTYPE_FXII = list(
      description        = "HAE subtype indicator: 1 = HAE with normal C1 inhibitor and factor XII mutation (HAE-FXII), 0 = other subtype.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HAE-C1INH-Type1/Type2; the source's reference subject per Table S10)",
      notes              = "ESTIMATED in the source's final ER model on both the baseline hazard (THETA(20)) and EC50 (THETA(21)); the source column HAETYPE codes 4 = HAE-FXII (Methods S2 $PK block). Neither coefficient is reported. The paper notes HAE-FXII and HAE-PLG subtypes had slightly lower model-predicted efficacy than HAE-C1INH-Type1/2 but with wide prediction intervals reflecting the limited number of such patients (Results Section 3.5)."
    ),
    CONMED_SCREENED = list(
      description        = "Concomitant and rescue medication indicators screened on the ER parameters: analgesics, anti-inflammatories, antihistamines, antibacterials, C1 esterase inhibitor, conestat alfa, icatibant acetate, and other rescue medicines.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (medication not used)",
      notes              = "SCREENED but NOT retained: all sixteen coefficients THETA(22) through THETA(37) are 0 FIX in the Methods S2 control stream, so none of these medications modifies the baseline hazard or EC50 in the final ER model. Recorded here to preserve the provenance of the source's covariate screen."
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "garadacimab", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "garadacimab", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "garadacimab", units = "mg", specimen = "plasma", verified = TRUE),
    cumhaz      = list(analyte = "Cumulative hazard of an HAE attack", units = NA_character_, specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 177L,
    n_studies      = 3L,
    age_range      = "12-73 years (median 41)",
    weight_range   = "43.3-153 kg (median 79.2)",
    sex_female_pct = NULL,
    race_ethnicity = NULL,
    disease_state  = "Patients with hereditary angioedema of all subtypes (HAE-C1INH-Type1, HAE-C1INH-Type2, and HAE-nC1INH including HAE-FXII and HAE-PLG), randomized to placebo or receiving at least one dose of garadacimab.",
    dose_range     = "Phase II subcutaneous 75, 200, or 600 mg once monthly (28 +/- 2 days); phase III subcutaneous 200 mg once monthly (30 +/- 4 days) after a loading dose of two 200 mg subcutaneous injections; plus placebo.",
    regions        = "multinational",
    biomarkers     = "Repeated time to investigator-confirmed HAE attack. Baseline attack rate during the run-in period ranged 0.13-3.11 attacks per week (median 0.61).",
    notes          = "ER analysis dataset: 177 unique patients with HAE pooled from the phase II study (NCT03712228) and the two phase III studies (pivotal VANGUARD NCT04656418 and open-label extension NCT04739059); the two phase I healthy-volunteer studies are excluded by the control stream (IGNORE STUDY.EQ.1001 and STUDY.EQ.1003). One patient from the pivotal phase III study was excluded from the final ER model because of overlapping comorbidities that may have contributed to attacks or whose symptoms may have been mistaken for HAE attacks (Section 2.3; IGNORE ID==167). Reference subject for the reported structural estimates (Table S10 and the Figure 3 caption): 41 years of age, body weight 79 kg, baseline monthly attack rate 2.9, female sex, non-Chinese and non-Japanese, HAE-C1INH-Type1/2, and no prior HAE treatment."
  )

  ini({
    # ------------------------------------------------------------------
    # PK structural parameters and covariate effects -- identical to the
    # companion PopPK model, Table S4 FINAL column. The Methods S2 ER
    # control stream re-implements the same two-compartment system in
    # $DES and drives the hazard from A(2)/(V2I/1000). See
    # modellib('Garcia_2025_garadacimab') for the full annotation,
    # including the Table S4 mixed-scale warning (continuous covariates
    # are tabulated as power EXPONENTS, categorical covariates as
    # back-transformed FOLD CHANGES).
    # ------------------------------------------------------------------
    lcl <- log(0.00664); label("Clearance at reference covariates (CL, L/h)")                    # Table S4 final model CL = 0.00664 L/h
    lvc <- log(2.37);    label("Central volume of distribution at reference weight (V2, L)")     # Table S4 final model V2 = 2.37 L
    lq  <- log(0.00685); label("Intercompartmental clearance at reference weight (Q, L/h)")      # Table S4 final model Q = 0.00685 L/h
    lvp <- log(1.41);    label("Peripheral volume of distribution at reference weight (V3, L)")  # Table S4 final model V3 = 1.41 L
    lka <- log(0.00824); label("First-order subcutaneous absorption rate constant (ka, 1/h)")    # Table S4 final model ka = 0.00824 1/h

    logitfdepot <- log(0.387 / (1 - 0.387)); label("Logit of absolute subcutaneous bioavailability (unitless; F1 = 0.387)")  # Table S4 final model F1 = 0.387

    e_wt_cl <- 1.16;  label("Power exponent of body weight on CL and Q, centred at 70 kg (unitless)")  # Table S4 final model "Weight effect on CL and Q" = 1.16
    e_wt_vc <- 0.843; label("Power exponent of body weight on Vc and Vp, centred at 70 kg (unitless)") # Table S4 final model "Weight effect on V2 and V3" = 0.843

    e_race_japanese_cl <- log(1.27); label("Log fold change in CL for Japanese vs non-Japanese subjects")       # Table S4 final model 1.27-fold
    e_race_chinese_cl  <- log(1.02); label("Log fold change in CL for Chinese vs non-Chinese subjects")         # Table S4 final model 1.02-fold
    e_dis_hae_cl       <- log(1.05); label("Log fold change in CL for patients with HAE vs healthy volunteers") # Table S4 final model 1.05-fold

    e_creat_cl <- -0.0343; label("Power exponent of baseline serum creatinine on CL, centred at 0.75 mg/dL (unitless)") # Table S4 final model -0.0343
    e_alt_cl   <- -0.0773; label("Power exponent of baseline ALT on CL, centred at 25 U/L (unitless)")                  # Table S4 final model -0.0773
    e_tbili_cl <- -0.136;  label("Power exponent of baseline total bilirubin on CL, centred at 8 umol/L (unitless)")    # Table S4 final model -0.136

    # Table S5 final model, in lower-triangle order: IIV-CL = 0.175 (CV% 43.8);
    # V2-CL covariance = 0.263 (Corr 0.743); IIV-V2 = 0.717 (CV% 102);
    # F1-CL covariance = 0.154 (Corr 0.614); F1-V2 covariance = 0.174
    # (Corr 0.342); IIV-F1 = 0.359 (logit scale). Comments must sit OUTSIDE the
    # c() -- rxode2 fails to parse a trailing comment inside a multi-line omega
    # c() block.
    etalcl + etalvc + etalogitfdepot ~ c(
      0.175,
      0.263, 0.717,
      0.154, 0.174, 0.359
    )

    # ------------------------------------------------------------------
    # ER (HAE attack hazard) structural parameters -- Table S8 FINAL
    # column. Table S8 reports these BACK-TRANSFORMED onto the natural
    # scale, which is confirmed against the Methods S2 control stream's
    # initial estimates: BHAZTV = EXP(THETA(1)) with THETA(1) = -5.41
    # gives 0.00447 versus the tabulated base-model 0.00446, and
    # EC50TV = EXP(THETA(2)) with THETA(2) = 6.09 gives 440 versus the
    # tabulated base-model 308 (same scale, re-estimated).
    #
    # Cross-check on the baseline hazard: 0.00440 /h * 24 h/day * 28 days
    # = 2.96 attacks per 4 weeks, matching the Figure 3 reference
    # subject's stated "baseline monthly attack rate of 2.9".
    # ------------------------------------------------------------------
    lhaz_base <- log(0.00440); label("Baseline (run-in) HAE attack hazard (h0, 1/h)")  # Table S8 final model baseline HAE attack hazard = 0.00440 /h (95% CI 0.00381-0.00510)

    # Constant inhibitory on-treatment (study/placebo) effect. Table S8
    # reports the MULTIPLICATIVE hazard factor 0.728 = exp(THETA(3));
    # the canonical `iplac` is the inhibited FRACTION, so it is entered
    # here as 1 - 0.728 and applied in model() as (1 - iplac * ON_TREATMENT).
    iplac <- 0.272; label("Constant inhibitory on-treatment (placebo/study) effect on the baseline HAE attack hazard (fraction)")  # Table S8 final model placebo effect = 0.728 multiplicative; 1 - 0.728 = 0.272

    lec50 <- log(303); label("Garadacimab concentration producing half-maximal inhibition of the HAE attack hazard (EC50, ng/mL)")  # Table S8 final model EC50 = 303 ng/mL (95% CI 160-533)

    # Imax and the Hill coefficient were FIXED, not estimated. Table S8
    # reports "Imax 1 (fixed)" and "Hill coefficient 1.00 (fixed)", and
    # the main text states the final ER model "had a fixed Hill
    # coefficient of 1, a fixed Imax value of 1" (Results Section 3.5).
    # The Methods S2 control stream implements Imax on the logit scale as
    # EMAX = 1/(1+EXP(-THETA(4))) with THETA(4) = 7 FIX, i.e. 0.9990889 --
    # the control-stream value is used here because an Imax of exactly 1
    # would drive the hazard to exactly zero at high concentration, which
    # the source's own likelihood (LAMBDA + 1e-6) guards against.
    logitimax <- fixed(7);      label("Logit of the maximum fractional inhibition of the HAE attack hazard (Imax = 0.999089, unitless)")  # Methods S2 $THETA(4) = 7 FIX -> Imax = 1/(1+exp(-7)) = 0.999089; Table S8 rounds this to "1 (fixed)"
    lhill     <- fixed(log(1)); label("Hill coefficient of the inhibitory hazard sigmoid (unitless)")                                      # Methods S2 $THETA(5) = 0 FIX -> HILL = exp(0) = 1; Table S8 "Hill coefficient 1.00 (fixed)"

    # ------------------------------------------------------------------
    # ER interindividual variability -- Table S9 FINAL column, a full 2x2
    # OMEGA BLOCK on the baseline hazard and EC50, in that order
    # (Methods S2: $OMEGA BLOCK(2) over ETA(1) BHAZ, ETA(2) EC50).
    # Variances are on the log scale. The IIV on EC50 is extremely large
    # (variance 4.95, CV% 1180), so individual EC50 spans orders of
    # magnitude; this is what makes the POPULATION-MEAN relative risk much
    # higher than the typical-value relative risk at any given exposure,
    # and it is why the paper's therapeutic thresholds are defined on
    # population-mean predictions.
    #
    # Note: 0.0445/sqrt(0.206*4.95) = 0.0441 versus the tabulated
    # correlation 0.0407 (and the base-model row round-trips to 0.049
    # versus a tabulated 0.00756). The covariance point estimates are used
    # here as printed; the correlation column appears to be inconsistently
    # rounded. The discrepancy is immaterial -- the correlation is
    # negligible on either reading -- and the resulting matrix is
    # positive-definite (eigenvalues 4.9504, 0.2056).
    # ------------------------------------------------------------------
    # Table S9 final model, in lower-triangle order:
    #   IIV-baseline HAE attack hazard = 0.206  (CV% 47.9; 95% CI 0.142-0.289)
    #   baseline-hazard-EC50 covariance= 0.0445
    #   IIV-EC50                       = 4.95   (CV% 1.18e3; 95% CI 3.20-7.37)
    etalhaz_base + etalec50 ~ c(
      0.206,
      0.0445, 4.95
    )

    # The source fits this model with the RTTE event-density likelihood
    # (Methods S2 $ERROR: Y = 2*LAMBDA - 2*DV*LOG(LAMBDA)), so there is no
    # observation-error model to translate. This placeholder additive
    # residual is attached to the survivor-function output so the nlmixr2
    # likelihood machinery accepts the model for forward simulation. It is
    # NOT from the source -- see the vignette's Assumptions and deviations.
    addSd <- 0.001; label("Placeholder additive residual error on the survivor-function output `sur` (unitless); not from the source")
  })

  model({
    # --- embedded PK (identical to modellib('Garcia_2025_garadacimab')) ---
    cl <- exp(lcl + etalcl) *
      (WT / 70)^e_wt_cl *
      exp(e_race_japanese_cl * RACE_JAPANESE +
            e_race_chinese_cl * RACE_CHINESE +
            e_dis_hae_cl * DIS_HAE) *
      (CREAT / 0.75)^e_creat_cl *
      (ALT / 25)^e_alt_cl *
      (TBILI / 8)^e_tbili_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q  <- exp(lq) * (WT / 70)^e_wt_cl
    vp <- exp(lvp) * (WT / 70)^e_wt_vc
    ka <- exp(lka)

    fdepot <- expit(logitfdepot + etalogitfdepot)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    Cc <- central / vc * 1000

    # --- HAE attack hazard ---------------------------------------------
    haz_base <- exp(lhaz_base + etalhaz_base)
    ec50     <- exp(lec50 + etalec50)
    imax     <- expit(logitimax)
    hill     <- exp(lhill)

    # Numerical guard only; see the companion model file. Harmless here
    # because hill is fixed at 1, but kept so the expression stays safe if
    # a user re-estimates the Hill coefficient.
    cer <- max(Cc, 0)

    # Fraction of the baseline hazard remaining after the garadacimab
    # effect (Methods S2: f_drug = log(1 - Imax*Cp^gamma/(EC50^gamma + Cp^gamma))).
    drug_effect <- 1 - imax * cer^hill / (ec50^hill + cer^hill)

    # Constant on-treatment (study/placebo) effect: 1 during the untreated
    # run-in period, 1 - iplac = 0.728 thereafter.
    plac_effect <- 1 - iplac * ON_TREATMENT

    hazard <- haz_base * plac_effect * drug_effect

    # Relative risk of an HAE attack versus the subject's own untreated
    # run-in period, which is the quantity the paper reports in Figure 5,
    # Figure S8, and Table S11 (there as a population MEAN over the IIV,
    # not the typical-value trajectory computed here).
    rr <- plac_effect * drug_effect

    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0

    # Probability of remaining attack-free since the start of the
    # integration window.
    sur <- exp(-cumhaz)

    sur ~ add(addSd)
  })
}
