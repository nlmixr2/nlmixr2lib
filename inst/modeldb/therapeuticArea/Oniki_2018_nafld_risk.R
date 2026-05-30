Oniki_2018_nafld_risk <- function() {
  description <- paste0(
    "Population BMI-driven risk-prediction model for nonalcoholic fatty ",
    "liver disease (NAFLD) prevalence in 341 elderly Japanese health-",
    "screening participants (Oniki 2018). The logit of the probability ",
    "of NAFLD is a fixed baseline floor (-5) plus a sigmoidal-Emax ",
    "function of (BMI - 17) with Hill exponent 3.43 (Eqs. 2-4 of the ",
    "source paper). Covariates female sex, HDL-C and LDL-C act on the ",
    "Emax of the logit (logit_max); the PNPLA3 rs738409 heterozygote ",
    "and homozygote indicators and HbA1c act on the (BMI50 - 17) ",
    "half-saturation offset. BMI is clamped to [17, 30] kg/m^2 before ",
    "entering the sigmoidal term (s011 dataset column BMI_A). No ",
    "between-subject random effect is estimated (ETA1 FIX 0) and the ",
    "Bernoulli-likelihood fit (METHOD=COND LAPLACE LIKELIHOOD) reports ",
    "no $SIGMA residual. The fit was performed with NONMEM 7.2.0 ",
    "(s011 control stream); MINIMIZATION SUCCESSFUL HOWEVER, PROBLEMS ",
    "OCCURRED WITH THE MINIMIZATION, with R and S matrices ",
    "algorithmically singular; OBJ 1193.995. The published 95% CIs ",
    "are bootstrap-derived (984/1000 successful runs) per the source ",
    "supplement. Companion BMI prediction model in Oniki_2018_bmi."
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
    "Companion BMI model in Oniki_2018_bmi."
  )
  vignette <- "Oniki_2018_BMI_NAFLD"
  units <- list(
    time          = "year",
    dosing        = "n/a (population NAFLD-risk prediction model; no drug input)",
    concentration = "p_nafld (probability of NAFLD, 0-1; also logit_nafld)"
  )

  covariateData <- list(
    BMI = list(
      description        = "Body mass index for this subject at the time of the risk evaluation (kg/m^2). Clamped inside model() to the interval [17, 30] kg/m^2 before entering the sigmoidal logit-of-NAFLD function (the source NONMEM dataset reports the clamp via the BMI_A derived column: BMI < 17 -> 17, BMI > 30 -> 30). Time-varying with each per-subject visit in the longitudinal screening dataset.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centring offset: 17 kg/m^2 (the lower clamp). The sigmoidal half-saturation offset (BMI50 - 17) and the BMI-driven term (BMI - 17) both use 17 as the floor. The model is undefined for BMI < 17 because (BMI - 17) < 0 would raise to the non-integer Hill power 3.43; the in-model clamp guarantees BMI - 17 >= 0. Upper clamp at 30 reflects the source dataset's small obese-stratum prevalence (1.2% per source Discussion paragraph 5).",
      source_name        = "BMI_A"
    ),
    SEXF = list(
      description        = "Sex indicator; 1 = female, 0 = male. The Oniki 2018 dataset GENDER column already encodes 0 = male / 1 = female, matching the canonical SEXF orientation without inversion.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Used as a logit-additive shift on logit_max (B = e_sexf_lmax when SEXF = 1, B = 0 when SEXF = 0); female subjects have a higher Emax of the logit-of-NAFLD curve than males at the same BMI and other covariates (Oniki 2018 Eq. 3 / Figure 2b).",
      source_name        = "GENDER"
    ),
    PNPLA3_CG = list(
      description        = "Indicator for the PNPLA3 rs738409 c.444C>G (I148M) C/G heterozygote genotype; 1 = subject carries the C/G genotype, 0 = subject does not carry the C/G genotype. Paired with PNPLA3_GG to encode the three-level rs738409 genotype with two binary indicators (C/C is the reference when both are 0). Time-fixed per subject (germline genotype).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (C/C wild-type when PNPLA3_GG is also 0)",
      notes              = "Derive from the source PNPLA3 three-level column (0 = C/C, 1 = C/G, 2 = G/G) as PNPLA3_CG = as.integer(PNPLA3 == 1). Used together with PNPLA3_GG as multiplicative scalars on the (BMI50 - 17) half-saturation offset (Oniki 2018 Eq. 4 / Figure 2c). The C/G factor (0.761) is closer to 1 than the G/G factor (0.592), consistent with an additive allele-dose effect on (BMI50 - 17).",
      source_name        = "PNPLA3"
    ),
    PNPLA3_GG = list(
      description        = "Indicator for the PNPLA3 rs738409 c.444C>G (I148M) G/G homozygote genotype; 1 = subject carries the G/G genotype, 0 = subject does not carry the G/G genotype. Paired with PNPLA3_CG to encode the three-level rs738409 genotype with two binary indicators (C/C is the reference when both are 0). Time-fixed per subject (germline genotype).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (C/C wild-type when PNPLA3_CG is also 0)",
      notes              = "Derive from the source PNPLA3 three-level column (0 = C/C, 1 = C/G, 2 = G/G) as PNPLA3_GG = as.integer(PNPLA3 == 2). See PNPLA3_CG for the joint usage and reference category.",
      source_name        = "PNPLA3"
    ),
    HBA1C = list(
      description        = "Glycated hemoglobin (HbA1c, %). Routine clinical lab measurement on whole blood; National Glycohemoglobin Standardization Program (NGSP) units.",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used via a power form (HBA1C / 5.88)^e_hba1c_bmi50 on the (BMI50 - 17) half-saturation offset (Oniki 2018 Eq. 4 / Figure 2c). Centring value 5.88% is the dataset baseline mean (pooled across DsbA-L genotypes per Table 1: weighted mean of 5.80, 5.88, 5.91 is approximately 5.83-5.88).",
      source_name        = "HbA1c"
    ),
    HDLC = list(
      description        = "Serum high-density lipoprotein cholesterol (HDL-C, mg/dL). Routine clinical lipid-panel lab measurement.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used as a linear-deviation effect on logit_max: (HDLC - 69.4) * e_hdlc_lmax (Oniki 2018 Eq. 3 / Figure 2b). Centring value 69.4 mg/dL is the dataset baseline mean (pooled across DsbA-L genotypes per Table 1: weighted mean of 69.7, 68.8, 71.2 is approximately 69.4). Lower HDLC corresponds to higher Emax of the logit (the regression coefficient e_hdlc_lmax = -0.0603 is negative).",
      source_name        = "HDL"
    ),
    LDLC = list(
      description        = "Serum low-density lipoprotein cholesterol (LDL-C, mg/dL). Routine clinical lipid-panel lab measurement.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used as a linear-deviation effect on logit_max: (LDLC - 120) * e_ldlc_lmax (Oniki 2018 Eq. 3 / Figure 2b). Centring value 120 mg/dL approximates the dataset baseline mean (pooled across DsbA-L genotypes per Table 1: weighted mean of 125.9, 125.5, 121.0 is approximately 125 mg/dL; the source NONMEM stream centres on 120 mg/dL exactly, which is the upper bound of the optimal LDLC clinical reference range).",
      source_name        = "LDL"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 341L,
    n_studies      = 1L,
    n_observations = "2015 NAFLD records across 342 subject-ids in the NONMEM dataset (Oniki 2018 s011 .lst NO. OF DATA RECS); 341 unique subjects after excluding those with habitual alcohol intake or hepatitis B/C virus positivity (Oniki 2018 Methods, Subjects and study protocol).",
    age_range      = "Elderly Japanese cohort; baseline age mean 67.7 years (SD ~5.9) pooled across DsbA-L genotypes (Oniki 2018 Table 1).",
    age_median     = NA_character_,
    weight_range   = NA_character_,
    weight_median  = NA_character_,
    sex_female_pct = 100 * (80 + 56 + 9) / (192 + 129 + 20),
    race_ethnicity = c(Asian = 100),
    disease_state  = paste0(
      "General elderly Japanese cohort; baseline NAFLD prevalence ",
      "14.1-25.0% across DsbA-L genotype strata (Oniki 2018 Table 1). ",
      "Diabetes prevalence 11.5-14.7%; hypertension prevalence ",
      "40.6-45.0%; dyslipidemia prevalence 35.0-50.0%."
    ),
    dose_range     = "n/a (no drug input; population disease-risk model)",
    regions        = "Japan (Kumamoto)",
    notes          = paste0(
      "Retrospective longitudinal observation, 5.5 +/- 1.1 years of ",
      "follow-up. NAFLD was diagnosed by hepatic ultrasonography ",
      "scanning per the Japanese practical guidelines (four criteria: ",
      "diffuse hyperechoic echotexture, increased echo vs kidneys, ",
      "vascular blurring, deep attenuation; Oniki 2018 Methods, ",
      "Measurements). The NAFLD dataset is the longitudinal record of ",
      "FLD-positive vs FLD-negative annual ultrasonography readings."
    )
  )

  ini({
    # ============================================================
    # Structural-parameter values are from the s011 NONMEM control
    # stream's FINAL PARAMETER ESTIMATE block (in agreement with
    # Eqs. 2-4 of the published Oniki 2018 paper). The fit reported
    # MINIMIZATION SUCCESSFUL HOWEVER, PROBLEMS OCCURRED WITH THE
    # MINIMIZATION; R and S matrices were algorithmically singular,
    # so the covariance step was unobtainable. The published 95% CIs
    # are bootstrap-derived (984/1000 successful runs per Supplement
    # s009). OBJ 1193.995. See vignette Assumptions and deviations
    # for the covariance-step caveat.
    # ============================================================

    # ----- Baseline floor (FIXED in source) -----
    # NONMEM: $THETA -5.00E+00 FIX ; THETA(1) BASE
    # When (BMI - 17) -> 0 the sigmoidal-Emax contribution -> 0 and
    # logit(P(NAFLD)) -> base_logit = -5, i.e. P -> expit(-5) ~ 0.007.
    base_logit         <- fixed(-5)        ; label("Baseline logit floor of P(NAFLD) for BMI <= 17 (fixed; unitless logit)")  # s011 .lst THETA(1) BASE; FIXED in $THETA

    # ----- Sigmoidal-Emax structural parameters (Eqs. 2, 4) -----
    # NONMEM variable LGT50 in s011 corresponds to (BMI50 - 17) per
    # Eq. 4 of the paper. The variable name "LGT50" in the source
    # control stream is misleading (it is not a logit-units quantity);
    # we use the clearer bmi50_m17 name below to reflect that it is
    # the (BMI50 - 17) half-saturation offset.
    bmi50_m17          <- 6.42             ; label("(BMI50 - 17) baseline at reference covariates (kg/m^2)")  # s011 .lst THETA(2) LGT50
    logit_max_base     <- 4.17             ; label("Baseline Emax of the logit-of-NAFLD sigmoid at reference covariates (unitless logit)")  # s011 .lst THETA(3) LGTMAX
    hill_bmi           <- 3.43             ; label("Hill / gamma exponent of the BMI-driven sigmoid (unitless)")  # s011 .lst THETA(4) GAMMA

    # ----- Covariate-effect parameters (Eqs. 3, 4) -----
    # PNPLA3 multiplicative factors on (BMI50 - 17): PNPLA3 C/G and G/G
    # both reduce (BMI50 - 17) (leftward-shift the sigmoid).
    e_pnpla3_cg_bmi50  <- 0.761            ; label("Multiplicative factor on (BMI50 - 17) for PNPLA3 C/G heterozygote (vs C/C reference; unitless)")  # s011 .lst THETA(5) PNPLA3=1
    e_pnpla3_gg_bmi50  <- 0.592            ; label("Multiplicative factor on (BMI50 - 17) for PNPLA3 G/G homozygote (vs C/C reference; unitless)")    # s011 .lst THETA(6) PNPLA3=2

    # Sex additive shift on logit_max.
    e_sexf_lmax        <- 1.02             ; label("Logit-additive shift on logit_max for female sex (vs male reference; unitless logit)")  # s011 .lst THETA(7) GENDER

    # HDL-C linear-deviation slope on logit_max, centred on 69.4 mg/dL.
    e_hdlc_lmax        <- -0.0603          ; label("Slope on logit_max per (HDLC - 69.4) mg/dL (unitless logit per mg/dL)")  # s011 .lst THETA(8) HDL

    # HbA1c power exponent on (BMI50 - 17), centred on 5.88%.
    e_hba1c_bmi50      <- -3.34            ; label("Power exponent for (HBA1C / 5.88) on (BMI50 - 17) (unitless)")  # s011 .lst THETA(9) HbA1c

    # LDL-C linear-deviation slope on logit_max, centred on 120 mg/dL.
    e_ldlc_lmax        <- 0.00922          ; label("Slope on logit_max per (LDLC - 120) mg/dL (unitless logit per mg/dL)")  # s011 .lst THETA(10) LDL

    # ----- Between-subject variability -----
    # NONMEM: $OMEGA 0 FIX ; ETA(1). The source paper does not estimate
    # IIV on the NAFLD-risk logit; ETA1 is constrained to zero. The eta
    # is named to pair with base_logit (NONMEM: LGT = BASE + TVLGT +
    # ETA(1); the eta is an additive shift on the baseline logit floor).
    etabase_logit ~ fixed(0)  # s011 .lst $OMEGA 0 FIX ; ETA(1)

    # ----- Residual error (placeholder) -----
    # The NONMEM fit is METHOD=COND LAPLACE LIKELIHOOD on a Bernoulli
    # outcome with no $SIGMA. nlmixr2 / rxode2 does not natively express
    # the Bernoulli likelihood through the ini()/model() observation
    # syntax in this batch, so the deterministic typical-value
    # probability p_nafld is the declared model output and a tiny
    # placeholder additive residual is attached to satisfy the
    # observation declaration. This deviation does not change the
    # typical-value prediction (the source-paper-published probability
    # of NAFLD) and is called out in the vignette Assumptions and
    # deviations section. Mirrors the Hansson_2013c_sunitinib.R and
    # Schoemaker_2018_levetiracetam.R placeholder-residual pattern.
    addSd_p_nafld <- fixed(0.001)  ; label("Placeholder additive residual SD on typical-value p_nafld; the source NAFLD likelihood is Bernoulli (no source residual)")  # not from source; see Assumptions and deviations
  })

  model({
    # ----- PNPLA3 multiplicative factor on (BMI50 - 17) -----
    # NONMEM: A = 1 if PNPLA3 == 0; THETA(5) if PNPLA3 == 1; THETA(6)
    # if PNPLA3 == 2. With the two paired binary indicators PNPLA3_CG
    # and PNPLA3_GG (mutually exclusive: at most one is 1), the
    # IF/ELSE branches collapse to the closed-form expression below.
    # When both indicators are 0 (C/C reference) pnpla3_mult = 1.
    pnpla3_mult <- (1 - PNPLA3_CG - PNPLA3_GG) * 1 +
                   PNPLA3_CG * e_pnpla3_cg_bmi50 +
                   PNPLA3_GG * e_pnpla3_gg_bmi50

    # ----- (BMI50 - 17) individual offset (Oniki 2018 Eq. 4) -----
    # NONMEM: LGT50 = THETA(2) * A * ((HbA1c / 5.88) ** THETA(9))
    bmi50_m17_i <- bmi50_m17 * pnpla3_mult * (HBA1C / 5.88)^e_hba1c_bmi50

    # ----- logit_max individual (Oniki 2018 Eq. 3) -----
    # NONMEM: LGTMAX = THETA(3) + (HDL - 69.4) * THETA(8) + B + (LDL - 120) * THETA(10)
    # where B = 0 if GENDER == 0 (male), B = THETA(7) if GENDER == 1 (female).
    logit_max_i <- logit_max_base +
                   (HDLC - 69.4) * e_hdlc_lmax +
                   SEXF * e_sexf_lmax +
                   (LDLC - 120) * e_ldlc_lmax

    # ----- Clamp BMI to [17, 30] kg/m^2 (source dataset column BMI_A) -----
    # NONMEM: BMI < 17 -> BMI_A = 17; BMI > 30 -> BMI_A = 30.
    # The (BMI - 17) factor is raised to a non-integer Hill power 3.43,
    # so the floor at 17 guarantees the base is non-negative; the
    # ceiling at 30 reflects the source-dataset's small >30 stratum and
    # avoids extrapolating the sigmoidal saturation beyond the data.
    # rxode2 has no min()/max() in model(); the equivalent floor and
    # ceiling are encoded with sign-based selectors (Hamuro 2017 pattern).
    bmi_floored <- (BMI <= 17) * 17 + (BMI >  17) * BMI
    bmi_c       <- (bmi_floored >= 30) * 30 + (bmi_floored <  30) * bmi_floored

    # ----- Sigmoidal-Emax in (BMI_A - 17) on the logit (Oniki 2018 Eq. 2) -----
    # NONMEM: TVLGT = LGTMAX * (BMI_A - 17)^GAMMA / (LGT50^GAMMA + (BMI_A - 17)^GAMMA)
    bmi_term       <- (bmi_c - 17)^hill_bmi
    bmi50_term     <- bmi50_m17_i^hill_bmi
    tvlgt          <- logit_max_i * bmi_term / (bmi50_term + bmi_term)

    # ----- Logit and probability of NAFLD (Oniki 2018 Eq. 2 / Figure 2) -----
    # NONMEM: LGT = BASE + TVLGT + ETA(1); P2 = EXP(LGT) / (1 + EXP(LGT)).
    # With etabase_logit ~ fixed(0) the ETA(1) term is identically zero.
    logit_nafld <- base_logit + tvlgt + etabase_logit
    p_nafld     <- expit(logit_nafld)

    # ----- Observation: typical-value probability of NAFLD -----
    # Deterministic typical-value declaration with placeholder additive
    # residual (see ini() comment). Downstream callers can either use
    # p_nafld directly (continuous probability) or sample Bernoulli
    # outcomes externally via rbinom(1, 1, p_nafld) on rxSolve output.
    p_nafld ~ add(addSd_p_nafld)
  })
}
