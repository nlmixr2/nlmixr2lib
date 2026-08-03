CarlssonPetri_2018_semaglutide <- function() {
  description <- "One-compartment population PK model for once-weekly subcutaneous semaglutide (GLP-1 receptor agonist) in adults with type 2 diabetes, pooled across five SUSTAIN phase III trials (Carlsson Petri 2018)."
  reference   <- "Carlsson Petri KC, Ingwersen SH, Flint A, Zacho J, Overgaard RV. Semaglutide s.c. once-weekly in type 2 diabetes: a population pharmacokinetic analysis. Diabetes Therapy. 2018;9(4):1533-1547. doi:10.1007/s13300-018-0458-5"
  vignette    <- "CarlssonPetri_2018_semaglutide"
  units       <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot   = list(analyte = "semaglutide", units = "nmol", specimen = "administration site", verified = FALSE),
    central = list(analyte = "semaglutide", units = "nmol", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL/F with reference weight 85 kg (Table S3 exponent 0.774; Methods 'reference subject profile ... body weight of 85 kg'). Time-fixed at baseline. No covariate effect on V/F (Discussion: 'Covariate effects were not included for V/F').",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female)",
      notes              = "Reference subject profile is female (Methods). Male-vs-female CL/F ratio 1.04 (Table S3) applied as e_sexm_cl^(1 - SEXF): multiplier 1.04 when male (SEXF=0) and 1 when female (SEXF=1).",
      source_name        = "SEX (male=1)"
    ),
    AGE = list(
      description        = "Baseline subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Encoded categorically per the paper as <65 (reference), 65-74, and >=75 years (Methods). The two indicators AGE_65_74 = (AGE >= 65 & AGE < 75) and AGE_GE_75 = (AGE >= 75) are derived inline in model() from the continuous AGE column; effects on CL/F are 0.988 for AGE_65_74 and 0.961 for AGE_GE_75 (Table S3).",
      source_name        = "AGE"
    ),
    RACE_BLACK = list(
      description        = "1 = Black or African American, 0 = other race",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White, with the small 'American Indian or Alaska Native' [n=2] and 'Unknown' [n=41] subgroups merged into the White reference per Methods 'Groups that contained < 20 subjects were merged with the largest covariate group')",
      notes              = "Black-vs-White CL/F ratio 0.974 (Table S3). Reference includes merged AIAN + Unknown subgroups.",
      source_name        = "RACE"
    ),
    RACE_ASIAN = list(
      description        = "1 = Asian, 0 = other race",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White, with merged AIAN + Unknown subgroups per RACE_BLACK notes)",
      notes              = "Asian-vs-White CL/F ratio 0.989 (Table S3).",
      source_name        = "RACE"
    ),
    RACE_HISPANIC = list(
      description        = "1 = Hispanic or Latino, 0 = non-Hispanic",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Non-Hispanic or Non-Latino)",
      notes              = "The source paper treats Hispanic/Latino status as an ethnicity covariate orthogonal to race (Methods) rather than as a race indicator. The canonical RACE_HISPANIC (1 = Hispanic/Latino, 0 = non-Hispanic) is reused here because the indicator semantics are identical; the ethnicity-vs-race distinction affects population classification rather than the covariate encoding. Hispanic-vs-non-Hispanic CL/F ratio 1.06 (Table S3). SUSTAIN-Japan subjects had no ethnicity reported (Table S1) and are treated as non-Hispanic reference per the paper's analysis.",
      source_name        = "ETHNIC"
    ),
    INJSITE_ARM = list(
      description        = "1 = subject's dominant SC injection site is upper arm, 0 = abdomen (reference)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (abdomen)",
      notes              = "Per-subject indicator: 'the most frequently used injection site for an individual patient was used as the covariate value' (Methods). Upper-arm-vs-abdomen CL/F ratio 1.08 (Table S3).",
      source_name        = "INJSITE"
    ),
    INJSITE_THIGH = list(
      description        = "1 = subject's dominant SC injection site is thigh, 0 = abdomen (reference)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (abdomen)",
      notes              = "Per-subject indicator: 'the most frequently used injection site for an individual patient was used as the covariate value' (Methods). Thigh-vs-abdomen CL/F ratio 1.04 (Table S3).",
      source_name        = "INJSITE"
    ),
    RENALIMP_MILD = list(
      description        = "1 = mild renal impairment (eGFR 60-89 mL/min/1.73 m^2), 0 = normal or non-mild category",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal renal function; eGFR >= 90 mL/min/1.73 m^2)",
      notes              = "eGFR-based classification per Methods. Mild-vs-normal CL/F ratio 0.948 (Table S3). Paired with RENALIMP_MOD and RENALIMP_SEV.",
      source_name        = "RENAL"
    ),
    RENALIMP_MOD = list(
      description        = "1 = moderate renal impairment (eGFR 30-59 mL/min/1.73 m^2), 0 = other category",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal renal function; eGFR >= 90 mL/min/1.73 m^2)",
      notes              = "eGFR-based classification per Methods. Moderate-vs-normal CL/F ratio 0.955 (Table S3).",
      source_name        = "RENAL"
    ),
    RENALIMP_SEV = list(
      description        = "1 = severe renal impairment (eGFR < 30 mL/min/1.73 m^2), 0 = other category",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal renal function; eGFR >= 90 mL/min/1.73 m^2)",
      notes              = "eGFR-based classification per Methods. Severe-vs-normal CL/F ratio 0.920 (Table S3). Discussion notes the severe stratum came only from SUSTAIN 6 (n=33) and might be confounded by trial effect; a dedicated renal-impairment clinical pharmacology trial (Marbury 2017) found no clinically relevant effect.",
      source_name        = "RENAL"
    ),
    DOSE_SEMAGLUTIDE_MG = list(
      description        = "Per-subject maintenance semaglutide dose",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Values 0.5 or 1.0 mg once weekly (Methods). The paper tested a binary indicator for 0.5-vs-1.0 mg maintenance dose on CL/F to assess dose proportionality; the effect ratio is 1.00 (Table S3), confirming dose proportionality across 0.5 and 1.0 mg. The binary indicator DOSE_LOW_MAINT = (DOSE_SEMAGLUTIDE_MG < 0.75) is derived inline in model() and the CL/F multiplier is fixed at 1.00.",
      source_name        = "DOSE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1612L,                                                              # Carlsson Petri 2018 Results 'Demographics'
    n_studies      = 5L,                                                                 # SUSTAIN 1, 2, 3, 6, and SUSTAIN-Japan (Methods 'Data Sources')
    age_range      = "20-86 years (Japanese subjects >=20; global >=18)",                # Carlsson Petri 2018 Results
    age_mean       = "57 years (SD 10.6)",                                               # Table 2
    weight_range   = "39.7-198.3 kg",                                                    # Table 2
    weight_mean    = "86.2 kg (SD 22.5)",                                                # Table 2
    bmi_range      = "16.3-72.8 kg/m^2",                                                 # Table 2
    bmi_mean       = "31.1 kg/m^2 (SD 7.1)",                                             # Table 2
    sex_female_pct = 42.5,                                                               # Table 2 (685 / 1612)
    race_ethnicity = c(White = 52.0, Asian = 40.8, Black_or_African_American = 4.5,
                       American_Indian_or_Alaska_Native = 0.1, Unknown = 2.5),          # Table 2
    ethnicity_hispanic_pct = 15.0,                                                       # Table 2 (241 / 1612)
    disease_state  = "Adults with type 2 diabetes mellitus (HbA1c >= 7.0% at entry)",   # Methods
    diabetes_duration_mean = "8.1 years (SD 6.6; range 0-48.9)",                         # Table 2
    hba1c_mean     = "8.2% (SD 1; range 5.9-13.1)",                                      # Table 2
    dose_range     = "0.5 or 1.0 mg semaglutide once weekly SC after 4 weeks of 0.25 mg then 4 weeks of 0.5 mg escalation (Methods).",
    regions        = "Multinational (SUSTAIN 1-3 and 6 global; SUSTAIN-Japan is Japan-only).",
    trials         = c("NCT02054897", "NCT01930188", "NCT01885208",
                       "NCT01720446", "NCT02207374"),
    renal_function = "Normal (61.8%), mild (33.1%), moderate (3%), and severe (2.0%) renal impairment; severe only in SUSTAIN 6 (Table 2).",
    n_observations = 6781L,                                                              # Results 'Demographics' (mean 4.2 samples per subject)
    notes          = "Reference subject profile for the full-model typical values: non-Hispanic or non-Latino white female <65 years, body weight 85 kg (pre-specified as approximate median), normal renal function, dosed in abdomen with 1.0 mg semaglutide once weekly (Methods). Table 2 gives global demographics; Table S1 (ESM) gives per-trial breakdowns. American Indian or Alaska Native (n=2) and Unknown (n=41) subjects were merged with the White group for the covariate analysis (Table 2 footnote a)."
  )

  ini({
    # Structural parameters at reference covariates (Table S3, full population PK model).
    # Reference subject: non-Hispanic/Latino white female, <65 y, WT = 85 kg, normal renal
    # function, abdomen injection, 1.0 mg once-weekly maintenance dose.
    lka <- fixed(log(0.0286)); label("SC absorption rate constant (1/h)")                             # Table S3 (ka = 0.0286 h^-1, FIXED per Methods 'Semaglutide Absorption rate constant ... set to 0.0286 h^-1 based on data from clinical pharmacology trials with richly sampled PK profiles')
    lcl <- log(0.0478);        label("Apparent clearance at reference covariates (L/h)")             # Table S3 (CL/F = 0.0478 L/h, 95% CI 0.0468-0.0488, RSE 1.06%)
    lvc <- log(12.2);          label("Apparent volume of distribution (L)")                              # Table S3 (V/F = 12.2 L, 95% CI 12.1-12.4, RSE 0.487%)

    # Body weight power effect on CL/F, reference 85 kg (Methods 'reference subject profile ...
    # body weight of 85 kg').
    e_wt_cl <- 0.774;          label("Body weight power exponent on CL/F (unitless)")                    # Table S3 (0.774, 95% CI 0.724-0.823, RSE 3.27%)

    # Categorical covariate effect ratios on CL/F.  The paper's full-model equation applies each
    # ratio as theta^INDICATOR (Methods 'Exponents used for categorical covariate relations are
    # indicator variables assigned the value 1 for the actual category and else 0'); the ratio
    # itself is what Table S3 reports.
    e_sexm_cl        <- 1.04;  label("Male-vs-female CL/F ratio (theta^(1 - SEXF))")                     # Table S3 (Sex - male = 1.04, 95% CI 1.02-1.06)
    e_age_65_74_cl   <- 0.988; label("Age 65-74 vs <65 CL/F ratio (theta^AGE_65_74)")                    # Table S3 (Age 65-74 y = 0.988, 95% CI 0.966-1.01)
    e_age_ge_75_cl   <- 0.961; label("Age >=75 vs <65 CL/F ratio (theta^AGE_GE_75)")                     # Table S3 (Age >74 y = 0.961, 95% CI 0.916-1.01)
    e_dose_low_cl    <- fixed(1.00); label("0.5 mg vs 1.0 mg maintenance dose CL/F ratio (theta^DOSE_LOW_MAINT)")  # Table S3 (Maintenance dose 0.5 mg = 1.00, 95% CI 0.984-1.02; encoded as fixed(1.00) to preserve fidelity to the published estimate and confirm dose proportionality)
    e_race_black_cl  <- 0.974; label("Black-vs-White CL/F ratio (theta^RACE_BLACK)")                     # Table S3 (Race - black = 0.974, 95% CI 0.912-1.04)
    e_race_asian_cl  <- 0.989; label("Asian-vs-White CL/F ratio (theta^RACE_ASIAN)")                     # Table S3 (Race - Asian = 0.989, 95% CI 0.965-1.01)
    e_hisp_cl        <- 1.06;  label("Hispanic/Latino-vs-non-Hispanic CL/F ratio (theta^RACE_HISPANIC)") # Table S3 (Ethnicity - Hispanic or Latino = 1.06, 95% CI 1.03-1.10)
    e_site_thigh_cl  <- 1.04;  label("Thigh-vs-abdomen CL/F ratio (theta^INJSITE_THIGH)")                # Table S3 (Injection site - thigh = 1.04, 95% CI 0.993-1.08)
    e_site_arm_cl    <- 1.08;  label("Upper arm-vs-abdomen CL/F ratio (theta^INJSITE_ARM)")              # Table S3 (Injection site - upper arm = 1.08, 95% CI 1.03-1.12)
    e_renal_mild_cl  <- 0.948; label("Mild renal impairment vs normal CL/F ratio (theta^RENALIMP_MILD)") # Table S3 (Renal - Mild impairment = 0.948, 95% CI 0.930-0.965)
    e_renal_mod_cl   <- 0.955; label("Moderate renal impairment vs normal CL/F ratio (theta^RENALIMP_MOD)")  # Table S3 (Renal - Moderate impairment = 0.955, 95% CI 0.900-1.01)
    e_renal_sev_cl   <- 0.920; label("Severe renal impairment vs normal CL/F ratio (theta^RENALIMP_SEV)")    # Table S3 (Renal - Severe impairment = 0.920, 95% CI 0.846-0.995)

    # IIV as %CV (log-normal): omega^2 = log(1 + CV^2).  CL/F IIV 12.9% and V/F IIV 37.3% per
    # Table S3 (full model); paper does not report a CL-V covariance.
    etalcl ~ log(1 + 0.129^2)                                                                             # Table S3 (IIV CL/F = 12.9% CV, shrinkage 25.0%)
    etalvc ~ log(1 + 0.373^2)                                                                             # Table S3 (IIV V/F  = 37.3% CV, shrinkage 58.1% -- 'sparse sampling ... relatively high shrinkage for the volume', Discussion)

    # Proportional residual error on untransformed concentrations (Methods 'The model was
    # estimated on un-transformed concentration values, and a proportional error model was used
    # to describe the residual variability').
    propSd <- 0.238; label("Proportional residual error (fraction)")                                      # Table S3 (Proportional residual error = 23.8% CV, shrinkage 8.0%)
  })

  model({
    # Derived age-band indicators from continuous AGE (paper's categorical encoding <65 / 65-74
    # / >=75 y with <65 as reference).
    age_65_74_ind <- 0 + (AGE >= 65 & AGE < 75)
    age_ge_75_ind <- 0 + (AGE >= 75)

    # Derived low-dose indicator from the per-subject maintenance dose (paper's binary 0.5-vs-
    # 1.0 mg indicator with 1.0 mg as reference).
    dose_low_ind <- 0 + (DOSE_SEMAGLUTIDE_MG < 0.75)

    # Individual PK parameters (paper's CL/F full-model equation, applied as the product of
    # multiplicative covariate ratios raised to the corresponding indicator, plus a body-weight
    # power term; Methods 'The model was parameterized as' / Table S3).
    cl <- exp(lcl + etalcl) * (WT / 85)^e_wt_cl *
      e_sexm_cl^(1 - SEXF) *
      e_age_65_74_cl^age_65_74_ind *
      e_age_ge_75_cl^age_ge_75_ind *
      e_dose_low_cl^dose_low_ind *
      e_race_black_cl^RACE_BLACK *
      e_race_asian_cl^RACE_ASIAN *
      e_hisp_cl^RACE_HISPANIC *
      e_site_thigh_cl^INJSITE_THIGH *
      e_site_arm_cl^INJSITE_ARM *
      e_renal_mild_cl^RENALIMP_MILD *
      e_renal_mod_cl^RENALIMP_MOD *
      e_renal_sev_cl^RENALIMP_SEV
    vc <- exp(lvc + etalvc)
    ka <- exp(lka)

    kel <- cl / vc

    # One-compartment disposition with first-order SC absorption; F is folded into apparent
    # CL/F and V/F (Methods 'a one-compartment model with first-order absorption and elimination
    # ... parameterized ... in terms of ka, CL/F, and V/F').
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Concentration in nmol/L when dose is in nmol and vc is in L (semaglutide MW 4113.6 g/mol;
    # 1 mg = 1000/4113.6 = 243.1 nmol -- matches Overgaard 2019 semaglutide unit convention).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
