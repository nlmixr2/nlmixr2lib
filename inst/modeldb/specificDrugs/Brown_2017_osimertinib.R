Brown_2017_osimertinib <- function() {
  description <- "Joint two-compartment population PK model for osimertinib (AZD9291) and its active metabolite AZ5104 in advanced non-small cell lung cancer (NSCLC) patients pooled with healthy volunteers (Brown 2017). First-order oral absorption into a parent (osimertinib) compartment is followed by a second compartment (AZ5104) in series; the fraction of parent eliminated as AZ5104 is fixed at 0.25 per the publication. Body weight (allometric on parent CL/F and Vc/F and on AZ5104 CL/F), serum albumin (power on parent Vc/F), healthy-volunteer disease state (linear factor on both parent and AZ5104 CL/F), and ethnicity (Chinese, Japanese, Asian-other, and non-Asian non-Caucasian linear factors on AZ5104 CL/F) were retained as significant covariates."
  reference <- "Brown K, Comisar C, Witjes H, Maringwa J, de Greef R, Vishwanathan K, Cantarini M, Cox E. Population pharmacokinetics and exposure-response of osimertinib in patients with non-small cell lung cancer. Br J Clin Pharmacol. 2017;83(6):1216-1226. doi:10.1111/bcp.13223"
  vignette <- "Brown_2017_osimertinib"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline; reported in kg).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effects on parent CL/F (exponent 0.56), parent Vc/F (exponent 0.65), and AZ5104 CL/F (exponent 0.99); reference body weight 62 kg per Brown 2017 Table 1 baseline median.",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Baseline serum albumin concentration.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on parent Vc/F (exponent 1.33); reference albumin 39 g/L per Brown 2017 Table 1 baseline median.",
      source_name        = "Albumin"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-volunteer cohort indicator (1 = Study 5 healthy volunteer, 0 = NSCLC patient pooled from AURA / AURA2).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (NSCLC patient; reference cohort is the pooled AURA + AURA2 advanced-NSCLC population).",
      notes              = "Linear effects (1 + 0.44 * DIS_HEALTHY) on parent CL/F and (1 + 1.25 * DIS_HEALTHY) on AZ5104 CL/F. The healthy-volunteer cohort is Study 5 (D5160C00005, 32 subjects). The reference NSCLC values for clearance are the published typical values 14.2 L/h (parent) and 31.5 L/h (AZ5104); HV subjects therefore have ~44 percent higher parent CL/F and ~125 percent higher AZ5104 CL/F (35 percent and 55 percent lower AUCss, respectively).",
      source_name        = "POP"
    ),
    RACE_CHINESE = list(
      description        = "Chinese-heritage race indicator (1 = Chinese, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Chinese; the paper-defined reference category for race is Caucasian).",
      notes              = "Linear additive effect (1 + 0.17 * RACE_CHINESE) on AZ5104 CL/F; ~14.5 percent decrease in AZ5104 AUCss vs Caucasian. Brown 2017 carries Chinese, Japanese, Asian-other, and non-Asian-non-Caucasian as four mutually exclusive race indicators with Caucasian as the reference; subjects with missing race are encoded as 0 for every indicator (i.e., treated as the Caucasian reference).",
      source_name        = "Ethnic Asian Chinese (Brown 2017 Table 2)"
    ),
    RACE_JAPANESE = list(
      description        = "Japanese-heritage race indicator (1 = Japanese, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Japanese; the paper-defined reference category for race is Caucasian).",
      notes              = "Linear additive effect (1 + 0.20 * RACE_JAPANESE) on AZ5104 CL/F; ~16.7 percent decrease in AZ5104 AUCss vs Caucasian.",
      source_name        = "Ethnic Asian Japanese (Brown 2017 Table 2)"
    ),
    RACE_ASIAN_OTH = list(
      description        = "Asian-other race indicator (1 = Asian heritage other than Chinese or Japanese, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian-other; the paper-defined reference category for race is Caucasian).",
      notes              = "Linear additive effect (1 + 0.21 * RACE_ASIAN_OTH) on AZ5104 CL/F; ~17.4 percent decrease in AZ5104 AUCss vs Caucasian. Maps onto Brown 2017's 'Asian (not Japanese or Chinese)' category. Dominant reference cohort is Caucasian (not Chinese).",
      source_name        = "Ethnic Asian other (Brown 2017 Table 2)"
    ),
    RACE_OTHER = list(
      description        = "Other-race indicator (1 = race category Other, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the paper-defined reference category for race is Caucasian).",
      notes              = "Linear additive effect (1 + 0.10 * RACE_OTHER) on AZ5104 CL/F; ~9.1 percent decrease in AZ5104 AUCss vs Caucasian. Maps onto Brown 2017's non-Asian non-Caucasian 'Other' category.",
      source_name        = "Ethnic non-Asian non-Caucasian (Brown 2017 Table 2)"
    )
  )

  population <- list(
    n_subjects     = 780,
    n_studies      = 3,
    age_range      = "21-89 years",
    age_median     = "61 years",
    weight_range   = "33-122 kg",
    weight_median  = "62 kg",
    sex_female_pct = 36.9,
    race_ethnicity = c(
      White            = 24.2,
      Asian_other      = 24.2,
      Japanese         = 18.8,
      Chinese          = 15.3,
      Other            = 6.4,
      Missing          = 11.2
    ),
    disease_state  = "Advanced EGFR-mutation-positive non-small cell lung cancer (NSCLC; n = 748, 95.9 percent) pooled with healthy adult volunteers (n = 32, 4.1 percent). NSCLC patients enrolled in AURA (D5160C00001, phase I/II) or AURA2 (D5160C00002, phase II); HV cohort enrolled in Study 5 (D5160C00005, phase I).",
    dose_range     = "Once-daily oral osimertinib 20-240 mg in NSCLC patients (predominantly capsule in AURA escalation/expansion, film-coated tablet in AURA extension and AURA2); single 20 mg oral dose in healthy volunteers (capsule, phase 1 tablet, or oral solution; one fed cohort). 80 mg once daily was the approved dose at the time of analysis.",
    regions        = "Multiregional (AURA / AURA2 enrolled at sites in North America, Europe, and Asia; Study 5 was a single-centre HV study).",
    n_observations = "21 930 plasma concentrations from 780 individuals included in the final popPK; sub-LLOQ concentrations were excluded. LLOQ = 0.05 nmol/L (osimertinib) and 0.0515 nmol/L (AZ5104).",
    notes          = "Demographic counts and ranges reproduced from Brown 2017 Table 1 (continuous covariates: mean, SD, median, range, n, missing; categorical covariates: counts and percentages). Race percentages do not sum to 100 because a Missing category is reported separately; subjects with missing race are encoded as 0 for every RACE_* indicator in this implementation, i.e., treated as the Caucasian reference."
  )

  ini({
    # Structural parameters - reference values are NSCLC patient typical
    # values at 62 kg / 39 g/L albumin (Brown 2017 Table 1 medians).

    lka         <- log(0.24);   label("First-order oral absorption rate constant (1/h)")                  # Brown 2017 Table 2: ka = 0.24 1/h
    lcl         <- log(14.2);   label("Apparent osimertinib clearance, CL_parent/F, in NSCLC at 62 kg (L/h)")  # Brown 2017 Table 2: CLparent/F = 14.2 L/h
    lvc         <- log(986);    label("Apparent osimertinib central volume, V_parent/F, in NSCLC at 62 kg / 39 g/L albumin (L)")  # Brown 2017 Table 2: Vparent/F = 986 L
    lcl_az5104  <- log(31.5);   label("Apparent AZ5104 clearance, CL_metabolite/F, in NSCLC Caucasian at 62 kg (L/h)")  # Brown 2017 Table 2: CLmetabolite/F = 31.5 L/h
    lvc_az5104  <- log(207);    label("Apparent AZ5104 central volume, V_metabolite/F (L)")               # Brown 2017 Table 2: Vmetabolite/F = 207 L

    # Continuous-covariate power-form effects (Brown 2017 Methods Eq.
    # P_j = theta_k * (X_ij / median(X_j))^theta_i; estimated, not fixed).
    e_wt_cl         <- 0.56;   label("Power exponent for body weight on CL_parent (unitless)")            # Brown 2017 Table 2: 0.56
    e_wt_vc         <- 0.65;   label("Power exponent for body weight on V_parent (unitless)")             # Brown 2017 Table 2: 0.65
    e_alb_vc        <- 1.33;   label("Power exponent for albumin on V_parent (unitless)")                 # Brown 2017 Table 2: 1.33
    e_wt_cl_az5104  <- 0.99;   label("Power exponent for body weight on CL_metabolite (unitless)")        # Brown 2017 Table 2: 0.99

    # Categorical-covariate linear-form effects (Brown 2017 Methods Eq.
    # P_j = theta_k * (1 + theta_i)^X_ij with X_ij in {0,1}).
    e_healthy_cl              <- 0.44;   label("Linear coefficient for HV vs NSCLC on CL_parent (unitless)")    # Brown 2017 Table 2: 0.44
    e_healthy_cl_az5104       <- 1.25;   label("Linear coefficient for HV vs NSCLC on CL_metabolite (unitless)")  # Brown 2017 Table 2: 1.25
    e_race_chinese_cl_az5104 <- 0.17;   label("Linear coefficient for Chinese vs Caucasian on CL_metabolite (unitless)")  # Brown 2017 Table 2: 0.17
    e_race_japanese_cl_az5104 <- 0.20;  label("Linear coefficient for Japanese vs Caucasian on CL_metabolite (unitless)") # Brown 2017 Table 2: 0.20
    e_race_asian_oth_cl_az5104 <- 0.21; label("Linear coefficient for Asian-other vs Caucasian on CL_metabolite (unitless)")  # Brown 2017 Table 2: 0.21
    e_race_other_cl_az5104   <- 0.10;   label("Linear coefficient for non-Asian non-Caucasian Other vs Caucasian on CL_metabolite (unitless)")  # Brown 2017 Table 2: 0.10

    # IIV. Brown 2017 Table 2 reports omega values (square root of the
    # log-scale variance; the table's IIV percent column equals omega *
    # 100 to two significant figures). Variances used by ini() are the
    # squares of the reported omega values. Off-diagonal covariance:
    # rho * omega_CL_parent * omega_CL_metabolite = 0.90 * 0.46 * 0.52
    # per Brown 2017 Table 2 footnote a (Correlation r = 0.90).
    etalcl + etalcl_az5104 ~ c(0.2116, 0.21528, 0.2704)                     # Brown 2017 Table 2: omega_CL_parent = 0.46, omega_CL_metabolite = 0.52, r = 0.90
    etalka      ~ 0.7921                                                    # Brown 2017 Table 2: omega_ka = 0.89
    etalvc      ~ 0.2704                                                    # Brown 2017 Table 2: omega_V_parent = 0.52
    etalvc_az5104 ~ 0.3844                                                  # Brown 2017 Table 2: omega_V_metabolite = 0.62

    # Residual error. Brown 2017 reports a single combined proportional
    # plus additive error model in Table 2; the same value-pair is
    # applied to both observed analytes (osimertinib, AZ5104). The
    # additive-error units in Table 2 are misprinted as 'g l-1' in the
    # published table but the LLOQ and all reported concentrations are
    # in nmol/L, so 0.105 is interpreted as nmol/L (consistent with
    # ~2 x LLOQ for osimertinib). Converted to mg/L using the published
    # molecular weights for osimertinib (499.62 g/mol) and AZ5104
    # (485.59 g/mol; N-desmethyl-osimertinib).
    propSd        <- 0.244;      label("Proportional residual error on osimertinib (fraction)")                                            # Brown 2017 Table 2: 24.4%
    addSd         <- 5.246e-5;   label("Additive residual error on osimertinib (mg/L) = 0.105 nmol/L * 499.62 g/mol / 1e6")                # Brown 2017 Table 2: 0.105 nmol/L (units interpreted from context; see vignette Errata)
    propSd_az5104 <- 0.244;      label("Proportional residual error on AZ5104 (fraction; shared with osimertinib in Brown 2017 Table 2)")  # Brown 2017 Table 2: 24.4%
    addSd_az5104  <- 5.099e-5;   label("Additive residual error on AZ5104 (mg/L) = 0.105 nmol/L * 485.59 g/mol / 1e6")                     # Brown 2017 Table 2: 0.105 nmol/L; converted via AZ5104 MW
  })

  model({
    # Constants. Osimertinib molecular weight 499.62 g/mol (PubChem CID
    # 71496458); AZ5104 (N-desmethyl-osimertinib, C27H31N7O2) molecular
    # weight 485.59 g/mol. Stoichiometric correction (mw_az5104 /
    # mw_parent) converts the parent-mass elimination flux into the
    # AZ5104 mass formation flux so the central / metabolite
    # compartments and observed Cc / Cc_az5104 are in mass-per-volume
    # units (mg in compartment, mg/L in observation). The mg/L
    # observations correspond to nmol/L in Brown 2017 via factor
    # (1e6 / molecular weight); the validation vignette converts.
    mw_parent <- 499.62
    mw_az5104 <- 485.59
    # Fraction of parent CL/F appearing as AZ5104 - fixed at 0.25 per
    # Brown 2017 Methods. The publication notes this value is arbitrary
    # because fmet, V_AZ5104/F, and CL_AZ5104/F are confounded in a
    # joint parent / metabolite model; simulated osimertinib and AZ5104
    # concentrations are insensitive to fmet provided the same value is
    # used at simulation as at fit.
    fmet <- 0.25

    # Reference covariate values - Brown 2017 Table 1 medians.
    ref_wt  <- 62
    ref_alb <- 39

    # Individual parameters.
    ka <- exp(lka + etalka)

    cl <- exp(lcl + etalcl) *
          (WT / ref_wt)^e_wt_cl *
          (1 + e_healthy_cl * DIS_HEALTHY)

    vc <- exp(lvc + etalvc) *
          (WT / ref_wt)^e_wt_vc *
          (ALB / ref_alb)^e_alb_vc

    cl_az5104 <- exp(lcl_az5104 + etalcl_az5104) *
                 (WT / ref_wt)^e_wt_cl_az5104 *
                 (1 + e_healthy_cl_az5104 * DIS_HEALTHY) *
                 (1 + e_race_chinese_cl_az5104   * RACE_CHINESE) *
                 (1 + e_race_japanese_cl_az5104  * RACE_JAPANESE) *
                 (1 + e_race_asian_oth_cl_az5104 * RACE_ASIAN_OTH) *
                 (1 + e_race_other_cl_az5104     * RACE_OTHER)

    vc_az5104 <- exp(lvc_az5104 + etalvc_az5104)

    # ODE system (amounts in mg; 1:1 molar parent -> AZ5104 stoichiometry
    # rendered in mass via mw_az5104 / mw_parent).
    d/dt(depot)          <- -ka * depot
    d/dt(central)        <-  ka * depot - (cl / vc) * central
    d/dt(central_az5104) <-  fmet * (cl / vc) * central * (mw_az5104 / mw_parent) -
                             (cl_az5104 / vc_az5104) * central_az5104

    # Observations (mg/L = ug/mL; multiply by 1e6 / molecular weight to
    # get the nmol/L scale that Brown 2017 reports).
    Cc        <- central / vc
    Cc_az5104 <- central_az5104 / vc_az5104

    Cc        ~ prop(propSd) + add(addSd)
    Cc_az5104 ~ prop(propSd_az5104) + add(addSd_az5104)
  })
}
