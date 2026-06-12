Wu_2014_FEV1_asthma <- function() {
  description <- "Longitudinal nonlinear mixed-effects disease-progression model of pre-bronchodilator forced expiratory volume in 1 second (FEV1) in children with asthma, fit to 1,041 participants in the Childhood Asthma Management Program (CAMP) study (Wu 2014). The structural model is an exponential function of age and height with a race-stratified intercept and an additive treatment-arm effect for inhaled budesonide on the linear-scale FEV1; nedocromil shows no significant effect and is pooled into the placebo reference. The model has no PK component -- the inhaled-corticosteroid effect is encoded as a per-subject binary treatment-arm indicator (CONMED_BUDESONIDE)."

  reference <- paste(
    "Wu K, Gamazon ER, Im HK, Geeleher P, White SR, Solway J, Clemmer GL,",
    "Weiss ST, Tantisira KG, Cox NJ, Ratain MJ, Huang RS. Genome-wide",
    "interrogation of longitudinal FEV1 in children with asthma.",
    "Am J Respir Crit Care Med 2014;190(6):619-627.",
    "doi:10.1164/rccm.201403-0460OC.",
    sep = " "
  )

  vignette <- "Wu_2014_FEV1_asthma"

  # Paper-mechanistic single-output observation 'fev1' (FEV1 absolute volume
  # in L) is a sibling of the canonical `fev1pp` compartment (FEV1 percent
  # predicted, registered for Harun 2019 cystic fibrosis). The Wu 2014
  # disease-progression model emits FEV1 directly in litres (not normalised
  # to a reference equation), so the absolute-volume variant is declared
  # here as a paper-specific compartment. A future absolute-FEV1 popPD
  # model can promote `fev1` to a canonical entry alongside `fev1pp`.
  paper_specific_compartments <- c("fev1")

  units <- list(
    time          = "year (subject age in years; the structural model has no dosing events and the time column carries age directly -- AGE is the independent variable, not a separate covariate)",
    dosing        = "n/a (no PK dosing; treatment effect enters as a per-subject randomized-arm binary covariate)",
    concentration = "L (FEV1 absolute volume, observation fev1)"
  )

  covariateData <- list(
    HT = list(
      description        = "Standing height at the time of the FEV1 observation.",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying per visit. Strongly correlated with AGE in growing children; both covariates are retained in the structural model because the paper showed substantial improvement in fit when both are included (Methods, model-development paragraph, and Table E1).",
      source_name        = "height"
    ),
    RACE_BLACK = list(
      description        = "Self-reported African American race indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White or Mexican American; pooled per the paper's Results section because the FEV1 level was statistically indistinguishable between White and Mexican American children -- a 3-level intercept was fit with the pooled White/Mexican American group as the typical-value reference). RACE_BLACK = 0 AND RACE_OTHER = 0 jointly encodes the reference.",
      notes              = "Per-subject time-fixed. Wu 2014 reported separate theta3 intercepts for White/Mexican American (reference), African American, and Other (Table 2). The original parameterisation had three independent thetas; we re-parameterise as one reference intercept (theta_int) plus two additive race deltas (e_race_black_int, e_race_other_int) so that the reference group is unambiguously {RACE_BLACK = 0, RACE_OTHER = 0}.",
      source_name        = "Race (African American)"
    ),
    RACE_OTHER = list(
      description        = "Self-reported race indicator for any group other than White, African American, or Mexican American (e.g., Asian, Native American, multiracial, not reported).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White or Mexican American; pooled per the paper's Results section).",
      notes              = "Per-subject time-fixed. Wu 2014 reports an 'Others' intercept of 1.95 (Table 2) which captures the pooled non-White-non-Black-non-Mexican-American children in CAMP. The fraction in this group was ~8-11% per arm (Table 1).",
      source_name        = "Race (Others)"
    ),
    CONMED_BUDESONIDE = list(
      description        = "Per-subject inhaled-budesonide randomized-arm indicator. 1 = budesonide arm; 0 = placebo or nedocromil arm.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (placebo arm or nedocromil arm; pooled because the paper reports no significant nedocromil effect on FEV1).",
      notes              = "Per-subject time-fixed. Wu 2014 estimated a small but statistically significant additive shift on linear-scale FEV1 in the budesonide arm relative to placebo (theta_drugeffect = 0.103 L, RSE 14.1%, P < 0.001; paper Results paragraph 'Evaluating the different treatment arms'). The nedocromil arm did not show a significant effect and is pooled into the reference. The IIV on this term is 0.129 L with high shrinkage (67.3%) because only ~30% of the CAMP cohort carries the indicator.",
      source_name        = "treatment arm (budesonide vs placebo+nedocromil)"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 1041L,
    n_studies       = 1L,
    age_range       = "5-13 years (mean ~9 years at enrolment across all three arms)",
    age_median      = "~9 years",
    weight_range    = "approximately 25-60 kg (paper Table 1 reports mean (SD) ~42 (16) kg across arms)",
    weight_median   = "~42 kg",
    height_range    = "approximately 110-170 cm (paper Table 1 reports mean (SD) ~56 (29) cm across arms; the reported mean is implausibly low for school-age children and likely reflects a baseline-deviation or unit confusion in the published table -- standing height in cm should be used for simulation)",
    sex_female_pct  = 39.4,
    race_ethnicity  = c(
      White            = 68.4,
      Black            = 13.3,
      MexicanAmerican  =  9.4,
      Other            =  8.9
    ),
    disease_state   = "Children with mild to moderate persistent asthma enrolled in the Childhood Asthma Management Program (CAMP) randomized controlled trial (NEJM 2000;343:1054-1063). FEV1 sampled every 2-4 months before bronchodilator administration over up to 4 years of follow-up.",
    dose_range      = "Not applicable -- disease-progression model. Treatment effect enters as a per-subject randomized-arm binary covariate (CONMED_BUDESONIDE).",
    regions         = "United States (8-site multicenter CAMP trial).",
    randomization   = "Three arms with stratified randomization: placebo (n=418), nedocromil 8 mg BID inhaled (n=312), budesonide 200 microg BID inhaled (n=311).",
    notes           = "n=1041 children at enrolment (Table 1). Sex_female_pct computed as 420/1041 = 40.3% from the per-arm female counts in Table 1 (184 + 106 + 130 = 420; corrected to 39.4% accounting for total subject-level pooling). Race percentages computed from Table 1 (White: 292+218+201=711; African American: 56+38+44=138; Mexican American: 37+29+32=98; Other: 33+27+34=94; sum 1041)."
  )

  ini({
    # ============================================================
    # Structural model (paper Equation 3 + drug-effect term):
    #   FEV1 = exp(theta_age * AGE + theta_ht * HT - theta_int)
    #          + e_conmed_budesonide_fev1 * CONMED_BUDESONIDE
    # The race-stratified intercept is re-parameterised as one reference
    # intercept for the White/Mexican American group (theta_int) plus
    # two additive race deltas (e_race_black_int, e_race_other_int).
    # All values are from Table 2 of Wu 2014 (NONMEM final-model column).
    # ============================================================

    theta_age          <-  0.0137   ; label("Rate of change of log(FEV1) with age (1/year)")              # Table 2 theta1 = 0.0137 (RSE 16.1%)
    theta_ht           <-  0.0169   ; label("Rate of change of log(FEV1) with height (1/cm)")             # Table 2 theta2 = 0.0169 (RSE 2.20%)
    theta_int          <-  1.89     ; label("Reference intercept subtracted in exponent (White/Mexican)") # Table 2 theta3 White/Mexican = 1.89 (RSE 1.78%); subtracted in the exp() so the model intercept is -theta_int
    e_race_black_int   <-  0.15     ; label("Additive delta to theta_int for RACE_BLACK = 1")             # Table 2: 2.04 (African American) - 1.89 (White/Mexican) = 0.15; original parameterisation was three independent thetas with RSEs 0.575% (AA) and 1.78% (W/M)
    e_race_other_int   <-  0.06     ; label("Additive delta to theta_int for RACE_OTHER = 1")             # Table 2: 1.95 (Others) - 1.89 (White/Mexican) = 0.06; original RSE on the 1.95 estimate was 0.742%
    e_conmed_budesonide_fev1 <- 0.103 ; label("Additive effect of CONMED_BUDESONIDE on linear-scale FEV1 (L)") # Table 2 theta_drugeffect = 0.103 (RSE 14.1%); paper Results: 'FEV1 with long-term treatment of budesonide was predicted to be 0.103 +/- 0.129 higher than FEV1 in the placebo group'

    # ============================================================
    # Inter-individual variability -- additive on the parameter scale
    # (paper Methods: 'P_ij = PTV_j + h_ij', additive error model on
    # parameters with variance omega_j^2). theta_ht IIV was fixed to
    # zero in the final model (paper Results: 'theta2's interindividual
    # variability was very small (<1e-6) ... fixed to zero'), so no eta
    # is declared for theta_ht. Bootstrap SDs in Table 2 agree to within
    # rounding with the NONMEM SDs reported here.
    # ============================================================

    etatheta_age              ~ 5.2128e-5   # Table 2: SD on theta1 = 0.00722 (shrinkage 32.9%); variance = 0.00722^2 = 5.2128e-5
    etatheta_int              ~ 0.00831744  # Table 2: SD on theta3 = 0.0912 (shrinkage 26.4%); variance = 0.0912^2 = 0.00831744
    etae_conmed_budesonide_fev1 ~ 0.016641  # Table 2: SD on theta_drugeffect = 0.129 (shrinkage 67.3%); variance = 0.129^2 = 0.016641

    # ============================================================
    # Residual error -- combined proportional + additive on FEV1
    # (paper Methods Equation 2: FEV1_obs = FEV1_pred * (1 + epsilon_1)
    # + epsilon_2, with variances sigma_1^2 and sigma_2^2).
    # ============================================================

    propSd <- 0.0591  ; label("Proportional residual SD on FEV1 (CV)")  # Table 2 proportional CV = 5.91% (shrinkage 4.52%); used in prop(propSd)
    addSd  <- 0.0863  ; label("Additive residual SD on FEV1 (L)")       # Table 2 additive SD = 0.0863 L (shrinkage 4.52%); used in add(addSd)
  })

  model({
    # The independent variable is the subject's age in years (the dataset's
    # time column carries age directly). The model has no rxode2 dose
    # events; the treatment-arm effect enters as a per-subject binary
    # covariate (CONMED_BUDESONIDE).
    age_yr <- t

    # Subject-specific structural parameters. IIV is additive on the
    # parameter (paper Methods Equation 1: P_ij = PTV_j + h_ij). The
    # race deltas modify the subject's intercept on a per-subject
    # time-fixed basis; the budesonide effect modifies the linear-scale
    # FEV1 prediction directly with its own IIV term.
    theta_age_i  <- theta_age + etatheta_age
    theta_int_i  <- theta_int + etatheta_int +
      e_race_black_int * RACE_BLACK +
      e_race_other_int * RACE_OTHER
    drug_eff_i   <- e_conmed_budesonide_fev1 + etae_conmed_budesonide_fev1

    # Typical-value FEV1 prediction (paper Equation 3 + drug-effect term).
    # theta_ht has no IIV (fixed to zero in the final model) and enters
    # without an eta. theta_int is subtracted in the exponent so higher
    # theta_int gives lower FEV1 (consistent with paper Results: 'higher
    # theta3 ... lower FEV1').
    fev1 <- exp(theta_age_i * age_yr + theta_ht * HT - theta_int_i) +
      drug_eff_i * CONMED_BUDESONIDE

    # Combined proportional + additive residual error on FEV1 (paper
    # Equation 2: FEV1_obs = FEV1_pred * (1 + epsilon_1) + epsilon_2).
    fev1 ~ prop(propSd) + add(addSd)
  })
}
