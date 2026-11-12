Kim_2023_DMD_6MWD_FFSOL <- function() {
  description <- "Natural-history multivariate disease-progression model of the 6-minute walk distance (6MWD) and the soleus muscle fat fraction (FFSOL, MR spectroscopy) versus age in ambulatory boys with Duchenne muscular dystrophy (DMD) from the ImagingDMD study (Kim 2023). 6MWD falls along a sigmoid Imax function of age (maximum fractional decrease fixed to 1) and FFSOL rises along a sigmoid Emax function of age whose maximum increase is bounded to (0, 2); the two endpoints are fitted jointly and linked through correlated random effects (the ages at half-maximal change of both endpoints and the 6MWD intercept). Baseline 6MWD and baseline FFSOL enter as power covariates whose exponents differ between boys younger and older than 7 years at baseline, and baseline corticosteroid use delays both the 6MWD decline and the FFSOL increase. No drug input; the model time is the subject's age in years."
  reference <- paste(
    "Kim S, Willcocks RJ, Daniels MJ, Morales JF, Yoon DY, Triplett WT,",
    "Barnard AM, Senesac CR, Lott DJ, Forbes SC, Walter GA, Sweeney HL,",
    "Vandenborne K. (2023). Multivariate modeling of magnetic resonance",
    "biomarkers and clinical outcome measures for Duchenne muscular dystrophy",
    "clinical trials. CPT Pharmacometrics Syst Pharmacol 12(10):1437-1449.",
    "doi:10.1002/psp4.13021.",
    "Parameter estimates: Table 1 ('6MWD + FFSOL' column); structure:",
    "Equations 1-5; covariate centring values and parameter transforms:",
    "Supplement 2 (Monolix model code).",
    sep = " "
  )
  vignette <- "Kim_2023_DMD_6MWD_FF"
  units <- list(
    time = "year (the subject's age; the model time IS age, not time since a dose)",
    dosing = "n/a (natural-history model, no drug input)",
    concentration = "m (walk_dist) and fraction 0-1 (ff_sol)"
  )

  covariateData <- list(
    WALK_DIST_BL = list(
      description = "Subject's own baseline 6-minute walk distance",
      units = "m",
      type = "continuous",
      reference_category = NULL,
      notes = "Measured at the (derived) baseline visit, the first visit at which 6MWD, FFVL and FFSOL were all available (Methods 'Data exclusion and formatting'). Enters as power terms whose exponent and centring value depend on the baseline-age group: (WALK_DIST_BL / 363)^beta for boys younger than 7 years at baseline and (WALK_DIST_BL / 364)^beta for older boys (Supplement 2: 'logtBSSMWTYOUNG = log(BSSMWTYOUNG/363)', 'logtBSSMWTOLD = log(BSSMWTOLD/364)'). Effects on the 6MWD intercept (old 0.68, young 0.9), the age at half-maximal 6MWD decrease (old 0.09), the FFSOL Hill coefficient (old -0.68) and, on the logit scale, the maximum FFSOL increase (young -1.32).",
      source_name = "BSSMWT, split into BSSMWTYOUNG / BSSMWTOLD"
    ),
    FF_SOL_BL = list(
      description = "Subject's own baseline soleus fat fraction (MR spectroscopy)",
      units = "fraction (0-1)",
      type = "continuous",
      reference_category = NULL,
      notes = "Same baseline visit as WALK_DIST_BL. Enters as power terms with group-specific centring: (FF_SOL_BL / 0.06)^beta for boys younger than 7 years at baseline and (FF_SOL_BL / 0.09)^beta for older boys (Supplement 2: 'logtBSFFSOLYOUNG = log(BSFFSOLYOUNG/0.06)', 'logtBSFFSOLOLD = log(BSFFSOLOLD/0.09)'). Effects on the age at half-maximal 6MWD decrease (young -0.09), the age at half-maximal FFSOL increase (old -0.07, young -0.16) and the FFSOL intercept (old 1.08, young 0.99).",
      source_name = "BSFFSOL, split into BSFFSOLYOUNG / BSFFSOLOLD"
    ),
    AGE_BASE = list(
      description = "Age at the subject's (derived) baseline visit, time-fixed",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only to assign the baseline-age group that selects which exponent and centring value apply to WALK_DIST_BL and FF_SOL_BL: young if AGE_BASE < 7 years, old otherwise (Methods: 'baseline scores of the measures to vary by young versus old baseline age (cutoff: 7 years old)'; Monolix column BAGECAT). The paper does not state the side of the boundary a boy of exactly 7.0 years falls on. Age itself was removed from the covariate list; the model time is the subject's current age.",
      source_name = "BAGE (BAGECAT for the group)"
    ),
    CONMED_STEROID = list(
      description = "Corticosteroid use at baseline (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Time-fixed steroid-use status at baseline (Table 1 abbreviations: 'STERUSE, steroid use status at baseline'). Enters as exp(beta * CONMED_STEROID) on the age at half-maximal 6MWD decrease (beta 0.13) and on the age at half-maximal FFSOL increase (beta 0.14).",
      source_name = "STERUSE"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 118L,
    n_studies = 1L,
    n_observations = "567 observations of each of 6MWD, FFVL and FFSOL (Results 'Data summary'); missing 6MWD from inability to walk handled as censored (M4).",
    age_range = "5-20 years at assessment; baseline 5-13 years",
    sex_female_pct = 0,
    disease_state = "Duchenne muscular dystrophy, ambulatory at the derived baseline visit (completed the 10-m walk/run in 45 s or less without assistance); boys still ambulatory after age 15 and those on experimental dystrophin-restoration therapies or with mutations likely to give non-zero dystrophin were excluded.",
    dose_range = "Not applicable -- natural-history model, no drug.",
    regions = "United States (ImagingDMD: University of Florida, Oregon Health & Science University, Children's Hospital of Philadelphia)",
    notes = "Figure 2f subgroups: no steroid + young 10, no steroid + old 8, steroid + young 27, steroid + old 73 (young = baseline age under 7 years). Annual MR and functional assessments over 1-9 years of follow-up. Fitted in Monolix 2021R1 (SAEM)."
  )

  ini({
    # Final 6MWD + FFSOL multivariate model, Kim 2023 Table 1 ('6MWD +
    # FFSOL' column). Supplement 2 <PARAMETER> values are the initial values
    # of that run, so only the transforms, covariate centring values and
    # the FIXED Imax are taken from it.

    # 6MWD sigmoid Imax (Kim 2023 Equation 4)
    le0_walk_dist <- log(355.4); label("6MWD extrapolated to age 0, S0_6MWD (m)") # Table 1 TV_S0,6MWD = 355.4 (RSE 1.05%)
    limax_walk_dist <- fixed(log(1)); label("Maximum fractional decrease in 6MWD, DPmax_6MWD (unitless)") # Table 1 TV_DPmax,6MWD = 1 (Fixed); Supplement 2 'Imax_pop = {value=1, method=FIXED}'
    lt50max_walk_dist <- log(10.96); label("Age at half of the maximum 6MWD decrease, DPT50_6MWD (years)") # Table 1 TV_DPT50,6MWD = 10.96 (RSE 5.31%)
    lhill_walk_dist <- log(52.92); label("Hill coefficient of the 6MWD age sigmoid, gamma_6MWD (unitless)") # Table 1 TV_gamma,6MWD = 52.92 (RSE 10.3%)

    # FFSOL sigmoid Emax (Kim 2023 Equation 5). Supplement 2 declares
    # 'Emax = {distribution=logitNormal, max=2}', so DPmax_FF is logit-normal
    # on (0, 2): logit(DPmax_FF / 2) carries the covariate and the eta.
    le0_ff_sol <- log(0.06); label("FFSOL extrapolated to age 0, S0_FF (fraction)") # Table 1 TV_S0,FF = 0.06 (RSE 3.46%)
    logitemax_ff_sol <- qlogis(0.49 / 2); label("Maximum increase in FFSOL, DPmax_FF, as logit(DPmax_FF / 2) (fraction)") # Table 1 TV_DPmax,FF = 0.49 (RSE 7.6%)
    lt50max_ff_sol <- log(12.16); label("Age at half of the maximum FFSOL increase, DPT50_FF (years)") # Table 1 TV_DPT50,FF = 12.16 (RSE 4.72%)
    lhill_ff_sol <- log(6.87); label("Hill coefficient of the FFSOL age sigmoid, gamma_FF (unitless)") # Table 1 TV_gamma,FF = 6.87 (RSE 5.25%)

    # Covariate effects (Kim 2023 Equation 3): power terms of baseline
    # scores centred on group medians, exp(beta) for steroid use; the
    # DPmax_FF effect is additive on its logit scale.
    e_conmed_steroid_t50max_walk_dist <- 0.13; label("Steroid use effect on DPT50_6MWD, exp(beta) (unitless)") # Table 1 beta_STERUSE_DPT50,6MWD = 0.13 (RSE 41.1%)
    e_ff_sol_bl_young_t50max_walk_dist <- -0.09; label("Power exponent of baseline FFSOL (young group) on DPT50_6MWD (unitless)") # Table 1 beta_BSFFYoung_DPT50,6MWD = -0.09 (RSE 61.9%)
    e_walk_dist_bl_old_t50max_walk_dist <- 0.09; label("Power exponent of baseline 6MWD (old group) on DPT50_6MWD (unitless)") # Table 1 beta_BS6MWDOld_DPT50,6MWD = 0.09 (RSE 54.4%)
    e_walk_dist_bl_old_e0_walk_dist <- 0.68; label("Power exponent of baseline 6MWD (old group) on S0_6MWD (unitless)") # Table 1 beta_BS6MWDOld_S0,6MWD = 0.68 (RSE 7.57%)
    e_walk_dist_bl_young_e0_walk_dist <- 0.9; label("Power exponent of baseline 6MWD (young group) on S0_6MWD (unitless)") # Table 1 beta_BS6MWDYoung_S0,6MWD = 0.9 (RSE 15.1%)
    e_walk_dist_bl_young_emax_ff_sol <- -1.32; label("Coefficient of log baseline 6MWD (young group) on logit(DPmax_FF / 2) (unitless)") # Table 1 beta_BS6MWDYoung_DPmax,FF = -1.32 (RSE 81.4%)
    e_conmed_steroid_t50max_ff_sol <- 0.14; label("Steroid use effect on DPT50_FF, exp(beta) (unitless)") # Table 1 beta_STERUSE_DPT50,FF = 0.14 (RSE 34.1%)
    e_ff_sol_bl_old_t50max_ff_sol <- -0.07; label("Power exponent of baseline FFSOL (old group) on DPT50_FF (unitless)") # Table 1 beta_BSFFOld_DPT50,FF = -0.07 (RSE 34.7%)
    e_ff_sol_bl_young_t50max_ff_sol <- -0.16; label("Power exponent of baseline FFSOL (young group) on DPT50_FF (unitless)") # Table 1 beta_BSFFYoung_DPT50,FF = -0.16 (RSE 29.3%)
    e_walk_dist_bl_old_hill_ff_sol <- -0.68; label("Power exponent of baseline 6MWD (old group) on gamma_FF (unitless)") # Table 1 beta_BS6MWDOld_gamma,FF = -0.68 (RSE 28.5%)
    e_ff_sol_bl_old_e0_ff_sol <- 1.08; label("Power exponent of baseline FFSOL (old group) on S0_FF (unitless)") # Table 1 beta_BSFFOld_S0,FF = 1.08 (RSE 5.29%)
    e_ff_sol_bl_young_e0_ff_sol <- 0.99; label("Power exponent of baseline FFSOL (young group) on S0_FF (unitless)") # Table 1 beta_BSFFYoung_S0,FF = 0.99 (RSE 4.98%)

    # Between-subject variability. Table 1 prints omega as the SD of the
    # random effect (log scale, logit scale for DPmax_FF); variances are
    # omega^2. The 3 x 3 block carries the three estimated correlations.
    etalhill_walk_dist ~ 0.45^2 # Table 1 omega_gamma,6MWD = 0.45 (RSE 19.7%)
    etalogitemax_ff_sol ~ 0.5^2 # Table 1 omega_DPmax,FF = 0.5 (RSE 15.1%)
    etalhill_ff_sol ~ 0.24^2 # Table 1 omega_gamma,FF = 0.24 (RSE 13.5%)
    etale0_ff_sol ~ 0.08^2 # Table 1 omega_S0,FF = 0.08 (RSE 34.8%)
    etalt50max_walk_dist + etale0_walk_dist + etalt50max_ff_sol ~ c(
      0.15^2,
      0.58 * 0.15 * 0.06, 0.06^2,
      0.9 * 0.15 * 0.13, 0.43 * 0.06 * 0.13, 0.13^2
    ) # Table 1 omega_DPT50,6MWD = 0.15, omega_S0,6MWD = 0.06, omega_DPT50,FF = 0.13; correlations DPT50,6MWD-S0,6MWD 0.58, DPT50,6MWD-DPT50,FF 0.9, S0,6MWD-DPT50,FF 0.43

    # Residual error: proportional for 6MWD, Monolix combined1 for FFSOL
    # (Supplement 2 'errorModel=proportional(b_1)', 'combined1(a_2, b_2)').
    propSd_walk_dist <- 0.12; label("Proportional residual SD, 6MWD (fraction)") # Table 1 Prop_6MWD = 0.12 (RSE 5.3%)
    addSd_ff_sol <- 0.01; label("Additive residual SD, FFSOL (fraction units)") # Table 1 Add_FF = 0.01 (RSE 21.6%)
    propSd_ff_sol <- 0.1; label("Proportional residual SD, FFSOL (fraction)") # Table 1 Prop_FF = 0.1 (RSE 12.4%)
  })
  model({
    # The model time is the subject's age in years (Monolix regressor
    # TIME = age at assessment).
    age <- t

    # Baseline-age group: young below 7 years at baseline. Each baseline
    # score carries a separate exponent and centring value per group; the
    # power term of the other group is 1.
    young <- 0
    if (AGE_BASE < 7) {
      young <- 1
    }
    old <- 1 - young

    e0_walk_dist <- exp(le0_walk_dist + etale0_walk_dist) *
      (WALK_DIST_BL / 364)^(e_walk_dist_bl_old_e0_walk_dist * old) *
      (WALK_DIST_BL / 363)^(e_walk_dist_bl_young_e0_walk_dist * young)
    imax_walk_dist <- exp(limax_walk_dist)
    t50max_walk_dist <- exp(lt50max_walk_dist + etalt50max_walk_dist) *
      exp(e_conmed_steroid_t50max_walk_dist * CONMED_STEROID) *
      (FF_SOL_BL / 0.06)^(e_ff_sol_bl_young_t50max_walk_dist * young) *
      (WALK_DIST_BL / 364)^(e_walk_dist_bl_old_t50max_walk_dist * old)
    hill_walk_dist <- exp(lhill_walk_dist + etalhill_walk_dist)

    e0_ff_sol <- exp(le0_ff_sol + etale0_ff_sol) *
      (FF_SOL_BL / 0.09)^(e_ff_sol_bl_old_e0_ff_sol * old) *
      (FF_SOL_BL / 0.06)^(e_ff_sol_bl_young_e0_ff_sol * young)
    emax_ff_sol <- 2 *
      expit(
        logitemax_ff_sol +
          e_walk_dist_bl_young_emax_ff_sol * young * log(WALK_DIST_BL / 363) +
          etalogitemax_ff_sol
      )
    t50max_ff_sol <- exp(lt50max_ff_sol + etalt50max_ff_sol) *
      exp(e_conmed_steroid_t50max_ff_sol * CONMED_STEROID) *
      (FF_SOL_BL / 0.09)^(e_ff_sol_bl_old_t50max_ff_sol * old) *
      (FF_SOL_BL / 0.06)^(e_ff_sol_bl_young_t50max_ff_sol * young)
    hill_ff_sol <- exp(lhill_ff_sol + etalhill_ff_sol) *
      (WALK_DIST_BL / 364)^(e_walk_dist_bl_old_hill_ff_sol * old)

    # Kim 2023 Equation 4: 6MWD = S0 * (1 - DPmax * Age^g / (DPT50^g + Age^g)).
    # The age fraction is written as 1 / (1 + (DPT50 / Age)^g): with the
    # steep Hill coefficient (about 53) the textbook form overflows to
    # Inf / Inf, whereas this form tends cleanly to 0 at young ages
    # (including age 0) and to 1 at old ages.
    frac_walk_dist <- 1 / (1 + (t50max_walk_dist / age)^hill_walk_dist)
    walk_dist <- e0_walk_dist * (1 - imax_walk_dist * frac_walk_dist)

    # Kim 2023 Equation 5: FF = S0 + DPmax * Age^g / (DPT50^g + Age^g).
    frac_ff_sol <- 1 / (1 + (t50max_ff_sol / age)^hill_ff_sol)
    ff_sol <- e0_ff_sol + emax_ff_sol * frac_ff_sol

    walk_dist ~ prop(propSd_walk_dist)
    ff_sol ~ add(addSd_ff_sol) + prop(propSd_ff_sol) + combined1()
  })
}
