Kim_2023_DMD_6MWD_FFVL <- function() {
  description <- "Natural-history multivariate disease-progression model of the 6-minute walk distance (6MWD) and the vastus lateralis muscle fat fraction (FFVL, MR spectroscopy) versus age in ambulatory boys with Duchenne muscular dystrophy (DMD) from the ImagingDMD study (Kim 2023). 6MWD falls along a sigmoid Imax function of age (maximum fractional decrease fixed to 1) and FFVL rises along a sigmoid Emax function of age; the two endpoints are fitted jointly and linked through correlated random effects (the ages at half-maximal change of both endpoints and the 6MWD intercept). Baseline 6MWD and baseline FFVL enter as power covariates whose exponents differ between boys younger and older than 7 years at baseline, and baseline corticosteroid use delays the 6MWD decline. No drug input; the model time is the subject's age in years."
  reference <- paste(
    "Kim S, Willcocks RJ, Daniels MJ, Morales JF, Yoon DY, Triplett WT,",
    "Barnard AM, Senesac CR, Lott DJ, Forbes SC, Walter GA, Sweeney HL,",
    "Vandenborne K. (2023). Multivariate modeling of magnetic resonance",
    "biomarkers and clinical outcome measures for Duchenne muscular dystrophy",
    "clinical trials. CPT Pharmacometrics Syst Pharmacol 12(10):1437-1449.",
    "doi:10.1002/psp4.13021.",
    "Parameter estimates: Table 1 ('6MWD + FFVL' column); structure:",
    "Equations 1-5; covariate centring values and parameter transforms:",
    "Supplement 2 (Monolix model code).",
    sep = " "
  )
  vignette <- "Kim_2023_DMD_6MWD_FF"
  units <- list(
    time = "year (the subject's age; the model time IS age, not time since a dose)",
    dosing = "n/a (natural-history model, no drug input)",
    concentration = "m (walk_dist) and fraction 0-1 (ff_vl)"
  )

  covariateData <- list(
    WALK_DIST_BL = list(
      description = "Subject's own baseline 6-minute walk distance",
      units = "m",
      type = "continuous",
      reference_category = NULL,
      notes = "Measured at the (derived) baseline visit, the first visit at which 6MWD, FFVL and FFSOL were all available (Methods 'Data exclusion and formatting'). Enters as power terms whose exponent and centring value depend on the baseline-age group: (WALK_DIST_BL / 363)^beta for boys younger than 7 years at baseline and (WALK_DIST_BL / 367)^beta for older boys (Supplement 2: 'logtBSSMWTYOUNG = log(BSSMWTYOUNG/363)', 'logtBSSMWTOLD = log(BSSMWTOLD/367)'). Effects on the 6MWD intercept (old 0.63, young 0.89), on the age at half-maximal FFVL increase (young 0.2) and on the FFVL Hill coefficient (young -0.99).",
      source_name = "BSSMWT, split into BSSMWTYOUNG / BSSMWTOLD"
    ),
    FF_VL_BL = list(
      description = "Subject's own baseline vastus lateralis fat fraction (MR spectroscopy)",
      units = "fraction (0-1)",
      type = "continuous",
      reference_category = NULL,
      notes = "Same baseline visit as WALK_DIST_BL. Enters as power terms with group-specific centring: (FF_VL_BL / 0.04)^beta for boys younger than 7 years at baseline and (FF_VL_BL / 0.17)^beta for older boys (Supplement 2: 'logtBSFFVLYOUNG = log(BSFFVLYOUNG/0.04)', 'logtBSFFVLOLD = log(BSFFVLOLD/0.17)'). Effects on the 6MWD intercept (old -0.04), the age at half-maximal 6MWD decrease (old -0.05, young -0.14), the 6MWD Hill coefficient (old 0.28), the age at half-maximal FFVL increase (old -0.12, young -0.17) and the FFVL Hill coefficient (young -0.32).",
      source_name = "BSFFVL, split into BSFFVLYOUNG / BSFFVLOLD"
    ),
    AGE_BASE = list(
      description = "Age at the subject's (derived) baseline visit, time-fixed",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only to assign the baseline-age group that selects which exponent and centring value apply to WALK_DIST_BL and FF_VL_BL: young if AGE_BASE < 7 years, old otherwise (Methods: 'baseline scores of the measures to vary by young versus old baseline age (cutoff: 7 years old)'; Monolix column BAGECAT). The paper does not state the side of the boundary a boy of exactly 7.0 years falls on. Age itself was removed from the covariate list; the model time is the subject's current age.",
      source_name = "BAGE (BAGECAT for the group)"
    ),
    CONMED_STEROID = list(
      description = "Corticosteroid use at baseline (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Time-fixed steroid-use status at baseline (Table 1 abbreviations: 'STERUSE, steroid use status at baseline'). Enters as exp(beta * CONMED_STEROID) on the age at half-maximal 6MWD decrease (beta 0.09, about a 9% delay).",
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
    # Final 6MWD + FFVL multivariate model, Kim 2023 Table 1 ('6MWD + FFVL'
    # column). Supplement 2 <PARAMETER> values are the initial values of
    # that run, so only the transforms, covariate centring values and the
    # FIXED Imax are taken from it.

    # 6MWD sigmoid Imax (Kim 2023 Equation 4)
    le0_walk_dist <- log(354.02); label("6MWD extrapolated to age 0, S0_6MWD (m)") # Table 1 TV_S0,6MWD = 354.02 (RSE 0.955%)
    limax_walk_dist <- fixed(log(1)); label("Maximum fractional decrease in 6MWD, DPmax_6MWD (unitless)") # Table 1 TV_DPmax,6MWD = 1 (Fixed); Supplement 2 'Imax_pop = {value=1, method=FIXED}'
    lt50max_walk_dist <- log(11.73); label("Age at half of the maximum 6MWD decrease, DPT50_6MWD (years)") # Table 1 TV_DPT50,6MWD = 11.73 (RSE 4.35%)
    lhill_walk_dist <- log(53.7); label("Hill coefficient of the 6MWD age sigmoid, gamma_6MWD (unitless)") # Table 1 TV_gamma,6MWD = 53.7 (RSE 10.5%)

    # FFVL sigmoid Emax (Kim 2023 Equation 5)
    le0_ff_vl <- log(0.032); label("FFVL extrapolated to age 0, S0_FF (fraction)") # Table 1 TV_S0,FF = 0.032 (RSE 9.04%)
    logitemax_ff_vl <- qlogis(0.72); label("Maximum increase in FFVL, DPmax_FF, on the logit scale (fraction)") # Table 1 TV_DPmax,FF = 0.72 (RSE 4.62%); Supplement 2 'Emax = {distribution=logitNormal}' on (0, 1)
    lt50max_ff_vl <- log(11.59); label("Age at half of the maximum FFVL increase, DPT50_FF (years)") # Table 1 TV_DPT50,FF = 11.59 (RSE 2.29%)
    lhill_ff_vl <- log(7.62); label("Hill coefficient of the FFVL age sigmoid, gamma_FF (unitless)") # Table 1 TV_gamma,FF = 7.62 (RSE 6.72%)

    # Covariate effects (Kim 2023 Equation 3): power terms of baseline
    # scores centred on group medians, exp(beta) for steroid use.
    e_conmed_steroid_t50max_walk_dist <- 0.09; label("Steroid use effect on DPT50_6MWD, exp(beta) (unitless)") # Table 1 beta_STERUSE_DPT50,6MWD = 0.09 (RSE 49.8%)
    e_ff_vl_bl_old_t50max_walk_dist <- -0.05; label("Power exponent of baseline FFVL (old group) on DPT50_6MWD (unitless)") # Table 1 beta_BSFFOld_DPT50,6MWD = -0.05 (RSE 53.3%)
    e_ff_vl_bl_young_t50max_walk_dist <- -0.14; label("Power exponent of baseline FFVL (young group) on DPT50_6MWD (unitless)") # Table 1 beta_BSFFYoung_DPT50,6MWD = -0.14 (RSE 23.4%)
    e_ff_vl_bl_old_hill_walk_dist <- 0.28; label("Power exponent of baseline FFVL (old group) on gamma_6MWD (unitless)") # Table 1 beta_BSFFOld_gamma,6MWD = 0.28 (RSE 40.9%)
    e_ff_vl_bl_old_e0_walk_dist <- -0.04; label("Power exponent of baseline FFVL (old group) on S0_6MWD (unitless)") # Table 1 beta_BSFFOld_S0,6MWD = -0.04 (RSE 41.1%)
    e_walk_dist_bl_old_e0_walk_dist <- 0.63; label("Power exponent of baseline 6MWD (old group) on S0_6MWD (unitless)") # Table 1 beta_BS6MWDOld_S0,6MWD = 0.63 (RSE 8.14%)
    e_walk_dist_bl_young_e0_walk_dist <- 0.89; label("Power exponent of baseline 6MWD (young group) on S0_6MWD (unitless)") # Table 1 beta_BS6MWDYoung_S0,6MWD = 0.89 (RSE 11.1%)
    e_ff_vl_bl_old_t50max_ff_vl <- -0.12; label("Power exponent of baseline FFVL (old group) on DPT50_FF (unitless)") # Table 1 beta_BSFFOld_DPT50,FF = -0.12 (RSE 20.6%)
    e_ff_vl_bl_young_t50max_ff_vl <- -0.17; label("Power exponent of baseline FFVL (young group) on DPT50_FF (unitless)") # Table 1 beta_BSFFYoung_DPT50,FF = -0.17 (RSE 25.3%)
    e_walk_dist_bl_young_t50max_ff_vl <- 0.2; label("Power exponent of baseline 6MWD (young group) on DPT50_FF (unitless)") # Table 1 beta_BS6MWDYoung_DPT50,FF = 0.2 (RSE 129%)
    e_ff_vl_bl_young_hill_ff_vl <- -0.32; label("Power exponent of baseline FFVL (young group) on gamma_FF (unitless)") # Table 1 beta_BSFFYoung_gamma,FF = -0.32 (RSE 30.9%)
    e_walk_dist_bl_young_hill_ff_vl <- -0.99; label("Power exponent of baseline 6MWD (young group) on gamma_FF (unitless)") # Table 1 beta_BS6MWDYoung_gamma,FF = -0.99 (RSE 69.3%)

    # Between-subject variability. Table 1 prints omega as the SD of the
    # random effect (log scale, logit scale for DPmax_FF); variances are
    # omega^2. The 3 x 3 block carries the three estimated correlations.
    etalhill_walk_dist ~ 0.51^2 # Table 1 omega_gamma,6MWD = 0.51 (RSE 19.4%)
    etalogitemax_ff_vl ~ 0.61^2 # Table 1 omega_DPmax,FF = 0.61 (RSE 18.7%)
    etalhill_ff_vl ~ 0.36^2 # Table 1 omega_gamma,FF = 0.36 (RSE 13.8%)
    etale0_ff_vl ~ 0.4^2 # Table 1 omega_S0,FF = 0.4 (RSE 21.3%)
    etalt50max_walk_dist + etale0_walk_dist + etalt50max_ff_vl ~ c(
      0.17^2,
      0.57 * 0.17 * 0.06, 0.06^2,
      0.81 * 0.17 * 0.17, 0.54 * 0.06 * 0.17, 0.17^2
    ) # Table 1 omega_DPT50,6MWD = 0.17, omega_S0,6MWD = 0.06, omega_DPT50,FF = 0.17; correlations DPT50,6MWD-S0,6MWD 0.57, DPT50,6MWD-DPT50,FF 0.81, S0,6MWD-DPT50,FF 0.54

    # Residual error: proportional for 6MWD, Monolix combined1 for FFVL
    # (Supplement 2 'errorModel=proportional(b1)', 'combined1(a2, b2)').
    propSd_walk_dist <- 0.11; label("Proportional residual SD, 6MWD (fraction)") # Table 1 Prop_6MWD = 0.11 (RSE 4.77%)
    addSd_ff_vl <- 0.01; label("Additive residual SD, FFVL (fraction units)") # Table 1 Add_FF = 0.01 (RSE 26.3%)
    propSd_ff_vl <- 0.12; label("Proportional residual SD, FFVL (fraction)") # Table 1 Prop_FF = 0.12 (RSE 9.53%)
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
      (FF_VL_BL / 0.17)^(e_ff_vl_bl_old_e0_walk_dist * old) *
      (WALK_DIST_BL / 367)^(e_walk_dist_bl_old_e0_walk_dist * old) *
      (WALK_DIST_BL / 363)^(e_walk_dist_bl_young_e0_walk_dist * young)
    imax_walk_dist <- exp(limax_walk_dist)
    t50max_walk_dist <- exp(lt50max_walk_dist + etalt50max_walk_dist) *
      exp(e_conmed_steroid_t50max_walk_dist * CONMED_STEROID) *
      (FF_VL_BL / 0.17)^(e_ff_vl_bl_old_t50max_walk_dist * old) *
      (FF_VL_BL / 0.04)^(e_ff_vl_bl_young_t50max_walk_dist * young)
    hill_walk_dist <- exp(lhill_walk_dist + etalhill_walk_dist) *
      (FF_VL_BL / 0.17)^(e_ff_vl_bl_old_hill_walk_dist * old)

    e0_ff_vl <- exp(le0_ff_vl + etale0_ff_vl)
    emax_ff_vl <- expit(logitemax_ff_vl + etalogitemax_ff_vl)
    t50max_ff_vl <- exp(lt50max_ff_vl + etalt50max_ff_vl) *
      (FF_VL_BL / 0.17)^(e_ff_vl_bl_old_t50max_ff_vl * old) *
      (FF_VL_BL / 0.04)^(e_ff_vl_bl_young_t50max_ff_vl * young) *
      (WALK_DIST_BL / 363)^(e_walk_dist_bl_young_t50max_ff_vl * young)
    hill_ff_vl <- exp(lhill_ff_vl + etalhill_ff_vl) *
      (FF_VL_BL / 0.04)^(e_ff_vl_bl_young_hill_ff_vl * young) *
      (WALK_DIST_BL / 363)^(e_walk_dist_bl_young_hill_ff_vl * young)

    # Kim 2023 Equation 4: 6MWD = S0 * (1 - DPmax * Age^g / (DPT50^g + Age^g)).
    # The age fraction is written as 1 / (1 + (DPT50 / Age)^g): with the
    # steep Hill coefficient (about 54) the textbook form overflows to
    # Inf / Inf, whereas this form tends cleanly to 0 at young ages
    # (including age 0) and to 1 at old ages.
    frac_walk_dist <- 1 / (1 + (t50max_walk_dist / age)^hill_walk_dist)
    walk_dist <- e0_walk_dist * (1 - imax_walk_dist * frac_walk_dist)

    # Kim 2023 Equation 5: FF = S0 + DPmax * Age^g / (DPT50^g + Age^g).
    frac_ff_vl <- 1 / (1 + (t50max_ff_vl / age)^hill_ff_vl)
    ff_vl <- e0_ff_vl + emax_ff_vl * frac_ff_vl

    walk_dist ~ prop(propSd_walk_dist)
    ff_vl ~ add(addSd_ff_vl) + prop(propSd_ff_vl) + combined1()
  })
}
