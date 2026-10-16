Luu_2022_antiVegf_mbma <- function() {
  description <- "MBMA. Model-based meta-analysis of best-corrected visual acuity (SCORE_BCVA, ETDRS letters) in neovascular age-related macular degeneration (nAMD) under intravitreal anti-VEGF therapy with ranibizumab, aflibercept, brolucizumab or abicipar, or under sham, fitted to study-arm-level aggregate data from 22 randomised trials (55 arms, 9548 subjects) (Luu 2022). Natural disease progression is an indirect-response model in which SCORE_BCVA declines from its baseline toward a steady state gss = kin/kout; each drug stimulates kin through C/(EC50 + C) with a maximal effect Emax_ss * (1 + dEmax0 * exp(-kEmax * t)) that is largest at treatment onset. Vitreous PK is one compartment (4 mL) with first-order elimination at a fixed, drug-specific half-life; vitreous concentration is converted to nM with the molecular weight of each drug, and EC50 is drug-specific. The drug in an arm is selected by the TRT_RANIBIZUMAB / TRT_AFLIBERCEPT / TRT_BROLUCIZUMAB / TRT_ABICIPAR indicators (all 0 = sham). Arm-average baseline SCORE_BCVA (centred at 56.4 letters) is a covariate on the kout half-life (power) and on Emax_ss (on the logit scale). Between-arm variability is on the kout half-life, logit Emax_ss and dEmax0. The residual is additive in letters with separate SDs for sham arms, for treated arms in the first year and for treated arms after the first year, each divided by sqrt(N_ARM). Suitable scope is arm-mean SCORE_BCVA trajectories, NOT individual patients. Bevacizumab arms were in the fitted data but have no published half-life and are not supported."
  reference <- "Luu KT, Seal J, Green M, Winskill C, Attar M. Effect of Anti-VEGF Therapy on the Disease Progression of Neovascular Age-Related Macular Degeneration: A Systematic Review and Model-Based Meta-Analysis. J Clin Pharmacol. 2022;62(5):594-608. doi:10.1002/jcph.2002. PMCID: PMC9305109."
  vignette <- "Luu_2022_antiVegf"
  units <- list(
    time = "week",
    dosing = "mg (intravitreal injection into the vitreous compartment)",
    concentration = "Cc in mg/L (vitreous); Cc_nM in nM drives the drug effect",
    response = "SCORE_BCVA (arm-mean ETDRS letters, 0-100), output bcva"
  )

  compartmentData <- list(
    central = list(
      analyte = "anti-VEGF drug of the arm (ranibizumab, aflibercept, brolucizumab or abicipar)",
      units = "mg",
      specimen = "vitreous",
      verified = TRUE
    ),
    bcva = list(analyte = "SCORE_BCVA", units = "ETDRS letters", specimen = "not applicable", verified = TRUE)
  )

  # Between-arm (MBMA) random effects; each pairs with the typical value it
  # names (lthalf_kout, logitemaxss, ldemax0).
  paper_specific_etas <- c("eta_arm_thalf_kout", "eta_arm_emaxss", "eta_arm_demax0")

  covariateData <- list(
    SCORE_BCVA = list(
      description = "Arm-average baseline best-corrected visual acuity of the trial arm",
      units = "ETDRS letters (0-100)",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed per arm. Sets the initial condition bcva(0) and enters as a covariate on the kout half-life (power, exponent 2.57) and on Emax_ss (linear in log(SCORE_BCVA/56.4) on the logit scale, slope 2.00) (Luu 2022 Table 1 'Base_HL_kout', 'Base_Emax'; Methods: covariates 'centered on a typical value, eg, the median of the study population' and 'tested in log-space with a linear model'). The centring value is not printed; 56.4 letters is the median of the 55 arm baselines in Luu 2022 Table 2.",
      source_name = "Average Baseline Best-Corrected Visual Acuity (Luu 2022 Table 2)"
    ),
    TRT_RANIBIZUMAB = list(
      description = "Trial-arm indicator: 1 = the arm received intravitreal ranibizumab, 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "MBMA arm-level indicator. Selects the ranibizumab vitreous half-life, molecular weight and EC50. At most one TRT_* indicator may be 1; all four 0 = sham arm (no drug, sham residual SD).",
      source_name = "drug of the arm (Luu 2022 Table 2 'Arms')"
    ),
    TRT_AFLIBERCEPT = list(
      description = "Trial-arm indicator: 1 = the arm received intravitreal aflibercept, 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "MBMA arm-level indicator. Selects the aflibercept vitreous half-life, molecular weight and EC50. At most one TRT_* indicator may be 1.",
      source_name = "drug of the arm (Luu 2022 Table 2 'Arms')"
    ),
    TRT_BROLUCIZUMAB = list(
      description = "Trial-arm indicator: 1 = the arm received intravitreal brolucizumab, 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "MBMA arm-level indicator. Selects the brolucizumab vitreous half-life, molecular weight and EC50. At most one TRT_* indicator may be 1.",
      source_name = "drug of the arm (Luu 2022 Table 2 'Arms')"
    ),
    TRT_ABICIPAR = list(
      description = "Trial-arm indicator: 1 = the arm received intravitreal abicipar pegol, 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "MBMA arm-level indicator. Selects the abicipar vitreous half-life, molecular weight and EC50. At most one TRT_* indicator may be 1.",
      source_name = "drug of the arm (Luu 2022 Table 2 'Arms')"
    ),
    N_ARM = list(
      description = "Number of subjects assessed in the arm-mean SCORE_BCVA data point",
      units = "participants",
      type = "count",
      reference_category = NULL,
      notes = "MBMA weight. Luu 2022 Methods: 'The residual error was weighted by the square root of the number of subjects assessed in each data point', encoded as residual SD = sigma / sqrt(N_ARM). Affects only the residual, not the typical-value prediction. Must be strictly positive.",
      source_name = "number of subjects assessed per data point (Luu 2022 Methods)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 9548L,
    n_studies = 22L,
    n_arms = 55L,
    n_observations = "500+ arm-mean BCVA observations (Luu 2022 Abstract; Table 2 'Data Points' per trial)",
    age_range = "collected as aggregate arm demographics but not reported",
    weight_range = "not reported",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported",
    disease_state = "Neovascular (wet) age-related macular degeneration; arm-average baseline BCVA 47.0-65.8 ETDRS letters (median 56.4) across the 55 arms of Luu 2022 Table 2.",
    dose_range = "Intravitreal ranibizumab 0.3-0.5 mg, aflibercept 0.5-2 mg, brolucizumab 0.5-6 mg, abicipar 1-2 mg (single dose, Q4 x 3 loading, Q4, Q8 or Q12); bevacizumab 1.25 mg Q4 (2 arms) and sham arms.",
    regions = "Multinational (published phase 1-3 trials plus Allergan internal abicipar data).",
    notes = "Trials (Luu 2022 Table 2): ANCHOR, BAMBOO, BRAMD, C-10-083, C-12-006, CANTREAT, CATT, CEDAR, CYPRESS, EOP1003, EOP1004, EXCITE, HAGA, HARRIER, HAWK, MARINA, MORI, PIER, REACH, SEQUOIA, VIEW 1, VIEW 2. Only monotherapy arms with nominal fixed dosing; pro re nata and treat-and-extend arms excluded. Study duration 4.33-104 weeks. Arm-level aggregate data only; no individual-level data."
  )

  ini({
    # Disease progression (Luu 2022 Table 1; estimated from sham data only)
    lgss <- log(35.4); label("Steady-state SCORE_BCVA under natural disease progression, gss = kin/kout (ETDRS letters)") # Table 1 'g ss' 35.4 letters (RSE 1.46%)
    lthalf_kout <- log(21.6); label("Half-life of SCORE_BCVA deterioration, ln2/kout, at baseline 56.4 letters (week)") # Table 1 'HL_kout' 21.6 weeks (RSE 22.4%)

    # Vitreous PK (Luu 2022 Eq 2, Table 1)
    lvvit <- fixed(log(0.004)); label("Vitreous volume (L)") # Methods: 'assuming the vitreous had a well-established volume of 4 mL'
    lthalf_ranibizumab <- fixed(log(1.29)); label("Intravitreal half-life of ranibizumab (week)") # Table 1 'HL_Rani' 1.29 weeks (Fixed)
    lthalf_brolucizumab <- fixed(log(0.629)); label("Intravitreal half-life of brolucizumab (week)") # Table 1 'HL_Bro' 0.629 weeks (Fixed)
    lthalf_aflibercept <- fixed(log(0.786)); label("Intravitreal half-life of aflibercept (week)") # Table 1 'HL_Aflib' 0.786 weeks (Fixed)
    lthalf_abicipar <- fixed(log(0.649)); label("Intravitreal half-life of abicipar (week)") # Table 1 'HL_Abic' 0.649 weeks (Fixed)

    # Drug effect (Luu 2022 Eq 1, Table 1)
    logitemaxss <- logit(0.842); label("Logit of the steady-state maximal drug effect Emax_ss on kin at baseline 56.4 letters (logit of unitless fraction)") # Table 1 'E max ss' 0.842 (RSE 3.24%)
    ldemax0 <- log(5.87); label("Fold change over Emax_ss at treatment onset, dEmax0 (unitless)") # Table 1 'Delta E max 0' 5.87 (RSE 15.8%)
    lthalf_kemax <- log(0.494); label("Half-life of the decline of dEmax0, ln2/kEmax (week)") # Table 1 'HL_kE max' 0.494 weeks (RSE 32.7%)
    lec50_ranibizumab <- log(14.7); label("Ranibizumab EC50 in the vitreous (nM)") # Table 1 'Ra_EC 50' 14.7 nM (RSE 8.43%)
    lec50_aflibercept <- log(1.12); label("Aflibercept EC50 in the vitreous (nM)") # Table 1 'Af_EC 50' 1.12 nM (RSE 87.1%)
    lec50_abicipar <- log(0.0308); label("Abicipar EC50 in the vitreous (nM)") # Table 1 'Ab_EC 50' 0.0308 nM (RSE 6.04%)
    lec50_brolucizumab <- log(0.0026); label("Brolucizumab EC50 in the vitreous (nM)") # Table 1 'Br_EC 50' 0.0026 nM (RSE 6.11%)

    # Baseline SCORE_BCVA covariate effects (Luu 2022 Table 1)
    e_score_bcva_thalf_kout <- 2.57; label("Power exponent of SCORE_BCVA/56.4 on the kout half-life (unitless)") # Table 1 'Base_HL_kout' 2.57 (RSE 23.9%)
    e_score_bcva_emaxss <- 2.00; label("Slope of logit(Emax_ss) on log(SCORE_BCVA/56.4) (unitless)") # Table 1 'Base_Emax' 2.00 (RSE 6.95%)

    # Between-arm variability (Luu 2022 Table 1). The percentages are read as
    # 100 x SD of eta on the transformed scale; variance = (pct/100)^2.
    eta_arm_thalf_kout ~ 0.9801 # Table 1 'Between-arm variability HL_kout' 99% (RSE 22.9%) -> SD 0.99 on log(HL_kout)
    eta_arm_emaxss ~ 0.00874225 # Table 1 'Between-arm variability E max' 9.35% (RSE 33.6%) -> SD 0.0935 on logit(Emax_ss)
    eta_arm_demax0 ~ 1.3924 # Table 1 'Between-arm variability Delta E max' 118% (RSE 26.0%) -> SD 1.18 on log(dEmax0)

    # Residual error for an arm of one subject; divided by sqrt(N_ARM) in model()
    addSd <- 5.2; label("Additive residual SD, treated arms during the first year (ETDRS letters, for N_ARM = 1)") # Table 1 'Sigma 1' 5.2 letters (RSE 28.1%)
    addSd_late <- 4.84; label("Additive residual SD, treated arms after the first year (ETDRS letters, for N_ARM = 1)") # Table 1 'Sigma 2' 4.84 letters (RSE 15.1%)
    addSd_sham <- 6.64; label("Additive residual SD, sham arms (ETDRS letters, for N_ARM = 1)") # Table 1 'Sigma 3' 6.64 letters (RSE 5.35%)
  })
  model({
    # Drug selection. trt is 1 for a treated arm and 0 for a sham arm.
    trt <- TRT_RANIBIZUMAB + TRT_AFLIBERCEPT + TRT_BROLUCIZUMAB + TRT_ABICIPAR
    kel <- log(2) * (TRT_RANIBIZUMAB / exp(lthalf_ranibizumab) +
      TRT_AFLIBERCEPT / exp(lthalf_aflibercept) +
      TRT_BROLUCIZUMAB / exp(lthalf_brolucizumab) +
      TRT_ABICIPAR / exp(lthalf_abicipar))
    ec50 <- TRT_RANIBIZUMAB * exp(lec50_ranibizumab) +
      TRT_AFLIBERCEPT * exp(lec50_aflibercept) +
      TRT_BROLUCIZUMAB * exp(lec50_brolucizumab) +
      TRT_ABICIPAR * exp(lec50_abicipar) +
      (1 - trt)
    # Molecular weights (kDa) are not printed in Luu 2022; taken from the
    # sources it cites: Lucentis label (ranibizumab ~48 kDa), Eylea label, ref 32
    # (aflibercept 115 kDa total including glycosylation), Beovu label
    # (brolucizumab ~26 kDa) and Callanan 2018, ref 16 (abicipar 34 kDa).
    # The (1 - trt) term keeps the sham-arm denominator nonzero; no drug is
    # present there.
    mw <- TRT_RANIBIZUMAB * 48 + TRT_AFLIBERCEPT * 115 + TRT_BROLUCIZUMAB * 26 +
      TRT_ABICIPAR * 34 + (1 - trt)

    # Disease progression; kin follows from gss = kin/kout
    gss <- exp(lgss)
    lbcva <- log(SCORE_BCVA / 56.4)
    thalf_kout <- exp(lthalf_kout + e_score_bcva_thalf_kout * lbcva + eta_arm_thalf_kout)
    kout <- log(2) / thalf_kout
    kin <- gss * kout

    # Drug effect parameters. Emax_ss is logit-normal with the covariate on the
    # logit scale; dEmax0 is log-normal.
    emaxss <- expit(logitemaxss + e_score_bcva_emaxss * lbcva + eta_arm_emaxss)
    demax0 <- exp(ldemax0 + eta_arm_demax0)
    kemax <- log(2) / exp(lthalf_kemax)

    # Vitreous PK (Eq 2); mg/L to nM is x 1000 / MW(kDa)
    vvit <- exp(lvvit)
    d/dt(central) <- -kel * central
    Cc <- central / vvit
    Cc_nM <- Cc * 1000 / mw

    # Indirect response on kin (Eq 1); time is weeks from the first dose
    emax_t <- emaxss * (1 + demax0 * exp(-kemax * time))
    bcva(0) <- SCORE_BCVA
    d/dt(bcva) <- kin * (1 + emax_t * Cc_nM / (ec50 + Cc_nM)) - kout * bcva
    bcva_cfb <- bcva - SCORE_BCVA

    # Residual SD by arm type and treatment year, weighted by sqrt(N_ARM)
    late <- 0
    if (time > 52) {
      late <- 1
    }
    sdUnit <- (1 - trt) * addSd_sham + trt * ((1 - late) * addSd + late * addSd_late)
    sdArm <- sdUnit / sqrt(N_ARM)
    bcva ~ add(sdArm)
  })
}
