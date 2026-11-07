Kosinski_2023_kidneyTransplant_egfr <- function() {
  description <- paste(
    "Longitudinal disease-progression model of estimated glomerular filtration",
    "rate (eGFR, MDRD4) in adult kidney-transplant recipients from the day of",
    "transplant to 6 years (Kosinski 2023, C-Path Transplant Therapeutics",
    "Consortium pooled data, 4868 recipients). Algebraic Chapman-Richards",
    "growth with stabilizing decay (CR-SD): egfr(t) = (rmin + amp) *",
    "(1 - exp(-k_recovery * t)) + amp * (exp(-k_decline * t) - 1), where",
    "k_recovery is the early graft-recovery rate, k_decline the long-term",
    "decline rate, rmin the plateau the long-term trajectory decays towards and",
    "amp the drop from the recovery peak (rmin + amp) to that plateau. Baseline",
    "donor, graft and regimen covariates (donor type, cold ischemia time, donor",
    "age and sex, recipient age, sex and race, prior kidney transplant,",
    "induction agent and calcineurin-inhibitor use) act proportionally.",
    "Combined proportional and time-decaying additive residual error. No drug",
    "dosing; t = 0 is the transplant."
  )
  reference <- paste(
    "Kosinski L, Frey E, Klein A, O'Doherty I, Romero K, Stegall M,",
    "Helantera I, Gaber AO, Fitzsimmons WE, Aggarwal V. Longitudinal estimated glomerular",
    "filtration rate (eGFR) modeling in long-term renal function to inform",
    "clinical trial design in kidney transplantation.",
    "Clin Transl Sci. 2023;16:1680-1690. doi:10.1111/cts.13579.",
    sep = " "
  )
  vignette <- "Kosinski_2023_kidneyTransplant_egfr"
  units <- list(
    time = "year",
    dosing = "(none; disease-progression model with no drug dosing, t = 0 is the kidney transplant)",
    concentration = "mL/min/1.73 m^2 (eGFR by the 4-variable MDRD equation with the 186 coefficient; the output is a renal-function measure, not a drug concentration)"
  )

  covariateData <- list(
    DONOR_DECEASED = list(
      description = "Deceased-donor graft indicator: 1 = graft from a deceased donor, 0 = living donor. Time-fixed.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "The paper's 'Living donor' effect on ka (2.77) is coded with deceased donor as the baseline, so the model applies it as (1 + e * (1 - DONOR_DECEASED)). Table S1: 51.3% deceased, 48.7% living.",
      source_name = "Donor type (Living donor)"
    ),
    COLD_ISCHEMIA_TIME = list(
      description = "Cold ischemia time of the transplanted kidney (hours). Time-fixed.",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = "Centred on 10 h in the ka equation (Table S2: 1 + theta_CIT * (CIT - 10)). Table S1: median 10.18 h overall (range 0-85), 19.50 h for deceased and 0.95 h for living donors.",
      source_name = "CIT"
    ),
    RACE_BLACK = list(
      description = "Recipient Black race indicator: 1 = Black, 0 = non-Black. Time-fixed.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Acts on ka (-0.32), ke (+3.97) and eGFRdelta (+0.11). Table S1: 12.3% Black. The MDRD4 eGFR itself includes a race term.",
      source_name = "Black race / Black recipients"
    ),
    SEXF = list(
      description = "Recipient sex: 1 = female, 0 = male. Time-fixed.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Acts on ka only ('Female recipients' 0.27). Table S1: 38.5% female.",
      source_name = "Female recipients"
    ),
    TX_KIDNEY_PRIOR = list(
      description = "Prior kidney transplant indicator: 1 = the recipient had at least one earlier kidney transplant, 0 = first kidney transplant. Time-fixed.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Acts on ka only ('Previous kidney transplant' 0.24). Table S1: 9.8% yes.",
      source_name = "Previous kidney transplant"
    ),
    INDUCTION_ATG = list(
      description = "Induction-agent indicator: 1 = anti-thymocyte globulin induction. Member of the mutually exclusive INDUCTION_* set; all zero = basiliximab or daclizumab.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Effect on ka 0.084 (95% CI -0.046 to 0.22, not significant but retained). Table S1: 31.9%.",
      source_name = "Induction therapy: Anti-thymocyte globulin"
    ),
    INDUCTION_ALEMTUZUMAB = list(
      description = "Induction-agent indicator: 1 = alemtuzumab induction. Member of the mutually exclusive INDUCTION_* set; all zero = basiliximab or daclizumab.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Effect on ka 0.5. Table S1: 4.4%.",
      source_name = "Induction therapy: Alemtuzumab"
    ),
    INDUCTION_OTHER = list(
      description = "Induction-agent indicator: 1 = another induction agent (not basiliximab, daclizumab, anti-thymocyte globulin or alemtuzumab). Member of the mutually exclusive INDUCTION_* set.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Effect on ka 0.21. Table S1: 8.7%.",
      source_name = "Induction therapy: Other"
    ),
    INDUCTION_NONE = list(
      description = "Induction-agent indicator: 1 = no induction therapy. Member of the mutually exclusive INDUCTION_* set.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Effect on ka -0.015 (95% CI -0.14 to 0.11, not significant but retained). Table S1: 7.8%. The reference (all INDUCTION_* = 0) is basiliximab or daclizumab (39.4%).",
      source_name = "Induction therapy: None"
    ),
    CONMED_CNI = list(
      description = "Calcineurin-inhibitor-based immunosuppression indicator: 1 = regimen contains tacrolimus or cyclosporine, 0 = CNI-free regimen (e.g. mTOR inhibitor or belatacept). Time-fixed (baseline regimen).",
      units = "(binary)",
      type = "binary",
      reference_category = 1,
      notes = "The paper's 'CNI not used' effect on ke (-0.91) is applied as (1 + e * (1 - CONMED_CNI)), so CNI use (81.8% of patients, Table S1) is the reference. Table S3: 804 of the 872 CNI-free patients came from the BENEFIT / BENEFIT-EXT belatacept trials.",
      source_name = "CNI not used"
    ),
    AGE_DONOR = list(
      description = "Kidney donor age at donation (years). Time-fixed.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Centred on 44 years (Table S2) in the ke, eGFRdelta and eGFRstable equations. Table S1: median 45, range 0-84.",
      source_name = "Donor age"
    ),
    SEXF_DONOR = list(
      description = "Kidney donor sex: 1 = female donor, 0 = male donor. Time-fixed.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Acts on eGFRdelta only ('Female donors' 0.091). Table S1: 51.7% female donors.",
      source_name = "Female donors"
    ),
    AGE = list(
      description = "Recipient age at transplant (years). Time-fixed.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Centred on 51 years (Table S2) in the eGFRdelta and eGFRstable equations. Table S1: median 51, range 18-84.",
      source_name = "Recipient age"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 4868L,
    n_studies = 5L,
    age_range = "18-84 years (recipients)",
    age_median = "51 years (recipients)",
    weight_range = "not reported",
    sex_female_pct = 38.5,
    race_ethnicity = c(Black = 12.3, NonBlack = 87.6),
    disease_state = "Adult kidney-transplant recipients (end-stage kidney disease), followed from transplant to 6 years",
    dose_range = "not applicable (no drug dosing); 81.8% on calcineurin-inhibitor-based immunosuppression",
    regions = "USA (Mayo Clinic Rochester, Houston Methodist), Finland (Helsinki University Hospital) and international (BENEFIT, BENEFIT-EXT)",
    n_observations = "143,517 eGFR values after data cleaning (80% training / 20% test split)",
    notes = "C-Path Transplant Therapeutics Consortium data: three single-centre observational cohorts (2000-2020) and the BENEFIT and BENEFIT-EXT belatacept trials. Demographics in Table S1 are for the 4904 patients before exclusions. Donor type 51.3% deceased; donor age median 45 years; CIT median 10.18 h. eGFR values above 120 mL/min/1.73 m^2 were treated as above the upper limit of quantification (M3-like likelihood)."
  )

  ini({
    # Typical values: Table 2 'Beta' column; reference subject = deceased
    # donor, CIT 10 h, non-Black male recipient aged 51, first transplant,
    # basiliximab/daclizumab induction, CNI-based regimen, male donor aged 44.
    lk_recovery_egfr <- log(44.7)
    label("Early graft-recovery rate constant k_recovery (paper's ka) (1/year)") # Table 2 ka = 44.7
    lk_decline_egfr <- log(0.0034)
    label("Long-term eGFR decline rate constant k_decline (paper's ke) (1/year)") # Table 2 ke = 0.0034
    lamp_decline_egfr <- log(46.1)
    label("Amplitude of the decline from the recovery peak to the plateau (paper's eGFRdelta) (mL/min/1.73 m^2)") # Table 2 eGFRdelta = 46.1
    lrmin_egfr <- log(8.85)
    label("Plateau the long-term eGFR trajectory decays to (paper's eGFRstable) (mL/min/1.73 m^2)") # Table 2 eGFRstable = 8.85

    # Covariate effects on k_recovery: proportional, k * (1 + beta * cov)
    e_cold_ischemia_time_k_recovery_egfr <- -0.012
    label("Fractional change in k_recovery per hour of cold ischemia time above 10 h") # Table 2 ka CIT = -0.012
    e_donor_living_k_recovery_egfr <- 2.77
    label("Fractional change in k_recovery for a living-donor graft") # Table 2 ka 'Living donor' = 2.77
    e_race_black_k_recovery_egfr <- -0.32
    label("Fractional change in k_recovery for a Black recipient") # Table 2 ka 'Black race' = -0.32
    e_sexf_k_recovery_egfr <- 0.27
    label("Fractional change in k_recovery for a female recipient") # Table 2 ka 'Female recipients' = 0.27
    e_tx_kidney_prior_k_recovery_egfr <- 0.24
    label("Fractional change in k_recovery for a recipient with a previous kidney transplant") # Table 2 ka 'Previous kidney transplant' = 0.24
    e_induction_atg_k_recovery_egfr <- 0.084
    label("Fractional change in k_recovery for anti-thymocyte globulin induction") # Table 2 ka induction 'Anti-thymocyte globulin' = 0.084
    e_induction_alemtuzumab_k_recovery_egfr <- 0.5
    label("Fractional change in k_recovery for alemtuzumab induction") # Table 2 ka induction 'Alemtuzumab' = 0.5
    e_induction_other_k_recovery_egfr <- 0.21
    label("Fractional change in k_recovery for other induction agents") # Table 2 ka induction 'Other' = 0.21
    e_induction_none_k_recovery_egfr <- -0.015
    label("Fractional change in k_recovery for no induction therapy") # Table 2 ka induction 'None' = -0.015

    # Covariate effects on k_decline
    e_cni_free_k_decline_egfr <- -0.91
    label("Fractional change in k_decline for a calcineurin-inhibitor-free regimen") # Table 2 ke 'CNI not used' = -0.91
    e_age_donor_k_decline_egfr <- 0.020
    label("Fractional change in k_decline per year of donor age above 44 years") # Table 2 ke 'Donor age' = 0.020
    e_race_black_k_decline_egfr <- 3.97
    label("Fractional change in k_decline for a Black recipient") # Table 2 ke 'Black recipients' = 3.97

    # Covariate effects on amp_decline (eGFRdelta)
    e_age_donor_amp_decline_egfr <- -0.011
    label("Fractional change in amp_decline per year of donor age above 44 years") # Table 2 eGFRdelta 'Donor age' = -0.011
    e_sexf_donor_amp_decline_egfr <- 0.091
    label("Fractional change in amp_decline for a female donor") # Table 2 eGFRdelta 'Female donors' = 0.091
    e_age_amp_decline_egfr <- -0.0024
    label("Fractional change in amp_decline per year of recipient age above 51 years") # Table 2 eGFRdelta 'Recipient age' = -0.0024
    e_race_black_amp_decline_egfr <- 0.11
    label("Fractional change in amp_decline for a Black recipient") # Table 2 eGFRdelta 'Black recipients' = 0.11

    # Covariate effects on rmin (eGFRstable)
    e_age_donor_rmin_egfr <- 0.0091
    label("Fractional change in rmin per year of donor age above 44 years") # Table 2 eGFRstable 'Donor age' = 0.0091
    e_age_rmin_egfr <- 0.0058
    label("Fractional change in rmin per year of recipient age above 51 years") # Table 2 eGFRstable 'Recipient age' = 0.0058

    # IIV: Table 2 'IIV' column, read as omega^2 variances of log-normal etas
    etalk_recovery_egfr ~ 0.89 # Table 2 ka IIV = 0.89
    etalk_decline_egfr ~ 15.6 # Table 2 ke IIV = 15.6
    etalamp_decline_egfr ~ 0.081 # Table 2 eGFRdelta IIV = 0.081
    etalrmin_egfr ~ 0.78 # Table 2 eGFRstable IIV = 0.78

    # Residual error (Eq. 2): y = yhat * (1 + eps_p) + eps_a / t
    propSd <- 0.229
    label("Proportional residual error (fraction)") # Table 2 eps_prop = 0.229
    addSd <- 0.021
    label("Additive residual error coefficient; the additive SD is addSd / t (mL/min/1.73 m^2 x year)") # Table 2 eps_add = 0.021
  })

  model({
    # Table S2 covariate equations (proportional, continuous covariates
    # centred on the cohort medians: CIT 10 h, donor age 44, recipient age 51).
    k_recovery_egfr <- exp(lk_recovery_egfr + etalk_recovery_egfr) *
      (1 + e_cold_ischemia_time_k_recovery_egfr * (COLD_ISCHEMIA_TIME - 10)) *
      (1 + e_donor_living_k_recovery_egfr * (1 - DONOR_DECEASED)) *
      (1 + e_race_black_k_recovery_egfr * RACE_BLACK) *
      (1 + e_sexf_k_recovery_egfr * SEXF) *
      (1 + e_tx_kidney_prior_k_recovery_egfr * TX_KIDNEY_PRIOR) *
      (1 +
        e_induction_atg_k_recovery_egfr * INDUCTION_ATG +
        e_induction_alemtuzumab_k_recovery_egfr * INDUCTION_ALEMTUZUMAB +
        e_induction_other_k_recovery_egfr * INDUCTION_OTHER +
        e_induction_none_k_recovery_egfr * INDUCTION_NONE)

    k_decline_egfr <- exp(lk_decline_egfr + etalk_decline_egfr) *
      (1 + e_cni_free_k_decline_egfr * (1 - CONMED_CNI)) *
      (1 + e_age_donor_k_decline_egfr * (AGE_DONOR - 44)) *
      (1 + e_race_black_k_decline_egfr * RACE_BLACK)

    amp_decline_egfr <- exp(lamp_decline_egfr + etalamp_decline_egfr) *
      (1 + e_age_donor_amp_decline_egfr * (AGE_DONOR - 44)) *
      (1 + e_sexf_donor_amp_decline_egfr * SEXF_DONOR) *
      (1 + e_age_amp_decline_egfr * (AGE - 51)) *
      (1 + e_race_black_amp_decline_egfr * RACE_BLACK)

    rmin_egfr <- exp(lrmin_egfr + etalrmin_egfr) *
      (1 + e_age_donor_rmin_egfr * (AGE_DONOR - 44)) *
      (1 + e_age_rmin_egfr * (AGE - 51))

    # Recovery peak (paper's eGFRgrowth = eGFRstable + eGFRdelta)
    rmax_egfr <- rmin_egfr + amp_decline_egfr

    # Table 1 Eq. 5 (CR-SD); time = years since transplant
    egfr <- rmax_egfr * (1 - exp(-k_recovery_egfr * time)) +
      amp_decline_egfr * (exp(-k_decline_egfr * time) - 1)

    # Eq. 2: the additive term is divided by time, so it is large early after
    # transplant and undefined at time 0 (the source data had no observation
    # at or before the transplant).
    addSd_t <- addSd / time
    egfr ~ prop(propSd) + add(addSd_t)
  })
}
