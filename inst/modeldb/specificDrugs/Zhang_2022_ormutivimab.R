Zhang_2022_ormutivimab <- function() {
  description <- "Time-dependent population pharmacodynamic emax model for rabies virus neutralizing antibody (RVNA) activity after rabies vaccination in healthy Chinese adults, with a categorical drug-product covariate that contrasts Ormutivimab (rHRIG, a recombinant human anti-rabies IgG1 monoclonal antibody) against plasma-derived human rabies immunoglobulin (HRIG) (Zhang 2022). Output Cc is neutralizing antibody activity in IU/mL measured by the rapid fluorescent focus inhibition test (RFFIT). The published Y1 two-compartment PK overlay for the passive-antibody component of the combined drug+vaccine groups (E = Y1 + Y2) is NOT included here because the seven structural PK constants (Ka, V1, V2, K10, K12, K21, C0) are not reported anywhere on disk; see the vignette's Assumptions and deviations section for the omitted-component audit trail."
  reference <- "Zhang J, Hao Y, Liu L, et al. Population Pharmacodynamic Analyses of Human Anti-Rabies Virus Monoclonal Antibody (Ormutivimab) in Healthy Adult Subjects. Vaccines (Basel). 2022;10(8):1218. doi:10.3390/vaccines10081218. PMID: 36016106."
  vignette <- "Zhang_2022_ormutivimab"
  units <- list(time = "day", dosing = "IU", concentration = "IU/mL")

  covariateData <- list(
    DRUG_ORMU = list(
      description        = "Indicator for Ormutivimab administration",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HRIG, plasma-derived human rabies immunoglobulin comparator; also used for placebo + vaccine subjects who received no passive antibody)",
      notes              = "Time-fixed per subject. 1 = subject received Ormutivimab (rHRIG, the recombinant human anti-rabies IgG drug). 0 = subject received plasma-derived HRIG. Modifies the typical-value emax and ET50 of the time-dependent vaccine-induced RVNA emax model. Zhang 2022 Eq. 11 (Results section 3.4): emax = Emax_HRIG + 0.143 for rHRIG; ET50 = ET50_HRIG - 3.8 for rHRIG. The paper reports only this binary contrast; placebo + vaccine subjects (no passive antibody) inherit the HRIG-arm typical emax/ET50 because their fitted-time-emax curve captures vaccine response alone and the binary covariate is structured as a shift from the HRIG baseline.",
      source_name        = "antibody type"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 240L,
    n_studies        = 1L,
    age_range        = "18 - 55 years; phase IIb mean (SD) 39.7 (8.1) years",
    age_median       = "39.7 years (mean)",
    weight_range     = "phase IIb mean (SD) 64.7 (11.2) kg",
    weight_median    = "64.7 kg (mean)",
    sex_female_pct   = 55.0,
    race_ethnicity   = c(Chinese = 100),
    disease_state    = "Healthy adults without prior rabies exposure or rabies immunization.",
    dose_range       = "Single IM injection on Day 0 of Ormutivimab 20 or 40 IU/kg, or HRIG 20 IU/kg, or placebo, all in combination with Vero-cell rabies vaccine on Days 0, 3, 7, 14, 28 (Essen regimen).",
    regions          = "China (Chinese CDC sites in Beijing-Chaoyang and Shanxi).",
    n_subjects_total = 300L,
    notes            = "Final-model parameters were estimated on the phase IIb cohort (n = 240; four arms of ~60 subjects each: placebo+vaccine, HRIG 20 IU/kg + vaccine, Ormutivimab 20 IU/kg + vaccine, Ormutivimab 40 IU/kg + vaccine; demographics per Zhang 2022 Table 2). The phase IIa cohort (n = 60; passive-antibody only, no vaccine) was used to develop the upstream two-compartment passive-antibody PK model that is NOT included in this nlmixr2lib model file (parameter values were not published; see vignette Errata). ClinicalTrials.gov NCT02559921."
  )

  ini({
    # -----------------------------------------------------------------------
    # Final-model PD parameters (Zhang 2022 Table 3, "Final Model Estimates"
    # column; bootstrap medians and 95% PIs are tabulated alongside but the
    # ini() values reproduce the point estimates verbatim). The structural
    # equation is the time-dependent emax model (Zhang 2022 Eq. 10):
    #
    #     Y2(t) = e0 + emax * t^hill / (ET50^hill + t^hill)
    #
    # with the typical-value covariate equations (Zhang 2022 Eq. 11):
    #
    #     emax = Emax_HRIG + 0       if HRIG
    #     emax = Emax_HRIG + 0.143   if rHRIG (Ormutivimab)
    #     ET50 = ET50_HRIG + 0       if HRIG
    #     ET50 = ET50_HRIG - 3.8     if rHRIG (Ormutivimab)
    #
    # IIV is log-normal on emax and hill (Zhang 2022 Eq. 1: P_i = P_TV *
    # exp(eta_i)); ET50 and e0 have no IIV in the final model. Residual
    # error is the combined proportional + additive form (Zhang 2022 Eq. 4:
    # Y_obs = Y_pred * (1 + eps1) + eps2).
    # -----------------------------------------------------------------------

    # Typical emax / ET50 (HRIG baseline) - log-transformed because emax and
    # ET50 are constrained positive; the additive HRIG -> rHRIG shift is
    # applied in linear space inside model() per the paper's equation.
    lemax <- log(3.6);    label("Typical maximum vaccine-induced RVNA effect emax (IU/mL) at the HRIG-arm baseline")  # Zhang 2022 Table 3: emax 3.6 IU/mL (RSE 15.1%)
    lET50 <- log(10.5);   label("Typical time to half-maximum vaccine-induced RVNA effect ET50 (day) at the HRIG-arm baseline")  # Zhang 2022 Table 3: ET50 10.5 day (RSE 8.4%)

    # Categorical-covariate shifts (additive on the linear-scale typical
    # value). theta1 raises emax in the rHRIG arm by 0.143 IU/mL; theta2
    # lowers ET50 in the rHRIG arm by 3.8 days. The combined values are
    #   emax(rHRIG) = 3.6 + 0.143 = 3.743 IU/mL
    #   ET50(rHRIG) = 10.5 - 3.8  = 6.7  day
    # which yield a higher and faster vaccine-induced antibody peak in the
    # Ormutivimab + vaccine arms (Zhang 2022 Discussion).
    e_drug_ormu_emax <-  0.143; label("Additive Ormutivimab-vs-HRIG shift in typical emax (IU/mL)")  # Zhang 2022 Table 3: theta1 0.143 (RSE 40.1%); Eq. 11
    e_drug_ormu_ET50 <- -3.8;   label("Additive Ormutivimab-vs-HRIG shift in typical ET50 (day)")    # Zhang 2022 Table 3: theta2 -3.8 (RSE 19.1%); Eq. 11

    # Hill / sigmoidicity exponent of the time-dependent emax model. The
    # paper does not separately partition hill by drug type, so this is a
    # shared typical value across HRIG and Ormutivimab arms.
    lhill <- log(7.66); label("Hill / sigmoidicity exponent of the time-dependent emax model (unitless)")  # Zhang 2022 Table 3: Gamma 7.66 (RSE 22.6%)

    # e0 baseline offset. NOT log-transformed because the fitted value is
    # negative (-3.19): Y2 at t = 0 equals e0, and the linear-scale emax
    # equation lets the model offset e0 below zero so that the predicted
    # RVNA curve reproduces the observed near-zero baseline before vaccine
    # antibody production starts and rises smoothly toward the asymptote
    # e0 + emax in late-time observations. The combined drug+vaccine model
    # E = Y1 + Y2 used the omitted Y1 passive-antibody component to lift
    # the early-time curve back into the physically observed positive
    # range; in the present Y2-only file e0 should be interpreted as a
    # model-internal offset rather than a physiological baseline activity.
    e0 <- -3.19; label("Y2 model offset at t = 0 (IU/mL; fitted, can be negative)")  # Zhang 2022 Table 3: e0 -3.19 (RSE 16.9%)

    # Inter-individual variability. Paper Eq. 1: P_i = P_TV * exp(eta_i),
    # so eta is on log(P) and the reported CV% maps to omega^2 via
    # omega^2 = log(CV^2 + 1).
    #   omega(emax) 9.0%   -> log(0.09^2  + 1) = 0.008068
    #   omega(hill) 56.1% -> log(0.561^2 + 1) = 0.273690
    etalEmax  ~ 0.008068  # Zhang 2022 Table 3: omega(emax)  =  9.0%; log(0.09^2 + 1)  = 0.008068
    etalhill ~ 0.273690  # Zhang 2022 Table 3: omega(Gamma) = 56.1%; log(0.561^2 + 1) = 0.273690

    # Residual error - combined proportional + additive (Zhang 2022 Table 3
    # row "Residual Error"; corresponds to Eq. 4: Y_obs = Y_pred * (1 + eps1)
    # + eps2). Reported on the natural (linear) RVNA scale (IU/mL).
    addSd  <- 0.245; label("Additive residual error on RVNA (IU/mL)")        # Zhang 2022 Table 3: sigma_additive    0.245 (RSE 3.8%)
    propSd <- 0.094; label("Proportional residual error on RVNA (fraction)") # Zhang 2022 Table 3: sigma_proportional 9.4% (RSE 53.9%)
  })

  model({
    # Drug-product typical-value shifts (Zhang 2022 Eq. 11). Applied on the
    # linear scale before IIV is multiplied in (paper Eq. 1: P_i = P_TV
    # * exp(eta_i)).
    emax_tv <- exp(lemax) + e_drug_ormu_emax * DRUG_ORMU
    ET50_tv <- exp(lET50) + e_drug_ormu_ET50 * DRUG_ORMU

    # Individual parameters - log-normal IIV on emax and hill; ET50 and e0
    # carry only the typical value (no eta in final model per Table 3).
    emax  <- emax_tv * exp(etalEmax)
    ET50  <- ET50_tv
    hill <- exp(lhill + etalhill)

    # Time-dependent vaccine-induced RVNA emax (Zhang 2022 Eq. 10). The
    # guard at t <= 0 avoids a NaN from t^hill at negative integration
    # times that some solvers visit during initialization (hill is
    # non-integer, so (-eps)^hill is undefined in R).
    if (t <= 0) {
      Cc <- e0
    } else {
      Cc <- e0 + emax * t^hill / (ET50^hill + t^hill)
    }

    # Combined residual error (Zhang 2022 Eq. 4 / Table 3).
    Cc ~ add(addSd) + prop(propSd)
  })
}
