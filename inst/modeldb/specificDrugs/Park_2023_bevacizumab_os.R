Park_2023_bevacizumab_os <- function() {
  description <- "Parametric log-logistic time-to-event exposure-response model for OVERALL SURVIVAL in adults with metastatic or recurrent non-squamous non-small-cell lung cancer receiving bevacizumab 15 mg/kg IV every 3 weeks (biosimilar SB8 or EU-sourced reference bevacizumab) with paclitaxel and carboplatin induction followed by bevacizumab maintenance (Park 2023, phase III study SB8-G31-NSCLC, n = 224). The baseline hazard is log-logistic, h0(t) = (gamma / lambda) * (t / lambda)^(gamma - 1) / (1 + (t / lambda)^gamma) with lambda in months, so the drug-free median survival equals lambda. The steady-state trough serum bevacizumab concentration CTROUGH (ug/mL) reduces the hazard by a constant fraction through an inhibitory sigmoid Emax term, h(t) = h0(t) * (1 - emax * CTROUGH^hill / (ec50^hill + CTROUGH^hill)), with the Hill coefficient fixed at 20 (a near all-or-none switch around ec50). Treatment (SB8 versus reference) was not a predictor and no other covariate was retained. There is no PK layer: CTROUGH is supplied as a time-fixed data column, the observed pre-dose concentration before cycle 7 (or cycle 5). Forward simulation exposes `hazard` (deaths per month), `sur` (survival probability) and `medtte` (the closed-form median survival time in months for the subject's CTROUGH)."
  reference <- paste(
    "Park S, Jung JA, Ju S, Lim HS.",
    "Exposure-response analysis using time-to-event data for bevacizumab",
    "biosimilar SB8 and the reference bevacizumab.",
    "Front Pharmacol. 2023;14:1288308. doi:10.3389/fphar.2023.1288308.",
    "Structural equations from the Supplementary Material NONMEM control",
    "file; parameter estimates from Table 3.",
    "Sister model from the same paper: modellib('Park_2023_bevacizumab_pfs').",
    sep = " "
  )
  vignette <- "Park_2023_bevacizumab"
  units <- list(
    time = "month",
    dosing = "n/a (no drug-dosing events; the exposure driver is supplied as the CTROUGH data covariate, in ug/mL)",
    concentration = "probability (the model output `sur` is the overall-survival probability, not a drug concentration)"
  )

  compartmentData <- list(
    cumhaz = list(
      analyte = "cumulative hazard of death",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Observed steady-state trough (pre-dose) serum bevacizumab concentration",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Required input, time-fixed per subject. TOTAL serum bevacizumab (SB8 or EU-sourced reference) MEASURED, not model-predicted: the pre-dose sample before cycle 7 (221 patients) or, failing that, before cycle 5 (3 patients), at least 100 days after the first dose of 15 mg/kg IV every 3 weeks (Park 2023 Methods). The same value is applied over the whole follow-up from randomisation, including the induction cycles before steady state was reached. Enters the hazard through the inhibitory sigmoid Emax form emax * CTROUGH^hill / (ec50^hill + CTROUGH^hill), not linearly. Park 2023 Table 1: median (Q1-Q3) 115.1 (87.8-144.6) ug/mL for SB8, 125.5 (97.0-153.6) for the reference product, 119.5 (92.0-150.5) pooled.",
      source_name = "Css,trough (CONC in the control stream)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Significant in the univariate step for OS (Park 2023 Table 2, dOFV -8.956) but not retained in the forward selection after the concentration effect was added (dOFV -5.319 against the 6.63 criterion)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Significant in the univariate step for OS (Park 2023 Table 2, dOFV -7.613) but not retained in the forward selection after the concentration effect was added (dOFV -4.415 against the 6.63 criterion)."
    ),
    AGE = list(
      description = "Age, screened as < 65 versus >= 65 years",
      units = "years",
      type = "continuous",
      notes = "Screened and not retained (Park 2023 Table 2, dOFV -0.129)."
    ),
    SEXF = list(
      description = "Sex",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained (Park 2023 Table 2, dOFV -1.727)."
    ),
    ADA_POS = list(
      description = "Overall anti-drug antibody status, by cycle 7 and by end of treatment",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained (Park 2023 Table 2, dOFV -0.164 by cycle 7 and -1.742 by end of treatment)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 224L,
    n_studies = 1L,
    age_range = "Adults; 155 (69.2%) < 65 years and 69 (30.8%) >= 65 years (Park 2023 Table 1).",
    weight_range = "Median 72.1 kg (Q1-Q3 64.2-82.0) (Park 2023 Table 1).",
    sex_female_pct = 41.5,
    race_ethnicity = "White 92.9%, non-white 7.1% (Park 2023 Table 1).",
    disease_state = "Metastatic or recurrent non-squamous non-small-cell lung cancer (93.3% adenocarcinoma; 99.6% stage IV), ECOG performance status 0 (33%) or 1 (67%).",
    dose_range = "Bevacizumab (SB8, n = 100, or EU-sourced reference bevacizumab, n = 124) 15 mg/kg IV every 3 weeks with paclitaxel 200 mg/m^2 and carboplatin AUC 6 for 4-6 induction cycles, then bevacizumab 15 mg/kg IV every 3 weeks as maintenance until progression, unacceptable toxicity, death or end of study.",
    regions = "Europe (EU 25.4%; Russia, Ukraine, Belarus, Georgia, Serbia and others) and East Asia (Korea, Thailand, Taiwan; 7.1%) (Park 2023 Table 1).",
    notes = "Subset of the phase III comparative study SB8-G31-NSCLC (EudraCT 2015-004026-34, NCT02754882; Reck 2020) restricted to patients with an available steady-state trough concentration. Overall survival was assessed from randomisation; the visual predictive check spans 24 months (Park 2023 Figure 3). Fitted in NONMEM 7.4 by FOCE-I."
  )

  ini({
    # Park 2023 Table 3, overall-survival rows. The log-logistic scale lambda
    # is in months (Figures 3 and 6 are drawn in months, and the drug-free
    # median OS of Figure 6, 12 months, equals lambda).
    llambda <- log(12.0); label("Log-logistic scale lambda of the baseline hazard (month); equals the drug-free median survival time") # Table 3 OS lambda = 12.0 (RSE 4.13%, 95% CI 11.028-12.972)
    lgamma <- log(3.67); label("Log-logistic shape gamma of the baseline hazard (unitless)") # Table 3 OS gamma = 3.67 (RSE 3.13%, 95% CI 3.445-3.895)
    lemax <- log(0.704); label("Maximum fractional reduction of the hazard by bevacizumab (fraction)") # Table 3 OS EMAX = 0.704 (RSE 3.98%, 95% CI 0.649-0.759)
    lec50 <- log(77.5); label("Steady-state trough concentration at half-maximal hazard reduction (ug/mL)") # Table 3 OS EC50 = 77.5 ug/mL (RSE 2.80%, 95% CI 73.247-81.753)
    lhill <- fixed(log(20)); label("Hill coefficient of the hazard-reduction sigmoid (unitless)") # Table 3 OS HILL = 20 with no RSE; Results: tested over 1-100 and fixed at 20

    # The control stream writes LAMDA = THETA(1)*EXP(ETA(1)) but Table 3
    # reports no variance and the Discussion states that typical predictions
    # were made, so the random effect is carried at zero.
    etallambda ~ fixed(0) # Supplementary control file LAMDA = THETA(1)*EXP(ETA(1)); OMEGA not reported
  })

  model({
    lambda <- exp(llambda + etallambda)
    gamma <- exp(lgamma)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Supplementary control file $PK: EFF = EMAX*(CONC**HILL)/((EC50**HILL)+(CONC**HILL))
    eff <- emax * CTROUGH^hill / (ec50^hill + CTROUGH^hill)

    # Supplementary control file $DES: the log-logistic hazard H1/H2 scaled by
    # (1 - EFF), with DEL = 1e-6 added to time so the hazard is finite at t = 0
    # when gamma < 1.
    del <- 1e-6
    tscaled <- (t + del) / lambda
    hazard <- (gamma / lambda) * tscaled^(gamma - 1) / (1 + tscaled^gamma) * (1 - eff)

    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0
    sur <- exp(-cumhaz)

    # Supplementary control file MEDTTE: solving
    # (1 + (t / lambda)^gamma)^-(1 - eff) = 0.5 for t.
    medtte <- lambda * (2^(1 / (1 - eff)) - 1)^(1 / gamma)
  })
}
