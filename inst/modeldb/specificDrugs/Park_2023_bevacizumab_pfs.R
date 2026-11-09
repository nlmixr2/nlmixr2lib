Park_2023_bevacizumab_pfs <- function() {
  description <- "Parametric log-logistic time-to-event exposure-response model for PROGRESSION-FREE SURVIVAL (blinded independent central review; progression or death) in adults with metastatic or recurrent non-squamous non-small-cell lung cancer receiving bevacizumab 15 mg/kg IV every 3 weeks (biosimilar SB8 or EU-sourced reference bevacizumab) with paclitaxel and carboplatin induction followed by bevacizumab maintenance (Park 2023, phase III study SB8-G31-NSCLC, n = 224). The baseline hazard is log-logistic, h0(t) = (gamma / lambda) * (t / lambda)^(gamma - 1) / (1 + (t / lambda)^gamma) with lambda in months, so the drug-free median progression-free survival equals lambda. The steady-state trough serum bevacizumab concentration CTROUGH (ug/mL) reduces the hazard by a constant fraction through an inhibitory sigmoid Emax term, h(t) = h0(t) * (1 - emax * CTROUGH^hill / (ec50^hill + CTROUGH^hill)), with the Hill coefficient fixed at 20 (a near all-or-none switch around ec50). Treatment (SB8 versus reference) was not a predictor and no other covariate was retained. There is no PK layer: CTROUGH is supplied as a time-fixed data column, the observed pre-dose concentration before cycle 7 (or cycle 5). Forward simulation exposes `hazard` (progression or death events per month), `sur` (progression-free survival probability) and `medtte` (the closed-form median progression-free survival time in months for the subject's CTROUGH). The control file published in the supplement is the overall-survival one; Park 2023 Equations 3-4 and Table 3 give the same log-logistic sigmoid-Emax structure for both endpoints, so this model shares the overall-survival equations with only the Table 3 estimates differing."
  reference <- paste(
    "Park S, Jung JA, Ju S, Lim HS.",
    "Exposure-response analysis using time-to-event data for bevacizumab",
    "biosimilar SB8 and the reference bevacizumab.",
    "Front Pharmacol. 2023;14:1288308. doi:10.3389/fphar.2023.1288308.",
    "Structural equations from the Supplementary Material NONMEM control",
    "file; parameter estimates from Table 3.",
    "Sister model from the same paper: modellib('Park_2023_bevacizumab_os').",
    sep = " "
  )
  vignette <- "Park_2023_bevacizumab"
  units <- list(
    time = "month",
    dosing = "n/a (no drug-dosing events; the exposure driver is supplied as the CTROUGH data covariate, in ug/mL)",
    concentration = "probability (the model output `sur` is the progression-free survival probability, not a drug concentration)"
  )

  compartmentData <- list(
    cumhaz = list(
      analyte = "cumulative hazard of progression or death",
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
      notes = "Screened and not retained for PFS (Park 2023 Table 2, dOFV -2.822)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened and not retained for PFS (Park 2023 Table 2, dOFV -3.315)."
    ),
    AGE = list(
      description = "Age, screened as < 65 versus >= 65 years",
      units = "years",
      type = "continuous",
      notes = "Screened and not retained (Park 2023 Table 2, dOFV -0.448 for PFS)."
    ),
    SEXF = list(
      description = "Sex",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained (Park 2023 Table 2, dOFV -2.985 for PFS)."
    ),
    ADA_POS = list(
      description = "Overall anti-drug antibody status, by cycle 7 and by end of treatment",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained (Park 2023 Table 2, dOFV -0.818 by cycle 7 and -0.722 by end of treatment, for PFS)."
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
    notes = "Subset of the phase III comparative study SB8-G31-NSCLC (EudraCT 2015-004026-34, NCT02754882; Reck 2020) restricted to patients with an available steady-state trough concentration. Progression-free survival was assessed by blinded independent central review; the visual predictive check spans 24 months (Park 2023 Figure 3). Fitted in NONMEM 7.4 by FOCE-I."
  )

  ini({
    # Park 2023 Table 3, progression-free-survival rows. The log-logistic scale
    # lambda is in months (Figures 3 and 6 are drawn in months, and the
    # drug-free median PFS of Figure 6, about 8.3 months, equals lambda).
    llambda <- log(8.28); label("Log-logistic scale lambda of the baseline hazard (month); equals the drug-free median progression-free survival time") # Table 3 PFS lambda = 8.28 (RSE 7.97%, 95% CI 6.986-9.574)
    lgamma <- log(3.94); label("Log-logistic shape gamma of the baseline hazard (unitless)") # Table 3 PFS gamma = 3.94 (RSE 9.19%, 95% CI 3.230-4.650)
    lemax <- log(0.451); label("Maximum fractional reduction of the hazard by bevacizumab (fraction)") # Table 3 PFS EMAX = 0.451 (RSE 28.82%, 95% CI 0.196-0.706)
    lec50 <- log(77.7); label("Steady-state trough concentration at half-maximal hazard reduction (ug/mL)") # Table 3 PFS EC50 = 77.7 ug/mL (RSE 5.08%, 95% CI 69.958-85.442)
    lhill <- fixed(log(20)); label("Hill coefficient of the hazard-reduction sigmoid (unitless)") # Table 3 PFS HILL = 20 with no RSE; Results: tested over 1-100 and fixed at 20

    # The control stream writes LAMDA = THETA(1)*EXP(ETA(1)) but Table 3
    # reports no variance and the Discussion states that typical predictions
    # were made, so the random effect is carried at zero.
    etallambda ~ fixed(0) # Supplementary (OS) control file LAMDA = THETA(1)*EXP(ETA(1)); OMEGA not reported
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
