Zufferey_2018_fondaparinux <- function() {
  description <- "Parametric time-to-event model for major bleeding after major orthopaedic surgery under fondaparinux thromboprophylaxis (POP-A-RIX 2.5 mg once daily and PROPICE 1.5 mg once daily pooled cohorts; n = 1393, 64 adjudicated bleeding events). The hazard is hz(t) = h0(t) * exp(beta1*SEX + beta2*AUCinf/8.5 + beta3*LBM/44), with gamma-shaped baseline h0(t) = theta1*theta2*(t-theta3)*exp(-theta2*(t-theta3)) for t > theta3 and 0 otherwise (lag time theta3 ~= 17.6 h, peak ~4 days post-surgery). AUCinf is derived inside the model from daily dose and clearance using the paper's PK equation CL = 0.34 * (CRCL/60)^0.485 * exp(eta) (lean-body-weight Cockcroft-Gault CrCl)."
  reference <- paste(
    "Zufferey PJ, Ollier E, Delavenne X, Laporte S, Mismetti P, Duffull SB.",
    "Incidence and risk factors of major bleeding following major orthopaedic surgery",
    "with fondaparinux thromboprophylaxis. A time-to-event analysis.",
    "Br J Clin Pharmacol. 2018;84(10):2242-2251.",
    "doi:10.1111/bcp.13663.",
    sep = " "
  )
  vignette <- "Zufferey_2018_fondaparinux"
  units <- list(
    time          = "hour",
    dosing        = "mg (per-subject daily dose carried as the DOSE covariate; no rxode2 dose events)",
    concentration = "probability (the model output `sur` is the survival probability for avoiding a major bleeding event; AUCinf is reported in mg*h/L as a derived internal quantity)"
  )

  covariateData <- list(
    DOSE = list(
      description        = "Per-subject daily fondaparinux dose (mg/day). Time-fixed within a regimen; in the pooled cohort POP-A-RIX subjects received 2.5 mg once daily and PROPICE subjects received 1.5 mg once daily.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Carried as a per-record data column (use case (b) of the canonical DOSE entry). Used inside model() as AUCinf = DOSE / cl to reproduce the paper's procedure 'For each subject, AUCinf of fondaparinux was calculated by dividing the daily dose by the predicted clearance' (page 5).",
      source_name        = "DOSE"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance computed with lean body weight as the body-size descriptor (CrCl_LBW); raw mL/min, NOT BSA-normalised.",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column `CrCl_LBW`. Reference value 60 mL/min in the paper's clearance equation CL = 0.34 * (CrCl_LBW/60)^0.485 * exp(b) (page 5). Raw Cockcroft-Gault scale, following the `CLCR` raw-CrCl precedent (Delattre 2010 amikacin); the canonical CRCL register entry documents that BSA-normalised vs raw is paper-dependent and must be recorded in per-model notes. The pooled cohort range is 15-264 mL/min (Table 1).",
      source_name        = "CrCl_LBW"
    ),
    LBM = list(
      description        = "Lean body mass (kg), computed in the paper by the Janmahasatian et al. formula (Clin Pharmacokinet 2005;44:1051-65).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column `LBW` (lean body weight); same biological quantity as the canonical `LBM`. Reference value 44 kg in the hazard centring `LBM/44`, equal to the pooled-cohort mean (Table 1). Pooled-cohort range 26-93 kg.",
      source_name        = "LBW"
    ),
    SEXF = list(
      description        = "Biological sex (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Source column `SEX` is encoded 1 = male, 0 = female; the canonical SEXF inverts the values. The hazard coefficient `e_sexf_bleed` is applied as `e_sexf_bleed * (1 - SEXF)` inside model() to preserve the paper's reported positive value (1.62) corresponding to a higher bleeding hazard in males; reference category is 0 (male) per the canonical convention. Pooled-cohort female fraction 74% (Table 1).",
      source_name        = "SEX"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1393L,
    n_studies      = 2L,
    age_range      = "24-101 years (mean 76)",
    age_median     = "76 years",
    weight_range   = "35-172 kg (mean 68)",
    weight_median  = "68 kg",
    sex_female_pct = 74,
    race_ethnicity = NULL,
    disease_state  = "Adults undergoing major orthopaedic surgery (primary or revision hip arthroplasty, primary or revision knee arthroplasty, or hip fracture surgery) receiving fondaparinux thromboprophylaxis. The pooled cohort combines POP-A-RIX (n = 957 with CrCl > 30 mL/min; 2.5 mg fondaparinux once daily) and PROPICE (n = 436 with moderate renal impairment, CrCl 20-50 mL/min; 1.5 mg fondaparinux once daily).",
    dose_range     = "1.5 mg or 2.5 mg subcutaneously once daily, first dose at least 6 h postoperatively; recommended thromboprophylaxis duration 5 weeks",
    regions        = "France (two multicentre prospective open-label cohorts)",
    notes          = "64 adjudicated major-bleeding events (4.6% of pooled cohort; 5.2% by day 11). LBM mean 44 kg (range 26-93). Cockcroft-Gault CrCl with lean body weight (CrCl_LBW) mean 41 mL/min (range 10-173). Surgical mix in pooled cohort: 38% hip arthroplasty, 27% knee arthroplasty, 35% hip fracture. See Tables 1-2 of the source for full baseline demographics and event tallies. Source: Zufferey et al. 2018 (POP-A-RIX NCT01063543; PROPICE NCT00555438)."
  )

  ini({
    # Pharmacokinetic clearance (page 5, equation for CL):
    #   CL = 0.34 * (CrCl_LBW / 60)^0.485 * exp(b),  b ~ N(0, omega2)
    # Units: CL in L/h (so AUCinf = daily_dose / CL has units mg*h/L,
    # matching Table 1 AUCinf mean 8.5 mg*h/L for daily doses 1.5-2.5 mg).
    lcl       <- log(0.34);    label("Typical clearance at CRCL = 60 mL/min (L/h)")          # page 5, CL equation: TVCL = 0.34 L/h
    e_crcl_cl <- fixed(0.485); label("Power exponent on (CRCL / 60) for clearance (unitless)")# page 5, CL equation: exponent 0.485 (reported without uncertainty)

    # Baseline hazard (Table 3, Full hazard model hz(t)):
    #   h0(t) = theta1 * theta2 * (t - theta3) * exp(-theta2 * (t - theta3))   for t > theta3
    #   h0(t) = 0                                                              for t <= theta3
    # Gamma-shaped: 0 at t = theta3, peaks at t = theta3 + 1/theta2 (~100 h
    # ~= 4 days post-surgery), peak value = theta1 / e ~= 4.5e-4 1/h.
    # Functional form derived by differentiating the published closed-form
    # cumulative hazard Hz(t) = -(theta1/theta2) * [(theta2*(t-theta3) + 1)
    # * exp(-theta2*(t-theta3)) - 1] * exp(phi) (page 7); the (t - theta3)
    # zeroing for t < theta3 carries the Table 3 footnote convention. Time
    # unit is hours (Table 3 footnote a: "t, time since the end of surgery
    # (hours)").
    base_haz1 <- 0.00122; label("Baseline hazard scale parameter, theta1 (1/h)")        # Table 3 full model: theta1 = 0.00122 (95% CI 0.00002-0.00242)
    base_haz2 <- 0.0122;  label("Baseline hazard early-rate parameter, theta2 (1/h)")   # Table 3 full model: theta2 = 0.0122 (95% CI 0.0086-0.0158)
    t_lag     <- 17.6;    label("Baseline hazard lag time, theta3 (h)")                 # Table 3 full model: theta3 = 17.6 h (95% CI 9.23-25.9); zero before this time

    # Hazard log-linear coefficients (Table 3, Full hazard model hz(t)):
    #   phi = beta1*SEX + beta2*(AUCinf/8.5) + beta3*(LBW/44)
    # Reference values: AUCinf 8.5 mg*h/L (pooled-cohort mean) and LBW
    # 44 kg (pooled-cohort mean), per Table 1 and Table 3 footnote a.
    e_sexf_bleed <- 1.62;  label("Hazard log-coefficient on (1 - SEXF) ie male sex (unitless)") # Table 3: beta1*SEX = 1.62 (95% CI 0.76-2.48); applied as e_sexf_bleed * (1 - SEXF) to preserve paper sign (SEX = 1 if male, SEX = 0 if female)
    e_auc_bleed  <- 0.975; label("Hazard log-coefficient on AUCinf/8.5 (unitless)")             # Table 3: beta2*AUCinf/8.5 = 0.975 (95% CI 0.497-1.453)
    e_lbm_bleed  <- -1.93; label("Hazard log-coefficient on LBM/44 (unitless)")                 # Table 3: beta3*LBW/44 = -1.93 (95% CI -3.74 to -0.12)

    # Inter-individual variability on log-CL (page 5: "log normal between
    # subject variability of 34 (CV%)"). On the log-normal scale,
    # omega^2 = log(1 + CV^2) = log(1 + 0.34^2) = log(1.1156) = 0.10936.
    etalcl ~ 0.10936  # page 5 CL equation: CV% = 34 -> omega^2 = log(1 + 0.34^2) = 0.10936

    # No IIV on the hazard parameters. The source paper does not report
    # any random effect on theta1, theta2, theta3, beta1, beta2, beta3
    # (Table 3 reports only median bootstrap point estimates and 95% CIs,
    # not omegas). The TTE model is a typical-value hazard.

    # No residual error. Output is a survival probability / hazard rather
    # than a measured concentration; the source NONMEM run used a
    # parametric survival likelihood. `hazard`, `cumhaz`, and `sur` are
    # exposed as derived outputs for forward simulation.
  })
  model({
    # 1. Individual clearance (page 5 CL equation, in L/h)
    cl      <- exp(lcl + etalcl) * (CRCL / 60)^e_crcl_cl

    # 2. Derived AUCinf for a single dose at the current daily dose level
    # (page 5 simulation procedure: "AUCinf of fondaparinux was
    # calculated by dividing the daily dose by the predicted
    # clearance"). Units: DOSE (mg) / CL (L/h) = mg*h/L (matches Table 1
    # pooled-cohort mean 8.5 mg*h/L).
    auc_inf <- DOSE / cl

    # 3. Baseline hazard with lag (derived by differentiating the page 7
    # closed-form cumulative hazard; consistent with Figure 2 showing
    # h0 rising from 0 at the lag time, peaking near day 4, then decaying):
    #   h0(t) = theta1 * theta2 * (t - t_lag) * exp(-theta2 * (t - t_lag))  for t > t_lag
    #   h0(t) = 0                                                           for t <= t_lag
    # The `(t > t_lag) * (t - t_lag)` factor evaluates to 0/1 in rxode2 so
    # we get the piecewise zeroing without an `ifelse` branch.
    lag_term <- (t > t_lag) * (t - t_lag)
    h0       <- base_haz1 * base_haz2 * lag_term * exp(-base_haz2 * lag_term)

    # 4. Covariate-adjusted log-hazard. SEXF (canonical: 1 = female)
    # is inverted to (1 - SEXF) to recover the paper's SEX = 1 if male
    # encoding so e_sexf_bleed = +1.62 keeps the published sign.
    phi    <- e_sexf_bleed * (1 - SEXF) +
              e_auc_bleed  * (auc_inf / 8.5) +
              e_lbm_bleed  * (LBM     / 44)
    hazard <- h0 * exp(phi)

    # 5. Cumulative hazard and survival. The Zufferey 2018 closed-form
    # Hz(t) is published on page 7, but we use the ODE form so the
    # source-trace stays one-to-one with the baseline-hazard equation
    # above (no algebraic re-derivation between paper and code).
    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0
    sur          <- exp(-cumhaz)
  })
}
