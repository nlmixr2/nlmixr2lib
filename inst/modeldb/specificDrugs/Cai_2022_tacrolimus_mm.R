Cai_2022_tacrolimus_mm <- function() {
  description <- "Empirical steady-state Michaelis-Menten population model for oral immediate-release tacrolimus in adult Chinese recipients early after first liver transplantation (Cai 2022 nonlinear MM final model). The model has no compartments and no dose events: it predicts the patient's total daily dose (mg/day) from the measured steady-state whole-blood trough supplied as the covariate CTROUGH, DD = Vm * C0 / ((10 / theta) * Km + C0) with theta = min(POD, 10) (Eq. 7). Vm = 6.62 mg/day; Km = 6.46 ng/mL with power effects of postoperative day, hematocrit and total bilirubin and scale effects of recipient CYP3A5 expresser status and triazole antifungal co-therapy. Residual error is on the daily-dose scale. Invert it to find the daily dose that achieves a target trough (the use the source paper makes of it)."
  reference <- "Cai XJ, Li RD, Li JH, Tao YF, Zhang QB, Shen CH, Zhang XF, Wang ZX, Jiao Z. Prospective population pharmacokinetic study of tacrolimus in adult recipients early after liver transplantation: A comparison of Michaelis-Menten and theory-based pharmacokinetic models. Front Pharmacol. 2022;13:1031969. doi:10.3389/fphar.2022.1031969"
  vignette <- "Cai_2022_tacrolimus"
  units <- list(
    time = "day",
    dosing = "mg/day (no dose events: the model output daily_dose is the predicted total daily dose)",
    concentration = "ng/mL (the input CTROUGH is the steady-state whole-blood trough)"
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Measured steady-state whole-blood tacrolimus trough concentration (C0)",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "(a) TOTAL concentration in WHOLE BLOOD (EMIT assay), not plasma. (b) Pre-morning-dose trough at steady state, taken after at least three doses at a constant dose rate (Cai 2022 Section 2.1). (c) Measured TDM value, not a model prediction; it is the independent variable of Eq. 7, which predicts the daily dose that sustains it. Table 1 median 5.96 ng/mL (range 0.30-23.11).",
      source_name = "C0"
    ),
    POD = list(
      description = "Postoperative day after liver transplantation",
      units = "days",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Enters twice (Section 2.5.2.2 Eq. 7 and Section 3.2.2.2 final-model equation): as the power effect (POD / 15)^0.277 on Km, and as the time factor 10 / theta multiplying Km in the dose equation, with theta = POD for 0 < POD <= 10 and 10 for POD > 10. POD must be > 0. The printed normaliser of the power term is 15 days, not the Table 1 median of 10 days; the printed equation is followed.",
      source_name = "POD"
    ),
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Power effect (HCT / 31)^1.16 on Km (Section 3.2.2.2 final-model equation; Table 1 median 31.3 percent, range 15.4-53.6).",
      source_name = "HCT"
    ),
    TBILI = list(
      description = "Total serum bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Power effect (TBIL / 58.25)^0.286 on Km (Section 3.2.2.2 final-model equation; 58.25 umol/L = Table 1 median, range 5.1-786.3).",
      source_name = "TBIL"
    ),
    CYP3A5_EXPR = list(
      description = "Recipient CYP3A5 expresser status (1 = *1/*1 or *1/*3; 0 = *3/*3)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A5*3/*3 nonexpresser)",
      notes = "Recipient genotype at rs776746. Km multiplied by 0.635 when CYP3A5_EXPR = 1 (Table 4 theta CYP3A5*1 = -0.365; Eq. 11 scale model).",
      source_name = "CYP3A5*1"
    ),
    CONMED_AZOLE = list(
      description = "Concomitant triazole antifungal (voriconazole or fluconazole) co-therapy",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no triazole antifungal)",
      notes = "Time-varying. Cai 2022 'TAF' (Table 4 row labelled theta VRCZ). Km multiplied by 3.15 when CONMED_AZOLE = 1 (Table 4 theta = 2.15; Eq. 11 scale model).",
      source_name = "TAF"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on Km and not significant (Section 3.2.2.2).",
      source_name = "WT"
    ),
    DOSE_TAC_MGD = list(
      description = "Total daily tacrolimus dose",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "The daily dose is the dependent variable of this model, not a covariate (Section 2.5.2: daily dose was screened only for the theory-based model).",
      source_name = "DD"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 176L,
    n_observations = 906L,
    n_studies = 1L,
    age_range = "18-74 years",
    age_median = "51 years (mean 50.68, SD 9.71)",
    weight_range = "40-104 kg",
    weight_median = "64 kg (mean 64.54, SD 11.85)",
    sex_female_pct = 14.8,
    race_ethnicity = "Chinese (single centre, Huashan Hospital, Shanghai)",
    disease_state = "Adult recipients of a first liver transplant from donation after cardiac death, early postoperative period (POD 2-72, median 10) until discharge; most recipients hepatitis B, C or E positive.",
    dose_range = "Immediate-release oral tacrolimus (Prograf) q12h, starting 0.5-1 mg q12h and titrated to whole-blood C0 8-12 ng/mL in the first 3 months; daily dose 0.25-8 mg/day (median 2.5).",
    regions = "China",
    notes = "Cai 2022 Table 1 (demographics) and Table 2 (genotypes). Same 176-patient, 906-trough dataset as the companion theory-based model Cai_2022_tacrolimus.R."
  )

  ini({
    lvmax <- log(6.62); label("Maximum steady-state daily dosing rate Vm (mg/day)") # Table 4 Vm = 6.62 mg/day (RSE 5.2%)
    lkm <- log(6.46); label("Steady-state trough at half of Vm, Km, at POD 15, HCT 31%, TBIL 58.25 umol/L, CYP3A5 nonexpresser, no triazole antifungal (ng/mL)") # Table 4 Km = 6.46 ng/mL (RSE 14%)

    e_pod_km <- 0.277; label("Power exponent of (POD / 15) on Km (unitless)") # Table 4 theta POD = 0.277 (RSE 22.6%)
    e_hct_km <- 1.16; label("Power exponent of (HCT / 31) on Km (unitless)") # Table 4 theta HCT = 1.16 (RSE 22.8%)
    e_tbili_km <- 0.286; label("Power exponent of (TBIL / 58.25) on Km (unitless)") # Table 4 theta TBIL = 0.286 (RSE 23%)
    e_cyp3a5_expr_km <- -0.365; label("Fractional change in Km for CYP3A5 expressers (unitless)") # Table 4 theta CYP3A5*1 = -0.365 (RSE 19.4%)
    e_conmed_azole_km <- 2.15; label("Fractional change in Km with triazole antifungal co-therapy (unitless)") # Table 4 theta VRCZ = 2.15 (RSE 42.3%)

    # Vm carries no IIV (fixed to 0, Section 3.2.2.1). Table 4 'omega Km' = 65.3
    # percent, read as CV percent of a log-normal eta: omega^2 = log(CV^2 + 1).
    etalkm ~ 0.3552 # Table 4 omega Km = 65.3 percent; log(0.653^2 + 1) = 0.3552

    propSd <- 0.204; label("Proportional residual error on the predicted daily dose (fraction)") # Table 4 delta1 = 20.4 percent (RSE 10.5%)
    addSd <- 0.568; label("Additive residual error on the predicted daily dose (mg/day)") # Table 4 delta2 = 0.568 mg/day (RSE 10.2%)
  })

  model({
    vmax <- exp(lvmax)
    km <- exp(lkm + etalkm) *
      (POD / 15)^e_pod_km *
      (HCT / 31)^e_hct_km *
      (TBILI / 58.25)^e_tbili_km *
      (1 + e_cyp3a5_expr_km * CYP3A5_EXPR) *
      (1 + e_conmed_azole_km * CONMED_AZOLE)

    # Eq. 7 time factor: theta = POD for 0 < POD <= 10, and 10 thereafter, so
    # the factor 10 / theta inflates Km in the first postoperative days.
    theta_pod <- POD
    if (POD > 10) theta_pod <- 10
    km_pod <- (10 / theta_pod) * km

    # Eq. 7: predicted total daily dose (mg/day) sustaining the trough CTROUGH.
    daily_dose <- vmax * CTROUGH / (km_pod + CTROUGH)

    daily_dose ~ add(addSd) + prop(propSd)
  })
}
