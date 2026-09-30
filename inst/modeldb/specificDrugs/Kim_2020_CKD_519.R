Kim_2020_CKD_519 <- function() {
  description <- paste0(
    "Sequential population PK/PD model for the cholesteryl ester transfer ",
    "protein (CETP) inhibitor CKD-519 in healthy adult men given 50-400 mg ",
    "once daily with a standard meal for 14 days (Kim 2020). PK: ",
    "three-compartment disposition with an Erlang absorption chain (dose ",
    "into four transit compartments at Ktr, then an absorption compartment ",
    "draining to central at Ka) and a relative bioavailability that falls ",
    "with dose (Bmax * (1 - DOSE / (BA50 + DOSE))) and, for the 50-200 mg ",
    "cohorts only, decays exponentially with time since the first dose. ",
    "CETP activity is a turnover state whose first-order loss is stimulated ",
    "by an Emax function of plasma CKD-519 and whose production rises with ",
    "time on study (placebo Emax-in-time effect). HDL-C and LDL-C are ",
    "turnover states whose first-order elimination is respectively ",
    "stimulated (HDL-C) and inhibited (LDL-C) by a sigmoid function of CETP ",
    "activity; LDL-C production also rises linearly with time (placebo)."
  )
  reference <- paste0(
    "Kim CO, Jeon S, Han S, Park MS, Yim DS. A Population Pharmacokinetic ",
    "and Pharmacodynamic Model of CKD-519. Pharmaceutics. 2020;12(6):573. ",
    "doi:10.3390/pharmaceutics12060573. PK: Table 2 and Equations 1-3; ",
    "CETP activity: Table 3 and Equations 4-6; HDL-C and LDL-C: Table 4 and ",
    "Equations 7-10; model structure: Figure 1."
  )
  vignette <- "Kim_2020_CKD_519"

  # CETP activity is carried as a turnover state; it is not (yet) a
  # registered canonical compartment, so it is whitelisted for this paper.
  paper_specific_compartments <- c("cetp")
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL",
    cetp = "pmol (CETP activity assay units)",
    hdl = "mg/dL",
    ldl = "mg/dL"
  )

  covariateData <- list(
    DOSE = list(
      description = "Assigned once-daily CKD-519 dose of the subject's cohort",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Per-subject assigned dose (use case (a) of the register entry), ",
        "time-fixed; 50, 100, 200 or 400 mg in the source study. Enters the ",
        "relative bioavailability as BA = Bmax * (1 - DOSE / (BA50 + DOSE)) ",
        "(Kim 2020 Equation 2 and Table 2). Must equal the per-dose amount ",
        "in the event table. Set to any value for placebo subjects (who ",
        "receive no drug, so F never acts); 0 is conventional."
      ),
      source_name = "DOSE"
    ),
    DOSE_HIGH = list(
      description = "Indicator for the 400 mg cohort",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (50, 100 and 200 mg cohorts, and placebo)",
      notes = paste0(
        "1 = 400 mg once-daily cohort. Selects the time-decay constant of ",
        "the bioavailability term FT = exp(-alpha * TIME): alpha1 = 0.002 ",
        "1/h is estimated for the 50, 100 and 200 mg cohorts and alpha2 is ",
        "fixed to 0 for the 400 mg cohort (Kim 2020 Table 2), because the ",
        "400 mg group showed no loss of exposure on repeated dosing. The ",
        "paper defines the switch by dose group only; the model is ",
        "undefined for doses between 200 and 400 mg or above 400 mg."
      ),
      source_name = "400 mg DOSE group"
    )
  )

  compartmentData <- list(
    transit1 = list(analyte = "CKD-519", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "CKD-519", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "CKD-519", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "CKD-519", units = "mg", specimen = "administration site", verified = TRUE),
    depot = list(analyte = "CKD-519", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "CKD-519", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "CKD-519", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "CKD-519", units = "mg", specimen = "plasma", verified = TRUE),
    cetp = list(
      analyte = "cholesteryl ester transfer protein activity",
      units = "pmol",
      specimen = "plasma",
      verified = TRUE
    ),
    hdl = list(analyte = "HDL cholesterol", units = "mg/dL", specimen = "serum", verified = FALSE),
    ldl = list(analyte = "LDL cholesterol", units = "mg/dL", specimen = "serum", verified = FALSE)
  )

  population <- list(
    species = "human",
    n_subjects = 32L,
    n_studies = 1L,
    n_observations = "1392 CKD-519 concentrations, 2064 CETP activity values, 656 HDL-C and LDL-C values",
    age_range = "19-47 years (mean 32.2)",
    weight_range = "58.3-79.5 kg (mean 68.7)",
    sex_female_pct = 0,
    race_ethnicity = "Korean (single-centre study at Severance Hospital, Seoul)",
    disease_state = "Healthy male volunteers",
    dose_range = "CKD-519 50, 100, 200 or 400 mg (6 subjects each) or placebo (8 subjects) once daily with a standard breakfast (700-800 kcal, 5-25% fat) for 14 days",
    regions = "South Korea",
    notes = paste0(
      "Randomised, double-blind, placebo-controlled multiple-ascending-dose ",
      "study (NCT02753504 / NCT03210649); subjects were hospitalised from ",
      "day 1 to day 21. Demographics in Kim 2020 Table 1. Age, weight, BMI, ",
      "MDRD creatinine clearance and ALT were screened; serum creatinine on ",
      "CL/F entered at forward selection but was removed at backward ",
      "elimination, so the final models carry no demographic covariates."
    )
  )

  ini({
    # ---- PK (Kim 2020 Table 2) ----
    lcl <- log(6.4); label("Apparent clearance CL/F (L/h)") # Table 2: CL/F = 6.4 L/h (RSE 10.3%)
    lvc <- log(11.4); label("Apparent central volume V/F (L)") # Table 2: V/F = 11.4 L (RSE 11.0%)
    lvp <- log(45.4); label("Apparent peripheral volume 1 V2/F (L)") # Table 2: V2/F = 45.4 L (RSE 19.0%)
    lvp2 <- log(1006); label("Apparent peripheral volume 2 V3/F (L)") # Table 2: V3/F = 1006.0 L (RSE 32.4%)
    lq <- log(2.6); label("Apparent intercompartmental clearance 1 Q2/F (L/h)") # Table 2: Q2/F = 2.6 L/h (RSE 11.6%)
    lq2 <- log(3.3); label("Apparent intercompartmental clearance 2 Q3/F (L/h)") # Table 2: Q3/F = 3.3 L/h (RSE 15.4%)
    lka <- log(1.09); label("Absorption rate constant Ka (1/h)") # Table 2: Ka = 1.09 1/h (RSE 3.5%)
    lktr <- log(1.10); label("Transit rate constant Ktr (1/h)") # Table 2: Ktr = 1.10 1/h (RSE 3.4%)
    lfdepot_max <- log(1.6); label("Bmax: relative bioavailability extrapolated to zero dose (unitless)") # Table 2: Bmax = 1.6 (RSE 30.7%)
    ld50_fdepot <- log(90.1); label("BA50: dose halving the relative bioavailability (mg)") # Table 2: BA50 = 90.1 mg (RSE 24.5%)
    alpha_fdepot <- 0.002; label("alpha1: decay rate of relative bioavailability with time, 50-200 mg cohorts (1/h)") # Table 2: alpha1 = 0.002 (RSE 14.4%)
    alpha_fdepot_high <- fixed(0); label("alpha2: decay rate of relative bioavailability with time, 400 mg cohort (1/h)") # Table 2: alpha2 = 0 FIX

    # ---- CETP activity (Kim 2020 Table 3) ----
    lrbase_cetp <- log(350); label("Baseline CETP activity (pmol)") # Table 3: CETPbase = 350.0 pmole (RSE 1.9%)
    lkin_cetp <- log(164); label("Baseline zero-order production rate of CETP activity Kin,base (pmol/h)") # Table 3: Kinbase = 164.0 pmole/h (RSE 0.9%)
    lpbo_emax_cetp <- log(9.6); label("Kmax: maximal fractional placebo (time) increase in CETP production (unitless)") # Table 3: Kmax = 9.6 (RSE 2.2%)
    lpbo_t50_cetp <- log(9700); label("K50: time to half-maximal placebo increase in CETP production (h)") # Table 3: K50 = 9700.0 h (RSE 17.3%)
    lemax <- log(18.2); label("Emax: maximal fractional stimulation of CETP-activity loss by CKD-519 (unitless)") # Table 3: Emax = 18.2 (RSE 0.8%)
    lec50 <- log(587); label("EC50: CKD-519 concentration for half-maximal stimulation of CETP-activity loss (ng/mL)") # Table 3: EC50 = 587.0 ng/mL (RSE 3.7%)

    # ---- HDL-C (Kim 2020 Table 4, HDL-C block) ----
    lrbase_hdl <- log(50.0); label("Baseline HDL-C (mg/dL)") # Table 4 HDL-C: RB = 50.0 mg/dL (RSE 13.2%)
    lksyn_hdl <- log(1.26); label("Zero-order production rate of HDL-C (mg/dL/h)") # Table 4 HDL-C: Ksyn = 1.26 (printed unit 'h-1'; see vignette) (RSE 15.2%)
    lhill_hdl <- log(2.1); label("gamma: sigmoidicity of the CETP-activity effect on HDL-C elimination (unitless)") # Table 4 HDL-C: gamma = 2.1 (RSE 5.3%)
    lemax_hdl <- log(1.5); label("Rmax: maximal fractional stimulation of HDL-C elimination by CETP activity (unitless)") # Table 4 HDL-C: Rmax = 1.5 (RSE 56.5%)
    lec50_hdl <- log(185); label("R50: CETP activity giving half-maximal effect on HDL-C elimination (pmol)") # Table 4 HDL-C: R50 = 185.0 pmole (RSE 24.8%)

    # ---- LDL-C (Kim 2020 Table 4, LDL-C block) ----
    lrbase_ldl <- log(97.4); label("Baseline LDL-C (mg/dL)") # Table 4 LDL-C: RB = 97.4 mg/dL (RSE 1.9%)
    lksyn_ldl <- log(0.435); label("Zero-order production rate of LDL-C at time zero (mg/dL/h)") # Table 4 LDL-C: Ksyn = 0.435 (printed unit 'h-1'; see vignette) (RSE 1.5%)
    lpbo_slope_ldl <- log(0.0009); label("beta: linear placebo (time) increase in LDL-C production (1/h)") # Table 4 LDL-C: beta = 0.0009 (RSE 0.8%)
    lhill_ldl <- log(2.2); label("gamma: sigmoidicity of the CETP-activity effect on LDL-C elimination (unitless)") # Table 4 LDL-C: gamma = 2.2 (RSE 1.3%)
    lemax_ldl <- log(0.8); label("Rmax: maximal fractional inhibition of LDL-C elimination by CETP activity (unitless)") # Table 4 LDL-C: Rmax = 0.8 (RSE 1.1%)
    lec50_ldl <- log(80.7); label("R50: CETP activity giving half-maximal effect on LDL-C elimination (pmol)") # Table 4 LDL-C: R50 = 80.7 pmole (RSE 3.5%)

    # ---- IIV: exponential; table reports CV%, converted as omega^2 = log(CV^2 + 1) ----
    etalcl ~ 0.024045 # Table 2: omega CL/F = 15.6% CV
    etalvp2 ~ 0.290848 # Table 2: omega V3/F = 58.1% CV
    etalfdepot_max ~ 0.076520 # Table 2: omega Bmax = 28.2% CV
    etalktr ~ 0.019686 # Table 2: omega Ktr = 14.1% CV
    etalkin_cetp ~ 0.038075 # Table 3: omega Kinbase = 19.7% CV
    etalemax ~ 0.144987 # Table 3: omega Emax = 39.5% CV
    etalrbase_hdl ~ 0.029155 # Table 4 HDL-C: omega RB = 17.2% CV
    etalrbase_ldl ~ 0.054200 # Table 4 LDL-C: omega RB = 23.6% CV
    etalpbo_slope_ldl ~ 0.667150 # Table 4 LDL-C: omega beta = 97.4% CV
    etalhill_ldl ~ 0.230382 # Table 4 LDL-C: omega GAM = 50.9% CV

    # ---- Residual error (sigma rows read as standard deviations) ----
    propSd <- 0.30; label("Proportional residual error, CKD-519 (fraction)") # Table 2: sigma prop = 0.30 (RSE 4.2%)
    addSd <- fixed(0.0001); label("Additive residual error, CKD-519 (ng/mL)") # Table 2: sigma add = 0.0001 FIX
    addSd_cetp <- 40.4; label("Additive residual error, CETP activity (pmol)") # Table 3: sigma add = 40.4 (RSE 4.2%)
    propSd_cetp <- 0.112; label("Proportional residual error, CETP activity (fraction)") # Table 3: sigma prop = 0.112 (RSE 5.7%)
    addSd_hdl <- 4.6; label("Additive residual error, HDL-C (mg/dL)") # Table 4 HDL-C: sigma add = 4.6 (RSE 6.2%)
    propSd_hdl <- fixed(0.0001); label("Proportional residual error, HDL-C (fraction)") # Table 4 HDL-C: sigma prop printed '0.0.0001 FIX', read as 0.0001 FIX
    propSd_ldl <- 0.07; label("Proportional residual error, LDL-C (fraction)") # Table 4 LDL-C: sigma prop = 0.07 (RSE 1.8%); sigma add = 0 FIX so omitted
  })

  model({
    # ---- PK individual parameters ----
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    vp2 <- exp(lvp2 + etalvp2)
    q <- exp(lq)
    q2 <- exp(lq2)
    ka <- exp(lka)
    ktr <- exp(lktr + etalktr)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Relative bioavailability (Equations 1-3): F = BA * FT, with
    # BA = Bmax * (1 - DOSE / (BA50 + DOSE)) and FT = exp(-alpha * TIME).
    # TIME is time since the first dose (h); alpha is alpha1 for the
    # 50-200 mg cohorts and alpha2 (fixed 0) for the 400 mg cohort.
    fdepot_max <- exp(lfdepot_max + etalfdepot_max)
    d50_fdepot <- exp(ld50_fdepot)
    ba <- fdepot_max * (1 - DOSE / (d50_fdepot + DOSE))
    alpha <- alpha_fdepot * (1 - DOSE_HIGH) + alpha_fdepot_high * DOSE_HIGH
    ft <- exp(-alpha * t)
    fdepot <- ba * ft

    # ---- CETP activity (Equations 4-6) ----
    # Printed Equation 4 omits the state from the loss term; Kout is defined
    # in the text as a first-order rate constant, so the loss is
    # Kout * CETP * (1 + Drug). Kout = Kin,base / CETPbase keeps the
    # pre-dose state at CETPbase.
    rbase_cetp <- exp(lrbase_cetp)
    kin_cetp <- exp(lkin_cetp + etalkin_cetp)
    kout_cetp <- kin_cetp / rbase_cetp
    pbo_emax_cetp <- exp(lpbo_emax_cetp)
    pbo_t50_cetp <- exp(lpbo_t50_cetp)
    emax <- exp(lemax + etalemax)
    ec50 <- exp(lec50)

    # ---- HDL-C and LDL-C (Equations 7-10) ----
    rbase_hdl <- exp(lrbase_hdl + etalrbase_hdl)
    ksyn_hdl <- exp(lksyn_hdl)
    hill_hdl <- exp(lhill_hdl)
    emax_hdl <- exp(lemax_hdl)
    ec50_hdl <- exp(lec50_hdl)

    rbase_ldl <- exp(lrbase_ldl + etalrbase_ldl)
    ksyn_ldl <- exp(lksyn_ldl)
    pbo_slope_ldl <- exp(lpbo_slope_ldl + etalpbo_slope_ldl)
    hill_ldl <- exp(lhill_ldl + etalhill_ldl)
    emax_ldl <- exp(lemax_ldl)
    ec50_ldl <- exp(lec50_ldl)

    # The first-order elimination constants are not tabulated; they follow
    # from a pre-dose steady state at the baseline CETP activity, where the
    # CETP-activity effect (Equation 9) is already non-zero.
    res0_hdl <- emax_hdl * rbase_cetp^hill_hdl / (rbase_cetp^hill_hdl + ec50_hdl^hill_hdl)
    res0_ldl <- emax_ldl * rbase_cetp^hill_ldl / (rbase_cetp^hill_ldl + ec50_ldl^hill_ldl)
    kdeg_hdl <- ksyn_hdl / (rbase_hdl * (1 + res0_hdl))
    kdeg_ldl <- ksyn_ldl / (rbase_ldl * (1 - res0_ldl))

    # ---- ODEs ----
    d/dt(transit1) <- -ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(depot) <- ktr * transit4 - ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    f(transit1) <- fdepot

    # Plasma CKD-519 (mg/L -> ng/mL)
    Cc <- 1000 * central / vc

    pbo_cetp <- pbo_emax_cetp * t / (pbo_t50_cetp + t)
    drug_cetp <- emax * Cc / (ec50 + Cc)
    d/dt(cetp) <- kin_cetp * (1 + pbo_cetp) - kout_cetp * (1 + drug_cetp) * cetp
    cetp(0) <- rbase_cetp

    res_hdl <- emax_hdl * cetp^hill_hdl / (cetp^hill_hdl + ec50_hdl^hill_hdl)
    res_ldl <- emax_ldl * cetp^hill_ldl / (cetp^hill_ldl + ec50_ldl^hill_ldl)
    pbo_ldl <- pbo_slope_ldl * t

    d/dt(hdl) <- ksyn_hdl - kdeg_hdl * (1 + res_hdl) * hdl
    hdl(0) <- rbase_hdl
    d/dt(ldl) <- ksyn_ldl * (1 + pbo_ldl) - kdeg_ldl * (1 - res_ldl) * ldl
    ldl(0) <- rbase_ldl

    Cc ~ add(addSd) + prop(propSd)
    cetp ~ add(addSd_cetp) + prop(propSd_cetp)
    hdl ~ add(addSd_hdl) + prop(propSd_hdl)
    ldl ~ prop(propSd_ldl)
  })
}
